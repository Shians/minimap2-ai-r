#include <Rcpp.h>
#include <stdint.h>  // Add this for integer types
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <zlib.h>
#include <chrono>
#include <sstream>

extern "C" {
    #include "minimap.h"
}

using namespace Rcpp;

// [[Rcpp::export]]
List align_sequences_cpp(
    const std::string& reference_file,
    const std::vector<std::string>& query_seqs,
    const std::vector<std::string>& query_names,
    const std::vector<std::string>& query_quals,
    const std::string& preset = "map-ont",
    int n_threads = 3
) {
    // Validate input lengths
    size_t n_seqs = query_seqs.size();
    if (n_seqs == 0) {
        stop("No query sequences provided");
    }
    
    if (query_names.size() != n_seqs) {
        stop("Number of query names (" + std::to_string(query_names.size()) + 
             ") does not match number of sequences (" + std::to_string(n_seqs) + ")");
    }
    
    if (query_quals.size() != n_seqs) {
        stop("Number of quality strings (" + std::to_string(query_quals.size()) + 
             ") does not match number of sequences (" + std::to_string(n_seqs) + ")");
    }
    
    Rcerr << "[minimap2] Initializing alignment with preset: " << preset << "\n";
    Rcerr << "[minimap2] Processing " << n_seqs << " sequences\n";
    
    // Initialize index options and mapping options
    mm_idxopt_t idx_opt;
    mm_mapopt_t map_opt;
    
    // Set preset parameters
    if (mm_set_opt(0, &idx_opt, &map_opt) < 0) {
        stop("Invalid preset option");
    }

    if (mm_set_opt(preset.c_str(), &idx_opt, &map_opt) < 0) {
        stop("Invalid preset option");
    }

    map_opt.flag |= MM_F_CIGAR;

    // Log indexing options
    Rcerr << "[minimap2] Index options:\n"
          << "  k-mer size (k): " << idx_opt.k << "\n"
          << "  window size (w): " << idx_opt.w << "\n"
          << "  bucket bits: " << idx_opt.bucket_bits << "\n"
          << "  mini batch size: " << idx_opt.mini_batch_size << "\n"
          << "  flag: " << std::hex << idx_opt.flag 
          << (idx_opt.flag & MM_I_HPC ? " (homopolymer-compressed)" : "")
          << (idx_opt.flag & MM_I_NO_SEQ ? " (no-seq)" : "")
          << (idx_opt.flag & MM_I_NO_NAME ? " (no-name)" : "")
          << std::dec << "\n";

    // Log mapping options
    Rcerr << "[minimap2] Mapping options:\n"
          << "  scoring: match=" << map_opt.a 
          << " mismatch=" << map_opt.b
          << " gap-open(Q)=" << map_opt.q 
          << " gap-extend(E)=" << map_opt.e << "\n"
          << "  min chain score: " << map_opt.min_chain_score << "\n"
          << "  min number of minimizers: " << map_opt.min_cnt << "\n"
          << "  bandwidth (dp): " << map_opt.bw << "\n"
          << "  bandwidth (long gaps): " << map_opt.bw_long << "\n"
          << "  max gap (ref): " << map_opt.max_gap_ref << "\n"
          << "  max gap (query): " << map_opt.max_gap << "\n"
          << "  flags: " << std::hex << map_opt.flag;

    // Print relevant mapping flags
    if (map_opt.flag & MM_F_CIGAR) Rcerr << " (cigar)";
    if (map_opt.flag & MM_F_OUT_SAM) Rcerr << " (sam)";
    if (map_opt.flag & MM_F_OUT_CS) Rcerr << " (cs)";
    if (map_opt.flag & MM_F_OUT_MD) Rcerr << " (md)";
    if (map_opt.flag & MM_F_SPLICE) Rcerr << " (splice)";
    if (map_opt.flag & MM_F_SR) Rcerr << " (short-read)";
    if (map_opt.flag & MM_F_FOR_ONLY) Rcerr << " (forward-only)";
    if (map_opt.flag & MM_F_REV_ONLY) Rcerr << " (reverse-only)";
    if (map_opt.flag & MM_F_SOFTCLIP) Rcerr << " (softclip)";
    Rcerr << std::dec << "\n";

    if (map_opt.best_n)
        Rcerr << "  max secondary alignments (best_n): " << map_opt.best_n << "\n";
    if (map_opt.max_qlen)
        Rcerr << "  max query length: " << map_opt.max_qlen << "\n";
    
    Rcerr << "[minimap2] Opening reference file: " << reference_file << "\n";
    
    // Open index reader
    mm_idx_reader_t *reader = mm_idx_reader_open(reference_file.c_str(), &idx_opt, 0);
    if (!reader) {
        stop("Failed to open reference file");
    }
    
    // Read the index
    mm_idx_t *mi = mm_idx_reader_read(reader, n_threads);
    mm_mapopt_update(&map_opt, mi);
    if (!mi) {
        mm_idx_reader_close(reader);
        stop("Failed to read index");
    }
    mm_idx_reader_close(reader);
    
    Rcerr << "[minimap2] Index loaded: " << mi->n_seq << " sequences, k=" << idx_opt.k 
          << ", w=" << idx_opt.w << "\n";
    
    // Initialize thread buffer
    mm_tbuf_t *tbuf = mm_tbuf_init();
    if (!tbuf) {
        mm_idx_destroy(mi);
        stop("Failed to initialize thread buffer");
    }
    
    // Store SAM output
    std::vector<std::string> sam_lines;
    
    // Store header lines
    sam_lines.push_back("@HD\tVN:1.6\tSO:unknown");
    for (int i = 0; i < (int)mi->n_seq; ++i) {
        sam_lines.push_back("@SQ\tSN:" + std::string(mi->seq[i].name) + 
                           "\tLN:" + std::to_string(mi->seq[i].len));
    }
    sam_lines.push_back("@PG\tID:minimap2\tPN:minimap2\tVN:2.24\tCL:minimap2");
    
    // Process sequences
    for (size_t seq_idx = 0; seq_idx < n_seqs; ++seq_idx) {
        const std::string& query_seq = query_seqs[seq_idx];
        std::string query_name = query_names[seq_idx];
        const std::string& query_qual = query_quals[seq_idx];
        
        // Strip leading '@' and everything after first whitespace from query name
        if (!query_name.empty() && query_name[0] == '@') {
            query_name = query_name.substr(1);
        }
        size_t space_pos = query_name.find_first_of(" \t");
        if (space_pos != std::string::npos) {
            query_name = query_name.substr(0, space_pos);
        }
        
        if (seq_idx > 0 && seq_idx % 1000 == 0) {
            Rcerr << "[minimap2] Processed " << seq_idx << " sequences\n";
            Rcpp::checkUserInterrupt();  // Allow user interruption
        }
        
        // Perform alignment
        int n_regs;
        mm_reg1_t *regs = mm_map(mi, 
                                query_seq.length(),
                                query_seq.c_str(),
                                &n_regs, 
                                tbuf,
                                &map_opt,
                                query_name.c_str());
        
        // Store alignment records
        if (n_regs == 0) {
            std::stringstream ss;
            ss << query_name << "\t4\t*\t0\t0\t*\t*\t0\t0\t"
               << query_seq << "\t" << query_qual;
            sam_lines.push_back(ss.str());
        } else {
            for (int i = 0; i < n_regs; ++i) {
                mm_reg1_t *r = &regs[i];
                
                // Get CIGAR string
                std::string cigar;
                if (r->p) {
                    uint32_t *cigar_array = r->p->cigar;
                    int n_cigar = r->p->n_cigar;
                    for (int j = 0; j < n_cigar; ++j) {
                        cigar += std::to_string(cigar_array[j] >> 4) + 
                                "MIDNSHP=XB"[cigar_array[j] & 0xf];
                    }
                }
                
                // Calculate flag
                int flag = 0;
                if (r->rev) flag |= 0x10;
                if (n_regs > 1) flag |= 0x100;
                
                // Print alignment
                std::stringstream ss;
                ss << query_name << "\t"
                   << flag << "\t"
                   << mi->seq[r->rid].name << "\t"
                   << r->rs + 1 << "\t"
                   << r->mapq << "\t"
                   << (cigar.empty() ? "*" : cigar) << "\t"
                   << "*\t0\t0\t"
                   << query_seq << "\t"
                   << query_qual << "\t"
                   << "NM:i:" << r->blen - r->mlen << "\t"
                   << "ms:i:" << r->score << "\t"
                   << "AS:i:" << r->score << "\t"
                   << "nn:i:" << r->mlen;
                sam_lines.push_back(ss.str());
            }
        }
        
        free(regs);
    }
    
    Rcerr << "[minimap2] Completed processing " << n_seqs << " sequences\n";
    
    // Cleanup
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);
    
    // Convert to R
    return wrap(sam_lines);
}

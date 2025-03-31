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

// Helper function to calculate SAM flags
uint16_t calculate_sam_flag(const mm_reg1_t* r) {
    uint16_t flag = 0;
    
    // 0x1   PAIRED        // not used in single-end alignment
    // 0x2   PROPER_PAIR   // set if proper_frag is true
    if (r->proper_frag) flag |= 0x2;
    
    // 0x4   UNMAP        // never set here as this is for mapped reads
    // 0x8   MUNMAP       // not used in single-end alignment
    
    // 0x10  REVERSE      // set if rev is true
    if (r->rev) flag |= 0x10;
    
    // 0x20  MREVERSE     // not used in single-end alignment
    
    // 0x40  READ1        // not used in single-end alignment
    // 0x80  READ2        // not used in single-end alignment
    
    // 0x100 SECONDARY    // set if not primary (parent != id)
    if (!r->sam_pri && (r->parent != r->id)) flag |= 0x100;
    
    // 0x200 QCFAIL       // not used
    
    // 0x400 DUP          // not used
    
    // 0x800 SUPPLEMENTARY // set if not primary and parent == id
    if (!r->sam_pri && (r->parent == r->id)) flag |= 0x800;
    
    return flag;
}

std::string generate_cigar(const mm_reg1_t *r, int query_length, const mm_mapopt_t &map_opt) {
    std::string cigar;
    
    // Determine clipping character based on flags
    char clip_char = 'S';  // default to soft clipping
    if (((!r->sam_pri) || (r->parent != r->id && (map_opt.flag & MM_F_SECONDARY_SEQ))) && 
        !(map_opt.flag & MM_F_SOFTCLIP)) {
        clip_char = 'H';
    }
    
    // Calculate clipping lengths based on strand orientation
    int start_clip, end_clip;
    if (r->rev) { // reverse strand
        start_clip = query_length - r->qe;
        end_clip = r->qs;
    } else { // forward strand
        start_clip = r->qs;
        end_clip = query_length - r->qe;
    }

    // Add clip at start if needed
    if (start_clip > 0) {
        cigar += std::to_string(start_clip) + clip_char;
    }
    
    // Add main alignment CIGAR
    if (r->p) {
        uint32_t *cigar_array = r->p->cigar;
        int n_cigar = r->p->n_cigar;
        for (int j = 0; j < n_cigar; ++j) {
            cigar += std::to_string(cigar_array[j] >> 4) + 
                    "MIDNSHP=XB"[cigar_array[j] & 0xf];
        }
    }
    
    // Add clip at end if needed
    if (end_clip > 0) {
        cigar += std::to_string(end_clip) + clip_char;
    }
    
    return cigar.empty() ? "*" : cigar;
}

// Add this helper function near the top of the file
std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.length(), 'N');
    for (size_t i = 0; i < seq.length(); ++i) {
        char base = seq[seq.length() - 1 - i];
        switch (base) {
            case 'A': rc[i] = 'T'; break;
            case 'T': rc[i] = 'A'; break;
            case 'G': rc[i] = 'C'; break;
            case 'C': rc[i] = 'G'; break;
            case 'a': rc[i] = 't'; break;
            case 't': rc[i] = 'a'; break;
            case 'g': rc[i] = 'c'; break;
            case 'c': rc[i] = 'g'; break;
            default: rc[i] = base;
        }
    }
    return rc;
}

// [[Rcpp::export]]
SEXP build_index_cpp(
    const std::string& reference_file,
    const std::string& preset,
    int n_threads
) {
    // Initialize index options and mapping options
    mm_idxopt_t idx_opt;
    mm_mapopt_t map_opt;
    // Set preset parameters
    if (mm_set_opt(0, &idx_opt, &map_opt) < 0) {
        stop("Failed to initialize options");
    }

    if (mm_set_opt(preset.c_str(), &idx_opt, &map_opt) < 0) {
        stop("Invalid preset option");
    }

    map_opt.flag |= MM_F_CIGAR;

    // Open index reader
    mm_idx_reader_t *reader = mm_idx_reader_open(reference_file.c_str(), &idx_opt, 0);
    if (!reader) {
        stop("Failed to open reference file");
    }
    
    // Read the index
    mm_idx_t *mi = mm_idx_reader_read(reader, n_threads);
    mm_mapopt_update(&map_opt, mi);
    mm_idx_reader_close(reader);
    
    if (!mi) {
        stop("Failed to read index");
    }

      // Create external pointer with custom deleter
    XPtr<mm_idx_t> idx_ptr(mi, true);
    return idx_ptr;
}

// [[Rcpp::export]]
std::vector<std::string> align_sequences_cpp(
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
    
    // Initialize thread buffer
    mm_tbuf_t *tbuf = mm_tbuf_init();
    if (!tbuf) {
        mm_idx_destroy(mi);
        stop("Failed to initialize thread buffer");
    }
    
    // Store SAM output
    std::vector<std::string> sam_lines;
    
    // Add header lines
    sam_lines.push_back("@HD\tVN:1.6\tSO:unknown");
    for (int i = 0; i < (int)mi->n_seq; ++i) {
        sam_lines.push_back("@SQ\tSN:" + std::string(mi->seq[i].name) + 
                           "\tLN:" + std::to_string(mi->seq[i].len));
    }
    sam_lines.push_back("@PG\tID:minimap2\tPN:minimap2\tVN:2.24\tCL:minimap2");
    
    // Process each sequence
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
            Rcpp::checkUserInterrupt();
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
        
        if (n_regs == 0) {
            // Store unmapped read
            std::stringstream ss;
            ss << query_name << "\t4\t*\t0\t0\t*\t*\t0\t0\t"
               << query_seq << "\t" << query_qual;
            sam_lines.push_back(ss.str());
        } else {
            // Store mapped reads
            for (int i = 0; i < n_regs; ++i) {
                mm_reg1_t *r = &regs[i];
                
                // Get CIGAR string and determine if hard clipping is used
                std::string cigar = generate_cigar(r, query_seq.length(), map_opt);
                bool is_hard_clipped = ((!r->sam_pri) || (r->parent != r->id && (map_opt.flag & MM_F_SECONDARY_SEQ))) && 
                                     !(map_opt.flag & MM_F_SOFTCLIP);

                // Trim sequence and quality if hard clipping is used
                std::string output_seq = query_seq;
                std::string output_qual = query_qual;
                if (is_hard_clipped && r->qe > r->qs) {
                    output_seq = query_seq.substr(r->qs, r->qe - r->qs);
                    if (!query_qual.empty()) {
                        output_qual = query_qual.substr(r->qs, r->qe - r->qs);
                    }
                }

                // Reverse complement sequence and reverse quality string if reverse strand
                if (r->rev) {
                    output_seq = reverse_complement(output_seq);
                    if (!output_qual.empty()) {
                        std::reverse(output_qual.begin(), output_qual.end());
                    }
                }

                // Calculate flag
                uint16_t flag = calculate_sam_flag(r);
                
                // Format alignment record
                std::stringstream ss;
                ss << query_name << "\t"
                   << flag << "\t"
                   << mi->seq[r->rid].name << "\t"
                   << r->rs + 1 << "\t"
                   << r->mapq << "\t"
                   << cigar << "\t"
                   << "*\t0\t0\t"
                   << output_seq << "\t"  // Use trimmed sequence
                   << output_qual << "\t"  // Use trimmed quality
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
    
    return sam_lines;
}

// [[Rcpp::export]]
std::vector<std::string> aligner_map_cpp(
    SEXP mi_ptr,
    const std::vector<std::string>& query_seqs,
    const std::vector<std::string>& query_names,
    const std::vector<std::string>& query_quals
) {
    // Cast SEXP to mm_idx_t*
    mm_idx_t* mi = Rcpp::as<XPtr<mm_idx_t>>(mi_ptr).get();

    mm_mapopt_t map_opt;
    mm_idxopt_t idx_opt;
    if (mm_set_opt(0, &idx_opt, &map_opt) < 0) {
        stop("Invalid preset option");
    }
    if (mm_set_opt("map-ont", &idx_opt, &map_opt) < 0) {
        stop("Invalid preset option");
    }
    map_opt.flag |= MM_F_CIGAR;

    mm_mapopt_update(&map_opt, mi);

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

    // Initialize thread buffer
    mm_tbuf_t *tbuf = mm_tbuf_init();
    if (!tbuf) {
        stop("Failed to initialize thread buffer");
    }
    
    // Store SAM output
    std::vector<std::string> sam_lines;
    
    // Add header lines
    sam_lines.push_back("@HD\tVN:1.6\tSO:unknown");
    for (int i = 0; i < (int)mi->n_seq; ++i) {
        sam_lines.push_back("@SQ\tSN:" + std::string(mi->seq[i].name) + 
                           "\tLN:" + std::to_string(mi->seq[i].len));
    }
    sam_lines.push_back("@PG\tID:minimap2\tPN:minimap2\tVN:2.24\tCL:minimap2");
    
    // Process each sequence
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
            Rcpp::checkUserInterrupt();
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
        
        if (n_regs == 0) {
            // Store unmapped read
            std::stringstream ss;
            ss << query_name << "\t4\t*\t0\t0\t*\t*\t0\t0\t"
               << query_seq << "\t" << query_qual;
            sam_lines.push_back(ss.str());
        } else {
            // Store mapped reads
            for (int i = 0; i < n_regs; ++i) {
                mm_reg1_t *r = &regs[i];
                
                // Get CIGAR string and determine if hard clipping is used
                std::string cigar = generate_cigar(r, query_seq.length(), map_opt);
                bool is_hard_clipped = ((!r->sam_pri) || (r->parent != r->id && (map_opt.flag & MM_F_SECONDARY_SEQ))) && 
                                     !(map_opt.flag & MM_F_SOFTCLIP);

                // Trim sequence and quality if hard clipping is used
                std::string output_seq = query_seq;
                std::string output_qual = query_qual;
                if (is_hard_clipped && r->qe > r->qs) {
                    output_seq = query_seq.substr(r->qs, r->qe - r->qs);
                    if (!query_qual.empty()) {
                        output_qual = query_qual.substr(r->qs, r->qe - r->qs);
                    }
                }

                // Reverse complement sequence and reverse quality string if reverse strand
                if (r->rev) {
                    output_seq = reverse_complement(output_seq);
                    if (!output_qual.empty()) {
                        std::reverse(output_qual.begin(), output_qual.end());
                    }
                }

                // Calculate flag
                uint16_t flag = calculate_sam_flag(r);
                
                // Format alignment record
                std::stringstream ss;
                ss << query_name << "\t"
                   << flag << "\t"
                   << mi->seq[r->rid].name << "\t"
                   << r->rs + 1 << "\t"
                   << r->mapq << "\t"
                   << cigar << "\t"
                   << "*\t0\t0\t"
                   << output_seq << "\t"
                   << output_qual << "\t"
                   << "NM:i:" << r->blen - r->mlen << "\t"
                   << "ms:i:" << r->score << "\t"
                   << "AS:i:" << r->score << "\t"
                   << "nn:i:" << r->mlen;
                sam_lines.push_back(ss.str());
            }
        }
        
        free(regs);
    }
    
    // Cleanup
    mm_tbuf_destroy(tbuf);
    
    return sam_lines;
}

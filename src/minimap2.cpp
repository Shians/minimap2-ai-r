#include <Rcpp.h>
#include <stdint.h>  // Add this for integer types
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <zlib.h>

extern "C" {
    #include "minimap.h"
}

using namespace Rcpp;

// [[Rcpp::export]]
SEXP build_index_cpp(
    const std::string& fn_idx_in,
    const std::string& preset,
    int k,
    int w,
    int min_cnt,
    int min_chain_score,
    int bw,
    int n_threads
) {
    mm_idxopt_t idx_opt;
    mm_mapopt_t map_opt;
    mm_idx_reader_t* reader;
    mm_idx_t* idx;

    // Initialize options with preset or defaults
    if (mm_set_opt(preset.empty() ? 0 : preset.c_str(), &idx_opt, &map_opt) < 0) {
        stop("Failed to set options");
    }

    // Override with user-specified parameters
    idx_opt.k = k;
    idx_opt.w = w;
    map_opt.min_cnt = min_cnt;
    map_opt.min_chain_score = min_chain_score;
    map_opt.bw = bw;

    // Open index reader
    reader = mm_idx_reader_open(fn_idx_in.c_str(), &idx_opt, 0);
    if (!reader) {
        stop("Failed to open input file");
    }

    // Read the index
    idx = mm_idx_reader_read(reader, n_threads);
    if (!idx) {
        mm_idx_reader_close(reader);
        stop("Failed to build index");
    }

    Rcout << "Index built successfully" << std::endl;
    Rcout << "Number of sequences: " << idx->n_seq << std::endl;

    mm_idx_seq_t* seq = idx->seq;
    for (int i = 0; i < idx->n_seq; i++) {
        Rcout << "Sequence " << i << ": " << seq[i].name << std::endl;
        Rcout << "Sequence length: " << seq[i].len << std::endl;
        Rcout << "Sequence offset: " << seq[i].offset << std::endl;
    }

    // Close the reader
    mm_idx_reader_close(reader);

    // Create external pointer with custom deleter
    XPtr<mm_idx_t> idx_ptr(idx, true);
    return idx_ptr;
}

// [[Rcpp::export]]
List align_sequences_cpp(SEXP idx_ptr,
                    const std::string& query_seq,
                    const std::string& query_name,
                    bool cs,
                    bool MD) {
    // Get index from external pointer
    XPtr<mm_idx_t> idx(idx_ptr);
    
    // Initialize thread buffer
    mm_tbuf_t *tbuf = mm_tbuf_init();
    if (!tbuf) {
        stop("Failed to initialize thread buffer");
    }
    
    // Set mapping options
    mm_idxopt_t idx_opt;
    mm_mapopt_t map_opt;
    
    // Initialize options with preset or defaults
    if (mm_set_opt(0, &idx_opt, &map_opt) < 0) {
        stop("Failed to set options");
    }
    
    // Set alignment flags
    map_opt.flag |= MM_F_CIGAR; // perform alignment
    // if (cs) map_opt.flag |= MM_F_OUT_CS;
    // if (MD) map_opt.flag |= MM_F_OUT_MD;
    
    try {
        // Perform alignment
        int n_regs;

        mm_idx_t *alignment_index = idx.get();

        mm_reg1_t *regs = mm_map(
            alignment_index,
            query_seq.length(),
            query_seq.c_str(),
            &n_regs,
            tbuf,
            &map_opt,
            query_name.c_str()
        );

        Rcout << "Number of alignments: " << n_regs << std::endl;
        
        // Process results
        List results;
        for (int i = 0; i < n_regs; ++i) {
            List hit;
            hit["query_name"] = String(query_name);
            hit["query_start"] = static_cast<int32_t>(regs[i].qs);
            hit["query_end"] = static_cast<int32_t>(regs[i].qe);
            hit["target_name"] = String(idx->seq[regs[i].rid].name);
            hit["target_start"] = static_cast<int32_t>(regs[i].rs);
            hit["target_end"] = static_cast<int32_t>(regs[i].re);
            hit["strand"] = static_cast<int32_t>(regs[i].rev ? -1 : 1);
            hit["mapq"] = static_cast<int32_t>(regs[i].mapq);
            hit["score"] = static_cast<int32_t>(regs[i].score);
            
            // Process CIGAR string
            if (regs[i].p) {
                std::string cigar;
                uint32_t *cigar_array = regs[i].p->cigar;
                int n_cigar = regs[i].p->n_cigar;
                for (int j = 0; j < n_cigar; ++j) {
                    cigar += std::to_string(cigar_array[j] >> 4) + 
                            "MIDNSHP=XB"[cigar_array[j] & 0xf];
                }
                hit["cigar"] = String(cigar);
                
                // Generate CS tag if requested
                if (cs) {
                    char *cs_str = nullptr;
                    int cs_max_len = 0;
                    int cs_len = mm_gen_cs(NULL, &cs_str, &cs_max_len, idx.get(), 
                                         &regs[i], query_seq.c_str(), 0);
                    if (cs_len > 0) {
                        hit["cs"] = std::string(cs_str, cs_len);
                    }
                    free(cs_str);
                }
                
                // Generate MD tag if requested
                if (MD) {
                    char *md_str = nullptr;
                    int md_max_len = 0;
                    int md_len = mm_gen_MD(NULL, &md_str, &md_max_len, idx.get(), 
                                         &regs[i], query_seq.c_str());
                    if (md_len > 0) {
                        hit["MD"] = std::string(md_str, md_len);
                    }
                    free(md_str);
                }
            }
            
            results.push_back(hit);
        }
        
        free(regs);  // Free the array of regions
        mm_tbuf_destroy(tbuf);  // Clean up thread buffer
        
        return results;
    } catch (...) {
        mm_tbuf_destroy(tbuf);  // Clean up in case of error
        throw;
    }
}

// [[Rcpp::export]]
void destroy_index_cpp(SEXP idx_ptr) {
    // The XPtr destructor will automatically call mm_idx_destroy
    XPtr<mm_idx_t> idx(idx_ptr);
}


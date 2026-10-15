#include <Rcpp.h>
#include <stdint.h>
#include <stdlib.h>
#include <algorithm>
#include <cctype>
#include <sstream>
#include <string>
#include <vector>

extern "C" {
    #include "minimap.h"
}

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Resource ownership
//
// Everything minimap2 allocates is released by a destructor, so an R error,
// a user interrupt (Rcpp::checkUserInterrupt throws) or any C++ exception
// unwinds without leaking.
// ---------------------------------------------------------------------------

// An index together with the mapping options derived from the preset it was
// built with. The R object returned by build_index() wraps one of these.
struct Minimap2Index {
    mm_idx_t *mi;
    mm_mapopt_t map_opt;
    std::string preset;

    Minimap2Index() : mi(nullptr) {}
    ~Minimap2Index() {
        if (mi) mm_idx_destroy(mi);
    }
    Minimap2Index(const Minimap2Index &) = delete;
    Minimap2Index &operator=(const Minimap2Index &) = delete;
};

// XPtr's default finalizer calls `delete`, which runs ~Minimap2Index and
// hence mm_idx_destroy().
typedef XPtr<Minimap2Index> Minimap2IndexPtr;

struct IndexReaderGuard {
    mm_idx_reader_t *reader;
    explicit IndexReaderGuard(mm_idx_reader_t *r) : reader(r) {}
    ~IndexReaderGuard() { if (reader) mm_idx_reader_close(reader); }
};

struct ThreadBufferGuard {
    mm_tbuf_t *tbuf;
    explicit ThreadBufferGuard(mm_tbuf_t *b) : tbuf(b) {}
    ~ThreadBufferGuard() { if (tbuf) mm_tbuf_destroy(tbuf); }
};

// mm_map() returns a malloc'ed array whose elements each own a malloc'ed
// `p` (CIGAR and DP details).
struct RegsGuard {
    mm_reg1_t *regs;
    int n_regs;
    RegsGuard(mm_reg1_t *r, int n) : regs(r), n_regs(n) {}
    ~RegsGuard() {
        if (!regs) return;
        for (int i = 0; i < n_regs; ++i) free(regs[i].p);
        free(regs);
    }
};

// Return the index behind an R external pointer, or stop() if the pointer
// is not an index or has been invalidated (e.g. by saveRDS()/readRDS() or
// by restoring a workspace, both of which leave a NULL address).
static Minimap2Index *get_index(SEXP idx_ptr) {
    if (TYPEOF(idx_ptr) != EXTPTRSXP) {
        stop("idx must be a minimap2 index created by build_index()");
    }
    Minimap2Index *index = static_cast<Minimap2Index *>(R_ExternalPtrAddr(idx_ptr));
    if (!index || !index->mi) {
        stop("index is no longer valid (it may have been saved and reloaded); "
             "rebuild it with build_index()");
    }
    return index;
}

// ---------------------------------------------------------------------------
// SAM formatting, following minimap2's mm_write_sam3() in format.c for
// single-end reads.
// ---------------------------------------------------------------------------

static int sam_flag(const mm_reg1_t *r) {
    int flag = 0;
    if (r->rev) flag |= 0x10;
    if (r->parent != r->id) flag |= 0x100;
    else if (!r->sam_pri) flag |= 0x800;
    // 0x2 (proper pair) only applies to multi-segment reads, never here.
    return flag;
}

// Supplementary alignments, and secondary alignments when SEQ is printed for
// them, are hard-clipped unless soft clipping was requested.
static bool uses_hard_clip(int flag, int64_t opt_flag) {
    return ((flag & 0x800) || ((flag & 0x100) && (opt_flag & MM_F_SECONDARY_SEQ))) &&
           !(opt_flag & MM_F_SOFTCLIP);
}

static std::string sam_cigar(const mm_reg1_t *r, int flag, int query_length, int64_t opt_flag) {
    if (!r->p) return "*";

    char clip_char = uses_hard_clip(flag, opt_flag) ? 'H' : 'S';
    int start_clip = r->rev ? query_length - r->qe : r->qs;
    int end_clip = r->rev ? r->qs : query_length - r->qe;

    std::string cigar;
    if (start_clip > 0) cigar += std::to_string(start_clip) + clip_char;
    for (uint32_t k = 0; k < r->p->n_cigar; ++k) {
        uint32_t op = r->p->cigar[k];
        cigar += std::to_string(op >> 4) + MM_CIGAR_STR[op & 0xf];
    }
    if (end_clip > 0) cigar += std::to_string(end_clip) + clip_char;
    return cigar.empty() ? "*" : cigar;
}

static char complement(char base) {
    switch (base) {
        case 'A': return 'T'; case 'T': return 'A';
        case 'G': return 'C'; case 'C': return 'G';
        case 'a': return 't'; case 't': return 'a';
        case 'g': return 'c'; case 'c': return 'g';
        default: return base;
    }
}

// Write seq[start, start + len), reversed (and complemented) when requested.
static std::string sam_field(const std::string &seq, int start, int len, bool rev, bool comp) {
    std::string out = seq.substr(start, len);
    if (rev) {
        std::reverse(out.begin(), out.end());
        if (comp) std::transform(out.begin(), out.end(), out.begin(), complement);
    }
    return out;
}

// SEQ and QUAL columns. `has_qual` is false when QUAL is absent ("*" or ""),
// in which case QUAL is written as "*" and never sliced.
static std::string sam_seq_qual(const mm_reg1_t *r, int flag, const std::string &seq,
                                const std::string &qual, bool has_qual, int64_t opt_flag) {
    int l_seq = (int)seq.length();
    int start = 0, len = l_seq;
    if (r && !((flag & 0x900) == 0 || (opt_flag & MM_F_SOFTCLIP))) {
        if ((flag & 0x100) && !(opt_flag & MM_F_SECONDARY_SEQ)) return "*\t*";
        start = r->qs;
        len = r->qe - r->qs;
    }
    bool rev = r && r->rev;
    std::string out = sam_field(seq, start, len, rev, true);
    out += '\t';
    out += has_qual ? sam_field(qual, start, len, rev, false) : "*";
    return out;
}

// NM/ms/AS/nn optional tags (dp_max0/dp_score/n_ambi live on r->p).
static std::string sam_extra_tags(const mm_reg1_t *r) {
    if (!r->p) return "";
    std::stringstream ss;
    ss << "NM:i:" << (r->blen - r->mlen + (int)r->p->n_ambi) << "\t"
       << "ms:i:" << r->p->dp_max0 << "\t"
       << "AS:i:" << r->p->dp_score << "\t"
       << "nn:i:" << (int)r->p->n_ambi;
    return ss.str();
}

static std::string clean_query_name(std::string name) {
    if (!name.empty() && name[0] == '@') name = name.substr(1);
    size_t space_pos = name.find_first_of(" \t");
    if (space_pos != std::string::npos) name = name.substr(0, space_pos);
    return name;
}

static void sam_header(std::vector<std::string> &sam_lines, const Minimap2Index &index) {
    const mm_idx_t *mi = index.mi;
    sam_lines.push_back("@HD\tVN:1.6\tSO:unsorted\tGO:query");
    for (uint32_t i = 0; i < mi->n_seq; ++i) {
        sam_lines.push_back("@SQ\tSN:" + std::string(mi->seq[i].name) +
                            "\tLN:" + std::to_string(mi->seq[i].len));
    }
    sam_lines.push_back(std::string("@PG\tID:minimap2\tPN:minimap2\tVN:") + MM_VERSION +
                        "\tCL:minimap2 -a -x " + index.preset);
}

// ---------------------------------------------------------------------------
// Exported functions
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
SEXP build_index_cpp(
    const std::string& reference_file,
    const std::string& preset,
    int n_threads
) {
    Minimap2IndexPtr index(new Minimap2Index(), true);
    mm_idxopt_t idx_opt;

    mm_set_opt(0, &idx_opt, &index->map_opt);
    if (mm_set_opt(preset.c_str(), &idx_opt, &index->map_opt) < 0) {
        stop("Invalid preset: '" + preset + "'. See the minimap2 manual (-x) for valid presets.");
    }
    index->preset = preset;
    index->map_opt.flag |= MM_F_OUT_SAM | MM_F_CIGAR;
    // Load the whole reference as a single index part; only one part is
    // ever mapped against.
    idx_opt.batch_size = 0x7fffffffffffffffULL;

    {
        IndexReaderGuard reader(mm_idx_reader_open(reference_file.c_str(), &idx_opt, 0));
        if (!reader.reader) stop("Failed to open reference file: " + reference_file);
        index->mi = mm_idx_reader_read(reader.reader, n_threads);
    }
    if (!index->mi) stop("Failed to read index from: " + reference_file);

    mm_mapopt_update(&index->map_opt, index->mi);
    return index;
}

// [[Rcpp::export]]
std::vector<std::string> aligner_map_cpp(
    SEXP mi_ptr,
    const std::vector<std::string>& query_seqs,
    const std::vector<std::string>& query_names,
    const std::vector<std::string>& query_quals
) {
    const Minimap2Index *index = get_index(mi_ptr);
    const mm_mapopt_t &map_opt = index->map_opt;

    size_t n_seqs = query_seqs.size();
    if (n_seqs == 0) stop("No query sequences provided");
    if (query_names.size() != n_seqs) {
        stop("Number of query names (" + std::to_string(query_names.size()) +
             ") does not match number of sequences (" + std::to_string(n_seqs) + ")");
    }
    if (query_quals.size() != n_seqs) {
        stop("Number of quality strings (" + std::to_string(query_quals.size()) +
             ") does not match number of sequences (" + std::to_string(n_seqs) + ")");
    }

    ThreadBufferGuard tbuf(mm_tbuf_init());
    if (!tbuf.tbuf) stop("Failed to initialize thread buffer");

    std::vector<std::string> sam_lines;
    sam_header(sam_lines, *index);

    for (size_t seq_idx = 0; seq_idx < n_seqs; ++seq_idx) {
        if (seq_idx > 0 && seq_idx % 1000 == 0) {
            Rcerr << "[minimap2] Processed " << seq_idx << " sequences\n";
            Rcpp::checkUserInterrupt();
        }

        const std::string &query_seq = query_seqs[seq_idx];
        const std::string &query_qual = query_quals[seq_idx];
        std::string query_name = clean_query_name(query_names[seq_idx]);
        bool has_qual = !(query_qual.empty() || query_qual == "*");
        if (has_qual && query_qual.length() != query_seq.length()) {
            stop("Quality string length (" + std::to_string(query_qual.length()) +
                 ") does not match sequence length (" + std::to_string(query_seq.length()) +
                 ") for query '" + query_name + "'");
        }

        int n_regs = 0;
        mm_reg1_t *regs_raw = mm_map(index->mi, (int)query_seq.length(), query_seq.c_str(),
                                     &n_regs, tbuf.tbuf, &map_opt, query_name.c_str());
        RegsGuard regs(regs_raw, n_regs);

        if (n_regs == 0) {
            sam_lines.push_back(query_name + "\t4\t*\t0\t0\t*\t*\t0\t0\t" +
                                sam_seq_qual(nullptr, 4, query_seq, query_qual, has_qual, map_opt.flag));
            continue;
        }

        for (int i = 0; i < n_regs; ++i) {
            const mm_reg1_t *r = &regs.regs[i];
            int flag = sam_flag(r);
            std::stringstream ss;
            ss << query_name << "\t"
               << flag << "\t"
               << index->mi->seq[r->rid].name << "\t"
               << r->rs + 1 << "\t"
               << r->mapq << "\t"
               << sam_cigar(r, flag, (int)query_seq.length(), map_opt.flag) << "\t"
               << "*\t0\t0\t"
               << sam_seq_qual(r, flag, query_seq, query_qual, has_qual, map_opt.flag) << "\t"
               << sam_extra_tags(r);
            sam_lines.push_back(ss.str());
        }
    }

    return sam_lines;
}

// Free an index now rather than waiting for garbage collection. The R object
// becomes invalid, and later use gives an R error.
// [[Rcpp::export]]
void index_destroy_cpp(SEXP mi_ptr) {
    if (TYPEOF(mi_ptr) != EXTPTRSXP) return;
    Minimap2Index *index = static_cast<Minimap2Index *>(R_ExternalPtrAddr(mi_ptr));
    if (!index) return;
    R_ClearExternalPtr(mi_ptr);
    delete index;
}

// [[Rcpp::export]]
List index_info_cpp(SEXP mi_ptr) {
    Minimap2Index *index = nullptr;
    if (TYPEOF(mi_ptr) == EXTPTRSXP) {
        index = static_cast<Minimap2Index *>(R_ExternalPtrAddr(mi_ptr));
    }
    if (!index || !index->mi) return List::create(_["valid"] = false);

    const mm_idx_t *mi = index->mi;
    uint64_t total_len = 0;
    for (uint32_t i = 0; i < mi->n_seq; ++i) total_len += mi->seq[i].len;
    return List::create(
        _["valid"] = true,
        _["preset"] = index->preset,
        _["n_seq"] = (int)mi->n_seq,
        _["total_length"] = (double)total_len,
        _["k"] = mi->k,
        _["w"] = mi->w,
        _["splice"] = (bool)(index->map_opt.flag & MM_F_SPLICE)
    );
}

// [[Rcpp::export]]
int calculate_effective_seq_size_cpp(const std::string& cigar) {
    if (cigar == "*") return 0;

    int length = 0;
    int current_num = 0;

    for (char c : cigar) {
        if (std::isdigit(c)) {
            current_num = current_num * 10 + (c - '0');
        } else {
            switch (c) {
                // Operations that consume query sequence
                case 'M':
                case 'I':
                case 'S':
                case '=':
                case 'X':
                    length += current_num;
                    break;
                // Operations that don't consume query sequence
                case 'D':
                case 'N':
                case 'H':
                case 'P':
                    break;
                default:
                    Rcpp::warning("Unknown CIGAR operation: %c", c);
            }
            current_num = 0;
        }
    }

    return length;
}

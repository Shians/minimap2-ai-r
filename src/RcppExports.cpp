// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// build_index_cpp
SEXP build_index_cpp(const std::string& fn_idx_in, const std::string& preset, int k, int w, int min_cnt, int min_chain_score, int bw, int n_threads);
RcppExport SEXP _minimap2_build_index_cpp(SEXP fn_idx_inSEXP, SEXP presetSEXP, SEXP kSEXP, SEXP wSEXP, SEXP min_cntSEXP, SEXP min_chain_scoreSEXP, SEXP bwSEXP, SEXP n_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type fn_idx_in(fn_idx_inSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type preset(presetSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type min_cnt(min_cntSEXP);
    Rcpp::traits::input_parameter< int >::type min_chain_score(min_chain_scoreSEXP);
    Rcpp::traits::input_parameter< int >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< int >::type n_threads(n_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(build_index_cpp(fn_idx_in, preset, k, w, min_cnt, min_chain_score, bw, n_threads));
    return rcpp_result_gen;
END_RCPP
}
// align_sequences_cpp
List align_sequences_cpp(SEXP idx_ptr, const std::string& query_seq, const std::string& query_name, bool cs, bool MD);
RcppExport SEXP _minimap2_align_sequences_cpp(SEXP idx_ptrSEXP, SEXP query_seqSEXP, SEXP query_nameSEXP, SEXP csSEXP, SEXP MDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type idx_ptr(idx_ptrSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type query_seq(query_seqSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type query_name(query_nameSEXP);
    Rcpp::traits::input_parameter< bool >::type cs(csSEXP);
    Rcpp::traits::input_parameter< bool >::type MD(MDSEXP);
    rcpp_result_gen = Rcpp::wrap(align_sequences_cpp(idx_ptr, query_seq, query_name, cs, MD));
    return rcpp_result_gen;
END_RCPP
}
// destroy_index_cpp
void destroy_index_cpp(SEXP idx_ptr);
RcppExport SEXP _minimap2_destroy_index_cpp(SEXP idx_ptrSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type idx_ptr(idx_ptrSEXP);
    destroy_index_cpp(idx_ptr);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_minimap2_build_index_cpp", (DL_FUNC) &_minimap2_build_index_cpp, 8},
    {"_minimap2_align_sequences_cpp", (DL_FUNC) &_minimap2_align_sequences_cpp, 5},
    {"_minimap2_destroy_index_cpp", (DL_FUNC) &_minimap2_destroy_index_cpp, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_minimap2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

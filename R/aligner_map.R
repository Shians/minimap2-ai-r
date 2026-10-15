#' Map sequences using a pre-built minimap2 index
#'
#' Maps each query against `idx` with the preset the index was built with.
#' Mapping runs on a single thread.
#'
#' @param idx A `minimap2_index` object created by [build_index()].
#' @param query_seqs Character vector of query sequences.
#' @param query_names Character vector of query names. Defaults to `"seq1"`,
#'   `"seq2"`, and so on.
#' @param query_quals Character vector of Phred+33 quality strings, each the
#'   same length as its sequence. `"*"` or `""` marks a query without
#'   qualities, for which QUAL is written as `*`. Defaults to no qualities.
#'
#' @return An [AlignmentRecords][AlignmentRecords-class] object holding the
#'   SAM header and one record per alignment (unmapped queries get one record
#'   with flag 4).
#' @export
aligner_map <- function(idx, query_seqs, query_names = NULL, query_quals = NULL) {
    if (!inherits(idx, "minimap2_index")) {
        stop("idx must be a minimap2 index created by build_index()")
    }

    if (!is.character(query_seqs)) {
        stop("query_seqs must be a character vector")
    }

    n_seqs <- length(query_seqs)

    if (is.null(query_names)) {
        query_names <- paste0("seq", seq_len(n_seqs))
    } else if (length(query_names) != n_seqs) {
        stop("Length of query_names must match length of query_seqs")
    }

    if (is.null(query_quals)) {
        query_quals <- rep("*", n_seqs)
    } else if (length(query_quals) != n_seqs) {
        stop("Length of query_quals must match length of query_seqs")
    }

    sam_to_records(aligner_map_cpp(idx, query_seqs, query_names, query_quals))
}

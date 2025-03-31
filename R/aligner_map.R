#' Map sequences using a pre-built minimap2 index
#'
#' @param idx An external pointer to a minimap2 index created by build_index()
#' @param query_seqs Character vector of query sequences
#' @param query_names Character vector of query names (optional)
#' @param query_quals Character vector of quality scores (optional)
#'
#' @return Character vector of SAM format alignment records
#' @export
aligner_map <- function(idx, query_seqs, query_names = NULL, query_quals = NULL) {
    # Input validation
    if (!is.character(query_seqs)) {
        stop("query_seqs must be a character vector")
    }
    
    n_seqs <- length(query_seqs)
    
    # Handle default query names
    if (is.null(query_names)) {
        query_names <- paste0("seq", seq_len(n_seqs))
    } else if (length(query_names) != n_seqs) {
        stop("Length of query_names must match length of query_seqs")
    }
    
    # Handle default quality scores
    if (is.null(query_quals)) {
        query_quals <- rep("*", n_seqs)
    } else if (length(query_quals) != n_seqs) {
        stop("Length of query_quals must match length of query_seqs")
    }
    
    # Call the C++ function
    sam_to_records(aligner_map_cpp(idx, query_seqs, query_names, query_quals))
} 
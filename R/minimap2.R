#' R Interface to Minimap2
#' @name minimap2
#' @keywords internal
"_PACKAGE"

#' @importFrom Rcpp sourceCpp
#' @useDynLib minimap2, .registration = TRUE
NULL

#' Calculate Effective Sequence Size from CIGAR String
#'
#' This function calculates the effective size of a sequence based on its CIGAR string.
#' The effective size includes all operations that consume the query sequence (M, I, S, =, X).
#'
#' @param cigar A character string containing the CIGAR string
#' @return An integer representing the effective sequence size
#' @examples
#' calculate_effective_seq_size("10M2I3M") # returns 15
#' calculate_effective_seq_size("5S10M5S") # returns 20
#' calculate_effective_seq_size("10M2D3M") # returns 13
#'
#' @export
calculate_effective_seq_size <- function(cigar) {
    if (!is.character(cigar)) {
        stop("CIGAR must be a character string")
    }
    calculate_effective_seq_size_cpp(cigar)
}

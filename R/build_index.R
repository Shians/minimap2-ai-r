#' Build a minimap2 index from a reference sequence
#'
#' Indexes a reference and records the preset's mapping options with it, so
#' that [aligner_map()] maps with the same preset the index was built for.
#'
#' The index lives in C memory and cannot be saved: an index restored with
#' [readRDS()] or from a saved workspace is invalid, and [aligner_map()]
#' signals an error for it. Rebuild the index in each session instead.
#'
#' @param reference_file Path to the reference sequence file (FASTA, optionally
#'   gzipped).
#' @param preset Character string naming a minimap2 preset (the value of
#'   `minimap2 -x`). Default is `"map-ont"`. Any preset accepted by the
#'   vendored minimap2 is valid; common ones include:
#'   \itemize{
#'     \item `"map-ont"`: Oxford Nanopore genomic reads
#'     \item `"map-pb"`: PacBio CLR genomic reads
#'     \item `"map-hifi"`, `"lr:hq"`: PacBio HiFi and Nanopore Q20 genomic reads
#'     \item `"sr"`: short genomic reads
#'     \item `"splice"`, `"splice:hq"`: long-read spliced alignment
#'     \item `"asm5"`, `"asm10"`, `"asm20"`: assembly to reference mapping
#'   }
#'   See the minimap2 manual for the full list.
#' @param n_threads Number of threads used to build the index. Default is 3.
#'   Mapping with [aligner_map()] is single-threaded regardless.
#'
#' @return A `minimap2_index` object: an external pointer to the index and its
#'   mapping options.
#' @export
#'
#' @examples
#' \dontrun{
#' # Build index for Oxford Nanopore reads
#' idx <- build_index("reference.fa", preset = "map-ont")
#'
#' # Build index for spliced alignment of long RNA reads with 8 threads
#' idx <- build_index("reference.fa", preset = "splice", n_threads = 8)
#' }
build_index <- function(reference_file, preset = "map-ont", n_threads = 3) {
    if (!is.character(reference_file) || length(reference_file) != 1) {
        stop("reference_file must be a single character string")
    }

    if (!file.exists(reference_file)) {
        stop("reference file does not exist: ", reference_file)
    }

    if (!is.character(preset) || length(preset) != 1 || is.na(preset)) {
        stop("preset must be a single character string")
    }

    if (!is.numeric(n_threads) || length(n_threads) != 1 ||
        n_threads < 1 || n_threads != round(n_threads)) {
        stop("n_threads must be a positive integer")
    }

    # minimap2's own mm_set_opt() validates the preset.
    idx <- build_index_cpp(reference_file, preset, as.integer(n_threads))
    class(idx) <- "minimap2_index"
    idx
}

#' @export
print.minimap2_index <- function(x, ...) {
    info <- index_info_cpp(x)
    if (!isTRUE(info$valid)) {
        cat("<minimap2_index: invalid; rebuild it with build_index()>\n")
    } else {
        cat("<minimap2_index>\n")
        cat("  preset:    ", info$preset, "\n", sep = "")
        cat("  sequences: ", info$n_seq, " (", format(info$total_length, big.mark = ","),
            " bp)\n", sep = "")
        cat("  k, w:      ", info$k, ", ", info$w, "\n", sep = "")
    }
    invisible(x)
}

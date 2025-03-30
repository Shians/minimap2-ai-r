#' R Interface to Minimap2
#' @name minimap2
#' @keywords internal
"_PACKAGE"

#' @importFrom Rcpp sourceCpp
#' @useDynLib minimap2, .registration = TRUE
NULL

#' Create a Minimap2 Aligner Object
#'
#' @param reference_file Path to the reference sequence file (FASTA/FASTQ)
#' @param preset Preset parameters for different types of sequences (e.g. "lr:hq", "splice:hq") see details
#' @param k K-mer length
#' @param w Minimizer window length
#' @param min_cnt Minimum number of minimizers
#' @param min_chain_score Minimum chaining score
#' @param bw Band width
#' @param n_threads Number of threads to use
#'
#' @details
#' The `preset` parameter is used to set the preset parameters for different types of sequences.
#' The available presets are:
#'  - "map-pb": PacBio CLR genomic reads (~15% error rate)
#'  - "map-ont": Oxford Nanopore reads (~20% error rate)
#'  - "map-hifi": PacBio HiFi/CCS genomic reads (>99% accuracy)
#'  - "lr:hq": Nanopore Q20 genomic reads
#'  - "sr": Short Illumina reads (~1% error rate)
#'  - "splice": Long-read spliced alignment (RNA-seq)
#'  - "splice:hq": Long-read spliced alignment (RNA-seq) with high quality
#'  - "asm5": Assembly to reference mapping (<5% divergence)
#'  - "asm10": Assembly to reference mapping (<10% divergence)
#'  - "ava-pb": PacBio all-vs-all overlap mapping
#'  - "ava-ont": ONT all-vs-all overlap mapping
#'
#' @return A Minimap2Aligner object
#' @export
Minimap2Aligner <- function(reference_file,
                          preset,
                          k = 15,
                          w = 10,
                          min_cnt = 3,
                          min_chain_score = 40,
                          bw = 500,
                          n_threads = 3) {
    if (!file.exists(reference_file)) {
        stop("Reference file not found: ", reference_file)
    }

    idx_ptr <- build_index_cpp(reference_file,
                          preset = preset,
                          k = k,
                          w = w,
                          min_cnt = min_cnt,
                          min_chain_score = min_chain_score,
                          bw = bw,
                          n_threads = n_threads)

    structure(list(
        idx_ptr = idx_ptr,
        reference_file = reference_file
    ), class = "Minimap2Aligner")
}

#' Align Sequences Using Minimap2
#'
#' @param aligner A Minimap2Aligner object
#' @param query_seq Query sequence to align
#' @param query_name Name of the query sequence
#' @param cs Whether to output the cs tag
#' @param MD Whether to output the MD tag
#' @return A list of alignments
#' @export
align_sequences <- function(aligner,
                          query_seq,
                          query_name = "",
                          cs = FALSE,
                          MD = FALSE) {
    if (!inherits(aligner, "Minimap2Aligner")) {
        stop("First argument must be a Minimap2Aligner object")
    }

    align_sequences_cpp(aligner$idx_ptr,
                   query_seq = query_seq,
                   query_name = query_name,
                   cs = cs,
                   MD = MD)
}

#' Clean Up Minimap2 Aligner Object
#'
#' @param aligner A Minimap2Aligner object
#' @export
destroy_aligner <- function(aligner) {
    if (!inherits(aligner, "Minimap2Aligner")) {
        stop("Argument must be a Minimap2Aligner object")
    }
    destroy_index_cpp(aligner$idx_ptr)
}

#' Print Minimap2 Aligner Object
#'
#' @param x A Minimap2Aligner object
#' @param ... Additional arguments passed to print
#' @export
print.Minimap2Aligner <- function(x, ...) {
    cat("Minimap2 Aligner\n")
    cat("Reference file:", x$reference_file, "\n")
    invisible(x)
}

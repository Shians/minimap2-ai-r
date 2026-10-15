#' Align sequences using minimap2
#'
#' Builds an index of `reference_file` with [build_index()], maps the queries
#' with [aligner_map()], and frees the index. To map several batches against
#' one reference, build the index once and call [aligner_map()] directly.
#'
#' @param reference_file Path to the reference sequence file (FASTA)
#' @param query_seqs Character vector of query sequences
#' @param query_names Character vector of query sequence names
#' @param query_quals Character vector of quality strings
#' @param preset Character string naming a minimap2 preset (the value of
#'   `minimap2 -x`), validated by minimap2 itself. Default is `"map-ont"`.
#'   See [build_index()] for common presets.
#' @param n_threads Number of threads used to build the index. Mapping is
#'   single-threaded.
#'
#' @return An AlignmentRecords object containing the alignment results
#'
#' @examples
#' \dontrun{
#' # Align ONT reads
#' records <- align_sequences(
#'   reference_file = "reference.fa",
#'   query_seqs = c("ACGT...", "TGCA..."),
#'   query_names = c("read1", "read2"),
#'   query_quals = c("!!!!", "!!!!"),
#'   preset = "map-ont"
#' )
#' 
#' # Check alignment results
#' n_records(records)
#' n_mapped(records)
#' head(records(records))
#' }
#'
#' @importFrom methods new
#' @export
align_sequences <- function(
    reference_file,
    query_seqs,
    query_names,
    query_quals,
    preset = "map-ont",
    n_threads = 3
) {
    # Input validation
    if (!file.exists(reference_file)) {
        stop("Reference file not found: ", reference_file)
    }
    
    if (length(query_seqs) == 0) {
        stop("No query sequences provided")
    }
    
    if (length(query_seqs) != length(query_names)) {
        stop("Number of query names (", length(query_names), 
             ") does not match number of sequences (", length(query_seqs), ")")
    }
    
    if (length(query_seqs) != length(query_quals)) {
        stop("Number of quality strings (", length(query_quals), 
             ") does not match number of sequences (", length(query_seqs), ")")
    }
    
    if (!is.numeric(n_threads) || length(n_threads) != 1 || n_threads < 1) {
        stop("Number of threads must be positive")
    }

    idx <- build_index(reference_file, preset = preset, n_threads = n_threads)
    on.exit(index_destroy_cpp(idx), add = TRUE)

    aligner_map(idx, query_seqs, query_names = query_names, query_quals = query_quals)
}

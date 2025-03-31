#' Align sequences using minimap2
#'
#' @param reference_file Path to the reference sequence file (FASTA)
#' @param query_seqs Character vector of query sequences
#' @param query_names Character vector of query sequence names
#' @param query_quals Character vector of quality strings
#' @param preset Preset parameters for different sequence types:
#'   \itemize{
#'     \item "map-pb": PacBio CLR genomic reads
#'     \item "map-ont": Oxford Nanopore reads
#'     \item "map-hifi": PacBio HiFi/CCS genomic reads
#'     \item "lr:hq": Nanopore Q20 genomic reads
#'     \item "sr": Short Illumina reads
#'     \item "splice": Long-read spliced alignment
#'     \item "splice:hq": High-quality long-read spliced alignment
#'     \item "asm5": Assembly to reference mapping (<5% divergence)
#'     \item "asm10": Assembly to reference mapping (<10% divergence)
#'   }
#' @param n_threads Number of threads to use for alignment
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
    
    # Check for valid preset
    valid_presets <- c(
        "map-pb", "map-ont", "map-hifi", "lr:hq", "sr",
        "splice", "splice:hq", "asm5", "asm10"
    )
    if (!preset %in% valid_presets) {
        stop("Invalid preset: ", preset, "\nValid presets are: ",
             paste(valid_presets, collapse = ", "))
    }
    
    # Check thread count
    if (n_threads < 1) {
        stop("Number of threads must be positive")
    }
    
    # Perform alignment
    sam_lines <- align_sequences_cpp(
        reference_file = reference_file,
        query_seqs = query_seqs,
        query_names = query_names,
        query_quals = query_quals,
        preset = preset,
        n_threads = n_threads
    )
    
    # Convert to AlignmentRecords object
    records <- sam_to_records(sam_lines)
    
    # Return results
    return(records)
} 
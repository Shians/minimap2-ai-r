#' Build a minimap2 index from a reference sequence
#'
#' @param reference_file Path to the reference sequence file (FASTA format)
#' @param preset Character string specifying the preset to use. Default is "map-ont".
#'   Available presets include:
#'   \itemize{
#'     \item "map-ont" - Oxford Nanopore read mapping
#'     \item "map-pb" - PacBio read mapping
#'     \item "map-hifi" - PacBio HiFi read mapping
#'     \item "sr" - Short-read mapping
#'     \item "asm5" - Long assembly to reference mapping (>95% identity)
#'     \item "asm10" - Long assembly to reference mapping (>85% identity)
#'     \item "splice" - Long-read spliced alignment
#'   }
#' @param n_threads Number of threads to use for index construction. Default is 3.
#'
#' @return An external pointer to the minimap2 index
#' @export
#'
#' @examples
#' \dontrun{
#' # Build index for Oxford Nanopore reads
#' idx <- build_index("reference.fa", preset = "map-ont")
#' 
#' # Build index for PacBio reads with 8 threads
#' idx <- build_index("reference.fa", preset = "map-pb", n_threads = 8)
#' }
build_index <- function(reference_file, preset = "map-ont", n_threads = 3) {
    # Input validation
    if (!is.character(reference_file) || length(reference_file) != 1) {
        stop("reference_file must be a single character string")
    }
    
    if (!file.exists(reference_file)) {
        stop("reference file does not exist: ", reference_file)
    }
    
    if (!is.character(preset) || length(preset) != 1) {
        stop("preset must be a single character string")
    }
    
    if (!is.numeric(n_threads) || length(n_threads) != 1 || 
        n_threads < 1 || n_threads != round(n_threads)) {
        stop("n_threads must be a positive integer")
    }
    
    # Valid presets
    valid_presets <- c("map-ont", "map-pb", "map-hifi", "sr", 
                      "asm5", "asm10", "splice")
    if (!preset %in% valid_presets) {
        stop("Invalid preset. Must be one of: ", 
             paste(valid_presets, collapse = ", "))
    }
    
    # Call the C++ function
    build_index_cpp(reference_file, preset, n_threads)
} 
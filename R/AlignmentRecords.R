#' AlignmentRecords Class
#'
#' @slot header Character vector storing SAM header lines
#' @slot records tibble containing alignment records
#' @importFrom tibble tibble
#' @importFrom methods new
#' @export
setClass("AlignmentRecords",
         slots = list(
             header = "character",
             records = "tbl_df"
         ))

#' Get header from alignment records
#'
#' @param x An object with alignment records
#' @return Character vector of header lines
#' @export
setGeneric("header", function(x) standardGeneric("header"))

#' @rdname header
#' @export
setMethod("header", "AlignmentRecords", function(x) {
    x@header
})

#' Get records from alignment records
#'
#' @param x An object with alignment records
#' @return tibble of alignment records
#' @export
setGeneric("records", function(x) standardGeneric("records"))

#' @rdname records
#' @export
setMethod("records", "AlignmentRecords", function(x) {
    x@records
})

#' Get number of records
#'
#' @param x An object with alignment records
#' @return Integer number of records
#' @export
setGeneric("n_records", function(x) standardGeneric("n_records"))

#' @rdname n_records
#' @export
setMethod("n_records", "AlignmentRecords", function(x) {
    nrow(x@records)
})

#' Get number of mapped records
#'
#' @param x An object with alignment records
#' @return Integer number of mapped records
#' @export
setGeneric("n_mapped", function(x) standardGeneric("n_mapped"))

#' @rdname n_mapped
#' @export
setMethod("n_mapped", "AlignmentRecords", function(x) {
    sum(!(x@records$flag & 4))  # Count records where unmapped bit (4) is not set
})

#' Get number of unmapped records
#'
#' @param x An object with alignment records
#' @return Integer number of unmapped records
#' @export
setGeneric("n_unmapped", function(x) standardGeneric("n_unmapped"))

#' @rdname n_unmapped
#' @export
setMethod("n_unmapped", "AlignmentRecords", function(x) {
    sum(x@records$flag & 4)  # Count records where unmapped bit (4) is set
})

#' Show method for AlignmentRecords
#'
#' @param object AlignmentRecords object
#' @export
setMethod("show", "AlignmentRecords", function(object) {
    cat("AlignmentRecords object\n")
    cat("Header lines:", length(header(object)), "\n")
    cat("Total records:", n_records(object), "\n")
    cat("Mapped records:", n_mapped(object), "\n")
    cat("Unmapped records:", n_unmapped(object), "\n")
    if (n_records(object) > 0) {
        print(head(records(object), n = 5))
        if (n_records(object) > 5) {
            cat("... with", n_records(object) - 5, "more records\n")
        }
    }
})

#' Create a new AlignmentRecords object
#'
#' @param header Character vector of SAM header lines
#' @param records tibble containing alignment records
#' @return An AlignmentRecords object
#' @export
AlignmentRecords <- function(header = character(), records = tibble::tibble(
    qname = character(),
    flag = integer(),
    rname = character(),
    pos = integer(),
    mapq = integer(),
    cigar = character(),
    rnext = character(),
    pnext = integer(),
    tlen = integer(),
    seq = character(),
    qual = character(),
    tags = character()
)) {
    methods::new("AlignmentRecords",
        header = header,
        records = records)
}

#' Convert SAM output to AlignmentRecords
#'
#' @param sam_lines Character vector of SAM format lines
#' @return AlignmentRecords object
#' @importFrom stringr str_detect str_split_fixed
#' @importFrom tibble as_tibble
#' @export
sam_to_records <- function(sam_lines) {
    # Split into header and alignment lines
    header_lines <- sam_lines[stringr::str_detect(sam_lines, "^@")]
    aln_lines <- sam_lines[!stringr::str_detect(sam_lines, "^@")]

    if (length(aln_lines) == 0) {
        return(AlignmentRecords(header = header_lines))
    }

    # Parse alignment lines
    fields <- stringr::str_split_fixed(aln_lines, "\t", 12)  # Changed to 12 columns

    # Create records tibble with mandatory fields
    records <- tibble::tibble(
        qname = fields[,1],
        flag = as.integer(fields[,2]),
        rname = fields[,3],
        pos = as.integer(fields[,4]),
        mapq = as.integer(fields[,5]),
        cigar = fields[,6],
        rnext = fields[,7],
        pnext = as.integer(fields[,8]),
        tlen = as.integer(fields[,9]),
        seq = fields[,10],
        qual = fields[,11],
        tags = fields[,12]
    )

    AlignmentRecords(header = header_lines, records = records)
}

random_dna <- function(n) {
    paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}

reverse_complement <- function(seq) {
    paste(rev(strsplit(chartr("ACGTacgt", "TGCAtgca", seq), "")[[1]]), collapse = "")
}

write_fasta <- function(seqs, path = tempfile(fileext = ".fa")) {
    writeLines(paste0(">", names(seqs), "\n", seqs), path)
    path
}

# Two random 20 kb contigs and a chimeric read joining a[1001:5000] to
# b[5001:7000]. Mapped with map-ont, the b part is a supplementary alignment
# with qs close to 4000 (> 1).
chimera_fixture <- function() {
    set.seed(1)
    contigs <- c(a = random_dna(20000), b = random_dna(20000))
    read <- paste0(substr(contigs[["a"]], 1001, 5000), substr(contigs[["b"]], 5001, 7000))
    list(ref = write_fasta(contigs), read = read)
}

flag_set <- function(flag, bit) bitwAnd(flag, bit) != 0

# Length of the leading hard clip of a CIGAR string (0 if none).
hard_clip_start <- function(cigar) {
    m <- regmatches(cigar, regexpr("^[0-9]+H", cigar))
    if (length(m) == 0) 0L else as.integer(sub("H", "", m))
}

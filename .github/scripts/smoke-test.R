#!/usr/bin/env Rscript
# Compiles the package against whatever is currently vendored in
# src/minimap2/ and runs a minimal end-to-end alignment. Used by the
# minimap2 sync workflow to verify a new upstream release didn't break
# the Rcpp wrapper (e.g. a changed struct field or function signature)
# before it gets pushed and tagged.
options(warn = 1)

pkgbuild::compile_dll(".", quiet = FALSE)
pkgload::load_all(".", quiet = TRUE, compile = FALSE)

ref <- tempfile(fileext = ".fa")
writeLines(c(">chr1", paste(rep("ACGTACGTAC", 20), collapse = "")), ref)

idx <- build_index(ref, preset = "map-ont", n_threads = 1)
stopifnot(!is.null(idx))

res <- align_sequences(
    reference_file = ref,
    query_seqs = paste(rep("ACGTACGTAC", 20), collapse = ""),
    query_names = "r1",
    query_quals = strrep("!", 200),
    preset = "map-ont",
    n_threads = 1
)

stopifnot(n_records(res) > 0)
stopifnot(n_mapped(res) > 0)
stopifnot(n_records(res) == n_mapped(res) + n_unmapped(res))

mapped_idx <- res@records$flag == 0
stopifnot(any(mapped_idx))
stopifnot(calculate_effective_seq_size(res@records$cigar[mapped_idx][1]) > 0)

idx2 <- build_index(ref, preset = "map-ont", n_threads = 1)
res2 <- aligner_map(idx2, query_seqs = "ACGTACGTAC")
stopifnot(n_records(res2) > 0)

cat("Smoke test passed: n_records =", n_records(res), " n_mapped =", n_mapped(res), "\n")

# --- Index object -----------------------------------------------------------

test_that("build_index returns a classed index that prints", {
    fx <- chimera_fixture()
    idx <- build_index(fx$ref, n_threads = 1)

    expect_s3_class(idx, "minimap2_index")
    expect_output(print(idx), "preset:    map-ont")
    expect_output(print(idx), "sequences: 2")
})

test_that("invalid presets are rejected by minimap2's option parser", {
    fx <- chimera_fixture()
    expect_error(build_index(fx$ref, preset = "no-such-preset"), "Invalid preset")
    expect_error(
        align_sequences(fx$ref, fx$read, "r1", "*", preset = "no-such-preset"),
        "Invalid preset"
    )
})

test_that("presets missing from the old hard-coded whitelist are accepted", {
    fx <- chimera_fixture()
    expect_s3_class(build_index(fx$ref, preset = "asm20", n_threads = 1), "minimap2_index")
    expect_s3_class(build_index(fx$ref, preset = "lr:hq", n_threads = 1), "minimap2_index")
})

test_that("an index restored with readRDS gives an R error, not a crash", {
    fx <- chimera_fixture()
    idx <- build_index(fx$ref, n_threads = 1)
    path <- tempfile(fileext = ".rds")
    saveRDS(idx, path)
    restored <- readRDS(path)

    expect_error(aligner_map(restored, fx$read), "no longer valid")
    expect_output(print(restored), "invalid")
    # The original index is unaffected.
    expect_gt(n_mapped(aligner_map(idx, fx$read)), 0)
})

test_that("aligner_map rejects objects that are not an index", {
    expect_error(aligner_map("ref.fa", "ACGT"), "created by build_index")
})

# --- Preset is honoured at mapping time --------------------------------------

test_that("a splice index maps a spliced read with an N operation", {
    set.seed(2)
    exon1 <- random_dna(600)
    exon2 <- random_dna(600)
    intron <- paste0("GT", random_dna(1996), "AG")
    ref <- write_fasta(c(gene = paste0(random_dna(1000), exon1, intron, exon2, random_dna(1000))))
    read <- paste0(exon1, exon2)

    spliced <- records(aligner_map(build_index(ref, preset = "splice", n_threads = 1), read))
    genomic <- records(aligner_map(build_index(ref, preset = "map-ont", n_threads = 1), read))

    primary <- spliced[!flag_set(spliced$flag, 0x900), ]
    expect_equal(nrow(primary), 1)
    expect_match(primary$cigar, "^600M2000N600M$")
    expect_false(any(grepl("N", genomic$cigar)))
})

# --- SAM fidelity ------------------------------------------------------------

test_that("default qualities work for a hard-clipped supplementary with qs > 1", {
    fx <- chimera_fixture()
    idx <- build_index(fx$ref, n_threads = 1)

    recs <- records(aligner_map(idx, fx$read))
    supp <- recs[flag_set(recs$flag, 0x800), ]

    expect_equal(nrow(recs), 2)
    expect_equal(nrow(supp), 1)
    expect_equal(supp$rname, "b")
    expect_match(supp$cigar, "^[0-9]+H[0-9]+M$")
    clip <- hard_clip_start(supp$cigar)
    expect_gt(clip, 1)
    expect_equal(supp$seq, substr(fx$read, clip + 1, nchar(fx$read)))
    expect_true(all(recs$qual == "*"))
})

test_that("empty quality strings are treated as absent", {
    fx <- chimera_fixture()
    idx <- build_index(fx$ref, n_threads = 1)
    recs <- records(aligner_map(idx, fx$read, query_quals = ""))
    expect_true(all(recs$qual == "*"))
})

test_that("real qualities are clipped and reversed with the sequence", {
    fx <- chimera_fixture()
    idx <- build_index(fx$ref, n_threads = 1)
    qual <- paste(sample(c("!", "5", "I"), nchar(fx$read), replace = TRUE), collapse = "")
    rc_read <- reverse_complement(fx$read)
    rc_qual <- paste(rev(strsplit(qual, "")[[1]]), collapse = "")

    recs <- records(aligner_map(idx, rc_read, query_quals = rc_qual))
    primary <- recs[!flag_set(recs$flag, 0x900), ]
    supp <- recs[flag_set(recs$flag, 0x800), ]

    # Both alignments are on the reverse strand, so SAM shows the original read.
    expect_true(all(flag_set(recs$flag, 0x10)))
    expect_equal(primary$seq, fx$read)
    expect_equal(primary$qual, qual)
    clip <- hard_clip_start(supp$cigar)
    expect_gt(clip, 1)
    expect_equal(supp$seq, substr(fx$read, clip + 1, nchar(fx$read)))
    expect_equal(supp$qual, substr(qual, clip + 1, nchar(qual)))
})

test_that("a quality string of the wrong length is an error", {
    fx <- chimera_fixture()
    idx <- build_index(fx$ref, n_threads = 1)
    expect_error(aligner_map(idx, fx$read, query_quals = "IIII"), "Quality string length")
})

test_that("secondary alignments are soft-clipped with SEQ and QUAL '*'", {
    set.seed(3)
    repeat_unit <- random_dna(5000)
    ref <- write_fasta(c(
        c1 = paste0(random_dna(3000), repeat_unit, random_dna(3000)),
        c2 = paste0(random_dna(3000), repeat_unit, random_dna(3000))
    ))
    read <- paste0(random_dna(200), substr(repeat_unit, 1001, 4000))
    qual <- strrep("I", nchar(read))

    recs <- records(aligner_map(build_index(ref, n_threads = 1), read, query_quals = qual))
    secondary <- recs[flag_set(recs$flag, 0x100), ]

    expect_gte(nrow(secondary), 1)
    expect_true(all(secondary$seq == "*"))
    expect_true(all(secondary$qual == "*"))
    expect_false(any(grepl("H", secondary$cigar)))
})

test_that("single-end records never carry the proper-pair flag", {
    fx <- chimera_fixture()
    idx <- build_index(fx$ref, n_threads = 1)
    recs <- records(aligner_map(idx, c(fx$read, random_dna(500))))
    expect_false(any(flag_set(recs$flag, 0x2)))
    expect_false(any(flag_set(recs$flag, 0x1)))
})

test_that("unmapped reads keep their sequence with flag 4", {
    fx <- chimera_fixture()
    idx <- build_index(fx$ref, n_threads = 1)
    unmapped_read <- random_dna(500)
    recs <- records(aligner_map(idx, unmapped_read, "u1"))
    expect_equal(recs$flag, 4L)
    expect_equal(recs$seq, unmapped_read)
    expect_equal(recs$qual, "*")
})

test_that("@PG reports the vendored minimap2 version", {
    fx <- chimera_fixture()
    hdr <- header(aligner_map(build_index(fx$ref, n_threads = 1), fx$read))
    pg <- grep("^@PG", hdr, value = TRUE)

    expect_length(pg, 1)
    expect_no_match(pg, "VN:2\\.24")
    expect_match(pg, "\tVN:2\\.[0-9]+-r[0-9]+")
    pkg_version <- as.character(utils::packageVersion("minimap2"))
    expect_match(pg, paste0("\tVN:", pkg_version), fixed = TRUE)
})

# --- align_sequences reuses the index + map path -----------------------------

test_that("align_sequences matches build_index + aligner_map", {
    fx <- chimera_fixture()
    one_shot <- align_sequences(fx$ref, fx$read, "r1", "*", n_threads = 1)
    two_step <- aligner_map(build_index(fx$ref, n_threads = 1), fx$read, "r1")

    expect_s4_class(one_shot, "AlignmentRecords")
    expect_identical(records(one_shot), records(two_step))
    expect_identical(header(one_shot), header(two_step))
})

test_that("align_sequences honours the preset", {
    fx <- chimera_fixture()
    one_shot <- align_sequences(fx$ref, fx$read, "r1", "*", preset = "map-hifi", n_threads = 1)
    expect_match(grep("^@PG", header(one_shot), value = TRUE), "-x map-hifi", fixed = TRUE)
})

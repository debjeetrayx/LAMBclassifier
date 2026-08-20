# The specification is a frozen scientific artefact. These tests assert that the
# shipped signatures, universe and thresholds are exactly the specified ones.

FROZEN_SIGNATURE_SHA256 <- "8626709901e967f04edfaaa6034e8d67c088da4ab6aa8e4e30e5a0d548989ffa"
FROZEN_UNIVERSE_SHA256  <- "81454b9a99665ed0ff17892eb66db028d67123d7b1aeb5779e39752694c3c7c0"
FROZEN_UNIVERSE_SIZE    <- 13468L
FROZEN_SPEC_VERSION     <- "1.0.0"

test_that("shipped specification matches the frozen 1.0.0 identifiers", {
    spec <- lamb_specification()
    expect_identical(spec$version, FROZEN_SPEC_VERSION)
    expect_identical(spec$signature_sha256, FROZEN_SIGNATURE_SHA256)
    expect_identical(spec$universe_sha256, FROZEN_UNIVERSE_SHA256)
    expect_identical(as.integer(spec$universe_size), FROZEN_UNIVERSE_SIZE)
})

test_that("METHODS.md documents the hashes the package actually computes", {
    # The documented hash must match what the package computes, so published
    # metrics always describe the classifier that is shipped.
    path <- system.file("METHODS.md", package = "LAMBclassifier")
    skip_if(path == "", "METHODS.md not installed")

    methods_text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    spec <- lamb_specification()

    expect_true(
        grepl(spec$signature_sha256, methods_text, fixed = TRUE),
        info = "Signature hash in METHODS.md does not match lamb_specification()."
    )
    expect_true(
        grepl(spec$universe_sha256, methods_text, fixed = TRUE),
        info = "Universe hash in METHODS.md does not match lamb_specification()."
    )
})

test_that("signature composition is exactly as frozen", {
    expected <- c(chondroid = 9L, squamous = 9L, basal_like = 14L,
                  luminal = 15L, apocrine = 14L, malo_her2 = 9L)
    observed <- vapply(default_gene_sets(), length, integer(1))
    expect_identical(observed[names(expected)], expected)
})

test_that("signature membership is as specified", {
    gs <- default_gene_sets()
    expect_true("SPDEF" %in% gs$luminal)
    expect_true("KRT5" %in% gs$squamous)
    expect_false("KRT5" %in% gs$basal_like)
    expect_false("AQP3" %in% gs$apocrine)
})

test_that("hash is sensitive to any signature change, in content or order", {
    base <- lamb_specification()

    reordered <- default_gene_sets()
    reordered$apocrine <- rev(reordered$apocrine)
    expect_false(identical(
        lamb_specification(gene_sets = reordered)$signature_sha256,
        base$signature_sha256
    ))

    added <- default_gene_sets()
    added$basal_like <- c(added$basal_like, "KRT5")
    expect_false(identical(
        lamb_specification(gene_sets = added)$signature_sha256,
        base$signature_sha256
    ))

    # universe is independent of signature edits
    expect_identical(
        lamb_specification(gene_sets = added)$universe_sha256,
        base$universe_sha256
    )
})

test_that("calibration thresholds are frozen", {
    cal <- lamb_calibration()
    expect_identical(cal$min_margin, 0.02)
    expect_identical(cal$min_gene_coverage, 0.6)
})

test_that("every signature gene is present in the reference universe", {
    universe <- lamb_universe()
    outside <- sort(unlist(lapply(default_gene_sets(), setdiff, y = universe),
                           use.names = FALSE))
    expect_identical(
        outside, character(0),
        info = paste("Signature genes outside the universe:",
                     paste(outside, collapse = ", "))
    )
})

test_that("every declared signature gene is also scorable", {
    universe <- lamb_universe()
    declared <- vapply(default_gene_sets(), length, integer(1))
    scorable <- vapply(default_gene_sets(),
                       function(g) length(intersect(g, universe)), integer(1))
    expect_identical(scorable, declared)
})

test_that("ROPN1B and its sibling ROPN1 are both in the universe", {
    universe <- lamb_universe()
    expect_true("ROPN1" %in% universe)
    expect_true("ROPN1B" %in% universe)
    expect_true("ROPN1B" %in% default_gene_sets()$basal_like)
})

test_that("specification version is documented in METHODS.md", {
    path <- system.file("METHODS.md", package = "LAMBclassifier")
    skip_if(path == "", "METHODS.md not installed")
    methods_text <- paste(readLines(path, warn = FALSE), collapse = "\n")
    expect_true(
        grepl(lamb_specification()$version, methods_text, fixed = TRUE),
        info = "Specification version in METHODS.md does not match the package."
    )
})

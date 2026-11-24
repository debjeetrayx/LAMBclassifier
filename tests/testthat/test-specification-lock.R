# Shipped signatures, universe and thresholds must equal the specification.

FROZEN_SIGNATURE_SHA256 <- "fe7c5dc714d7c6e72b1b70042076b6d067164d18f8b95648a299fe4c9be0336d"
FROZEN_UNIVERSE_SHA256  <- "81454b9a99665ed0ff17892eb66db028d67123d7b1aeb5779e39752694c3c7c0"
FROZEN_UNIVERSE_SIZE    <- 13468L
FROZEN_SPEC_VERSION     <- "2.0.0"

test_that("shipped specification matches the frozen 2.0.0 identifiers", {
    spec <- lamb_specification()
    expect_identical(spec$version, FROZEN_SPEC_VERSION)
    expect_identical(spec$signature_sha256, FROZEN_SIGNATURE_SHA256)
    expect_identical(spec$universe_sha256, FROZEN_UNIVERSE_SHA256)
    expect_identical(as.integer(spec$universe_size), FROZEN_UNIVERSE_SIZE)
})

test_that("METHODS.md documents the hashes the package actually computes", {
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
                  luminal = 15L, apocrine = 14L, malo = 15L)
    observed <- vapply(default_gene_sets(), length, integer(1))
    expect_identical(observed[names(expected)], expected)
})

test_that("signature membership is as specified", {
    gs <- default_gene_sets()
    expect_true("SPDEF" %in% gs$luminal)
    expect_true("KRT5" %in% gs$squamous)
    expect_false("KRT5" %in% gs$basal_like)
    expect_false("AQP3" %in% gs$apocrine)
    expect_length(intersect(gs$malo, DEFAULT_HER2_AMPLICON_GENES), 0L)
    expect_identical(gs$malo, c("TBX10", "BMP7", "TDO2", "SCGB2A1", "PHGDH", "ACY3", "RASD1", "CEACAM5",
                                "MAPK8IP1", "NDUFA4L2", "ACOX2", "DEFB1", "RNF145", "SCGB2A2", "DUSP6"))
    expect_identical(names(gs), c("chondroid", "squamous", "basal_like", "luminal", "apocrine", "malo"))
    expect_identical(DEFAULT_MALO_HER2_GENES, DEFAULT_MALO_GENES)   # deprecated alias
})

test_that("1.0.0 class gene lists are unchanged in 2.0.0", {
    gs <- default_gene_sets()
    expect_identical(gs$chondroid, c("SP7", "ACAN", "COL2A1", "PRELP", "TIMP4", "IBSP", "DLX5", "SOX11", "BMP2"))
    expect_identical(gs$squamous, c("BNC1", "ZBED2", "TP63", "POU3F1", "FOXQ1", "KRT6A", "KRT5", "KRT14", "DSG3"))
    expect_identical(gs$basal_like, c("FOXC1", "SOX10", "ELF5", "NFIB", "SOX11", "SFRP1", "BCL11A", "GABRP",
                                      "ID4", "VGLL1", "KLF5", "ROPN1", "ROPN1B", "EN1"))
    expect_identical(gs$luminal, c("ESR1", "GATA3", "MYB", "PGR", "RERG", "ANKRD30A", "STC2", "TFF1", "CA12",
                                   "INPP4B", "AGR3", "ANXA9", "BCL2", "FOXA1", "SPDEF"))
    expect_identical(gs$apocrine, c("KYNU", "SRD5A1", "CLDN8", "TSPAN8", "HPD", "ENPP3", "CUX2", "FKBP5",
                                    "UGT2B28", "CLCA2", "HPGD", "FOXA1", "AR", "SERHL2"))
})

test_that("HER2 amplicon rule is as specified", {
    cal <- lamb_calibration()
    expect_identical(cal$her2_amplicon$genes, c("ERBB2", "GRB7", "MIEN1", "PGAP3", "STARD3"))
    expect_identical(cal$her2_amplicon$genes, DEFAULT_HER2_AMPLICON_GENES)
    expect_identical(cal$her2_amplicon$q, 9)
    expect_identical(cal$her2_amplicon$min_high, 3L)
    expect_identical(cal$her2_amplicon$min_measured, 3L)
    expect_true(all(cal$her2_amplicon$genes %in% lamb_universe()))
    expect_identical(lamb_specification()$her2_amplicon$q, 9)
})

test_that("the amplicon module is covered by the signature hash", {
    cal <- lamb_calibration()
    cal$her2_amplicon$genes <- rev(cal$her2_amplicon$genes)
    expect_false(identical(lamb_specification(calibration = cal)$signature_sha256,
                           lamb_specification()$signature_sha256))
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

test_that("locked specification hashes are deterministic and sensitive", {
    a <- lamb_specification()
    b <- lamb_specification()
    expect_identical(a$signature_sha256, b$signature_sha256)
    expect_identical(a$universe_sha256, b$universe_sha256)
    expect_length(a$signature_sha256, 1)
    expect_length(a$universe_sha256, 1)

    changed <- default_gene_sets()
    changed$luminal <- rev(changed$luminal)
    c <- lamb_specification(gene_sets = changed)
    expect_false(identical(a$signature_sha256, c$signature_sha256))
    expect_identical(a$universe_sha256, c$universe_sha256)
})

test_that("binary metrics retain all four cells and Wilson intervals", {
    x <- lamb_binary_metrics(
        predicted = c("yes", "yes", "no", "no"),
        reference = c(TRUE, FALSE, TRUE, FALSE),
        positive_predicted = "yes", positive_reference = TRUE
    )
    expect_equal(unname(unlist(x[c("TP", "FP", "FN", "TN")])), rep(1, 4))
    expect_equal(x$Sensitivity, 0.5)
    expect_equal(x$Specificity, 0.5)
    expect_lt(x$Sensitivity_Lower, x$Sensitivity)
    expect_gt(x$Sensitivity_Upper, x$Sensitivity)
})

test_that("cross-classification normalization and Cramer's V are correct", {
    pred <- c("A", "A", "B", "B")
    ref <- c("X", "X", "Y", "Y")
    tab <- lamb_cross_tab(pred, ref, normalize = "prediction")
    expect_equal(rowSums(tab), c(A = 1, B = 1))
    expect_equal(lamb_cramers_v(pred, ref), 1)
})

test_that("result export includes locked provenance", {
    set.seed(20260821)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), 3000)))
    x <- matrix(
        seq_len(length(genes)), ncol = 1,
        dimnames = list(genes, "sample1")
    )
    result <- LAMBclassifier(x, run_pam50 = FALSE)
    path <- tempfile(fileext = ".csv")
    export_results(result, path, include_scores = FALSE)
    stem <- sub("\\.[^.]+$", "", path)
    expect_true(file.exists(path))
    expect_true(file.exists(paste0(stem, "_provenance.tsv")))
    expect_true(file.exists(paste0(stem, "_signatures.tsv")))
    exported <- read.csv(path, check.names = FALSE)
    expect_true(all(c("Classifier_Version", "Signature_SHA256", "Universe_SHA256") %in% names(exported)))
})

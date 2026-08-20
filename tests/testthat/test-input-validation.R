# Bad input must fail loudly rather than produce a confident wrong answer.

make_expr <- function(n_extra = 3000, n_samples = 2, seed = 20260820) {
    set.seed(seed)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), n_extra)))
    matrix(
        stats::runif(length(genes) * n_samples),
        nrow = length(genes),
        dimnames = list(genes, paste0("S", seq_len(n_samples)))
    )
}

test_that("missing rownames is an error", {
    x <- make_expr()
    rownames(x) <- NULL
    expect_error(lamb_score(x), "rownames")
})

test_that("empty dimensions are errors", {
    # A zero-row matrix has no rownames, so the rowname guard fires first.
    expect_error(
        lamb_score(matrix(numeric(0), nrow = 0, ncol = 2,
                          dimnames = list(character(0), c("A", "B")))),
        "rownames"
    )
    expect_error(
        lamb_score(matrix(numeric(0), nrow = 0, ncol = 2,
                          dimnames = list(NULL, c("A", "B")))),
        "rownames"
    )
    expect_error(
        lamb_score(matrix(numeric(0), nrow = 2, ncol = 0,
                          dimnames = list(c("G1", "G2"), character(0)))),
        "no samples"
    )
})

test_that("non-numeric data is rejected with a clear message", {
    x <- make_expr()
    storage.mode(x) <- "character"
    expect_error(lamb_score(x), "must be numeric")
})

test_that("an all-missing matrix is rejected rather than returning -Inf", {
    x <- make_expr()
    x[] <- NA_real_
    expect_error(lamb_score(x), "no non-missing values")
})

test_that("duplicated gene symbols warn and are de-duplicated, not double-counted", {
    x <- make_expr()
    rownames(x)[2] <- rownames(x)[1]
    expect_warning(scores <- lamb_score(x), "duplicated gene symbol")

    deduped <- x[!duplicated(rownames(x)), , drop = FALSE]
    expect_equal(scores, lamb_score(deduped), ignore_attr = TRUE)
})

test_that("untransformed data is detected and log2-transformed with a warning", {
    x <- make_expr() * 1000
    expect_warning(lamb_score(x), "looks untransformed")
})

test_that("classifier propagates the same input guards", {
    x <- make_expr()
    rownames(x) <- NULL
    expect_error(LAMBclassifier(x, run_pam50 = FALSE), "rownames")
})

test_that("data frames are accepted as well as matrices", {
    x <- make_expr()
    from_matrix <- LAMBclassifier(x, run_pam50 = FALSE)
    from_df <- LAMBclassifier(as.data.frame(x), run_pam50 = FALSE)
    expect_equal(from_matrix$Pred_Class, from_df$Pred_Class)
})

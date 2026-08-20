# A call implies every exclusion made to reach it. A gate that could not be
# evaluated must not be treated as though it had returned "no".

universe_matrix <- function(extra = 3000, seed = 20260821) {
    set.seed(seed)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), extra)))
    matrix(stats::runif(length(genes)), ncol = 1,
           dimnames = list(genes, "S1"))
}

drop_signature <- function(x, sig) {
    x[setdiff(rownames(x), default_gene_sets()[[sig]]), , drop = FALSE]
}

test_that("a complete matrix still yields a real call", {
    res <- suppressMessages(suppressWarnings(
        LAMBclassifier(universe_matrix(), run_pam50 = FALSE)
    ))
    expect_false(res$Pred_Class == "NC")
    expect_equal(res$Unscored_Signatures, "")
    expect_false(grepl("incomplete_gates", res$Flag))
})

test_that("dropping any required signature forces NC", {
    for (sig in c("chondroid", "squamous", "basal_like", "luminal", "apocrine")) {
        x <- drop_signature(universe_matrix(), sig)
        res <- suppressMessages(suppressWarnings(
            LAMBclassifier(x, run_pam50 = FALSE)
        ))
        expect_identical(res$Pred_Class, "NC", info = sig)
        expect_true(grepl("incomplete_gates", res$Flag), info = sig)
        expect_true(grepl(sig, res$Unscored_Signatures, fixed = TRUE), info = sig)
        expect_true(is.na(res$Margin), info = sig)
    }
})

test_that("dropping only MA-Lo genes does not force NC", {
    # MA-Lo refines a call inside gate 3; it excludes no class.
    x <- drop_signature(universe_matrix(), "malo_her2")
    res <- suppressMessages(suppressWarnings(
        LAMBclassifier(x, run_pam50 = FALSE)
    ))
    expect_false(res$Pred_Class == "NC")
    expect_false(grepl("incomplete_gates", res$Flag))
})

test_that("a panel with only luminal-lineage genes cannot produce a confident call", {
    gs <- default_gene_sets()
    keep <- unique(c(gs$luminal, gs$apocrine, sample(lamb_universe(), 500)))
    set.seed(7)
    x <- matrix(stats::runif(length(keep)), ncol = 1, dimnames = list(keep, "partial"))
    x[gs$luminal, 1] <- 5
    res <- suppressMessages(suppressWarnings(
        LAMBclassifier(x, run_pam50 = FALSE)
    ))
    expect_identical(res$Pred_Class, "NC")
    expect_true(grepl("incomplete_gates", res$Flag))
})

test_that("scoring is refused below 100 reference-universe genes", {
    gs <- default_gene_sets()
    keep <- unique(c(gs$luminal, gs$apocrine))
    set.seed(3)
    x <- matrix(stats::runif(length(keep)), ncol = 1, dimnames = list(keep, "tiny"))
    expect_error(lamb_score(x), "locked universe")
    expect_error(LAMBclassifier(x, run_pam50 = FALSE), "locked universe")
})

test_that("low coverage in a non-winning signature is still flagged", {
    x <- universe_matrix()
    gs <- default_gene_sets()
    # keep squamous scorable (3 genes) but well under 60% coverage
    x <- x[setdiff(rownames(x), gs$squamous[4:length(gs$squamous)]), , drop = FALSE]
    x[gs$luminal, 1] <- 5   # force an LHS-side call, not squamous
    res <- suppressMessages(suppressWarnings(
        LAMBclassifier(x, run_pam50 = FALSE)
    ))
    expect_false(res$Pred_Class == "Squamous")
    expect_true(grepl("low_coverage", res$Flag))
})

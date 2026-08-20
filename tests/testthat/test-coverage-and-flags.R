# Coverage handling and quality flags are what make targeted-panel use honest.

full_matrix <- function(n_extra = 3000, seed = 20260820) {
    set.seed(seed)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), n_extra)))
    matrix(stats::runif(length(genes)), ncol = 1, dimnames = list(genes, "S1"))
}

test_that("coverage attribute reports requested, found and fraction per signature", {
    scores <- suppressMessages(lamb_score(full_matrix()))
    cov <- attr(scores, "coverage")
    expect_true(all(c("Signature", "N_Requested", "N_Found", "Coverage") %in% names(cov)))
    expect_setequal(cov$Signature, names(default_gene_sets()))
    expect_true(all(cov$Coverage >= 0 & cov$Coverage <= 1))
    expect_true(all(cov$N_Found <= cov$N_Requested))
})

test_that("a signature below min_genes_per_sig is left unscored rather than guessed", {
    x <- full_matrix()
    # strip all but one chondroid gene
    gs <- default_gene_sets()
    drop <- setdiff(gs$chondroid, gs$chondroid[1])
    x <- x[setdiff(rownames(x), drop), , drop = FALSE]

    scores <- suppressMessages(suppressWarnings(
        lamb_score(x, min_genes_per_sig = 3)
    ))
    expect_true(is.na(scores[, "chondroid"]))
    expect_false(is.na(scores[, "luminal"]))
})

test_that("close calls are flagged when the margin is below the threshold", {
    # Construct near-tied luminal and apocrine scores by leaving data unboosted
    # and searching seeds for a small margin.
    found <- FALSE
    for (seed in 1:60) {
        res <- suppressMessages(suppressWarnings(
            LAMBclassifier(full_matrix(seed = seed), run_pam50 = FALSE)
        ))
        if (!is.na(res$Margin) && res$Margin < lamb_calibration()$min_margin) {
            expect_true(grepl("close_call", res$Flag))
            found <- TRUE
            break
        }
    }
    skip_if_not(found, "no close call arose in the searched seeds")
})

test_that("margin threshold governs the close_call flag consistently", {
    res <- suppressMessages(suppressWarnings(
        LAMBclassifier(full_matrix(), run_pam50 = FALSE)
    ))
    if (!is.na(res$Margin)) {
        expect_equal(
            grepl("close_call", res$Flag),
            res$Margin < lamb_calibration()$min_margin
        )
    }
})

test_that("targeted panels with enough universe genes still produce calls", {
    gs <- default_gene_sets()
    set.seed(3)
    panel <- unique(c(unlist(gs), sample(lamb_universe(), 400)))
    x <- matrix(stats::runif(length(panel)), ncol = 1,
                dimnames = list(panel, "spatial1"))
    res <- suppressMessages(suppressWarnings(
        LAMBclassifier(x, run_pam50 = FALSE)
    ))
    expect_equal(nrow(res), 1L)
    expect_true(res$Pred_Class %in% c(LAMB_CLASS_ORDER, "NC"))
})

test_that("a matrix sharing no genes with the universe is refused", {
    set.seed(5)
    x <- matrix(stats::runif(300), ncol = 1,
                dimnames = list(paste0("NOSUCHGENE", seq_len(300)), "S1"))
    expect_error(LAMBclassifier(x, run_pam50 = FALSE), "locked universe")
})

test_that("class order and colours are consistent and complete", {
    # Six substantive classes plus the "NC" non-call.
    expect_length(LAMB_CLASS_ORDER, 7L)
    expect_true("NC" %in% LAMB_CLASS_ORDER)
    expect_length(setdiff(LAMB_CLASS_ORDER, "NC"), 6L)
    expect_true(all(LAMB_CLASS_ORDER %in% names(LAMB_CLASS_COLORS)))
    expect_false(any(duplicated(LAMB_CLASS_ORDER)))
})

test_that("scoring is invariant to adding unrelated genes", {
    x <- full_matrix()
    set.seed(11)
    extra <- matrix(stats::runif(50), ncol = 1,
                    dimnames = list(paste0("UNRELATED", 1:50), "S1"))
    padded <- rbind(x, extra)
    a <- suppressMessages(suppressWarnings(LAMBclassifier(x, run_pam50 = FALSE)))
    b <- suppressMessages(suppressWarnings(LAMBclassifier(padded, run_pam50 = FALSE)))
    expect_identical(a$Pred_Class, b$Pred_Class)
})

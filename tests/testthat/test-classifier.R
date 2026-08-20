test_that("classification is invariant to processing a sample alone", {
    set.seed(20260820)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), 3000)))
    x <- matrix(
        stats::runif(length(genes) * 3), nrow = length(genes),
        dimnames = list(genes, paste0("S", 1:3))
    )
    together <- LAMBclassifier(x, run_pam50 = FALSE)
    alone <- LAMBclassifier(x[, "S2", drop = FALSE], run_pam50 = FALSE)
    cols <- c(
        "Pred_Class", "Margin", "Metaplastic_Score", "Score_Chondroid",
        "Score_Squamous", "Score_Basal_Like", "Score_LHS",
        "Score_Apocrine", "Score_MALO_HER2"
    )
    expect_equal(
        as.data.frame(together[together$Sample == "S2", cols]),
        as.data.frame(alone[, cols]),
        ignore_attr = TRUE
    )
})

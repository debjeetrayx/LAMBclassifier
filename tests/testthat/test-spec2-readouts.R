# 101 universe genes holding 1..101: value r has percentile 101 - r, so
# values 92..101 are the top 9%.
amp_sample <- function(module_values) {
    module <- DEFAULT_HER2_AMPLICON_GENES
    filler <- setdiff(lamb_universe(), c(module, DEFAULT_ALL_INPUT_GENES))[seq_len(101 - length(module))]
    genes <- c(module, filler)
    v <- setNames(numeric(length(genes)), genes)
    v[module] <- module_values
    v[filler] <- setdiff(1:101, module_values)
    v
}

amp_call <- function(v) {
    x <- matrix(v, ncol = 1, dimnames = list(names(v), "S1"))
    .lamb_s2_her2_amplicon(x)
}

test_that("HER2 amplicon reproduces a hand-computed example", {
    # Percentiles: ERBB2 0, GRB7 6, MIEN1 9 (on the cut-off), PGAP3 10, STARD3 100.
    v <- amp_sample(c(ERBB2 = 101, GRB7 = 95, MIEN1 = 92, PGAP3 = 91, STARD3 = 1))
    a <- amp_call(v)
    expect_identical(a$HER2_Amplicon, "high")
    expect_identical(a$N_High, 3L)
    expect_identical(a$N_Measured, 5L)

    # MIEN1 one rank lower: percentile 10, two genes reach the cut-off.
    v2 <- amp_sample(c(ERBB2 = 101, GRB7 = 95, MIEN1 = 90, PGAP3 = 91, STARD3 = 1))
    a2 <- amp_call(v2)
    expect_identical(a2$HER2_Amplicon, "not high")
    expect_identical(a2$N_High, 2L)
})

test_that("the percentile is taken within universe genes, without missing values", {
    v <- amp_sample(c(ERBB2 = 101, GRB7 = 95, MIEN1 = 92, PGAP3 = 91, STARD3 = 1))
    # Counted, a non-universe row would push MIEN1 to percentile 9.9.
    out <- c(v, NOT_A_UNIVERSE_GENE = 1000)
    expect_identical(amp_call(out)$HER2_Amplicon, "high")
    # A missing universe value is excluded from n, so n stays 101.
    na_gene <- setdiff(lamb_universe(), names(v))[1]
    expect_identical(amp_call(c(v, setNames(NA_real_, na_gene)))$HER2_Amplicon, "high")
    # Measured with a top value instead, it pushes MIEN1 to percentile 9.9.
    expect_identical(amp_call(c(v, setNames(200, na_gene)))$N_High, 2L)
})

test_that("ties take the average rank", {
    # MIEN1 tied with a filler at 93: average rank 92.5, percentile 8.5.
    v <- amp_sample(c(ERBB2 = 101, GRB7 = 95, MIEN1 = 93, PGAP3 = 91, STARD3 = 1))
    v[names(v)[v == 92]] <- 93
    a <- amp_call(v)
    expect_identical(a$N_High, 3L)
    expect_identical(a$HER2_Amplicon, "high")
})

test_that("fewer than 3 measured amplicon genes is not assessable", {
    v <- amp_sample(c(ERBB2 = 101, GRB7 = 100, MIEN1 = 99, PGAP3 = 98, STARD3 = 97))
    v[c("MIEN1", "PGAP3", "STARD3")] <- NA
    expect_identical(amp_call(v)$HER2_Amplicon, "not assessable")
    expect_identical(amp_call(v)$N_Measured, 2L)
    gone <- v[!names(v) %in% c("PGAP3", "STARD3")]
    gone["MIEN1"] <- 99
    expect_identical(amp_call(gone)$HER2_Amplicon, "high")   # 3 measured, all in the top 9%
})

# ---- end to end through LAMBclassifier ----

s2_matrix <- function(boost, seed = 20261120) {
    set.seed(seed)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), 3000)))
    x <- matrix(stats::runif(length(genes), 0, 1), ncol = 1, dimnames = list(genes, "S1"))
    gs <- c(default_gene_sets(), list(amplicon = DEFAULT_HER2_AMPLICON_GENES))
    for (sig in names(boost)) x[intersect(gs[[sig]], rownames(x)), 1] <- boost[[sig]]
    x
}

run_s2 <- function(x, ...) suppressMessages(suppressWarnings(LAMBclassifier(x, run_pam50 = FALSE, ...)))

test_that("HER2 amplicon is reported on every class and never changes it", {
    for (sig in c("chondroid", "squamous", "basal_like", "luminal", "apocrine", "malo")) {
        boost <- list(5)
        names(boost) <- sig
        if (sig == "malo") boost <- list(apocrine = 4, malo = 5)
        base <- run_s2(s2_matrix(boost))
        amp <- run_s2(s2_matrix(c(boost, list(amplicon = 10))))
        expect_identical(base$HER2_Amplicon, "not high", info = sig)
        expect_identical(amp$HER2_Amplicon, "high", info = sig)
        expect_identical(amp$Pred_Class, base$Pred_Class, info = sig)
    }
})

test_that("the deprecated malo_her2_genes argument still runs, with a warning", {
    x <- s2_matrix(list(apocrine = 4, malo = 5))
    expect_warning(res <- suppressMessages(LAMBclassifier(x, run_pam50 = FALSE, malo_her2_genes = DEFAULT_MALO_GENES)),
                   "deprecated")
    expect_identical(res$Pred_Class, run_s2(x)$Pred_Class)
})

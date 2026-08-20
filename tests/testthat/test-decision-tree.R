# Each gate of the decision tree must be reachable and must route correctly.

# Build a sample where the named signatures are top-ranked. Because singscore
# ranks within a sample, boosting a signature's genes drives its score up.
synth <- function(boost, n_extra = 3000, seed = 20260820) {
    set.seed(seed)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), n_extra)))
    x <- matrix(stats::runif(length(genes), 0, 1), ncol = 1,
                dimnames = list(genes, "S1"))
    gs <- default_gene_sets()
    for (i in seq_along(boost)) {
        sig <- names(boost)[i]
        x[intersect(gs[[sig]], rownames(x)), 1] <- boost[[i]]
    }
    x
}

call_of <- function(x) {
    suppressMessages(suppressWarnings(
        LAMBclassifier(x, run_pam50 = FALSE)$Pred_Class
    ))
}

test_that("gate 1 routes to Chondroid when chondroid leads", {
    expect_identical(call_of(synth(list(chondroid = 10))), "Chondroid")
})

test_that("gate 1 routes to Squamous when squamous leads", {
    expect_identical(call_of(synth(list(squamous = 10))), "Squamous")
})

test_that("gate 1 picks the higher of the two metaplastic signatures", {
    expect_identical(call_of(synth(list(chondroid = 10, squamous = 8))), "Chondroid")
    expect_identical(call_of(synth(list(chondroid = 8, squamous = 10))), "Squamous")
})

test_that("gate 2 routes to Basal-like only when it beats both luminal scores", {
    expect_identical(call_of(synth(list(basal_like = 10))), "Basal-like")
})

test_that("gate 3 separates LHS from Apocrine", {
    expect_identical(call_of(synth(list(luminal = 10))), "LHS")
    expect_identical(call_of(synth(list(apocrine = 10))), "Apocrine")
})

test_that("gate 4 routes to MA-Lo (HER2) when the amplicon signature leads", {
    # MA-Lo is only reachable inside gate 3, so the sample must first clear the
    # metaplastic and basal-like gates. Boosting a luminal-lineage signature
    # alongside malo_her2 does that, with malo_her2 highest.
    expect_identical(
        call_of(synth(list(apocrine = 8, malo_her2 = 10))),
        "MA-Lo (HER2)"
    )
})

test_that("metaplastic gate takes precedence over basal-like", {
    # Both boosted; step 1 must win because it is evaluated first.
    expect_identical(call_of(synth(list(squamous = 10, basal_like = 9))), "Squamous")
})

test_that("all six substantive classes are reachable", {
    calls <- c(
        call_of(synth(list(chondroid = 10))),
        call_of(synth(list(squamous = 10))),
        call_of(synth(list(basal_like = 10))),
        call_of(synth(list(luminal = 10))),
        call_of(synth(list(apocrine = 10))),
        call_of(synth(list(apocrine = 8, malo_her2 = 10)))
    )
    # LAMB_CLASS_ORDER also carries "NC", which is a non-call rather than a class.
    expect_setequal(unique(calls), setdiff(LAMB_CLASS_ORDER, "NC"))
})

test_that("margins are non-negative and finite for confident calls", {
    res <- suppressMessages(suppressWarnings(
        LAMBclassifier(synth(list(basal_like = 10)), run_pam50 = FALSE)
    ))
    expect_true(is.finite(res$Margin))
    expect_gte(res$Margin, 0)
})

test_that("gene row order does not affect the call", {
    x <- synth(list(apocrine = 10))
    set.seed(99)
    shuffled <- x[sample(nrow(x)), , drop = FALSE]
    expect_identical(call_of(x), call_of(shuffled))
})

test_that("results carry provenance columns for every sample", {
    # Same gene set, two columns, each pushed towards a different class.
    base <- synth(list())
    x <- cbind(base, base)
    colnames(x) <- c("A", "B")
    gs <- default_gene_sets()
    x[gs$luminal, "A"] <- 10
    x[gs$basal_like, "B"] <- 10
    res <- suppressMessages(suppressWarnings(LAMBclassifier(x, run_pam50 = FALSE)))
    expect_equal(nrow(res), 2L)
    expect_true(all(res$Signature_SHA256 == lamb_specification()$signature_sha256))
    expect_true(all(res$Classifier_Version == lamb_specification()$version))
})

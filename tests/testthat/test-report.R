pdf_pages <- function(path) {
    raw <- readBin(path, "raw", file.info(path)$size)
    txt <- rawToChar(raw[raw != as.raw(0)])
    m <- regmatches(txt, regexpr("/Type /Pages[^>]*/Count [0-9]+", txt, useBytes = TRUE))
    as.integer(sub(".*/Count ", "", m))
}

sim_results <- function(n = 4, seed = 11) {
    set.seed(seed)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), 2000)))
    x <- matrix(stats::rnorm(length(genes) * n, 6, 2), nrow = length(genes),
                dimnames = list(genes, paste0("S", seq_len(n))))
    sigs <- list(DEFAULT_LHS_GENES, DEFAULT_BASAL_LIKE_GENES, DEFAULT_SQUAMOUS_EPITHELIAL_GENES,
                 DEFAULT_MALO_GENES)
    for (i in seq_len(n)) {
        g <- intersect(sigs[[(i - 1) %% length(sigs) + 1]], genes)
        x[g, i] <- x[g, i] + 6
    }
    suppressWarnings(LAMBclassifier(x, run_pam50 = FALSE))
}

test_that(".lamb_confidence gives the three plain-language levels", {
    r <- data.frame(Pred_Class = c("LHS", "Basal-like", "NC", "Squamous"),
                    Margin = c(0.3, 0.01, NA, NA), Flag = c("", "close_call", "incomplete_gates", ""))
    expect_equal(.lamb_confidence(r), c("clear", "close call", "not classifiable", "not classifiable"))
    expect_equal(.lamb_confidence(r[0, ]), character())
})

test_that("multi-sample report has a cover page plus one page per sample", {
    res <- sim_results(4)
    f <- tempfile(fileext = ".pdf")
    expect_invisible(lamb_report(res, f, title = "Test cohort"))
    expect_true(file.exists(f))
    expect_equal(pdf_pages(f), nrow(res) + 1)
})

test_that("a single sample gives exactly one page", {
    res <- sim_results(1)
    f <- tempfile(fileext = ".pdf")
    lamb_report(res, f)
    expect_equal(pdf_pages(f), 1)
})

test_that("NC samples, NA scores, missing PAM50 and stripped attributes all render", {
    res <- sim_results(3)
    res$Pred_Class[2] <- "NC"
    res$Margin[2] <- NA
    res$Score_Squamous[2] <- NA
    res$Unscored_Signatures[2] <- "squamous"
    res$Score_MALO[3] <- NA
    expect_false("PAM50" %in% names(res))
    f <- tempfile(fileext = ".pdf")
    lamb_report(res, f)
    expect_equal(pdf_pages(f), 4)

    bare <- as.data.frame(res)[, c("Sample", "Pred_Class", "Margin")]   # no attributes, no scores
    f2 <- tempfile(fileext = ".pdf")
    lamb_report(bare, f2)
    expect_equal(pdf_pages(f2), 4)

    res$PAM50 <- c("LumA", NA, "Basal")
    f3 <- tempfile(fileext = ".pdf")
    lamb_report(res, f3)
    expect_equal(pdf_pages(f3), 4)
})

test_that("the decision path agrees with the call", {
    res <- sim_results(4)
    for (i in seq_len(nrow(res))) {
        p <- .lamb_path(as.list(as.data.frame(res)[i, ]))
        decided <- vapply(p$steps, `[[`, "", "status") == "decided"
        if (res$Pred_Class[i] == "NC") {
            expect_false(any(decided))
        } else {
            expect_equal(sum(decided), 1)
            expect_match(p$steps[[which(decided)]]$outcome, res$Pred_Class[i], fixed = TRUE)
        }
    }
    expect_true(all(c(LAMB_CLASS_ORDER) %in% names(.lamb_class_text)))
    expect_match(.lamb_class_text[["Squamous"]], "not a histological diagnosis", ignore.case = TRUE)
    expect_match(.lamb_class_text[["Chondroid"]], "not a histological diagnosis", ignore.case = TRUE)
    expect_match(.lamb_class_text[["Apocrine"]], "Not an ER/PR/AR test result")
    expect_match(.lamb_class_text[["MA-Lo"]], "Not an ER/PR/AR or HER2 test result")
    expect_false("MA-Lo (HER2)" %in% names(.lamb_class_text))
})

test_that("bad input is refused", {
    expect_error(lamb_report(data.frame(x = 1), tempfile()), "LAMBclassifier")
})

test_that("cluster report pages hold up to 12 clusters each", {
    res <- sim_results(14)
    res$Sample <- paste("cluster", seq_len(nrow(res)))
    calls <- data.frame(Cluster = c(res$Sample, "cluster 15"), Cells = c(rep(500, 14), 20),
                        Cell_type = c(rep("epithelial", 14), "not epithelial"),
                        Epithelial_share = 0.9, LAMB_class = c(res$Pred_Class, "LHS"), Note = "")
    f <- tempfile(fileext = ".pdf")
    lamb_report_cells(calls, res, f)
    expect_equal(pdf_pages(f), 2)
    f1 <- tempfile(fileext = ".pdf")
    lamb_report_cells(calls[1:3, ], res[1:3, ], f1)
    expect_equal(pdf_pages(f1), 1)
})

test_that("low gene coverage is named on reports and in the notes", {
    r <- data.frame(Sample = "s1", Pred_Class = "LHS", Margin = 0.3, Flag = "low_coverage")
    attr(r, "coverage") <- data.frame(Signature = c("luminal", "squamous", "malo"),
                                      N_Requested = c(15, 9, 15), N_Found = c(15, 4, 8),
                                      Coverage = c(1, 4 / 9, 8 / 15))
    txt <- .lamb_low_coverage_text(r)
    expect_match(txt, "Squamous 4/9")
    expect_match(txt, "MA-Lo 8/15")
    expect_false(grepl("LHS", txt))
    expect_match(.lamb_notes("LHS", "low_coverage", txt), "Few genes measured in this data for Squamous 4/9")
    attr(r, "coverage")$Coverage <- 1
    expect_identical(.lamb_low_coverage_text(r), "")
    f <- tempfile(fileext = ".pdf")
    expect_silent(lamb_report(r, f))
    expect_true(file.exists(f))
})

test_that("the drawn steps show only reached steps, at most two scores each", {
    res <- sim_results(4)
    for (i in seq_len(nrow(res))) {
        m <- .lamb_rpt_steps(as.list(as.data.frame(res)[i, ]))
        expect_true(all(vapply(m$steps, function(s) length(s$rows) <= 2, TRUE)))
        decides <- vapply(m$steps, `[[`, TRUE, "decides")
        if (res$Pred_Class[i] != "NC") expect_equal(sum(decides), 1)
    }
    malo <- list(Sample = "s", Pred_Class = "MA-Lo", Margin = 0.2, Score_Chondroid = -0.1,
                 Score_Squamous = -0.1, Score_Basal_Like = 0, Score_LHS = 0.05, Score_Apocrine = 0.1,
                 Score_MALO = 0.3)
    m <- .lamb_rpt_steps(malo)
    last <- m$steps[[length(m$steps)]]
    expect_true(last$decides)
    expect_match(last$a, "Yes: MA-Lo", fixed = TRUE)
    expect_equal(vapply(last$rows, `[[`, "", "text"), c("MA-Lo", "Apocrine"))
    # step 3 shows the lean, never a plain "LHS" or "Apocrine" that contradicts the MA-Lo class
    s3 <- m$steps[[length(m$steps) - 1]]
    expect_identical(s3$a, "Apocrine-leaning")
    expect_false(any(vapply(s3$rows, `[[`, TRUE, "win")))
    malo$Score_LHS <- 0.095
    expect_identical(m2 <- .lamb_rpt_steps(malo)$steps[[length(m$steps) - 1]]$a, "Balanced between LHS and Apocrine")
    nc <- list(Sample = "s", Pred_Class = "NC", Margin = NA, Score_Squamous = NA, Score_Chondroid = 0,
               Score_Basal_Like = 0, Score_LHS = 0, Score_Apocrine = 0, Unscored_Signatures = "squamous")
    m <- .lamb_rpt_steps(nc)
    expect_length(m$steps, 1)
    expect_match(m$steps[[1]]$a, "Squamous not scored")
    expect_match(m$rest, "not reached")
})

test_that("chondroid coverage below 100% is warned about on targeted panels only", {
    r <- data.frame(Sample = "s1", Pred_Class = "LHS", Margin = 0.3, Flag = "")
    attr(r, "coverage") <- data.frame(Signature = c("chondroid", "luminal"), N_Requested = c(9, 15),
                                      N_Found = c(8, 15), Coverage = c(8 / 9, 1))
    attr(r, "universe_coverage") <- 12000 / length(lamb_universe())   # whole transcriptome
    expect_identical(.lamb_s2_chondroid_text(r), "")
    attr(r, "universe_coverage") <- 400 / length(lamb_universe())     # targeted panel
    expect_identical(.lamb_s2_chondroid_text(r), "Chondroid: partly measured on this panel.")
    # 8/9 is above the 60% minimum, so it is not a low-coverage signature
    expect_identical(.lamb_low_coverage_text(r), "")
    attr(r, "coverage")$N_Found[1] <- 9
    expect_identical(.lamb_s2_chondroid_text(r), "")
    attr(r, "coverage")$N_Found[1] <- 0
    expect_identical(.lamb_s2_chondroid_text(r), "Chondroid: not measured in this data, so it cannot be called.")
})

test_that("the HER2 amplicon read-out renders on bulk and cluster pages", {
    res <- sim_results(4)
    expect_true(all(res$HER2_Amplicon %in% c("high", "not high", "not assessable")))
    r <- as.list(as.data.frame(res)[1, ])
    r$HER2_Amplicon <- "high"; r$HER2_Amplicon_N_High <- 4; r$HER2_Amplicon_N_Measured <- 5
    h <- .lamb_s2_her2_lines(r, res)
    expect_identical(h$word, "High")
    expect_match(h$txt, "4 of 5 measured amplicon genes rank in this sample's top 9%", fixed = TRUE)
    r$HER2_Amplicon <- "not assessable"
    expect_identical(.lamb_s2_her2_header(r, res)$txt, "HER2 amplicon: not assessed on this platform.")
    expect_identical(.lamb_s2_mini_text(list(HER2_Amplicon = "not assessable")),
                     "HER2 amplicon: not assessed on this platform.")
    expect_identical(.lamb_s2_mini_text(list(HER2_Amplicon = "not high")), "")
    expect_identical(.lamb_s2_mini_text(list(HER2_Amplicon = "high")), "HER2 amplicon high")
    hres <- data.frame(Sample = "h", Pred_Class = "LHS", Margin = 0.3, HER2_Amplicon = "high",
                       HER2_Amplicon_N_High = 4, HER2_Amplicon_N_Measured = 5)
    f <- tempfile(fileext = ".pdf")
    lamb_report(hres, f)
    expect_equal(pdf_pages(f), 1)

})

test_that("MA-Lo lean follows the higher of the LHS and Apocrine scores, for MA-Lo only", {
    r <- data.frame(Pred_Class = c("MA-Lo", "MA-Lo", "LHS", "MA-Lo"), Score_LHS = c(0.2, -0.1, 0.3, NA),
                    Score_Apocrine = c(0.1, 0.05, 0.1, 0.2))
    expect_identical(.lamb_malo_lean(r), c("LHS-leaning", "apocrine-leaning", "", ""))
})

test_that("MA-Lo lean tie follows the specification margin", {
    r <- data.frame(Pred_Class = "MA-Lo", Score_LHS = 0.14, Score_Apocrine = 0.10)
    expect_identical(.lamb_malo_lean(r), "LHS-leaning")
    expect_identical(.lamb_malo_lean(r, tie = 0.05), "balanced")
    attr(r, "specification") <- list(min_margin = 0.05)
    expect_identical(.lamb_malo_lean(r), "balanced")
})

test_that("MA-Lo lean is balanced when LHS and Apocrine differ by less than 0.02", {
    r <- data.frame(Pred_Class = "MA-Lo", Score_LHS = c(0.1123, 0.1123, 0.1377),
                    Score_Apocrine = c(0.1144, 0.1100, 0.0919))
    expect_identical(.lamb_malo_lean(r), c("balanced", "balanced", "LHS-leaning"))
    # the lean is the step-3 answer in a cluster panel, so the extra line carries HER2 only
    expect_identical(.lamb_s2_mini_text(cbind(r[1, ], HER2_Amplicon = "high")), "HER2 amplicon high")
})

test_that("step rows name the compared group, never a class that contradicts the call", {
    malo <- list(Sample = "s", Pred_Class = "MA-Lo", Margin = 0.2, Score_Chondroid = -0.1,
                 Score_Squamous = -0.1, Score_Basal_Like = 0, Score_LHS = 0.05, Score_Apocrine = 0.1,
                 Score_MALO = 0.3)
    m <- .lamb_rpt_steps(malo)
    expect_equal(vapply(m$steps[[1]]$rows, `[[`, "", "text"), c("Strongest metaplastic", "Strongest mammary"))
    expect_equal(vapply(m$steps[[2]]$rows, `[[`, "", "text"), c("Basal-like", "Strongest luminal"))
})

test_that("close-call wording names the pair without implying a probability", {
    r <- list(Sample = "s", Pred_Class = "MA-Lo", Margin = 0.01, Score_Chondroid = -0.1,
              Score_Squamous = -0.1, Score_Basal_Like = 0, Score_LHS = 0.2, Score_Apocrine = 0.1,
              Score_MALO = 0.21)
    expect_identical(.lamb_close_beat(r), "MA-Lo beat LHS")
    expect_identical(.lamb_close_beat(list(Pred_Class = "NC")), "The winning score beat the runner-up")
    expect_identical(.LAMB_LEGEND, c("wins the step", "does not win", "Further right = higher score."))
})

test_that("marker ranks say what they are and print below 1 as <1", {
    expect_identical(.lamb_gp_num(c(0, 0.4, 0.6, 47.5, 100)), c("<1", "<1", "1", "48", "100"))
    expect_identical(.lamb_gp_subtitle(10104),
                     "Percentile among the 10,104 genes measured in this sample (100 = highest). RNA, not IHC.")
    expect_identical(.lamb_gp_subtitle(NA),
                     "Percentile among the genes measured in this sample (100 = highest). RNA, not IHC.")
})

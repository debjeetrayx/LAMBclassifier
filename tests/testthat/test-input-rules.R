# Rules documented in vignettes/input-rules.Rmd; each test names its rule.

rules_genes <- function() unique(c(unlist(default_gene_sets()), utils::head(lamb_universe(), 600)))

rules_log2 <- function(n = 4, seed = 7) {
    set.seed(seed)
    g <- rules_genes()
    matrix(stats::rnorm(length(g) * n, 7, 2), length(g), n, dimnames = list(g, paste0("S", seq_len(n))))
}

quiet_run <- function(f, ...) {
    log <- character()
    r <- withCallingHandlers(lamb_run(f, out_dir = tempfile("lamb_"), run_pam50 = FALSE, ...),
                             message = function(m) { log <<- c(log, conditionMessage(m)); invokeRestart("muffleMessage") })
    r$log <- log
    r
}

test_that("rule: data centred per gene is refused", {
    x <- rules_log2(5)
    f <- tempfile(fileext = ".tsv")
    utils::write.table(x - rowMeans(x), f, sep = "\t", quote = FALSE, col.names = NA)
    expect_error(quiet_run(f), "centred per gene")
    utils::write.table(x - apply(x, 1, stats::median), f, sep = "\t", quote = FALSE, col.names = NA)
    expect_error(quiet_run(f), "centred per gene")
})

test_that("rule: z-scoring each sample is accepted (it keeps the order of genes)", {
    x <- rules_log2(4)
    f <- tempfile(fileext = ".tsv")
    utils::write.table(x, f, sep = "\t", quote = FALSE, col.names = NA)
    a <- quiet_run(f)
    utils::write.table(scale(x), f, sep = "\t", quote = FALSE, col.names = NA)
    b <- quiet_run(f)
    expect_identical(a$calls$LAMB_class, b$calls$LAMB_class)
})

test_that("rule: thousands separators are read as one number", {
    expect_equal(.lamb_num(c("1,234.5", "12", "1,234,567")), c(1234.5, 12, 1234567))
    expect_equal(.lamb_num(c("1.234,5", "7,25"), dec = ","), c(1234.5, 7.25))
    expect_equal(.lamb_num("1 234"), 1234)
    x <- round(2^rules_log2(3) * 10)
    f1 <- tempfile(fileext = ".csv"); f2 <- tempfile(fileext = ".csv")
    utils::write.csv(x, f1)
    utils::write.csv(format(x, big.mark = ","), f2)
    expect_identical(quiet_run(f1)$calls$LAMB_class, quiet_run(f2)$calls$LAMB_class)
})

test_that("rule: summary rows are dropped and negative log data is never logged again", {
    x <- rules_log2(3)
    x[1:50, ] <- -3                       # log values below -1
    y <- rbind(x, Total = colSums(abs(x)))  # a stray row pushes the maximum above 100
    f <- tempfile(fileext = ".tsv")
    utils::write.table(y, f, sep = "\t", quote = FALSE, col.names = NA)
    r <- expect_silent(quiet_run(f))
    expect_true(any(grepl("Dropped summary row", r$log)))
    expect_false(any(grepl("applied log2", r$log)))
    f2 <- tempfile(fileext = ".tsv")
    utils::write.table(x, f2, sep = "\t", quote = FALSE, col.names = NA)
    expect_identical(r$calls$LAMB_class, quiet_run(f2)$calls$LAMB_class)
})

test_that("rule: values = declares what the numbers are, skipping detection", {
    x <- round(2^rules_log2(3) * 10)       # whole numbers that are not counts
    f <- tempfile(fileext = ".tsv")
    utils::write.table(x, f, sep = "\t", quote = FALSE, col.names = NA)
    auto <- quiet_run(f)
    expect_true(any(grepl("raw counts", auto$log)))
    tpm <- quiet_run(f, values = "tpm")
    expect_true(any(grepl("declared as tpm", tpm$log)))
    expect_false(any(grepl("raw counts", tpm$log)))
    expect_error(lamb_run(f, values = "percent"), "should be one of")
})

test_that("rule: a GEO series matrix and non-upper-case names are refused with a reason", {
    f <- tempfile(fileext = ".txt")
    writeLines(c("!Series_title\t\"x\"", "!series_matrix_table_begin", "\"ID_REF\"\t\"GSM1\"", "\"1007_s_at\"\t7.1"), f)
    expect_error(quiet_run(f), "GEO series matrix")
    x <- rules_log2(2)
    rownames(x) <- paste0(substr(rownames(x), 1, 1), tolower(substring(rownames(x), 2)))
    f2 <- tempfile(fileext = ".tsv")
    utils::write.table(x, f2, sep = "\t", quote = FALSE, col.names = NA)
    expect_error(quiet_run(f2), "not in upper case")
})

test_that("rule: missing values are treated like unmeasured genes, per sample", {
    x <- rules_log2(3)
    lhs <- intersect(DEFAULT_LHS_GENES, rownames(x))
    clean <- suppressWarnings(suppressMessages(LAMBclassifier(x, run_pam50 = FALSE)))
    x[lhs[1], 1] <- NA                                     # one LHS gene missing in sample 1
    x[lhs, 3] <- NA                                        # every LHS gene missing in sample 3
    x[setdiff(rownames(x), unlist(default_gene_sets()))[1], 2] <- NA   # not a signature gene
    r <- suppressWarnings(suppressMessages(LAMBclassifier(x, run_pam50 = FALSE)))
    expect_false(r$Pred_Class[1] == "NC")                  # scored on its measured genes
    expect_identical(r$Pred_Class[2], clean$Pred_Class[2])   # gap outside signatures: same class
    # ranked without one of ~670 genes: scores move by about one rank step (1/670)
    expect_equal(r$Score_Chondroid[2], clean$Score_Chondroid[2], tolerance = 5e-3)
    expect_identical(r$Pred_Class[3], "NC")                # a required signature has no genes left
    expect_equal(unname(attr(r, "missing_values")[1:2]), c(1, 1))
    expect_lt(attr(r, "sample_coverage")[1, "luminal"], 1)
    f <- tempfile(fileext = ".tsv")
    utils::write.table(x, f, sep = "\t", quote = FALSE, col.names = NA, na = "")
    out <- quiet_run(f)
    expect_match(out$calls$Note[1], "1 signature gene value\\(s\\) are missing in this sample; it was scored on its measured genes")
    expect_match(out$calls$Note[3], "Not classifiable: .* missing in this sample")
})

test_that("rule: samples without missing values score exactly as before", {
    x <- rules_log2(4)
    a <- suppressWarnings(suppressMessages(LAMBclassifier(x, run_pam50 = FALSE)))
    y <- x; y[1, 4] <- NA
    b <- suppressWarnings(suppressMessages(LAMBclassifier(y, run_pam50 = FALSE)))
    sc <- grep("^Score_|^Margin$|^Pred_Class$|^Flag$", names(a), value = TRUE)
    expect_identical(as.data.frame(a)[1:3, sc], as.data.frame(b)[1:3, sc])
})

test_that("rule: duplicated gene names keep the highest-mean row (non-counts) or are summed (counts)", {
    x <- rules_log2(2)
    g <- rownames(x)[1]
    d <- rbind(x, x[g, , drop = FALSE] - 5)
    rownames(d)[nrow(d)] <- g
    out <- .lamb_symbols(d, function(...) NULL, counts = FALSE)
    expect_equal(unname(out[g, ]), unname(x[g, ]))
    cnt <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("A", "A"), c("s1", "s2")))
    expect_equal(unname(.lamb_symbols(cnt, function(...) NULL, counts = TRUE)["A", ]), c(3, 7))
})

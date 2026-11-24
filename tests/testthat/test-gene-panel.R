gp_pdf_pages <- function(path) {
    raw <- readBin(path, "raw", file.info(path)$size)
    txt <- rawToChar(raw[raw != as.raw(0)])
    m <- regmatches(txt, regexpr("/Type /Pages[^>]*/Count [0-9]+", txt, useBytes = TRUE))
    as.integer(sub(".*/Count ", "", m))
}

gp_sim_results <- function(n = 4, seed = 11) {
    set.seed(seed)
    genes <- unique(c(DEFAULT_ALL_INPUT_GENES, sample(lamb_universe(), 2000)))
    x <- matrix(stats::rnorm(length(genes) * n, 6, 2), nrow = length(genes),
                dimnames = list(genes, paste0("S", seq_len(n))))
    sigs <- list(DEFAULT_LHS_GENES, DEFAULT_BASAL_LIKE_GENES, DEFAULT_MALO_GENES)
    for (i in seq_len(n)) {
        g <- intersect(sigs[[(i - 1) %% length(sigs) + 1]], genes)
        x[g, i] <- x[g, i] + 6
    }
    suppressWarnings(LAMBclassifier(x, run_pam50 = FALSE))
}

test_that("percentile helper: 100 = highest, ties averaged, NA and non-universe genes excluded", {
    x <- matrix(c(5, 1, 3, 3, NA, 9,
                  1, 2, 3, 4, 5, 6),
                ncol = 2, dimnames = list(c("A", "B", "C", "D", "E", "Z"), c("s1", "s2")))
    uni <- c("A", "B", "C", "D", "E")   # Z is outside the universe
    gp <- .lamb_gp_percentiles(x, c("A", "B", "C", "D", "E", "Z", "Q"), universe = uni)
    # s1: measured A=5, B=1, C=3, D=3 (n = 4); ranks 4, 1, 2.5, 2.5 -> 100 * (rank - 1) / 3
    expect_equal(unname(gp$pct["s1", c("A", "B", "C", "D")]), c(100, 0, 50, 50))
    expect_true(is.na(gp$pct["s1", "E"]))
    expect_true(all(is.na(gp$pct[, "Z"])))
    expect_true(all(is.na(gp$pct[, "Q"])))
    # s2: A..E = 1..5 (n = 5), Z ignored although it is the highest value
    expect_equal(unname(gp$pct["s2", c("A", "B", "C", "D", "E")]), c(0, 25, 50, 75, 100))
    expect_equal(gp$n, c(s1 = 4L, s2 = 5L))
})

test_that("LAMBclassifier attaches gene percentiles and the ranked-gene count", {
    res <- gp_sim_results(3)
    pct <- attr(res, "gene_percentiles")
    expect_true(is.matrix(pct))
    expect_identical(rownames(pct), res$Sample)
    want <- unique(c(.LAMB_GP_MARKERS, unlist(default_gene_sets()), DEFAULT_HER2_AMPLICON_GENES))
    expect_setequal(colnames(pct), want)
    expect_identical(colnames(pct)[1:9], .LAMB_GP_MARKERS)
    v <- pct[!is.na(pct)]
    expect_true(all(v >= 0 & v <= 100))
    expect_equal(v, round(v, 1))
    n <- attr(res, "n_genes_ranked")
    expect_identical(names(n), res$Sample)
    expect_true(all(n > 1000))
})

test_that("lamb_run writes the percentile CSV and LAMB_genes.pdf", {
    out <- tempfile()
    f <- system.file("extdata", "example_bulk.csv", package = "LAMBclassifier")
    skip_if(!nzchar(f))
    r <- suppressWarnings(suppressMessages(lamb_run(f, out_dir = out, run_pam50 = FALSE)))
    csv <- utils::read.csv(file.path(out, "technical", "LAMB_gene_percentiles.csv"), check.names = FALSE)
    expect_identical(names(csv)[1], "Sample")
    expect_identical(csv$Sample, r$details$Sample)
    expect_identical(names(csv)[2:10], .LAMB_GP_MARKERS)
    expect_equal(gp_pdf_pages(file.path(out, "LAMB_genes.pdf")), nrow(r$calls))
    expect_match(paste(readLines(file.path(out, "technical", "run_log.txt")), collapse = " "), "LAMB_genes.pdf")
})

test_that("lamb_gene_report: one page per sample, subset, and genes not measured", {
    res <- gp_sim_results(4)
    f <- tempfile(fileext = ".pdf")
    expect_invisible(lamb_gene_report(res, f))
    expect_equal(gp_pdf_pages(f), 4)
    f2 <- tempfile(fileext = ".pdf")
    lamb_gene_report(res, f2, samples = res$Sample[2:3])
    expect_equal(gp_pdf_pages(f2), 2)

    # a whole group and some markers not measured: the page still renders
    pct <- attr(res, "gene_percentiles")
    pct[, intersect(colnames(pct), c(DEFAULT_CHONDROID_GENES, "ESR1", "PGR"))] <- NA
    attr(res, "gene_percentiles") <- pct
    f3 <- tempfile(fileext = ".pdf")
    lamb_gene_report(res, f3)
    expect_equal(gp_pdf_pages(f3), 4)
    f4 <- tempfile(fileext = ".pdf")
    lamb_report(res, f4)
    expect_equal(gp_pdf_pages(f4), 5)

    expect_error(lamb_gene_report(as.data.frame(res)[, c("Sample", "Pred_Class")], tempfile()),
                 "gene_percentiles")
})

test_that("marker strip lists measured markers only and is omitted when none is", {
    res <- gp_sim_results(1)
    s <- res$Sample[1]
    pct <- attr(res, "gene_percentiles")
    pct[, .LAMB_GP_MARKERS] <- NA
    attr(res, "gene_percentiles") <- pct
    f <- tempfile(fileext = ".pdf")
    grDevices::pdf(f, width = .LAMB_PAGE[["w"]], height = .LAMB_PAGE[["h"]])
    .lamb_new_page()
    expect_identical(.lamb_gp_strip(res, s, 0.6, 4, 5), 5)        # nothing drawn
    pct[, c("ESR1", "AR")] <- c(90, 10)
    attr(res, "gene_percentiles") <- pct
    y <- .lamb_gp_strip(res, s, 0.6, 4, 5)
    cap <- .lamb_text_height(.lamb_gp_subtitle(.lamb_gp_row(res, s)$n), 3.4, 0.72)   # wrapped caption
    grDevices::dev.off()
    expect_equal(y, 5 - 0.1 - cap - 0.12 - 2 * 0.2)              # two rows
})

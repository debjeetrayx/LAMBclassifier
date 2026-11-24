rb_genes <- function(n = 1200) {
    gl <- get_gene_length()
    keep <- !is.na(gl$SYMBOL) & nzchar(gl$SYMBOL) & !duplicated(gl$SYMBOL) & !duplicated(gl$gene_id)
    gl <- gl[keep, ]
    list(symbol = utils::head(gl$SYMBOL, n), ensembl = utils::head(gl$gene_id, n))
}

# read quietly, keeping the say() lines for inspection
rb_read <- function(path) {
    msgs <- character()
    r <- .lamb_read_bulk(path, say = function(...) msgs <<- c(msgs, paste0(...)))
    r$msgs <- msgs
    r
}

rb_dir <- function() {
    d <- tempfile("lamb_read_bulk_")
    dir.create(d)
    d
}

rb_log2 <- function(n_genes = 1200, n_samples = 4, seed = 1) {
    set.seed(seed)
    g <- rb_genes(n_genes)$symbol
    matrix(round(stats::rnorm(n_genes * n_samples, 8, 2), 6), n_genes,
           dimnames = list(g, paste0("GSM", 100 + seq_len(n_samples))))
}

write_salmon <- function(path, names, tpm, with_tpm = TRUE) {
    df <- data.frame(Name = names, Length = 1500L, EffectiveLength = 1350.5,
                     TPM = tpm, NumReads = round(tpm * 10, 3), check.names = FALSE)
    if (!with_tpm) df$TPM <- NULL
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_star <- function(path, genes, unstranded, fwd, rev) {
    fwd <- rep_len(fwd, length(genes)); rev <- rep_len(rev, length(genes))
    df <- data.frame(c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous", genes),
                     c(10, 20, 30, 40, unstranded), c(11, 21, 31, 41, fwd), c(12, 22, 32, 42, rev))
    dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
    utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}


# ---- plain tables ----

test_that("a GSE-shaped tab table (blank top-left cell) round-trips identically to read.table", {
    m <- rb_log2()
    f <- file.path(rb_dir(), "matrix.txt")
    lines <- c(paste(c("", colnames(m)), collapse = "\t"),
               vapply(seq_len(nrow(m)), function(i) paste(c(rownames(m)[i], format(m[i, ], digits = 15)), collapse = "\t"), ""))
    writeLines(lines, f)
    old <- as.matrix(utils::read.table(f, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE,
                                       quote = "\"", comment.char = ""))
    r <- rb_read(f)
    expect_identical(r$expr, old)
    expect_equal(r$units, "unknown")
    expect_equal(r$platform, "table")
    expect_equal(r$source_platform, "table")
})

test_that("R write.table style (header one field short) and gzip are read", {
    m <- rb_log2()
    d <- rb_dir()
    f <- file.path(d, "m.tsv.gz")
    con <- gzfile(f, "w")
    utils::write.table(m, con, sep = "\t", quote = FALSE)
    close(con)
    r <- rb_read(f)
    expect_equal(r$expr, m)
})

test_that("CSV with quoted headers and a named gene column", {
    m <- rb_log2()
    f <- file.path(rb_dir(), "m.csv")
    utils::write.csv(data.frame(Gene = rownames(m), m, check.names = FALSE), f, row.names = FALSE)
    expect_true(grepl('^"Gene","GSM101"', readLines(f, 1)))
    r <- rb_read(f)
    expect_equal(r$expr, m)
})

test_that("European CSV: semicolon separator and decimal comma", {
    m <- rb_log2()
    f <- file.path(rb_dir(), "euro.csv")
    utils::write.csv2(data.frame(Gene = rownames(m), m, check.names = FALSE), f, row.names = FALSE)
    r <- rb_read(f)
    expect_equal(r$expr, m)
    expect_true(any(grepl("decimal comma", r$msgs)))
})

test_that("integer count tables come back as a double matrix", {
    g <- rb_genes()$symbol
    m <- matrix(as.integer(stats::rpois(length(g) * 3, 50)), length(g), dimnames = list(g, c("A", "B", "C")))
    f <- file.path(rb_dir(), "counts.tsv")
    utils::write.table(data.frame(gene = g, m), f, sep = "\t", quote = FALSE, row.names = FALSE)
    r <- rb_read(f)
    expect_type(r$expr, "double")
    expect_equal(r$expr, m * 1)
})


# ---- annotation columns ----

test_that("extra annotation columns: symbol column chosen, others named and dropped", {
    m <- rb_log2()
    g <- rb_genes()
    df <- data.frame(gene_id = g$ensembl, symbol = g$symbol, description = "some protein",
                     entrez = seq_along(g$symbol), m, check.names = FALSE)
    f <- file.path(rb_dir(), "annot.tsv")
    utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
    r <- rb_read(f)
    expect_equal(dim(r$expr), dim(m))
    expect_identical(rownames(r$expr), g$symbol)
    expect_identical(colnames(r$expr), colnames(m))
    expect_true(any(grepl("'symbol'.*gene_id.*description", r$msgs)))
    expect_true(any(grepl("annotation column.*entrez", r$msgs)))
})

test_that("Ensembl-only gene column is accepted when there are no symbols", {
    m <- rb_log2()
    g <- rb_genes()
    df <- data.frame(description = "x", id = paste0(g$ensembl, ".7"), m, check.names = FALSE)
    f <- file.path(rb_dir(), "ens.csv")
    utils::write.csv(df, f, row.names = FALSE)
    r <- rb_read(f)
    expect_identical(rownames(r$expr), paste0(g$ensembl, ".7"))
})

test_that("probe IDs (no gene-like column) are refused in plain language", {
    m <- rb_log2()
    rownames(m) <- paste0(seq_len(nrow(m)), "_at")
    f <- file.path(rb_dir(), "probes.tsv")
    utils::write.table(data.frame(ID = rownames(m), m), f, sep = "\t", quote = FALSE, row.names = FALSE)
    expect_error(rb_read(f), "Could not find a column of gene names.*probe IDs")
})


# ---- transposed ----

test_that("samples in rows are transposed", {
    m <- rb_log2()
    df <- data.frame(Sample = colnames(m), ER_status = "pos", t(m), check.names = FALSE)
    f <- file.path(rb_dir(), "transposed.csv")
    utils::write.csv(df, f, row.names = FALSE)
    r <- rb_read(f)
    expect_equal(r$expr, m)
    expect_true(any(grepl("transposed", r$msgs)))
    expect_true(any(grepl("ER_status", r$msgs)))
})


# ---- Salmon / kallisto / RSEM ----

test_that("Salmon quant.sf gives ONE sample (TPM), named after its folder", {
    g <- rb_genes()$ensembl
    tpm <- stats::runif(length(g), 0, 200)
    f <- file.path(rb_dir(), "patient07", "quant.sf")
    write_salmon(f, g, tpm)
    r <- rb_read(f)
    expect_equal(ncol(r$expr), 1)
    expect_identical(colnames(r$expr), "patient07")
    expect_equal(unname(r$expr[, 1]), tpm)
    expect_equal(r$units, "tpm")
    expect_equal(r$platform, "salmon")
})

test_that("Salmon without a TPM column is refused", {
    g <- rb_genes()$symbol
    f <- file.path(rb_dir(), "s1", "quant.sf")
    write_salmon(f, g, stats::runif(length(g)), with_tpm = FALSE)
    expect_error(rb_read(f), "This looks like a Salmon file but has no TPM column.", fixed = TRUE)
})

test_that("transcript-level Salmon without gene field is refused, not guessed", {
    tx <- sprintf("ENST%011d.1", seq_len(1500))
    f <- file.path(rb_dir(), "s1", "quant.sf")
    write_salmon(f, tx, stats::runif(1500))
    expect_error(rb_read(f), "one value per transcript.*tx2gene")
})

test_that("GENCODE pipe-format transcript names are summed to genes", {
    g <- rb_genes()
    gi <- rep(seq_along(g$symbol), each = 2)
    tx <- sprintf("ENST%011d.1|%s.5|OTTHUMG1|OTTHUMT1|%s-20%d|%s|1500|protein_coding|",
                  seq_along(gi), g$ensembl[gi], g$symbol[gi], rep(1:2, length(g$symbol)), g$symbol[gi])
    tpm <- rep(c(1, 2), length(g$symbol))
    f <- file.path(rb_dir(), "s9", "quant.sf")
    write_salmon(f, tx, tpm)
    r <- rb_read(f)
    expect_equal(nrow(r$expr), length(g$symbol))
    expect_identical(rownames(r$expr), g$symbol)
    expect_true(all(r$expr[, 1] == 3))
    expect_true(any(grepl("summed", r$msgs)))
})

test_that("kallisto abundance.tsv gives one sample of TPM", {
    g <- rb_genes()$symbol
    tpm <- stats::runif(length(g), 0, 100)
    d <- file.path(rb_dir(), "pt3")
    dir.create(d)
    utils::write.table(data.frame(target_id = g, length = 1000, eff_length = 850.2, est_counts = tpm * 3, tpm = tpm),
                       file.path(d, "abundance.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    r <- rb_read(file.path(d, "abundance.tsv"))
    expect_identical(colnames(r$expr), "pt3")
    expect_equal(unname(r$expr[, 1]), tpm)
    expect_equal(c(r$units, r$platform), c("tpm", "kallisto"))
})

test_that("RSEM genes.results uses TPM", {
    g <- rb_genes()$ensembl
    tpm <- stats::runif(length(g), 0, 100)
    f <- file.path(rb_dir(), "TUM12.genes.results")
    utils::write.table(data.frame(gene_id = g, `transcript_id(s)` = "ENST1,ENST2", length = 1000, effective_length = 900,
                                  expected_count = tpm * 5, TPM = tpm, FPKM = tpm * 2, check.names = FALSE),
                       f, sep = "\t", quote = FALSE, row.names = FALSE)
    r <- rb_read(f)
    expect_identical(colnames(r$expr), "TUM12")
    expect_equal(unname(r$expr[, 1]), tpm)
    expect_equal(c(r$units, r$platform), c("tpm", "rsem"))
})


# ---- featureCounts / STAR / HTSeq ----

test_that("featureCounts output: comment line, annotation columns ignored, BAM paths cleaned", {
    g <- rb_genes()$ensembl
    cnt <- matrix(stats::rpois(length(g) * 2, 30), ncol = 2)
    f <- file.path(rb_dir(), "counts.txt")
    writeLines('# Program:featureCounts v2.0.6; Command:"featureCounts" "-a" "gtf"', f)
    df <- data.frame(Geneid = g, Chr = "chr1;chr1", Start = "100;200", End = "150;250", Strand = "+;+",
                     Length = 1234L, cnt, check.names = FALSE)
    names(df)[7:8] <- c("/data/bam/P1Aligned.sortedByCoord.out.bam", "/data/bam/P2.bam")
    suppressWarnings(utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE))
    r <- rb_read(f)
    expect_identical(colnames(r$expr), c("P1", "P2"))
    expect_equal(unname(r$expr), cnt * 1)
    expect_equal(c(r$units, r$platform), c("counts", "featurecounts"))
})

test_that("featureCounts .summary file is refused", {
    f <- file.path(rb_dir(), "counts.txt.summary")
    writeLines(c("Status\tP1.bam", "Assigned\t100", "Unassigned_NoFeatures\t5"), f)
    expect_error(rb_read(f), "featureCounts summary")
})

test_that("STAR ReadsPerGene: summary rows dropped, column with most counts chosen and named", {
    g <- rb_genes()$ensembl
    rev <- stats::rpois(length(g), 100)
    f <- file.path(rb_dir(), "P5ReadsPerGene.out.tab")
    write_star(f, g, unstranded = rev - 1, fwd = rpois(length(g), 2), rev = rev)
    r <- rb_read(f)
    expect_equal(dim(r$expr), c(length(g), 1))
    expect_identical(colnames(r$expr), "P5")
    expect_equal(unname(r$expr[, 1]), rev * 1)
    expect_false(any(grepl("^N_", rownames(r$expr))))
    expect_true(any(grepl("column 4 \\(read 2.*reverse", r$msgs)))
    expect_equal(c(r$units, r$platform), c("counts", "star"))
})

test_that("HTSeq-count (no header, __ summary rows) gives one sample", {
    g <- rb_genes()$symbol
    f <- file.path(rb_dir(), "P8.counts")
    utils::write.table(data.frame(c(g, "__no_feature", "__ambiguous"), c(seq_along(g), 5, 6)), f,
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    r <- rb_read(f)
    expect_identical(colnames(r$expr), "P8")
    expect_equal(nrow(r$expr), length(g))
    expect_equal(r$platform, "htseq")
})


# ---- folders ----

test_that("folder of Salmon sample folders combines by gene (intersection)", {
    g <- rb_genes()$symbol
    d <- rb_dir()
    write_salmon(file.path(d, "P1", "quant.sf"), g, seq_along(g))
    write_salmon(file.path(d, "P2", "quant.sf"), g[-1], seq_along(g[-1]))
    write_salmon(file.path(d, "P3", "quant.sf"), rev(g), seq_along(g))
    r <- rb_read(d)
    expect_identical(colnames(r$expr), c("P1", "P2", "P3"))
    expect_equal(nrow(r$expr), length(g) - 1)
    expect_equal(unname(r$expr[g[5], ]), c(5, 4, length(g) - 4))
    expect_equal(c(r$units, r$platform, r$source_platform), c("tpm", "multi_file", "salmon"))
    expect_true(any(grepl("kept 1,199 genes", r$msgs)))
})

test_that("folder of STAR files at top level is combined; junction files ignored", {
    g <- rb_genes()$ensembl
    d <- rb_dir()
    for (p in c("A", "B")) write_star(file.path(d, paste0(p, "_ReadsPerGene.out.tab")), g, rpois(length(g), 50), 1, 1)
    writeLines("chr1\t1\t2", file.path(d, "A_SJ.out.tab"))
    r <- rb_read(d)
    expect_identical(colnames(r$expr), c("A", "B"))
    expect_equal(r$source_platform, "star")
})

test_that("folder mixing formats is refused", {
    g <- rb_genes()$symbol
    d <- rb_dir()
    write_salmon(file.path(d, "P1", "quant.sf"), g, seq_along(g))
    write_star(file.path(d, "P2ReadsPerGene.out.tab"), g, rpois(length(g), 5), 1, 1)
    expect_error(rb_read(d), "mixes different kinds of files")
})

test_that("empty folder is refused", {
    expect_error(rb_read(rb_dir()), "No expression files")
})


# ---- Excel and RDS ----

test_that(".xlsx first sheet is read like a table", {
    skip_if_not_installed("readxl")
    skip_if_not_installed("writexl")
    m <- rb_log2()
    f <- file.path(rb_dir(), "cohort.xlsx")
    writexl::write_xlsx(list(expr = data.frame(Gene = rownames(m), m, check.names = FALSE),
                             notes = data.frame(a = 1)), f)
    r <- rb_read(f)
    expect_equal(r$expr, m)
    expect_equal(r$platform, "xlsx")
    expect_true(any(grepl("first sheet \\('expr'\\)", r$msgs)))
})

test_that(".rds matrix, data.frame and transposed matrix", {
    m <- rb_log2()
    d <- rb_dir()
    saveRDS(m, file.path(d, "m.rds"))
    saveRDS(as.data.frame(m), file.path(d, "df.rds"))
    saveRDS(t(m), file.path(d, "t.rds"))
    saveRDS(data.frame(gene = rownames(m), m, check.names = FALSE), file.path(d, "df2.rds"))
    for (f in c("m.rds", "df.rds", "t.rds", "df2.rds")) {
        r <- rb_read(file.path(d, f))
        expect_equal(r$expr, m, info = f)
        expect_equal(r$platform, "rds", info = f)
    }
})


# ---- value checks ----

test_that("missing values are kept and counted", {
    m <- rb_log2()
    m[1:3, 1] <- NA
    f <- file.path(rb_dir(), "na.tsv")
    utils::write.table(m, f, sep = "\t", quote = FALSE, na = "")
    r <- rb_read(f)
    expect_equal(sum(is.na(r$expr)), 3)
    expect_true(any(grepl("^3 value\\(s\\) are missing", r$msgs)))
})

test_that("all-negative values are refused", {
    m <- -abs(rb_log2())
    f <- file.path(rb_dir(), "neg.tsv")
    utils::write.table(m, f, sep = "\t", quote = FALSE)
    expect_error(rb_read(f), "All values are negative")
})

test_that("fewer than 1,000 genes warns but does not stop", {
    m <- rb_log2(n_genes = 200)
    f <- file.path(rb_dir(), "small.tsv")
    utils::write.table(m, f, sep = "\t", quote = FALSE)
    r <- rb_read(f)
    expect_equal(nrow(r$expr), 200)
    expect_true(any(grepl("only 200 genes.*whole-transcriptome", r$msgs)))
})

test_that("missing file and raw sequencing files get plain errors", {
    expect_error(.lamb_read_bulk(file.path(tempdir(), "nope.tsv")), "not found")
    f <- file.path(rb_dir(), "x.bam")
    writeLines("x", f)
    expect_error(rb_read(f), "raw sequencing")
})

test_that("say defaults to message()", {
    m <- rb_log2(n_genes = 50)
    f <- file.path(rb_dir(), "tiny.tsv")
    utils::write.table(m, f, sep = "\t", quote = FALSE)
    expect_message(.lamb_read_bulk(f), "only 50 genes")
})

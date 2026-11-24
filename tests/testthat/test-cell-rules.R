# Rules documented in vignettes/input-rules.Rmd.

cr_counts <- function(n = 120, seed = 3) {
    set.seed(seed)
    g <- unique(c(unlist(default_gene_sets()), utils::head(lamb_universe(), 400), "EPCAM", "CDH1", "KRT8", "KRT18"))
    m <- matrix(stats::rpois(length(g) * n, 2), length(g), n, dimnames = list(g, sprintf("AAAC%04d-1", seq_len(n))))
    Matrix::Matrix(m, sparse = TRUE)
}

test_that("rule: Matrix Market files under other names are read like a 10x folder, with a metadata.csv column as groups", {
    m <- cr_counts()
    d <- tempfile("atlas_"); dir.create(d)
    Matrix::writeMM(m, file.path(d, "count_matrix_sparse.mtx"))
    writeLines(rownames(m), file.path(d, "count_matrix_genes.tsv"))          # one column: symbols
    writeLines(colnames(m), file.path(d, "count_matrix_barcodes.tsv"))
    meta <- data.frame(barcode = colnames(m), celltype_major = rep(c("Cancer Epithelial", "Luminal"), each = 60),
                       nCount_RNA = seq_len(ncol(m)))
    utils::write.csv(meta, file.path(d, "metadata.csv"), row.names = FALSE)
    expect_identical(.lamb_detect_cells(d), "10x_mtx")
    x <- .lamb_read_cells(d, "10x_mtx", say = function(...) NULL)
    expect_equal(dim(x$counts), dim(m))
    expect_true("celltype_major" %in% names(x$meta))
    expect_false("nCount_RNA" %in% names(x$meta))            # an ID-like / continuous column is not offered
    expect_error(suppressMessages(lamb_run(d, out_dir = tempfile())), "celltype_major")   # no clusters: message names the column
    r <- suppressWarnings(suppressMessages(lamb_run(d, out_dir = tempfile(), clusters = "celltype_major")))
    expect_setequal(r$calls$Cluster, c("Cancer Epithelial", "Luminal"))
    expect_error(suppressMessages(lamb_run(d, out_dir = tempfile(), clusters = "celltype")), "neither a file nor a cell annotation column")
})

test_that("rule: a single-cell count matrix saved as a text table is refused, not scored cell by cell", {
    m <- as.matrix(cr_counts(600))
    m[m < 4] <- 0   # about 86% zeros, as in droplet data
    f <- tempfile(fileext = ".tsv")
    utils::write.table(m, f, sep = "\t", quote = FALSE, col.names = NA)
    expect_error(suppressMessages(lamb_run(f, out_dir = tempfile(), run_pam50 = FALSE)), "single-cell count matrix")
})

test_that("rule: sparse non-negative data is never mistaken for centred data", {
    x <- matrix(c(rep(0, 800), stats::rpois(400, 5)), 400, 3, dimnames = list(utils::head(lamb_universe(), 400), c("a", "b", "c")))
    expect_silent(.lamb_refuse_centred(x))
    expect_error(.lamb_refuse_centred(x - rowMeans(x)), "centred per gene")
})

test_that("rule: the cell-type gate leaves out immune clusters and keeps named tumour clusters", {
    set.seed(5)
    g <- unique(c(unlist(default_gene_sets()), utils::head(lamb_universe(), 400), unlist(.LAMB_GATE_MODULES)))
    m <- matrix(stats::rpois(length(g) * 3, 1), length(g), 3, dimnames = list(g, c("tumour", "macrophage", "chondroid")))
    m[.LAMB_GATE_MODULES$epithelial, "tumour"] <- 400
    m[c(.LAMB_GATE_MODULES$immune, .LAMB_GATE_MODULES$epithelial[1:2]), "macrophage"] <- c(rep(400, 12), 50, 50)  # some ambient epithelial RNA
    m[.LAMB_GATE_MODULES$fibroblast, "chondroid"] <- 300; m[c("COL2A1", "ACAN", "SP7"), "chondroid"] <- 500
    gate <- .lamb_gate(m, function(...) NULL)
    expect_identical(unname(gate), c("epithelial", "immune", "fibroblast-like"))
})

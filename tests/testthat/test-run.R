sim_genes <- function() {
    unique(c(unlist(default_gene_sets()), utils::head(lamb_universe(), 400),
             "EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "PTPRC", "CD3E", "CD2", "LYZ", "CD68"))
}

write_10x <- function(dir, m, xenium = FALSE, clusters = NULL) {
    mdir <- file.path(dir, if (xenium) "cell_feature_matrix" else "filtered_feature_bc_matrix")
    dir.create(mdir, recursive = TRUE)
    Matrix::writeMM(m, file.path(mdir, "matrix.mtx"))
    feat <- data.frame(id = rownames(m), name = rownames(m), type = "Gene Expression")
    if (xenium) feat <- rbind(feat, data.frame(id = "NegControlProbe_1", name = "NegControlProbe_1", type = "Negative Control Probe"))
    if (xenium) m <- rbind(m, Matrix::Matrix(0, 1, ncol(m), sparse = TRUE))
    if (xenium) Matrix::writeMM(m, file.path(mdir, "matrix.mtx"))
    utils::write.table(feat, file.path(mdir, "features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    writeLines(colnames(m), file.path(mdir, "barcodes.tsv"))
    if (!is.null(clusters)) {
        cd <- file.path(dir, "analysis", "clustering", "gene_expression_graphclust")
        dir.create(cd, recursive = TRUE)
        utils::write.csv(data.frame(Barcode = names(clusters), Cluster = clusters), file.path(cd, "clusters.csv"), row.names = FALSE)
    }
    if (xenium) {
        con <- gzfile(file.path(dir, "cells.csv.gz"), "w")
        utils::write.csv(data.frame(cell_id = colnames(m), x_centroid = stats::runif(ncol(m)) * 1000,
                                    y_centroid = stats::runif(ncol(m)) * 1000), con, row.names = FALSE)
        close(con)
    }
}

# two groups of cells: an epithelial luminal-like group and a non-epithelial group
sim_cells <- function(n = 120) {
    set.seed(7)
    g <- sim_genes()
    m <- matrix(stats::rpois(length(g) * n, 0.3), nrow = length(g), dimnames = list(g, paste0("c", seq_len(n), "-1")))
    epi <- seq_len(n / 2)
    m[intersect(c(DEFAULT_LHS_GENES, "EPCAM", "CDH1"), g), epi] <- m[intersect(c(DEFAULT_LHS_GENES, "EPCAM", "CDH1"), g), epi] + 5
    imm <- c("PTPRC", "CD3E", "CD2", "LYZ", "CD68"); m[imm, -epi] <- m[imm, -epi] + 5   # an immune module
    m[intersect(.LAMB_EPITHELIAL, g), -epi] <- 0   # no epithelial gene in the non-epithelial group
    list(m = Matrix::Matrix(m, sparse = TRUE),
         clusters = setNames(ifelse(seq_len(n) %in% epi, 1, 2), colnames(m)))
}

test_that("bulk: raw counts are detected, converted and give the same calls as LAMBclassifier on log2 TPM", {
    gl <- get_gene_length()
    g <- intersect(sim_genes(), gl$SYMBOL)
    set.seed(3)
    counts <- matrix(stats::rpois(length(g) * 3, 200), nrow = length(g), dimnames = list(g, paste0("S", 1:3)))
    f <- tempfile(fileext = ".csv")
    utils::write.csv(counts, f)
    out <- tempfile()
    r <- suppressWarnings(suppressMessages(lamb_run(f, out_dir = out, run_pam50 = FALSE)))
    ref <- suppressWarnings(LAMBclassifier(log2(counts_to_tpm(counts, gl) + 1), run_pam50 = FALSE))
    expect_equal(r$calls$LAMB_class, ref$Pred_Class)
    expect_true(all(file.exists(file.path(out, c("LAMB_calls.csv", "LAMB_report.pdf", "technical/run_log.txt")))))
    expect_match(paste(readLines(file.path(out, "technical/run_log.txt")), collapse = " "), "raw counts")
})

test_that("bulk: Ensembl IDs with versions are mapped to symbols", {
    gl <- get_gene_length()
    g <- intersect(sim_genes(), gl$SYMBOL)
    x <- matrix(stats::runif(length(g) * 2, 0, 10), nrow = length(g), dimnames = list(g, c("A", "B")))
    rownames(x) <- paste0(gl$gene_id[match(g, gl$SYMBOL)], ".7")
    y <- .lamb_symbols(x, function(...) NULL, counts = FALSE)
    expect_true(all(rownames(y) %in% g))
    expect_gt(nrow(y), 0.9 * length(g))
})

test_that("counts_to_tpm works with the bundled gene-length table", {
    gl <- get_gene_length()
    g <- utils::head(gl$SYMBOL[nzchar(gl$SYMBOL)], 200)
    counts <- matrix(10, nrow = length(g), ncol = 2, dimnames = list(g, c("a", "b")))
    tpm <- suppressWarnings(counts_to_tpm(counts, gl))
    expect_equal(unname(colSums(tpm)), c(1e6, 1e6))
})

test_that("single-cell folder: instrument clusters are used and non-epithelial clusters are marked", {
    s <- sim_cells()
    d <- tempfile(); write_10x(d, s$m, clusters = s$clusters)
    out <- tempfile()
    r <- suppressWarnings(suppressMessages(lamb_run(d, out_dir = out)))
    expect_equal(nrow(r$calls), 2)
    expect_equal(r$calls$Cell_type[r$calls$Cluster == "cluster 1"], "epithelial")
    expect_equal(r$calls$Cell_type[r$calls$Cluster == "cluster 2"], "immune")
    expect_false(file.exists(file.path(out, "LAMB_xenium_explorer.csv")))
})

test_that("Xenium folder: control features dropped, explorer file has every cell once", {
    s <- sim_cells()
    d <- tempfile(); write_10x(d, s$m, xenium = TRUE, clusters = s$clusters)
    out <- tempfile()
    r <- suppressWarnings(suppressMessages(lamb_run(d, out_dir = out)))
    ex <- utils::read.csv(file.path(out, "LAMB_xenium_explorer.csv"))
    expect_equal(nrow(ex), ncol(s$m))
    expect_false(anyDuplicated(ex$cell_id) > 0)
    expect_setequal(names(ex), c("cell_id", "group", "color"))
    expect_true(all(ex$group[ex$cell_id %in% names(s$clusters)[s$clusters == 2]] == "Not epithelial"))
    expect_match(paste(readLines(file.path(out, "technical/run_log.txt")), collapse = " "), "1 control / codeword")
})

test_that("small clusters are reported but not classified; own cluster file is honoured", {
    s <- sim_cells()
    own <- s$clusters; own[1:10] <- 3
    f <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(barcode = names(own), group = own), f, row.names = FALSE)
    d <- tempfile(); write_10x(d, s$m, clusters = s$clusters)
    r <- suppressWarnings(suppressMessages(lamb_run(d, out_dir = tempfile(), clusters = f, min_cells = 20)))
    expect_equal(r$calls$LAMB_class[r$calls$Cluster == "3"], "Too few cells")
})

test_that("Visium: the run log counts spots, not cells, and spots the instrument left unclustered are named so", {
    s <- sim_cells()
    d <- tempfile(); outs <- file.path(d, "outs")
    write_10x(outs, s$m, clusters = s$clusters[-1])   # Space Ranger clusters skip some spots
    dir.create(file.path(outs, "spatial"))
    n <- ncol(s$m)
    utils::write.csv(data.frame(barcode = colnames(s$m), in_tissue = 1L, array_row = seq_len(n), array_col = seq_len(n),
                                pxl_row_in_fullres = seq_len(n) * 100, pxl_col_in_fullres = seq_len(n) * 200),
                     file.path(outs, "spatial", "tissue_positions.csv"), row.names = FALSE)
    writeLines('{"spot_diameter_fullres": 177.5}', file.path(outs, "spatial", "scalefactors_json.json"))
    out <- tempfile()
    r <- suppressWarnings(suppressMessages(lamb_run(d, out_dir = out)))
    log <- readLines(file.path(out, "technical/run_log.txt"))
    expect_true("1 spots have no cluster and are left out." %in% log)
    expect_true(any(startsWith(log, paste0("Spots: ", n - 1, "; clusters: 2"))))
    expect_true(any(grepl("share of spots detecting one of", log)))
    expect_false(any(grepl("cells have no cluster", log)))
})

test_that("a panel without any epithelial marker leaves Epithelial_share empty, not 100%", {
    s <- sim_cells()
    m <- s$m[setdiff(rownames(s$m), .LAMB_EPITHELIAL), ]
    d <- tempfile(); write_10x(d, m, clusters = s$clusters)
    out <- tempfile()
    r <- suppressWarnings(suppressMessages(lamb_run(d, out_dir = out)))
    expect_true(all(is.na(r$calls$Epithelial_share)))
    calls <- utils::read.csv(file.path(out, "LAMB_calls.csv"), colClasses = "character")
    expect_true(all(calls$Epithelial_share == ""))
    log <- readLines(file.path(out, "technical/run_log.txt"))
    expect_true("No epithelial marker on this panel: Epithelial_share is left empty." %in% log)
    expect_false(any(grepl("detecting one of \\.$", log)))
})

test_that("a folder without clusters stops with a plain instruction", {
    s <- sim_cells()
    d <- tempfile(); write_10x(d, s$m)
    expect_error(suppressMessages(lamb_run(d, out_dir = tempfile())), "clusters =")
})

test_that("slide-level class is the majority class, 'mixed' below half, percentages add to 100", {
    calls <- data.frame(Cluster = c("a", "b", "c", "d"), Cells = c(60, 40, 30, 500),
                        Cell_type = c("epithelial", "epithelial", "epithelial", "not epithelial"),
                        LAMB_class = c("Basal-like", "LHS", "Too few cells", "Basal-like"), stringsAsFactors = FALSE)
    s <- .lamb_slide_summary(calls, function(...) NULL)
    expect_equal(s$call, "Basal-like")
    expect_equal(s$table[["Basal-like (%)"]], 60)
    expect_equal(sum(unlist(s$table[grepl("%", names(s$table))])), 100)
    calls$Cells[1:2] <- c(45, 40); calls$LAMB_class[3] <- "Apocrine"
    expect_match(.lamb_slide_summary(calls, function(...) NULL)$call, "^mixed: Basal-like")
})

test_that("Xenium with an H&E: cells are placed through the alignment and the slide table is written", {
    skip_if_not_installed("png")
    s <- sim_cells()
    d <- tempfile(); write_10x(d, s$m, xenium = TRUE, clusters = s$clusters)
    img <- tempfile(fileext = ".png")
    png::writePNG(array(0.9, c(300, 300, 3)), img)          # H&E exported at 1/4 of a 1200 px image
    al <- tempfile(fileext = ".csv")
    utils::write.table(matrix(c(1 / 0.2125, 0, 0, 0, 1 / 0.2125, 0, 0, 0, 1), 3, byrow = TRUE), al,
                       sep = ",", row.names = FALSE, col.names = FALSE)   # H&E pixel = 1 um
    out <- tempfile()
    r <- suppressWarnings(suppressMessages(lamb_run(d, out_dir = out, he_image = img, he_alignment = al, he_downsample = 4)))
    expect_true(file.exists(file.path(out, "LAMB_slide.csv")))
    expect_match(paste(readLines(file.path(out, "technical/run_log.txt")), collapse = " "), "100% of cells fall on the image")
    expect_error(suppressMessages(lamb_run(d, out_dir = tempfile(), he_image = img, he_alignment = al)), "he_downsample")
    expect_error(suppressMessages(lamb_run(d, out_dir = tempfile(), he_image = img)), "he_alignment")
})

test_that("a 'mixed' slide call names every class at 10% or more, ties included", {
    calls <- data.frame(Cluster = paste("c", 1:3), Cells = c(100, 100, 100),
                        Cell_type = "epithelial", LAMB_class = c("LHS", "Basal-like", "MA-Lo"))
    s <- .lamb_slide_summary(calls, function(...) NULL)$call
    expect_match(s, "^mixed: ")
    for (k in c("LHS", "Basal-like", "MA-Lo")) expect_match(s, k)
})

test_that("bulk: rounded TPM summing to one million per sample is read as TPM, not counts", {
    gl <- get_gene_length()
    g <- intersect(sim_genes(), gl$SYMBOL)
    set.seed(5)
    tpm <- matrix(stats::rgamma(length(g) * 3, 0.5, 0.01), nrow = length(g), dimnames = list(g, paste0("S", 1:3)))
    tpm <- round(sweep(tpm, 2, colSums(tpm), "/") * 1e6)
    expect_false(.lamb_is_counts(tpm))
    f <- tempfile(fileext = ".csv"); utils::write.csv(tpm, f)
    out <- tempfile()
    r <- suppressWarnings(suppressMessages(lamb_run(f, out_dir = out, run_pam50 = FALSE)))
    ref <- suppressWarnings(LAMBclassifier(log2(tpm + 1), run_pam50 = FALSE))
    expect_equal(r$calls$LAMB_class, ref$Pred_Class)
    expect_match(paste(readLines(file.path(out, "technical/run_log.txt")), collapse = " "), "rounded TPM")
    counts <- matrix(stats::rpois(length(g) * 3, c(200, 400, 800)), nrow = length(g), byrow = TRUE)
    expect_true(.lamb_is_counts(counts))
})

test_that("bulk: Ensembl IDs that mostly fail to map stop the run with examples", {
    gl <- get_gene_length()
    g <- intersect(sim_genes(), gl$SYMBOL)
    x <- matrix(stats::runif(length(g) * 2, 0, 10), nrow = length(g), dimnames = list(g, c("A", "B")))
    rownames(x) <- sprintf("ENSG%011d", seq_along(g) + 90000000)
    expect_error(.lamb_symbols(x, function(...) NULL, counts = FALSE), "Ensembl IDs match human genes")
})

test_that("bulk: mostly unrecognised gene names stop the run", {
    x <- matrix(stats::runif(600 * 2, 0, 10), nrow = 600,
                dimnames = list(c(utils::head(lamb_universe(), 200), paste0("probe_", 1:400)), c("A", "B")))
    expect_error(.lamb_check_names(x), "not recognised as human gene symbols")
    expect_silent(.lamb_check_names(x[1:200, ]))
})

test_that("NC and low-coverage notes point at the data, not the platform", {
    n <- .lamb_notes(c("NC", "LHS"), c("", "low_coverage"), "Squamous 4/9")
    expect_match(n[1], "too few genes in this data")
    expect_match(n[2], "in this data for Squamous 4/9")
    expect_false(any(grepl("platform", n)))
})

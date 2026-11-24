# Synthetic fixtures laid out as each vendor documents its output.

quiet <- function(...) NULL

rc_matrix <- function(genes = c("EPCAM", "KRT8", "AR", "ESR1", "TP63", "PTPRC"), n = 6, prefix = "c") {
    set.seed(11)
    m <- matrix(stats::rpois(length(genes) * n, 3), nrow = length(genes),
                dimnames = list(genes, paste0(prefix, seq_len(n), "-1")))
    Matrix::Matrix(m, sparse = TRUE)
}

rc_clusters_file <- function(dir, sub, ids, cl) {
    d <- file.path(dir, "analysis", "clustering", sub)
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(data.frame(Barcode = ids, Cluster = cl), file.path(d, "clusters.csv"), row.names = FALSE)
}

# 10x mtx folder; controls = extra non-"Gene Expression" features
rc_write_mtx <- function(mdir, m, controls = 0, gz = FALSE) {
    dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
    feat <- data.frame(id = paste0("ENSG", seq_len(nrow(m))), name = rownames(m), type = "Gene Expression")
    if (controls) {
        feat <- rbind(feat, data.frame(id = paste0("NegControlProbe_", seq_len(controls)),
                                       name = paste0("NegControlProbe_", seq_len(controls)), type = "Negative Control Probe"))
        m <- rbind(m, Matrix::Matrix(1, controls, ncol(m), sparse = TRUE))
    }
    Matrix::writeMM(methods::as(m, "CsparseMatrix"), file.path(mdir, "matrix.mtx"))
    utils::write.table(feat, file.path(mdir, "features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    writeLines(colnames(m), file.path(mdir, "barcodes.tsv"))
    if (gz) for (f in c("matrix.mtx", "features.tsv", "barcodes.tsv")) {
        p <- file.path(mdir, f)
        writeLines(readLines(p), gzfile(paste0(p, ".gz")))
        unlink(p)
    }
}

rc_write_xenium <- function(dir, m) {
    rc_write_mtx(file.path(dir, "cell_feature_matrix"), m, controls = 2, gz = TRUE)
    con <- gzfile(file.path(dir, "cells.csv.gz"), "w")
    utils::write.csv(data.frame(cell_id = colnames(m), x_centroid = seq_len(ncol(m)) * 10.5,
                                y_centroid = seq_len(ncol(m)) * 2.25), con, row.names = FALSE)
    close(con)
    writeLines('{"run_name": "x", "pixel_size": 0.2125, "panel_name": "p"}', file.path(dir, "experiment.xenium"))
    rc_clusters_file(dir, "gene_expression_graphclust", colnames(m), rep(1:2, length.out = ncol(m)))
}

# ---- HDF5 fixtures (rhdf5) ----

h5_attr <- function(f, path, name, value) {
    fid <- rhdf5::H5Fopen(f)
    on.exit(rhdf5::H5Fclose(fid))
    oid <- rhdf5::H5Oopen(fid, path)
    rhdf5::h5writeAttribute(value, oid, name)
    rhdf5::H5Oclose(oid)
}

h5_groups <- function(f, groups) for (g in groups) rhdf5::h5createGroup(f, g)

rc_write_10x_h5 <- function(f, m, controls = 0) {
    m <- methods::as(m, "CsparseMatrix")
    name <- rownames(m); type <- rep("Gene Expression", nrow(m))
    if (controls) {
        m <- rbind(m, Matrix::Matrix(1, controls, ncol(m), sparse = TRUE))
        name <- c(name, paste0("NegControlProbe_", seq_len(controls))); type <- c(type, rep("Negative Control Probe", controls))
    }
    m <- methods::as(m, "CsparseMatrix")
    rhdf5::h5createFile(f)
    h5_groups(f, c("matrix", "matrix/features"))
    rhdf5::h5write(as.integer(m@x), f, "matrix/data")
    rhdf5::h5write(as.integer(m@i), f, "matrix/indices")
    rhdf5::h5write(as.integer(m@p), f, "matrix/indptr")
    rhdf5::h5write(as.integer(dim(m)), f, "matrix/shape")
    rhdf5::h5write(colnames(m), f, "matrix/barcodes")
    rhdf5::h5write(name, f, "matrix/features/name")
    rhdf5::h5write(paste0("ENSG", seq_along(name)), f, "matrix/features/id")
    rhdf5::h5write(type, f, "matrix/features/feature_type")
    rhdf5::h5closeAll()
}

# AnnData group for a sparse matrix of cells x genes, given genes x cells `m`
h5ad_sparse <- function(f, path, m, fmt = "csr") {
    m <- methods::as(m, "CsparseMatrix")          # genes x cells, CSC == cells x genes CSR
    if (fmt == "csc") m <- methods::as(Matrix::t(m), "CsparseMatrix")
    rhdf5::h5createGroup(f, path)
    rhdf5::h5write(as.numeric(m@x), f, paste0(path, "/data"))
    rhdf5::h5write(as.integer(m@i), f, paste0(path, "/indices"))
    rhdf5::h5write(as.integer(m@p), f, paste0(path, "/indptr"))
    shape <- if (fmt == "csr") rev(dim(m)) else dim(m)   # (n_obs, n_var)
    h5_attr(f, path, "encoding-type", paste0(fmt, "_matrix"))
    h5_attr(f, path, "shape", as.integer(shape))
}

h5ad_df <- function(f, path, index, cols = list(), cat_new = list(), cat_old = list()) {
    rhdf5::h5createGroup(f, path)
    rhdf5::h5write(index, f, paste0(path, "/_index"))
    h5_attr(f, path, "_index", "_index")
    h5_attr(f, path, "encoding-type", "dataframe")
    for (k in names(cols)) rhdf5::h5write(cols[[k]], f, paste0(path, "/", k))
    for (k in names(cat_new)) {            # anndata >= 0.8
        v <- factor(cat_new[[k]])
        rhdf5::h5createGroup(f, paste0(path, "/", k))
        rhdf5::h5write(levels(v), f, paste0(path, "/", k, "/categories"))
        codes <- as.integer(v) - 1L; codes[is.na(codes)] <- -1L
        rhdf5::h5write(codes, f, paste0(path, "/", k, "/codes"))
        h5_attr(f, paste0(path, "/", k), "encoding-type", "categorical")
    }
    if (length(cat_old)) rhdf5::h5createGroup(f, paste0(path, "/__categories"))
    for (k in names(cat_old)) {            # anndata 0.7
        v <- factor(cat_old[[k]])
        rhdf5::h5write(levels(v), f, paste0(path, "/__categories/", k))
        rhdf5::h5write(as.integer(v) - 1L, f, paste0(path, "/", k))
    }
}

lognorm <- function(m) { x <- methods::as(m, "CsparseMatrix"); x@x <- log1p(x@x / 7); x }

# ---- detection ----

test_that("detection picks the platform from the files present", {
    m <- rc_matrix()
    d <- tempfile(); rc_write_mtx(file.path(d, "outs", "filtered_feature_bc_matrix"), m)
    expect_equal(.lamb_detect_cells(d), "10x_mtx")
    expect_equal(.lamb_detect_cells(file.path(d, "outs", "filtered_feature_bc_matrix")), "10x_mtx")
    x <- tempfile(); rc_write_xenium(x, m)
    expect_equal(.lamb_detect_cells(x), "xenium")
    e <- tempfile(); dir.create(e)
    expect_null(.lamb_detect_cells(e))
    # a folder merely *named* like a vendor output is not detected
    v <- file.path(tempfile(), "visium_outs"); dir.create(file.path(v, "spatial"), recursive = TRUE)
    expect_null(.lamb_detect_cells(v))
    b <- tempfile(fileext = ".csv"); utils::write.csv(as.matrix(m), b)
    expect_null(.lamb_detect_cells(b))
    r <- tempfile(fileext = ".rds"); saveRDS(as.matrix(m), r)
    expect_null(.lamb_detect_cells(r))
})

# ---- 10x mtx and Xenium: identical to the frozen reader ----
test_that("10x mtx: same output as the original .lamb_load_cells", {
    m <- rc_matrix()
    d <- tempfile(); rc_write_mtx(file.path(d, "filtered_feature_bc_matrix"), m)
    rc_clusters_file(d, "gene_expression_graphclust", colnames(m), c(1, 1, 2, 2, 3, 3))
    old <- .lamb_load_cells(d, "single_cell", quiet)
    new <- .lamb_read_cells(d, .lamb_detect_cells(d), quiet)
    expect_identical(.lamb_detect_cells(d), "10x_mtx")
    expect_identical(new$counts, old$counts)
    expect_identical(new$clusters, old$clusters)
    expect_identical(new$xy, old$xy)
    expect_identical(new$pixel_size, old$pixel_size)
    expect_equal(new$unit, "cell"); expect_equal(new$platform, "10x_mtx")
    expect_equal(dim(new$counts), dim(m))
    expect_equal(as.matrix(new$counts), as.matrix(m))
})

test_that("Xenium: same output as the original .lamb_load_cells; controls dropped; xy in microns", {
    m <- rc_matrix()
    d <- tempfile(); rc_write_xenium(d, m)
    msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
    old <- .lamb_load_cells(d, "xenium", quiet)
    new <- .lamb_read_cells(d, "xenium", say)
    expect_identical(.lamb_detect_cells(d), "xenium")
    for (k in c("counts", "clusters", "xy", "pixel_size")) expect_identical(new[[k]], old[[k]])
    expect_equal(rownames(new$counts), rownames(m))
    expect_equal(new$xy$x, seq_len(ncol(m)) * 10.5)
    expect_equal(new$clusters[["c3-1"]], "cluster 1")
    expect_true(any(grepl("2 control / codeword features dropped", msgs)))
    expect_true(any(grepl("6 cells, 6 genes, 2 groups, positions in microns", msgs)))
})

test_that("duplicated gene symbols are summed (sparse-safe)", {
    m <- rc_matrix(genes = c("KRT8", "AR", "KRT8", "ESR1"))
    d <- tempfile(); rc_write_mtx(file.path(d, "filtered_feature_bc_matrix"), m)
    new <- .lamb_read_cells(d, "10x_mtx", quiet)
    expect_equal(rownames(new$counts), c("KRT8", "AR", "ESR1"))
    expect_equal(as.matrix(new$counts), rowsum(as.matrix(m), rownames(m), reorder = FALSE))
    expect_s4_class(new$counts, "dgCMatrix")
})

test_that("older Cell Ranger clusters (graphclust/) are found when gene_expression_graphclust/ is absent", {
    m <- rc_matrix()
    d <- tempfile(); rc_write_mtx(file.path(d, "outs", "filtered_feature_bc_matrix"), m)
    rc_clusters_file(file.path(d, "outs"), "graphclust", colnames(m), c(1, 1, 1, 2, 2, 2))
    new <- .lamb_read_cells(d, "10x_mtx", quiet)
    expect_equal(unname(new$clusters[colnames(m)]), rep(c("cluster 1", "cluster 2"), each = 3))
})

# ---- 10x .h5 ----

for (backend in c("rhdf5", "hdf5r")) {
    test_that(paste0("10x .h5 (", backend, "): matrix, symbols, feature types; Xenium .h5 accepted"), {
        skip_if_not_installed("rhdf5")
        skip_if_not_installed(backend)
        withr_opt <- options(LAMBclassifier.h5_backend = backend); on.exit(options(withr_opt))
        m <- rc_matrix()
        d <- tempfile(); dir.create(file.path(d, "outs"), recursive = TRUE)
        f <- file.path(d, "outs", "filtered_feature_bc_matrix.h5")
        rc_write_10x_h5(f, m, controls = 1)
        rc_clusters_file(file.path(d, "outs"), "gene_expression_graphclust", colnames(m), rep(1:3, 2))
        expect_equal(.lamb_detect_cells(f), "10x_h5")
        expect_equal(.lamb_detect_cells(d), "10x_h5")
        a <- .lamb_read_cells(f, "10x_h5", quiet)
        b <- .lamb_read_cells(d, "10x_h5", quiet)
        expect_equal(as.matrix(a$counts), as.matrix(m))
        expect_identical(a$counts, b$counts)
        expect_equal(a$clusters[["c4-1"]], "cluster 1")
        # Xenium folder that only ships cell_feature_matrix.h5
        x <- tempfile(); rc_write_xenium(x, m); unlink(file.path(x, "cell_feature_matrix"), recursive = TRUE)
        rc_write_10x_h5(file.path(x, "cell_feature_matrix.h5"), m, controls = 2)
        expect_equal(.lamb_detect_cells(x), "xenium")
        xr <- .lamb_read_cells(x, "xenium", quiet)
        expect_equal(as.matrix(xr$counts), as.matrix(m))
        expect_equal(nrow(xr$xy), ncol(m))
    })
}

# ---- AnnData ----

for (backend in c("rhdf5", "hdf5r")) {
    test_that(paste0("h5ad (", backend, "): counts layer preferred over log-normalised X; new categoricals"), {
        skip_if_not_installed("rhdf5")
        skip_if_not_installed(backend)
        withr_opt <- options(LAMBclassifier.h5_backend = backend); on.exit(options(withr_opt))
        m <- rc_matrix()
        f <- tempfile(fileext = ".h5ad"); rhdf5::h5createFile(f)
        h5ad_sparse(f, "X", lognorm(m))
        rhdf5::h5createGroup(f, "layers"); h5ad_sparse(f, "layers/counts", m)
        h5ad_df(f, "obs", colnames(m), cat_new = list(leiden = c("0", "0", "1", "1", NA, "2")))
        h5ad_df(f, "var", rownames(m))
        rhdf5::h5closeAll()
        expect_equal(.lamb_detect_cells(f), "h5ad")
        msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
        r <- .lamb_read_cells(f, "h5ad", say)
        expect_equal(as.matrix(r$counts), as.matrix(m))
        expect_equal(unname(r$clusters), c("0", "0", "1", "1", NA, "2"))
        expect_equal(names(r$clusters), colnames(m))
        expect_true(any(grepl("layers/counts", msgs)))
        expect_true(any(grepl("obs column 'leiden'", msgs)))
    })

    test_that(paste0("h5ad (", backend, "): dense counts X, Ensembl index with feature_name, old categoricals"), {
        skip_if_not_installed("rhdf5")
        skip_if_not_installed(backend)
        withr_opt <- options(LAMBclassifier.h5_backend = backend); on.exit(options(withr_opt))
        m <- rc_matrix()
        f <- tempfile(fileext = ".h5ad"); rhdf5::h5createFile(f)
        rhdf5::h5write(as.matrix(m) * 1.0, f, "X")      # R genes x cells -> HDF5 (n_obs, n_var)
        h5ad_df(f, "obs", colnames(m), cat_old = list(louvain = c("a", "b", "a", "b", "a", "b")))
        h5ad_df(f, "var", paste0("ENSG0000000000", seq_len(nrow(m))), cat_new = list(feature_name = rownames(m)))
        rhdf5::h5closeAll()
        r <- .lamb_read_cells(f, "h5ad", quiet)
        expect_equal(rownames(r$counts), rownames(m))
        expect_equal(as.matrix(r$counts), as.matrix(m), ignore_attr = TRUE)
        expect_equal(unname(r$clusters), c("a", "b", "a", "b", "a", "b"))
    })

    test_that(paste0("h5ad (", backend, "): raw/X in CSC with raw/var; no cluster column gives NULL"), {
        skip_if_not_installed("rhdf5")
        skip_if_not_installed(backend)
        withr_opt <- options(LAMBclassifier.h5_backend = backend); on.exit(options(withr_opt))
        m <- rc_matrix()
        f <- tempfile(fileext = ".h5ad"); rhdf5::h5createFile(f)
        h5ad_sparse(f, "X", lognorm(m[1:4, ]))                # X: HVG subset, normalised
        rhdf5::h5createGroup(f, "raw"); h5ad_sparse(f, "raw/X", m, fmt = "csc")
        h5ad_df(f, "raw/var", rownames(m))
        h5ad_df(f, "obs", colnames(m), cols = list(n_genes = rep(5L, ncol(m))))
        h5ad_df(f, "var", rownames(m)[1:4])
        rhdf5::h5closeAll()
        msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
        r <- .lamb_read_cells(f, "h5ad", say)
        expect_equal(as.matrix(r$counts), as.matrix(m))
        expect_null(r$clusters)
        expect_true(any(grepl("Supply them with clusters = 'my_clusters.csv'", msgs)))
    })
}

test_that("h5ad with only log-normalised values stops with a plain message", {
    skip_if_not_installed("rhdf5")
    m <- rc_matrix()
    f <- tempfile(fileext = ".h5ad"); rhdf5::h5createFile(f)
    h5ad_sparse(f, "X", lognorm(m))
    rhdf5::h5createGroup(f, "layers"); h5ad_sparse(f, "layers/spliced", lognorm(m))
    h5ad_df(f, "obs", colnames(m)); h5ad_df(f, "var", rownames(m))
    rhdf5::h5closeAll()
    expect_error(.lamb_read_cells(f, "h5ad", quiet), "needs the raw counts")
    expect_error(.lamb_read_cells(f, "h5ad", quiet), "other layers, not read: spliced")
})

# anndata >= 0.10 (as written by anndata 0.10.9): obs index named 'barcode',
# nullable-integer / nullable-boolean columns (values + mask), counts layer in CSC.
h5ad_nullable <- function(f, path, values, mask, type) {
    rhdf5::h5createGroup(f, path)
    rhdf5::h5write(values, f, paste0(path, "/values"))
    rhdf5::h5write(mask, f, paste0(path, "/mask"))
    h5_attr(f, path, "encoding-type", type)
}

for (backend in c("rhdf5", "hdf5r")) {
    test_that(paste0("h5ad (", backend, "): anndata 0.10 nullable obs columns, CSC counts layer, duplicated cell names"), {
        skip_if_not_installed("rhdf5")
        skip_if_not_installed(backend)
        withr_opt <- options(LAMBclassifier.h5_backend = backend); on.exit(options(withr_opt))
        m <- rc_matrix(n = 40)
        cells <- colnames(m); cells[40] <- cells[1]
        f <- tempfile(fileext = ".h5ad"); rhdf5::h5createFile(f)
        h5ad_sparse(f, "X", lognorm(m), fmt = "csc")
        rhdf5::h5createGroup(f, "layers"); h5ad_sparse(f, "layers/counts", m, fmt = "csc")
        rhdf5::h5createGroup(f, "obs")
        rhdf5::h5write(cells, f, "obs/barcode")
        h5_attr(f, "obs", "_index", "barcode"); h5_attr(f, "obs", "encoding-type", "dataframe")
        h5ad_nullable(f, "obs/n_umi", seq_len(40), c(TRUE, rep(FALSE, 39)), "nullable-integer")
        h5ad_nullable(f, "obs/flag", rep(c(TRUE, FALSE), 20), c(TRUE, rep(FALSE, 39)), "nullable-boolean")
        h5ad_df(f, "var", rownames(m))
        rhdf5::h5closeAll()
        msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
        r <- .lamb_read_cells(f, "h5ad", say)
        expect_equal(as.matrix(r$counts), as.matrix(m), ignore_attr = TRUE)
        expect_equal(colnames(r$counts), make.unique(cells))
        expect_true(any(grepl("1 duplicated cell names", msgs)))
        expect_equal(r$meta$flag, c(NA, rep(c("FALSE", "TRUE"), 19), "FALSE"))
        expect_false("n_umi" %in% names(r$meta))       # one value per cell: not a group column
        g <- .lamb_groups("flag", NULL, colnames(r$counts), quiet, "h5ad", r$meta)
        expect_equal(g[1:3], c(NA, "FALSE", "TRUE"))
    })

    test_that(paste0("h5ad (", backend, "): raw/var with gene IDs only takes symbols from var, matched by ID"), {
        skip_if_not_installed("rhdf5")
        skip_if_not_installed(backend)
        withr_opt <- options(LAMBclassifier.h5_backend = backend); on.exit(options(withr_opt))
        m <- rc_matrix()
        ids <- paste0("ENSG9900000000", seq_len(nrow(m)))
        f <- tempfile(fileext = ".h5ad"); rhdf5::h5createFile(f)
        h5ad_sparse(f, "X", lognorm(m[c(2, 4, 5), ]))         # X: subset of genes, normalised
        rhdf5::h5createGroup(f, "raw"); h5ad_sparse(f, "raw/X", m)
        h5ad_df(f, "raw/var", ids)                            # index only
        h5ad_df(f, "var", ids[c(2, 4, 5)], cat_new = list(feature_name = rownames(m)[c(2, 4, 5)]))
        h5ad_df(f, "obs", colnames(m))
        rhdf5::h5closeAll()
        msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
        r <- .lamb_read_cells(f, "h5ad", say)
        expect_equal(rownames(r$counts), rownames(m)[c(2, 4, 5)])   # the fake IDs map to nothing
        expect_equal(as.matrix(r$counts), as.matrix(m)[c(2, 4, 5), ], ignore_attr = TRUE)
        expect_true(any(grepl("var/feature_name, matched to raw/var by gene ID", msgs)))
    })
}

test_that("a 'no clusters' hint is not printed when the user names the clusters", {
    skip_if_not_installed("rhdf5")
    m <- rc_matrix()
    f <- tempfile(fileext = ".h5ad"); rhdf5::h5createFile(f)
    h5ad_sparse(f, "X", m)
    h5ad_df(f, "obs", colnames(m)); h5ad_df(f, "var", rownames(m))
    rhdf5::h5closeAll()
    msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
    .lamb_read_cells(f, "h5ad", say, user_clusters = TRUE)
    expect_false(any(grepl("no cell clusters", msgs)))
    .lamb_read_cells(f, "h5ad", say)
    expect_true(any(grepl("no cell clusters", msgs)))
})

# ---- Visium ----

rc_write_visium <- function(d, m, header = TRUE, h5 = FALSE) {
    outs <- file.path(d, "outs"); sp <- file.path(outs, "spatial")
    dir.create(sp, recursive = TRUE)
    if (h5) rc_write_10x_h5(file.path(outs, "filtered_feature_bc_matrix.h5"), m) else
        rc_write_mtx(file.path(outs, "filtered_feature_bc_matrix"), m, gz = TRUE)
    n <- ncol(m)
    pos <- data.frame(barcode = colnames(m), in_tissue = c(rep(1L, n - 1), 0L), array_row = seq_len(n),
                      array_col = seq_len(n), pxl_row_in_fullres = seq_len(n) * 100, pxl_col_in_fullres = seq_len(n) * 200)
    if (header) utils::write.csv(pos, file.path(sp, "tissue_positions.csv"), row.names = FALSE) else
        utils::write.table(pos, file.path(sp, "tissue_positions_list.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
    writeLines('{"tissue_hires_scalef": 0.17, "tissue_lowres_scalef": 0.05, "fiducial_diameter_fullres": 160.0, "spot_diameter_fullres": 110.0}',
               file.path(sp, "scalefactors_json.json"))
    rc_clusters_file(outs, "graphclust", colnames(m), rep(1:2, length.out = n))
}

test_that("Visium (newer layout): in-tissue spots only, pixels converted with 55 um / spot_diameter_fullres", {
    m <- rc_matrix()
    d <- tempfile(); rc_write_visium(d, m)
    expect_equal(.lamb_detect_cells(d), "visium")
    expect_equal(.lamb_detect_cells(file.path(d, "outs")), "visium")
    r <- .lamb_read_cells(d, "visium", quiet)
    expect_equal(r$unit, "spot")
    expect_equal(ncol(r$counts), ncol(m) - 1)
    expect_equal(r$pixel_size, 55 / 110)
    expect_equal(r$xy$x, (1:5) * 200 * 0.5)       # pxl_col -> x
    expect_equal(r$xy$y, (1:5) * 100 * 0.5)       # pxl_row -> y
    expect_equal(r$clusters[["c2-1"]], "cluster 2")
})

test_that("Visium (older layout): headerless tissue_positions_list.csv and an .h5 matrix", {
    skip_if_not_installed("rhdf5")
    m <- rc_matrix()
    d <- tempfile(); rc_write_visium(d, m, header = FALSE, h5 = TRUE)
    expect_equal(.lamb_detect_cells(d), "visium")
    r <- .lamb_read_cells(d, "visium", quiet)
    expect_equal(ncol(r$counts), ncol(m) - 1)
    expect_equal(r$xy$cell_id, colnames(m)[1:5])
    expect_equal(r$xy$x[2], 400 * 0.5)
})

test_that("Visium as downloaded from 10x: <sample>_filtered_feature_bc_matrix.h5 next to spatial/ and analysis/", {
    skip_if_not_installed("rhdf5")
    m <- rc_matrix()
    d <- tempfile(); rc_write_visium(d, m, header = FALSE, h5 = TRUE)
    flat <- file.path(d, "outs")
    h5 <- file.path(flat, "V1_Breast_Cancer_Block_A_Section_1_filtered_feature_bc_matrix.h5")
    file.rename(file.path(flat, "filtered_feature_bc_matrix.h5"), h5)
    expect_equal(.lamb_detect_cells(flat), "visium")
    expect_equal(.lamb_detect_cells(h5), "visium")
    for (p in c(flat, h5)) {
        r <- .lamb_read_cells(p, "visium", quiet)
        expect_equal(r$unit, "spot")
        expect_equal(nrow(r$xy), ncol(m) - 1)
        expect_equal(r$clusters[["c2-1"]], "cluster 2")
    }
    # a raw matrix beside it is not taken for the filtered one
    file.copy(h5, file.path(flat, "V1_raw_feature_bc_matrix.h5"))
    expect_equal(.lamb_find_h5(flat), h5)
})

test_that("Visium: the matrix .h5 or filtered_feature_bc_matrix/ given on its own is read with its positions", {
    skip_if_not_installed("rhdf5")
    m <- rc_matrix()
    d <- tempfile(); rc_write_visium(d, m, h5 = TRUE)
    h5 <- file.path(d, "outs", "filtered_feature_bc_matrix.h5")
    expect_equal(.lamb_detect_cells(h5), "visium")
    r <- .lamb_read_cells(h5, "visium", quiet)
    expect_equal(r$unit, "spot")
    expect_equal(r$xy$x[2], 400 * 0.5)
    d2 <- tempfile(); rc_write_visium(d2, m)
    md <- file.path(d2, "outs", "filtered_feature_bc_matrix")
    expect_equal(.lamb_detect_cells(md), "visium")
    expect_equal(nrow(.lamb_read_cells(md, "visium", quiet)$xy), ncol(m) - 1)
    # a plain Cell Ranger .h5 (no spatial/ beside it) stays single-cell
    d3 <- tempfile(); dir.create(d3); rc_write_10x_h5(file.path(d3, "filtered_feature_bc_matrix.h5"), m)
    expect_equal(.lamb_detect_cells(file.path(d3, "filtered_feature_bc_matrix.h5")), "10x_h5")
})

# ---- Visium HD ----

# Space Ranger 3 bin folder: filtered matrix (.h5 and folder), parquet positions
# with in_tissue = 0 rows, scalefactors and analysis/clustering.
rc_write_hd_bin <- function(root, bin, m) {
    bd <- file.path(root, "binned_outputs", bin)
    rc_write_mtx(file.path(bd, "filtered_feature_bc_matrix"), m, gz = TRUE)
    if (requireNamespace("rhdf5", quietly = TRUE)) rc_write_10x_h5(file.path(bd, "filtered_feature_bc_matrix.h5"), m)
    dir.create(file.path(bd, "spatial"))
    n <- ncol(m)
    arrow::write_parquet(data.frame(barcode = c(colnames(m), "s_off-1"), in_tissue = c(rep(1L, n), 0L),
                                    array_row = 0:n, array_col = 0:n, pxl_row_in_fullres = (0:n) * 10,
                                    pxl_col_in_fullres = (0:n) * 20),
                         file.path(bd, "spatial", "tissue_positions.parquet"))
    um <- as.numeric(sub("^square_0*([0-9]+)um$", "\\1", bin))
    writeLines(sprintf('{"spot_diameter_fullres": %s, "bin_size_um": %s, "microns_per_pixel": 0.3251}', um / 0.3251, um),
               file.path(bd, "spatial", "scalefactors_json.json"))
    rc_clusters_file(bd, "gene_expression_graphclust", colnames(m)[-1], rep(1:2, length.out = n - 1))
    bd
}

test_that("Visium HD: only square_016um kept, it is used and said; 2 um is never picked by itself", {
    skip_if_not_installed("arrow")
    d <- tempfile()
    rc_write_hd_bin(d, "square_002um", rc_matrix(n = 3, prefix = "s_002um_"))
    rc_write_hd_bin(d, "square_016um", rc_matrix(prefix = "s_016um_"))
    expect_equal(.lamb_detect_cells(d), "visium_hd")
    msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
    r <- .lamb_read_cells(d, "visium_hd", say)
    expect_equal(ncol(r$counts), 6)
    expect_true(any(grepl("no usable square_008um folder; using square_016um", msgs)))
    expect_equal(r$xy$x[2], 20 * 0.3251)
    expect_equal(ncol(.lamb_read_cells(d, "visium_hd", quiet, hd_bin = "square_002um")$counts), 3)
    expect_error(.lamb_read_cells(d, "visium_hd", quiet, hd_bin = "square_004um"), "Bins present: square_002um, square_016um")
})

test_that("Visium HD: a bin's matrix .h5 or folder given on its own reads that bin with positions", {
    skip_if_not_installed("arrow")
    skip_if_not_installed("rhdf5")
    d <- tempfile()
    bd <- rc_write_hd_bin(d, "square_008um", rc_matrix(prefix = "s_008um_"))
    for (p in c(file.path(bd, "filtered_feature_bc_matrix.h5"), file.path(bd, "filtered_feature_bc_matrix"))) {
        expect_equal(.lamb_detect_cells(p), "visium_hd")
        r <- .lamb_read_cells(p, "visium_hd", quiet)
        expect_equal(r$unit, "bin")
        expect_equal(nrow(r$xy), 6)
        expect_equal(length(r$clusters), 5)
    }
})

test_that("Visium HD: square_008um used by default, parquet positions in microns, unit = bin", {
    skip_if_not_installed("arrow")
    m <- rc_matrix(prefix = "s_008um_0000")
    d <- tempfile()
    for (b in c("square_002um", "square_008um")) {
        bd <- file.path(d, "outs", "binned_outputs", b)
        mm <- if (b == "square_008um") m else rc_matrix(n = 3, prefix = "s_002um_")
        rc_write_mtx(file.path(bd, "filtered_feature_bc_matrix"), mm, gz = TRUE)
        dir.create(file.path(bd, "spatial"))
        n <- ncol(mm)
        arrow::write_parquet(data.frame(barcode = colnames(mm), in_tissue = rep(1L, n), array_row = seq_len(n),
                                        array_col = seq_len(n), pxl_row_in_fullres = seq_len(n) * 10,
                                        pxl_col_in_fullres = seq_len(n) * 20),
                             file.path(bd, "spatial", "tissue_positions.parquet"))
        writeLines('{"spot_diameter_fullres": 29.1, "bin_size_um": 8.0, "microns_per_pixel": 0.2738, "tissue_hires_scalef": 0.08}',
                   file.path(bd, "spatial", "scalefactors_json.json"))
        rc_clusters_file(bd, "gene_expression_graphclust", colnames(mm), rep(1:2, length.out = n))
    }
    expect_equal(.lamb_detect_cells(d), "visium_hd")
    r <- .lamb_read_cells(d, "visium_hd", quiet)
    expect_equal(r$unit, "bin")
    expect_equal(ncol(r$counts), ncol(m))
    expect_equal(colnames(r$counts), colnames(m))
    expect_equal(r$xy$x, seq_len(ncol(m)) * 20 * 0.2738)
    expect_equal(r$clusters[[colnames(m)[2]]], "cluster 2")
    r2 <- .lamb_read_cells(d, "visium_hd", quiet, hd_bin = "square_002um")
    expect_equal(ncol(r2$counts), 3)
})

# ---- CosMx ----

test_that("CosMx flat files: controls and background dropped, ids fov_cellID, pixels to microns, no groups", {
    d <- tempfile(); dir.create(d)
    ex <- data.frame(fov = c(1, 1, 1, 2, 2), cell_ID = c(0, 1, 2, 1, 2),
                     EPCAM = c(9, 1, 0, 3, 4), KRT8 = c(9, 2, 2, 0, 1), NegPrb1 = 1, NegPrb2 = 0, SystemControl1 = 5)
    con <- gzfile(file.path(d, "Run5_S1_exprMat_file.csv.gz"), "w"); utils::write.csv(ex, con, row.names = FALSE); close(con)
    md <- data.frame(fov = c(1, 1, 2, 2), cell_ID = c(1, 2, 1, 2), Area = 100,
                     CenterX_global_px = c(1000, 2000, 3000, 4000), CenterY_global_px = c(10, 20, 30, 40))
    utils::write.csv(md, file.path(d, "Run5_S1_metadata_file.csv"), row.names = FALSE)
    expect_equal(.lamb_detect_cells(d), "cosmx")
    msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
    r <- .lamb_read_cells(d, "cosmx", say)
    expect_equal(rownames(r$counts), c("EPCAM", "KRT8"))
    expect_equal(colnames(r$counts), c("1_1", "1_2", "2_1", "2_2"))
    expect_equal(unname(as.matrix(r$counts)["EPCAM", ]), c(1, 0, 3, 4))
    expect_equal(r$xy$x, c(1000, 2000, 3000, 4000) * 0.12028)
    expect_null(r$clusters)
    expect_true(any(msgs == "This CosMx run has no cell clusters. Supply them with clusters = 'my_clusters.csv' (two columns: cell id, cluster)."))
    expect_true(any(grepl("3 negative-probe / system-control", msgs)))
})

test_that("CosMx: each exprMat file is paired with the metadata of its own run; a folder of several runs asks for one", {
    d <- tempfile(); dir.create(d)
    for (run in c("Lung9_Rep1", "Lung5_Rep2")) {
        k <- if (run == "Lung9_Rep1") 1 else 2
        ex <- data.frame(fov = c(1, 1, 2), cell_ID = c(0, 1, 1), EPCAM = c(9, k, 3), KRT8 = 1, NegPrb1 = 0)
        utils::write.csv(ex, file.path(d, paste0(run, "_exprMat_file.csv")), row.names = FALSE)
        # as exported: fov quoted, a centre in global pixels
        writeLines(c('"fov","cell_ID","Area","CenterX_global_px","CenterY_global_px"',
                     paste0('"1",1,100,', 1000 * k, ',10'), paste0('"2",1,100,', 2000 * k, ',20')),
                   file.path(d, paste0(run, "_metadata_file.csv")))
    }
    expect_equal(.lamb_detect_cells(d), "cosmx")
    expect_error(.lamb_read_cells(d, "cosmx", quiet), "holds 2 CosMx runs.*Give one of these files")
    f <- file.path(d, "Lung5_Rep2_exprMat_file.csv")
    expect_equal(.lamb_detect_cells(f), "cosmx")
    r <- .lamb_read_cells(f, "cosmx", quiet)
    expect_equal(colnames(r$counts), c("1_1", "2_1"))
    expect_equal(unname(as.matrix(r$counts)["EPCAM", ]), c(2, 3))
    expect_equal(r$xy$x, c(2000, 4000) * 0.12028)
    # an exprMat file without the metadata of its run is not taken for CosMx
    unlink(file.path(d, "Lung9_Rep1_metadata_file.csv"))
    expect_null(.lamb_detect_cells(file.path(d, "Lung9_Rep1_exprMat_file.csv")))
})

# ---- MERSCOPE ----

test_that("MERSCOPE files with a run prefix (<run>_cell_by_gene.csv) as public exports name them", {
    d <- tempfile(); dir.create(d)
    for (reg in c("1017_region_0_", "1017_region_1_")) {
        k <- if (reg == "1017_region_0_") 1 else 5
        writeLines(c("cell,EPCAM,KRT8,Blank-0,Blank-1", paste(0:2, c(3, 0, 1) * k, 1, 0, 1, sep = ",")),
                   file.path(d, paste0(reg, "cell_by_gene.csv")))
        # first header empty, fov column empty, as in the export
        writeLines(c(",fov,volume,center_x,center_y,min_x,min_y,max_x,max_y,barcodeCount",
                     paste0(0:2, ",,500.1,", c(-2.2, 20.5, 30.5) * k, ",", c(1.5, 2.5, 3.5), ",0,0,1,1,")),
                   file.path(d, paste0(reg, "cell_metadata.csv")))
    }
    expect_equal(.lamb_detect_cells(d), "merscope")
    expect_error(.lamb_read_cells(d, "merscope", quiet), "holds 2 MERSCOPE runs")
    f <- file.path(d, "1017_region_1_cell_by_gene.csv")
    expect_equal(.lamb_detect_cells(f), "merscope")
    r <- .lamb_read_cells(f, "merscope", quiet)
    expect_equal(colnames(r$counts), c("0", "1", "2"))
    expect_equal(rownames(r$counts), c("EPCAM", "KRT8"))
    expect_equal(unname(as.matrix(r$counts)["EPCAM", ]), c(15, 0, 5))
    expect_equal(r$xy$x, c(-2.2, 20.5, 30.5) * 5)
    # one run alone in its folder is found from the folder
    d2 <- tempfile(); dir.create(d2)
    file.copy(file.path(d, c("1017_region_0_cell_by_gene.csv", "1017_region_0_cell_metadata.csv")), d2)
    expect_equal(.lamb_detect_cells(d2), "merscope")
    expect_equal(ncol(.lamb_read_cells(d2, "merscope", quiet)$counts), 3)
})

test_that("MERSCOPE: Blank-* dropped, 19-digit cell ids kept exact, centres already microns", {
    d <- tempfile(); dir.create(d)
    ids <- c("1234567890123456789", "1234567890123456790", "1234567890123456791")
    writeLines(c("cell,EPCAM,KRT8,Blank-0001,Blank-0002",
                 paste(ids, c(3, 0, 1), c(1, 1, 0), 0, 1, sep = ",")), file.path(d, "cell_by_gene.csv"))
    writeLines(c(",fov,volume,center_x,center_y,min_x",
                 paste(ids, 0, 500, c(10.5, 20.5, 30.5), c(1.5, 2.5, 3.5), 0, sep = ",")), file.path(d, "cell_metadata.csv"))
    expect_equal(.lamb_detect_cells(d), "merscope")
    expect_equal(.lamb_detect_cells(file.path(d, "cell_by_gene.csv")), "merscope")
    r <- .lamb_read_cells(d, "merscope", quiet)
    expect_equal(colnames(r$counts), ids)
    expect_equal(rownames(r$counts), c("EPCAM", "KRT8"))
    expect_equal(r$xy$cell_id, ids)
    expect_equal(r$xy$x, c(10.5, 20.5, 30.5))
    expect_null(r$clusters)
    expect_match(.lamb_no_clusters_text("merscope"), "This MERSCOPE run has no cell clusters")
})

# ---- Seurat / SingleCellExperiment ----

test_that("Seurat and SingleCellExperiment .rds: detected and read as before", {
    m <- rc_matrix(n = 8)
    if (requireNamespace("SeuratObject", quietly = TRUE)) {
        s <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = m))
        SeuratObject::Idents(s) <- factor(rep(c("A", "B"), 4))
        f <- tempfile(fileext = ".rds"); saveRDS(s, f)
        expect_equal(.lamb_detect_cells(f), "seurat")
        old <- .lamb_load_cells(f, "single_cell", quiet)
        new <- .lamb_read_cells(f, "seurat", quiet)
        expect_identical(new$counts, old$counts)
        expect_identical(new$clusters, old$clusters)
        expect_null(new$xy)
    }
    if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        sce <- SingleCellExperiment::SingleCellExperiment(list(counts = m))
        SingleCellExperiment::colLabels(sce) <- rep(c("x", "y"), 4)
        f <- tempfile(fileext = ".rds"); saveRDS(sce, f)
        expect_equal(.lamb_detect_cells(f), "sce")
        new <- .lamb_read_cells(f, "sce", quiet)
        expect_equal(unname(new$clusters), rep(c("x", "y"), 4))
        expect_equal(as.matrix(new$counts), as.matrix(m))
    }
    if (!requireNamespace("SeuratObject", quietly = TRUE) && !requireNamespace("SingleCellExperiment", quietly = TRUE))
        skip("neither SeuratObject nor SingleCellExperiment installed")
})

test_that("Seurat v5 merged object: split counts layers are joined, every cell kept", {
    skip_if_not_installed("SeuratObject")
    a <- rc_matrix(n = 6, prefix = "a"); b <- rc_matrix(genes = c("EPCAM", "KRT8", "AR", "GATA3"), n = 4, prefix = "b")
    s <- suppressWarnings(merge(SeuratObject::CreateSeuratObject(a, project = "s1"),
                                SeuratObject::CreateSeuratObject(b, project = "s2")))
    s$cell_type <- rep(c("x", "y"), 5)
    expect_gt(length(SeuratObject::Layers(s[["RNA"]], search = "counts")), 1)
    f <- tempfile(fileext = ".rds"); saveRDS(s, f)
    msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
    r <- .lamb_read_cells(f, "seurat", say)
    expect_equal(ncol(r$counts), 10)
    expect_equal(colnames(r$counts), colnames(s))
    expect_equal(sort(rownames(r$counts)), sort(union(rownames(a), rownames(b))))
    expect_equal(unname(r$counts["GATA3", colnames(b)]), unname(b["GATA3", ]))
    expect_equal(sum(r$counts), sum(a) + sum(b))
    expect_true(any(grepl("2 counts layers", msgs)))
    # default identities are the sample names: not taken as clusters
    expect_null(r$clusters)
    expect_true(any(grepl("e.g. clusters = \"cell_type\"", msgs)))
    expect_equal(.lamb_groups("cell_type", r$clusters, colnames(r$counts), quiet, "seurat", r$meta), rep(c("x", "y"), 5))
})

test_that("Seurat: v3 Assay read; SCT/integrated default assay skipped for RNA; normalised counts refused", {
    skip_if_not_installed("SeuratObject")
    m <- rc_matrix(n = 8)
    old <- options(Seurat.object.assay.version = "v3"); on.exit(options(old))
    s <- suppressWarnings(SeuratObject::CreateSeuratObject(m))
    expect_s4_class(s[["RNA"]], "Assay")
    SeuratObject::Idents(s) <- factor(rep(c("A", "B"), 4))
    f <- tempfile(fileext = ".rds"); saveRDS(s, f)
    r <- .lamb_read_cells(f, "seurat", quiet)
    expect_equal(as.matrix(r$counts), as.matrix(m))
    expect_equal(unname(r$clusters), rep(c("A", "B"), 4))
    # an integrated assay with data only, set as default
    integ <- SeuratObject::CreateAssayObject(data = lognorm(m))
    suppressWarnings(s[["integrated"]] <- integ)
    SeuratObject::DefaultAssay(s) <- "integrated"
    saveRDS(s, f)
    msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
    r <- .lamb_read_cells(f, "seurat", say)
    expect_equal(as.matrix(r$counts), as.matrix(m))
    expect_true(any(grepl("assay 'RNA' \\(default assay 'integrated' skipped\\)", msgs)))
    # log-normalised values in the counts slot
    s2 <- suppressWarnings(SeuratObject::CreateSeuratObject(lognorm(m)))
    saveRDS(s2, f)
    expect_error(.lamb_read_cells(f, "seurat", quiet), "normalised values")
})

test_that("SingleCellExperiment: colData cluster column, Ensembl row names with rowData Symbol, no counts assay", {
    skip_if_not_installed("SingleCellExperiment")
    m <- rc_matrix(n = 8)
    ids <- paste0("ENSG9900000000", seq_len(nrow(m)))
    cm <- m; rownames(cm) <- ids
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = cm))
    SummarizedExperiment::rowData(sce)$Symbol <- rownames(m)
    SummarizedExperiment::colData(sce)$cluster <- rep(c("k1", "k2"), 4)
    f <- tempfile(fileext = ".rds"); saveRDS(sce, f)
    msgs <- character(); say <- function(...) msgs <<- c(msgs, paste0(...))
    r <- .lamb_read_cells(f, "sce", say)
    expect_equal(rownames(r$counts), rownames(m))
    expect_equal(unname(r$clusters), rep(c("k1", "k2"), 4))
    expect_true(any(grepl("rowData column Symbol", msgs)))
    expect_true(any(grepl("colData column 'cluster'", msgs)))
    lc <- SingleCellExperiment::SingleCellExperiment(list(logcounts = lognorm(m)))
    saveRDS(lc, f)
    expect_error(.lamb_read_cells(f, "sce", quiet), "no 'counts' assay \\(assays: logcounts\\)")
})

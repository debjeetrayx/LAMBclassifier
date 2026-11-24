# ---- Single-cell and spatial input readers ----
# Entry points, used by lamb_run():
#
#   .lamb_detect_cells(input)                 -> platform string, or NULL
#   .lamb_read_cells(input, platform, say)    -> list(counts, clusters, xy,
#                                                     pixel_size, platform, unit)
#
# Platforms: "10x_mtx", "10x_h5", "h5ad", "seurat", "sce", "xenium", "visium",
# "visium_hd", "cosmx", "merscope". NULL from detection means "not single-cell
# or spatial": the bulk path takes over.
#
# Returned list (same shapes for every platform):
#   counts     dgCMatrix, genes (HGNC symbols) x cells / spots / bins, raw counts
#   clusters   named character vector (cell id -> group) or NULL (none found)
#   xy         data.frame(cell_id, x, y) in MICRONS, or NULL (no positions)
#   pixel_size microns per pixel of the platform's image coordinates, or NULL
#   platform   the platform string
#   unit       "cell", "spot" (Visium) or "bin" (Visium HD)
#
# Microns per pixel: Visium 55 / spot_diameter_fullres (a spot is 55 um);
# Visium HD microns_per_pixel, else bin size / spot_diameter_fullres; CosMx
# 0.12028. Xenium and MERSCOPE coordinates are already microns.

.LAMB_PLATFORM_NAMES <- c(
    "10x_mtx" = "10x Cell Ranger", "10x_h5" = "10x Cell Ranger (.h5)", h5ad = "AnnData (.h5ad)",
    seurat = "Seurat", sce = "SingleCellExperiment", xenium = "Xenium", visium = "Visium",
    visium_hd = "Visium HD", cosmx = "CosMx", merscope = "MERSCOPE")

.LAMB_VISIUM_SPOT_UM <- 55
.LAMB_COSMX_UM_PER_PX <- 0.12028


# ---- detection ----

.lamb_detect_cells <- function(input) {
    if (!file.exists(input)) return(NULL)
    if (!dir.exists(input)) return(.lamb_detect_cells_file(input))
    d <- input
    has <- function(...) any(file.exists(file.path(d, c(...))))
    # Xenium first: its cell_feature_matrix/ would otherwise look like plain 10x.
    if (has("experiment.xenium", "cells.csv.gz", "cells.parquet", "cell_feature_matrix",
            "cell_feature_matrix.h5")) return("xenium")
    if (!is.null(.lamb_find_hd_bin(d, NULL))) return("visium_hd")
    if (!is.null(.lamb_find_visium_root(d))) return("visium")
    if (length(.lamb_find_cosmx(d))) return("cosmx")
    if (length(.lamb_find_merscope(d))) return("merscope")
    if (!is.null(.lamb_find_mtx(d))) return("10x_mtx")
    if (!is.null(.lamb_find_h5(d))) return("10x_h5")
    NULL
}

.lamb_detect_cells_file <- function(f) {
    if (grepl("\\.rds$", f, ignore.case = TRUE)) {
        obj <- readRDS(f)
        if (inherits(obj, "Seurat")) return("seurat")
        if (inherits(obj, "SingleCellExperiment")) return("sce")
        return(NULL)
    }
    if (grepl("\\.h5ad$", f, ignore.case = TRUE)) return("h5ad")
    if (grepl("\\.h5$", f, ignore.case = TRUE)) {
        # look inside when an HDF5 reader is present: an AnnData file saved as .h5
        top <- tryCatch(.lamb_h5_open(f)$ls(), error = function(e) NULL)
        if (!is.null(top) && all(c("X", "obs", "var") %in% top)) return("h5ad")
        return("10x_h5")
    }
    if (grepl("_exprMat_file\\.csv(\\.gz)?$", f) && length(.lamb_find_cosmx(dirname(f)))) return("cosmx")
    if (basename(f) %in% c("cell_by_gene.csv", "cell_by_gene.csv.gz") &&
        length(.lamb_find_merscope(dirname(f)))) return("merscope")
    NULL
}

.lamb_find_mtx <- function(dir) {
    for (d in c(dir, file.path(dir, "cell_feature_matrix"), file.path(dir, "filtered_feature_bc_matrix"),
                file.path(dir, "outs", "filtered_feature_bc_matrix"))) {
        if (any(file.exists(file.path(d, c("matrix.mtx.gz", "matrix.mtx"))))) return(d)
    }
    if (!is.null(.lamb_mtx_triplet(dir))) return(dir)
    NULL
}

# Matrix Market triplets under other names (e.g. count_matrix_sparse.mtx);
# accepted only when exactly one of each file is present.
.lamb_mtx_triplet <- function(dir) {
    if (!dir.exists(dir)) return(NULL)
    f <- list.files(dir, full.names = TRUE)
    f <- f[!startsWith(basename(f), "._")]
    one <- function(re) { x <- f[grepl(re, basename(f), ignore.case = TRUE)]; if (length(x) == 1) x else NULL }
    mtx <- one("\\.mtx(\\.gz)?$")
    feat <- one("(genes|features)\\.tsv(\\.gz)?$")
    bc <- one("barcodes\\.tsv(\\.gz)?$")
    if (is.null(mtx) || is.null(feat) || is.null(bc)) return(NULL)
    list(matrix = mtx, features = feat, barcodes = bc)
}

.lamb_find_h5 <- function(dir) {
    f <- file.path(c(dir, file.path(dir, "outs")),
                   rep(c("filtered_feature_bc_matrix.h5", "cell_feature_matrix.h5"), each = 2))
    f <- f[file.exists(f) & !dir.exists(f)]
    if (length(f)) f[1] else NULL
}

.lamb_find_visium_root <- function(dir) {
    for (r in c(dir, file.path(dir, "outs"))) {
        pos <- file.path(r, "spatial", c("tissue_positions.csv", "tissue_positions_list.csv"))
        if (any(file.exists(pos)) && (!is.null(.lamb_find_mtx(r)) || !is.null(.lamb_find_h5(r)))) return(r)
    }
    NULL
}

# bin = NULL: any square_* bin folder counts (detection); otherwise that bin.
.lamb_find_hd_bin <- function(dir, bin = "square_008um") {
    ok <- function(b) file.exists(file.path(b, "spatial", "tissue_positions.parquet")) &&
        (!is.null(.lamb_find_mtx(b)) || !is.null(.lamb_find_h5(b)))
    if (ok(dir)) return(dir)
    for (r in c(file.path(dir, "binned_outputs"), file.path(dir, "outs", "binned_outputs"))) {
        if (!dir.exists(r)) next
        bins <- if (is.null(bin)) list.files(r, pattern = "^square_") else bin
        for (b in file.path(r, bins)) if (ok(b)) return(b)
    }
    NULL
}

.lamb_find_cosmx <- function(dir) {
    e <- list.files(dir, pattern = "_exprMat_file\\.csv(\\.gz)?$", full.names = TRUE)
    m <- list.files(dir, pattern = "_metadata_file\\.csv(\\.gz)?$", full.names = TRUE)
    if (length(e) && length(m)) c(expr = e[1], meta = m[1]) else character()
}

.lamb_find_merscope <- function(dir) {
    e <- file.path(dir, c("cell_by_gene.csv", "cell_by_gene.csv.gz"))
    m <- file.path(dir, c("cell_metadata.csv", "cell_metadata.csv.gz"))
    if (any(file.exists(e)) && any(file.exists(m))) c(expr = e[file.exists(e)][1], meta = m[file.exists(m)][1]) else character()
}


# ---- reading: dispatcher ----

.lamb_read_cells <- function(input, platform, say = message, hd_bin = "square_008um") {
    platform <- match.arg(platform, names(.LAMB_PLATFORM_NAMES))
    out <- switch(platform,
                  seurat = , sce = .lamb_read_rds_cells(input, platform, say),
                  "10x_mtx" = , "10x_h5" = .lamb_read_10x(input, platform, say),
                  xenium = .lamb_read_10x(input, "xenium", say),
                  h5ad = .lamb_read_h5ad(input, say),
                  visium = .lamb_read_visium(input, say),
                  visium_hd = .lamb_read_visium_hd(input, hd_bin, say),
                  cosmx = .lamb_read_cosmx(input, say),
                  merscope = .lamb_read_merscope(input, say))
    out$platform <- platform
    if (is.null(out$unit)) out$unit <- "cell"
    out <- out[c("counts", "clusters", "xy", "pixel_size", "platform", "unit", "meta")]
    names(out) <- c("counts", "clusters", "xy", "pixel_size", "platform", "unit", "meta")
    out$meta <- .lamb_meta_columns(out$meta, colnames(out$counts))
    units <- c(cell = "cells", spot = "spots", bin = "bins")[[out$unit]]
    grp <- if (is.null(out$clusters)) "no groups" else
        paste0(length(unique(out$clusters)), " groups")
    pos <- if (is.null(out$xy)) "no positions" else "positions in microns"
    say("Loaded ", .LAMB_PLATFORM_NAMES[[platform]], ": ", format(ncol(out$counts), big.mark = ","), " ",
        units, ", ", format(nrow(out$counts), big.mark = ","), " genes, ", grp, ", ", pos, ".")
    if (is.null(out$clusters)) say(.lamb_no_clusters_text(platform, names(out$meta)))
    out
}

.lamb_no_clusters_text <- function(platform, meta_cols = NULL) {
    what <- c(cell = "cell", spot = "spot", bin = "bin")[[if (platform == "visium") "spot" else
        if (platform == "visium_hd") "bin" else "cell"]]
    txt <- paste0("This ", .LAMB_PLATFORM_NAMES[[platform]], " run has no ", what, " clusters. ")
    if (length(meta_cols)) {
        # Columns named like cell labels first.
        pref <- grepl("cell.?type|cluster|annot|label|level|ident|subtype|lineage", meta_cols, ignore.case = TRUE) &
            !grepl("ontology", meta_cols, ignore.case = TRUE)
        meta_cols <- c(meta_cols[pref], meta_cols[!pref])
        paste0(txt, "Name one of its cell annotation columns, e.g. clusters = \"", meta_cols[1], "\" (available: ",
               paste(utils::head(meta_cols, 12), collapse = ", "), if (length(meta_cols) > 12) ", ..." else "",
               "), or supply a file: clusters = 'my_clusters.csv' (two columns: ", what, " id, cluster).")
    } else {
        paste0(txt, "Supply them with clusters = 'my_clusters.csv' (two columns: ", what, " id, cluster).")
    }
}

# Usable group columns: 2-500 distinct values and at most one per ten cells
# (a repeating label, not an ID or a continuous measurement).
.lamb_meta_columns <- function(meta, cells) {
    if (is.null(meta) || !nrow(meta) || !ncol(meta)) return(NULL)
    meta <- as.data.frame(meta, stringsAsFactors = FALSE)
    if (is.null(rownames(meta)) || !all(cells %in% rownames(meta))) return(NULL)
    meta <- meta[cells, , drop = FALSE]
    keep <- vapply(meta, function(v) {
        (is.character(v) || is.factor(v) || is.logical(v) || (is.numeric(v) && all(v == round(v), na.rm = TRUE))) &&
            length(unique(stats::na.omit(v))) %in% 2:min(500, max(2, floor(length(cells) / 10)))
    }, TRUE)
    if (!any(keep)) return(NULL)
    meta[, keep, drop = FALSE]
}

# Cell metadata table next to a matrix; first column = cell barcode.
.lamb_read_meta_file <- function(root) {
    f <- list.files(root, pattern = "(^|_)(cell_)?metadata\\.(csv|tsv)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
    if (!length(f)) return(NULL)
    t <- utils::read.delim(f[1], sep = if (grepl("\\.csv", f[1])) "," else "\t", stringsAsFactors = FALSE,
                           check.names = FALSE)
    if (anyDuplicated(t[[1]])) return(NULL)
    rownames(t) <- as.character(t[[1]])
    t[-1]
}


# ---- Seurat / SingleCellExperiment ----

.lamb_read_rds_cells <- function(input, platform, say) {
    obj <- readRDS(input)
    if (inherits(obj, "Seurat")) {
        if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Install SeuratObject to read a Seurat object.", call. = FALSE)
        say("Seurat object: raw counts of the default assay; groups = active identities.")
        counts <- SeuratObject::LayerData(obj, layer = "counts")
        return(list(counts = methods::as(counts, "CsparseMatrix"),
                    clusters = stats::setNames(as.character(SeuratObject::Idents(obj)), colnames(obj)), xy = NULL, meta = obj[[]]))
    }
    if (inherits(obj, "SingleCellExperiment")) {
        if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("Install SingleCellExperiment.", call. = FALSE)
        say("SingleCellExperiment: counts assay; groups = colLabels().")
        lab <- SingleCellExperiment::colLabels(obj)
        return(list(counts = methods::as(SummarizedExperiment::assay(obj, "counts"), "CsparseMatrix"),
                    clusters = if (is.null(lab)) NULL else stats::setNames(as.character(lab), colnames(obj)), xy = NULL,
                    meta = as.data.frame(SummarizedExperiment::colData(obj))))
    }
    stop("The .rds file holds a ", class(obj)[1], "; expected a Seurat or SingleCellExperiment object, ",
         "or a genes x samples matrix for bulk.", call. = FALSE)
}


# ---- 10x matrices (mtx folder or .h5) ----

# rowsum(reorder = FALSE) for a sparse matrix, which base rowsum() refuses.
.lamb_sum_dup_rows <- function(m, sym) {
    u <- unique(sym)
    g <- Matrix::sparseMatrix(i = match(sym, u), j = seq_along(sym), x = 1, dims = c(length(u), length(sym)))
    out <- methods::as(g %*% m, "CsparseMatrix")
    dimnames(out) <- list(u, colnames(m))
    out
}

# Gene Expression features only; duplicated symbols summed.
.lamb_tidy_features <- function(m, sym, type, bc, say) {
    m <- methods::as(m, "CsparseMatrix")
    keep <- if (!is.null(type)) type == "Gene Expression" else rep(TRUE, length(sym))
    say("Matrix: ", sum(keep), " gene features kept (", sum(!keep), " control / codeword features dropped).")
    m <- m[keep, , drop = FALSE]
    sym <- sym[keep]
    if (anyDuplicated(sym)) m <- .lamb_sum_dup_rows(m, sym) else rownames(m) <- sym
    if (is.null(rownames(m))) rownames(m) <- sym
    colnames(m) <- bc
    m
}

.lamb_read_mtx_dir <- function(mdir, say) {
    rd <- function(f) { p <- file.path(mdir, f); if (file.exists(paste0(p, ".gz"))) gzfile(paste0(p, ".gz")) else p }
    if (any(file.exists(file.path(mdir, c("matrix.mtx.gz", "matrix.mtx"))))) {
        m <- Matrix::readMM(rd("matrix.mtx"))
        ff <- if (any(file.exists(file.path(mdir, c("features.tsv", "features.tsv.gz"))))) "features.tsv" else "genes.tsv"
        feat <- utils::read.delim(rd(ff), header = FALSE, stringsAsFactors = FALSE)
        bc <- utils::read.delim(rd("barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)[[1]]
    } else {
        tr <- .lamb_mtx_triplet(mdir)
        say("Matrix Market files: ", basename(tr$matrix), ", ", basename(tr$features), ", ", basename(tr$barcodes), ".")
        m <- Matrix::readMM(tr$matrix)
        feat <- utils::read.delim(tr$features, header = FALSE, stringsAsFactors = FALSE)
        bc <- utils::read.delim(tr$barcodes, header = FALSE, stringsAsFactors = FALSE)[[1]]
    }
    sym <- feat[[min(2, ncol(feat))]]   # 10x: id, symbol, type; other exports: symbol only
    .lamb_tidy_features(m, sym, if (ncol(feat) >= 3) feat[[3]] else NULL, bc, say)
}

# 10x HDF5 matrix: Cell Ranger v3+ (/matrix/...) or v2 (/<genome>/...).
.lamb_read_10x_h5_file <- function(f, say) {
    h <- .lamb_h5_open(f)
    paths <- h$ls()
    base <- if ("matrix" %in% paths) "matrix" else {
        g <- sub("/data$", "", grep("^[^/]+/data$", paths, value = TRUE))
        if (!length(g)) stop("The file ", f, " is not a 10x matrix .h5 (no /matrix group).", call. = FALSE)
        g[1]
    }
    p <- function(x) paste0(base, "/", x)
    shape <- as.numeric(h$read(p("shape")))
    m <- Matrix::sparseMatrix(i = as.numeric(h$read(p("indices"))) + 1, p = as.numeric(h$read(p("indptr"))),
                              x = as.numeric(h$read(p("data"))), dims = shape)
    if (p("features/name") %in% paths) {
        sym <- as.character(h$read(p("features/name")))
        type <- if (p("features/feature_type") %in% paths) as.character(h$read(p("features/feature_type"))) else NULL
    } else {
        sym <- as.character(h$read(p("gene_names")))
        type <- NULL
    }
    .lamb_tidy_features(m, sym, type, as.character(h$read(p("barcodes"))), say)
}

# Matrix from a folder (mtx preferred, else .h5) or from an .h5 file.
# Returns list(m, root): root is the folder the instrument's analysis/ sits next to.
.lamb_read_10x_matrix <- function(input, say) {
    if (!dir.exists(input)) return(list(m = .lamb_read_10x_h5_file(input, say), root = dirname(input)))
    mdir <- .lamb_find_mtx(input)
    if (!is.null(mdir)) {
        root <- if (basename(mdir) %in% c("cell_feature_matrix", "filtered_feature_bc_matrix")) dirname(mdir) else mdir
        return(list(m = .lamb_read_mtx_dir(mdir, say), root = root))
    }
    h5 <- .lamb_find_h5(input)
    if (is.null(h5)) stop("No matrix.mtx(.gz) found under ", input, call. = FALSE)
    say("Matrix file: ", basename(h5), ".")
    list(m = .lamb_read_10x_h5_file(h5, say), root = dirname(h5))
}

# graphclust/ is the older Cell Ranger / Space Ranger name, tried second.
.lamb_find_clusters <- function(root, say) {
    for (sub in c("gene_expression_graphclust", "graphclust")) {
        f <- Filter(file.exists, file.path(c(root, file.path(root, "outs"), dirname(root)),
                                           "analysis", "clustering", sub, "clusters.csv"))
        if (length(f)) {
            t <- utils::read.csv(f[1], stringsAsFactors = FALSE)
            say("Groups: the instrument's graph-based clusters (", f[1], ").")
            return(stats::setNames(paste0("cluster ", t[[2]]), t[[1]]))
        }
    }
    NULL
}

# 10x mtx / 10x h5 / Xenium.
.lamb_read_10x <- function(input, platform, say) {
    mat <- .lamb_read_10x_matrix(input, say)
    root <- mat$root
    cl <- .lamb_find_clusters(root, say)
    xy <- NULL
    if (platform == "xenium") {
        cf <- file.path(root, "cells.csv.gz")
        pq <- file.path(root, "cells.parquet")
        if (file.exists(cf)) {
            ct <- utils::read.csv(gzfile(cf), stringsAsFactors = FALSE)
            xy <- data.frame(cell_id = ct$cell_id, x = ct$x_centroid, y = ct$y_centroid, stringsAsFactors = FALSE)
        } else if (file.exists(pq) && requireNamespace("arrow", quietly = TRUE)) {
            ct <- as.data.frame(arrow::read_parquet(pq))
            xy <- data.frame(cell_id = as.character(ct$cell_id), x = ct$x_centroid, y = ct$y_centroid,
                             stringsAsFactors = FALSE)
        } else {
            say("cells.csv.gz not found: no tissue map.")
        }
        if (!is.null(xy)) say("Positions: cell centroids from Xenium, already in microns.")
    }
    px <- 0.2125
    ex <- file.path(root, "experiment.xenium")
    if (file.exists(ex)) {
        j <- paste(readLines(ex, warn = FALSE), collapse = " ")
        v <- regmatches(j, regexpr('"pixel_size"[^0-9]*[0-9.]+', j))
        if (length(v)) px <- as.numeric(sub('.*[^0-9.]', "", v))
    }
    list(counts = mat$m, clusters = cl, xy = xy, pixel_size = px, unit = "cell", meta = .lamb_read_meta_file(root))
}


# ---- Visium and Visium HD (Space Ranger) ----

.lamb_json_number <- function(file, key) {
    if (!file.exists(file)) return(NA_real_)
    j <- paste(readLines(file, warn = FALSE), collapse = " ")
    v <- regmatches(j, regexpr(paste0('"', key, '"\\s*:\\s*[-0-9.eE+]+'), j))
    if (length(v)) as.numeric(sub('.*:\\s*', "", v)) else NA_real_
}

.lamb_spatial_xy <- function(pos, um_per_px) {
    data.frame(cell_id = as.character(pos[[1]]),
               x = as.numeric(pos[[6]]) * um_per_px, y = as.numeric(pos[[5]]) * um_per_px,
               stringsAsFactors = FALSE)
}

# Positions table columns (both layouts): barcode, in_tissue, array_row,
# array_col, pxl_row_in_fullres, pxl_col_in_fullres.
.lamb_keep_in_tissue <- function(m, pos, unit, say) {
    pos <- pos[as.integer(pos[[2]]) == 1, , drop = FALSE]
    keep <- colnames(m) %in% as.character(pos[[1]])
    if (any(!keep)) say(sum(!keep), " ", unit, "s outside the tissue (in_tissue = 0) left out.")
    list(m = m[, keep, drop = FALSE], pos = pos)
}

.lamb_read_visium <- function(input, say) {
    root <- .lamb_find_visium_root(input)
    if (is.null(root)) stop("No Visium spatial/tissue_positions file with a filtered matrix found under ", input, call. = FALSE)
    m <- .lamb_read_10x_matrix(root, say)$m
    sp <- file.path(root, "spatial")
    new <- file.path(sp, "tissue_positions.csv")
    pos <- if (file.exists(new)) utils::read.csv(new, stringsAsFactors = FALSE) else
        utils::read.csv(file.path(sp, "tissue_positions_list.csv"), header = FALSE, stringsAsFactors = FALSE)
    k <- .lamb_keep_in_tissue(m, pos, "spot", say)
    d <- .lamb_json_number(file.path(sp, "scalefactors_json.json"), "spot_diameter_fullres")
    xy <- NULL
    um <- NULL
    if (is.na(d) || d <= 0) {
        say("scalefactors_json.json (spot_diameter_fullres) not found: spot positions cannot be put in microns, no tissue map.")
    } else {
        um <- .LAMB_VISIUM_SPOT_UM / d
        xy <- .lamb_spatial_xy(k$pos, um)
        say("Positions: spot centres converted to microns (55 um spot = ", signif(d, 5),
            " full-resolution pixels, so ", signif(um, 4), " um per pixel).")
    }
    list(counts = k$m, clusters = .lamb_find_clusters(root, say), xy = xy, pixel_size = um, unit = "spot")
}

.lamb_read_visium_hd <- function(input, bin, say) {
    bdir <- .lamb_find_hd_bin(input, bin)
    if (is.null(bdir)) {
        have <- list.files(c(file.path(input, "binned_outputs"), file.path(input, "outs", "binned_outputs")), pattern = "^square_")
        stop("No Visium HD bin folder '", bin, "' with a filtered matrix and spatial/tissue_positions.parquet under ",
             input, if (length(have)) paste0(". Bins present: ", paste(have, collapse = ", ")), call. = FALSE)
    }
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("Reading Visium HD positions (tissue_positions.parquet) needs the arrow package. ",
             "Install it with install.packages('arrow') and run again.", call. = FALSE)
    }
    say("Visium HD: using the ", basename(bdir), " bins.")
    m <- .lamb_read_10x_matrix(bdir, say)$m
    pos <- as.data.frame(arrow::read_parquet(file.path(bdir, "spatial", "tissue_positions.parquet")))
    k <- .lamb_keep_in_tissue(m, pos, "bin", say)
    sf <- file.path(bdir, "spatial", "scalefactors_json.json")
    um <- .lamb_json_number(sf, "microns_per_pixel")
    how <- "microns_per_pixel in scalefactors_json.json"
    if (is.na(um)) {
        size <- .lamb_json_number(sf, "bin_size_um")
        if (is.na(size)) size <- suppressWarnings(as.numeric(sub("^square_0*([0-9]+)um$", "\\1", basename(bdir))))
        um <- size / .lamb_json_number(sf, "spot_diameter_fullres")
        how <- paste0(size, " um bin / spot_diameter_fullres")
    }
    xy <- NULL
    if (is.na(um) || um <= 0) {
        um <- NULL
        say("scalefactors_json.json not usable: bin positions cannot be put in microns, no tissue map.")
    } else {
        xy <- .lamb_spatial_xy(k$pos, um)
        say("Positions: bin centres converted to microns (", signif(um, 4), " um per pixel, from ", how, ").")
    }
    list(counts = k$m, clusters = .lamb_find_clusters(bdir, say), xy = xy, pixel_size = um, unit = "bin")
}


# ---- CosMx and MERSCOPE (flat files) ----

# CSV reader: data.table when present (fast on large flat files), else utils.
# id_cols are read as text (MERSCOPE cell ids exceed double precision).
.lamb_read_flat <- function(path, id_cols = 1) {
    gz <- grepl("\\.gz$", path)
    # fread() needs R.utils for .gz files; base R reads them without it.
    if (requireNamespace("data.table", quietly = TRUE) && (!gz || requireNamespace("R.utils", quietly = TRUE))) {
        hdr <- names(data.table::fread(path, nrows = 0))
        cc <- if (length(id_cols)) stats::setNames(rep("character", length(id_cols)), hdr[id_cols]) else NULL
        return(as.data.frame(data.table::fread(path, colClasses = cc, check.names = FALSE, showProgress = FALSE)))
    }
    con <- if (grepl("\\.gz$", path)) gzfile(path) else path
    hdr <- names(utils::read.csv(con, nrows = 1, check.names = FALSE))
    cc <- rep(NA_character_, length(hdr)); cc[id_cols] <- "character"
    con <- if (grepl("\\.gz$", path)) gzfile(path) else path
    utils::read.csv(con, colClasses = cc, check.names = FALSE, stringsAsFactors = FALSE)
}

.lamb_flat_counts <- function(tab, gene_cols, ids) {
    x <- as.matrix(tab[, gene_cols, drop = FALSE])
    storage.mode(x) <- "double"
    m <- .lamb_dgc(Matrix::t(.lamb_dgc(x)))
    dimnames(m) <- list(gene_cols, ids)
    m
}

.lamb_read_cosmx <- function(input, say) {
    f <- .lamb_find_cosmx(if (dir.exists(input)) input else dirname(input))
    if (!length(f)) stop("No CosMx *_exprMat_file.csv and *_metadata_file.csv found in ", input, call. = FALSE)
    ex <- .lamb_read_flat(f[["expr"]], id_cols = integer())
    if (!all(c("fov", "cell_ID") %in% names(ex))) stop(basename(f[["expr"]]), " needs the columns fov and cell_ID.", call. = FALSE)
    bg <- ex$cell_ID == 0
    if (any(bg)) say(sum(bg), " rows with cell_ID 0 (transcripts outside any cell) left out.")
    ex <- ex[!bg, , drop = FALSE]
    other <- setdiff(names(ex), c("fov", "cell_ID", "cell"))
    ctrl <- grepl("^(NegPrb|Negative|SystemControl|FalseCode)", other)
    genes <- other[!ctrl]
    say("CosMx: ", length(genes), " genes kept (", sum(ctrl), " negative-probe / system-control columns dropped).")
    ids <- paste(ex$fov, ex$cell_ID, sep = "_")
    m <- .lamb_flat_counts(ex, genes, ids)
    if (anyDuplicated(genes)) m <- .lamb_sum_dup_rows(m, genes)
    say("CosMx cell id = <fov>_<cell_ID> (for example ", ids[1], ").")

    md <- .lamb_read_flat(f[["meta"]], id_cols = integer())
    xy <- NULL
    if (all(c("fov", "cell_ID", "CenterX_global_px", "CenterY_global_px") %in% names(md))) {
        xy <- data.frame(cell_id = paste(md$fov, md$cell_ID, sep = "_"),
                         x = md$CenterX_global_px * .LAMB_COSMX_UM_PER_PX,
                         y = md$CenterY_global_px * .LAMB_COSMX_UM_PER_PX, stringsAsFactors = FALSE)
        xy <- xy[xy$cell_id %in% ids, , drop = FALSE]
        say("Positions: CenterX/Y_global_px converted to microns (", .LAMB_COSMX_UM_PER_PX, " um per CosMx pixel).")
    } else {
        say(basename(f[["meta"]]), " has no CenterX_global_px / CenterY_global_px: no tissue map.")
    }
    list(counts = m, clusters = NULL, xy = xy, pixel_size = .LAMB_COSMX_UM_PER_PX, unit = "cell")
}

.lamb_read_merscope <- function(input, say) {
    f <- .lamb_find_merscope(if (dir.exists(input)) input else dirname(input))
    if (!length(f)) stop("No MERSCOPE cell_by_gene.csv and cell_metadata.csv found in ", input, call. = FALSE)
    ex <- .lamb_read_flat(f[["expr"]], id_cols = 1)
    ids <- as.character(ex[[1]])
    other <- names(ex)[-1]
    blank <- grepl("^Blank-", other)
    genes <- other[!blank]
    say("MERSCOPE: ", length(genes), " genes kept (", sum(blank), " Blank-* control columns dropped).")
    m <- .lamb_flat_counts(ex, genes, ids)
    if (anyDuplicated(genes)) m <- .lamb_sum_dup_rows(m, genes)

    md <- .lamb_read_flat(f[["meta"]], id_cols = 1)
    xy <- NULL
    if (all(c("center_x", "center_y") %in% names(md))) {
        xy <- data.frame(cell_id = as.character(md[[1]]), x = as.numeric(md$center_x), y = as.numeric(md$center_y),
                         stringsAsFactors = FALSE)
        xy <- xy[xy$cell_id %in% ids, , drop = FALSE]
        say("Positions: center_x / center_y from cell_metadata.csv, already in microns.")
    } else {
        say("cell_metadata.csv has no center_x / center_y: no tissue map.")
    }
    list(counts = m, clusters = NULL, xy = xy, pixel_size = NULL, unit = "cell")
}


# ---- HDF5 access (rhdf5, else hdf5r) ----

.lamb_h5_backend <- function() {
    opt <- getOption("LAMBclassifier.h5_backend")
    for (b in c(opt, "rhdf5", "hdf5r")) if (requireNamespace(b, quietly = TRUE)) return(b)
    stop("Reading .h5 / .h5ad files needs an HDF5 reader. Install one with ",
         "install.packages('BiocManager'); BiocManager::install('rhdf5')  ",
         "(or install.packages('hdf5r')) and run again.", call. = FALSE)
}

# Small uniform interface: ls() all object paths (no leading '/'), read(path),
# attrs(path) as a named list. 1-D results come back as plain vectors.
.lamb_h5_open <- function(f) {
    flat <- function(x) {
        if (inherits(x, "integer64")) x <- as.numeric(x)
        if (is.factor(x)) x <- as.character(x)
        if (is.data.frame(x) || length(dim(x)) > 1) x else as.vector(x)
    }
    if (.lamb_h5_backend() == "rhdf5") {
        list(
            ls = function() {
                l <- rhdf5::h5ls(f, recursive = TRUE, all = FALSE)
                sub("^/+", "", file.path(sub("^/$", "", l$group), l$name))
            },
            read = function(p) flat(rhdf5::h5read(f, p, bit64conversion = "double")),
            attrs = function(p) lapply(rhdf5::h5readAttributes(f, p), flat))
    } else {
        with_h <- function(fun) {
            h <- hdf5r::H5File$new(f, mode = "r")
            on.exit(h$close_all())
            fun(h)
        }
        list(
            ls = function() with_h(function(h) h$ls(recursive = TRUE)$name),
            read = function(p) with_h(function(h) flat(h[[p]]$read())),
            attrs = function(p) with_h(function(h) {
                o <- h[[p]]
                n <- hdf5r::h5attr_names(o)
                stats::setNames(lapply(n, function(a) flat(hdf5r::h5attr(o, a))), n)
            }))
    }
}


# ---- AnnData (.h5ad) ----

.lamb_is_whole <- function(v) {
    v <- v[!is.na(v)]
    if (length(v) > 1e6) v <- v[seq_len(1e6)]
    length(v) > 0 && all(v >= 0) && all(abs(v - round(v)) < 1e-6)
}

# A matrix element of an h5ad (dense dataset or csr/csc group) -> dgCMatrix genes x cells.
.lamb_h5ad_matrix <- function(h, paths, p) {
    if (paste0(p, "/data") %in% paths) {
        a <- h$attrs(p)
        fmt <- a[["encoding-type"]]
        if (is.null(fmt)) fmt <- paste0(a[["h5sparse_format"]], "_matrix")
        shape <- as.numeric(if (!is.null(a[["shape"]])) a[["shape"]] else a[["h5sparse_shape"]])
        x <- as.numeric(h$read(paste0(p, "/data")))
        i <- as.numeric(h$read(paste0(p, "/indices"))) + 1
        ptr <- as.numeric(h$read(paste0(p, "/indptr")))
        # CSR of cells x genes is CSC of genes x cells
        if (fmt == "csr_matrix") return(Matrix::sparseMatrix(i = i, p = ptr, x = x, dims = rev(shape)))
        if (fmt == "csc_matrix") return(Matrix::t(Matrix::sparseMatrix(i = i, p = ptr, x = x, dims = shape)))
        stop("Unknown sparse matrix format in ", p, ": ", fmt, call. = FALSE)
    }
    x <- h$read(p)   # HDF5 (cells, genes) in C order reads in R as genes x cells
    storage.mode(x) <- "double"
    .lamb_dgc(x)
}

# A data frame element (obs / var): its index and named columns.
.lamb_h5ad_index <- function(h, paths, grp) {
    if (!grp %in% paths) return(NULL)
    if (!any(startsWith(paths, paste0(grp, "/")))) {   # anndata < 0.7: compound table
        tab <- h$read(grp)
        return(as.character(tab[[if ("index" %in% names(tab)) "index" else 1]]))
    }
    idx <- h$attrs(grp)[["_index"]]
    if (is.null(idx)) idx <- "_index"
    as.character(h$read(paste0(grp, "/", idx)))
}

.lamb_h5ad_column <- function(h, paths, grp, col) {
    p <- paste0(grp, "/", col)
    if (!grp %in% paths) return(NULL)
    if (!any(startsWith(paths, paste0(grp, "/")))) {   # anndata < 0.7: compound table
        tab <- h$read(grp)
        return(if (col %in% names(tab)) as.character(tab[[col]]) else NULL)
    }
    if (!p %in% paths) return(NULL)
    if (paste0(p, "/codes") %in% paths) {   # anndata >= 0.8 categorical
        lev <- as.character(h$read(paste0(p, "/categories")))
        codes <- as.integer(h$read(paste0(p, "/codes")))
    } else {
        v <- h$read(p)
        cp <- paste0(grp, "/__categories/", col)   # anndata 0.7 categorical
        if (!cp %in% paths) return(as.character(v))
        lev <- as.character(h$read(cp))
        codes <- as.integer(v)
    }
    out <- rep(NA_character_, length(codes))
    ok <- codes >= 0
    out[ok] <- lev[codes[ok] + 1]
    out
}

.lamb_ensembl_share <- function(x) mean(grepl("^ENS[A-Z]*G[0-9]+", x))

.lamb_read_h5ad <- function(input, say) {
    h <- .lamb_h5_open(input)
    paths <- h$ls()
    cand <- c("layers/counts", "raw/X", "X")
    cand <- cand[cand %in% paths]
    if (!length(cand)) stop("The .h5ad file has no expression matrix (X).", call. = FALSE)
    m <- NULL
    for (p in cand) {
        x <- .lamb_h5ad_matrix(h, paths, p)
        if (.lamb_is_whole(x@x)) { m <- x; src <- p; break }
    }
    if (is.null(m)) {
        stop("This .h5ad holds only normalised (non-whole-number) values in ", paste(cand, collapse = ", "),
             ". LAMB adds up raw counts per cluster, so it needs the raw counts. ",
             "Save them in the file (adata.layers['counts'] = raw counts, or adata.raw) and run again.", call. = FALSE)
    }
    say("AnnData: raw counts from ", src, ".")
    var <- if (src == "raw/X") "raw/var" else "var"
    cells <- .lamb_h5ad_index(h, paths, "obs")

    # gene symbols: var index or the first symbol-like column
    opts <- c("feature_name", "gene_symbols", "gene_symbol", "gene_name", "gene_names", "symbol", "gene_ids")
    sym <- .lamb_h5ad_index(h, paths, var)
    from <- paste0(var, " index")
    if (is.null(sym) || .lamb_ensembl_share(sym) > 0.5) {
        for (o in opts) {
            v <- .lamb_h5ad_column(h, paths, var, o)
            if (!is.null(v) && .lamb_ensembl_share(v) <= 0.5) { sym <- v; from <- paste0(var, "/", o); break }
        }
    }
    if (length(sym) != nrow(m)) stop("Gene names in ", var, " do not match the matrix size.", call. = FALSE)
    if (.lamb_ensembl_share(sym) > 0.5) {
        gl <- get_gene_length()
        s2 <- gl$SYMBOL[match(sub("\\.[0-9]+$", "", sym), gl$gene_id)]
        keep <- !is.na(s2) & nzchar(s2)
        say("Gene IDs are Ensembl: mapped ", sum(keep), " of ", length(sym), " to symbols.")
        m <- m[keep, , drop = FALSE]
        sym <- s2[keep]
    } else {
        say("Gene symbols: ", from, ".")
    }
    sym <- as.character(sym)
    if (anyDuplicated(sym)) {
        say(sum(duplicated(sym)), " duplicated gene names: counts summed.")
        m <- .lamb_sum_dup_rows(m, sym)
    } else {
        rownames(m) <- sym
    }
    colnames(m) <- cells

    cl <- NULL
    for (k in c("leiden", "louvain", "cluster", "clusters", "seurat_clusters")) {
        v <- .lamb_h5ad_column(h, paths, "obs", k)
        if (!is.null(v)) {
            cl <- stats::setNames(v, cells)
            say("Groups: the obs column '", k, "' (", length(unique(stats::na.omit(v))), " groups).")
            break
        }
    }
    # Every obs column, so a cell-type column can be named with clusters = "...".
    obs_cols <- unique(sub("/.*", "", sub("^obs/", "", grep("^obs/", paths, value = TRUE))))
    obs_cols <- setdiff(obs_cols, c("_index", "__categories", h$attrs("obs")[["_index"]]))
    meta <- NULL
    if (length(obs_cols)) {
        cols <- lapply(obs_cols, function(k) tryCatch(.lamb_h5ad_column(h, paths, "obs", k), error = function(e) NULL))
        ok <- vapply(cols, function(v) length(v) == length(cells), TRUE)
        if (any(ok)) {
            meta <- as.data.frame(stats::setNames(cols[ok], obs_cols[ok]), stringsAsFactors = FALSE, check.names = FALSE)
            rownames(meta) <- make.unique(cells)
        }
    }
    list(counts = .lamb_dgc(m), clusters = cl, xy = NULL, pixel_size = NULL, unit = "cell", meta = meta)
}

# Any matrix -> dgCMatrix (never a symmetric / diagonal / pattern class).
.lamb_dgc <- function(x) {
    methods::as(methods::as(methods::as(x, "CsparseMatrix"), "generalMatrix"), "dMatrix")
}

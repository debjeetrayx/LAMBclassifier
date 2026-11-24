# Frozen earlier single-cell loader: regression reference for R/read_cells.R. Do not edit.

.lamb_load_cells <- function(input, data_type, say) {
    if (!dir.exists(input)) {
        obj <- readRDS(input)
        if (inherits(obj, "Seurat")) {
            if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Install SeuratObject to read a Seurat object.", call. = FALSE)
            say("Seurat object: raw counts of the default assay; groups = active identities.")
            counts <- SeuratObject::LayerData(obj, layer = "counts")
            return(list(counts = methods::as(counts, "CsparseMatrix"),
                        clusters = stats::setNames(as.character(SeuratObject::Idents(obj)), colnames(obj)), xy = NULL))
        }
        if (inherits(obj, "SingleCellExperiment")) {
            if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) stop("Install SingleCellExperiment.", call. = FALSE)
            say("SingleCellExperiment: counts assay; groups = colLabels().")
            lab <- SingleCellExperiment::colLabels(obj)
            return(list(counts = methods::as(SummarizedExperiment::assay(obj, "counts"), "CsparseMatrix"),
                        clusters = if (is.null(lab)) NULL else stats::setNames(as.character(lab), colnames(obj)), xy = NULL))
        }
        stop("The .rds file holds a ", class(obj)[1], "; expected a Seurat or SingleCellExperiment object, ",
             "or a genes x samples matrix for bulk.", call. = FALSE)
    }
    mdir <- .lamb_find_mtx(input)
    if (is.null(mdir)) stop("No matrix.mtx(.gz) found under ", input, call. = FALSE)
    rd <- function(f) { p <- file.path(mdir, f); if (file.exists(paste0(p, ".gz"))) gzfile(paste0(p, ".gz")) else p }
    m <- Matrix::readMM(rd("matrix.mtx"))
    feat <- utils::read.delim(rd("features.tsv"), header = FALSE, stringsAsFactors = FALSE)
    bc <- utils::read.delim(rd("barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)[[1]]
    m <- methods::as(m, "CsparseMatrix")
    keep <- if (ncol(feat) >= 3) feat[[3]] == "Gene Expression" else rep(TRUE, nrow(feat))
    say("Matrix: ", sum(keep), " gene features kept (", sum(!keep), " control / codeword features dropped).")
    m <- m[keep, , drop = FALSE]
    sym <- feat[[2]][keep]
    if (anyDuplicated(sym)) m <- methods::as(rowsum(m, sym, reorder = FALSE), "CsparseMatrix") else rownames(m) <- sym
    if (is.null(rownames(m))) rownames(m) <- sym
    colnames(m) <- bc

    root <- if (basename(mdir) %in% c("cell_feature_matrix", "filtered_feature_bc_matrix")) dirname(mdir) else mdir
    cl_file <- Filter(file.exists, file.path(c(root, file.path(root, "outs"), dirname(root)),
                                             "analysis", "clustering", "gene_expression_graphclust", "clusters.csv"))
    cl <- NULL
    if (length(cl_file)) {
        t <- utils::read.csv(cl_file[1], stringsAsFactors = FALSE)
        cl <- stats::setNames(paste0("cluster ", t[[2]]), t[[1]])
        say("Groups: the instrument's graph-based clusters (", cl_file[1], ").")
    }
    xy <- NULL
    if (data_type == "xenium") {
        cf <- file.path(root, "cells.csv.gz")
        if (file.exists(cf)) {
            ct <- utils::read.csv(gzfile(cf), stringsAsFactors = FALSE)
            xy <- data.frame(cell_id = ct$cell_id, x = ct$x_centroid, y = ct$y_centroid, stringsAsFactors = FALSE)
        } else {
            say("cells.csv.gz not found: no tissue map.")
        }
    }
    px <- 0.2125
    ex <- file.path(root, "experiment.xenium")
    if (file.exists(ex)) {
        j <- paste(readLines(ex, warn = FALSE), collapse = " ")
        v <- regmatches(j, regexpr('"pixel_size"[^0-9]*[0-9.]+', j))
        if (length(v)) px <- as.numeric(sub('.*[^0-9.]', "", v))
    }
    list(counts = m, clusters = cl, xy = xy, pixel_size = px)
}


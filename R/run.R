#' Run LAMBclassifier on a file or folder in one command
#'
#' @description
#' Point it at your data and it does the rest: it recognises the data type,
#' prepares the expression values, classifies, and writes a small fixed set of
#' result files.
#'
#' The simplest way: call \code{lamb_run()} with nothing in the brackets. It asks
#' whether your data is a file or a folder, opens your computer's own file
#' picker, runs, and opens the report. For a desktop icon that does the same
#' without opening R, run \code{\link{lamb_install_launcher}()} once.
#'
#' From a terminal:
#' \preformatted{Rscript -e 'LAMBclassifier::lamb_run("path/to/data")'}
#'
#' @details
#' \strong{What it accepts}
#' \describe{
#'   \item{Bulk, one patient or many}{A table (\code{.csv}, \code{.tsv},
#'     \code{.txt}, gzipped or not, any of \code{,} \code{;} tab separators,
#'     decimal commas, \code{.xlsx}) with genes in rows or in columns (HGNC
#'     symbols or Ensembl IDs; extra annotation columns are ignored); a single
#'     quantifier file (Salmon \code{quant.sf}, kallisto \code{abundance.tsv},
#'     RSEM \code{genes.results}, featureCounts, STAR \code{ReadsPerGene.out.tab},
#'     HTSeq-count); a folder of such files, one per patient; or an \code{.rds}
#'     matrix. Raw counts, TPM / FPKM and log2 values are recognised and
#'     converted automatically.}
#'   \item{Single-cell}{A 10x Cell Ranger folder or \code{.h5} file, an
#'     \code{.h5ad} (AnnData) file with raw counts, or an \code{.rds} file holding
#'     a Seurat or SingleCellExperiment object.}
#'   \item{Spatial}{The output folder of Xenium, Visium (Space Ranger), Visium HD,
#'     CosMx (flat files) or MERSCOPE. CosMx and MERSCOPE carry no clusters:
#'     supply them with \code{clusters}.}
#' }
#'
#' \strong{Single-cell and spatial: groups of cells.} LAMB scores groups of cells
#' (pseudobulk), never single cells. By default the groups are the graph-based
#' clusters the instrument already wrote
#' (\code{analysis/clustering/gene_expression_graphclust/clusters.csv}); for a
#' Seurat object, its active identities; for a SingleCellExperiment, its
#' \code{colLabels()}. Supply \code{clusters} to use your own.
#'
#' \strong{Cell-type gate.} The classes describe breast tumour epithelium. Within
#' each cluster's pseudobulk, gene modules for epithelium, immune cells,
#' vessels, adipocytes and fibroblasts are compared (mean within-cluster
#' percentile rank of each module's 3 highest-ranked genes; 3 measured genes needed). A cluster
#' is classified unless another module ranks above the epithelial one; such
#' clusters are named by what dominates (immune, vascular, adipocyte,
#' fibroblast-like), left out of the summary, and their class is not shown.
#' Metaplastic tumour cells can lose epithelial genes (chondroid, neural-like):
#' name those clusters with \code{tumour_clusters} to classify them.
#'
#' \strong{Output} (in \code{out_dir}):
#' \describe{
#'   \item{\code{LAMB_calls.csv}}{One row per sample (bulk) or per cluster
#'     (single-cell, spatial): the class, the HER2 amplicon read-out ("high",
#'     "not high" or "not assessable"; it never changes the class), how sure
#'     the call is, and a plain-language note.}
#'   \item{\code{LAMB_report.pdf}}{Bulk: one page per sample (see
#'     \code{\link{lamb_report}}) showing the class, what it means, how sure the
#'     call is, gene coverage and the decision path; a cohort overview first
#'     when there is more than one sample.}
#'   \item{\code{LAMB_genes.pdf}}{Bulk: one page per sample; single-cell and
#'     spatial: one page per classified cluster. Every gene LAMB scores as its
#'     rank within the sample (100 = highest), grouped by signature (see
#'     \code{\link{lamb_gene_report}}).}
#'   \item{\code{LAMB_slide.csv}}{Single-cell and spatial: the slide-level class
#'     (the class holding most classified epithelial cells; "mixed" when no class
#'     reaches half) and the percentage of epithelial cells in each class.}
#'   \item{\code{LAMB_summary.pdf}}{Single-cell and spatial: the slide-level
#'     class and the percentages; for spatial data a tissue map coloured by
#'     class, drawn on the H&E for Xenium when \code{he_image} and
#'     \code{he_alignment} are given.}
#'   \item{\code{LAMB_clusters.pdf}}{Single-cell and spatial: why each
#'     epithelial cluster got its class.}
#'   \item{\code{LAMB_xenium_explorer.csv}}{Xenium only: load in Xenium Explorer
#'     (Cells > Cell groups) to see the classes on the tissue. Every cell is in
#'     exactly one group.}
#'   \item{\code{technical/}}{Scores, signature coverage, provenance hashes and a
#'     run log, for whoever checks or reproduces the run.}
#' }
#'
#' @param input Path to a file or folder (see "What it accepts"). Leave it out
#'   to choose with a file picker.
#' @param out_dir Folder for the results. Default: \code{LAMB_results_<name>}
#'   next to the input.
#' @param clusters Optional. Your own cell groups: the name of a cell
#'   annotation column in the data (an \code{.h5ad} obs column, Seurat
#'   \code{meta.data}, SingleCellExperiment \code{colData}, or a
#'   \code{metadata.csv} next to a matrix), for example \code{"cell_type"}; a
#'   CSV file with two columns (barcode, group); or a named vector (names = cell
#'   barcodes).
#' @param data_type \code{"auto"} (default), \code{"bulk"}, \code{"single_cell"}
#'   or \code{"spatial"}. Leave on \code{"auto"}: the platform is recognised from
#'   the files present.
#' @param tumour_clusters Optional. Names of clusters to classify as tumour
#'   whatever their markers, for metaplastic tumour cells that have lost
#'   epithelial genes (single-cell and spatial only).
#' @param min_cells Clusters with fewer cells are reported but not classified
#'   (default 50).
#' @param values Bulk only. What the numbers are: \code{"auto"} (default)
#'   decides from the values and says so in the run log; \code{"counts"},
#'   \code{"tpm"} or \code{"log2"} declare it. Declare it when whole numbers are
#'   not raw counts (for example rounded TPM or intensities): automatic
#'   detection would treat them as counts and normalise them by gene length.
#'   See \code{vignette("input-rules")}.
#' @param run_pam50 Also run PAM50 on bulk data. \code{"auto"} (default) runs it
#'   only with 20 or more samples, because PAM50 centres each gene across the
#'   samples run together; \code{TRUE} / \code{FALSE} force it. Ignored for
#'   single-cell and spatial data.
#' @param he_image Optional, Xenium only. The post-run H&E image: a PNG / JPEG, or
#'   a TIFF / OME-TIFF that GDAL can read (for a JPEG2000-compressed OME-TIFF,
#'   export a smaller PNG, see \code{he_downsample}).
#' @param he_alignment Required with \code{he_image}: the \code{matrix.csv} that
#'   Xenium Explorer writes when the H&E is aligned (Image alignment > export).
#'   It maps full-resolution H&E pixels to Xenium pixels.
#' @param he_downsample Only for a PNG / JPEG export: how many times smaller it is
#'   than the full-resolution H&E that was aligned (for example 16). TIFFs are
#'   measured automatically.
#'
#' @return Invisibly, a list with \code{calls} (the table written to
#'   \code{LAMB_calls.csv}), \code{details} (full classifier output) and
#'   \code{out_dir}.
#'
#' @examples
#' # The bundled public example; results go to a temporary folder.
#' res <- lamb_run(system.file("extdata", "example_bulk.csv", package = "LAMBclassifier"),
#'                 out_dir = tempfile("lamb_"))
#' res$calls
#' list.files(res$out_dir)
#'
#' \dontrun{
#' # Your own data: results go to LAMB_results_<name> next to the input.
#' lamb_run("expression_counts.csv")
#' lamb_run("xenium_output_folder")
#' lamb_run("xenium_output_folder", he_image = "HE.ome.tif",
#'          he_alignment = "HE_alignment_files/matrix.csv")
#' lamb_run("cellranger_sample/outs")
#' lamb_run("my_seurat_object.rds")
#' }
#' @export
lamb_run <- function(input,
                     out_dir = NULL,
                     clusters = NULL,
                     data_type = c("auto", "bulk", "single_cell", "spatial"),
                     min_cells = 50,
                     run_pam50 = "auto",
                     values = c("auto", "counts", "tpm", "log2"),
                     he_image = NULL,
                     he_alignment = NULL,
                     he_downsample = NULL,
                     tumour_clusters = NULL) {
    data_type <- match.arg(data_type)
    values <- match.arg(values)
    picked <- missing(input) || is.null(input)
    if (picked) {
        if (!.lamb_can_ask()) {
            stop("No input given, and this R session cannot show a file picker (it is not interactive).\n",
                 "Give the path instead, for example: lamb_run(\"my_data.csv\")", call. = FALSE)
        }
        kind <- .lamb_ask_kind()
        input <- if (is.null(kind)) NULL else lamb_choose_input(kind)
        if (is.null(input)) {
            message("Cancelled: nothing was run.")
            return(invisible(NULL))
        }
    }
    if (!file.exists(input)) stop("Input not found: ", input, call. = FALSE)
    input <- normalizePath(input)
    log <- character()
    say <- function(...) {
        msg <- paste0(...)
        message(msg)
        log <<- c(log, msg)
    }

    platform <- if (data_type == "bulk") NULL else .lamb_detect_cells(input)
    if (is.null(platform)) {
        if (data_type %in% c("single_cell", "spatial")) {
            stop("Could not recognise ", input, " as single-cell or spatial data. Expected a 10x folder or .h5, ",
                 "an .h5ad file, a Seurat / SingleCellExperiment .rds, or a Xenium, Visium, Visium HD, ",
                 "CosMx or MERSCOPE output folder.", call. = FALSE)
        }
        data_type <- "bulk"
    } else if (data_type == "auto") {
        data_type <- if (platform %in% c("xenium", "visium", "visium_hd", "cosmx", "merscope")) "spatial" else "single_cell"
    }
    say("Input: ", input)
    say("Data type: ", data_type, if (!is.null(platform)) paste0(" (", platform, ")"))

    if (is.null(out_dir)) {
        stem <- sub("\\.(csv|tsv|txt|rds)(\\.gz)?$", "", basename(input), ignore.case = TRUE)
        out_dir <- file.path(dirname(input), paste0("LAMB_results_", stem))
    }
    dir.create(file.path(out_dir, "technical"), recursive = TRUE, showWarnings = FALSE)

    res <- if (data_type == "bulk") {
        .lamb_run_bulk(input, out_dir, run_pam50, say, values)
    } else {
        he <- list(image = he_image, alignment = he_alignment, downsample = he_downsample)
        .lamb_run_cells(input, out_dir, clusters, platform, min_cells, say, he, tumour_clusters)
    }

    spec <- lamb_specification()
    writeLines(c(log, "",
                 paste("LAMBclassifier package", as.character(utils::packageVersion("LAMBclassifier"))),
                 paste("Specification", spec$version),
                 paste("Signature SHA-256", spec$signature_sha256),
                 paste("Universe SHA-256", spec$universe_sha256),
                 paste("R", R.version.string),
                 paste("Run (UTC)", format(Sys.time(), tz = "UTC", usetz = TRUE))),
               file.path(out_dir, "technical", "run_log.txt"))
    message("Results written to: ", normalizePath(out_dir))
    if (picked || .lamb_can_ask() && identical(Sys.getenv("LAMB_LAUNCHER"), "1")) {
        main <- file.path(out_dir, if (data_type == "bulk") "LAMB_report.pdf" else "LAMB_summary.pdf")
        lamb_open(main)
        lamb_open(out_dir)
    }
    invisible(c(res, list(out_dir = normalizePath(out_dir))))
}


# ---- bulk ----

.lamb_run_bulk <- function(input, out_dir, run_pam50, say, values = "auto") {
    rd <- .lamb_read_bulk(input, say)
    if (values != "auto") {
        say("Values declared as ", values, " (values = \"", values, "\"); automatic detection skipped.")
        rd$units <- values
    }
    # Summary rows some tools and spreadsheets append: never genes.
    summ <- grepl("^(total|sum|all|__.*|n_(unmapped|multimapping|nofeature|ambiguous))$", rownames(rd$expr), ignore.case = TRUE)
    if (any(summ)) {
        say("Dropped summary row(s) that are not genes: ", paste(rownames(rd$expr)[summ], collapse = ", "), ".")
        rd$expr <- rd$expr[!summ, , drop = FALSE]
    }
    x <- rd$expr
    # A single-cell matrix saved as a table must not be scored cell by cell.
    if (ncol(x) >= 500 && mean(x == 0, na.rm = TRUE) > 0.7 && all(x == round(x), na.rm = TRUE)) {
        stop("This looks like a single-cell count matrix (", format(ncol(x), big.mark = ","), " columns, ",
             round(100 * mean(x == 0, na.rm = TRUE)), "% zeros) saved as a text table. LAMB scores groups of cells, ",
             "not single cells. Give it as a 10x folder, an .h5ad or a Seurat / SingleCellExperiment .rds file, ",
             "with cell clusters.", call. = FALSE)
    }
    counts <- if (rd$units == "unknown") .lamb_is_counts(x) else rd$units == "counts"
    x <- .lamb_check_names(.lamb_symbols(x, say, counts = counts))
    x <- .lamb_bulk_scale(x, say, units = rd$units)
    .lamb_refuse_centred(x)
    say("Samples: ", ncol(x), "; genes: ", nrow(x))

    if (identical(run_pam50, "auto")) {
        run_pam50 <- ncol(x) >= 20
        say(if (run_pam50) "PAM50 also run (20 or more samples)." else
            "PAM50 not run: it centres genes across the cohort and needs 20 or more samples.")
    } else if (isTRUE(run_pam50) && ncol(x) < 20) {
        say("Warning: PAM50 with fewer than 20 samples depends on which samples are run together.")
        if (ncol(x) < 2) { run_pam50 <- FALSE; say("PAM50 cannot run on one sample: skipped.") }
    }

    det <- LAMBclassifier(x, run_pam50 = isTRUE(run_pam50))
    calls <- data.frame(Sample = det$Sample, LAMB_class = det$Pred_Class,
                        HER2_Amplicon = det$HER2_Amplicon,
                        MA_Lo_Lean = .lamb_malo_lean(det),
                        Confidence = .lamb_confidence(det),
                        Note = .lamb_notes(det$Pred_Class, det$Flag, .lamb_low_coverage_text(det)), stringsAsFactors = FALSE)
    miss <- .lamb_missing_sig(x)[calls$Sample]
    nc <- calls$LAMB_class == "NC" & miss > 0
    ok <- calls$LAMB_class != "NC" & miss > 0
    calls$Note[nc] <- paste0("Not classifiable: ", miss[nc], " signature gene value(s) are missing in this sample, ",
                             "leaving a required signature with too few genes. Fill or re-export the missing values.")
    calls$Note[ok] <- paste0(miss[ok], " signature gene value(s) are missing in this sample; it was scored on its ",
                             "measured genes", ifelse(grepl("low_coverage", det$Flag[ok]),
                                                      ", and some signatures fell below 60% of their genes: read with care.", "."))
    if (any(miss > 0)) say(sum(miss > 0), " sample(s) have missing values in signature genes: scored on their measured ",
                           "genes (", sum(nc), " could not be classified). See Note in LAMB_calls.csv.")
    .lamb_write_common(det, calls, out_dir)
    lamb_report(det, file.path(out_dir, "LAMB_report.pdf"))
    say("Report: LAMB_report.pdf (one page per sample", if (nrow(calls) > 1) ", after a cohort overview", ").")
    lamb_gene_report(det, file.path(out_dir, "LAMB_genes.pdf"))
    say("Gene detail: LAMB_genes.pdf (one page per sample, every scored gene as its rank within the sample).")
    list(calls = calls, details = det)
}

# Per-gene centring changes within-sample gene order, which is what LAMB scores.
.lamb_refuse_centred <- function(x) {
    if (ncol(x) < 3 || min(x, na.rm = TRUE) >= 0) return(invisible())   # centred data always has negatives
    tol <- 1e-6 * max(1, max(abs(x), na.rm = TRUE))
    near0 <- function(v) mean(abs(v) < tol, na.rm = TRUE)
    if (near0(rowMeans(x, na.rm = TRUE)) > 0.8 || near0(apply(x, 1, stats::median, na.rm = TRUE)) > 0.8) {
        stop("These values look centred per gene: almost every gene averages 0 across the samples. ",
             "Centring changes the order of genes within each sample, which is what LAMB scores, so the ",
             "calls would be wrong. Give the values before centring (counts, TPM, or log2 expression).",
             call. = FALSE)
    }
    invisible()
}

# Per sample: signature genes holding a missing value (unlike a gene the platform lacks).
.lamb_missing_sig <- function(x) {
    sig <- intersect(unique(unlist(default_gene_sets())), rownames(x))
    if (!length(sig)) return(stats::setNames(integer(ncol(x)), colnames(x)))
    colSums(is.na(x[sig, , drop = FALSE]))
}

.lamb_is_counts <- function(x) {
    v <- x[!is.na(x)]
    length(v) > 0 && all(v >= 0) && all(abs(v - round(v)) < 1e-8) && max(v) > 100 && !.lamb_sums_like_tpm(x)
}

# Every sample summing to about one million is the mark of TPM; read counts never total alike.
.lamb_sums_like_tpm <- function(x) {
    s <- colSums(x, na.rm = TRUE)
    all(s > 5e5 & s < 1.01e6) && (length(s) == 1 || max(s) / min(s) < 1.02)
}

.lamb_bulk_scale <- function(x, say, units = "unknown") {
    mx <- max(x, na.rm = TRUE)
    if (units == "counts") {
        say("Raw counts: converted to TPM, then log2(TPM + 1).")
        return(log2(counts_to_tpm(x, get_gene_length()) + 1))
    }
    if (units %in% c("tpm", "fpkm")) {
        say("Values are ", toupper(units), ": applied log2(x + 1).")
        return(log2(x + 1))
    }
    if (units == "log2") {
        say("Values are log2 already: used as they are.")
        return(x)
    }
    if (!.lamb_is_counts(x) && all(abs(x - round(x)) < 1e-8, na.rm = TRUE) && .lamb_sums_like_tpm(x)) {
        say("Whole numbers that add up to about one million per sample: read as rounded TPM, applied log2(x + 1).")
        return(log2(x + 1))
    }
    if (.lamb_is_counts(x)) {
        say("Values look like raw counts: converted to TPM, then log2(TPM + 1).")
        x <- log2(counts_to_tpm(x, get_gene_length()) + 1)
    } else if (mx > 100 && min(x, na.rm = TRUE) >= 0) {
        say("Values look like TPM / FPKM (not log): applied log2(x + 1).")
        x <- log2(x + 1)
    } else if (mx > 100) {
        # Negative values mean log scale already; log2(x + 1) would make them NA.
        say("Values include negatives, so they are taken as log-transformed already and used as they are ",
            "(the largest value, ", signif(mx, 3), ", is unusually high for log data: check for a total or summary row).")
    } else {
        say("Values look log-transformed already: used as they are.")
    }
    x
}

# A file whose names are mostly not genes would otherwise look like a small targeted panel.
.lamb_check_names <- function(x) {
    known <- rownames(x) %in% c(get_gene_length()$SYMBOL, lamb_universe())
    if (mean(known) < 0.5 && sum(rownames(x) %in% lamb_universe()) < 1000) {
        stop("Most row names are not recognised as human gene symbols or Ensembl gene IDs (",
             sum(!known), " of ", length(known), ", for example ",
             paste(utils::head(rownames(x)[!known], 3), collapse = ", "), "). ",
             "Give HGNC symbols (ESR1, ERBB2) or Ensembl gene IDs (ENSG...) in the first column.", call. = FALSE)
    }
    invisible(x)
}

# Map Ensembl IDs to HGNC symbols with the bundled table; collapse duplicated
# symbols (sum for counts, the row with the highest mean otherwise).
.lamb_symbols <- function(x, say, counts) {
    ids <- sub("\\.[0-9]+$", "", rownames(x))
    if (mean(grepl("^ENSG[0-9]+$", ids)) > 0.5) {
        gl <- get_gene_length()
        sym <- gl$SYMBOL[match(ids, gl$gene_id)]
        keep <- !is.na(sym) & nzchar(sym)
        say("Gene IDs are Ensembl: mapped ", sum(keep), " of ", length(ids), " to symbols.")
        if (mean(keep) < 0.5 && sum(sym[keep] %in% lamb_universe()) < 1000) {
            stop("Only ", sum(keep), " of ", length(ids), " Ensembl IDs match human genes (for example ",
                 paste(utils::head(rownames(x)[!keep], 3), collapse = ", "), " do not). ",
                 "Check that the IDs are human Ensembl gene IDs (ENSG...), or give HGNC symbols instead.", call. = FALSE)
        }
        x <- x[keep, , drop = FALSE]
        rownames(x) <- sym[keep]
    }
    if (anyDuplicated(rownames(x))) {
        n <- sum(duplicated(rownames(x)))
        if (counts) {
            x <- rowsum(x, rownames(x), reorder = FALSE)
            say(n, " duplicated gene names: counts summed.")
        } else {
            o <- order(rowMeans(x, na.rm = TRUE), decreasing = TRUE)
            x <- x[o, , drop = FALSE]
            x <- x[!duplicated(rownames(x)), , drop = FALSE]
            say(n, " duplicated gene names: kept the row with the highest mean.")
        }
    }
    x
}


# ---- single-cell and spatial ----

.LAMB_EPITHELIAL <- c("EPCAM", "CDH1", "KRT5", "KRT7", "KRT8", "KRT14", "KRT17", "KRT18", "KRT19",
                      "CLDN4", "CLDN7", "ELF3", "GRHL2", "ESRP1", "ERBB3", "SPINT2")

# Cell-type gate. Top 3 genes per module, because a module mixes lineages
# (T, B, myeloid) and most of its genes are off in any one cell type.
.LAMB_GATE_MODULES <- list(
    epithelial = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "KRT7", "CLDN4", "CLDN7", "ELF3", "GRHL2", "ESRP1"),
    immune     = c("PTPRC", "CD3E", "CD2", "CD79A", "MS4A1", "JCHAIN", "LYZ", "CD68", "CD14", "FCER1G", "TYROBP", "NKG7"),
    vascular   = c("PECAM1", "VWF", "CDH5", "CLDN5", "EMCN", "RGS5", "PDGFRB"),
    adipocyte  = c("ADIPOQ", "PLIN1", "LEP", "CIDEA"),
    fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "FAP"))

.lamb_gate <- function(pb, say) {
    pr <- apply(pb, 2, function(v) rank(v, ties.method = "average") / length(v))
    if (is.null(dim(pr))) pr <- matrix(pr, ncol = 1, dimnames = list(rownames(pb), colnames(pb)))
    used <- lapply(.LAMB_GATE_MODULES, intersect, rownames(pb))
    used <- used[lengths(used) >= 3]
    if (!"epithelial" %in% names(used)) {
        say("Fewer than 3 epithelial gate genes on this panel: every cluster is classified (no cell-type gate).")
        return(stats::setNames(rep("epithelial", ncol(pb)), colnames(pb)))
    }
    top3 <- function(g) apply(pr[g, , drop = FALSE], 2, function(v) mean(sort(v, decreasing = TRUE)[1:3]))
    sc <- vapply(used, top3, numeric(ncol(pb)))
    if (is.null(dim(sc))) sc <- matrix(sc, nrow = 1, dimnames = list(colnames(pb), names(used)))
    top <- colnames(sc)[apply(sc, 1, which.max)]
    say("Cell-type gate: a cluster is classified unless an ", paste(setdiff(names(used), "epithelial"), collapse = ", "),
        " gene module ranks above the epithelial module in its pseudobulk (modules with 3+ measured genes).")
    stats::setNames(c(epithelial = "epithelial", immune = "immune", vascular = "vascular",
                      adipocyte = "adipocyte", fibroblast = "fibroblast-like")[top], colnames(pb))
}

.lamb_run_cells <- function(input, out_dir, clusters, platform, min_cells, say, he = list(), tumour_clusters = NULL) {
    loaded <- .lamb_read_cells(input, platform, say, user_clusters = !is.null(clusters))
    counts <- loaded$counts
    if (identical(loaded$unit, "spot")) {
        say("Visium spots hold several cells each, often of different types: a class describes a cluster ",
            "of spots, and the cell-type gate reads 'spots' where it says 'cells'.")
    }
    units <- c(cell = "cells", spot = "spots", bin = "bins")[[loaded$unit]]
    grp <- .lamb_groups(clusters, loaded$clusters, colnames(counts), say, platform, loaded$meta)
    keep <- !is.na(grp)
    if (!any(keep)) stop("No ", sub("s$", "", units), " has a cluster. Check that the cluster barcodes match the data.", call. = FALSE)
    if (any(!keep)) say(sum(!keep), " ", units, " have no cluster and are left out.")
    counts <- counts[, keep, drop = FALSE]
    grp <- factor(as.character(grp[keep]))
    say(toupper(substr(units, 1, 1)), substring(units, 2), ": ", ncol(counts), "; clusters: ", nlevels(grp), "; genes: ", nrow(counts))

    epi_genes <- intersect(.LAMB_EPITHELIAL, rownames(counts))
    if (length(epi_genes)) {
        cell_epi <- Matrix::colSums(counts[epi_genes, , drop = FALSE] > 0) > 0
        say("Epithelial_share (reported, not used to decide): share of ", units, " detecting one of ",
            paste(epi_genes, collapse = ", "), ".")
    } else {
        cell_epi <- rep(NA, ncol(counts))
        say("No epithelial marker on this panel: Epithelial_share is left empty.")
    }

    n_cells <- as.vector(table(grp)[levels(grp)])
    epi_share <- as.vector(tapply(cell_epi, grp, mean)[levels(grp)])
    pb <- vapply(levels(grp), function(g) Matrix::rowSums(counts[, grp == g, drop = FALSE]), numeric(nrow(counts)))
    if (is.null(dim(pb))) pb <- matrix(pb, ncol = 1, dimnames = list(rownames(counts), levels(grp)))
    big <- n_cells >= min_cells
    if (!any(big)) {
        stop("No cluster has ", min_cells, " or more cells (largest: ", max(n_cells), "), so nothing can be classified. ",
             "Check that the cluster file's cell ids match this data (", sum(!is.na(grp)), " cells matched), ",
             "or lower min_cells.", call. = FALSE)
    }
    log_cpm <- log2(t(t(pb[, big, drop = FALSE]) / colSums(pb[, big, drop = FALSE])) * 1e6 + 1)
    say("Pseudobulk per cluster: summed counts, log2(CPM + 1); clusters under ", min_cells, " cells are not classified.")

    det <- LAMBclassifier(log_cpm, run_pam50 = FALSE)
    cls <- stats::setNames(rep("Too few cells", nlevels(grp)), levels(grp))
    cls[det$Sample] <- det$Pred_Class
    flag <- stats::setNames(rep("", nlevels(grp)), levels(grp))
    flag[det$Sample] <- det$Flag
    amp <- stats::setNames(rep("", nlevels(grp)), levels(grp))
    amp[det$Sample] <- det$HER2_Amplicon
    lean <- stats::setNames(rep("", nlevels(grp)), levels(grp))
    lean[det$Sample] <- .lamb_malo_lean(det)
    gate <- .lamb_gate(pb, say)
    if (length(tumour_clusters)) {
        bad <- setdiff(tumour_clusters, levels(grp))
        if (length(bad)) stop("tumour_clusters names cluster(s) not in the data: ", paste(bad, collapse = ", "), call. = FALSE)
        gate[tumour_clusters] <- "epithelial"
        say("Classified as tumour whatever their markers (tumour_clusters): ", paste(tumour_clusters, collapse = ", "), ".")
    }
    epithelial <- unname(gate[levels(grp)] == "epithelial")
    calls <- data.frame(Cluster = levels(grp), Cells = n_cells,
                        Cell_type = unname(gate[levels(grp)]),
                        Epithelial_share = round(epi_share, 2),
                        LAMB_class = unname(cls),
                        HER2_Amplicon = unname(amp),
                        MA_Lo_Lean = unname(lean),
                        Confidence = "",
                        Note = .lamb_notes(unname(cls), unname(flag), .lamb_low_coverage_text(det)), stringsAsFactors = FALSE)
    conf <- stats::setNames(.lamb_confidence(det), det$Sample)
    calls$Confidence <- ifelse(calls$Cluster %in% names(conf), unname(conf[calls$Cluster]), "")
    calls$Note[!epithelial] <- paste0("Not tumour epithelium (", calls$Cell_type[!epithelial], " genes dominate): ",
                                      "the class does not apply. If this is tumour, name it with tumour_clusters =.")
    calls$Note[calls$LAMB_class == "Too few cells"] <- paste0("Fewer than ", min_cells, " cells: not classified.")
    # Never show an epithelial class on non-epithelial clusters.
    calls$LAMB_class[!epithelial] <- "Not epithelial"
    calls$HER2_Amplicon[!epithelial] <- ""
    calls$MA_Lo_Lean[!epithelial] <- ""
    calls$Confidence[!epithelial] <- ""
    calls <- calls[order(calls$Cell_type != "epithelial", -calls$Cells), ]
    .lamb_write_common(det, calls, out_dir)

    cell_class <- cls[as.character(grp)]
    cell_class[!epithelial[match(as.character(grp), levels(grp))]] <- "Not epithelial"
    names(cell_class) <- colnames(counts)
    slide <- .lamb_slide_summary(calls, say)
    utils::write.csv(slide$table, file.path(out_dir, "LAMB_slide.csv"), row.names = FALSE)
    xy <- NULL
    img <- NULL
    if (!is.null(loaded$xy)) {
        xy <- loaded$xy[match(colnames(counts), loaded$xy$cell_id), ]
        if (platform == "xenium") {
            .lamb_write_explorer(cell_class, file.path(out_dir, "LAMB_xenium_explorer.csv"))
            if (!is.null(he$image)) {
                img <- .lamb_he_layer(he, xy, loaded$pixel_size, say)
                xy[, c("x", "y")] <- img$xy
            }
        }
    }
    if (!is.null(he$image) && platform != "xenium") {
        say("he_image is used for Xenium only (it needs Xenium Explorer's alignment file): ignored.")
    }
    .lamb_pdf_cells(calls, cell_class, xy, slide, img, file.path(out_dir, "LAMB_summary.pdf"), say)
    lamb_report_cells(calls, det, file.path(out_dir, "LAMB_clusters.pdf"))
    say("Reports: LAMB_summary.pdf (slide overview) and LAMB_clusters.pdf (why each cluster got its class).")
    shown <- calls$Cluster[calls$Cell_type == "epithelial" & calls$Cluster %in% det$Sample]
    if (length(shown)) {
        lamb_gene_report(det, file.path(out_dir, "LAMB_genes.pdf"), samples = shown)
        say("Gene detail: LAMB_genes.pdf (one page per classified cluster, every scored gene as its rank within the cluster).")
    }
    list(calls = calls, slide = slide$table, details = det)
}

.lamb_groups <- function(user, found, cells, say, platform = NULL, meta = NULL) {
    if (!is.null(user)) {
        if (is.character(user) && length(user) == 1 && !file.exists(user)) {
            if (!is.null(meta) && user %in% names(meta)) {
                v <- meta[[user]]
                say("Groups: the cell annotation column '", user, "' (", length(unique(stats::na.omit(v))), " groups).")
                return(as.character(v[match(cells, rownames(meta))]))
            }
            stop("clusters = \"", user, "\" is neither a file nor a cell annotation column of this data",
                 if (length(names(meta))) paste0(" (available: ", paste(utils::head(names(meta), 12), collapse = ", "), ")") else "",
                 ".", call. = FALSE)
        }
        if (is.character(user) && length(user) == 1) {
            t <- utils::read.csv(user, stringsAsFactors = FALSE)
            user <- stats::setNames(as.character(t[[2]]), t[[1]])
            say("Groups: your cluster file.")
        } else {
            say("Groups: the clusters you supplied.")
        }
        if (is.null(names(user))) stop("`clusters` needs cell barcodes as names.", call. = FALSE)
        return(unname(user[cells]))
    }
    if (is.null(found)) {
        stop(if (is.null(platform)) "No clusters found. Supply them with clusters = 'my_clusters.csv' (two columns: cell id, cluster)."
             else .lamb_no_clusters_text(platform, names(meta)), call. = FALSE)
    }
    unname(found[cells])
}


# ---- shared output ----

.lamb_notes <- function(cls, flag, low_cov = "") {
    note <- rep("", length(cls))
    flag[is.na(flag)] <- ""
    note[grepl("close_call", flag)] <- "Close call: two classes scored almost the same."
    lc <- grepl("low_coverage", flag)
    note[lc] <- paste(note[lc], if (nzchar(low_cov)) paste0("Few genes measured in this data for ", low_cov, ".")
                      else "Some signature genes are missing in this data.")
    note[cls == "NC"] <- "Not classifiable: a required signature had too few genes in this data to be scored."
    trimws(note)
}

.lamb_write_common <- function(det, calls, out_dir) {
    utils::write.csv(calls, file.path(out_dir, "LAMB_calls.csv"), row.names = FALSE, na = "")
    suppressMessages(export_results(det, file.path(out_dir, "technical", "LAMB_details.csv")))
    utils::write.table(attr(det, "coverage"), file.path(out_dir, "technical", "signature_coverage.tsv"),
                       sep = "\t", quote = FALSE, row.names = FALSE)
    pct <- attr(det, "gene_percentiles")
    if (!is.null(pct)) {
        gp <- data.frame(id = rownames(pct), pct, check.names = FALSE, stringsAsFactors = FALSE)
        names(gp)[1] <- if ("Cluster" %in% names(calls)) "Cluster" else "Sample"
        utils::write.csv(gp, file.path(out_dir, "technical", "LAMB_gene_percentiles.csv"), row.names = FALSE)
    }
}

.lamb_colours <- function() {
    c(LAMB_CLASS_COLORS, "Not epithelial" = "#d9d9d9", "Too few cells" = "#f0f0f0")
}

.lamb_write_explorer <- function(cell_class, path) {
    out <- data.frame(cell_id = names(cell_class), group = unname(cell_class),
                      color = unname(.lamb_colours()[cell_class]), stringsAsFactors = FALSE)
    out$color[is.na(out$color)] <- "#888888"
    stopifnot(!anyDuplicated(out$cell_id))   # Xenium Explorer refuses a cell in two groups
    utils::write.csv(out, path, row.names = FALSE)
}

.lamb_footer <- function() {
    graphics::mtext("Research use only: expression-defined labels, not a clinical diagnosis.",
                    side = 1, line = 0.5, adj = 0, cex = 0.7, col = "#777777", outer = TRUE)
}

# Slide-level class: most epithelial cells; "mixed" when none reaches half.
.lamb_slide_summary <- function(calls, say) {
    epi <- calls[calls$Cell_type == "epithelial" & calls$LAMB_class %in% LAMB_CLASS_ORDER, ]
    lv <- intersect(LAMB_CLASS_ORDER, unique(epi$LAMB_class))
    n <- tapply(epi$Cells, factor(epi$LAMB_class, levels = LAMB_CLASS_ORDER), sum)
    n[is.na(n)] <- 0
    pct <- if (sum(n)) 100 * n / sum(n) else n
    top <- names(which.max(pct))
    call <- if (!sum(n)) "none (no classified epithelial cluster)" else if (max(pct) >= 50) top else {
        o <- order(-pct); o <- o[pct[o] >= 10 | seq_along(o) <= 2]   # every class at 10% or more, at least two
        paste0("mixed: ", paste0(names(pct)[o], " ", round(pct[o]), "%", collapse = " / "))
    }
    say("Slide-level class: ", call, " (class with most classified epithelial cells; 'mixed' below 50%, free choice).")
    tab <- data.frame(Slide_class = call, Epithelial_cells_classified = sum(n),
                      Epithelial_clusters_classified = nrow(epi), stringsAsFactors = FALSE, check.names = FALSE)
    for (k in LAMB_CLASS_ORDER[LAMB_CLASS_ORDER != "NC"]) tab[[paste0(k, " (%)")]] <- round(unname(pct[k]), 1)
    tab[["NC (%)"]] <- round(unname(pct["NC"]), 1)
    list(table = tab, call = call, pct = pct[pct > 0])
}

.lamb_pdf_cells <- function(calls, cell_class, xy, slide, img, path, say) {
    pct <- slide$pct
    grDevices::pdf(path, width = 11.69, height = 8.27, bg = "white")
    on.exit(grDevices::dev.off())
    sub <- paste0(format(sum(calls$Cells[calls$Cell_type == "epithelial"]), big.mark = ","),
                  " epithelial cells in ", sum(calls$Cell_type == "epithelial"),
                  " clusters; every cluster in LAMB_calls.csv.")
    if (!is.null(xy)) {
        if (!is.null(img)) {
            graphics::layout(matrix(1:3, 1), widths = c(1.5, 1.5, 0.8))
        } else {
            graphics::layout(matrix(1:2, 1), widths = c(2, 1))
        }
        graphics::par(oma = c(2.5, 1.5, 5, 1.5), mar = c(1, 1, 2, 1))
        set.seed(1)
        idx <- which(!is.na(xy$x))
        if (length(idx) > 80000) idx <- sort(sample(idx, 80000))
        # Larger dots for small samples so they stay visible.
        grow <- min(6, sqrt(50000 / max(length(idx), 1)))
        cc <- cell_class[idx]
        col <- .lamb_colours()[cc]
        col[is.na(col)] <- "#888888"
        rx <- range(xy$x[idx]); ry <- range(xy$y[idx])
        if (!is.null(img)) {
            .lamb_draw_he(img, rx, ry); graphics::title("H&E", cex.main = 1, font.main = 1)
            .lamb_draw_he(img, rx, ry)
            keep <- cc != "Not epithelial"
            graphics::points(xy$x[idx][keep], xy$y[idx][keep], pch = 16, cex = max(0.22, min(1, 0.22 * grow)),
                             col = grDevices::adjustcolor(col[keep], 0.8))
            graphics::title("LAMB classes on the H&E", cex.main = 1, font.main = 1)
        } else {
            o <- order(cc != "Not epithelial")
            graphics::plot(xy$x[idx][o], xy$y[idx][o], pch = 16, cex = max(0.15, min(0.9, 0.15 * grow)), col = col[o], asp = 1, axes = FALSE,
                           xlab = "", ylab = "", ylim = rev(ry))
            graphics::title("LAMB classes on the tissue", cex.main = 1, font.main = 1)
        }
        leg <- c(names(pct), "Not epithelial")
        graphics::legend("topleft", legend = leg, fill = .lamb_colours()[leg], border = NA, bty = "n", cex = 0.8,
                         bg = grDevices::adjustcolor("white", 0.7))
    } else {
        graphics::par(oma = c(2.5, 1.5, 5, 1.5))
    }
    k <- rev(names(pct)); v <- rev(unname(pct))
    pad <- max(0, 8 - length(v))   # empty slots below keep bars slim and top-aligned
    k <- c(rep("", pad), k); v <- c(rep(NA, pad), v)
    bcol <- .lamb_colours()[k]; bcol[is.na(bcol)] <- "#888888"
    graphics::par(mar = c(2, 7, 2, 3.5))
    b <- graphics::barplot(v, horiz = TRUE, names.arg = k, las = 1, col = bcol, border = NA, xlim = c(0, 100),
                           axes = FALSE, cex.names = 0.85, space = 0.6)
    if (any(!is.na(v))) {
        graphics::text(v[!is.na(v)], b[!is.na(v)], paste0(round(v[!is.na(v)], 1), "%"), pos = 4, cex = 0.85, xpd = NA)
    } else {
        graphics::text(50, mean(b), "No tumour-epithelium cluster\nwas classified", cex = 0.85, col = "#666666")
    }
    graphics::title("% of epithelial cells", cex.main = 1, font.main = 1)
    graphics::mtext(paste("Slide class:", sub("^mixed: (.*)$", "mixed (\\1)", slide$call)), side = 3, line = 2.8, adj = 0, cex = 1.35, font = 2, outer = TRUE)
    graphics::mtext("LAMB classes of the epithelial cells", side = 3, line = 1.4, adj = 0, cex = 0.9, outer = TRUE)
    graphics::mtext(sub, side = 3, line = 0.3, adj = 0, cex = 0.7, col = "#666666", outer = TRUE)
    .lamb_footer()
}


# ---- Individual gene expression: within-sample percentile ranks ----

.LAMB_GP_MARKERS <- c("ESR1", "PGR", "AR", "ERBB2", "FOXA1", "GATA3", "KRT5", "TP63", "MKI67")

.LAMB_GP_SUBTITLE <- "Rank within this sample (100 = highest of %s genes measured). RNA, not IHC."

# Percentile among the sample's measured universe genes, 100 = highest, ties
# averaged (the HER2 amplicon ranking, read from the top). Genes outside the
# universe or missing are NA.
.lamb_gp_percentiles <- function(data, genes, universe = lamb_universe()) {
    if (!is.matrix(data)) data <- as.matrix(data)
    data <- data[!duplicated(rownames(data)), , drop = FALSE]   # first occurrence, as in lamb_score()
    mat <- data[intersect(universe, rownames(data)), , drop = FALSE]
    genes <- unique(as.character(genes))
    ids <- colnames(mat) %||% as.character(seq_len(ncol(mat)))
    pct <- matrix(NA_real_, ncol(mat), length(genes), dimnames = list(ids, genes))
    n <- stats::setNames(integer(ncol(mat)), ids)
    for (j in seq_len(ncol(mat))) {
        v <- mat[, j]
        names(v) <- rownames(mat)
        v <- v[!is.na(v)]
        n[j] <- length(v)
        if (length(v) < 2L) next
        p <- 100 * (rank(v, ties.method = "average") - 1) / (length(v) - 1)
        hit <- intersect(genes, names(v))
        pct[j, hit] <- p[hit]
    }
    list(pct = round(pct, 1), n = n)
}

.lamb_gp_subtitle <- function(n) {
    sprintf(.LAMB_GP_SUBTITLE, if (is.na(n)) "the" else format(n, big.mark = ","))
}

# One sample's percentiles (named, NA = not measured) and N, or NULL.
.lamb_gp_row <- function(results, sample) {
    pct <- attr(results, "gene_percentiles")
    if (is.null(pct) || !sample %in% rownames(pct)) return(NULL)
    n <- attr(results, "n_genes_ranked")
    list(pct = pct[sample, ], n = if (is.null(n)) NA_integer_ else unname(n[sample]))
}

# Gene, bar 0-100, value. Returns the y below the row.
.lamb_gp_bar <- function(gene, value, x0, x1, y, lab_w, cex = 0.8, font = 1, col = .LAMB_INK, lwd = 4) {
    num_w <- graphics::strwidth("100", units = "user", cex = cex) + 0.08
    bx0 <- x0 + lab_w; bx1 <- x1 - num_w
    graphics::text(x0, y, gene, adj = c(0, 0.5), cex = cex, font = font, col = col)
    graphics::segments(bx0, y, bx1, y, col = "#e8e8e8", lwd = lwd, lend = 1)
    graphics::segments(bx0, y, bx0 + (bx1 - bx0) * value / 100, y, col = "#7d7d7d", lwd = lwd, lend = 1)
    graphics::text(x1, y, round(value), adj = c(1, 0.5), cex = cex, col = col)
}

# Marker strip for the bulk page; measured markers only, nothing when none.
.lamb_gp_strip <- function(results, sample, x0, x1, ytop) {
    g <- .lamb_gp_row(results, sample)
    if (is.null(g)) return(ytop)
    v <- g$pct[intersect(.LAMB_GP_MARKERS, names(g$pct))]
    v <- v[!is.na(v)]
    if (!length(v)) return(ytop)
    rh <- min(0.2, (ytop - 0.42 - 0.85) / max(length(v) - 1, 1))   # stay clear of the footer
    if (rh < 0.13) return(ytop)
    .lamb_rpt_label(x0, ytop, "Marker genes")
    sub <- .lamb_gp_subtitle(g$n)
    cex <- .lamb_fit_cex(sub, x1 - x0 + 0.4, 0.72, min_cex = 0.5)   # may run into the gutter, not the bars
    graphics::text(x0, ytop - 0.2, sub, adj = c(0, 0.5), cex = cex, col = .LAMB_MUTED)
    y <- ytop - 0.42
    for (k in names(v)) {
        .lamb_gp_bar(k, v[[k]], x0, x1, y, lab_w = 0.7, cex = 0.8, lwd = 4)
        y <- y - rh
    }
    y
}


# ---- gene-detail report ----

.LAMB_GP_COLUMNS <- list(c("chondroid", "squamous"), "basal_like", "luminal", "apocrine", c("malo", "her2"))

#' Write a gene-detail PDF: one page per sample or cluster
#'
#' @description
#' Shows the expression of every gene LAMB scores, as its within-sample
#' percentile rank: the rank among the sample's measured reference-universe
#' genes, 100 = highest, 0 = lowest. One A4 landscape page per sample, with
#' each class signature and the HER2 amplicon module as a group, genes sorted
#' by rank. The groups behind the predicted class are in bold.
#'
#' Ranks are within one sample, so they read the same on any platform; they
#' are RNA levels, not IHC results. Genes not measured are left out.
#'
#' @param results Output of [LAMBclassifier()], with its
#'   \code{"gene_percentiles"} attribute.
#' @param file Path of the PDF to write.
#' @param samples Optional. The samples (or clusters) to include, in order.
#'   Default: all.
#'
#' @return Invisibly, \code{file}.
#'
#' @examples
#' \dontrun{
#' res <- LAMBclassifier(expr_log2)
#' lamb_gene_report(res, "LAMB_genes.pdf")
#' }
#' @export
lamb_gene_report <- function(results, file, samples = NULL) {
    if (!is.data.frame(results) || !all(c("Sample", "Pred_Class") %in% names(results))) {
        stop("`results` must be the output of LAMBclassifier() (columns Sample and Pred_Class).", call. = FALSE)
    }
    if (is.null(attr(results, "gene_percentiles"))) {
        stop("`results` has no \"gene_percentiles\" attribute: use the LAMBclassifier() output as it is ",
             "(subsetting or joining with dplyr drops it).", call. = FALSE)
    }
    rows <- as.data.frame(results, stringsAsFactors = FALSE)
    if (!is.null(samples)) rows <- rows[match(as.character(samples), rows$Sample), , drop = FALSE]
    rows <- rows[!is.na(rows$Sample), , drop = FALSE]
    if (!nrow(rows)) stop("No sample to report.", call. = FALSE)
    sigs <- attr(results, "signatures") %||% default_gene_sets()
    amp <- attr(results, "specification")$her2_amplicon$genes %||% DEFAULT_HER2_AMPLICON_GENES
    groups <- c(sigs, list(her2 = amp))

    grDevices::pdf(file, width = .LAMB_PAGE[["w"]], height = .LAMB_PAGE[["h"]], bg = "white",
                   title = "LAMB gene detail", useDingbats = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (i in seq_len(nrow(rows))) {
        .lamb_new_page()
        .lamb_gp_page(as.list(rows[i, , drop = FALSE]), results, groups)
        .lamb_footer_page(results, i, nrow(rows))
    }
    invisible(file)
}

.lamb_gp_page <- function(r, results, groups) {
    W <- .LAMB_PAGE[["w"]]; H <- .LAMB_PAGE[["h"]]; x0 <- .LAMB_X0
    sample <- as.character(r[["Sample"]])
    g <- .lamb_gp_row(results, sample)
    graphics::text(x0, H - 0.45, "Gene expression detail", adj = c(0, 0.5), cex = 0.8, col = .LAMB_MUTED)
    graphics::text(x0, H - 0.9, sample, adj = c(0, 0.5), cex = 1.6, font = 2)
    cls <- as.character(r[["Pred_Class"]])
    if (!length(cls) || is.na(cls)) cls <- "NC"
    chip <- .lamb_chip(x0, H - 1.2, if (cls == "NC") "Not classifiable" else cls, cex = 1.3, pad = 0.1)
    tx <- x0 + chip[["w"]] + 0.3
    cy <- H - 1.2 - chip[["h"]] / 2
    graphics::text(tx, cy + 0.1, .lamb_gp_subtitle(if (is.null(g)) NA else g$n), adj = c(0, 0.5),
                   cex = 0.85, col = .LAMB_MUTED)
    graphics::text(tx, cy - 0.12, "Genes of each LAMB program, highest first. Bold: the program(s) behind this class.",
                   adj = c(0, 0.5), cex = 0.85, col = .LAMB_MUTED)
    if (is.null(g)) return(invisible())

    key <- c("Chondroid" = "chondroid", "Squamous" = "squamous", "Basal-like" = "basal_like",
             "LHS" = "luminal", "Apocrine" = "apocrine", "MA-Lo" = "malo")[cls]
    amp_high <- identical(as.character(r[["HER2_Amplicon"]]), "high")
    label <- c(.LAMB_SIG_LABEL, her2 = "HER2 amplicon")

    ncol <- length(.LAMB_GP_COLUMNS); gap <- 0.3
    cw <- (W - 2 * x0 - (ncol - 1) * gap) / ncol
    ytop <- H - 1.2 - chip[["h"]] - 0.4
    rh <- 0.21
    for (ci in seq_len(ncol)) {
        gx0 <- x0 + (ci - 1) * (cw + gap); gx1 <- gx0 + cw
        y <- ytop
        for (k in .LAMB_GP_COLUMNS[[ci]]) {
            genes <- groups[[k]]
            if (is.null(genes)) next
            v <- g$pct[intersect(genes, names(g$pct))]
            v <- sort(v[!is.na(v)], decreasing = TRUE)
            if (!length(v)) next
            on <- (!is.na(key) && k == key) || (k == "her2" && amp_high)
            col <- if (on) .LAMB_INK else .LAMB_MUTED
            graphics::text(gx0, y, label[[k]], adj = c(0, 0.5), cex = 1, font = if (on) 2 else 1, col = col)
            graphics::segments(gx0, y - 0.13, gx1, y - 0.13, col = if (on) .LAMB_INK else .LAMB_FAINT, lwd = 0.6)
            y <- y - 0.34
            for (gene in names(v)) {
                .lamb_gp_bar(gene, v[[gene]], gx0, gx1, y, lab_w = 0.85, cex = 0.8, col = col, lwd = 3.5)
                y <- y - rh
            }
            y <- y - 0.3
        }
    }
    invisible()
}

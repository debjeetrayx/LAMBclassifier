# =============================================================================
# LAMBclassifier: Visualization Functions
# =============================================================================

utils::globalVariables(c(
    "Mean_SQ", "Squamous_Score", "Mean_SS", "Sensory_Score",
    "Pred_Class", "Score", "score_axis",
    "Score_Sensory", "Score_Basal_Like", "Score_LHS", "Score_Apocrine", "Score_MALO_HER2",
    "Score_Squamous", "Score_Chondroid",
    "Y_Expr", "X_Expr",
    "scaled", "Density", "reference", "prediction", "Reference", "Prediction",
    "Count", "parent", "x_start", "y_start", "x_end", "y_end",
    "nudge_val", "full_label", "hjust_val",
    "display_label", "fill_col", "nudge_y", "v_just",
    "meta_signal", "max_mammary", "max_lum_apo", "edge_type"
))

#' Plot Decision Boundaries
#'
#' @description Creates scatter plots showing sample positions at each decision boundary gate.
#' @param results Data frame from LAMBclassifier().
#' @param metadata Optional metadata data frame.
#' @param output_prefix File path prefix for output files.
#' @return Invisible NULL.
#' @export
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
plot_decision_boundaries <- function(results, metadata = NULL, output_prefix) {
    if (!is.null(metadata)) {
        plot_data <- dplyr::left_join(results, metadata, by = "Sample")
        color_var <- if ("subtype" %in% colnames(metadata)) "subtype" else "Pred_Class"
    } else {
        plot_data <- results
        color_var <- "Pred_Class"
    }

    my_theme <- theme_minimal() +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))

    get_limits <- function(x, y) {
        val_range <- range(c(x, y), na.rm = TRUE)
        if (is.infinite(val_range[1])) return(c(-1, 1))
        padding <- diff(val_range) * 0.05
        if (padding == 0) padding <- 0.5
        c(val_range[1] - padding, val_range[2] + padding)
    }

    plot_data$meta_signal <- pmax(plot_data$Score_Chondroid, plot_data$Score_Squamous)
    plot_data$max_mammary <- pmax(plot_data$Score_Basal_Like, plot_data$Score_LHS, plot_data$Score_Apocrine)
    plot_data$max_lum_apo <- pmax(plot_data$Score_LHS, plot_data$Score_Apocrine)

    # 1. Metaplastic Gate
    lims1 <- get_limits(plot_data$max_mammary, plot_data$meta_signal)
    p1 <- ggplot(plot_data, aes(x = max_mammary, y = meta_signal)) +
        geom_point(aes(color = .data[[color_var]]), alpha = 0.7, size = 2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
        labs(title = "1. Metaplastic Gate", x = "max(Basal, LHS, Apocrine)", y = "max(Chondroid, Squamous)") +
        coord_fixed(xlim = lims1, ylim = lims1) +
        my_theme +
        scale_color_manual(values = LAB_CLASS_COLORS)

    # 2. Basal-like Gate
    dat2 <- plot_data[!plot_data$Pred_Class %in% c("Chondroid", "Squamous"), ]
    if (nrow(dat2) > 0) {
        lims2 <- get_limits(dat2$max_lum_apo, dat2$Score_Basal_Like)
        p2 <- ggplot(dat2, aes(x = max_lum_apo, y = Score_Basal_Like)) +
            geom_point(aes(color = .data[[color_var]]), alpha = 0.7, size = 2) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
            labs(title = "2. Basal-like Gate", x = "max(LHS, Apocrine)", y = "Basal-like Score") +
            coord_fixed(xlim = lims2, ylim = lims2) +
            my_theme +
            scale_color_manual(values = LAB_CLASS_COLORS)
    } else {
        p2 <- ggplot() + theme_void()
    }

    # 3. Receptor Gate
    dat3 <- plot_data[plot_data$Pred_Class %in% c("LHS", "MA-Lo (HER2)", "Apocrine"), ]
    if (nrow(dat3) > 0) {
        lims3 <- get_limits(dat3$Score_Apocrine, dat3$Score_LHS)
        p3 <- ggplot(dat3, aes(x = Score_Apocrine, y = Score_LHS)) +
            geom_point(aes(color = .data[[color_var]]), alpha = 0.7, size = 2) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 1) +
            labs(title = "3. Receptor Gate", x = "Apocrine Score", y = "LHS Score") +
            coord_fixed(xlim = lims3, ylim = lims3) +
            my_theme +
            scale_color_manual(values = LAB_CLASS_COLORS)
    } else {
        p3 <- ggplot() + theme_void()
    }

    pdf_file <- paste0(output_prefix, "_Decision_Boundaries.pdf")
    png_file <- paste0(output_prefix, "_Decision_Boundaries.png")

    pdf(pdf_file, width = 12, height = 4)
    gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
    dev.off()

    png(png_file, width = 1200, height = 400, res = 100)
    gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
    dev.off()

    message("Saved: ", pdf_file, " and ", png_file)
    invisible(NULL)
}

#' Plot Lineage Tree
#'
#' @description Creates a publication-grade visualization of the classification hierarchy.
#' @param results Data frame from LAMBclassifier().
#' @param metadata Optional metadata data frame.
#' @param output_prefix File path prefix for output files.
#' @return Invisible NULL.
#' @export
#' @import ggplot2
plot_lineage_tree <- function(results, metadata = NULL, output_prefix) {
    plot_data <- if (!is.null(metadata)) dplyr::left_join(results, metadata, by = "Sample") else results
    total_n <- nrow(plot_data)

    cnt <- function(cls) sum(plot_data$Pred_Class == cls, na.rm = TRUE)
    pct <- function(n) sprintf("%.1f%%", 100 * n / total_n)

    n_chondroid  <- cnt("Chondroid")
    n_squamous   <- cnt("Squamous")
    n_basal_like <- cnt("Basal-like")
    n_lhs        <- cnt("LHS")
    n_malow      <- cnt("MA-Lo (HER2)")
    n_apocrine   <- cnt("Apocrine")

    n_meta_node <- n_chondroid + n_squamous
    n_diff_lum  <- n_lhs + n_malow + n_apocrine
    n_lasp_node <- n_basal_like + n_diff_lum

    nodes <- tibble::tibble(
        id = c("MaSC", "Metaplastic_Branch", "Chondroid", "Squamous",
               "LASP_Anchor", "Basal_like", "Differentiated_Luminal",
               "LHS", "MA_low", "Apocrine"),
        title = c("Mammary Stem Cell\n(MaSC)",
                  "Metaplastic Lineage", "Chondroid", "Squamous",
                  "LASP Anchor\n(Luminal Progenitor)",
                  "Basal-like\n(Progenitor State)",
                  "Differentiated Luminal\nLineage",
                  "LHS\n(Luminal ER+)",
                  "MA-Lo (HER2)\n(Amplicon expression)",
                  "Apocrine\n(Molecular Apocrine)"),
        stat = c(paste0("Root (N=", total_n, ")"),
                 paste0(pct(n_meta_node), " (n=", n_meta_node, ")"),
                 paste0(pct(n_chondroid), " (n=", n_chondroid, ")"),
                 paste0(pct(n_squamous), " (n=", n_squamous, ")"),
                 paste0("Bipotent Anchor (n=", n_lasp_node, ")"),
                 paste0(pct(n_basal_like), " (n=", n_basal_like, ")"),
                 paste0("Committed (n=", n_diff_lum, ")"),
                 paste0(pct(n_lhs), " (n=", n_lhs, ")"),
                 paste0(pct(n_malow), " (n=", n_malow, ")"),
                 paste0(pct(n_apocrine), " (n=", n_apocrine, ")")),
        x = c(0, 2.0, 4.2, 4.2, 2.0, 4.2, 4.2, 6.4, 6.4, 6.4),
        y = c(0, 3.0, 4.0, 2.0, -2.0, 0.0, -3.0, -1.7, -3.0, -4.3),
        parent = c(NA, "MaSC", "Metaplastic_Branch", "Metaplastic_Branch",
                   "MaSC", "LASP_Anchor", "LASP_Anchor",
                   "Differentiated_Luminal", "Differentiated_Luminal", "Differentiated_Luminal")
    )

    nodes$fill_col <- sapply(nodes$id, function(id) {
        if (id == "Chondroid") return(LAB_CLASS_COLORS["Chondroid"])
        if (id == "Squamous") return(LAB_CLASS_COLORS["Squamous"])
        if (id == "Basal_like") return(LAB_CLASS_COLORS["Basal-like"])
        if (id == "LHS") return(LAB_CLASS_COLORS["LHS"])
        if (id == "MA_low") return(LAB_CLASS_COLORS["MA-Lo (HER2)"])
        if (id == "Apocrine") return(LAB_CLASS_COLORS["Apocrine"])
        if (id == "LASP_Anchor") return("#FFF8DC")
        if (id == "MaSC") return("#F4F6F7")
        return("#FAFAFA")
    })

    nodes$border_col <- sapply(nodes$id, function(id) {
        if (id %in% c("Chondroid", "Squamous", "Basal_like", "LHS", "MA_low", "Apocrine")) return("gray20")
        if (id == "LASP_Anchor") return("#B8860B")
        return("gray50")
    })

    nodes$text_col <- sapply(nodes$id, function(id) {
        if (id %in% c("LHS", "Basal_like", "Chondroid", "Apocrine")) return("white")
        return("black")
    })

    w <- 1.45
    h <- 0.70
    nodes$xmin <- nodes$x - w / 2
    nodes$xmax <- nodes$x + w / 2
    nodes$ymin <- nodes$y - h / 2
    nodes$ymax <- nodes$y + h / 2

    edges <- dplyr::filter(nodes, !is.na(parent))
    edges <- dplyr::select(edges, from = parent, to = id, x_end = x, y_end = y)
    edges <- dplyr::left_join(edges, dplyr::select(nodes, id, x, y), by = c("from" = "id"))
    edges <- dplyr::rename(edges, x_start = x, y_start = y)
    edges$x_start_clean <- edges$x_start + w / 2
    edges$x_end_clean   <- edges$x_end - w / 2
    edges$edge_type     <- ifelse(edges$to == "MA_low", "dashed", "solid")

    n_nc         <- cnt("NC")
    sub_title_text <- "Proposed differentiation hierarchy (including ERBB2-amplicon expression)"
    if (n_nc > 0) sub_title_text <- paste0(sub_title_text, sprintf(" | Unclassified (NC): %d (%.1f%%)", n_nc, 100 * n_nc / total_n))

    p <- ggplot() +
        geom_segment(data = edges, aes(x = x_start_clean, y = y_start, xend = x_end_clean, yend = y_end, linetype = edge_type),
                     color = "gray40", linewidth = 0.8, arrow = arrow(length = unit(0.20, "cm"), type = "closed")) +
        scale_linetype_identity() +
        geom_rect(data = nodes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill_col, color = border_col),
                  linewidth = 0.75) +
        scale_fill_identity() + scale_color_identity() +
        geom_text(data = nodes, aes(x = x, y = y + 0.09, label = title, color = text_col),
                  fontface = "bold", size = 2.9, lineheight = 0.85) +
        geom_text(data = nodes, aes(x = x, y = y - 0.20, label = stat, color = text_col),
                  fontface = "italic", size = 2.3) +
        theme_void() +
        theme(plot.margin = margin(15, 15, 15, 15),
              plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
              plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30")) +
        coord_cartesian(xlim = c(-1.0, 7.5), ylim = c(-5.2, 4.8)) +
        labs(title = "LAMBclassifier Lineage Tree", subtitle = sub_title_text)

    ggsave(paste0(output_prefix, "_Lineage_Tree.pdf"), p, width = 12, height = 7.5)
    ggsave(paste0(output_prefix, "_Lineage_Tree.png"), p, width = 12, height = 7.5, dpi = 300)
    message("Saved publication lineage tree: ", paste0(output_prefix, "_Lineage_Tree.pdf"))
    invisible(NULL)
}

#' Plot Comprehensive Expression Heatmap
#'
#' @description Generates an annotated expression heatmap of signature genes organized by LAMB class.
#' @param expr_mat Expression matrix (log2 TPM).
#' @param results Data frame from LAMBclassifier().
#' @param metadata Optional metadata data frame.
#' @param genes Optional gene vector.
#' @param output_prefix Output file path prefix.
#' @param show_names Logical, whether to show sample names.
#' @param cluster_cols Logical, whether to cluster columns.
#' @param scale_rows Logical, whether to center/scale rows.
#' @param cluster_rows_by_class Logical, whether to order genes by signature class.
#' @param lab_column Optional name of comparison column in metadata.
#' @param add_top_n_markers Reserved compatibility argument. Additional markers
#'   are not added by the locked implementation.
#' @return Invisible NULL.
#' @export
#' @import pheatmap
plot_comprehensive_heatmap <- function(expr_mat,
                                       results,
                                       metadata = NULL,
                                       genes = NULL,
                                       output_prefix,
                                       show_names = FALSE,
                                       cluster_cols = FALSE,
                                       scale_rows = FALSE,
                                       cluster_rows_by_class = FALSE,
                                       lab_column = NULL,
                                       add_top_n_markers = 0) {
    common_samples <- intersect(colnames(expr_mat), results$Sample)
    expr_mat <- expr_mat[, common_samples]
    saved_sigs <- attr(results, "signatures")
    results <- results[match(common_samples, results$Sample), ]

    if (!is.null(metadata)) {
        common_cols <- setdiff(intersect(names(results), names(metadata)), "Sample")
        if (length(common_cols) > 0) metadata <- metadata[, !names(metadata) %in% common_cols, drop = FALSE]
        results <- dplyr::left_join(results, metadata, by = "Sample")
        results <- results[match(common_samples, results$Sample), ]
    }
    attr(results, "signatures") <- saved_sigs

    if (is.null(genes)) {
        if (!is.null(attr(results, "signatures"))) {
            genes <- unique(unlist(attr(results, "signatures")))
        } else {
            genes <- DEFAULT_ALL_INPUT_GENES
        }
    }

    genes <- intersect(genes, rownames(expr_mat))
    if (length(genes) < 2) {
        warning("Fewer than 2 genes available for heatmap.")
        return(NULL)
    }

    class_order <- LAB_CLASS_ORDER[LAB_CLASS_ORDER %in% results$Pred_Class]
    results$Pred_Class <- factor(results$Pred_Class, levels = class_order)
    ord <- order(results$Pred_Class)

    row_gaps <- NULL
    annot_row <- NULL
    sigs <- attr(results, "signatures")
    if (is.null(sigs)) sigs <- default_gene_sets()

    if (cluster_rows_by_class) {
        sig_order <- c("chondroid", "squamous", "basal_like", "luminal", "apocrine", "malo_her2")
        sig_label_map <- c(
            "chondroid"  = "Chondroid",
            "squamous"   = "Squamous",
            "basal_like" = "Basal-like",
            "luminal"    = "LHS",
            "apocrine"   = "Apocrine",
            "malo_her2"  = "MA-Lo (HER2)"
        )
        gene_sig_map <- data.frame(gene = character(), signature = character(), sig_label = character(), stringsAsFactors = FALSE)
        for (sig_name in names(sigs)) {
            sig_genes <- intersect(sigs[[sig_name]], genes)
            if (length(sig_genes) > 0) {
                s_name <- tolower(sig_name)
                s_lbl <- if (s_name %in% names(sig_label_map)) sig_label_map[[s_name]] else sig_name
                gene_sig_map <- rbind(gene_sig_map, data.frame(gene = sig_genes, signature = s_name, sig_label = s_lbl, stringsAsFactors = FALSE))
            }
        }
        gene_sig_map <- gene_sig_map[!duplicated(gene_sig_map$gene), ]
        missing_genes <- setdiff(genes, gene_sig_map$gene)
        if (length(missing_genes) > 0) {
            gene_sig_map <- rbind(gene_sig_map, data.frame(gene = missing_genes, signature = "other", sig_label = "Other", stringsAsFactors = FALSE))
        }
        sig_order_full <- c(sig_order, "other")
        gene_sig_map$signature <- factor(gene_sig_map$signature, levels = sig_order_full)
        gene_sig_map <- gene_sig_map[order(gene_sig_map$signature), ]
        genes <- gene_sig_map$gene
        sig_counts <- table(gene_sig_map$signature)
        sig_counts <- sig_counts[sig_counts > 0]
        row_gaps <- cumsum(sig_counts[-length(sig_counts)])

        annot_row <- data.frame(
            Signature = factor(gene_sig_map$sig_label, levels = c("Chondroid", "Squamous", "Basal-like", "LHS", "Apocrine", "MA-Lo (HER2)", "Other")),
            row.names = gene_sig_map$gene
        )
    }

    heat_mat <- as.matrix(expr_mat[genes, ord])
    annot_df <- data.frame(Pred_Class = results$Pred_Class[ord], row.names = common_samples[ord])
    if ("PAM50" %in% colnames(results)) annot_df$PAM50 <- results$PAM50[ord]

    ann_colors <- list(Pred_Class = LAB_CLASS_COLORS[names(LAB_CLASS_COLORS) %in% unique(annot_df$Pred_Class)])
    if ("PAM50" %in% colnames(annot_df)) ann_colors$PAM50 <- PAM50_COLORS[names(PAM50_COLORS) %in% unique(annot_df$PAM50)]
    if (!is.null(annot_row)) {
        ann_colors$Signature <- LAB_CLASS_COLORS[names(LAB_CLASS_COLORS) %in% unique(annot_row$Signature)]
    }

    class_counts <- table(results$Pred_Class[ord])
    gaps <- cumsum(class_counts[-length(class_counts)])

    if (scale_rows) {
        mat <- heat_mat
        for (i in 1:3) {
            mat <- t(scale(t(mat), center = TRUE, scale = FALSE))
            mat <- scale(mat, center = TRUE, scale = FALSE)
        }
        heat_display <- mat
        lim <- quantile(abs(heat_display), 0.99, na.rm = TRUE)
        if (is.na(lim) || lim == 0) lim <- 1
        breaks <- seq(-lim, lim, length.out = 101)
        colors <- colorRampPalette(c("green3", "black", "red3"))(100)
        title_suffix <- "(Centered)"
        file_suffix <- "_Heatmap_Centered.pdf"
    } else {
        heat_display <- heat_mat
        q_breaks <- quantile(heat_mat, probs = seq(0, 1, length.out = 101), na.rm = TRUE)
        breaks <- unique(q_breaks)
        if (length(breaks) < 10) breaks <- seq(min(heat_mat, na.rm = TRUE), max(heat_mat, na.rm = TRUE), length.out = 101)
        colors <- colorRampPalette(c("white", "firebrick3"))(length(breaks) - 1)
        title_suffix <- "(Log2 TPM)"
        file_suffix <- "_Heatmap_TPM.pdf"
    }

    if (cluster_cols) {
        title_suffix <- paste(title_suffix, "- Unsupervised")
        file_suffix <- sub(".pdf$", "_Unsupervised.pdf", file_suffix)
    }

    pdf_height <- max(8, nrow(heat_mat) * 0.15)
    pdf_width <- max(10, ncol(heat_mat) * 0.06)
    do_cluster_rows <- if (cluster_rows_by_class && !is.null(row_gaps)) FALSE else TRUE

    ph <- pheatmap::pheatmap(
        heat_display,
        annotation_col = annot_df,
        annotation_row = annot_row,
        annotation_colors = ann_colors,
        show_colnames = show_names,
        show_rownames = TRUE,
        cluster_rows = do_cluster_rows,
        cluster_cols = cluster_cols,
        treeheight_row = if (do_cluster_rows) 50 else 0,
        treeheight_col = if (cluster_cols) 50 else 0,
        gaps_col = if (cluster_cols) NULL else gaps,
        gaps_row = if (!do_cluster_rows && !is.null(row_gaps)) row_gaps else NULL,
        color = colors,
        breaks = breaks,
        border_color = "grey30",
        fontsize_row = 8,
        main = paste("LAMBclassifier Expression Heatmap", title_suffix),
        silent = TRUE,
        width = pdf_width,
        height = pdf_height
    )

    pdf(file = paste0(output_prefix, file_suffix), width = pdf_width, height = pdf_height)
    grid::grid.draw(ph$gtable)
    dev.off()

    png_suffix <- sub(".pdf$", ".png", file_suffix)
    png(filename = paste0(output_prefix, png_suffix), width = pdf_width, height = pdf_height, units = "in", res = 300)
    grid::grid.draw(ph$gtable)
    dev.off()

    message("Saved: ", output_prefix, file_suffix, " and ", png_suffix)
    invisible(NULL)
}

#' Plot Expression Scatter
#'
#' @description Creates a 2D scatter plot of two genes, colored by LAMB class.
#' @param norm_data Expression matrix.
#' @param results Data frame from LAMBclassifier().
#' @param g1 X-axis gene symbol.
#' @param g2 Y-axis gene symbol.
#' @param output_prefix File path prefix.
#' @export
plot_expression_scatter <- function(norm_data, results, g1, g2, output_prefix) {
    if (!g1 %in% rownames(norm_data) || !g2 %in% rownames(norm_data)) {
        warning("Gene not found in expression data.")
        return(NULL)
    }

    plot_data <- tibble::tibble(
        Sample = colnames(norm_data),
        X_Expr = as.numeric(norm_data[g1, ]),
        Y_Expr = as.numeric(norm_data[g2, ])
    )
    plot_data <- dplyr::left_join(plot_data, results[, c("Sample", "Pred_Class")], by = "Sample")

    p <- ggplot(plot_data, aes(x = Y_Expr, y = X_Expr, color = Pred_Class)) +
        geom_point(size = 1.8, alpha = 0.8) +
        labs(title = paste(g1, "vs", g2), x = paste(g2, "(Log2 TPM)"), y = paste(g1, "(Log2 TPM)"), color = "Class") +
        scale_color_manual(values = LAB_CLASS_COLORS) +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5), aspect.ratio = 1)

    filename <- paste0(output_prefix, "_Scatter_", g1, "_vs_", g2, ".pdf")
    ggsave(filename, p, width = 6, height = 5)
    ggsave(sub(".pdf$", ".png", filename), p, width = 6, height = 5, dpi = 300)
    message("Saved Scatter Plot: ", filename)
}

#' Plot Prediction--Reference Cross-classification
#'
#' @description Generates a cross-classification heatmap between LAMB
#'   predictions and external labels. The reference labels are not assumed to
#'   be ground truth or to map one-to-one to LAMB classes.
#' @param results Data frame containing prediction and reference columns.
#' @param truth_col Name of the column containing external reference labels.
#' @param output_prefix File path prefix.
#' @param title_suffix Optional title suffix.
#' @export
#' @import ggplot2
plot_confusion_matrix <- function(results, truth_col, output_prefix, title_suffix = "") {
    if (!truth_col %in% colnames(results)) return(NULL)
    dat <- results[!is.na(results$Pred_Class) & !is.na(results[[truth_col]]), ]
    if (nrow(dat) == 0) return(NULL)

    curr_df <- as.data.frame(table(Prediction = dat$Pred_Class, Reference = dat[[truth_col]]))
    colnames(curr_df) <- c("Prediction", "Reference", "Count")

    p <- ggplot(curr_df, aes(x = Reference, y = Prediction, fill = Count)) +
        geom_tile(color = "white") +
        geom_text(aes(label = Count), color = "black", size = 4) +
        scale_fill_gradient(low = "white", high = "steelblue") +
        labs(title = paste("Prediction--reference cross-classification:", title_suffix), subtitle = paste("N =", sum(curr_df$Count)),
             x = paste("Reference (", truth_col, ")", sep = ""), y = "LAMB Prediction") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", hjust = 0.5))

    filename <- paste0(output_prefix, "_Confusion_", truth_col, ".pdf")
    ggsave(filename, p, width = 6, height = 5)
    ggsave(sub(".pdf$", ".png", filename), p, width = 6, height = 5, dpi = 300)
    message("Saved prediction--reference cross-classification: ", filename)
}

#' Generate Standard Heatmaps
#'
#' @description Generates Supervised TPM, Supervised Centered, and Unsupervised Centered heatmaps.
#' @param expr_mat Expression matrix (log2 TPM).
#' @param results Data frame from LAMBclassifier().
#' @param metadata Optional metadata data frame.
#' @param genes Optional gene vector.
#' @param output_prefix Output file path prefix.
#' @param lab_column Optional metadata comparison column name.
#' @param show_names Logical, whether to show sample names.
#' @param add_top_n_markers Reserved compatibility argument. Additional markers
#'   are not added by the locked implementation.
#' @export
generate_standard_heatmaps <- function(expr_mat, results, metadata = NULL, genes = NULL, output_prefix, lab_column = NULL, show_names = FALSE, add_top_n_markers = 0) {
    # 1. Supervised TPM
    plot_comprehensive_heatmap(
        expr_mat = expr_mat, results = results, metadata = metadata,
        genes = genes, output_prefix = output_prefix, scale_rows = FALSE,
        cluster_cols = FALSE, lab_column = lab_column, show_names = show_names
    )

    # 2. Supervised Centered (Signature Grouped)
    plot_comprehensive_heatmap(
        expr_mat = expr_mat, results = results, metadata = metadata,
        genes = genes, output_prefix = output_prefix, scale_rows = TRUE,
        cluster_cols = FALSE, cluster_rows_by_class = TRUE, lab_column = lab_column, show_names = show_names
    )

    # 3. Unsupervised Centered
    plot_comprehensive_heatmap(
        expr_mat = expr_mat, results = results, metadata = metadata,
        genes = genes, output_prefix = output_prefix, scale_rows = TRUE,
        cluster_cols = TRUE, lab_column = lab_column, show_names = show_names
    )
}

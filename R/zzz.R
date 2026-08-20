# =============================================================================
# LAMBclassifier: Package-Level Declarations
# =============================================================================

#' LAMBclassifier: Luminal-Apocrine-Metaplastic-Basal Breast Cancer Classifier
#'
#' A breast cancer sample classifier based on lineage-aware mammary breast
#' classification. Evaluates single-sample rank enrichment across six developmental
#' endpoints.
#'
#' @docType package
#' @name LAMBclassifier-package
#' @keywords internal
#' @importFrom stats median quantile model.matrix sd var predict density runif
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom utils head read.table write.csv
"_PACKAGE"

# Suppress R CMD check notes for NSE variables used in dplyr/ggplot2
# These are column names in data frames accessed via non-standard evaluation
utils::globalVariables(c(
    # dplyr/tidyr column references
    "Parameters", "Class", "Jaccard", "Freq", "Total", "Pct",
    "PAM50", "LAB", "Pred_Class", "Sample",

    # ggplot2 aesthetics
    "Score_Squamous", "Score_Chondroid",
    "Score_MALO_HER2", "Score_LHS", "Score_Sensory", "Score_Basal_Like",
    "Axis_Metaplastic", "Axis_Basal_Like", "Axis_Luminal",

    # Lineage tree coordinates
    "x", "y", "xend", "yend", "node_x", "node_y",
    "x_start", "y_start", "x_end", "y_end",
    "x_start_clean", "x_end_clean", "xmin", "xmax", "ymin", "ymax",
    "border_col", "text_col", "stat", "plot_score_lhs", "title",
    "nudge_val", "full_label", "hjust_val",
    "display_label", "fill_col", "nudge_y", "v_just",
    "Y_Expr", "X_Expr", "scaled", "Density",
    "reference", "prediction", "Reference", "Prediction", "Count", "parent",
    "id", "label", "hjust", "plot_score_apo",

    # Internal package data
    "gene.length", "pam50.robust",

    # Pivot operations
    "."
))

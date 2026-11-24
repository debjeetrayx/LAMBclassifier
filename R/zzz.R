#' LAMBclassifier: Luminal-Apocrine-Metaplastic-Basal Breast Cancer Classifier
#'
#' Assigns breast cancer expression profiles to six expression-defined classes
#' (LHS, Apocrine, MA-Lo, Basal-like, Chondroid, Squamous) with single-sample
#' singscore rank enrichment and a fixed decision tree, and reports an
#' ERBB2-amplicon read-out beside every class. Research use only. Start with
#' \code{\link{lamb_run}}.
#'
#' @docType package
#' @name LAMBclassifier-package
#' @keywords internal
#' @importFrom stats median quantile model.matrix sd var predict density runif
#' @importFrom grDevices colorRampPalette dev.off pdf png
#' @importFrom utils head read.table write.csv
"_PACKAGE"

# Columns used through non-standard evaluation (dplyr, ggplot2).
utils::globalVariables(c(
    "Parameters", "Class", "Jaccard", "Freq", "Total", "Pct",
    "PAM50", "LAB", "Pred_Class", "Sample",

    "Score_Squamous", "Score_Chondroid",
    "Score_MALO", "Score_LHS", "Score_Basal_Like",
    "Axis_Metaplastic", "Axis_Basal_Like", "Axis_Luminal",

    "x", "y", "xend", "yend", "node_x", "node_y",
    "x_start", "y_start", "x_end", "y_end",
    "x_start_clean", "x_end_clean", "xmin", "xmax", "ymin", "ymax",
    "border_col", "text_col", "stat", "plot_score_lhs", "title",
    "nudge_val", "full_label", "hjust_val",
    "display_label", "fill_col", "nudge_y", "v_just",
    "Y_Expr", "X_Expr", "scaled", "Density",
    "reference", "prediction", "Reference", "Prediction", "Count", "parent",
    "id", "label", "hjust", "plot_score_apo",

    "gene.length", "pam50.robust",

    "."
))

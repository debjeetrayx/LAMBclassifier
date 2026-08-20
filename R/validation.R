# Reproducibility and validation helpers for a locked LAMB specification.

#' Return the Locked LAMB Classifier Specification
#'
#' Produces deterministic SHA-256 identifiers for the ordered signature genes
#' and reference universe. These identifiers should accompany every exported
#' result so that an analysis cannot silently mix classifier versions.
#'
#' @param gene_sets Named list of signature gene vectors.
#' @param calibration Calibration list returned by [lamb_calibration()].
#' @return A list containing the classifier version, hashes, signatures,
#'   calibration constants, and universe size.
#' @export
lamb_specification <- function(gene_sets = default_gene_sets(),
                               calibration = lamb_calibration()) {
    if (is.null(names(gene_sets)) || any(!nzchar(names(gene_sets)))) {
        stop("`gene_sets` must be a named list.")
    }

    signature_text <- paste(
        vapply(names(gene_sets), function(nm) {
            paste0(nm, "=", paste(as.character(gene_sets[[nm]]), collapse = ","))
        }, character(1)),
        collapse = "\n"
    )
    universe <- calibration$universe

    list(
        version = as.character(calibration$version),
        signature_sha256 = digest::digest(signature_text, algo = "sha256", serialize = FALSE),
        universe_sha256 = digest::digest(paste(universe, collapse = "\n"),
                                         algo = "sha256", serialize = FALSE),
        universe_size = length(universe),
        signatures = gene_sets,
        min_margin = calibration$min_margin,
        min_gene_coverage = calibration$min_gene_coverage
    )
}


.wilson_interval <- function(x, n, conf_level = 0.95) {
    if (!is.finite(n) || n <= 0) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    p <- x / n
    denominator <- 1 + z^2 / n
    centre <- (p + z^2 / (2 * n)) / denominator
    half_width <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator
    c(estimate = p, lower = centre - half_width, upper = centre + half_width)
}


#' Binary Performance with Wilson Confidence Intervals
#'
#' Computes a complete confusion matrix and standard binary metrics after
#' removing rows with missing predictions or reference labels. The positive
#' prediction and reference values must be declared explicitly.
#'
#' @param predicted Predicted labels.
#' @param reference Independent reference labels.
#' @param positive_predicted Value(s) counted as a positive prediction.
#' @param positive_reference Value(s) counted as a positive reference.
#' @param conf_level Confidence level for Wilson intervals.
#' @return A one-row data frame with counts and metrics.
#' @export
lamb_binary_metrics <- function(predicted, reference,
                                positive_predicted,
                                positive_reference,
                                conf_level = 0.95) {
    if (length(predicted) != length(reference)) {
        stop("`predicted` and `reference` must have equal length.")
    }
    keep <- !is.na(predicted) & !is.na(reference)
    pred <- predicted[keep] %in% positive_predicted
    ref <- reference[keep] %in% positive_reference

    tp <- sum(pred & ref)
    fp <- sum(pred & !ref)
    fn <- sum(!pred & ref)
    tn <- sum(!pred & !ref)

    sensitivity <- .wilson_interval(tp, tp + fn, conf_level)
    specificity <- .wilson_interval(tn, tn + fp, conf_level)
    ppv <- .wilson_interval(tp, tp + fp, conf_level)
    npv <- .wilson_interval(tn, tn + fn, conf_level)

    data.frame(
        N = length(ref), TP = tp, FP = fp, FN = fn, TN = tn,
        Sensitivity = sensitivity["estimate"], Sensitivity_Lower = sensitivity["lower"],
        Sensitivity_Upper = sensitivity["upper"],
        Specificity = specificity["estimate"], Specificity_Lower = specificity["lower"],
        Specificity_Upper = specificity["upper"],
        PPV = ppv["estimate"], PPV_Lower = ppv["lower"], PPV_Upper = ppv["upper"],
        NPV = npv["estimate"], NPV_Lower = npv["lower"], NPV_Upper = npv["upper"],
        Balanced_Accuracy = mean(c(sensitivity["estimate"], specificity["estimate"]), na.rm = TRUE),
        Conf_Level = conf_level,
        row.names = NULL,
        check.names = FALSE
    )
}


#' Cross-tabulate LAMB and an External Reference
#'
#' This function describes association between taxonomies. It does not assume
#' that the external reference is ground truth or that its categories have a
#' one-to-one mapping to LAMB.
#'
#' @param predicted LAMB labels.
#' @param reference External labels.
#' @param normalize Return counts, row proportions (`"prediction"`), or column
#'   proportions (`"reference"`).
#' @return A matrix.
#' @export
lamb_cross_tab <- function(predicted, reference,
                           normalize = c("none", "prediction", "reference")) {
    normalize <- match.arg(normalize)
    keep <- !is.na(predicted) & !is.na(reference)
    tab <- table(Prediction = predicted[keep], Reference = reference[keep])
    if (normalize == "prediction") return(prop.table(tab, margin = 1))
    if (normalize == "reference") return(prop.table(tab, margin = 2))
    tab
}


#' Cramer's V for Cross-classification Association
#'
#' @param predicted Predicted labels.
#' @param reference External labels.
#' @return Cramer's V, or `NA` when either variable has fewer than two levels.
#' @export
lamb_cramers_v <- function(predicted, reference) {
    tab <- lamb_cross_tab(predicted, reference, normalize = "none")
    if (nrow(tab) < 2L || ncol(tab) < 2L || sum(tab) == 0L) return(NA_real_)
    chi <- suppressWarnings(stats::chisq.test(tab, correct = FALSE)$statistic)
    unname(sqrt(chi / (sum(tab) * min(nrow(tab) - 1L, ncol(tab) - 1L))))
}

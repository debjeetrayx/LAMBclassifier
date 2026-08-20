#' @title Reference Gene Universe
#' @description The fixed gene space in which all signature scores are computed.
#' @return Character vector of HGNC symbols.
#' @export
lamb_universe <- function() LAMB_CALIBRATION$universe

#' @title Calibration Constants
#' @description Decision thresholds used by \code{\link{LAMBclassifier}}.
#' @return A list of calibration constants.
#' @export
lamb_calibration <- function() LAMB_CALIBRATION


#' @title Score LAMB Signatures
#'
#' @description
#' Computes a singscore rank-enrichment score for each signature, independently
#' for every sample.
#'
#' @details
#' Input must be log2-transformed. Genes outside the reference universe are
#' dropped before scoring, which keeps scores comparable across datasets.
#'
#' @param data Gene x sample matrix, rows named with HGNC gene symbols,
#'   log2-transformed.
#' @param gene_sets Named list of signature gene vectors.
#' @param min_genes_per_sig Minimum number of present genes required to score a
#'   signature; otherwise that signature is assigned `NA`.
#' @param check_coverage Warn when a signature has too few genes present.
#'
#' @return A samples x signatures matrix, with a \code{"coverage"} attribute
#'   giving genes requested, found and the fraction per signature, and a
#'   \code{"universe_coverage"} attribute.
#'
#' @importFrom singscore rankGenes simpleScore
#' @export
lamb_score <- function(data,
                       gene_sets = default_gene_sets(),
                       min_genes_per_sig = 3,
                       check_coverage = TRUE) {
    if (!is.matrix(data)) data <- as.matrix(data)
    if (is.null(rownames(data))) stop("`data` must have gene symbols as rownames.")
    if (ncol(data) < 1L) stop("`data` has no samples.")
    if (nrow(data) < 1L) stop("`data` has no genes.")
    if (!is.numeric(data)) stop("`data` must be numeric; found ", typeof(data), ".")
    if (all(is.na(data))) stop("`data` contains no non-missing values.")

    dup <- unique(rownames(data)[duplicated(rownames(data))])
    if (length(dup)) {
        warning("`data` has ", length(dup), " duplicated gene symbol(s), e.g. ",
                paste(utils::head(dup, 3), collapse = ", "),
                ". Only the first occurrence of each is scored; ",
                "collapse duplicates before scoring to control which row is used.",
                call. = FALSE)
        data <- data[!duplicated(rownames(data)), , drop = FALSE]
    }

    if (max(data, na.rm = TRUE) > 100) {
        warning("Values above 100 detected - `data` looks untransformed. ",
                "Applying log2(x + 1).", call. = FALSE)
        data <- log2(data + 1)
    }

    universe <- lamb_universe()
    keep <- intersect(universe, rownames(data))

    # Scoring happens strictly inside the locked universe. Falling back to
    # arbitrary supplied genes would change the ranking space while the output
    # still carried the locked universe hash, so it is refused instead.
    if (length(keep) < 100L) {
        stop("Only ", length(keep), " of the ", length(universe),
             " reference-universe genes are present. LAMBclassifier scores ",
             "within the locked universe and will not score outside it. ",
             "At least 100 universe genes are required.", call. = FALSE)
    }
    if (length(keep) < 1000L) {
        message("Targeted panel: ", length(keep),
                " reference-universe genes present. Check per-signature ",
                "coverage before interpreting calls.")
    }

    mat <- data[keep, , drop = FALSE]
    ranked <- singscore::rankGenes(mat)

    cov_rows <- list()
    score_cols <- list()
    for (sig in names(gene_sets)) {
        requested <- unique(gene_sets[[sig]])
        found <- intersect(requested, rownames(mat))

        cov_rows[[sig]] <- data.frame(
            Signature   = sig,
            N_Requested = length(requested),
            N_Found     = length(found),
            Coverage    = round(length(found) / length(requested), 3),
            stringsAsFactors = FALSE
        )

        if (length(found) < min_genes_per_sig) {
            warning("Signature '", sig, "' has only ", length(found), " gene(s) present ",
                    "(minimum required: ", min_genes_per_sig, "). Deactivating signature (score = NA).",
                    call. = FALSE)
            score_cols[[sig]] <- rep(NA_real_, ncol(mat))
            next
        }
        score_cols[[sig]] <- singscore::simpleScore(ranked, upSet = found)$TotalScore
    }

    scores <- do.call(cbind, score_cols)
    rownames(scores) <- colnames(mat)
    colnames(scores) <- names(gene_sets)

    coverage <- do.call(rbind, cov_rows)
    rownames(coverage) <- NULL

    if (check_coverage) {
        low <- coverage$Signature[coverage$Coverage < LAMB_CALIBRATION$min_gene_coverage]
        if (length(low)) {
            warning("Low gene coverage for: ", paste(low, collapse = ", "),
                    ". Affected calls will be flagged.", call. = FALSE)
        }
    }

    attr(scores, "coverage") <- coverage
    attr(scores, "universe_coverage") <- length(keep) / length(universe)
    scores
}


#' @title Default LAMB Gene Sets
#' @description The six signatures used by \code{\link{LAMBclassifier}}.
#' @return Named list of character vectors.
#' @export
default_gene_sets <- function() {
    list(
        chondroid  = DEFAULT_CHONDROID_GENES,
        squamous   = DEFAULT_SQUAMOUS_EPITHELIAL_GENES,
        basal_like = DEFAULT_BASAL_LIKE_GENES,
        luminal    = DEFAULT_LHS_GENES,
        apocrine   = DEFAULT_APOCRINE_GENES,
        malo_her2  = DEFAULT_MALO_HER2_GENES
    )
}

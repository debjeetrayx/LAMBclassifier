#' Luminal-Apocrine-Metaplastic-Basal Classifier (Bulk, Single-Cell & Spatial)
#'
#' @description
#' Assigns each sample to one of six breast cancer subtypes: LHS,
#' MA-Lo (HER2), Apocrine, Basal-like, Chondroid or Squamous.
#'
#' @details
#' Samples are classified via a binary lineage tree followed by a separate
#' ERBB2-amplicon expression check:
#' \enumerate{
#'   \item Metaplastic (Chondroid vs Squamous) vs committed mammary epithelium.
#'   \item Step 2: LASP Anchor (Luminal Progenitor) binary division into Basal-like vs Luminal Lineage.
#'   \item Among luminal samples, binary decision between LHS vs Apocrine lineage.
#' }
#' Luminal samples are separately evaluated for an ERBB2 17q12 amplicon
#' expression signal (\code{MA-Lo (HER2)}). This expression label is not a
#' direct DNA-amplification assay or clinical HER2 diagnosis.
#'
#' \code{Flag} reports quality concerns and does not affect the class:
#' \describe{
#'   \item{close_call}{Two subtypes scored almost equally (margin < min_margin).}
#'   \item{low_coverage}{Too few signature genes present to score reliably.}
#' }
#'
#' @param data Gene x sample matrix, rows named with HGNC gene symbols,
#'   log2-transformed.
#' @param chondroid_genes,squamous_genes,basal_like_genes Gene vectors for the
#'   corresponding signatures.
#' @param lhs_genes,apocrine_genes,malo_her2_genes Gene vectors for the
#'   corresponding signatures.
#' @param min_genes_per_sig Minimum required present genes for a signature to be evaluated (default: 3).
#'   If a signature has fewer present genes in the input dataset, its score is set to NA and its branch check is skipped.
#' @param run_pam50 Logical. Also run PAM50 subtyping and append a
#'   \code{PAM50} column.
#' @param pam50_data PAM50 centroids. \code{NULL} uses the bundled centroids.
#' @param gene_length_info Gene lengths for PAM50 mapping. \code{NULL} uses the
#'   bundled data.
#' @param calibration Calibration constants. \code{NULL} uses the bundled set;
#'   supply your own only for a platform needing its own fit, such as a
#'   targeted panel.
#'
#' @return A tibble with one row per sample: \code{Sample}, \code{Pred_Class},
#'   \code{Margin}, \code{Flag}, the metaplastic score, and
#'   the six signature scores. Signatures and per-signature gene coverage are
#'   attached as attributes.
#'
#' @export
#' @importFrom singscore rankGenes simpleScore
#' @importFrom genefu molecular.subtyping
LAMBclassifier <- function(data,
                           chondroid_genes   = DEFAULT_CHONDROID_GENES,
                           squamous_genes    = DEFAULT_SQUAMOUS_EPITHELIAL_GENES,
                           basal_like_genes  = DEFAULT_BASAL_LIKE_GENES,
                           lhs_genes         = DEFAULT_LHS_GENES,
                           apocrine_genes    = DEFAULT_APOCRINE_GENES,
                           malo_her2_genes   = DEFAULT_MALO_HER2_GENES,
                           min_genes_per_sig = 3,
                           run_pam50         = TRUE,
                           pam50_data        = NULL,
                           gene_length_info  = NULL,
                           calibration       = NULL) {

    if (!is.matrix(data)) data <- as.matrix(data)
    if (run_pam50 && is.null(gene_length_info)) gene_length_info <- get_gene_length()
    if (run_pam50 && is.null(pam50_data)) pam50_data <- get_pam50_data()
    cal <- calibration %||% lamb_calibration()

    gene_sets <- list(
        chondroid  = chondroid_genes  %||% DEFAULT_CHONDROID_GENES,
        squamous   = squamous_genes   %||% DEFAULT_SQUAMOUS_EPITHELIAL_GENES,
        basal_like = basal_like_genes %||% DEFAULT_BASAL_LIKE_GENES,
        luminal    = lhs_genes        %||% DEFAULT_LHS_GENES,
        apocrine   = apocrine_genes   %||% DEFAULT_APOCRINE_GENES,
        malo_her2  = malo_her2_genes  %||% DEFAULT_MALO_HER2_GENES
    )
    specification <- lamb_specification(gene_sets = gene_sets, calibration = cal)

    scores   <- lamb_score(data, gene_sets = gene_sets, min_genes_per_sig = min_genes_per_sig)
    coverage <- attr(scores, "coverage")
    low_cov  <- coverage$Signature[coverage$Coverage < cal$min_gene_coverage]

    predictions <- rep(NA_character_, nrow(scores))
    margins     <- rep(NA_real_, nrow(scores))
    flags       <- rep("", nrow(scores))

    # ---- Gate evaluability --------------------------------------------------
    # A call names the class a sample belongs to *and* implies every exclusion
    # made on the way there. An unscored signature makes the corresponding gate
    # unevaluable, so the sample cannot be assigned rather than being allowed to
    # fall through as though the gate had returned "no".
    required <- c("chondroid", "squamous", "basal_like", "luminal", "apocrine")
    unscored <- vapply(
        required,
        function(sig) is.na(scores[, sig]),
        logical(nrow(scores))
    )
    if (is.null(dim(unscored))) {
        unscored <- matrix(unscored, nrow = nrow(scores),
                           dimnames = list(NULL, required))
    }
    unevaluable <- apply(unscored, 1L, any)
    missing_sigs <- apply(unscored, 1L, function(z) paste(required[z], collapse = ","))

    # ---- Step 1: metaplastic vs mammary (Uniform Threshold-Free Comparison) --
    meta_raw    <- pmax(scores[, "chondroid"], scores[, "squamous"])
    mam_signal  <- pmax(scores[, "basal_like"], scores[, "luminal"], scores[, "apocrine"])
    meta_signal <- meta_raw - mam_signal

    is_meta     <- !is.na(meta_raw) & !is.na(mam_signal) & (meta_raw > mam_signal)

    if (any(is_meta)) {
        chondroid_wins <- is_meta & (scores[, "chondroid"] > scores[, "squamous"])
        predictions[chondroid_wins] <- "Chondroid"
        predictions[is_meta & !chondroid_wins] <- "Squamous"
        margins[is_meta] <- meta_raw[is_meta] - mam_signal[is_meta]
    }

    # ---- Step 2: Basal-like Dual-Winner Check (Beats LHS and Apocrine) -------
    in_mammary  <- !is_meta
    is_basal    <- in_mammary & !is.na(scores[, "basal_like"]) &
        !is.na(scores[, "luminal"]) & !is.na(scores[, "apocrine"]) &
        (scores[, "basal_like"] > scores[, "luminal"]) &
        (scores[, "basal_like"] > scores[, "apocrine"])

    if (any(is_basal)) {
        predictions[is_basal] <- "Basal-like"
        max_lum_apo <- pmax(scores[is_basal, "luminal"], scores[is_basal, "apocrine"])
        margins[is_basal] <- scores[is_basal, "basal_like"] - max_lum_apo
    }

    # ---- Step 3: LHS vs Apocrine & ERBB2-amplicon expression check -----------
    sig_of     <- cal$signature_of_class
    in_luminal <- in_mammary & !is_basal

    if (any(in_luminal)) {
        idx  <- which(in_luminal)
        lum  <- scores[idx, "luminal"]
        apo  <- scores[idx, "apocrine"]
        malo <- scores[idx, "malo_her2"]

        usable <- !is.na(lum) & !is.na(apo)
        if (any(usable)) {
            # 1. Binary lineage decision: LHS vs Apocrine
            lum_wins <- usable & (lum > apo)
            predictions[idx[lum_wins]] <- "LHS"
            margins[idx[lum_wins]] <- lum[lum_wins] - apo[lum_wins]

            apo_wins <- usable & !lum_wins
            predictions[idx[apo_wins]] <- "Apocrine"
            margins[idx[apo_wins]] <- apo[apo_wins] - lum[apo_wins]

            # 2. Separate ERBB2-amplicon expression check (MA-Lo (HER2))
            max_lum_lineage <- pmax(lum, apo)
            is_her2_amp <- usable & !is.na(malo) & (malo > max_lum_lineage)
            if (any(is_her2_amp)) {
                predictions[idx[is_her2_amp]] <- "MA-Lo (HER2)"
                margins[idx[is_her2_amp]] <- malo[is_her2_amp] - max_lum_lineage[is_her2_amp]
            }
        }
        if (any(!usable)) predictions[idx[!usable]] <- "NC"
    }

    # Unevaluable samples are NC regardless of what the reachable gates produced.
    predictions[unevaluable] <- "NC"
    margins[unevaluable] <- NA_real_

    predictions[is.na(predictions)] <- "NC"

    # ---- Quality flags -------------------------------------------------------
    add_flag <- function(f, which, tag) {
        f[which] <- ifelse(nzchar(f[which]), paste0(f[which], ";", tag), tag)
        f
    }

    borderline <- !is.na(margins) &
        predictions != "NC" &
        margins < cal$min_margin
    flags <- add_flag(flags, borderline, "close_call")

    # A call depends on every signature on its decision path. Flag when any of
    # them is below the coverage threshold.
    if (length(low_cov)) {
        decision_relevant <- union(required, "malo_her2")
        any_low <- any(decision_relevant %in% low_cov)
        flags <- add_flag(flags, rep(any_low, length(flags)) & predictions != "NC",
                          "low_coverage")
    }

    # Record why a sample could not be assigned.
    flags <- add_flag(flags, unevaluable, "incomplete_gates")

    results <- tibble::tibble(
        Sample     = rownames(scores),
        Pred_Class = predictions,
        Margin     = round(margins, 4),
        Flag       = flags,
        Unscored_Signatures = missing_sigs,
        Classifier_Version = specification$version,
        Signature_SHA256   = specification$signature_sha256,
        Universe_SHA256    = specification$universe_sha256,

        Metaplastic_Score  = round(meta_signal, 4),

        Score_Chondroid  = round(scores[, "chondroid"], 4),
        Score_Squamous   = round(scores[, "squamous"], 4),
        Score_Basal_Like = round(scores[, "basal_like"], 4),
        Score_LHS        = round(scores[, "luminal"], 4),
        Score_Apocrine   = round(scores[, "apocrine"], 4),
        Score_MALO_HER2  = round(scores[, "malo_her2"], 4)
    )

    if (run_pam50 && !is.null(pam50_data)) {
        message("Running PAM50 classification...")
        pam50_res <- run_pam50_internal(data, pam50_data, gene_length_info)
        results <- dplyr::left_join(results, pam50_res, by = "Sample")
    }

    # set after any join: dplyr drops custom attributes
    attr(results, "signatures") <- gene_sets
    attr(results, "coverage") <- coverage
    attr(results, "calibration_version") <- cal$version
    attr(results, "specification") <- specification

    results
}


#' Run PAM50 Intrinsic Subtyping (Internal)
#' @keywords internal
#' @importFrom genefu molecular.subtyping
run_pam50_internal <- function(data, pam50_centroids, gene_length_info) {
    annots <- gene_length_info[
        gene_length_info$entrezid %in%
            pam50_centroids$centroids.map$EntrezGene.ID,
        c("SYMBOL", "entrezid")
    ]

    data_tbl <- tibble::as_tibble(data, rownames = "SYMBOL")
    annots <- dplyr::left_join(annots, data_tbl, by = "SYMBOL")
    annots <- unique(annots)
    rownames(annots) <- annots$SYMBOL

    dat <- t(annots[, -c(1:2), drop = FALSE])
    annot_matrix <- annots[, 1:2]
    colnames(annot_matrix) <- c("Gene.Symbol", "EntrezGene.ID")

    output <- genefu::molecular.subtyping(
        sbt.model = "pam50",
        data = dat,
        annot = annot_matrix,
        do.mapping = FALSE
    )

    tibble::tibble(
        Sample = names(output$subtype),
        PAM50 = as.character(output$subtype)
    )
}


#' Get gene length information from internal data
#' @keywords internal
get_gene_length <- function() {
    env <- parent.env(environment())
    if (exists("gene.length", envir = env, inherits = TRUE)) {
        return(get("gene.length", envir = env, inherits = TRUE))
    }
    tryCatch({
        pkg_env <- asNamespace("LAMBclassifier")
        if (exists("gene.length", envir = pkg_env, inherits = FALSE)) {
            return(get("gene.length", envir = pkg_env))
        }
    }, error = function(e) NULL)
    NULL
}

#' Get PAM50 centroids from internal data
#' @keywords internal
get_pam50_data <- function() {
    env <- parent.env(environment())
    if (exists("pam50.robust", envir = env, inherits = TRUE)) {
        return(get("pam50.robust", envir = env, inherits = TRUE))
    }
    tryCatch({
        pkg_env <- asNamespace("LAMBclassifier")
        if (exists("pam50.robust", envir = pkg_env, inherits = FALSE)) {
            return(get("pam50.robust", envir = pkg_env))
        }
    }, error = function(e) NULL)
    NULL
}

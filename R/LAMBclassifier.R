#' Luminal-Apocrine-Metaplastic-Basal Classifier (Bulk, Single-Cell & Spatial)
#'
#' @description
#' Assigns each sample to one of six breast cancer subtypes: LHS,
#' MA-Lo, Apocrine, Basal-like, Chondroid or Squamous.
#'
#' @details
#' Samples are classified via a binary lineage tree (specification 2.0.0):
#' \enumerate{
#'   \item Metaplastic (Chondroid vs Squamous) vs committed mammary epithelium.
#'   \item Basal-like vs luminal lineage.
#'   \item Among luminal samples, LHS vs Apocrine; the sample is called \code{MA-Lo} instead when the MA-Lo
#'     score exceeds both the LHS and the Apocrine scores.
#' }
#'
#' One read-out sits beside the class and never changes it:
#' \describe{
#'   \item{\code{HER2_Amplicon}}{An ERBB2 17q12 amplicon modifier, given for
#'     every sample whatever its class. Within the sample's reference-universe
#'     genes each gene gets a percentile from the top,
#'     \code{100 * (1 - (rank - 1) / (n - 1))} (0 = highest; ties averaged;
#'     missing values excluded from \code{n}). "high" when at least 3 of the 5
#'     module genes (\code{\link{DEFAULT_HER2_AMPLICON_GENES}}) have a
#'     percentile of at most \code{q} = 9; "not assessable" when fewer than 3
#'     module genes are measured; "not high" otherwise. An expression read-out,
#'     not a DNA-amplification assay or a clinical HER2 result.}
#' }
#'
#' \code{Flag} reports quality concerns and does not affect the class:
#' \describe{
#'   \item{close_call}{Two subtypes scored almost equally (margin < min_margin).}
#'   \item{low_coverage}{One of the six class signatures has under
#'     \code{min_gene_coverage} (0.60) of its genes measured in that sample.}
#'   \item{incomplete_gates}{A gate could not be evaluated, so the sample is
#'     \code{NC}; \code{Unscored_Signatures} names the signatures.}
#' }
#'
#' @param data Gene x sample matrix, rows named with HGNC gene symbols,
#'   log2-transformed.
#' @param chondroid_genes,squamous_genes,basal_like_genes Gene vectors for the
#'   corresponding signatures.
#' @param lhs_genes,apocrine_genes,malo_genes Gene vectors for the
#'   corresponding signatures.
#' @param malo_her2_genes Deprecated (specification 1.0.0 name). When given,
#'   it is used as \code{malo_genes}, with a warning.
#' @param min_genes_per_sig Minimum number of measured genes for a signature to
#'   be scored; below it the score is \code{NA} and dependent calls are NC.
#' @param run_pam50 Logical. Also run PAM50 subtyping and append a
#'   \code{PAM50} column.
#' @param pam50_data PAM50 centroids. \code{NULL} uses the bundled centroids.
#' @param gene_length_info Gene lengths for PAM50 mapping. \code{NULL} uses the
#'   bundled data.
#' @param calibration Calibration constants. \code{NULL} uses the bundled set
#'   (\code{\link{lamb_calibration}()}). Other constants make a different
#'   classifier.
#'
#' @return A tibble with one row per sample: \code{Sample}, \code{Pred_Class},
#'   \code{Margin}, \code{Flag}, \code{Unscored_Signatures},
#'   \code{Classifier_Version}, \code{Signature_SHA256},
#'   \code{Universe_SHA256}, \code{Metaplastic_Score}, the six class signature
#'   scores (\code{Score_Chondroid}, \code{Score_Squamous},
#'   \code{Score_Basal_Like}, \code{Score_LHS}, \code{Score_Apocrine},
#'   \code{Score_MALO}), \code{HER2_Amplicon} and the counts behind it
#'   (\code{HER2_Amplicon_N_High}, \code{HER2_Amplicon_N_Measured}), and
#'   \code{PAM50} when \code{run_pam50 = TRUE}. Signatures, per-signature gene
#'   coverage and the specification are attached as attributes, and:
#'   \describe{
#'     \item{\code{"gene_percentiles"}}{Samples x genes matrix: each gene's
#'       percentile rank among the sample's measured reference-universe genes
#'       (100 = highest, 0 = lowest, ties averaged, 1 decimal), for nine
#'       markers (ESR1, PGR, AR, ERBB2, FOXA1, GATA3, KRT5, TP63, MKI67), every
#'       signature gene and the HER2 amplicon genes. \code{NA} when not
#'       measured.}
#'     \item{\code{"n_genes_ranked"}}{Per sample, the number of measured
#'       universe genes the ranks are taken among.}
#'   }
#'
#' @export
#' @importFrom singscore rankGenes simpleScore
#' @importFrom genefu intrinsic.cluster.predict
LAMBclassifier <- function(data,
                           chondroid_genes   = DEFAULT_CHONDROID_GENES,
                           squamous_genes    = DEFAULT_SQUAMOUS_EPITHELIAL_GENES,
                           basal_like_genes  = DEFAULT_BASAL_LIKE_GENES,
                           lhs_genes         = DEFAULT_LHS_GENES,
                           apocrine_genes    = DEFAULT_APOCRINE_GENES,
                           malo_genes        = DEFAULT_MALO_GENES,
                           min_genes_per_sig = 3,
                           run_pam50         = TRUE,
                           pam50_data        = NULL,
                           gene_length_info  = NULL,
                           calibration       = NULL,
                           malo_her2_genes   = NULL) {

    if (!is.matrix(data)) data <- as.matrix(data)
    if (run_pam50 && is.null(gene_length_info)) gene_length_info <- get_gene_length()
    if (run_pam50 && is.null(pam50_data)) pam50_data <- get_pam50_data()
    cal <- calibration %||% lamb_calibration()
    if (!is.null(malo_her2_genes)) {
        warning("`malo_her2_genes` is deprecated (specification 1.0.0); use `malo_genes`. ",
                "It was used as the MA-Lo signature.", call. = FALSE)
        malo_genes <- malo_her2_genes
    }

    gene_sets <- list(
        chondroid  = chondroid_genes  %||% DEFAULT_CHONDROID_GENES,
        squamous   = squamous_genes   %||% DEFAULT_SQUAMOUS_EPITHELIAL_GENES,
        basal_like = basal_like_genes %||% DEFAULT_BASAL_LIKE_GENES,
        luminal    = lhs_genes        %||% DEFAULT_LHS_GENES,
        apocrine   = apocrine_genes   %||% DEFAULT_APOCRINE_GENES,
        malo       = malo_genes       %||% DEFAULT_MALO_GENES
    )
    specification <- lamb_specification(gene_sets = gene_sets, calibration = cal)

    scores   <- lamb_score(data, gene_sets = gene_sets, min_genes_per_sig = min_genes_per_sig)
    coverage <- attr(scores, "coverage")
    low_cov  <- coverage$Signature[coverage$Coverage < cal$min_gene_coverage]

    predictions <- rep(NA_character_, nrow(scores))
    margins     <- rep(NA_real_, nrow(scores))
    flags       <- rep("", nrow(scores))

    # ---- Gate evaluability ----
    # An unscored gate cannot be treated as "no": the sample becomes NC.
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

    # ---- Step 1: metaplastic vs mammary ----
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

    # ---- Step 2: Basal-like vs luminal ----
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

    # ---- Step 3: LHS vs Apocrine, then MA-Lo ----
    in_luminal <- in_mammary & !is_basal

    if (any(in_luminal)) {
        idx  <- which(in_luminal)
        lum  <- scores[idx, "luminal"]
        apo  <- scores[idx, "apocrine"]
        malo <- scores[idx, "malo"]

        usable <- !is.na(lum) & !is.na(apo)
        if (any(usable)) {
            lum_wins <- usable & (lum > apo)
            predictions[idx[lum_wins]] <- "LHS"
            margins[idx[lum_wins]] <- lum[lum_wins] - apo[lum_wins]

            apo_wins <- usable & !lum_wins
            predictions[idx[apo_wins]] <- "Apocrine"
            margins[idx[apo_wins]] <- apo[apo_wins] - lum[apo_wins]

            max_lum_lineage <- pmax(lum, apo)
            is_malo <- usable & !is.na(malo) & (malo > max_lum_lineage)
            if (any(is_malo)) {
                predictions[idx[is_malo]] <- "MA-Lo"
                margins[idx[is_malo]] <- malo[is_malo] - max_lum_lineage[is_malo]
            }
        }
        if (any(!usable)) predictions[idx[!usable]] <- "NC"
    }

    predictions[unevaluable] <- "NC"
    margins[unevaluable] <- NA_real_

    predictions[is.na(predictions)] <- "NC"

    # ---- Quality flags ----
    add_flag <- function(f, which, tag) {
        f[which] <- ifelse(nzchar(f[which]), paste0(f[which], ";", tag), tag)
        f
    }

    borderline <- !is.na(margins) &
        predictions != "NC" &
        margins < cal$min_margin
    flags <- add_flag(flags, borderline, "close_call")

    # Per sample: missing values can drop coverage below the platform's.
    decision_relevant <- union(required, "malo")
    sample_cov <- attr(scores, "sample_coverage")
    per_low <- if (is.null(sample_cov)) rep(length(low_cov) > 0 && any(decision_relevant %in% low_cov), nrow(scores)) else
        apply(sample_cov[, intersect(decision_relevant, colnames(sample_cov)), drop = FALSE] < cal$min_gene_coverage, 1, any)
    flags <- add_flag(flags, per_low & predictions != "NC", "low_coverage")

    flags <- add_flag(flags, unevaluable, "incomplete_gates")

    # ---- Read-outs beside the class (never change it) ----
    amp <- .lamb_s2_her2_amplicon(data, cal)
    if (!is.null(rownames(scores))) amp <- amp[match(rownames(scores), amp$Sample), , drop = FALSE]
    gp <- .lamb_gp_percentiles(data, c(.LAMB_GP_MARKERS, unlist(gene_sets), cal$her2_amplicon$genes), cal$universe)

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
        Score_MALO       = round(scores[, "malo"], 4),
        HER2_Amplicon    = amp$HER2_Amplicon,
        HER2_Amplicon_N_High     = amp$N_High,
        HER2_Amplicon_N_Measured = amp$N_Measured
    )

    if (run_pam50 && !is.null(pam50_data)) {
        message("Running PAM50 classification...")
        pam50_res <- run_pam50_internal(data, pam50_data, gene_length_info)
        results <- dplyr::left_join(results, pam50_res, by = "Sample")
    }

    # set after any join: dplyr drops custom attributes
    attr(results, "signatures") <- gene_sets
    attr(results, "coverage") <- coverage
    attr(results, "sample_coverage") <- attr(scores, "sample_coverage")
    attr(results, "missing_values") <- attr(scores, "missing_values")
    attr(results, "universe_coverage") <- attr(scores, "universe_coverage")
    attr(results, "calibration_version") <- cal$version
    attr(results, "specification") <- specification
    attr(results, "gene_percentiles") <- gp$pct[match(results$Sample, rownames(gp$pct)), , drop = FALSE]
    attr(results, "n_genes_ranked") <- gp$n[match(results$Sample, names(gp$n))]

    results
}


# ---- Specification 2.0.0 read-outs ----

# HER2 amplicon modifier, one row per sample. Validated as written: do not
# change the rule without re-validating.
.lamb_s2_her2_amplicon <- function(data, cal = lamb_calibration()) {
    spec <- cal$her2_amplicon
    if (is.null(spec)) stop("The calibration has no `her2_amplicon` entry (specification 2.0.0 or later needed).",
                            call. = FALSE)
    if (!is.matrix(data)) data <- as.matrix(data)
    data <- data[!duplicated(rownames(data)), , drop = FALSE]   # first occurrence, as in lamb_score()
    mat <- data[intersect(cal$universe, rownames(data)), , drop = FALSE]
    module <- as.character(spec$genes)
    one <- function(v) {
        v <- v[!is.na(v)]
        measured <- intersect(module, names(v))
        if (length(measured) < spec$min_measured || length(v) < 2L) {
            return(c(N_High = NA_integer_, N_Measured = length(measured)))
        }
        pct <- 100 * (1 - (rank(v, ties.method = "average") - 1) / (length(v) - 1))
        c(N_High = sum(pct[measured] <= spec$q), N_Measured = length(measured))
    }
    counts <- vapply(seq_len(ncol(mat)), function(j) {
        v <- mat[, j]
        names(v) <- rownames(mat)
        one(v)
    }, c(N_High = 0L, N_Measured = 0L))
    n_high <- as.integer(counts["N_High", ])
    n_meas <- as.integer(counts["N_Measured", ])
    status <- ifelse(n_meas < spec$min_measured | is.na(n_high), "not assessable",
                     ifelse(n_high >= spec$min_high, "high", "not high"))
    data.frame(Sample = colnames(mat) %||% as.character(seq_len(ncol(mat))),
               HER2_Amplicon = status, N_High = n_high, N_Measured = n_meas,
               stringsAsFactors = FALSE)
}

# PAM50 intrinsic subtyping with the bundled centroids.
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

    # Not molecular.subtyping(): it finds `pam50.robust` only via the search path.
    output <- genefu::intrinsic.cluster.predict(
        sbt.model = pam50_centroids,
        data = dat,
        annot = annot_matrix,
        do.mapping = FALSE
    )

    tibble::tibble(
        Sample = rownames(dat),
        PAM50 = as.character(output$subtype)
    )
}


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

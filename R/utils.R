# Helper functions, class colours and orderings, and the default gene
# signatures used by LAMBclassifier.

`%||%` <- function(x, y) if (is.null(x)) y else x


# =============================================================================
# Standard class orderings and colour palettes
# =============================================================================

#' @title LAMBclassifier Class Order
#' @description Display order for the six LAMB classes.
#' @export
LAMB_CLASS_ORDER <- c(
    "LHS", "MA-Lo (HER2)", "Apocrine",
    "Basal-like", "Chondroid", "Squamous", "NC"
)

#' @rdname LAMB_CLASS_ORDER
#' @export
LAB_CLASS_ORDER <- LAMB_CLASS_ORDER

#' @title LAMBclassifier Class Colours
#' @description Colour palette for the LAMB classes.
#' @export
LAMB_CLASS_COLORS <- c(
    "LHS"           = "#000080", # Navy Blue (Luminal Hormone Sensing)
    "MA-Lo (HER2)"  = "#6BAED6", # Light Blue (Differentiated Luminal / HER2)
    "Apocrine"      = "purple",  # Purple (Molecular Apocrine)
    "Basal-like"    = "#CA2529", # Red (Basal-like / Luminal Progenitor)
    "Chondroid"     = "forestgreen", # Green (Chondroid Metaplastic)
    "Squamous"      = "#E69F00", # Orange (Squamous Metaplastic)
    "NC"            = "#888888"  # Grey (Not Classified / Unclassified)
)

#' @rdname LAMB_CLASS_COLORS
#' @export
LAB_CLASS_COLORS <- LAMB_CLASS_COLORS

#' @title PAM50 Class Order
#' @description Standard display order for PAM50 intrinsic subtypes.
#' @export
PAM50_ORDER <- c("LumA", "LumB", "Her2", "Basal", "Normal")

#' @title PAM50 Class Colours
#' @description Colour palette for PAM50 intrinsic subtypes.
#' @export
PAM50_COLORS <- c(
    "LumA"   = "#08519c", # Dark Blue
    "LumB"   = "#6baed6", # Light Blue
    "Her2"   = "#e377c2", # Pink
    "Basal"  = "#d62728", # Red
    "Normal" = "#2ca02c" # Green
)


#' @title Calculate Jaccard Index
#'
#' @description Computes Jaccard similarity between two sets.
#'
#' @param set1 Character vector
#' @param set2 Character vector
#'
#' @return Numeric Jaccard index (0-1)
#'
#' @export
calculate_jaccard <- function(set1, set2) {
    if (length(set1) == 0 && length(set2) == 0) {
        return(1)
    }
    length(intersect(set1, set2)) / length(union(set1, set2))
}


#' @title Summary of LAMBclassifier Results
#'
#' @description Prints a formatted summary of classification results.
#'
#' @param results Data frame returned by LAMBclassifier()
#'
#' @return Invisible NULL, prints summary to console
#'
#' @examples
#' \dontrun{
#' results <- LAMBclassifier(data)
#' summary_lamb_results(results)
#' }
#'
#' @export
summary_lamb_results <- function(results) {
    cat("\n=== LAMBclassifier Results Summary ===\n\n")

    # Class distribution
    cat("Classification Distribution:\n")
    class_table <- dplyr::count(results, Pred_Class, name = "n")
    class_table <- dplyr::mutate(class_table, Pred_Class = factor(Pred_Class, levels = LAMB_CLASS_ORDER))
    class_table <- dplyr::arrange(class_table, Pred_Class)
    for (i in seq_len(nrow(class_table))) {
        cls <- as.character(class_table$Pred_Class[i])
        pct <- round(100 * class_table$n[i] / sum(class_table$n), 1)
        cat(sprintf("  %-25s: %4d (%5.1f%%)\n", cls, class_table$n[i], pct))
    }

    # PAM50 if available
    if ("PAM50" %in% colnames(results)) {
        cat("\nPAM50 Distribution:\n")
        pam50_table <- dplyr::count(results, PAM50, name = "n")
        pam50_table <- dplyr::mutate(pam50_table, PAM50 = factor(PAM50, levels = PAM50_ORDER))
        pam50_table <- dplyr::arrange(pam50_table, PAM50)
        for (i in seq_len(nrow(pam50_table))) {
            p <- as.character(pam50_table$PAM50[i])
            cat(sprintf("  %s: %d\n", p, pam50_table$n[i]))
        }
    }

    cat("\n")
    invisible(NULL)
}


#' @title Export Results to CSV
#'
#' @description Exports LAMBclassifier results to a CSV file. When the results
#'   contain the specification attached by [LAMBclassifier()], companion
#'   provenance and ordered-signature TSV files are written beside the CSV.
#'
#' @param results Data frame returned by LAMBclassifier()
#' @param file_path Output file path
#' @param include_scores Logical, whether to include scoring columns
#'
#' @return Invisible file path
#'
#' @export
export_results <- function(results, file_path, include_scores = TRUE) {
    specification <- attr(results, "specification", exact = TRUE)

    if (!include_scores) {
        keep_cols <- c(
            "Sample", "Pred_Class", "PAM50", "Classifier_Version",
            "Signature_SHA256", "Universe_SHA256"
        )
        keep_cols <- keep_cols[keep_cols %in% colnames(results)]
        results <- results[, keep_cols, drop = FALSE]
    }

    write.csv(results, file_path, row.names = FALSE)
    message("Results exported to: ", file_path)

    if (!is.null(specification)) {
        stem <- sub("\\.[^.]+$", "", file_path)
        provenance_path <- paste0(stem, "_provenance.tsv")
        signatures_path <- paste0(stem, "_signatures.tsv")

        provenance <- data.frame(
            Key = c(
                "classifier_version", "signature_sha256", "universe_sha256",
                "universe_size", "min_margin", "min_gene_coverage", "exported_utc"
            ),
            Value = c(
                specification$version, specification$signature_sha256,
                specification$universe_sha256, specification$universe_size,
                specification$min_margin, specification$min_gene_coverage,
                format(Sys.time(), tz = "UTC", usetz = TRUE)
            ),
            stringsAsFactors = FALSE
        )
        utils::write.table(
            provenance, provenance_path, sep = "\t", quote = FALSE,
            row.names = FALSE
        )

        signatures <- do.call(rbind, lapply(names(specification$signatures), function(nm) {
            genes <- as.character(specification$signatures[[nm]])
            data.frame(
                Signature = nm, Order = seq_along(genes), Gene = genes,
                stringsAsFactors = FALSE
            )
        }))
        utils::write.table(
            signatures, signatures_path, sep = "\t", quote = FALSE,
            row.names = FALSE
        )
        message("Specification exported to: ", provenance_path, " and ", signatures_path)
    }

    invisible(file_path)
}

# ---- Default gene signatures ------------------------------------------------

#' @title Default Gene Signature: Squamous Epithelial
#' @description Squamous metaplastic carcinoma signature.
#' @export
DEFAULT_SQUAMOUS_EPITHELIAL_GENES <- c(
    "BNC1", "ZBED2", "TP63", "POU3F1", "FOXQ1", "KRT6A",
    "KRT5", "KRT14", "DSG3"
)

#' @title Default Gene Signature: Chondroid
#' @description Chondroid metaplastic carcinoma signature.
#' @export
DEFAULT_CHONDROID_GENES <- c(
    "SP7", "ACAN", "COL2A1", "PRELP",
    "TIMP4", "IBSP", "DLX5", "SOX11", "BMP2"
)

#' @title Default Gene Signature: Metaplastic (Union)
#' @description Union of DEFAULT_CHONDROID_GENES and DEFAULT_SQUAMOUS_EPITHELIAL_GENES.
#' @export
DEFAULT_METAPLASTIC_GENES <- unique(c(
    DEFAULT_CHONDROID_GENES,
    DEFAULT_SQUAMOUS_EPITHELIAL_GENES
))

#' @title Default Gene Signature: Basal-like
#' @description Basal-like breast cancer signature.
#' @export
DEFAULT_BASAL_LIKE_GENES <- c(
    "FOXC1", "SOX10", "ELF5", "NFIB",
    "SOX11", "SFRP1", "BCL11A", "GABRP",
    "ID4", "VGLL1", "KLF5", "ROPN1", "ROPN1B", "EN1"
)

#' @title Default Gene Signature: Luminal Hormone Sensing (LHS)
#' @description Luminal hormone sensing (LHS) signature.
#' @export
DEFAULT_LHS_GENES <- c(
    "ESR1", "GATA3", "MYB", "PGR", "RERG",
    "ANKRD30A", "STC2", "TFF1", "CA12", "INPP4B",
    "AGR3", "ANXA9", "BCL2", "FOXA1", "SPDEF"
)

#' @title Default Gene Signature: Apocrine
#' @description Molecular apocrine carcinoma signature.
#' @export
DEFAULT_APOCRINE_GENES <- c(
    "KYNU", "SRD5A1", "CLDN8", "TSPAN8", "HPD",
    "ENPP3", "CUX2", "FKBP5", "UGT2B28", "CLCA2",
    "HPGD", "FOXA1", "AR", "SERHL2"
)

#' @title Default Gene Signature: MA-Lo (HER2)
#' @description ERBB2-amplicon expression and co-target signature (including
#'   FGFR4 and PRODH). This is not a direct DNA-amplification assay.
#' @export
DEFAULT_MALO_HER2_GENES <- c(
    "ERBB2", "GRB7", "MIEN1", "PGAP3",
    "STARD3", "TCAP", "PNMT", "FGFR4", "PRODH"
)

#' @title All Default LAMBclassifier Input Genes
#' @description Consolidated list of all genes across default signatures.
#' @export
DEFAULT_ALL_INPUT_GENES <- unique(c(
    DEFAULT_CHONDROID_GENES,
    DEFAULT_SQUAMOUS_EPITHELIAL_GENES,
    DEFAULT_BASAL_LIKE_GENES,
    DEFAULT_LHS_GENES,
    DEFAULT_APOCRINE_GENES,
    DEFAULT_MALO_HER2_GENES
))

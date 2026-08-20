# =============================================================================
# LAMBclassifier: Data Preparation Utilities
# =============================================================================
#
# LAMBclassifier expects a Gene × Sample expression matrix in log2-transformed
# space. This file provides helper functions to convert raw RNA-seq count
# matrices or linear-scale TPM files into the expected format.
#
# WHICH FORMAT DOES MY DATA NEED TO BE IN?
# -----------------------------------------
# | Input format              | Action needed                              |
# |---------------------------|--------------------------------------------|
# | Raw counts (STAR, RSEM)   | raw_counts = TRUE  (converts to TPM first) |
# | Linear TPM (already norm) | log2_transform = TRUE                      |
# | log2 TPM / log2 counts+1  | Neither flag needed (already correct)      |
#
# GENE IDENTIFIERS
# ----------------
# Row names must be HGNC gene symbols (e.g. "ESR1", "TP63"). Ensembl IDs are
# not directly supported. Convert with biomaRt or org.Hs.eg.db before calling
# these functions.
#
# =============================================================================

#' Prepare Expression Data for LAMBclassifier
#'
#' @description
#' Validates and optionally transforms a gene expression matrix into the
#' log2-scale format expected by \code{LAMBclassifier()}. Handles two common
#' preprocessing needs: (1) converting raw counts to TPM using gene lengths,
#' and (2) log2-transforming linear-scale TPM values.
#'
#' @param data Gene × Sample expression matrix (data frame or matrix).
#'   Rows must be HGNC gene symbols; columns must be sample IDs.
#' @param log2_transform Logical. Apply \code{log2(x + 1)} transformation.
#'   Set \code{TRUE} when input is in linear TPM scale (default \code{FALSE}).
#' @param raw_counts Logical. If \code{TRUE}, convert raw counts to TPM
#'   using gene lengths before log2 transformation. Requires
#'   \code{gene_length_info} (default \code{FALSE}).
#' @param gene_length_info Data frame with columns \code{SYMBOL} and
#'   \code{length} (in bp). Required when \code{raw_counts = TRUE}. The
#'   bundled \code{gene.length} object in \code{sysdata.rda} can be used.
#' @param id_type Character. Gene identifier type in row names. Only
#'   \code{"SYMBOL"} (HGNC) is currently supported (default \code{"SYMBOL"}).
#'
#' @examples
#' \dontrun{
#' # Load raw TPM data and log2 transform
#' expr <- read.table("tpm_matrix.txt", header = TRUE, row.names = 1)
#' data <- prepare_lamb_data(expr, log2_transform = TRUE)
#'
#' # Data already log2 transformed
#' data <- prepare_lamb_data(log2_expr, log2_transform = FALSE)
#' }
#'
#' @export
prepare_lamb_data <- function(data,
                           log2_transform = FALSE,
                           raw_counts = FALSE,
                           gene_length_info = NULL,
                           id_type = "SYMBOL") {
    # Validate input
    if (!is.data.frame(data) && !is.matrix(data)) {
        stop("Input data must be a data.frame or matrix")
    }

    # Convert to matrix if needed
    if (is.data.frame(data)) {
        data <- as.matrix(data)
    }

    # Check for gene symbols
    if (id_type != "SYMBOL") {
        if (is.null(gene_length_info)) {
            stop("Gene length info (sysdata.rda) required for ID mapping")
        }
        message("Note: ID mapping currently supports SYMBOL. Ensure input has Gene Symbols.")
    }

    # TPM conversion from raw counts
    if (raw_counts) {
        if (is.null(gene_length_info)) {
            stop("gene_length_info required for raw count to TPM conversion")
        }
        message("Transforming raw counts to TPM...")
        data <- counts_to_tpm(data, gene_length_info)
    }

    # Log2 transformation
    if (log2_transform) {
        message("Applying log2(x+1) transformation...")
        data <- log2(data + 1)
    }

    # Check for negative values
    if (any(data < 0, na.rm = TRUE)) {
        warning(
            "Negative values detected in expression data. This may indicate ",
            "data is already log2 transformed or contains errors."
        )
    }

    return(data)
}


#' @title Convert Raw Counts to TPM
#'
#' @description Converts raw count data to Transcripts Per Million (TPM).
#'
#' @param counts Gene x Sample count matrix
#' @param gene_length_info Data frame with SYMBOL and length columns
#'
#' @return TPM normalized expression matrix
#'
#' @keywords internal
counts_to_tpm <- function(counts, gene_length_info) {
    # Match genes
    common_genes <- intersect(rownames(counts), gene_length_info$SYMBOL)

    if (length(common_genes) < 100) {
        warning("Only ", length(common_genes), " genes matched. Check gene symbol format.")
    }

    counts <- counts[common_genes, , drop = FALSE]
    gene_lengths <- gene_length_info$length[match(common_genes, gene_length_info$SYMBOL)]

    # Calculate RPK (reads per kilobase)
    rpk <- counts / (gene_lengths / 1000)

    # Calculate TPM
    tpm <- t(t(rpk) / colSums(rpk) * 1e6)

    return(tpm)
}


#' @title Load Expression Data from File
#'
#' @description Convenience function to load expression data from common formats.
#'
#' @param file_path Path to expression data file
#' @param sep Field separator (default "\\t" for tab-delimited)
#' @param header Logical, whether file has header row (default TRUE)
#' @param row.names Column index for row names (default 1)
#' @param log2_transform Logical, whether to apply log2 transformation
#'
#' @return Expression matrix ready for LAMBclassifier
#'
#' @examples
#' \dontrun{
#' data <- load_expression_data("expression.txt", log2_transform = TRUE)
#' }
#'
#' @export
load_expression_data <- function(file_path,
                                 sep = "\t",
                                 header = TRUE,
                                 row.names = 1,
                                 log2_transform = FALSE) {
    if (!file.exists(file_path)) {
        stop("File not found: ", file_path)
    }

    data <- read.table(file_path,
        sep = sep,
        header = header,
        row.names = row.names,
        check.names = FALSE
    )

    return(prepare_lamb_data(data, log2_transform = log2_transform))
}

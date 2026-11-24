#' Prepare Expression Data for LAMBclassifier
#'
#' @description
#' Validates and optionally transforms a gene expression matrix into the
#' log2-scale format expected by \code{LAMBclassifier()}: raw counts are
#' converted to TPM, and linear TPM is log2-transformed. Row names must be HGNC
#' symbols.
#'
#' @param data Gene × Sample expression matrix (data frame or matrix).
#'   Rows must be HGNC gene symbols; columns must be sample IDs.
#' @param log2_transform Logical. Apply \code{log2(x + 1)} transformation.
#'   Set \code{TRUE} when input is in linear TPM scale (default \code{FALSE}).
#' @param raw_counts Logical. If \code{TRUE}, convert raw counts to TPM
#'   using gene lengths before log2 transformation. Requires
#'   \code{gene_length_info} (default \code{FALSE}).
#' @param gene_length_info Data frame with columns \code{SYMBOL} and
#'   \code{length} or \code{effective_length} (in bp). Required when
#'   \code{raw_counts = TRUE}. \code{\link{lamb_run}} uses the bundled gene
#'   lengths automatically.
#' @param id_type Character. Gene identifier type in row names. Only
#'   \code{"SYMBOL"} (HGNC) is supported; for Ensembl IDs use \code{\link{lamb_run}}.
#'
#' @return A numeric gene x sample matrix.
#'
#' @examples
#' \dontrun{
#' # Linear TPM
#' expr <- read.table("tpm_matrix.txt", header = TRUE, row.names = 1)
#' data <- prepare_lamb_data(expr, log2_transform = TRUE)
#'
#' # Already log2
#' data <- prepare_lamb_data(log2_expr, log2_transform = FALSE)
#' }
#'
#' @export
prepare_lamb_data <- function(data,
                           log2_transform = FALSE,
                           raw_counts = FALSE,
                           gene_length_info = NULL,
                           id_type = "SYMBOL") {
    if (!is.data.frame(data) && !is.matrix(data)) {
        stop("Input data must be a data.frame or matrix")
    }

    if (is.data.frame(data)) {
        data <- as.matrix(data)
    }

    if (id_type != "SYMBOL") {
        stop("prepare_lamb_data() takes HGNC gene symbols only. lamb_run() maps Ensembl IDs to symbols.", call. = FALSE)
    }

    if (raw_counts) {
        if (is.null(gene_length_info)) {
            stop("gene_length_info required for raw count to TPM conversion")
        }
        message("Transforming raw counts to TPM...")
        data <- counts_to_tpm(data, gene_length_info)
    }

    if (log2_transform) {
        message("Applying log2(x+1) transformation...")
        data <- log2(data + 1)
    }

    if (any(data < 0, na.rm = TRUE)) {
        warning(
            "Negative values detected in expression data. This may indicate ",
            "data is already log2 transformed or contains errors."
        )
    }

    return(data)
}


# Raw counts to TPM; gene_length_info needs SYMBOL and a length column.
counts_to_tpm <- function(counts, gene_length_info) {
    common_genes <- intersect(rownames(counts), gene_length_info$SYMBOL)

    if (length(common_genes) < 100) {
        warning("Only ", length(common_genes), " genes matched. Check gene symbol format.")
    }

    counts <- counts[common_genes, , drop = FALSE]
    len_col <- if ("length" %in% names(gene_length_info)) "length" else "effective_length"  # bundled table: effective_length
    gene_lengths <- gene_length_info[[len_col]][match(common_genes, gene_length_info$SYMBOL)]
    ok <- !is.na(gene_lengths) & gene_lengths > 0
    counts <- counts[ok, , drop = FALSE]
    gene_lengths <- gene_lengths[ok]

    rpk <- counts / (gene_lengths / 1000)

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

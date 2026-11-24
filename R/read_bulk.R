# ---- Bulk expression input reader ----
# Formats are recognised by header names, never by column position; ambiguous
# files stop with a message rather than a guess. Never maps IDs or rescales.
#
# Returns list(expr = numeric matrix (genes x samples; rownames = gene IDs as
# found, colnames = sample names), units = "counts" | "tpm" | "fpkm" | "log2" |
# "unknown", platform = "table" | "salmon" | "kallisto" | "featurecounts" |
# "star" | "htseq" | "rsem" | "xlsx" | "rds" | "multi_file", source_platform =
# the per-file platform for a folder, otherwise the same as platform).

.LAMB_ENSEMBL_GENE_RE <- "^ENS[A-Z]*G[0-9]{6,}(\\.[0-9]+)?(_PAR_Y)?$"
.LAMB_TRANSCRIPT_RE   <- "^(ENS[A-Z]*T[0-9]{6,}(\\.[0-9]+)?(_PAR_Y)?|[NX][MR]_[0-9]+(\\.[0-9]+)?)$"
.LAMB_TEXT_EXT_RE     <- "\\.(sf|tsv|txt|csv|tab|results|counts)(\\.gz)?$"

.lamb_read_bulk <- function(input, say = message) {
    if (!file.exists(input)) stop("File or folder not found: ", input, call. = FALSE)
    res <- if (dir.exists(input)) .lamb_read_bulk_folder(input, say) else .lamb_read_bulk_file(input, say)
    .lamb_bulk_final_checks(res, say)
}


# ---- one file ----

.lamb_read_bulk_file <- function(path, say) {
    b <- tolower(basename(path))
    if (grepl("\\.rds$", b)) return(.lamb_bulk_from_rds(path, say))
    if (grepl("\\.xlsx?$", b)) return(.lamb_bulk_from_xlsx(path, say))
    if (grepl("\\.h5$", b)) {
        stop("Binary .h5 quantification files (e.g. kallisto abundance.h5) cannot be read. ",
             "Give the text version instead (kallisto abundance.tsv).", call. = FALSE)
    }
    if (grepl("\\.(bam|sam|cram|fastq|fq|fastq\\.gz|fq\\.gz)$", b)) {
        stop("This is a raw sequencing/alignment file (", basename(path), "). LAMB needs a table of ",
             "gene expression values (e.g. Salmon quant.sf, featureCounts output, or a genes x samples table).",
             call. = FALSE)
    }
    .lamb_bulk_from_text(path, say)
}

# Name used for a single-sample file: the file name stem; for generic names
# (quant.sf, abundance.tsv, ReadsPerGene.out.tab) the folder holding the file.
.lamb_sample_stem <- function(path) {
    b <- basename(path)
    b <- sub("\\.gz$", "", b, ignore.case = TRUE)
    b <- sub("\\.(genes|isoforms)\\.results$", "", b, ignore.case = TRUE)
    b <- sub("[._]?ReadsPerGene\\.out\\.tab$", "", b)
    b <- sub("\\.(sf|tsv|txt|csv|tab|results|counts|xlsx|xls|rds)$", "", b, ignore.case = TRUE)
    if (tolower(b) %in% c("", "quant", "quant.genes", "abundance")) {
        b <- basename(dirname(normalizePath(path, mustWork = FALSE)))
    }
    b
}

.lamb_known_symbols <- local({
    cache <- NULL
    function() {
        if (is.null(cache)) {
            gl <- get_gene_length()
            s <- if (is.null(gl)) character() else unique(gl$SYMBOL)
            cache <<- s[!is.na(s) & nzchar(s)]
        }
        cache
    }
})

# Fraction of values that are known HGNC symbols, and that are Ensembl gene IDs.
.lamb_gene_score <- function(v) {
    v <- as.character(v)
    v <- v[!is.na(v) & nzchar(v)]
    if (!length(v)) return(c(symbol = 0, ensembl = 0))
    if (length(v) > 5000) v <- v[unique(round(seq(1, length(v), length.out = 5000)))]
    c(symbol = mean(v %in% .lamb_known_symbols()), ensembl = mean(grepl(.LAMB_ENSEMBL_GENE_RE, v)))
}

.lamb_looks_gene <- function(v, min = 0.25) max(.lamb_gene_score(v)) >= min

.lamb_examples <- function(v, n = 3) {
    v <- unique(as.character(v[!is.na(v) & nzchar(as.character(v))]))
    paste(utils::head(v, n), collapse = ", ")
}


# ---- text files ----

.lamb_split_line <- function(line, sep) {
    scan(text = line, what = "", sep = sep, quote = "\"", quiet = TRUE,
         na.strings = character(), strip.white = TRUE, blank.lines.skip = FALSE)
}

.lamb_sniff_sep <- function(lines) {
    lines <- utils::head(lines, 50)
    count <- function(ch) vapply(lines, function(l) sum(strsplit(l, "", fixed = TRUE)[[1]] == ch), 0)
    k_tab <- count("\t")
    if (all(k_tab >= 1)) return("\t")
    k_semi <- count(";")
    if (all(k_semi >= 1)) return(";")
    if (all(count(",") >= 1)) return(",")
    if (all(grepl("[ ]", lines))) return("")
    stop("Could not tell how the columns in this file are separated (expected tabs, commas or ",
         "semicolons). Save the table as tab-separated or comma-separated text.", call. = FALSE)
}

.lamb_sniff_dec <- function(lines, sep) {
    if (sep == ",") return(".")
    fields <- unlist(lapply(utils::head(lines, 50), function(l) {
        if (sep == "") strsplit(trimws(l), "[ \t]+")[[1]] else strsplit(l, sep, fixed = TRUE)[[1]]
    }))
    fields <- gsub("\"", "", trimws(fields))
    comma <- sum(grepl("^-?[0-9]+,[0-9]+([eE][-+]?[0-9]+)?$", fields))
    dot <- sum(grepl("^-?[0-9]*\\.[0-9]+([eE][-+]?[0-9]+)?$", fields))
    if (comma > dot) "," else "."
}

.lamb_num <- function(v, dec = ".") {
    v <- trimws(as.character(v))
    # Thousands separators, as spreadsheets export them: "1,234.5" / "1.234,5" / "1 234".
    v <- gsub("\xc2\xa0", " ", v, fixed = TRUE, useBytes = TRUE)   # no-break space (UTF-8 bytes), any locale
    v <- gsub("(?<=[0-9])[ '](?=[0-9]{3}\\b)", "", v, perl = TRUE)
    if (dec == ",") {
        th <- grepl("^[+-]?[0-9]{1,3}(\\.[0-9]{3})+(,[0-9]*)?$", v)
        v[th] <- gsub(".", "", v[th], fixed = TRUE)
        v <- gsub(",", ".", v, fixed = TRUE)
    } else {
        th <- grepl("^[+-]?[0-9]{1,3}(,[0-9]{3})+(\\.[0-9]*)?$", v)
        v[th] <- gsub(",", "", v[th], fixed = TRUE)
    }
    suppressWarnings(as.numeric(v))
}

.lamb_bulk_from_text <- function(path, say) {
    con <- file(path, "r")   # file() reads gzip transparently
    head_lines <- tryCatch(readLines(con, n = 200, warn = FALSE), finally = close(con))
    if (length(head_lines)) head_lines[1] <- sub("^\xef\xbb\xbf", "", head_lines[1], useBytes = TRUE)
    if (any(grepl("^!(Series_|series_matrix_table_begin|Sample_)", head_lines))) {
        stop(basename(path), " is a GEO series matrix as downloaded: it holds microarray probe IDs and ",
             "description lines, not one row per gene. Map the probes to gene symbols with the platform ",
             "annotation (GPL) and save a table with one row per gene, then run LAMB on that table.", call. = FALSE)
    }
    skip <- 0L
    while (skip < length(head_lines) &&
           (grepl("^#", head_lines[skip + 1]) || !nzchar(trimws(head_lines[skip + 1])))) skip <- skip + 1L
    body <- if (skip > 0) head_lines[-seq_len(skip)] else head_lines
    body <- body[nzchar(trimws(body))]
    if (!length(body)) stop("The file ", basename(path), " is empty.", call. = FALSE)
    if (length(body) == 1) stop("The file ", basename(path), " has only one line: no expression values found.", call. = FALSE)

    sep <- .lamb_sniff_sep(body)
    dec <- .lamb_sniff_dec(body[-1], sep)
    first <- .lamb_split_line(body[1], sep)
    stem <- .lamb_sample_stem(path)

    is_star <- length(first) == 4 && first[1] %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")
    headerless <- is_star || (length(first) >= 2 && !is.na(.lamb_num(first[2], dec)) &&
                              all(!is.na(.lamb_num(first[-1], dec))) &&
                              (.lamb_looks_gene(first[1], 1) || grepl("^__", first[1])))
    if (grepl("^status$", first[1], ignore.case = TRUE) && any(grepl("^Assigned", body[-1]))) {
        stop("This is a featureCounts summary file (", basename(path), "). Give the main featureCounts ",
             "output instead (the file without .summary).", call. = FALSE)
    }

    tab <- tryCatch(
        utils::read.table(path, sep = sep, dec = dec, header = FALSE, skip = skip + !headerless,
                          quote = "\"", comment.char = "", na.strings = c("NA", ""),
                          strip.white = sep != "\t", stringsAsFactors = FALSE, check.names = FALSE),
        error = function(e) {
            stop("Could not read ", basename(path), " as a table: rows have different numbers of columns ",
                 "or the separator is inconsistent (", conditionMessage(e), "). Check that gene names and ",
                 "sample names do not contain the column separator.", call. = FALSE)
        })
    if (sep == ";" || dec == ",") {
        say("Read ", basename(path), " with '", if (sep == "\t") "tab" else if (sep == "") "space" else sep,
            "' as separator and ", if (dec == ",") "decimal comma" else "decimal point", ".")
    }

    if (is_star) return(.lamb_bulk_star(tab, stem, path, say))
    if (headerless) {
        if (ncol(tab) != 2) {
            stop("The file ", basename(path), " has no header row, so the sample names are unknown. ",
                 "Add a first row with a name for each column.", call. = FALSE)
        }
        ids <- as.character(tab[[1]])
        drop <- grepl("^__", ids)
        platform <- "table"; units <- "unknown"
        if (any(drop)) {
            say("HTSeq-count file: dropped ", sum(drop), " summary rows (", paste(ids[drop], collapse = ", "), ").")
            platform <- "htseq"; units <- "counts"
        }
        x <- matrix(as.double(tab[[2]][!drop]), ncol = 1, dimnames = list(ids[!drop], stem))
        say("No header row: one sample, named after the file (", stem, ").")
        return(list(expr = x, units = units, platform = platform))
    }

    hdr <- .lamb_split_line(body[1], sep)
    if (length(hdr) == ncol(tab) - 1) hdr <- c("", hdr)
    if (length(hdr) != ncol(tab)) {
        stop("The header row of ", basename(path), " has ", length(hdr), " names but the rows have ",
             ncol(tab), " columns. Check the separator and that every column has a name.", call. = FALSE)
    }
    .lamb_bulk_from_df(as.list(tab), hdr, stem, say, dec = dec, platform = "table")
}

.lamb_bulk_star <- function(tab, stem, path, say) {
    ids <- as.character(tab[[1]])
    keep <- !grepl("^N_", ids)
    totals <- vapply(2:4, function(j) sum(as.numeric(tab[[j]][keep]), na.rm = TRUE), 0)
    labels <- c("unstranded", "read 1 on the gene strand (forward)", "read 2 on the gene strand (reverse)")
    j <- which.max(totals)
    say("STAR ReadsPerGene file: dropped the 4 summary rows (N_unmapped, N_multimapping, N_noFeature, ",
        "N_ambiguous); used column ", j + 1, " (", labels[j], "), the one with the most counts (",
        paste0(c("unstranded ", "forward ", "reverse "), format(totals, big.mark = ","), collapse = "; "), ").")
    x <- matrix(as.double(tab[[j + 1]][keep]), ncol = 1, dimnames = list(ids[keep], stem))
    list(expr = x, units = "counts", platform = "star")
}


# ---- Excel and RDS ----

.lamb_bulk_from_xlsx <- function(path, say) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
        stop("Reading Excel files needs the 'readxl' package. Install it with install.packages(\"readxl\"), ",
             "or save the sheet as a CSV file and use that instead.", call. = FALSE)
    }
    sheets <- readxl::excel_sheets(path)
    df <- readxl::read_excel(path, sheet = 1, col_names = TRUE, .name_repair = "minimal")
    say("Excel file: read the first sheet ('", sheets[1], "')",
        if (length(sheets) > 1) paste0("; ignored ", length(sheets) - 1, " other sheet(s)") else "", ".")
    hdr <- names(df)
    hdr[is.na(hdr)] <- ""
    cols <- lapply(as.list(df), function(v) if (inherits(v, c("POSIXt", "Date"))) as.character(v) else v)
    res <- .lamb_bulk_from_df(cols, hdr, .lamb_sample_stem(path), say, platform = "xlsx")
    if (res$platform == "table") res$platform <- "xlsx"
    res
}

.lamb_bulk_from_rds <- function(path, say) {
    obj <- readRDS(path)
    if (is.matrix(obj) && is.numeric(obj)) {
        if (is.null(rownames(obj)) && is.null(colnames(obj))) {
            stop("The matrix in ", basename(path), " has no gene names (row names). Save it with gene names ",
                 "as row names and sample names as column names.", call. = FALSE)
        }
        x <- .lamb_orient(obj, say)
        storage.mode(x) <- "double"
        if (is.null(colnames(x))) colnames(x) <- paste0("Sample_", seq_len(ncol(x)))
        return(list(expr = x, units = "unknown", platform = "rds"))
    }
    if (is.matrix(obj)) obj <- as.data.frame(obj, stringsAsFactors = FALSE)
    if (is.data.frame(obj)) {
        auto_rn <- .row_names_info(obj) < 0
        cols <- lapply(as.list(obj), function(v) if (is.factor(v)) as.character(v) else v)
        res <- .lamb_bulk_from_df(cols, names(obj), .lamb_sample_stem(path), say,
                                  row_ids = if (auto_rn) NULL else rownames(obj), platform = "rds")
        res$platform <- "rds"
        return(res)
    }
    stop("The .rds file holds a ", class(obj)[1], "; for bulk data it should hold a genes x samples ",
         "matrix or data frame.", call. = FALSE)
}

# Numeric matrix with dimnames: genes must be in rows; transpose when only the
# column names look like genes.
.lamb_orient <- function(x, say) {
    if (!is.null(rownames(x)) && .lamb_looks_gene(rownames(x))) return(x)
    if (!is.null(colnames(x)) && .lamb_looks_gene(colnames(x))) {
        say("Genes are in the columns and samples in the rows: the table was turned round (transposed).")
        return(t(x))
    }
    stop("Neither the row names nor the column names look like gene names (row names start ",
         .lamb_examples(rownames(x)), "). LAMB needs gene symbols (e.g. ESR1) or Ensembl gene IDs ",
         "(e.g. ENSG00000091831).", call. = FALSE)
}


# ---- a parsed table: quantifier formats, annotation columns, orientation ----

.lamb_has <- function(h, need) all(need %in% h)

.lamb_bulk_from_df <- function(cols, hdr, stem, say, dec = ".", row_ids = NULL, platform = "table") {
    hdr <- trimws(as.character(hdr))
    h <- tolower(hdr)
    col <- function(name) cols[[match(name, h)]]

    # ---- single-sample quantifier outputs, recognised by their header names ----
    if (.lamb_has(h, c("name", "effectivelength", "numreads"))) {
        if (!"tpm" %in% h) stop("This looks like a Salmon file but has no TPM column.", call. = FALSE)
        say("Salmon quantification file: used the TPM column (one sample: ", stem, ").")
        return(.lamb_tx_or_gene(col("name"), .lamb_num(col("tpm"), dec), stem, "Salmon", "salmon", say))
    }
    if (.lamb_has(h, c("target_id", "eff_length", "est_counts"))) {
        if (!"tpm" %in% h) stop("This looks like a kallisto file but has no tpm column.", call. = FALSE)
        say("kallisto quantification file: used the tpm column (one sample: ", stem, ").")
        return(.lamb_tx_or_gene(col("target_id"), .lamb_num(col("tpm"), dec), stem, "kallisto", "kallisto", say))
    }
    if (.lamb_has(h, c("gene_id", "expected_count", "effective_length"))) {
        if (!"tpm" %in% h) stop("This looks like an RSEM file but has no TPM column.", call. = FALSE)
        tpm <- .lamb_num(col("tpm"), dec)
        ids <- as.character(col("gene_id"))
        if (match("transcript_id", h, nomatch = 0) == 1) {
            say("RSEM isoforms file: summed transcript TPM into genes using the gene_id column.")
            ok <- !is.na(ids) & nzchar(ids)
            s <- rowsum(tpm[ok], ids[ok], reorder = FALSE)
            x <- matrix(s[, 1], ncol = 1, dimnames = list(rownames(s), stem))
        } else {
            say("RSEM gene results file: used the TPM column (one sample: ", stem, ").")
            x <- matrix(tpm, ncol = 1, dimnames = list(ids, stem))
        }
        return(list(expr = x, units = "tpm", platform = "rsem"))
    }
    fc_annot <- c("geneid", "chr", "start", "end", "strand", "length")
    if (.lamb_has(h, fc_annot)) {
        cnt <- which(!h %in% fc_annot)
        if (!length(cnt)) stop("This looks like a featureCounts file but has no count columns.", call. = FALSE)
        nm <- basename(hdr[cnt])
        nm <- sub("\\.(bam|sam|cram)$", "", nm, ignore.case = TRUE)
        nm <- sub("[._]?Aligned(\\.sortedByCoord)?\\.out$", "", nm)
        nm <- sub("[._]sorted$", "", nm)
        if (length(cnt) == 1 && !nzchar(nm)) nm <- stem
        x <- do.call(cbind, lapply(cols[cnt], function(v) if (is.numeric(v)) as.double(v) else .lamb_num(v, dec)))
        dimnames(x) <- list(as.character(col("geneid")), nm)
        say("featureCounts file: used ", length(cnt), " count column(s); ignored Chr, Start, End, Strand, Length.")
        return(list(expr = x, units = "counts", platform = "featurecounts"))
    }

    # ---- generic table ----
    n <- length(cols)
    if (!is.null(row_ids)) {
        cols <- c(list(as.character(row_ids)), cols)
        hdr <- c("(row names)", hdr); h <- tolower(hdr); n <- n + 1
    }
    blank_na <- vapply(seq_len(n), function(j) !nzchar(hdr[j]) && all(is.na(cols[[j]])), TRUE)
    if (any(blank_na)) {
        cols <- cols[!blank_na]; hdr <- hdr[!blank_na]; h <- h[!blank_na]; n <- length(cols)
    }
    is_num <- vapply(cols, function(v) is.numeric(v) || (is.logical(v) && all(is.na(v))), TRUE)
    junk <- 0
    for (j in which(!is_num)) {
        v <- cols[[j]]
        if (!is.character(v) && !is.logical(v)) next
        present <- !is.na(v) & nzchar(trimws(as.character(v)))
        if (!any(present)) next
        num <- .lamb_num(v, dec)
        if (mean(!is.na(num[present])) >= 0.9 && !.lamb_looks_gene(v)) {
            junk <- junk + sum(present & is.na(num))
            cols[[j]] <- num
            is_num[j] <- TRUE
        }
    }
    if (junk > 0) say(junk, " cell(s) held text where a number was expected; they are treated as missing.")

    annot_re <- "^(entrez([ ._]?gene)?([ ._]?id)?|gene[ ._]?length|length|eff(ective)?[ ._]?length|chr(om(osome)?)?|start|end|strand|gc([ ._]?content)?|width)$"
    annot <- is_num & grepl(annot_re, h)
    if (any(annot)) {
        say("Ignored annotation column(s) that are not samples: ", paste(hdr[annot], collapse = ", "), ".")
        is_num[annot] <- FALSE
    }
    if (!any(is_num)) stop("No column of numbers was found: this does not look like an expression table.", call. = FALSE)

    char_idx <- which(!is_num & !annot)
    scores <- vapply(char_idx, function(j) .lamb_gene_score(cols[[j]]), c(symbol = 0, ensembl = 0))
    gene <- NA_integer_
    if (length(char_idx)) {
        sym <- scores["symbol", ]; ens <- scores["ensembl", ]
        if (max(sym) >= 0.5 || (max(sym) >= 0.25 && max(sym) >= max(ens))) {
            gene <- char_idx[which.max(sym)]
        } else if (max(ens) >= 0.25) {
            gene <- char_idx[which.max(ens)]
        } else if (max(sym) >= 0.25) {
            gene <- char_idx[which.max(sym)]
        }
    }

    if (is.na(gene)) {
        if (.lamb_looks_gene(hdr[is_num])) return(.lamb_bulk_transposed(cols, hdr, is_num, char_idx, say, platform))
        what <- if (length(char_idx)) paste0("the first text column ('", hdr[char_idx[1]], "') starts ",
                                             .lamb_examples(cols[[char_idx[1]]]))
                else "there is no column of names"
        if (length(char_idx) && .lamb_looks_gene(toupper(as.character(cols[[char_idx[1]]])))) {
            stop("The gene names are not in upper case (", .lamb_examples(cols[[char_idx[1]]]), "). Names like ",
                 "'Gata3' are mouse genes, and LAMB is built for human tumours only. If these are human genes ",
                 "written in lower case, convert the names to upper case (GATA3) and run again.", call. = FALSE)
        }
        stop("Could not find a column of gene names: ", what, ". LAMB needs gene symbols (e.g. ESR1) or ",
             "Ensembl gene IDs (e.g. ENSG00000091831). If these are microarray probe IDs, convert them to ",
             "gene symbols first.", call. = FALSE)
    }

    other <- setdiff(char_idx, gene)
    if (length(other)) {
        lab <- ifelse(nzchar(hdr[other]), hdr[other], paste0("column ", other))
        say("Gene names taken from column '", if (nzchar(hdr[gene])) hdr[gene] else paste("column", gene),
            "'; ignored other text column(s): ", paste(lab, collapse = ", "), ".")
    }
    ids <- as.character(cols[[gene]])
    keep_rows <- !is.na(ids) & nzchar(trimws(ids))
    if (any(!keep_rows)) say(sum(!keep_rows), " row(s) without a gene name were dropped.")
    .lamb_warn_dates(ids, say)

    num_idx <- which(is_num)
    x <- do.call(cbind, lapply(cols[num_idx], as.double))
    if (any(!keep_rows)) x <- x[keep_rows, , drop = FALSE]
    dimnames(x) <- list(ids[keep_rows], .lamb_sample_names(hdr[num_idx], say))
    list(expr = x, units = "unknown", platform = platform)
}

.lamb_bulk_transposed <- function(cols, hdr, is_num, char_idx, say, platform) {
    num_idx <- which(is_num)
    ids <- NULL
    for (j in char_idx) {
        v <- as.character(cols[[j]])
        if (!anyNA(v) && !anyDuplicated(v)) { ids <- v; break }
    }
    nr <- length(cols[[num_idx[1]]])
    if (is.null(ids)) ids <- paste0("Sample_", seq_len(nr))
    other <- setdiff(char_idx, char_idx[vapply(char_idx, function(j) identical(as.character(cols[[j]]), ids), TRUE)])
    say("Genes are in the columns and samples in the rows: the table was turned round (transposed).",
        if (length(other)) paste0(" Ignored text column(s): ", paste(hdr[other], collapse = ", "), ".") else "")
    x <- do.call(cbind, lapply(cols[num_idx], as.double))
    dimnames(x) <- list(.lamb_sample_names(ids, say), hdr[num_idx])
    list(expr = t(x), units = "unknown", platform = platform)
}

.lamb_sample_names <- function(nm, say) {
    nm <- as.character(nm)
    blank <- is.na(nm) | !nzchar(nm)
    if (any(blank)) {
        nm[blank] <- paste0("Sample_", which(blank))
        say(sum(blank), " sample column(s) had no name and were called ", paste(nm[blank], collapse = ", "), ".")
    }
    if (anyDuplicated(nm)) {
        say("Some sample names were repeated (", .lamb_examples(nm[duplicated(nm)]),
            "); a suffix was added to keep them apart.")
        nm <- make.unique(nm, sep = "_")
    }
    nm
}

# Gene names that Excel turned into dates (SEPT1 -> 1-Sep, MARCH1 -> 1-Mar).
.lamb_warn_dates <- function(ids, say) {
    d <- grepl("^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$|^[0-9]{4}-[0-9]{2}-[0-9]{2}", ids)
    if (any(d)) say(sum(d), " gene name(s) look like dates (e.g. ", .lamb_examples(ids[d]),
                    "): Excel probably converted them. Those genes will not be matched.")
}

# Transcript-level Salmon/kallisto output: sum to genes only when the IDs
# carry the gene (GENCODE pipe format); never guess a transcript-to-gene map.
.lamb_tx_or_gene <- function(ids, tpm, stem, tool, platform, say) {
    ids <- as.character(ids)
    if (mean(grepl("|", ids, fixed = TRUE)) > 0.5) {
        parts <- strsplit(ids, "|", fixed = TRUE)
        ensg <- vapply(parts, function(p) { g <- grep("^ENS[A-Z]*G[0-9]", p, value = TRUE); if (length(g)) g[1] else NA_character_ }, "")
        sym <- vapply(parts, function(p) if (length(p) >= 6 && nzchar(p[6])) p[6] else NA_character_, "")
        gene <- if (mean(!is.na(sym)) > 0.9) sym else ensg
        if (all(is.na(gene))) {
            stop(tool, " transcript names have extra fields but no gene field could be found. Give a gene-level ",
                 "file, or summarise transcripts to genes with a tx2gene table (e.g. tximport).", call. = FALSE)
        }
        ok <- !is.na(gene)
        s <- rowsum(tpm[ok], gene[ok], reorder = FALSE)
        say(tool, " values are per transcript; the transcript names carry the gene (GENCODE format), so ",
            format(sum(ok), big.mark = ","), " transcripts were summed into ", format(nrow(s), big.mark = ","), " genes.")
        return(list(expr = matrix(s[, 1], ncol = 1, dimnames = list(rownames(s), stem)),
                    units = "tpm", platform = platform))
    }
    if (mean(grepl(.LAMB_TRANSCRIPT_RE, ids)) > 0.5) {
        stop(tool, " gives one value per transcript (IDs like ", .lamb_examples(ids, 2), "), but LAMB needs one ",
             "value per gene, and the transcript names do not say which gene they belong to. Give a gene-level ",
             "file instead (e.g. Salmon run with --geneMap, which writes quant.genes.sf), or summarise ",
             "transcripts to genes with a tx2gene table (e.g. tximport).", call. = FALSE)
    }
    if (!.lamb_looks_gene(ids)) {
        stop("The names in this ", tool, " file (", .lamb_examples(ids), ") look neither like gene symbols nor ",
             "Ensembl gene IDs.", call. = FALSE)
    }
    list(expr = matrix(tpm, ncol = 1, dimnames = list(ids, stem)), units = "tpm", platform = platform)
}


# ---- a folder of single-sample files ----

.lamb_bulk_folder_files <- function(dir) {
    top <- list.files(dir, full.names = TRUE)
    top <- top[!dir.exists(top)]
    top <- top[grepl(paste0(.LAMB_TEXT_EXT_RE, "|\\.(xlsx|xls)$"), basename(top), ignore.case = TRUE) &
               !grepl("^[.~]|SJ\\.out\\.tab|\\.summary|^Log\\.|\\.log$", basename(top))]
    prefer <- c("quant.genes.sf", "quant.sf", "abundance.tsv")
    sub <- unlist(lapply(list.dirs(dir, recursive = FALSE), function(d) {
        f <- list.files(d)
        hit <- intersect(prefer, f)
        if (length(hit)) return(file.path(d, hit[1]))
        hit <- grep("ReadsPerGene\\.out\\.tab$|\\.genes\\.results$", f, value = TRUE)
        if (length(hit) == 1) file.path(d, hit) else character()
    }))
    c(top, sub)
}

.lamb_read_bulk_folder <- function(dir, say) {
    files <- .lamb_bulk_folder_files(dir)
    if (!length(files)) {
        stop("No expression files were found in the folder ", dir, ". Expected one file per sample ",
             "(e.g. quant.sf, abundance.tsv, ReadsPerGene.out.tab, *.genes.results, featureCounts or a table).",
             call. = FALSE)
    }
    if (length(files) == 1) {
        say("Folder holds one expression file: ", basename(files))
        return(.lamb_read_bulk_file(files, say))
    }
    notes <- character()
    parts <- lapply(files, function(f) {
        local_say <- function(...) notes <<- c(notes, paste0(...))
        tryCatch(.lamb_read_bulk_file(f, local_say),
                 error = function(e) stop("While reading ", basename(dirname(f)), "/", basename(f), ": ",
                                          conditionMessage(e), call. = FALSE))
    })
    plat <- vapply(parts, `[[`, "", "platform")
    if (length(unique(plat)) > 1) {
        tb <- table(plat)
        stop("This folder mixes different kinds of files (", paste0(names(tb), ": ", tb, " file(s)", collapse = "; "),
             "). Put one kind of file per folder.", call. = FALSE)
    }
    units <- unique(vapply(parts, `[[`, "", "units"))
    if (length(units) > 1) stop("The files in this folder hold different kinds of values (",
                                paste(units, collapse = ", "), "). Put one kind per folder.", call. = FALSE)
    say("Folder of ", length(files), " ", plat[1], " files.")
    for (m in unique(notes)) say("  ", m)

    mats <- lapply(parts, function(p) {
        x <- p$expr
        if (anyDuplicated(rownames(x))) {
            if (units %in% c("counts", "tpm", "fpkm")) {
                x <- rowsum(x, rownames(x), reorder = FALSE)
            } else {
                x <- x[order(rowMeans(x, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
                x <- x[!duplicated(rownames(x)), , drop = FALSE]
            }
        }
        x
    })
    all_genes <- unique(unlist(lapply(mats, rownames)))
    common <- Reduce(intersect, lapply(mats, rownames))
    if (!length(common)) stop("The files in this folder share no gene names, so they cannot be combined. ",
                              "Check they were made with the same annotation.", call. = FALSE)
    x <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
    colnames(x) <- .lamb_sample_names(colnames(x), say)
    say("Combined ", length(files), " files into ", ncol(x), " samples: kept ", format(length(common), big.mark = ","),
        " genes found in every file (of ", format(length(all_genes), big.mark = ","), " genes in total).")
    if (length(common) < 0.5 * min(vapply(mats, nrow, 0))) {
        say("Warning: fewer than half of the genes are shared by all files; the files may use different annotations.")
    }
    list(expr = x, units = units, platform = "multi_file", source_platform = plat[1])
}


# ---- checks common to every input ----

.lamb_bulk_final_checks <- function(res, say) {
    x <- res$expr
    if (!is.matrix(x)) x <- as.matrix(x)
    if (!is.double(x)) storage.mode(x) <- "double"
    if (!nrow(x) || !ncol(x)) stop("No expression values were found in the input.", call. = FALSE)
    n_na <- sum(is.na(x))
    if (n_na == length(x)) stop("Every value in the input is missing (empty or NA).", call. = FALSE)
    if (n_na > 0) say(format(n_na, big.mark = ","), " value(s) are missing (empty or NA); they are kept as missing.")
    v <- x[is.finite(x)]
    if (length(v) && all(v < 0)) {
        stop("All values are negative. LAMB needs expression levels (counts, TPM, or log2 expression), ",
             "not values centred on zero such as fold changes or z-scores.", call. = FALSE)
    }
    if (length(v) && all(v == 0)) stop("All values are zero: there is no expression signal in this input.", call. = FALSE)
    if (nrow(x) < 1000) {
        say("Warning: only ", nrow(x), " genes were found. This probably is not a whole-transcriptome file; ",
            "results may be unreliable.")
    }
    res$expr <- x
    if (is.null(res$source_platform)) res$source_platform <- res$platform
    res[c("expr", "units", "platform", "source_platform")]
}

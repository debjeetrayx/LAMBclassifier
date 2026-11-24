# ---- Regenerate inst/extdata/example_bulk.csv and inst/extdata/example_patient.csv ----
# Source (public):
#   GEO series GSE20194 (MAQC-II breast cancer, Affymetrix HG-U133A, GPL96).
#   Popovici V, Chen W, Gallas BG, et al. Effect of training-sample size and
#   classification difficulty on the accuracy of genomic predictors.
#   Breast Cancer Res 2010;12:R5. doi:10.1186/bcr2468
#   (MAQC Consortium. Nat Biotechnol 2010;28:827-838. doi:10.1038/nbt.1665)
#
#   Series matrix: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20194/matrix/GSE20194_series_matrix.txt.gz
#   Platform:      https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/soft/GPL96_family.soft.gz
#
# Terms: GEO data are public. NCBI places no restrictions on their use or
#   distribution (https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html);
#   submitters may hold rights over the underlying data. Cite the paper above.
#
# Processing (no re-normalisation):
#   1. Series-matrix values as deposited (MAS5, log2).
#   2. Probe -> gene: first symbol of the GPL96 "Gene Symbol" column (text
#      before " /// "); probes without a symbol dropped; probes of the same
#      symbol averaged (arithmetic mean), kept to 15 significant digits.
#   3. Classify all 278 profiles with LAMBclassifier(run_pam50 = FALSE).
#   4. Keep six unflagged profiles, one per LAMB class (largest margin in
#      its class, preferring receptor status consistent with the class).
#   5. Keep only genes in lamb_universe(): scoring happens strictly inside the
#      locked universe, so this gives the SAME calls and scores; checked below.
#   6. Values written at 15 significant digits (no further rounding). Do not
#      round more: probe averages that differ only in late digits can become
#      exact ties, which changes singscore ranks and shifts a reported score.
#
# Run from the package root:  Rscript data-raw/example_bulk.R
# Set LAMB_GEO_CACHE to a folder already holding the two GEO files to avoid
# downloading (~110 MB).

devtools::load_all(".", quiet = TRUE)

cache <- Sys.getenv("LAMB_GEO_CACHE", file.path(tempdir(), "geo"))
dir.create(cache, showWarnings = FALSE, recursive = TRUE)
files <- c(
    series = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20194/matrix/GSE20194_series_matrix.txt.gz",
    gpl    = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/soft/GPL96_family.soft.gz"
)
local <- c(series = file.path(cache, "GSE20194_series_matrix.txt.gz"),
           gpl = if (file.exists(file.path(cache, "GPL96.soft.gz"))) file.path(cache, "GPL96.soft.gz")
                 else file.path(cache, "GPL96_family.soft.gz"))
for (k in names(files)) if (!file.exists(local[[k]])) {
    options(timeout = max(600, getOption("timeout")))
    utils::download.file(files[[k]], local[[k]], mode = "wb")
}

# ---- 1-2. probe -> gene ----
sm <- utils::read.delim(gzfile(local[["series"]]), comment.char = "!", check.names = FALSE)
gl <- readLines(gzfile(local[["gpl"]]))
tab <- utils::read.delim(text = gl[(which(gl == "!platform_table_begin") + 1):(which(gl == "!platform_table_end") - 1)],
                         check.names = FALSE, quote = "", stringsAsFactors = FALSE)
sym <- trimws(sub(" ///.*", "", tab[["Gene Symbol"]][match(sm[[1]], tab$ID)]))
ok <- !is.na(sym) & nzchar(sym)
v <- as.matrix(sm[ok, -1])
expr <- rowsum(v, sym[ok]) / as.vector(table(sym[ok])[sort(unique(sym[ok]))])
# 15 significant digits (write.table default).
expr <- signif(expr, 15)
stopifnot(nrow(expr) == 13237, ncol(expr) == 278)

# ---- 3. classify the full cohort ----
full <- as.data.frame(LAMBclassifier(expr, run_pam50 = FALSE))

# ---- 4. six profiles, one per class ----
# Fixed here so the files are stable (unflagged, largest margin per class).
ids <- c(LHS = "GSM505471", "MA-Lo" = "GSM505546", Apocrine = "GSM505331",
         "Basal-like" = "GSM505423", Squamous = "GSM505491", Chondroid = "GSM505558")
ref <- full[match(ids, full$Sample), ]
stopifnot(identical(ref$Pred_Class, names(ids)), all(ref$Flag == ""))

# ---- 5-6. universe genes only, full precision; write; verify identical ----
genes <- intersect(rownames(expr), lamb_universe())
small <- expr[genes, ids]

write_tab <- function(m, path) {
    out <- data.frame(Gene = rownames(m), m, check.names = FALSE)
    utils::write.csv(out, path, row.names = FALSE, quote = FALSE)
}
dir.create("inst/extdata", showWarnings = FALSE)
write_tab(small, "inst/extdata/example_bulk.csv")
write_tab(small[, ids[["LHS"]], drop = FALSE], "inst/extdata/example_patient.csv")

# Read the files back the way lamb_run() does and compare with the full-matrix run.
num <- grep("^Score_|^Margin$", names(ref), value = TRUE)
same <- function(a, b) identical(a$Pred_Class, b$Pred_Class) && identical(a$Flag, b$Flag) &&
    max(abs(as.matrix(a[num]) - as.matrix(b[num]))) == 0
back <- function(f) as.matrix(utils::read.table(f, sep = ",", header = TRUE, row.names = 1, check.names = FALSE))
chk <- as.data.frame(LAMBclassifier(back("inst/extdata/example_bulk.csv"), run_pam50 = FALSE))
one <- as.data.frame(LAMBclassifier(back("inst/extdata/example_patient.csv"), run_pam50 = FALSE))
stopifnot(same(chk, ref), same(one, ref[1, ]))
run <- lamb_run("inst/extdata/example_bulk.csv", out_dir = file.path(tempdir(), "chk"), run_pam50 = FALSE)
stopifnot(same(as.data.frame(run$details), ref))

print(ref[, c("Sample", "Pred_Class", "Margin", "Flag")])
cat(length(genes), "genes x", length(ids), "samples written; calls and scores identical to the full 13,237-gene run\n")

# ---- Regenerate inst/extdata/example_10x/  (a REAL public single-cell subset) ----
# Source (public):
#   GEO sample GSM5354516 (patient CID3941, ER+ primary breast cancer, 10x
#   Chromium 3'), part of series GSE176078.
#   Wu SZ, Al-Eryani G, Roden DL, et al. A single-cell and spatially resolved
#   atlas of human breast cancers. Nat Genet 2021;53:1334-1347.
#   doi:10.1038/s41588-021-00911-1 (open access, CC BY 4.0).
#
#   File: https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354516/suppl/GSM5354516_CID3941.tar.gz
#   (2.4 MB; raw UMI counts, gene symbols, barcodes, and the authors' per-cell
#   metadata.csv with cell-type annotations)
#
# Terms: GEO data are public. NCBI places no restrictions on their use or
#   distribution (https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html).
#   Cite Wu et al. 2021 when using this example.
#
# Subset (exact):
#   * Cells: 120 of the 196 "Cancer Epithelial" cells, 60 of the 266 "T-cells"
#     and all 55 "B-cells" (celltype_major in the authors' metadata), sampled
#     with set.seed(20211)  -> 235 cells.
#   * Genes: those in lamb_universe() plus the epithelial-check genes used by
#     lamb_run() (EPCAM, CDH1, KRT5/7/8/14/17/18/19, CLDN4/7, ELF3, GRHL2,
#     ESRP1, ERBB3, SPINT2) that are present. Raw UMI counts unchanged.
#   * Written in Cell Ranger (v3+) layout: matrix.mtx.gz, features.tsv.gz,
#     barcodes.tsv.gz. features.tsv.gz column 1 is the Ensembl ID where the
#     package's gene table maps the symbol, otherwise the symbol itself.
#   * The "clusters" in analysis/clustering/gene_expression_graphclust/
#     clusters.csv are the AUTHORS' cell-type labels (celltype_major), not a
#     Cell Ranger graph-based clustering; the file layout is Cell Ranger's so
#     that lamb_run() picks it up.
#
# Run from the package root:  Rscript data-raw/example_10x.R

devtools::load_all(".", quiet = TRUE)

url <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5354nnn/GSM5354516/suppl/GSM5354516_CID3941.tar.gz"
work <- file.path(tempdir(), "wu2021")
dir.create(work, showWarnings = FALSE)
tgz <- file.path(work, basename(url))
if (!file.exists(tgz)) utils::download.file(url, tgz, mode = "wb")
utils::untar(tgz, exdir = work)
src <- file.path(work, "CID3941")

m <- Matrix::readMM(file.path(src, "count_matrix_sparse.mtx"))
genes <- readLines(file.path(src, "count_matrix_genes.tsv"))
bc <- readLines(file.path(src, "count_matrix_barcodes.tsv"))
meta <- utils::read.csv(file.path(src, "metadata.csv"), row.names = 1, stringsAsFactors = FALSE)
dimnames(m) <- list(genes, bc)
m <- methods::as(m, "CsparseMatrix")
stopifnot(all(m@x == round(m@x)), all(bc %in% rownames(meta)))

set.seed(20211)
type <- meta[bc, "celltype_major"]
take <- c("Cancer Epithelial" = 120, "T-cells" = 60, "B-cells" = 55)
cells <- unlist(lapply(names(take), function(k) sort(sample(bc[type == k], take[[k]]))))

keep <- union(intersect(lamb_universe(), genes), intersect(.LAMB_EPITHELIAL, genes))
keep <- genes[genes %in% keep]          # original gene order
sub <- m[keep, cells]

out <- "inst/extdata/example_10x"
unlink(out, recursive = TRUE)
dir.create(file.path(out, "analysis", "clustering", "gene_expression_graphclust"), recursive = TRUE)
mtx <- file.path(out, "matrix.mtx")
Matrix::writeMM(sub, mtx)
gzip_file <- function(f) { con <- gzfile(paste0(f, ".gz"), "wb", compression = 9); writeLines(readLines(f), con); close(con); unlink(f) }
gzip_file(mtx)
gl <- get_gene_length()
ens <- gl$gene_id[match(keep, gl$SYMBOL)]
ens[is.na(ens) | !nzchar(ens)] <- keep[is.na(ens) | !nzchar(ens)]
con <- gzfile(file.path(out, "features.tsv.gz"), "wb"); writeLines(paste(ens, keep, "Gene Expression", sep = "\t"), con); close(con)
con <- gzfile(file.path(out, "barcodes.tsv.gz"), "wb"); writeLines(cells, con); close(con)
utils::write.csv(data.frame(Barcode = cells, Cluster = meta[cells, "celltype_major"]),
                 file.path(out, "analysis", "clustering", "gene_expression_graphclust", "clusters.csv"),
                 row.names = FALSE)

res <- lamb_run(out, out_dir = file.path(tempdir(), "chk10x"))
print(res$calls)
stopifnot(res$calls$LAMB_class[res$calls$Cluster == "cluster Cancer Epithelial"] == "LHS")
cat(length(cells), "cells x", length(keep), "genes; nonzero entries:", Matrix::nnzero(sub), "\n")

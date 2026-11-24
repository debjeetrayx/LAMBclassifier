# ---- Regenerate inst/extdata/xenium_sim/   (SIMULATED data) ----
# THIS IS SIMULATED DATA. No tissue was measured. It exists only so that
# lamb_run() can be tried on a Xenium-format folder without downloading a
# multi-gigabyte public bundle.
#
# Derived from public data: the expression profiles in
# inst/extdata/example_bulk.csv (GEO GSE20194, MAQC-II; see
# data-raw/example_bulk.R for accession, citation and terms). Run that script
# first.
#
# Recipe (set.seed(20194)):
#   * Panel: every LAMB signature gene and epithelial-check gene present in
#     example_bulk.csv, plus 300 further lamb_universe() genes drawn at random,
#     plus 5 negative-control features (dropped by the reader, as on Xenium).
#   * Three epithelial groups of 110 cells, each Poisson-sampled from one
#     public profile: GSM505471 (LHS), GSM505423 (Basal-like),
#     GSM505546 (MA-Lo). Per-gene rate proportional to 2^log2 value,
#     per-cell depth drawn from a log-normal around 300 transcripts.
#   * One non-epithelial group of 120 cells: the median of the six example
#     profiles with all epithelial-check genes and LAMB signature genes set
#     to near-zero, and the immune gate genes (PTPRC, CD3E, ...) set among the
#     highest, so the cell-type gate leaves it out as immune.
#   * Positions: each epithelial group occupies its own elliptical region of a
#     2 x 1.5 mm field; the non-epithelial cells fill the space between.
#   * Clusters (analysis/clustering/gene_expression_graphclust/clusters.csv)
#     are the simulated groups, numbered 1-4 like Xenium graph clusters.
#   * experiment.xenium carries "pixel_size": 0.2125 and states it is simulated;
#     SIMULATED_DATA.txt in the folder says so too. The folder name is short
#     because R package paths must stay under 100 bytes.
#
# Run from the package root:  Rscript data-raw/xenium_sim.R

devtools::load_all(".", quiet = TRUE)
set.seed(20194)

bulk <- as.matrix(utils::read.csv("inst/extdata/example_bulk.csv", row.names = 1, check.names = FALSE))
sig <- unique(unlist(default_gene_sets()))
epi <- intersect(.LAMB_EPITHELIAL, rownames(bulk))
imm <- intersect(.LAMB_GATE_MODULES$immune, rownames(bulk))   # the stand-in group is immune-like
amp <- intersect(DEFAULT_HER2_AMPLICON_GENES, rownames(bulk))   # so the HER2-amplicon modifier can be shown
core <- intersect(c(sig, epi, imm, amp), rownames(bulk))
extra <- sample(setdiff(intersect(lamb_universe(), rownames(bulk)), core), 300)
panel <- sort(c(core, extra))

groups <- list(
    list(id = "GSM505471", n = 110, centre = c(500, 450), axes = c(380, 300)),   # LHS
    list(id = "GSM505423", n = 110, centre = c(1500, 450), axes = c(380, 300)),  # Basal-like
    list(id = "GSM505546", n = 110, centre = c(1000, 1150), axes = c(450, 250))  # MA-Lo
)
stroma <- apply(bulk[panel, ], 1, stats::median)
stroma[intersect(c(sig, epi), panel)] <- -6            # ~1/2^20 of a typical gene
stroma[imm] <- max(stroma) + 2                         # immune genes among the highest, as in immune cells

draw_cells <- function(log2_profile, n) {
    rate <- 2^log2_profile
    rate <- rate / sum(rate)
    depth <- round(stats::rlnorm(n, log(300), 0.35))
    vapply(depth, function(d) stats::rpois(length(rate), d * rate), numeric(length(rate)))
}
in_ellipse <- function(xy, c0, ax) ((xy[, 1] - c0[1]) / ax[1])^2 + ((xy[, 2] - c0[2]) / ax[2])^2 <= 1
place <- function(n, accept) {
    out <- matrix(numeric(), 0, 2)
    while (nrow(out) < n) {
        p <- cbind(stats::runif(4 * n, 0, 2000), stats::runif(4 * n, 0, 1500))
        out <- rbind(out, p[accept(p), , drop = FALSE])
    }
    out[seq_len(n), , drop = FALSE]
}

counts <- list(); xy <- list(); cl <- integer()
for (k in seq_along(groups)) {
    g <- groups[[k]]
    counts[[k]] <- draw_cells(bulk[panel, g$id], g$n)
    xy[[k]] <- place(g$n, function(p) in_ellipse(p, g$centre, g$axes))
    cl <- c(cl, rep(k, g$n))
}
outside <- function(p) !Reduce(`|`, lapply(groups, function(g) in_ellipse(p, g$centre, g$axes * 1.1)))
counts[[4]] <- draw_cells(stroma, 120)
xy[[4]] <- place(120, outside)
cl <- c(cl, rep(4L, 120))
m <- do.call(cbind, counts)
xy <- do.call(rbind, xy)
rownames(m) <- panel

# 5 negative-control features with sparse background counts
neg <- c(paste0("NegControlProbe_0000", 1:3), paste0("NegControlCodeword_050", 0:1))
negm <- matrix(stats::rpois(length(neg) * ncol(m), 0.02), length(neg))
full <- rbind(m, negm)

# Xenium-style cell ids: 8 letters a-p + "-1"
ids <- unique(replicate(ncol(full) * 2, paste0(paste(sample(letters[1:16], 8, TRUE), collapse = ""), "-1")))
ids <- ids[seq_len(ncol(full))]
o <- order(xy[, 1], xy[, 2])           # instruments write cells in no class order
full <- full[, o]; xy <- xy[o, ]; cl <- cl[o]
colnames(full) <- ids

out <- "inst/extdata/xenium_sim"
unlink(out, recursive = TRUE)
dir.create(file.path(out, "cell_feature_matrix"), recursive = TRUE)
dir.create(file.path(out, "analysis", "clustering", "gene_expression_graphclust"), recursive = TRUE)
gz_lines <- function(x, f) { con <- gzfile(f, "wb", compression = 9); writeLines(x, con); close(con) }
tmp <- tempfile(fileext = ".mtx")
Matrix::writeMM(methods::as(Matrix::Matrix(full, sparse = TRUE), "CsparseMatrix"), tmp)
gz_lines(readLines(tmp), file.path(out, "cell_feature_matrix", "matrix.mtx.gz"))
gl <- get_gene_length()
ens <- gl$gene_id[match(panel, gl$SYMBOL)]
ens[is.na(ens) | !nzchar(ens)] <- panel[is.na(ens) | !nzchar(ens)]
gz_lines(c(paste(ens, panel, "Gene Expression", sep = "\t"),
           paste(neg, neg, ifelse(grepl("Probe", neg), "Negative Control Probe", "Negative Control Codeword"), sep = "\t")),
         file.path(out, "cell_feature_matrix", "features.tsv.gz"))
gz_lines(ids, file.path(out, "cell_feature_matrix", "barcodes.tsv.gz"))
cells <- data.frame(cell_id = ids, x_centroid = round(xy[, 1], 2), y_centroid = round(xy[, 2], 2),
                    transcript_counts = colSums(m[, o]), control_probe_counts = colSums(negm[, o]))
con <- gzfile(file.path(out, "cells.csv.gz"), "wb"); utils::write.csv(cells, con, row.names = FALSE); close(con)
utils::write.csv(data.frame(Barcode = ids, Cluster = cl),
                 file.path(out, "analysis", "clustering", "gene_expression_graphclust", "clusters.csv"), row.names = FALSE)
writeLines(c('{',
             '    "note": "SIMULATED example for LAMBclassifier. Not a real Xenium run; no tissue was measured. See data-raw/xenium_sim.R.",',
             '    "run_name": "LAMB simulated example",',
             '    "panel_name": "simulated LAMB example panel",',
             paste0('    "panel_num_targets_predesigned": ', length(panel), ','),
             '    "pixel_size": 0.2125',
             '}'), file.path(out, "experiment.xenium"))

writeLines(c("SIMULATED DATA - no tissue was measured.",
             "Cells were Poisson-sampled from public GEO GSE20194 bulk profiles",
             "(GSM505471 LHS, GSM505423 Basal-like, GSM505546 MA-Lo) plus a",
             "non-epithelial stand-in group, placed in separate regions of a 2 x 1.5 mm field.",
             "For trying LAMBclassifier::lamb_run() on a Xenium-format folder only.",
             "Recipe: data-raw/xenium_sim.R in the package source. See ../README.md."),
           file.path(out, "SIMULATED_DATA.txt"))

res <- lamb_run(out, out_dir = file.path(tempdir(), "chkxen"))
print(res$calls)
want <- c("cluster 1" = "LHS", "cluster 2" = "Basal-like", "cluster 3" = "MA-Lo")
got <- stats::setNames(res$calls$LAMB_class, res$calls$Cluster)
stopifnot(identical(unname(got[names(want)]), unname(want)),
          res$calls$Cell_type[res$calls$Cluster == "cluster 4"] == "immune")
cat(ncol(full), "simulated cells x", length(panel), "panel genes\n")

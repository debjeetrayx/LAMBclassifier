# Regenerates the example report pages in dev/report_examples/ from PUBLIC data
# only (GSE20194, MAQC-II; and the package's inst/extdata examples).
# Run from the package root:
#   LAMB_GSE20194=<path to GSE20194_matrix.txt> LAMB_EXAMPLE_OUT=dev/report_examples \
#     Rscript dev/report_examples/make_examples.R
# then turn each PDF page into a PNG. The committed PNGs were made with pdfseparate +
# macOS `sips` because pdftoppm on the build machine substitutes a regular face for the
# non-embedded Helvetica-Bold (bold text renders as plain), which misrepresents the pages.
suppressMessages(devtools::load_all(".", quiet = TRUE))
out <- Sys.getenv("LAMB_EXAMPLE_OUT", tempdir())
cache <- Sys.getenv("LAMB_EXAMPLE_CACHE", tempdir())
keep_attrs <- function(x, from) { for (a in c("coverage", "specification")) attr(x, a) <- attr(from, a); x }

gse <- Sys.getenv("LAMB_GSE20194", "GSE20194_matrix.txt")
rds <- file.path(cache, "gse20194_res.rds")
if (!file.exists(rds)) {
  x <- as.matrix(read.table(gse, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE,
                            quote = "", comment.char = ""))
  res <- suppressWarnings(LAMBclassifier(x, run_pam50 = TRUE))
  # NC: every squamous gene removed
  x_nc <- x[!rownames(x) %in% DEFAULT_SQUAMOUS_EPITHELIAL_GENES, c("GSM505327", "GSM505333")]
  res_nc <- suppressWarnings(LAMBclassifier(x_nc, run_pam50 = FALSE))
  # low coverage: keep only 4 squamous and 4 MA-Lo genes (still scored, below 60%)
  drop <- c(DEFAULT_SQUAMOUS_EPITHELIAL_GENES[-(1:4)], DEFAULT_MALO_GENES[-(1:4)])
  x_low <- x[!rownames(x) %in% drop, c("GSM505327", "GSM505335", "GSM505417")]
  res_low <- suppressWarnings(LAMBclassifier(x_low, run_pam50 = FALSE))
  saveRDS(list(res = res, nc = res_nc, low = res_low), rds)
}
d <- readRDS(rds)
res <- d$res
print(table(res$Pred_Class))

# cohort: cover + one page per sample
keep <- c("GSM505327", "GSM505333", "GSM505335", "GSM505331", "GSM505417", "GSM505343", "GSM505534")
sub <- keep_attrs(res[match(keep, res$Sample), ], res)
print(as.data.frame(sub[, c("Sample", "Pred_Class", "Margin", "Flag", "PAM50")]))
lamb_report(sub, file.path(out, "bulk_cohort.pdf"), title = "GSE20194 example cohort")
lamb_report(keep_attrs(res[1, ], res), file.path(out, "bulk_single.pdf"))
lamb_report(keep_attrs(d$nc[1, ], d$nc), file.path(out, "bulk_nc.pdf"))
print(as.data.frame(d$low[, c("Sample", "Pred_Class", "Margin", "Flag")]))
lamb_report(d$low, file.path(out, "bulk_lowcov.pdf"), title = "Low gene coverage example")

# cells: cluster table built from bulk samples standing in for clusters
set.seed(1)
det <- keep_attrs(res[1:14, ], res)
det$Sample <- paste("cluster", 1:14)
calls <- data.frame(Cluster = c(det$Sample, "cluster 15", "cluster 16"),
                    Cells = c(round(runif(14, 100, 5000)), 30, 800),
                    Cell_type = c(rep("epithelial", 15), "not epithelial"), Epithelial_share = 0.9,
                    LAMB_class = c(det$Pred_Class, "Too few cells", "LHS"), Note = "")
lamb_report_cells(calls, det, file.path(out, "cells_many.pdf"))
det1 <- keep_attrs(d$low[3, ], d$low); det1$Sample <- "cluster 1"
lamb_report_cells(calls[c(1, 16), ], det1, file.path(out, "cells_one.pdf"))

# slide overview from the simulated Xenium example
xen <- system.file("extdata", "xenium_sim", package = "LAMBclassifier")
if (!nzchar(xen)) xen <- "inst/extdata/xenium_sim"
xo <- file.path(tempdir(), "xen_run")
suppressMessages(lamb_run(xen, out_dir = xo))
file.copy(file.path(xo, "LAMB_summary.pdf"), file.path(out, "slide_summary.pdf"), overwrite = TRUE)

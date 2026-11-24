# ---- Update LAMB_CALIBRATION in R/sysdata.rda to classifier specification 2.0.0 ----
# Edits only the calibration list; gene.length and pam50.robust are written back
# unchanged, and the 13,468-gene universe is kept as it was (same hash).
#
#   version            "1.0.0" -> "2.0.0"
#   signature_of_class "MA-Lo (HER2)" = "malo_her2"  ->  "MA-Lo" = "malo"
#   her2_amplicon      new: module genes, q = 9, min_high = 3, min_measured = 3
#
# min_margin (0.02) and min_gene_coverage (0.60) are unchanged.
#
# Run from the package root:  Rscript data-raw/calibration_2.0.0.R

e <- new.env()
load("R/sysdata.rda", envir = e)
cal <- e$LAMB_CALIBRATION
universe_before <- cal$universe

cal$version <- "2.0.0"
cal$signature_of_class <- c("Chondroid" = "chondroid", "Squamous" = "squamous",
                            "Basal-like" = "basal_like", "LHS" = "luminal",
                            "Apocrine" = "apocrine", "MA-Lo" = "malo")
cal$her2_amplicon <- list(genes = c("ERBB2", "GRB7", "MIEN1", "PGAP3", "STARD3"),
                          q = 9, min_high = 3L, min_measured = 3L)

stopifnot(identical(cal$universe, universe_before),
          identical(cal$min_margin, 0.02), identical(cal$min_gene_coverage, 0.6))
e$LAMB_CALIBRATION <- cal
save(list = c("gene.length", "LAMB_CALIBRATION", "pam50.robust"), envir = e,
     file = "R/sysdata.rda", compress = "xz")

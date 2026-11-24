# Evidence for the 0.02 margin and the 0.60 coverage flag, from public SCAN-B data.
# Usage: Rscript threshold_evidence.R <folder with GSE96058_matrix.txt>
# Writes threshold_evidence.tsv next to this script's working directory.
suppressMessages(library(LAMBclassifier))
args <- commandArgs(trailingOnly = TRUE)
data_dir <- if (length(args)) args[1] else "."
qs <- c(0.5, 0.9, 0.95, 0.99)
rows <- list()

# 1. Technical repeatability: 136 SCAN-B tumours sequenced twice.
cl <- read.csv(system.file("validation", "GSE96058_calls.csv", package = "LAMBclassifier"))
rownames(cl) <- cl$Sample
rp <- grep("repl$", cl$Sample, value = TRUE)
pr <- sub("repl$", "", rp)
keep <- pr %in% cl$Sample
rp <- rp[keep]; pr <- pr[keep]
d_margin <- abs(cl[rp, "Margin"] - cl[pr, "Margin"])
d_lean <- abs((cl[rp, "Score_LHS"] - cl[rp, "Score_Apocrine"]) - (cl[pr, "Score_LHS"] - cl[pr, "Score_Apocrine"]))
for (q in qs) {
    rows[[length(rows) + 1]] <- data.frame(Analysis = "technical repeat", Quantity = "abs change in margin",
                                           Group = paste0("quantile ", q), N = length(rp), Value = unname(quantile(d_margin, q)))
    rows[[length(rows) + 1]] <- data.frame(Analysis = "technical repeat", Quantity = "abs change in LHS minus Apocrine",
                                           Group = paste0("quantile ", q), N = length(rp), Value = unname(quantile(d_lean, q)))
}

# 2. Gene loss: the same random genes removed from every sample, as a platform would lack them.
x <- utils::read.delim(file.path(data_dir, "GSE96058_matrix.txt"), row.names = 1, check.names = FALSE)
x <- as.matrix(x[, !grepl("repl$", colnames(x))])
set.seed(20261124)
x <- x[intersect(rownames(x), lamb_universe()), sample(ncol(x), 600)]
full <- suppressMessages(suppressWarnings(LAMBclassifier(x, run_pam50 = FALSE)))
sim <- list()
for (f in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) for (s in 1:3) {
    set.seed(1000 * f + s)
    drop <- sample(rownames(x), round(f * nrow(x)))
    r <- suppressMessages(suppressWarnings(LAMBclassifier(x[setdiff(rownames(x), drop), ], run_pam50 = FALSE)))
    sc <- attr(r, "sample_coverage")
    dec <- intersect(c("chondroid", "squamous", "basal_like", "luminal", "apocrine", "malo"), colnames(sc))
    sim[[length(sim) + 1]] <- data.frame(call = r$Pred_Class, ref = full$Pred_Class,
                                         min_cov = apply(as.matrix(sc[, dec, drop = FALSE]), 1, min),
                                         flagged = grepl("low_coverage", r$Flag))
}
sim <- do.call(rbind, sim)
# Three outcomes per call: same class, a different class, or refused (NC).
outcome <- ifelse(sim$call == "NC", "refused (NC)", ifelse(sim$call == sim$ref, "same class", "different class"))
band <- cut(sim$min_cov, c(0, 0.3, 0.5, 0.6, 0.7, 1.001), right = FALSE,
            labels = c("<0.30", "0.30-0.49", "0.50-0.59", "0.60-0.69", ">=0.70"))
n_note <- "600 tumours x 21 gene-loss runs"
for (b in levels(band)) for (o in c("same class", "different class", "refused (NC)")) {
    k <- band == b
    rows[[length(rows) + 1]] <- data.frame(Analysis = "gene loss", Quantity = paste("share", o),
                                           Group = paste("lowest class-signature coverage", b),
                                           N = sum(k), Value = mean(outcome[k] == o))
}
scored <- sim$call != "NC"
for (g in c(FALSE, TRUE)) {
    k <- scored & sim$flagged == g
    rows[[length(rows) + 1]] <- data.frame(Analysis = "gene loss", Quantity = "share same class (scored calls)",
                                           Group = if (g) "low_coverage flagged" else "not flagged",
                                           N = sum(k), Value = mean(sim$call[k] == sim$ref[k]))
}
out <- do.call(rbind, rows)
out$Value <- round(out$Value, 4)
out$Note <- ifelse(out$Analysis == "gene loss", n_note, "136 tumour pairs")
utils::write.table(out, "threshold_evidence.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)

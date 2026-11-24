# A name defined in two files is silent in R: the file collated last wins.
test_that("no top-level name is defined in two R/ files", {
    dir <- test_path("..", "..", "R")
    skip_if_not(dir.exists(dir), "R/ sources not available (installed package)")
    defs <- do.call(rbind, lapply(list.files(dir, pattern = "[.]R$", full.names = TRUE), function(f) {
        l <- readLines(f, warn = FALSE)
        n <- sub("^([.A-Za-z_][A-Za-z0-9._]*)\\s*<-.*", "\\1", grep("^[.A-Za-z_][A-Za-z0-9._]*\\s*<-", l, value = TRUE))
        if (length(n)) data.frame(name = n, file = basename(f)) else NULL
    }))
    dup <- defs[defs$name %in% defs$name[duplicated(defs$name)], ]
    expect_equal(nrow(dup), 0, info = paste(dup$name, dup$file, collapse = "; "))
})

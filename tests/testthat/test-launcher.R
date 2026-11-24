withr_tempdir <- function() {
    d <- tempfile("launcher_")
    dir.create(d)
    d
}

test_that("macOS/Linux launcher is written, executable, and calls this R", {
    dir <- withr_tempdir()
    path <- suppressMessages(lamb_install_launcher(dir = dir, platform = "unix"))
    expect_equal(basename(path), "LAMB classifier.command")
    txt <- readLines(path)
    expect_equal(txt[1], "#!/bin/bash")
    expect_false(any(grepl("{{RSCRIPT}}", txt, fixed = TRUE)))
    expect_true(any(grepl(file.path(R.home("bin"), "Rscript"), txt, fixed = TRUE)))
    expect_true(any(grepl("library(LAMBclassifier)", txt, fixed = TRUE)))
    expect_true(any(grepl("export LAMB_LAUNCHER=1", txt, fixed = TRUE)))
    expect_true(any(grepl("Press any key", txt)))          # window stays open
    skip_on_os("windows")
    mode <- as.integer(file.info(path)$mode)
    expect_true(bitwAnd(mode, strtoi("100", 8L)) > 0)       # owner execute bit
    expect_identical(unname(file.access(path, 1)), 0L)
})

test_that("Windows launcher is written with CRLF and accepts a dropped file", {
    dir <- withr_tempdir()
    path <- suppressMessages(lamb_install_launcher(dir = dir, platform = "windows"))
    expect_equal(basename(path), "LAMB classifier.bat")
    raw <- readBin(path, "raw", file.info(path)$size)
    expect_true(grepl("\r\n", rawToChar(raw), fixed = TRUE))
    txt <- readLines(path)
    expect_true(any(grepl("Rscript.exe", txt, fixed = TRUE)))
    expect_true(any(grepl("%~1", txt, fixed = TRUE)))
    expect_true(any(grepl("LAMB_INPUT", txt, fixed = TRUE)))
    expect_true(any(grepl("^pause", txt)))
    expect_false(any(grepl("{{RSCRIPT}}", txt, fixed = TRUE)))
})

test_that("an existing launcher is never overwritten silently", {
    dir <- withr_tempdir()
    path <- suppressMessages(lamb_install_launcher(dir = dir, platform = "unix"))
    writeLines("mine", path)
    expect_error(lamb_install_launcher(dir = dir, platform = "unix"), "already exists")
    expect_equal(readLines(path), "mine")
    suppressMessages(lamb_install_launcher(dir = dir, platform = "unix", overwrite = TRUE))
    expect_equal(readLines(path)[1], "#!/bin/bash")
})

test_that("missing folder is a plain error", {
    expect_error(lamb_install_launcher(dir = file.path(tempdir(), "no_such_folder_xyz")), "Folder not found")
})

test_that("the file picker refuses plainly when no one can answer it", {
    skip_if(interactive())
    old <- Sys.getenv("LAMB_LAUNCHER", NA)
    Sys.unsetenv("LAMB_LAUNCHER")
    on.exit(if (!is.na(old)) Sys.setenv(LAMB_LAUNCHER = old))
    expect_error(lamb_choose_input("file"), "cannot show a file picker")
    expect_error(lamb_choose_input("folder"), "cannot show a file picker")
})

test_that("lamb_open does nothing outside an interactive session", {
    skip_if(interactive())
    old <- Sys.getenv("LAMB_LAUNCHER", NA)
    Sys.unsetenv("LAMB_LAUNCHER")
    on.exit(if (!is.na(old)) Sys.setenv(LAMB_LAUNCHER = old))
    expect_false(lamb_open(tempdir()))
})

# ---- No-code launcher: file picker, opening results, desktop icon ----

# The desktop launcher sets LAMB_LAUNCHER=1: Rscript is never interactive().
.lamb_can_ask <- function() interactive() || identical(Sys.getenv("LAMB_LAUNCHER"), "1")

# Returns the chosen path, or NULL when cancelled.
lamb_choose_input <- function(type = c("file", "folder"),
                              caption = if (type == "file") "Choose your expression data" else
                                  "Choose a folder (10x, Xenium, or one file per patient)") {
    type <- match.arg(type)
    if (!.lamb_can_ask()) {
        stop("No input given, and this R session cannot show a file picker (it is not interactive).\n",
             "Give the path instead, for example: lamb_run(\"my_data.csv\")", call. = FALSE)
    }
    path <- .lamb_pick(type, caption)
    if (is.null(path) || !length(path) || is.na(path[1]) || !nzchar(path[1])) return(NULL)
    normalizePath(path[1], mustWork = FALSE)
}

.lamb_pick <- function(type, caption) {
    if (identical(Sys.getenv("RSTUDIO"), "1") && requireNamespace("rstudioapi", quietly = TRUE) &&
        isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
        return(if (type == "file") rstudioapi::selectFile(caption = caption) else
            rstudioapi::selectDirectory(caption = caption))
    }
    if (identical(Sys.info()[["sysname"]], "Darwin") && nzchar(Sys.which("osascript"))) {
        # Native Finder dialog; tcltk on macOS needs XQuartz, which most people lack.
        what <- if (type == "file") "choose file" else "choose folder"
        cmd <- sprintf("POSIX path of (%s with prompt \"%s\")", what, gsub("\"", "'", caption))
        out <- suppressWarnings(system2("osascript", c("-e", shQuote(cmd)), stdout = TRUE, stderr = FALSE))
        if (!is.null(attr(out, "status")) || !length(out)) return(NULL)   # cancelled
        return(sub("/$", "", out[1]))
    }
    if (.Platform$OS.type == "windows") {
        return(if (type == "file") utils::choose.files(caption = caption, multi = FALSE) else
            utils::choose.dir(caption = caption))
    }
    if (isTRUE(capabilities("tcltk")) && requireNamespace("tcltk", quietly = TRUE)) {
        out <- tryCatch(
            if (type == "file") tcltk::tk_choose.files(caption = caption, multi = FALSE) else
                tcltk::tk_choose.dir(caption = caption),
            error = function(e) e)
        if (!inherits(out, "error")) return(if (length(out)) out else NULL)
    }
    prompt <- paste0(caption, " - type or paste the full path (empty to cancel): ")
    ans <- if (interactive()) readline(prompt) else {
        cat(prompt)
        con <- file("stdin"); on.exit(close(con), add = TRUE)
        x <- readLines(con, n = 1, warn = FALSE)
        if (length(x)) x else ""
    }
    ans <- gsub("^[\"']|[\"']$", "", trimws(ans))
    if (nzchar(ans)) ans else NULL
}

# Native dialogs pick a file or a folder, not both, so ask first.
.lamb_ask_kind <- function() {
    q <- "Is your data a single file, or a folder (10x, Xenium, Visium, CosMx, MERSCOPE, or one file per patient)?"
    if (identical(Sys.info()[["sysname"]], "Darwin") && nzchar(Sys.which("osascript"))) {
        cmd <- sprintf("button returned of (display dialog \"%s\" buttons {\"Cancel\", \"A folder\", \"A file\"} default button \"A file\" with title \"LAMB classifier\")", q)
        out <- suppressWarnings(system2("osascript", c("-e", shQuote(cmd)), stdout = TRUE, stderr = FALSE))
        if (!is.null(attr(out, "status")) || !length(out)) return(NULL)
        return(if (out[1] == "A folder") "folder" else "file")
    }
    if (.Platform$OS.type == "windows") {
        a <- utils::winDialog("yesnocancel", paste(q, "\n\nYes = a file, No = a folder"))
        return(switch(a, YES = "file", NO = "folder", NULL))
    }
    if (interactive()) {
        a <- utils::menu(c("A file", "A folder"), title = q)
        return(switch(a + 1, NULL, "file", "folder"))
    }
    cat(q, " [f = file, d = folder, empty = cancel]: ")
    con <- file("stdin"); on.exit(close(con), add = TRUE)
    a <- tolower(trimws(readLines(con, n = 1, warn = FALSE)))
    if (!length(a) || !nzchar(a)) return(NULL)
    if (startsWith(a, "d")) "folder" else "file"
}

# No-op outside an interactive session or the launcher.
lamb_open <- function(path) {
    if (!.lamb_can_ask() || is.null(path) || !file.exists(path)) return(invisible(FALSE))
    path <- normalizePath(path)
    sys <- Sys.info()[["sysname"]]
    ok <- tryCatch({
        if (.Platform$OS.type == "windows") {
            shell.exec <- get("shell.exec", envir = baseenv())
            shell.exec(path)
        } else if (identical(sys, "Darwin")) {
            system2("open", shQuote(path), wait = FALSE)
        } else {
            system2("xdg-open", shQuote(path), wait = FALSE, stdout = FALSE, stderr = FALSE)
        }
        TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    invisible(ok)
}

#' Put a double-click LAMB classifier icon on the desktop
#'
#' @description
#' Writes a small launcher file so that LAMB can be run without typing any R:
#' double-click it, choose your data in the window that opens, and the results
#' are written next to the data (see \code{\link{lamb_run}}).
#'
#' \itemize{
#'   \item macOS and Linux: \file{LAMB classifier.command}, made executable.
#'     On macOS, double-click it in Finder (the first time, macOS may ask you
#'     to confirm: right-click > Open).
#'   \item Windows: \file{LAMB classifier.bat}. Double-click it, or drop a file
#'     or folder onto it to classify that directly.
#' }
#'
#' The launcher calls the R found at \code{file.path(R.home("bin"), "Rscript")}
#' when it is written. If R is updated or moved, run this function again. The
#' window stays open at the end, so any error can be read.
#'
#' An existing launcher is never replaced silently: in an interactive session
#' you are asked first; otherwise set \code{overwrite = TRUE}.
#'
#' @param dir Folder to write the launcher into. Default: your Desktop.
#' @param platform \code{"auto"} (default; this computer), \code{"unix"}
#'   (macOS / Linux \file{.command}) or \code{"windows"} (\file{.bat}).
#' @param overwrite Replace an existing launcher without asking.
#'
#' @return Invisibly, the path of the launcher written (or \code{NULL} when you
#'   chose not to replace an existing one).
#'
#' @examples
#' # Writes to a temporary folder here; by default it goes to the Desktop.
#' path <- lamb_install_launcher(dir = tempdir())
#' readLines(path)
#' @export
lamb_install_launcher <- function(dir = file.path(path.expand("~"), "Desktop"),
                                  platform = c("auto", "unix", "windows"),
                                  overwrite = FALSE) {
    platform <- match.arg(platform)
    if (platform == "auto") platform <- if (.Platform$OS.type == "windows") "windows" else "unix"
    if (!dir.exists(dir)) stop("Folder not found: ", dir, call. = FALSE)

    windows <- platform == "windows"
    template <- system.file("launcher", if (windows) "launcher.bat" else "launcher.command",
                            package = "LAMBclassifier")
    if (!nzchar(template)) stop("Launcher template missing from the installed package.", call. = FALSE)
    target <- file.path(dir, if (windows) "LAMB classifier.bat" else "LAMB classifier.command")

    if (file.exists(target) && !overwrite) {
        if (!interactive()) {
            stop("A launcher already exists: ", target,
                 "\nIt was left unchanged. Use overwrite = TRUE to replace it.", call. = FALSE)
        }
        ans <- utils::askYesNo(paste0("A LAMB launcher already exists at\n  ", target, "\nReplace it?"),
                               default = FALSE)
        if (!isTRUE(ans)) {
            message("Left the existing launcher unchanged.")
            return(invisible(NULL))
        }
    }

    rscript <- .lamb_rscript(windows)
    txt <- gsub("{{RSCRIPT}}", rscript, readLines(template, warn = FALSE), fixed = TRUE)
    if (windows) {
        con <- file(target, "wb")
        writeLines(txt, con, sep = "\r\n")
    } else {
        con <- file(target, "wb")
        writeLines(txt, con, sep = "\n")
    }
    close(con)
    if (!windows) Sys.chmod(target, mode = "0755", use_umask = FALSE)
    message("Launcher written: ", target)
    invisible(target)
}

.lamb_rscript <- function(windows) {
    exe <- file.path(R.home("bin"), if (windows) "Rscript.exe" else "Rscript")
    if (windows) {
        if (.Platform$OS.type == "windows") exe <- normalizePath(exe, winslash = "\\", mustWork = FALSE)
        else exe <- gsub("/", "\\\\", exe)
    }
    exe
}

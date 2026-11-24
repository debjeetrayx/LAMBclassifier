# ---- Clinician-facing report pages ----
# User coordinates are inches (origin bottom-left) on one full-page region.

# ---- plain-language text ----

# Wording follows the README ("The six classes", "Limits").
.lamb_class_text <- c(
    "LHS" = "Luminal hormone-sensing expression pattern (oestrogen-driven luminal cells). Not an ER/PR test result.",
    "Apocrine" = "Molecular apocrine: an androgen-receptor (AR / FOXA1) driven luminal expression pattern. Not an ER/PR/AR test result.",
    "MA-Lo" = "Oestrogen-receptor activity well below LHS and androgen-receptor activity well below Apocrine. Not an ER/PR/AR or HER2 test result.",
    "Basal-like" = "Expression pattern of basal / luminal-progenitor breast cells. Not an ER/PR/HER2 test result.",
    "Chondroid" = "Cartilage-like expression program. Not a histological diagnosis of metaplastic carcinoma; no external validation.",
    "Squamous" = "Squamous keratin expression program. Not a histological diagnosis of metaplastic carcinoma; no external validation.",
    "NC" = "A required signature could not be scored (its genes were not measured or are missing in this sample), so no class is given."
)

.LAMB_SIG_LABEL <- c(chondroid = "Chondroid", squamous = "Squamous", basal_like = "Basal-like",
                     luminal = "LHS", apocrine = "Apocrine", malo = "MA-Lo")

.LAMB_RESEARCH_LINE <- "Research use only: expression-defined labels, not a clinical diagnosis."


# ---- confidence ----

# "clear", "close call" or "not classifiable", one per row.
.lamb_confidence <- function(results) {
    n <- nrow(results)
    if (!n) return(character())
    cls <- as.character(results$Pred_Class)
    margin <- if ("Margin" %in% names(results)) suppressWarnings(as.numeric(results$Margin)) else rep(NA_real_, n)
    flag <- if ("Flag" %in% names(results)) as.character(results$Flag) else rep("", n)
    flag[is.na(flag)] <- ""
    min_margin <- .lamb_min_margin(results)
    out <- rep("clear", n)
    out[grepl("close_call", flag) | (!is.na(margin) & margin < min_margin)] <- "close call"
    out[is.na(cls) | cls == "NC" | !cls %in% setdiff(LAMB_CLASS_ORDER, "NC") | is.na(margin)] <- "not classifiable"
    out
}

.lamb_min_margin <- function(results) {
    spec <- attr(results, "specification")
    m <- if (!is.null(spec$min_margin)) spec$min_margin else lamb_calibration()$min_margin
    as.numeric(m)
}

# Signatures below the coverage minimum, e.g. "Squamous 4/9, MA-Lo 8/15".
.lamb_low_coverage_text <- function(results) {
    cov <- attr(results, "coverage")
    if (is.null(cov) || !nrow(cov)) return("")
    low <- cov[!is.na(cov$Coverage) & cov$Coverage < .lamb_min_coverage(results), , drop = FALSE]
    if (!nrow(low)) return("")
    lab <- .LAMB_SIG_LABEL[low$Signature]
    lab[is.na(lab)] <- low$Signature[is.na(lab)]
    paste0(lab, " ", low$N_Found, "/", low$N_Requested, collapse = ", ")
}

.lamb_min_coverage <- function(results) {
    spec <- attr(results, "specification")
    m <- if (!is.null(spec$min_gene_coverage)) spec$min_gene_coverage else lamb_calibration()$min_gene_coverage
    as.numeric(m)
}


# ---- decision path model (shared by bulk and cell pages) ----

.lamb_report_num <- function(r, col) {
    if (!col %in% names(r)) return(NA_real_)
    v <- suppressWarnings(as.numeric(r[[col]]))
    if (!length(v)) NA_real_ else v[1]
}

.lamb_fmt <- function(x) ifelse(is.na(x), "not scored", formatC(x, format = "f", digits = 3))

# The deciding step comes from Pred_Class, not recomputed, so it agrees with
# the call even when stored scores are rounded.
.lamb_path <- function(r) {
    sc <- c(chondroid = .lamb_report_num(r, "Score_Chondroid"), squamous = .lamb_report_num(r, "Score_Squamous"),
            basal_like = .lamb_report_num(r, "Score_Basal_Like"), luminal = .lamb_report_num(r, "Score_LHS"),
            apocrine = .lamb_report_num(r, "Score_Apocrine"), malo = .lamb_report_num(r, "Score_MALO"))
    cls <- as.character(r[["Pred_Class"]])
    if (!length(cls) || is.na(cls)) cls <- "NC"
    decided <- switch(cls, "Chondroid" = 1L, "Squamous" = 1L, "Basal-like" = 2L,
                      "LHS" = 3L, "Apocrine" = 3L, "MA-Lo" = 3L, NA_integer_)
    required <- c("chondroid", "squamous", "basal_like", "luminal", "apocrine")
    missing <- required[is.na(sc[required])]
    nc_step <- NA_integer_
    if (is.na(decided)) {
        nc_step <- if (length(missing) || any(is.na(sc[c("chondroid", "squamous")]))) 1L
                   else if (any(is.na(sc[c("basal_like", "luminal", "apocrine")]))) 2L else 1L
    }
    stop_at <- if (!is.na(decided)) decided else nc_step

    status <- function(k) {
        if (k < stop_at) "passed" else if (k == stop_at) (if (is.na(decided)) "not_evaluable" else "decided") else "not_reached"
    }
    row <- function(label, value, key, win = FALSE, short = label) {
        list(label = label, value = unname(value), key = key, win = win, short = short)
    }
    best <- function(keys) {
        v <- sc[keys]
        if (all(is.na(v))) return(list(key = keys[1], value = NA_real_))
        k <- names(v)[which.max(v)]
        list(key = k, value = unname(v[k]))
    }
    gap <- function(a, b) if (is.na(a) || is.na(b)) NA_real_ else abs(a - b)
    # the deciding gap is the recorded Margin, so page and table agree exactly
    margin <- .lamb_report_num(r, "Margin")
    dgap <- function(a, b) if (!is.na(margin)) margin else gap(a, b)
    not_scored <- function(keys) paste(.LAMB_SIG_LABEL[keys], collapse = ", ")

    # ---- Step 1 ----
    met <- best(c("chondroid", "squamous"))
    mam <- best(c("basal_like", "luminal", "apocrine"))
    st1 <- status(1L)
    met_wins <- cls %in% c("Chondroid", "Squamous")
    g1 <- list(caption = NULL, rows = list(
        row(paste0("Metaplastic (strongest: ", .LAMB_SIG_LABEL[met$key], ")"), met$value, met$key,
            win = st1 != "not_evaluable" && met_wins, short = "Metaplastic"),
        row(paste0("Mammary (strongest: ", .LAMB_SIG_LABEL[mam$key], ")"), mam$value, mam$key,
            win = st1 == "passed" || (st1 == "decided" && !met_wins), short = "Mammary")))
    groups1 <- list(g1)
    if (st1 == "not_evaluable") {
        groups1 <- list(list(caption = "Signatures needed for step 1:", rows = lapply(required, function(k)
            row(paste0(.LAMB_SIG_LABEL[[k]], " program"), sc[k], k, short = .LAMB_SIG_LABEL[[k]]))))
    }
    if (st1 == "decided") {
        groups1[[2]] <- list(caption = "Which metaplastic program is stronger?", rows = list(
            row("Chondroid program", sc["chondroid"], "chondroid", win = cls == "Chondroid", short = "Chondroid"),
            row("Squamous program", sc["squamous"], "squamous", win = cls == "Squamous", short = "Squamous")))
    }
    out1 <- switch(st1,
        decided = paste0("Yes: the metaplastic program is stronger (ahead by ", .lamb_fmt(dgap(met$value, mam$value)),
                         "). The ", tolower(cls), " program beats the ",
                         tolower(setdiff(c("Chondroid", "Squamous"), cls)), " program by ",
                         .lamb_fmt(gap(sc["chondroid"], sc["squamous"])), ". Class: ", cls, "."),
        passed = paste0("No: mammary epithelium is stronger (ahead by ", .lamb_fmt(gap(met$value, mam$value)),
                        "). Go to step 2."),
        not_evaluable = paste0("Cannot be evaluated: ", if (length(missing)) not_scored(missing) else "a required signature",
                               " could not be scored (genes not measured or missing in this sample), so no class is given."),
        not_reached = "Not reached.")
    short1 <- switch(st1,
        decided = paste0("Metaplastic: yes (+", .lamb_fmt(dgap(met$value, mam$value)), "), ", cls),
        passed = paste0("Metaplastic: no (mammary +", .lamb_fmt(gap(met$value, mam$value)), ")"),
        not_evaluable = "Metaplastic: cannot be evaluated",
        not_reached = "Metaplastic: not reached")

    # ---- Step 2 ----
    st2 <- status(2L)
    lum_best <- best(c("luminal", "apocrine"))
    g2 <- list(caption = NULL, rows = list(
        row("Basal-like program", sc["basal_like"], "basal_like", win = st2 == "decided", short = "Basal-like"),
        row("LHS (hormone-sensing) program", sc["luminal"], "luminal",
            win = st2 == "passed" && lum_best$key == "luminal", short = "LHS"),
        row("Apocrine program", sc["apocrine"], "apocrine",
            win = st2 == "passed" && lum_best$key == "apocrine", short = "Apocrine")))
    out2 <- switch(st2,
        decided = paste0("Yes: basal-like is stronger than both (ahead by ",
                         .lamb_fmt(dgap(sc["basal_like"], lum_best$value)), "). Class: Basal-like."),
        passed = paste0("No: the ", .LAMB_SIG_LABEL[lum_best$key], " program is stronger than basal-like (by ",
                        .lamb_fmt(gap(sc["basal_like"], lum_best$value)), "). Go to step 3."),
        not_evaluable = "Cannot be evaluated: a required signature could not be scored.",
        not_reached = if (!is.na(decided)) "Not reached: the class was decided at step 1." else "Not reached.")
    short2 <- switch(st2,
        decided = paste0("Basal-like: yes (+", .lamb_fmt(dgap(sc["basal_like"], lum_best$value)), ")"),
        passed = paste0("Basal-like: no (", .LAMB_SIG_LABEL[lum_best$key], " stronger)"),
        not_evaluable = "Basal-like: cannot be evaluated",
        not_reached = "Basal-like: not reached")

    # ---- Step 3 ----
    st3 <- status(3L)
    lhs_wins <- !is.na(sc["luminal"]) && !is.na(sc["apocrine"]) && sc["luminal"] > sc["apocrine"]
    if (st3 == "decided" && cls %in% c("LHS", "Apocrine")) lhs_wins <- cls == "LHS"
    lin <- if (lhs_wins) "LHS" else "Apocrine"
    is_malo <- cls == "MA-Lo"
    g3 <- list(caption = NULL, rows = list(
        row("LHS (hormone-sensing) program", sc["luminal"], "luminal", win = st3 == "decided" && lhs_wins, short = "LHS"),
        row("Apocrine program", sc["apocrine"], "apocrine", win = st3 == "decided" && !lhs_wins, short = "Apocrine")))
    g3b <- list(caption = "MA-Lo check: is the MA-Lo program stronger than both?", rows = list(
        row("MA-Lo program", sc["malo"], "malo", win = st3 == "decided" && is_malo,
            short = "MA-Lo"),
        row(paste0("Stronger luminal program (", lin, ")"), max(sc["luminal"], sc["apocrine"]),
            if (lhs_wins) "luminal" else "apocrine", win = st3 == "decided" && !is_malo, short = lin)))
    malo_txt <- if (is.na(sc["malo"])) {
        "The MA-Lo check was skipped: the MA-Lo signature could not be scored."
    } else if (is_malo) {
        paste0("The MA-Lo program is stronger than both (ahead by ",
               .lamb_fmt(dgap(sc["malo"], max(sc["luminal"], sc["apocrine"]))), "), so it overrides.")
    } else {
        "The MA-Lo program is not stronger, so it does not override."
    }
    out3 <- switch(st3,
        decided = paste0(lin, " beats ", setdiff(c("LHS", "apocrine"), lin)[1], " (by ",
                         .lamb_fmt(if (is_malo) gap(sc["luminal"], sc["apocrine"]) else dgap(sc["luminal"], sc["apocrine"])), "). ", malo_txt, " Class: ", cls, "."),
        not_reached = if (!is.na(decided)) paste0("Not reached: the class was decided at step ", decided, ".") else "Not reached.",
        "Not reached.")
    short3 <- switch(st3,
        decided = paste0(lin, " vs apocrine: ", lin, if (is_malo) "; MA-Lo overrides" else ""),
        not_reached = "LHS vs apocrine: not reached",
        "LHS vs apocrine: not reached")
    if (st3 == "decided" && lin == "Apocrine") {
        out3 <- paste0("Apocrine beats LHS (by ", .lamb_fmt(if (is_malo) gap(sc["luminal"], sc["apocrine"]) else dgap(sc["luminal"], sc["apocrine"])), "). ",
                       malo_txt, " Class: ", cls, ".")
        short3 <- paste0("LHS vs apocrine: Apocrine", if (is_malo) "; MA-Lo overrides" else "")
    }

    runner <- switch(as.character(stop_at),
        "1" = if (met_wins) paste0("mammary epithelium (", .LAMB_SIG_LABEL[mam$key], ")") else NA_character_,
        "2" = .LAMB_SIG_LABEL[[lum_best$key]],
        "3" = if (is_malo) lin else setdiff(c("LHS", "Apocrine"), lin)[1],
        NA_character_)
    list(
        cls = cls, scores = sc, missing = missing, runner_up = runner,
        steps = list(
            list(n = 1L, status = st1, outcome = out1, short = short1, groups = groups1,
                 question = "Is a metaplastic expression program (chondroid or squamous) stronger than every mammary program?"),
            list(n = 2L, status = st2, outcome = out2, short = short2, groups = list(g2),
                 question = "Is the basal-like program stronger than both luminal programs (LHS and apocrine)?"),
            list(n = 3L, status = st3, outcome = out3, short = short3, groups = list(g3, g3b),
                 question = "Which luminal program is stronger: LHS (hormone-sensing) or apocrine? Then the MA-Lo check.")
        )
    )
}


# ---- drawing primitives (user units = inches) ----

.LAMB_PAGE <- c(w = 11.69, h = 8.27)
.LAMB_INK <- "#1a1a1a"
.LAMB_MUTED <- "#5f5f5f"
.LAMB_FAINT <- "#b5b5b5"

.lamb_new_page <- function() {
    graphics::par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", xpd = NA)
    graphics::plot.new()
    graphics::plot.window(xlim = c(0, .LAMB_PAGE[["w"]]), ylim = c(0, .LAMB_PAGE[["h"]]))
}

.lamb_lh <- function(cex) 0.2 * cex * 1.18   # line height in inches at pointsize 12

.lamb_wrap <- function(txt, width, cex = 1, font = 1) {
    out <- character()
    for (para in strsplit(txt, "\n", fixed = TRUE)[[1]]) {
        words <- strsplit(para, " +")[[1]]
        words <- words[nzchar(words)]
        line <- ""
        for (w in words) {
            cand <- if (nzchar(line)) paste(line, w) else w
            if (nzchar(line) && graphics::strwidth(cand, units = "user", cex = cex, font = font) > width) {
                out <- c(out, line)
                line <- w
            } else {
                line <- cand
            }
        }
        out <- c(out, line)
    }
    out
}

# Returns the y just below the block.
.lamb_text <- function(x, ytop, txt, width, cex = 1, font = 1, col = .LAMB_INK, adj = 0, max_lines = Inf) {
    lines <- .lamb_wrap(txt, width, cex, font)
    if (length(lines) > max_lines) {
        lines <- lines[seq_len(max_lines)]
        lines[max_lines] <- paste0(sub("[[:space:]]*\\S+$", "", lines[max_lines]), " ...")
    }
    lh <- .lamb_lh(cex)
    for (i in seq_along(lines)) {
        graphics::text(x, ytop - (i - 0.5) * lh, lines[i], adj = c(adj, 0.5), cex = cex, font = font, col = col)
    }
    ytop - length(lines) * lh
}

.lamb_fit_cex <- function(txt, width, cex, font = 1, min_cex = 0.55) {
    w <- graphics::strwidth(txt, units = "user", cex = cex, font = font)
    if (w <= width) cex else max(min_cex, cex * width / w)
}

.lamb_text_height <- function(txt, width, cex = 1, font = 1) {
    length(.lamb_wrap(txt, width, cex, font)) * .lamb_lh(cex)
}

.lamb_on_colour <- function(col) {
    rgb <- grDevices::col2rgb(col) / 255
    lum <- sum(c(0.2126, 0.7152, 0.0722) * rgb)
    if (lum > 0.55) .LAMB_INK else "white"
}

.lamb_class_col <- function(cls) {
    col <- unname(LAMB_CLASS_COLORS[cls])
    col[is.na(col)] <- "#888888"
    col
}

# Coloured chip carrying the class name (meaning never depends on colour alone).
.lamb_chip <- function(x, ytop, label, cex = 2.4, pad = 0.12) {
    col <- .lamb_class_col(label)
    w <- graphics::strwidth(label, units = "user", cex = cex, font = 2) + 2 * pad
    h <- .lamb_lh(cex) + pad
    graphics::rect(x, ytop - h, x + w, ytop, col = col, border = NA)
    graphics::text(x + pad, ytop - h / 2, label, adj = c(0, 0.5), cex = cex, font = 2, col = .lamb_on_colour(col))
    c(w = w, h = h)
}

.lamb_footer_page <- function(results, page, n_pages) {
    spec <- attr(results, "specification")
    if (is.null(spec)) spec <- lamb_specification()
    pkg <- tryCatch(as.character(utils::packageVersion("LAMBclassifier")), error = function(e) "unknown")
    meta <- paste0("LAMBclassifier ", pkg, "  |  specification ", spec$version,
                   "  |  signatures ", substr(spec$signature_sha256, 1, 12),
                   "  |  ", format(Sys.Date(), "%Y-%m-%d"),
                   "  |  page ", page, " of ", n_pages)
    W <- .LAMB_PAGE[["w"]]
    graphics::segments(.LAMB_X0, 0.62, W - .LAMB_X0, 0.62, col = .LAMB_FAINT, lwd = 0.6)
    graphics::text(.LAMB_X0, 0.45, .LAMB_RESEARCH_LINE, adj = c(0, 0.5), cex = 0.78, font = 2, col = .LAMB_MUTED)
    graphics::text(.LAMB_X0, 0.25, meta, adj = c(0, 0.5), cex = 0.7, col = .LAMB_MUTED)
}


# ---- small shared helpers ----

.LAMB_GAUGE_MAX <- 0.6

.lamb_gauge_x <- function(m, x0, x1) {
    x0 + (x1 - x0) * sqrt(pmin(pmax(m, 0), .LAMB_GAUGE_MAX) / .LAMB_GAUGE_MAX)
}

# "MA-Lo beat LHS", or the generic form when the runner-up is unknown.
.lamb_close_beat <- function(r) {
    p <- .lamb_path(r)
    if (is.na(p$runner_up)) return("The winning score beat the runner-up")
    paste0(p$cls, " beat ", p$runner_up)
}

.lamb_unscored_text <- function(r) {
    u <- if ("Unscored_Signatures" %in% names(r)) as.character(r[["Unscored_Signatures"]]) else ""
    if (is.na(u) || !nzchar(u)) return("a required signature")
    keys <- strsplit(u, ",", fixed = TRUE)[[1]]
    lab <- .LAMB_SIG_LABEL[keys]
    lab[is.na(lab)] <- keys[is.na(lab)]
    paste(lab, collapse = ", ")
}

.lamb_axis_range <- function(sc) {
    v <- sc[!is.na(sc)]
    lo <- min(c(-0.5, v - 0.05)); hi <- max(c(0.5, v + 0.05))
    c(lo, hi)
}

# ---- compact decision steps (what the pages draw) ----

# Only the steps reached, each with at most its two deciding scores.
.lamb_rpt_steps <- function(r) {
    p <- .lamb_path(r)
    sc <- p$scores; cls <- p$cls
    lab <- .LAMB_SIG_LABEL
    st <- vapply(p$steps, `[[`, "", "status")
    top <- function(keys) {
        v <- sc[keys]
        if (all(is.na(v))) keys[1] else names(v)[which.max(v)]
    }
    row <- function(text, key, win) list(text = text, key = key, value = unname(sc[key]), win = win)
    # q/a: bulk page wording; qs/as: the shorter cluster-panel wording
    step <- function(n, q, qs, a, rows = list(), decides = FALSE, as = sub(":.*", "", a)) {
        list(n = n, q = q, qs = qs, a = a, as = as, rows = rows, decides = decides)
    }
    out <- list()
    if (st[1] == "not_evaluable") {
        miss <- if (length(p$missing)) paste(lab[p$missing], collapse = ", ") else .lamb_unscored_text(r)
        out[[1]] <- step("1", "Can the gene programs be compared?", "Scored?",
                         paste0("No: ", miss, " not scored"), decides = TRUE,
                         as = paste(miss, "not scored"))
    } else if (st[1] %in% c("passed", "decided")) {
        yes <- st[1] == "decided"
        met <- if (yes && cls %in% c("Chondroid", "Squamous")) tolower(cls) else top(c("chondroid", "squamous"))
        mam <- top(c("basal_like", "luminal", "apocrine"))
        out[[1]] <- step("1", "Is a metaplastic program stronger than the mammary ones?", "Metaplastic?",
                         if (yes) paste0("Yes: ", cls) else "No",
                         list(row("Strongest metaplastic", met, yes),
                              row("Strongest mammary", mam, !yes)), yes)
    }
    if (st[2] %in% c("passed", "decided")) {
        yes <- st[2] == "decided"
        lum <- top(c("luminal", "apocrine"))
        out[[length(out) + 1]] <- step("2", "Is basal-like stronger than the luminal programs?", "Basal-like?",
                                       if (yes) "Yes: Basal-like" else "No",
                                       list(row("Basal-like", "basal_like", yes),
                                            row("Strongest luminal", lum, !yes)), yes)
    }
    if (st[3] == "decided") {
        malo <- cls == "MA-Lo"
        lhs <- if (cls %in% c("LHS", "Apocrine")) cls == "LHS" else
            (!is.na(sc["luminal"]) && !is.na(sc["apocrine"]) && sc["luminal"] > sc["apocrine"])
        lin <- if (lhs) "luminal" else "apocrine"
        lt <- if (malo) .lamb_lean_text(.lamb_malo_lean(r)[1]) else ""
        a3 <- if (malo && nzchar(lt)) paste0(toupper(substr(lt, 1, 1)), substring(lt, 2)) else lab[[lin]]
        out[[length(out) + 1]] <- step("3", "Which luminal program is stronger?", "Luminal type:", a3,
                                       list(row("LHS", "luminal", lhs && !malo), row("Apocrine", "apocrine", !lhs && !malo)), !malo)
        out[[length(out) + 1]] <- if (is.na(sc["malo"])) {
            step("", "MA-Lo check: is the MA-Lo program stronger?", "MA-Lo?", "Skipped: not scored",
                 as = "not scored")
        } else if (malo) {
            step("", "MA-Lo check: is MA-Lo stronger than both?", "MA-Lo?", "Yes: MA-Lo",
                 list(row("MA-Lo", "malo", TRUE), row(lab[[lin]], lin, FALSE)), TRUE)
        } else {
            step("", "MA-Lo check: is the MA-Lo program stronger?", "MA-Lo?", "No")
        }
    }
    skipped <- which(st == "not_reached")
    rest <- if (!length(skipped)) "" else paste0(
        if (length(skipped) > 1) paste0("Steps ", paste(skipped, collapse = " and ")) else paste("Step", skipped),
        if (cls == "NC") " not reached." else " not needed.")
    list(cls = cls, steps = out, rest = rest, any_score = any(!is.na(sc)))
}


# ---- drawing primitives for the compact pages ----

# Winner filled, runner-up hollow.
.LAMB_LEGEND <- c("wins the step", "does not win", "Further right = higher score.")
.lamb_rpt_dot <- function(x, y, win, decides, cls, cex = 1.6) {
    if (win) {
        graphics::points(x, y, pch = 21, cex = cex, lwd = 0.8, col = .LAMB_INK,
                         bg = if (decides) .lamb_class_col(cls) else "#555555")
    } else {
        graphics::points(x, y, pch = 21, cex = cex * 0.9, lwd = 1.2, col = "#7a7a7a", bg = "white")
    }
}

.lamb_rpt_label <- function(x, y, txt, cex = 0.8, font = 1, col = .LAMB_MUTED) {
    graphics::text(x, y, txt, adj = c(0, 0.5), cex = cex, font = font, col = col)
}


# ---- bulk page ----

.LAMB_X0 <- 0.6

# MA-Lo lies between LHS and Apocrine; its lean is whichever of the two scores is higher,
# or "balanced" when they differ by less than the class margin.
.lamb_malo_lean <- function(r, tie = .lamb_min_margin(r)) {
    cls <- as.character(r[["Pred_Class"]]); lhs <- suppressWarnings(as.numeric(r[["Score_LHS"]]))
    apo <- suppressWarnings(as.numeric(r[["Score_Apocrine"]]))
    if (is.null(lhs) || is.null(apo) || !length(cls)) return(rep("", length(cls)))
    out <- ifelse(abs(lhs - apo) < tie, "balanced",
                  ifelse(lhs > apo, "LHS-leaning", "apocrine-leaning"))
    out[is.na(cls) | cls != "MA-Lo" | is.na(lhs) | is.na(apo)] <- ""
    out
}

.lamb_lean_text <- function(lean) if (identical(lean, "balanced")) "balanced between LHS and Apocrine" else lean

.lamb_rpt_header <- function(r, title, idx, n, results = NULL) {
    W <- .LAMB_PAGE[["w"]]; H <- .LAMB_PAGE[["h"]]; x0 <- .LAMB_X0
    graphics::text(x0, H - 0.45, title, adj = c(0, 0.5), cex = 0.8, col = .LAMB_MUTED)
    if (n > 1) graphics::text(W - x0, H - 0.45, paste0("Sample ", idx, " of ", n), adj = c(1, 0.5),
                              cex = 0.8, col = .LAMB_MUTED)
    graphics::text(x0, H - 0.9, as.character(r[["Sample"]]), adj = c(0, 0.5), cex = 1.6, font = 2)

    cls <- as.character(r[["Pred_Class"]])
    if (!length(cls) || is.na(cls)) cls <- "NC"
    label <- if (cls == "NC") "Not classifiable" else cls
    ctop <- H - 1.3
    chip <- .lamb_chip(x0, ctop, label, cex = 2, pad = 0.14)
    ccy <- ctop - chip[["h"]] / 2

    tx <- x0 + chip[["w"]] + 0.35
    tw <- min(6.2, W - x0 - tx)
    ctext <- .lamb_class_text[[cls]]
    if (cls == "NC") {
        u <- .lamb_unscored_text(r)
        ctext <- paste0(toupper(substr(u, 1, 1)), substring(u, 2),
                        " could not be scored (genes not measured or missing in this sample), so no class is given.")
    }
    lean <- .lamb_malo_lean(r)
    if (length(lean) && nzchar(lean[1])) {
        lt <- .lamb_lean_text(lean[1])
        lt <- paste0(toupper(substr(lt, 1, 1)), substring(lt, 2),
                     switch(lean[1], "LHS-leaning" = " (LHS above Apocrine)",
                            "apocrine-leaning" = " (Apocrine above LHS)", ""))
        ctext <- paste0(lt, ". ", ctext)
    }
    meaning <- .lamb_wrap(ctext, tw, 1)
    pam <- if ("PAM50" %in% names(r)) as.character(r[["PAM50"]]) else NA_character_
    has_pam <- length(pam) && !is.na(pam) && nzchar(pam)
    amp <- .lamb_s2_her2_header(r, results)
    amp_h <- if (is.null(amp)) 0 else 0.1 + .lamb_text_height(amp$txt, tw, 0.85, amp$font)
    block <- length(meaning) * .lamb_lh(1) + amp_h + if (has_pam) 0.08 + .lamb_lh(0.8) else 0
    y <- .lamb_text(tx, ccy + block / 2, ctext, tw, cex = 1)
    if (!is.null(amp)) y <- .lamb_text(tx, y - 0.1, amp$txt, tw, cex = 0.85, font = amp$font, col = amp$col)
    if (has_pam) {
        .lamb_text(tx, y - 0.08, paste0("PAM50 (different taxonomy, for reference): ", pam), tw,
                   cex = 0.8, col = .LAMB_MUTED)
    }
    min(ctop - chip[["h"]], ccy - block / 2) - 0.45
}

.lamb_rpt_confidence <- function(r, conf, min_margin, x0, x1, ytop) {
    w <- x1 - x0
    .lamb_rpt_label(x0, ytop, "How clear-cut")
    word <- switch(conf, "clear" = "Clear", "close call" = "Close call", "Not classifiable")
    graphics::text(x0, ytop - 0.42, word, adj = c(0, 0.5), cex = 1.7, font = 2)
    margin <- .lamb_report_num(r, "Margin")
    msg <- switch(conf,
        "clear" = paste0("The winning score beat the runner-up by at least the ", min_margin, " margin (here ",
                         .lamb_fmt(margin), ")."),
        "close call" = paste0(.lamb_close_beat(r), " by less than the ", min_margin, " margin (here ",
                              .lamb_fmt(margin), "). Treat the class as uncertain."),
        { u <- .lamb_unscored_text(r)
          paste0(toupper(substr(u, 1, 1)), substring(u, 2), " could not be scored: its genes were not measured or are missing in this sample.") })
    y <- .lamb_text(x0, ytop - 0.72, msg, w, cex = 0.9) - 0.25
    if (!is.na(margin)) {
        # slim margin bar on a square-root scale so the small minimum stays visible
        bx0 <- x0; bx1 <- x0 + min(w, 3)
        xt <- .lamb_gauge_x(min_margin, bx0, bx1)
        xm <- .lamb_gauge_x(margin, bx0, bx1)
        graphics::segments(bx0, y, bx1, y, col = "#e6e6e6", lwd = 5, lend = 1)
        graphics::segments(bx0, y, xm, y, col = "#6b6b6b", lwd = 5, lend = 1)
        graphics::segments(xt, y - 0.07, xt, y + 0.07, col = .LAMB_INK, lwd = 1)
        graphics::text(bx0, y - 0.17, "close call", adj = c(0, 0.5), cex = 0.62, col = .LAMB_MUTED)
        graphics::text(bx1, y - 0.17, "clear", adj = c(1, 0.5), cex = 0.62, col = .LAMB_MUTED)
        y <- y - 0.3
    }
    y
}

# Drawn only when there is something the reader needs: returns ytop otherwise.
.lamb_rpt_coverage <- function(r, coverage, conf, x0, x1, ytop) {
    low <- .lamb_low_coverage_text(r)
    chon <- .lamb_s2_chondroid_text(r)
    if (!nzchar(low) && !nzchar(chon)) return(ytop + 0.35)
    .lamb_rpt_label(x0, ytop, "Genes measured")
    y <- ytop - 0.22
    if (nzchar(low)) {
        txt <- paste0(if (grepl("(^| )0/", low)) "Genes measured, below " else "Few genes measured, below ",
                      round(100 * .lamb_min_coverage(r)), "%: ", low,
                      ".", if (conf != "not classifiable") " Read the class with care.")
        y <- .lamb_text(x0, y, txt, x1 - x0, cex = 0.9, font = 2) - 0.06
    }
    if (nzchar(chon)) y <- .lamb_text(x0, y, chon, x1 - x0, cex = 0.9)
    y
}


# ---- specification 2.0.0 read-outs on the pages ----

.lamb_s2_chondroid_text <- function(results) {
    cov <- attr(results, "coverage")
    if (is.null(cov) || !nrow(cov) || !"chondroid" %in% cov$Signature) return("")
    k <- cov[cov$Signature == "chondroid", , drop = FALSE][1, ]
    if (is.na(k$N_Found) || k$N_Found >= k$N_Requested) return("")
    if (k$N_Found == 0) return("Chondroid: not measured in this data, so it cannot be called.")
    uc <- attr(results, "universe_coverage")
    if (is.null(uc) || uc * length(lamb_universe()) >= 1000) return("")   # warn on targeted panels only
    "Chondroid: partly measured on this panel."
}

.lamb_s2_q <- function(results) {
    spec <- attr(results, "specification")
    q <- spec$her2_amplicon$q %||% lamb_calibration()$her2_amplicon$q
    as.numeric(q)
}

.lamb_s2_her2_lines <- function(r, results) {
    if (!"HER2_Amplicon" %in% names(r)) return(NULL)
    st <- as.character(r[["HER2_Amplicon"]])
    if (!length(st) || is.na(st)) return(NULL)
    q <- .lamb_s2_q(results)
    hi <- .lamb_report_num(r, "HER2_Amplicon_N_High")
    ms <- .lamb_report_num(r, "HER2_Amplicon_N_Measured")
    amp <- lamb_calibration()$her2_amplicon
    counts <- if (!is.na(hi) && !is.na(ms)) paste0(hi, " of ", ms, " measured amplicon genes rank in this sample's top ",
                                                   q, "% (", amp$min_high, " needed).") else ""
    switch(st,
        "high" = list(word = "High", txt = counts),
        "not high" = list(word = "Not high", txt = counts),
        "not assessable" = list(word = "Not assessed", txt = ""),
        list(word = st, txt = ""))
}

.lamb_s2_her2_header <- function(r, results) {
    h <- .lamb_s2_her2_lines(r, results)
    if (is.null(h)) return(NULL)
    high <- h$word == "High"
    if (h$word == "Not assessed") {
        return(list(txt = "HER2 amplicon: not assessed on this platform.", font = 1, col = .LAMB_MUTED))
    }
    txt <- paste0("HER2 amplicon: ", tolower(h$word), ". ", h$txt,
                  if (high) " Expression only, not HER2 IHC/ISH; it does not change the class." else "")
    list(txt = txt, font = if (high) 2 else 1, col = if (high) .LAMB_INK else .LAMB_MUTED)
}


.lamb_rpt_path <- function(m, x0, x1, ytop) {
    .lamb_rpt_label(x0, ytop, if (m$cls == "NC") "Why no class" else "Why this class")
    y <- ytop - 0.32
    if (!m$any_score && m$cls != "NC") {
        .lamb_text(x0, y, "The scores behind this call were not kept with these results.", x1 - x0,
                   cex = 0.9, col = .LAMB_MUTED)
        return(invisible(y))
    }
    vals <- unlist(lapply(m$steps, function(s) vapply(s$rows, `[[`, 0, "value")))
    rng <- .lamb_axis_range(vals)
    lab_w <- 2.05
    ax0 <- x0 + 0.35 + lab_w + 0.2; ax1 <- x1 - 0.15
    sx <- function(v) ax0 + (ax1 - ax0) * (v - rng[1]) / diff(rng)
    rh <- 0.27
    for (s in m$steps) {
        cy <- y - 0.13
        if (nzchar(s$n)) {
            graphics::symbols(x0 + 0.12, cy, circles = 0.11, inches = FALSE, add = TRUE,
                              bg = "#efefef", fg = NA)
            graphics::text(x0 + 0.12, cy, s$n, cex = 0.78, font = 2, col = .LAMB_INK)
        }
        qx <- x0 + 0.35
        graphics::text(qx, cy, s$q, adj = c(0, 0.5), cex = 0.95, col = .LAMB_INK)
        ax <- qx + graphics::strwidth(s$q, units = "user", cex = 0.95) + 0.15
        graphics::text(ax, cy, s$a, adj = c(0, 0.5), cex = 0.95, font = 2, col = .LAMB_INK)
        y <- y - 0.3
        for (rw in s$rows) {
            ry <- y - rh / 2
            graphics::text(qx + lab_w, ry, rw$text, adj = c(1, 0.5), cex = 0.82,
                           font = if (rw$win) 2 else 1, col = if (rw$win) .LAMB_INK else .LAMB_MUTED)
            graphics::segments(ax0, ry, ax1, ry, col = "#ededed", lwd = 1)
            if (is.na(rw$value)) {
                graphics::text(ax0, ry, "not scored", adj = c(0, 0.5), cex = 0.75, font = 3, col = .LAMB_MUTED)
            } else {
                .lamb_rpt_dot(sx(rw$value), ry, rw$win, s$decides, m$cls)
            }
            y <- y - rh
        }
        y <- y - if (length(s$rows)) 0.22 else 0.1
    }
    if (nzchar(m$rest)) {
        graphics::text(x0 + 0.35, y - 0.1, m$rest, adj = c(0, 0.5), cex = 0.85, col = .LAMB_FAINT)
        y <- y - 0.35
    }
    if (any(vapply(m$steps, function(s) length(s$rows) > 0, TRUE))) {
        ly <- y - 0.15
        lx <- x0 + 0.35
        for (k in seq_along(.LAMB_LEGEND)) {
            if (k < 3) .lamb_rpt_dot(lx + 0.05, ly, k == 1, FALSE, m$cls, cex = 0.9)
            tx <- lx + if (k < 3) 0.15 else 0
            graphics::text(tx, ly, .LAMB_LEGEND[[k]], adj = c(0, 0.5), cex = 0.72, col = .LAMB_MUTED)
            lx <- tx + graphics::strwidth(.LAMB_LEGEND[[k]], units = "user", cex = 0.72) + 0.35
        }
        y <- ly - 0.15
    }
    invisible(y)
}


# ---- cover page ----

.lamb_rpt_cover <- function(results, title, conf) {
    W <- .LAMB_PAGE[["w"]]; H <- .LAMB_PAGE[["h"]]; x0 <- .LAMB_X0
    graphics::text(x0, H - 0.75, title, adj = c(0, 0.5), cex = 1.7, font = 2)
    graphics::text(x0, H - 1.1, paste0(nrow(results), " samples. One page per sample follows, in input order."),
                   adj = c(0, 0.5), cex = 0.95, col = .LAMB_MUTED)

    cls <- as.character(results$Pred_Class)
    cls[is.na(cls)] <- "NC"
    lv <- LAMB_CLASS_ORDER
    counts <- table(factor(cls, levels = lv))
    ytop <- H - 1.75
    .lamb_rpt_label(x0, ytop, "Classes", cex = 0.85)
    bx0 <- x0 + 1.65; bx1 <- 6.2
    rh <- 0.4
    mx <- max(counts, 1)
    for (i in seq_along(lv)) {
        cy <- ytop - 0.35 - (i - 0.5) * rh
        k <- counts[[i]]
        graphics::text(bx0 - 0.15, cy, if (lv[i] == "NC") "Not classifiable" else lv[i], adj = c(1, 0.5),
                       cex = 0.9, col = if (k) .LAMB_INK else .LAMB_FAINT)
        xe <- bx0 + (bx1 - bx0) * k / mx
        if (k > 0) graphics::rect(bx0, cy - rh * 0.3, xe, cy + rh * 0.3, col = .lamb_class_col(lv[i]), border = NA)
        graphics::text(xe + 0.1, cy, if (k) paste0(k, "  (", round(100 * k / sum(counts)), "%)") else "0",
                       adj = c(0, 0.5), cex = 0.85, col = if (k) .LAMB_INK else .LAMB_FAINT)
    }

    cx <- 7.3
    .lamb_rpt_label(cx, ytop, "How clear-cut", cex = 0.85)
    cf <- table(factor(conf, levels = c("clear", "close call", "not classifiable")))
    nm <- c("Clear", "Close call", "Not classifiable")
    for (i in seq_along(cf)) {
        cy <- ytop - 0.35 - (i - 0.5) * rh
        graphics::text(cx, cy, nm[i], adj = c(0, 0.5), cex = 0.95, col = if (cf[[i]]) .LAMB_INK else .LAMB_FAINT)
        graphics::text(cx + 1.75, cy, cf[[i]], adj = c(1, 0.5), cex = 0.95, font = 2,
                       col = if (cf[[i]]) .LAMB_INK else .LAMB_FAINT)
    }
    y <- ytop - 0.35 - 3 * rh - 0.3
    low <- .lamb_low_coverage_text(results)
    if (nzchar(low)) {
        y <- .lamb_text(cx, y, paste0("Few genes measured for ", low, ": read every call with care."),
                        W - x0 - cx, cex = 0.85, font = 2) - 0.1
    }
    chon <- .lamb_s2_chondroid_text(results)
    if (nzchar(chon)) y <- .lamb_text(cx, y, chon, W - x0 - cx, cex = 0.85) - 0.1
    if ("HER2_Amplicon" %in% names(results)) {
        amp <- as.character(results$HER2_Amplicon)
        .lamb_text(cx, y - 0.1, paste0("HER2 amplicon high: ", sum(amp == "high", na.rm = TRUE), " of ", length(amp),
                                       " samples", if (any(amp == "not assessable", na.rm = TRUE))
                                           paste0(" (", sum(amp == "not assessable", na.rm = TRUE), " not assessed)") else "",
                                       "."), W - x0 - cx, cex = 0.85)
    }

    hy <- ytop - 0.35 - length(lv) * rh - 0.55
    .lamb_rpt_label(x0, hy, "How to read each page", cex = 0.85)
    help <- c("Each page gives the class (the coloured box), how clear-cut the call is, and the steps that led to it.",
              paste0("Clear: the winning score beat the runner-up by at least the ", .lamb_min_margin(results),
                     " margin; otherwise a close call. Not classifiable: a needed signature could not be measured."),
              "HER2 amplicon is reported beside every class and never changes it.",
              "Classes are expression patterns, not ER/PR/HER2 test results or a histological diagnosis.")
    y <- hy - 0.25
    for (h in help) y <- .lamb_text(x0, y, h, W - 2 * x0, cex = 0.85, col = .LAMB_INK) - 0.06
}


# ---- exported report ----

#' Write a one-page-per-sample PDF report for clinicians
#'
#' @description
#' Turns the output of [LAMBclassifier()] into a PDF that can be read without
#' training in the classifier: one A4 landscape page per sample, preceded by a
#' cohort cover page when there is more than one sample.
#'
#' @details
#' Each sample page shows the class and its meaning, how sure the call is
#' (margin against the 0.02 minimum), low gene coverage when it matters, the
#' HER2 amplicon read-out, the decision steps that were reached, and a strip
#' of nine marker genes (ESR1, PGR, AR, ERBB2, FOXA1, GATA3, KRT5, TP63,
#' MKI67) as within-sample percentile ranks. Every page carries a
#' research-use footer with the package and specification versions.
#'
#' Keep the attributes of the classifier output: subsetting a tibble with
#' \code{[} or joining it with dplyr drops the \code{"coverage"},
#' \code{"gene_percentiles"} and \code{"specification"} attributes. Without
#' them the page leaves out coverage and the marker strip and uses the bundled
#' specification.
#'
#' @param results Output of [LAMBclassifier()] for bulk samples.
#' @param file Path of the PDF to write.
#' @param title Optional report title. Default: "LAMB classification report".
#'
#' @return Invisibly, \code{file}.
#'
#' @examples
#' \dontrun{
#' res <- LAMBclassifier(expr_log2, run_pam50 = TRUE)
#' lamb_report(res, "LAMB_report.pdf")
#' }
#' @export
lamb_report <- function(results, file, title = NULL) {
    if (!is.data.frame(results) || !all(c("Sample", "Pred_Class") %in% names(results))) {
        stop("`results` must be the output of LAMBclassifier() (columns Sample and Pred_Class).", call. = FALSE)
    }
    if (!nrow(results)) stop("`results` has no samples.", call. = FALSE)
    title <- if (is.null(title)) "LAMB classification report" else as.character(title)
    conf <- .lamb_confidence(results)
    min_margin <- .lamb_min_margin(results)
    coverage <- attr(results, "coverage")
    rows <- as.data.frame(results, stringsAsFactors = FALSE)
    n <- nrow(rows)
    n_pages <- n + (n > 1)
    W <- .LAMB_PAGE[["w"]]; x0 <- .LAMB_X0

    grDevices::pdf(file, width = .LAMB_PAGE[["w"]], height = .LAMB_PAGE[["h"]], bg = "white",
                   title = title, useDingbats = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    page <- 0L
    if (n > 1) {
        .lamb_new_page()
        page <- page + 1L
        .lamb_rpt_cover(results, title, conf)
        .lamb_footer_page(results, page, n_pages)
    }
    for (i in seq_len(n)) {
        r <- as.list(rows[i, , drop = FALSE])
        .lamb_new_page()
        page <- page + 1L
        top <- .lamb_rpt_header(r, title, i, n, results)
        lx1 <- 4.0
        # NC: the chip and header already say it, so no confidence block
        y <- if (identical(as.character(r[["Pred_Class"]]), "NC")) top + 0.35 else
            .lamb_rpt_confidence(r, conf[i], min_margin, x0, lx1, top)
        y <- .lamb_rpt_coverage(results, coverage, conf[i], x0, lx1, y - 0.35)
        .lamb_gp_strip(results, as.character(r[["Sample"]]), x0, lx1, y - 0.4)
        .lamb_rpt_path(.lamb_rpt_steps(r), 4.75, W - x0, top)
        .lamb_footer_page(results, page, n_pages)
    }
    invisible(file)
}


# ---- single-cell / spatial page ----

.lamb_rpt_mini <- function(m, cl, conf, x0, y0, w, h, extra = "") {
    x1 <- x0 + w; ytop <- y0 + h
    graphics::rect(x0, y0, x1, ytop, border = "#dddddd", lwd = 0.7)
    pad <- 0.12
    cx <- x0 + pad; cw <- w - 2 * pad
    cells <- paste0(format(cl$Cells, big.mark = ","), " cells")
    cells_w <- graphics::strwidth(cells, units = "user", cex = 0.68)
    nm <- .lamb_wrap(as.character(cl$Cluster), cw - cells_w - 0.1, 0.85, 2)[1]
    graphics::text(cx, ytop - 0.2, nm, adj = c(0, 0.5), cex = 0.85, font = 2)
    graphics::text(x1 - pad, ytop - 0.2, cells, adj = c(1, 0.5), cex = 0.68, col = .LAMB_MUTED)
    cls <- m$cls
    label <- if (cls == "NC") "Not classifiable" else cls
    chip <- .lamb_chip(cx, ytop - 0.36, label, cex = 0.9, pad = 0.06)
    word <- switch(conf, "clear" = "Clear", "close call" = "Close call", "")
    if (nzchar(word)) graphics::text(cx + chip[["w"]] + 0.1, ytop - 0.36 - chip[["h"]] / 2, word,
                                     adj = c(0, 0.5), cex = 0.72, font = 2)
    y <- ytop - 0.36 - chip[["h"]] - 0.1

    vals <- unlist(lapply(m$steps, function(s) vapply(s$rows, `[[`, 0, "value")))
    rng <- .lamb_axis_range(vals)
    ax0 <- x0 + w * 0.5; ax1 <- x1 - pad - 0.04
    sx <- function(v) ax0 + (ax1 - ax0) * (v - rng[1]) / diff(rng)
    n_lines <- length(m$steps) + nzchar(m$rest) + nzchar(extra)
    lh <- min(0.24, (y - y0 - 0.06) / max(n_lines, 1))
    tcex <- 0.6
    for (s in m$steps) {
        cy <- y - lh / 2
        txt <- paste(s$qs, s$as)
        avail <- (if (length(s$rows)) ax0 - 0.08 else x1 - pad) - cx
        fc <- .lamb_fit_cex(txt, avail, tcex, font = if (s$decides) 2 else 1, min_cex = 0.45)
        graphics::text(cx, cy, txt, adj = c(0, 0.5), cex = fc, font = if (s$decides) 2 else 1, col = .LAMB_INK)
        if (length(s$rows)) {
            graphics::segments(ax0, cy, ax1, cy, col = "#ededed", lwd = 0.8)
            for (rw in s$rows[order(vapply(s$rows, `[[`, TRUE, "win"))]) {   # winner drawn last, on top
                if (!is.na(rw$value)) .lamb_rpt_dot(sx(rw$value), cy, rw$win, s$decides, cls, cex = 1.05)
            }
        }
        y <- y - lh
    }
    if (nzchar(m$rest)) {
        graphics::text(cx, y - lh / 2, m$rest, adj = c(0, 0.5), cex = tcex, col = .LAMB_FAINT)
        y <- y - lh
    }
    if (nzchar(extra)) {
        graphics::text(cx, y - lh / 2, extra, adj = c(0, 0.5), cex = .lamb_fit_cex(extra, w - 2 * pad, tcex, min_cex = 0.45),
                       col = .LAMB_INK)
    }
}

.lamb_s2_mini_text <- function(d) {
    amp <- if ("HER2_Amplicon" %in% names(d)) as.character(d[["HER2_Amplicon"]]) else NA_character_
    her2 <- if (!length(amp) || is.na(amp) || amp == "not high") "" else
        if (amp == "not assessable") "HER2 amplicon: not assessed on this platform." else paste("HER2 amplicon", amp)
    her2   # the MA-Lo lean is already the step-3 answer in the panel
}

# Single-cell / spatial: one panel per epithelial cluster, up to 12 per page.
# `calls`: Cluster, Cells, Cell_type, Epithelial_share, LAMB_class, Note.
# `details`: LAMBclassifier() output keyed by Sample = cluster name.
lamb_report_cells <- function(calls, details, file, title = NULL) {
    need <- c("Cluster", "Cells", "Cell_type", "LAMB_class")
    if (!all(need %in% names(calls))) stop("`calls` needs columns ", paste(need, collapse = ", "), call. = FALSE)
    title <- if (is.null(title)) "LAMB classes per cell cluster" else as.character(title)
    calls <- as.data.frame(calls, stringsAsFactors = FALSE)
    det <- as.data.frame(details, stringsAsFactors = FALSE)
    shown <- calls[calls$Cell_type == "epithelial" & as.character(calls$Cluster) %in% det$Sample, , drop = FALSE]
    others <- calls[!calls$Cluster %in% shown$Cluster, , drop = FALSE]
    conf_all <- .lamb_confidence(details)
    names(conf_all) <- det$Sample

    per <- 12L
    n_pages <- max(1L, ceiling(nrow(shown) / per))
    W <- .LAMB_PAGE[["w"]]; H <- .LAMB_PAGE[["h"]]; x0 <- .LAMB_X0
    grDevices::pdf(file, width = W, height = H, bg = "white", title = title, useDingbats = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    for (p in seq_len(n_pages)) {
        .lamb_new_page()
        graphics::text(x0, H - 0.55, title, adj = c(0, 0.5), cex = 1.4, font = 2)
        sub <- paste0(nrow(shown), if (nrow(shown) == 1) " epithelial cluster, classified from its pooled cells. "
                      else " epithelial clusters, each classified from its pooled cells. ",
                      "Filled dot = wins the step; further right = higher score.")
        yb <- .lamb_text(x0, H - 0.75, sub, W - 2 * x0, cex = 0.78, col = .LAMB_MUTED)
        low <- .lamb_low_coverage_text(details)
        if (nzchar(low)) {
            yb <- .lamb_text(x0, yb - 0.05, paste0("Few genes measured in this data for ", low,
                                                   ": read every cluster with care."),
                             W - 2 * x0, cex = 0.78, font = 2)
        }
        chon <- .lamb_s2_chondroid_text(details)
        if (nzchar(chon)) yb <- .lamb_text(x0, yb - 0.05, chon, W - 2 * x0, cex = 0.78)
        idx <- shown[seq_len(nrow(shown)) > (p - 1) * per & seq_len(nrow(shown)) <= p * per, , drop = FALSE]
        gtop <- yb - 0.2; gbot <- 1.1
        cols <- 4; rws <- 3; gap <- 0.18
        pw <- (W - 2 * x0 - (cols - 1) * gap) / cols
        ph <- (gtop - gbot - (rws - 1) * gap) / rws
        for (j in seq_len(nrow(idx))) {
            cl <- idx[j, ]
            d <- as.list(det[match(as.character(cl$Cluster), det$Sample), , drop = FALSE])
            ci <- (j - 1) %% cols; ri <- (j - 1) %/% cols
            .lamb_rpt_mini(.lamb_rpt_steps(d), cl, unname(conf_all[as.character(cl$Cluster)]),
                           x0 + ci * (pw + gap), gtop - (ri + 1) * ph - ri * gap, pw, ph,
                           extra = .lamb_s2_mini_text(d))
        }
        if (!nrow(idx)) .lamb_text(x0, gtop, "No epithelial cluster was classified.", W - 2 * x0, cex = 1)
        if (p == n_pages && nrow(others)) {
            ne <- sum(others$Cell_type != "epithelial")
            tf <- sum(others$Cell_type == "epithelial")
            parts <- c(if (ne) paste0(ne, " non-epithelial cluster(s) (the classes do not apply)"),
                       if (tf) paste0(tf, " epithelial cluster(s) not classified (too few cells)"))
            .lamb_text(x0, 0.98, paste0("Not shown: ", paste(parts, collapse = " and "), ". ",
                                        "All clusters are listed in LAMB_calls.csv."),
                       W - 2 * x0, cex = 0.7, col = .LAMB_MUTED)
        }
        .lamb_footer_page(details, p, n_pages)
    }
    invisible(file)
}

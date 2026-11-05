# src/bdhs2022/tables_figures.R
# Formatting of saved results into manuscript tables, and the figures.
# Nothing is estimated here.

fmt <- function(x, d = 2) {
  x <- round(as.numeric(x), d)
  x[!is.na(x) & x == 0] <- 0                 # avoid "-0.00"
  ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
}
fmt_ci <- function(est, lo, hi, d = 2) {
  sep <- ifelse(!is.na(lo) & lo < 0, " to ", "–")
  ifelse(is.na(lo), fmt(est, d), paste0(fmt(est, d), " (", fmt(lo, d), sep, fmt(hi, d), ")"))
}
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", fmt(p, 3)))

# Readable labels for variables and model terms
var_labels <- c(age_group = "Age (years)", age_group_alt = "Age (years)",
                education = "Education", working = "Currently working",
                residence = "Place of residence", wealth = "Household wealth tertile",
                birth_order = "Birth order", trimester = "Trimester of pregnancy",
                hypertension = "Hypertension", gest_weeks = "Duration of pregnancy (weeks)",
                bmi = "Pre-pregnancy BMI (kg/m²)", age = "Age (years)",
                age_c = "Age (per year)", bmi_measured = "Measured BMI (kg/m²)")
reference_level <- c(age_group = "<25", age_group_alt = "<30", residence = "Rural",
                     wealth = "Lowest", trimester = "1st", education = "Primary or below",
                     hypertension = "No")

split_term <- function(term) {
  v <- vapply(term, function(t) {
    hit <- names(var_labels)[startsWith(t, names(var_labels))]
    if (length(hit)) hit[which.max(nchar(hit))] else t
  }, "")
  data.frame(variable = unname(v), level = substring(term, nchar(v) + 1))
}
term_label <- function(term) {
  s <- split_term(term)
  lab <- ifelse(s$variable %in% names(var_labels), var_labels[s$variable], s$variable)
  ifelse(s$level == "", lab, paste0(lab, ": ", s$level))
}

# ---- Tables ---------------------------------------------------------------------
# Results come in long format with a `scale` column (released, unconverted);
# the tables place the two scales side by side.

scale_label <- c(released = "Released scale", unconverted = "Unconverted scale")
scale_of_outcome <- function(o) ifelse(grepl("_alt$", o), "unconverted", "released")

# Table 1 (manuscript): characteristics by household wealth tertile, and GDM
# on the released and unconverted scales, in one table.
# Variable rows carry the p values; category rows are indented below them.
# tw: table1_by_wealth; t1_cat / t1_cont: GDM by characteristic (both scales).
table1_wealth <- function(tw, t1_cat, t1_cont) {
  groups <- c("Total", "Lowest", "Middle", "Highest")
  t1_cat$scale <- scale_of_outcome(t1_cat$outcome)
  t1_cont$scale <- scale_of_outcome(t1_cont$outcome)
  npct <- function(n, pct) paste0(n, " (", fmt(pct, 1), ")")
  wcell <- function(r, g) {
    x <- r[r$group == g, ]
    if (!nrow(x)) return("–")
    if (is.na(x$n)) paste0(fmt(x$value, 1), " (", fmt(x$sd, 1), ")") else npct(x$n, x$value)
  }
  gcat <- function(v, l, s) {
    x <- t1_cat[t1_cat$variable == v & t1_cat$level == l & t1_cat$scale == s, ]
    npct(x$n_gdm, x$pct_gdm)
  }
  gp <- function(v, s) fmt_p(t1_cat$p_value[t1_cat$variable == v & t1_cat$scale == s][1])
  first <- tw[tw$variable == tw$variable[1], ]
  ns <- vapply(groups, function(g) sum(first$n[first$group == g]), 0)
  gt <- function(s) { x <- tw[tw$variable == paste0("gdm_", s) & tw$group == "Total", ]; npct(x$n, x$value) }

  out <- list(c("All women, n", ns[["Total"]], ns[["Lowest"]], ns[["Middle"]], ns[["Highest"]], "",
                gt("released"), "", gt("unconverted"), ""))
  for (v in unique(t1_cat$variable)) {
    r <- tw[tw$variable == v, ]
    base <- t1_cat[t1_cat$variable == v & t1_cat$scale == "released", ]
    wp <- if (nrow(r)) fmt_p(r$p_value[1]) else ""
    out[[length(out) + 1]] <- c(paste0(var_labels[[v]], ", n (%)"), rep("", 4), wp,
                                "", gp(v, "released"), "", gp(v, "unconverted"))
    for (i in seq_len(nrow(base))) {
      l <- base$level[i]
      w <- if (nrow(r)) vapply(groups[-1], function(g) wcell(r[r$level == l, ], g), "") else rep("–", 3)
      out[[length(out) + 1]] <- c(paste0("  ", l), npct(base$n[i], base$pct_of_sample[i]), w, "",
                                  gcat(v, l, "released"), "", gcat(v, l, "unconverted"), "")
    }
  }
  for (v in unique(t1_cont$variable)) {
    r <- tw[tw$variable == v, ]
    gc <- function(s) { x <- t1_cont[t1_cont$variable == v & t1_cont$scale == s & t1_cont$group == "GDM", ]
      c(paste0(fmt(x$mean, 1), " (", fmt(x$sd, 1), ")"), fmt_p(x$p_value)) }
    out[[length(out) + 1]] <- c(paste0(var_labels[[v]], ", mean (SD)"), vapply(groups, function(g) wcell(r, g), ""),
                                fmt_p(r$p_value[1]), gc("released"), gc("unconverted"))
  }
  out <- as.data.frame(do.call(rbind, out))
  names(out) <- c("Characteristic", groups, "p (wealth)", "GDM, released scale", "p (released)",
                  "GDM, unconverted scale", "p (unconverted)")
  rownames(out) <- NULL
  out
}

# Effects table with reference rows for factors; one column (and optionally a
# p value column) per scale x model x measure.
effects_table <- function(eff, with_p = TRUE) {
  if (is.null(eff$scale)) eff$scale <- "released"
  eff$column <- paste(scale_label[eff$scale], eff$model, eff$measure, sep = ": ")
  s <- split_term(unique(eff$term))
  rows <- do.call(rbind, lapply(unique(s$variable), function(v) {
    lv <- s$level[s$variable == v]
    if (v %in% names(reference_level)) {
      data.frame(variable = v, level = c(reference_level[[v]], lv), term = c(NA, paste0(v, lv)))
    } else data.frame(variable = v, level = "per unit", term = v)
  }))
  out <- data.frame(Characteristic = var_labels[rows$variable], Category = rows$level)
  for (cc in unique(eff$column)) {
    e <- eff[eff$column == cc, ]
    out[[paste0(cc, " (95% CI)")]] <- vapply(rows$term, function(t) {
      if (is.na(t)) return("1.00")
      r <- e[e$term == t, ]
      if (nrow(r)) fmt_ci(r$estimate, r$ci_lower, r$ci_upper) else ""
    }, "")
    if (with_p) out[[paste0(cc, " p")]] <- vapply(rows$term, function(t) {
      r <- e[e$term %in% t, ]
      if (nrow(r)) fmt_p(r$p_value) else ""
    }, "")
  }
  # reference rows are blank for models that do not contain the variable
  for (cc in grep("(95% CI)", names(out), value = TRUE, fixed = TRUE)) {
    has <- tapply(out[[cc]] != "" & out[[cc]] != "1.00", rows$variable, any)
    out[[cc]][is.na(rows$term) & !has[rows$variable]] <- ""
  }
  out
}

# Inequality indices and decomposition. Compact (main text): contribution and
# percent with 95% CI per scale. Full (supplement): adds elasticity and the
# concentration index of each covariate.
inequality_table <- function(ci_boot, decomp, full = FALSE) {
  idx_names <- c(erreygers = "Erreygers corrected concentration index",
                 wagstaff = "Wagstaff normalised index", relative = "Relative concentration index")
  block <- function(s) {
    ci <- ci_boot[ci_boot$scale == s, ]
    ci <- ci[match(names(idx_names), ci$index), ]
    dec <- decomp[decomp$scale == s, ]
    top <- data.frame(Component = unname(idx_names), Elasticity = "", `Concentration index of covariate` = "",
                      `Contribution (95% CI)` = fmt_ci(ci$estimate, ci$ci_lower, ci$ci_upper, 3),
                      `Percent (95% CI)` = c("100.0", "", ""), check.names = FALSE)
    body <- data.frame(
      Component = ifelse(dec$term == "residual", "Residual (unexplained)", term_label(dec$term)),
      Elasticity = fmt(dec$elasticity, 3),
      `Concentration index of covariate` = fmt(dec$concentration_index, 3),
      `Contribution (95% CI)` = fmt_ci(dec$contribution, dec$contribution_ci_lower, dec$contribution_ci_upper, 3),
      `Percent (95% CI)` = fmt_ci(dec$percent, dec$percent_ci_lower, dec$percent_ci_upper, 1),
      check.names = FALSE)
    b <- rbind(top, body)
    if (!full) b <- b[, c("Component", "Contribution (95% CI)", "Percent (95% CI)")]
    names(b)[-1] <- paste0(scale_label[[s]], ": ", names(b)[-1])
    b
  }
  blocks <- lapply(names(scale_label), block)
  out <- blocks[[1]]
  for (b in blocks[-1]) out <- cbind(out, b[match(out$Component, b$Component), -1, drop = FALSE])
  rownames(out) <- NULL
  out
}

# Sensitivity analyses: one panel per scale
sensitivity_table <- function(prev, eff, ci) {
  show <- c("wealthHighest", "residenceUrban", "trimester3rd")
  out <- data.frame(Scale = scale_label[prev$scale], Analysis = paste(prev$id, prev$scenario), n = prev$n,
                    `Prevalence, % (95% CI)` = fmt_ci(100 * prev$prevalence, 100 * prev$ci_lower,
                                                      100 * prev$ci_upper, 1), check.names = FALSE)
  key <- paste(prev$scale, prev$id)
  out$`Erreygers index (95% CI)` <- vapply(key, function(k) {
    r <- ci[paste(ci$scale, ci$id) == k, ]
    if (nrow(r)) fmt_ci(r$estimate, r$ci_lower, r$ci_upper, 3) else ""
  }, "")
  for (t in show) {
    out[[paste0("aPR ", term_label(t), " (95% CI)")]] <- vapply(key, function(k) {
      r <- eff[paste(eff$scale, eff$id) == k & eff$term == t, ]
      if (nrow(r)) fmt_ci(r$estimate, r$ci_lower, r$ci_upper) else ""
    }, "")
  }
  out <- out[order(match(prev$scale, names(scale_label)), as.numeric(sub("S", "", prev$id))), ]
  rownames(out) <- NULL
  out
}


# ---- Figures ----------------------------------------------------------------------
# Journal style: Times (serif) type, black and grey only, white background,
# L-shaped axes, panel letters in bold. Released scale: black, solid line,
# filled circles; unconverted scale: grey, dashed line, open circles.

ink <- "black"
grey <- "grey50"
scale_cols <- c(released = "black", unconverted = "grey45")
scale_pch <- c(released = 19, unconverted = 21)
scale_lty <- c(released = 1, unconverted = 2)

# Draws a figure twice at full-page width (190 mm): vector PDF with embedded
# Times New Roman, and PNG at 1000 dpi. width and height set the aspect ratio.
save_figure <- function(file, width, height, draw, pointsize = 9) {
  w <- 190 / 25.4
  h <- w * height / width
  style <- function() par(family = "Times New Roman", cex.axis = 0.85, cex.lab = 0.95, mgp = c(2.3, 0.6, 0), tcl = -0.3,
                          las = 1, bty = "l", col.axis = "black", col.lab = "black", fg = "black")
  png(file, width = w, height = h, units = "in", res = 1000, type = "cairo", family = "Times New Roman",
      pointsize = pointsize, bg = "white")
  style(); draw(); dev.off()
  cairo_pdf(sub(".png", ".pdf", file, fixed = TRUE), width = w, height = h, family = "Times New Roman", pointsize = pointsize, bg = "white")
  style(); draw(); invisible(dev.off())
}
panel_letter <- function(l, title = "") mtext(trimws(paste0("(", l, ") ", title)), side = 3, adj = 0, font = 2, line = 0.6, cex = 0.95)
comma <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

# Box-and-arrow diagram in the style of a journal flow chart: a main column of
# boxes on the left joined by vertical arrows, and exclusion boxes on the right
# reached by horizontal arrows. Text is left-aligned inside each box.
flow_chart <- function(main, side, file, wrap_main = 62, wrap_side = 50) {
  wrapl <- function(s, w) unlist(lapply(strsplit(s, "\n", fixed = TRUE)[[1]], function(p) if (p == "") "" else strwrap(p, w)))
  m_lines <- lapply(main, wrapl, w = wrap_main)
  s_lines <- lapply(side, wrapl, w = wrap_side)
  pad <- 0.9
  m_h <- vapply(m_lines, length, 0) + pad
  s_h <- vapply(s_lines, length, 0) + pad
  gap <- pmax(2.4, c(s_h[seq_len(length(main) - 1)] + 1.0))
  total <- sum(m_h) + sum(gap) + 0.6
  unit_px <- 56
  save_figure(file, 2400, round(total * unit_px), pointsize = 11, draw = function() {
  par(mar = c(0, 0, 0, 0))
  plot.new(); plot.window(c(0, 1), c(-total, 0), xaxs = "i", yaxs = "i")
  mx <- c(0.04, 0.50); sx <- c(0.60, 0.97); ax <- 0.20
  y <- -0.3
  for (i in seq_along(main)) {
    top <- y; bottom <- y - m_h[i]
    rect(mx[1], bottom, mx[2], top, border = "black", lwd = 1)
    text(mx[1] + 0.015, top - pad / 2 - (seq_along(m_lines[[i]]) - 0.5), m_lines[[i]], adj = c(0, 0.5), cex = 0.78)
    if (i < length(main)) {
      next_top <- bottom - gap[i]
      arrows(ax, bottom - 0.25, ax, next_top + 0.25, length = 0.06, lwd = 1)
      if (length(s_lines[[i]])) {
        mid <- (bottom + next_top) / 2
        arrows(ax + 0.03, mid, sx[1] - 0.03, mid, length = 0.06, lwd = 1)
        st <- mid + s_h[i] / 2
        rect(sx[1], mid - s_h[i] / 2, sx[2], st, border = "black", lwd = 1)
        text(sx[1] + 0.015, st - pad / 2 - (seq_along(s_lines[[i]]) - 0.5), s_lines[[i]], adj = c(0, 0.5), cex = 0.78)
      }
      y <- next_top
    }
  }
  })
}

# Selection of the study sample (Fig. 1). Steps that exclude nobody are
# merged into the box above; n_alt: outcome sample on the unconverted scale.
plot_flow <- function(flow, file, n_alt = NA) {
  n <- setNames(flow$n_remaining, flow$step)
  x <- setNames(flow$n_excluded, flow$step)
  get <- function(k) if (k %in% names(n)) n[[k]] else NA
  ex <- function(k) if (k %in% names(x)) x[[k]] else 0
  main <- c(
    sprintf("Ever-married women aged 15–49 years interviewed in BDHS 2022\n(n = %s)", comma(get("Women interviewed"))),
    sprintf("Currently pregnant women\n(n = %s)", comma(get("Currently pregnant"))),
    sprintf("Pregnant women in households selected for biomarker testing\n(n = %s)",
            comma(get("In household selected for biomarker testing"))),
    sprintf("Analytic sample: aged 18 years or older, fasted at least 8 hours, valid fasting plasma glucose\n(n = %s)",
            comma(get("Fasted at least 8 hours (analytic sample)"))),
    paste0(sprintf("Outcome and model sample, released scale\n(n = %s; no missing covariates)",
                   comma(get("Complete covariate data (model sample)"))),
           if (!is.na(n_alt)) sprintf("\n\nOutcome and model sample, unconverted scale\n(n = %s; no missing covariates)", comma(n_alt)) else ""))
  side <- c(
    sprintf("Excluded %s women who were not currently pregnant", comma(ex("Currently pregnant"))),
    sprintf("Excluded %s women in households not selected for biomarker testing",
            comma(ex("In household selected for biomarker testing"))),
    paste0(sprintf("Excluded %s women aged under 18 years", comma(ex("Aged 18 years or older"))),
           sprintf("\n\nExcluded %s women without a valid fasting plasma glucose value",
                   comma(ex("Valid fasting plasma glucose value"))),
           if (ex("Fasted at least 8 hours (analytic sample)") > 0)
             sprintf("\n\nExcluded %s women who fasted less than 8 hours", comma(ex("Fasted at least 8 hours (analytic sample)"))) else ""),
    paste0(sprintf("Excluded %s women with self-reported pre-existing diabetes",
                   comma(ex("No self-reported pre-existing diabetes"))),
           sprintf("\n\nExcluded %s women with fasting plasma glucose of 7.0 mmol/L or above on the released scale",
                   comma(ex("FPG below 7.0 mmol/L (outcome sample)"))),
           if (!is.na(n_alt)) sprintf("\n\nExcluded %s women with fasting plasma glucose of 7.0 mmol/L or above on the unconverted scale",
                                      comma(get("No self-reported pre-existing diabetes") - n_alt)) else ""))
  flow_chart(main, side, file)
}

# Two panels: concentration curves and decomposition (absolute contributions,
# because percentage shares are not informative when the index is close to zero)
plot_inequality <- function(curve, decomp, file) {
  save_figure(file, 3300, 1600, function() {
  layout(matrix(1:2, 1), widths = c(1, 1.4))
  par(mar = c(4, 4, 2.2, 1))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i",
       xlab = "Cumulative share of women, poorest to richest", ylab = "Cumulative share of GDM cases")
  abline(0, 1, lty = 3, col = "grey30")
  for (s in names(scale_cols)) {
    r <- curve[curve$scale == s, ]
    lines(r$cum_women, r$cum_gdm, lwd = 1.6, col = scale_cols[[s]], lty = scale_lty[[s]])
  }
  legend("topleft", c(scale_label[names(scale_cols)], "Line of equality"), lty = c(scale_lty, 3),
         lwd = c(1.6, 1.6, 1), col = c(scale_cols, "grey30"), bty = "n", cex = 0.8)
  panel_letter("A")
  ref <- decomp[decomp$scale == "released", ]
  ord <- c("residual", ref$term[ref$term != "residual"][order(ref$contribution[ref$term != "residual"])])
  n <- length(ord)
  par(mar = c(4, 13, 2.6, 1))
  plot(NA, xlim = range(c(decomp$contribution_ci_lower, decomp$contribution_ci_upper, 0), na.rm = TRUE),
       ylim = c(0.5, n + 0.5), yaxt = "n", ylab = "", xlab = "Contribution to the Erreygers index (95% CI)")
  abline(v = 0, lty = 3, col = "grey30")
  off <- c(released = 0.15, unconverted = -0.15)
  for (s in names(scale_cols)) {
    r <- decomp[decomp$scale == s, ]
    y <- match(r$term, ord) + off[[s]]
    segments(r$contribution_ci_lower, y, r$contribution_ci_upper, y, col = scale_cols[[s]], lwd = 1.2)
    points(r$contribution, y, pch = scale_pch[[s]], col = scale_cols[[s]], bg = "white", cex = 0.9)
  }
  axis(2, at = seq_len(n), labels = ifelse(ord == "residual", "Residual (unexplained)", term_label(ord)),
       cex.axis = 0.8, tick = FALSE)
  legend("right", scale_label[names(scale_cols)], pch = scale_pch, col = scale_cols, pt.bg = "white",
         lwd = 1.2, bty = "n", cex = 0.8)
  panel_letter("B")
  })
}

plot_calibration <- function(cal, file) {
  save_figure(file, 3200, 1650, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 2.2, 1))
  for (s in names(scale_cols)) {
    r <- cal[cal$scale == s, ]
    plot(r$predicted, r$observed, xlim = c(0, 1), ylim = c(0, 1), pch = 19, col = "black", cex = 0.8,
         xaxs = "i", yaxs = "i", xlab = "Mean predicted probability (decile of predicted risk)",
         ylab = "Observed weighted proportion with GDM")
    segments(r$predicted, r$ci_lower, r$predicted, r$ci_upper, col = "black")
    abline(0, 1, lty = 3, col = "grey30")
    panel_letter(LETTERS[match(s, names(scale_cols))], scale_label[[s]])
  }
  })
}

# Weighted FPG distributions: empirical cumulative distribution in pregnant
# and non-pregnant women of the same age range, with the two thresholds, on
# each glucose scale.
fpg_groups <- function(women) {
  bio <- biomarker_women(women)
  preg <- women$v213 %in% 1
  age_range <- range(women$age[bio & preg])
  list(pregnant = bio & preg,
       non_pregnant = bio & !preg & women$age >= age_range[1] & women$age <= age_range[2])
}

plot_fpg_ecdf <- function(women, file) {
  groups <- fpg_groups(women)
  g1 <- groups$pregnant
  g2 <- groups$non_pregnant
  wecdf <- function(x, w) { ok <- !is.na(x); x <- x[ok]; w <- w[ok]; o <- order(x)
    list(x = x[o], y = cumsum(w[o]) / sum(w)) }
  save_figure(file, 3200, 1500, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 2.2, 1))
  for (v in c("fpg", "fpg_alt")) {
    s <- if (v == "fpg") "released" else "unconverted"
    a <- wecdf(women[[v]][g1], women$weight[g1]); b <- wecdf(women[[v]][g2], women$weight[g2])
    plot(a$x, a$y, type = "s", lwd = 1.6, col = "black", xlim = c(3, 10), ylim = c(0, 1), yaxs = "i",
         xlab = "Fasting plasma glucose (mmol/L)", ylab = "Cumulative weighted proportion")
    lines(b$x, b$y, type = "s", lwd = 1.6, col = "grey50", lty = 2)
    abline(v = c(5.1, 7.0), lty = 3, col = "grey30")
    legend("bottomright", c("Pregnant women", "Non-pregnant women, same age range"), lwd = 1.6, lty = c(1, 2),
           col = c("black", "grey50"), bty = "n", cex = 0.8)
    panel_letter(if (s == "released") "A" else "B", scale_label[[s]])
  }
  })
}

# Adjusted PRs across fasting thresholds on the unconverted scale (one panel
# per exposure), with the released-scale threshold and the 5.1 mmol/L threshold marked.
plot_threshold_sweep <- function(sweep, file) {
  save_figure(file, 3200, 1500, function() {
  par(mfrow = c(1, 2), mar = c(4, 4, 2.2, 1))
  panels <- c(wealthHighest = "Highest vs lowest wealth tertile", residenceUrban = "Urban vs rural residence")
  for (tm in names(panels)) {
    s <- sweep[sweep$term == tm, ]
    plot(s$threshold_unconverted, s$estimate, log = "y", pch = 19, col = "black", cex = 0.8,
         ylim = range(c(sweep$ci_lower, sweep$ci_upper, 1)),
         xlab = "FPG threshold, unconverted scale (mmol/L)", ylab = "Adjusted prevalence ratio (95% CI)")
    segments(s$threshold_unconverted, s$ci_lower, s$threshold_unconverted, s$ci_upper, col = "black")
    abline(h = 1, lty = 3, col = "grey30")
    abline(v = c(82 / 18, 5.1), lty = 2, col = "grey50")
    panel_letter(LETTERS[match(tm, names(panels))], panels[[tm]])
  }
  })
}

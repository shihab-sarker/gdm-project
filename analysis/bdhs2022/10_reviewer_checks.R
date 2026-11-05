# 10_reviewer_checks.R
# Analyses requested by the external reviews.
# (1) Scale-independent gradient: survey-weighted linear regression of log FPG.
#     The released value is the unconverted value x 1.11 (then rounded), so log
#     FPG differs between scales by a constant: coefficients other than the
#     intercept agree up to the rounding of sb267g.
# (2) Threshold sweep: adjusted PRs for the highest wealth tertile and urban
#     residence with GDM defined on the unconverted scale at thresholds from
#     4.4 to 5.4 mmol/L. sb267g is rounded, so released >= 5.1 is exactly
#     sb267 >= 82 mg/dL (unconverted 4.56); this point is added to the grid.

source(here::here("src", "common", "setup.R"))
women  <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))
x <- cfg$covariates
base <- sample[sample$pre_dm %in% 0 & !is.na(sample$fpg) & !is.na(sample$fpg_alt), ]

# ---- (1) log-FPG regression ----------------------------------------------------
fit_log <- function(d, y, label) {
  d$log_fpg <- log(d[[y]])
  m <- svyglm(reformulate(x, "log_fpg"), design = svy_domain(d, women))
  ci <- confint(m, ddf = degf(m$survey.design))
  data.frame(scale = label, n = nrow(d), term = names(coef(m)), coefficient = coef(m),
             pct_difference = 100 * (exp(coef(m)) - 1),
             ci_lower = 100 * (exp(ci[, 1]) - 1), ci_upper = 100 * (exp(ci[, 2]) - 1),
             p_value = summary(m)$coefficients[, 4], row.names = NULL)
}
log_fpg <- do.call(rbind, lapply(names(scales), function(s) {
  fit_log(base, if (s == "released") "fpg" else "fpg_alt", s)
}))
print(log_fpg)
save_table(log_fpg, "reviewer_log_fpg")

# ---- (2) threshold sweep -------------------------------------------------------
sweep_rows <- base[base$fpg_alt < cfg$fpg_diabetes, ]
thresholds <- sort(c(round(seq(4.4, 5.4, by = 0.1), 1), 82 / 18))
sweep <- do.call(rbind, lapply(thresholds, function(t) {
  d <- sweep_rows
  d$gdm_t <- as.numeric(d$fpg_alt >= t)
  e <- effect_table(fit_pr(svy_domain(d, women), x, "gdm_t"), "PR", "adjusted")
  e <- e[e$term %in% c("wealthHighest", "residenceUrban"), c("term", "estimate", "ci_lower", "ci_upper", "p_value")]
  cbind(threshold_unconverted = t, released_equivalent = t * 1.11, n = nrow(d), events = sum(d$gdm_t),
        prevalence = weighted.mean(d$gdm_t, d$weight), e, row.names = NULL)
}))
print(sweep)
save_table(sweep, "reviewer_threshold_sweep")
plot_threshold_sweep(sweep, file.path(dir_figures, "Figure_SA2.png"))

# ---- (3) early-pregnancy glucose ------------------------------------------------
by_month <- do.call(rbind, lapply(names(scales), function(s) {
  cbind(scale = s, fpg_by_month(base, women, if (s == "released") "fpg" else "fpg_alt", cfg))
}))
print(by_month)
save_table(by_month, "reviewer_fpg_by_month")
fasting <- fasting_by_trimester(base, women)
print(fasting)
save_table(fasting, "reviewer_fasting_by_trimester")

# ---- (4) FPG by trimester and slope per completed month ---------------------------
by_trimester <- do.call(rbind, lapply(names(scales), function(s) {
  d <- base
  d$fpg_x <- d[[if (s == "released") "fpg" else "fpg_alt"]]
  des <- svy_domain(d, women)
  do.call(rbind, lapply(levels(d$trimester), function(t) {
    q <- coef(svyquantile(~fpg_x, subset(des, trimester == t), c(0.25, 0.5, 0.75), ci = FALSE))
    data.frame(scale = s, trimester = t, n = sum(d$trimester %in% t), q25 = q[1], median = q[2], q75 = q[3], row.names = NULL)
  }))
}))
print(by_trimester)
save_table(by_trimester, "reviewer_fpg_by_trimester")

# Percentage change in FPG per completed month: unadjusted, and adjusted for
# the covariates of the log-FPG model (trimester replaced by month)
slope <- do.call(rbind, lapply(names(scales), function(s) {
  d <- base
  d$log_fpg <- log(d[[if (s == "released") "fpg" else "fpg_alt"]])
  d$month <- pmin(d$v214, 9)
  des <- svy_domain(d, women)
  do.call(rbind, lapply(c("unadjusted", "adjusted"), function(a) {
    rhs <- if (a == "unadjusted") "month" else c("month", setdiff(x, "trimester"))
    m <- svyglm(reformulate(rhs, "log_fpg"), design = des)
    ci <- confint(m, ddf = degf(m$survey.design))["month", ]
    data.frame(scale = s, model = a, n = nrow(d), pct_per_month = 100 * (exp(coef(m)[["month"]]) - 1),
               ci_lower = 100 * (exp(ci[[1]]) - 1), ci_upper = 100 * (exp(ci[[2]]) - 1),
               p_value = summary(m)$coefficients["month", 4])
  }))
}))
print(slope)
save_table(slope, "reviewer_fpg_slope_by_month")

# Women under 8 hours since eating by the computed times, although flagged as fasted
under8 <- function(d, label) data.frame(group = label, n = nrow(d), n_with_times = sum(!is.na(d$hours_ate)),
                                        n_under_8 = sum(d$hours_ate < 8, na.rm = TRUE))
fasting_times <- rbind(under8(women[biomarker_women(women), ], "All women tested, flagged fasted"),
                       under8(sample, "Analytic sample"), under8(base, "Women without known diabetes"))
print(fasting_times)
save_table(fasting_times, "reviewer_fasting_times")

# ---- (5) relative concentration index of continuous FPG -------------------------------
boot <- bootstrap_weights(sample, women, cfg)
rows <- sample$pre_dm %in% 0 & !is.na(sample$fpg) & !is.na(sample$fpg_alt)
ci_fpg <- do.call(rbind, lapply(names(scales), function(s) {
  y <- sample[[if (s == "released") "fpg" else "fpg_alt"]][rows]
  r <- sample$wealth_score[rows]
  reps <- apply(boot$weights[rows, ], 2, function(w) { k <- w > 0; concentration_index(y[k], r[k], w[k]) })
  cbind(scale = s, n = length(y), bootstrap_summary(concentration_index(y, r, sample$weight[rows]), reps, boot))
}))
print(ci_fpg)
save_table(ci_fpg, "reviewer_concentration_fpg")

# ---- (6) women beyond the x-axis of Fig. SA3 ----------------------------------------
groups <- fpg_groups(women)
above_10 <- do.call(rbind, lapply(c(released = "fpg", unconverted = "fpg_alt"), function(v) {
  do.call(rbind, lapply(names(groups), function(g) {
    x <- women[[v]][groups[[g]]]
    data.frame(scale = if (v == "fpg") "released" else "unconverted", group = g,
               n = sum(!is.na(x)), n_above_10 = sum(x > 10, na.rm = TRUE))
  }))
}))
print(above_10)
save_table(above_10, "reviewer_fpg_above_10")

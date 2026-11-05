# src/common/validity_checks.R
# Checks on the glucose outcome and on selection into the analytic sample.

# Women with a valid fasting glucose in the biomarker subsample (any
# pregnancy status): household selected, aged 18+, fasted >= 8 h.
biomarker_women <- function(women) {
  women$shbpbg %in% 1 & women$age >= 18 & women$fasted & !is.na(women$fpg)
}

# Weighted FPG distribution in pregnant and non-pregnant women of the same
# biomarker subsample. If the high prevalence reflected pregnancy, FPG would
# be lower, not higher, in non-pregnant women (FPG falls in early pregnancy;
# Mills et al. 1998).
fpg_distribution <- function(women, var = "fpg") {
  women$fpg_x <- women[[var]]
  bio <- biomarker_women(women)
  preg <- women$v213 %in% 1
  age_range <- range(women$age[bio & preg])
  groups <- list(
    "Pregnant" = bio & preg,
    "Non-pregnant" = bio & !preg,
    "Non-pregnant, same age range" = bio & !preg & women$age >= age_range[1] & women$age <= age_range[2])
  des <- svydesign(ids = ~v021, strata = ~v022, weights = ~weight, data = women)
  do.call(rbind, lapply(names(groups), function(g) {
    d <- subset(des, groups[[g]])
    q <- coef(svyquantile(~fpg_x, d, c(0.25, 0.5, 0.75), ci = FALSE))
    p51 <- svymean(~I(fpg_x >= 5.1), d)
    p70 <- svymean(~I(fpg_x >= 7.0), d)
    data.frame(scale = var, group = g, n = sum(groups[[g]]), q25 = q[1], median = q[2], q75 = q[3],
               pct_fpg_5.1_or_more = 100 * coef(p51)[2], se_5.1 = 100 * SE(p51)[2],
               pct_fpg_7.0_or_more = 100 * coef(p70)[2], se_7.0 = 100 * SE(p70)[2],
               row.names = NULL)
  }))
}

# Digit preference: share of values whose first decimal is 0 or 5
# (expected 20% if decimals were uniform), and the range of values.
glucose_heaping <- function(women) {
  bio <- biomarker_women(women)
  preg <- women$v213 %in% 1
  groups <- list("Pregnant" = bio & preg, "Non-pregnant" = bio & !preg)
  do.call(rbind, lapply(names(groups), function(g) {
    x <- women$fpg[groups[[g]]]
    data.frame(group = g, n = length(x),
               share_decimal_0_or_5 = mean(round(x * 10) %% 10 %in% c(0, 5)),
               expected_if_uniform = 0.2, min = min(x), max = max(x))
  }))
}

# Prevalence by completed month of pregnancy (months 9 and above pooled)
prevalence_by_month <- function(sample, women, outcome = "gdm") {
  sample$gest_month <- factor(pmin(sample$v214, 9))
  cbind(outcome = outcome, prevalence_by(svy_domain(sample, women), outcome, "gest_month"))
}

# Timing of the test by socioeconomic position: weighted distribution of
# trimester (completed months 1-3, 4-6, 7+) by residence and wealth tertile,
# with the Rao-Scott test, among the women of the main outcome sample.
timing_by_ses <- function(sample, women) {
  des <- subset(svy_domain(sample, women), !is.na(gdm))
  do.call(rbind, lapply(c("residence", "wealth"), function(v) {
    tab <- svytable(reformulate(c(v, "trimester")), des)
    n <- table(des$variables[[v]], des$variables$trimester)
    pct <- 100 * prop.table(tab, 1)
    p <- suppressWarnings(svychisq(reformulate(c(v, "trimester")), des, statistic = "F")$p.value)
    mean_m <- svyby(~gest_month, reformulate(v), des, svymean)
    data.frame(variable = v, level = rownames(pct), n = as.numeric(rowSums(n)),
               first = pct[, "1st"], second = pct[, "2nd"], third = pct[, "3rd"],
               mean_completed_months = mean_m$gest_month, p_value = as.numeric(p), row.names = NULL)
  }))
}

# Replication of the BDHS 2022 final report (NIPORT & ICF 2024, Table 14.6.1):
# women aged 18+ in biomarker households with a valid glucose value,
# household weight hv005. Diabetes = FPG >= 7.0 mmol/L or taking medication
# for diabetes (sb240 = 1). Computed on both glucose scales; the report
# gives 16.8% for women in 2022 (and 10% in 2017-18).
dhs_report_replication <- function(pr) {
  ok <- pr$hv104 %in% 2 & pr$hv105 >= 18 & pr$hv105 < 95 & pr$shbpbg %in% 1 &
    !pr$sb267 %in% c(NA, 994:996) & !pr$sb267g %in% c(NA, 994:996)
  p <- pr[ok, ]
  w <- p$hv005 / 1e6
  med <- p$sb240 %in% 1
  scales <- list(sb267g = p$sb267g, sb267_div18 = p$sb267 / 18)
  do.call(rbind, lapply(names(scales), function(s) {
    fpg <- scales[[s]]
    pct <- function(x) 100 * weighted.mean(x, w)
    data.frame(scale = s, n = nrow(p), weighted_n = sum(w),
               diabetes = pct(fpg >= 7 | med),
               fpg_7_or_more_no_medication = pct(fpg >= 7 & !med),
               fpg_6.1_6.9_no_medication = pct(fpg >= 6.1 & fpg < 7 & !med),
               fpg_3.9_6.0_no_medication = pct(fpg >= 3.9 & fpg < 6.1 & !med),
               diabetes_urban = 100 * weighted.mean((fpg >= 7 | med)[p$hv025 %in% 1], w[p$hv025 %in% 1]),
               diabetes_rural = 100 * weighted.mean((fpg >= 7 | med)[p$hv025 %in% 2], w[p$hv025 %in% 2]),
               report_2022 = 16.8, report_2017_18 = 10)
  }))
}

# Women included and excluded at each step of the participant flow
# (weighted means / proportions and standardised differences). Descriptive
# only: small excluded groups cannot show the absence of selection bias.
selection_comparison <- function(women, cfg) {
  selected <- women$v213 %in% 1 & women$shbpbg %in% 1
  adult    <- selected & women$age >= 18
  glucose  <- adult & !is.na(women$fpg)
  analytic <- glucose & women$fasted
  outcome  <- analytic & !is.na(women$gdm)
  steps <- list(
    "Aged under 18 (vs 18+)"            = list(pool = selected, kept = adult),
    "No valid glucose value"            = list(pool = adult, kept = glucose),
    "Pre-existing diabetes or FPG >= 7" = list(pool = analytic, kept = outcome))

  women$wealth_quintile <- factor(women$v190)
  des <- svydesign(ids = ~v021, strata = ~v022, weights = ~weight, data = women)
  describe <- function(d) {
    a <- svymean(~age, d, na.rm = TRUE)
    rows <- data.frame(variable = "age", level = "mean", estimate = coef(a),
                       var = coef(svyvar(~age, d, na.rm = TRUE)))
    for (v in c("residence", "education", "wealth_quintile")) {
      m <- svymean(reformulate(v), d, na.rm = TRUE)
      p <- coef(m)
      rows <- rbind(rows, data.frame(variable = v, level = sub(v, "", names(p)),
                                     estimate = p, var = p * (1 - p)))
    }
    rows
  }
  do.call(rbind, lapply(names(steps), function(s) {
    st <- steps[[s]]
    inc <- st$pool & st$kept
    exc <- st$pool & !st$kept
    if (sum(exc) < 2) return(NULL)
    di <- describe(subset(des, inc))
    de <- describe(subset(des, exc))
    data.frame(step = s, n_included = sum(inc), n_excluded = sum(exc),
               variable = di$variable, level = di$level,
               included = di$estimate, excluded = de$estimate,
               std_difference = (di$estimate - de$estimate) / sqrt((di$var + de$var) / 2),
               row.names = NULL)
  }))
}

# Each woman classified on both glucose scales (women without known diabetes
# and with a valid glucose value): unweighted n and survey-weighted % per cell.
paired_classification <- function(sample, women, cfg) {
  band <- function(x) cut(x, c(-Inf, cfg$fpg_gdm, cfg$fpg_diabetes, Inf), right = FALSE,
                          labels = c("Normal", "GDM range", "Diabetes range"))
  d <- sample[sample$pre_dm %in% 0 & !is.na(sample$fpg) & !is.na(sample$fpg_alt), ]
  d$released <- band(d$fpg)
  d$unconverted <- band(d$fpg_alt)
  d$cell <- paste(d$released, d$unconverted, sep = " / ")
  des <- svy_domain(d, women)
  cells <- as.data.frame(table(released = d$released, unconverted = d$unconverted), responseName = "n")
  cells <- cells[cells$n > 0, ]
  do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
    key <- paste(cells$released[i], cells$unconverted[i], sep = " / ")
    des$variables$in_cell <- as.numeric(des$variables$cell %in% key)
    p <- svyciprop(~in_cell, des, method = "logit")
    data.frame(cells[i, ], pct = 100 * as.numeric(p), ci_lower = 100 * confint(p)[1],
               ci_upper = 100 * confint(p)[2], row.names = NULL)
  }))
}

# The same women on both scales, so that the two scales differ only in the
# threshold: (a) women below 7.0 mmol/L on both scales; (b) all women without
# known diabetes, with the diabetes range counted as GDM.
common_sample <- function(sample, women, cfg) {
  des <- svy_domain(sample, women)
  sets <- list(
    "Below 7.0 mmol/L on both scales" = list(rows = !is.na(sample$gdm) & !is.na(sample$gdm_alt),
                                             outcome = c(released = "gdm", unconverted = "gdm_alt")),
    "Diabetes range counted as GDM" = list(rows = !is.na(sample$gdm_incl_dm) & !is.na(sample$gdm_incl_dm_alt),
                                           outcome = c(released = "gdm_incl_dm", unconverted = "gdm_incl_dm_alt")))
  do.call(rbind, lapply(names(sets), function(k) {
    st <- sets[[k]]
    d <- subset(des, st$rows)
    do.call(rbind, lapply(names(st$outcome), function(s) {
      o <- st$outcome[[s]]
      e <- effect_table(fit_pr(d, cfg$covariates, o), "PR", "adjusted")
      e <- e[e$term %in% c("wealthHighest", "residenceUrban"), c("term", "estimate", "ci_lower", "ci_upper", "p_value")]
      p <- weighted_prevalence(d, o)[, c("n", "events", "prevalence", "ci_lower", "ci_upper")]
      names(p)[4:5] <- c("prevalence_ci_lower", "prevalence_ci_upper")
      cbind(sample = k, scale = s, p, e, row.names = NULL)
    }))
  }))
}

# FPG by completed month of pregnancy on one scale: weighted median and IQR,
# and the weighted share within 0.3 mmol/L of the 5.1 mmol/L threshold.
fpg_by_month <- function(d, women, var, cfg) {
  d$fpg_x <- d[[var]]
  d$gest_month <- pmin(d$v214, 9)
  d$near <- as.numeric(round(abs(d$fpg_x - cfg$fpg_gdm), 6) <= 0.3)
  des <- svy_domain(d, women)
  do.call(rbind, lapply(sort(unique(d$gest_month)), function(m) {
    s <- subset(des, gest_month == m)
    q <- coef(svyquantile(~fpg_x, s, c(0.25, 0.5, 0.75), ci = FALSE))
    data.frame(month = m, n = sum(d$gest_month == m), q25 = q[1], median = q[2], q75 = q[3],
               pct_within_0.3 = 100 * coef(svymean(~near, s)), row.names = NULL)
  }))
}

# Fasting by approximate trimester: hours from the last food and the last
# drink other than plain water to the test (weighted median and IQR).
fasting_by_trimester <- function(d, women) {
  des <- svy_domain(d, women)
  quart <- function(s, v) {
    x <- subset(s, !is.na(s$variables[[v]]))
    coef(svyquantile(reformulate(v), x, c(0.25, 0.5, 0.75), ci = FALSE))
  }
  do.call(rbind, lapply(levels(d$trimester), function(t) {
    s <- subset(des, trimester == t)
    a <- quart(s, "hours_ate")
    g <- d$trimester %in% t
    out <- data.frame(trimester = t, n = sum(g), ate_q25 = a[1], ate_median = a[2], ate_q75 = a[3],
                      n_ate_under_8 = sum(d$hours_ate[g] < 8, na.rm = TRUE), n_drank_time = sum(!is.na(d$hours_drank[g])),
                      pct_drank_never = 100 * coef(svymean(~as.numeric(drank_never), s)),
                      pct_drank_within_8h = 100 * coef(svymean(~as.numeric(drank_8h), s)),
                      n_drank_within_8h = sum(d$drank_8h[g]), row.names = NULL)
    h <- if (out$n_drank_time > 1) quart(s, "hours_drank") else rep(NA, 3)
    cbind(out, drank_q25 = h[1], drank_median = h[2], drank_q75 = h[3])
  }))
}

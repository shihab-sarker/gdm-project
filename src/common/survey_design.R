# src/common/survey_design.R
# Complex survey design of BDHS 2022 (two-stage stratified cluster sample):
# PSU = v021, strata = v022, weight = v005 / 1e6.
#
# Domain (subpopulation) estimation. The analytic sample is a domain of the
# full sample, so variances must be computed over all PSUs of the survey,
# including PSUs with no woman in the analytic sample (Korn & Graubard 1999;
# Lumley 2004). Such PSUs contribute zero to every domain total; only the
# fact that they exist in their stratum matters. The design is therefore
# built from the analytic rows plus one placeholder row per PSU that has no
# analytic woman, and then restricted to the analytic rows. This gives the
# same standard errors as a design on all 30,078 women while keeping the
# bootstrap replicate-weight matrix small.

domain_frame <- function(sample, women) {
  others <- women[!women$v021 %in% sample$v021, c("v021", "v022", "weight")]
  others <- others[!duplicated(others$v021), ]
  frame <- dplyr::bind_rows(sample, others)
  frame$in_domain <- c(rep(TRUE, nrow(sample)), rep(FALSE, nrow(others)))
  frame
}

# Linearisation (Taylor series) design for the domain
svy_domain <- function(sample, women) {
  frame <- domain_frame(sample, women)
  des <- svydesign(ids = ~v021, strata = ~v022, weights = ~weight, data = frame)
  subset(des, in_domain)
}

# Design on the analytic rows only (no domain), used for comparison and by
# diagnostics that need the design variables on every row.
svy_restricted <- function(sample) {
  svydesign(ids = ~v021, strata = ~v022, weights = ~weight, data = sample)
}

# Rao-Wu rescaled bootstrap (Rao & Wu 1988; Preston 2009): PSUs are resampled
# within strata, n_h - 1 out of n_h, and the weights rescaled. Returns the
# replicate weights of the analytic rows (one column per replicate).
# The draws are those of as.svrepdesign(type = "subbootstrap") in the same random
# number order, so results match it exactly; PSU counts use tabulate(), because
# table() on every replicate of every stratum took almost all of its run time.
bootstrap_weights <- function(sample, women, cfg) {
  frame <- domain_frame(sample, women)
  des <- svydesign(ids = ~v021, strata = ~v022, weights = ~weight, data = frame)
  B <- cfg$bootstrap_B
  psu <- des$cluster[, 1]
  first <- !duplicated(psu)
  upsu <- psu[first]
  ustrata <- as.character(des$strata[, 1])[first]
  multiplier <- matrix(NA_real_, length(upsu), B)
  set.seed(cfg$seed)
  for (s in unique(ustrata)) {
    rows <- which(ustrata == s)
    n <- length(rows)
    drawn <- sample(upsu[rows])
    multiplier[rows, ] <- vapply(seq_len(B), function(b) {
      tabulate(match(sample(drawn, n - 1, replace = TRUE), upsu[rows]), n)
    }, numeric(n)) * n / (n - 1)
  }
  rep_w <- multiplier[match(psu, upsu), , drop = FALSE] * (1 / des$prob)
  list(weights = rep_w[frame$in_domain, , drop = FALSE], scale = 1 / (B - 1), rscales = rep(1, B))
}

# Bootstrap standard error and percentile 95% CI from replicate estimates
bootstrap_summary <- function(estimate, replicates, boot) {
  ok <- !is.na(replicates)
  se <- sqrt(svrVar(replicates[ok], boot$scale, boot$rscales[ok], mse = FALSE, coef = estimate))
  data.frame(estimate = estimate, boot_se = as.numeric(se),
             ci_lower = unname(quantile(replicates[ok], 0.025)),
             ci_upper = unname(quantile(replicates[ok], 0.975)),
             B = length(replicates), failed = sum(!ok))
}

# Weighted prevalence with a logit-transformed 95% CI (svyciprop)
weighted_prevalence <- function(design, var) {
  d <- subset(design, !is.na(design$variables[[var]]))
  y <- d$variables[[var]][weights(d) > 0]
  n <- length(y)
  if (sum(y) %in% c(0, n)) {   # all or no women with the outcome: no CI
    return(data.frame(n = n, events = sum(y), prevalence = mean(y),
                      se = NA, ci_lower = NA, ci_upper = NA))
  }
  p <- svyciprop(reformulate(var), d, method = "logit")
  data.frame(n = n, events = sum(y), prevalence = as.numeric(p),
             se = as.numeric(SE(p)), ci_lower = confint(p)[1], ci_upper = confint(p)[2])
}

prevalence_by <- function(design, var, by) {
  lev <- levels(factor(design$variables[[by]]))
  do.call(rbind, lapply(lev, function(l) {
    cbind(group = by, level = l,
          weighted_prevalence(subset(design, design$variables[[by]] %in% l), var))
  }))
}

# Design effect of the prevalence (Kish 1965): variance under the survey
# design divided by the variance under simple random sampling of the same n.
# The weighting component is Kish's 1 + CV^2 of the weights; the remainder is
# attributed to clustering and stratification.
design_effect <- function(design, sample, var) {
  m <- svymean(reformulate(var), design, na.rm = TRUE, deff = "replace")
  ok <- !is.na(sample[[var]])
  w <- sample$weight[ok]
  deff_weights <- sum(ok) * sum(w^2) / sum(w)^2
  data.frame(n = sum(ok), prevalence = as.numeric(coef(m)), se = as.numeric(SE(m)),
             deff = as.numeric(deff(m)), deff_weighting = deff_weights,
             deff_clustering = as.numeric(deff(m)) / deff_weights,
             effective_n = sum(ok) / as.numeric(deff(m)))
}

# Number of PSUs per stratum in the full sample and in the analytic domain
psu_structure <- function(sample, women) {
  count <- function(d) tapply(d$v021, factor(d$v022, levels = sort(unique(women$v022))),
                              function(x) length(unique(x)))
  count0 <- function(d) { x <- count(d); x[is.na(x)] <- 0; x }
  summ <- function(x, scope) data.frame(scope = scope, strata = length(x), min = min(x),
                                        median = median(x), max = max(x),
                                        strata_with_one_psu = sum(x == 1),
                                        strata_with_no_psu = sum(x == 0))
  rbind(summ(count0(women), "full sample"), summ(count0(sample), "analytic sample"))
}

# Standard errors from the domain design versus a design built on the
# analytic rows only (prevalence and logistic-model coefficients).
design_comparison <- function(sample, women, cfg) {
  f <- reformulate(cfg$covariates, "gdm")
  d_dom <- svy_domain(sample, women)
  d_res <- svy_restricted(sample)
  p_dom <- svymean(~gdm, d_dom, na.rm = TRUE)
  p_res <- svymean(~gdm, d_res, na.rm = TRUE)
  m_dom <- svyglm(f, d_dom, family = quasibinomial())
  m_res <- svyglm(f, d_res, family = quasibinomial())
  out <- data.frame(term = c("prevalence", names(coef(m_dom))),
                    se_domain_design = c(SE(p_dom), SE(m_dom)),
                    se_restricted_design = c(SE(p_res), SE(m_res)))
  out$pct_difference <- 100 * (out$se_domain_design / out$se_restricted_design - 1)
  out
}

# ---- Table 1: characteristics of women with and without GDM -----------------
# Categorical: unweighted n, weighted column % and weighted % with GDM, with
# the Rao-Scott corrected chi-square test (F form; Rao & Scott 1984).
# Continuous: weighted mean and SD, design-based t test.

describe_categorical <- function(design, vars, outcome = "gdm") {
  d <- subset(design, !is.na(design$variables[[outcome]]))
  do.call(rbind, lapply(vars, function(v) {
    x <- d$variables[[v]]
    y <- d$variables[[outcome]]
    tab_w <- svytable(reformulate(c(v, outcome)), d)
    tab_n <- table(x, y)
    # Rao-Scott test; flag tables with an expected count below 5
    expected <- outer(rowSums(tab_n), colSums(tab_n)) / sum(tab_n)
    test <- suppressWarnings(svychisq(reformulate(c(v, outcome)), d, statistic = "F"))
    data.frame(outcome = outcome, variable = v, level = rownames(tab_w),
               n = as.numeric(rowSums(tab_n)), n_gdm = as.numeric(tab_n[, "1"]),
               pct_of_sample = as.numeric(100 * rowSums(tab_w) / sum(tab_w)),
               pct_gdm = as.numeric(100 * tab_w[, "1"] / rowSums(tab_w)),
               n_missing = sum(is.na(x)), p_value = as.numeric(test$p.value),
               expected_count_below_5 = any(expected < 5))
  }))
}

describe_continuous <- function(design, vars, outcome = "gdm") {
  d <- subset(design, !is.na(design$variables[[outcome]]))
  d$variables$group <- factor(ifelse(d$variables[[outcome]] == 1, "GDM", "No GDM"))
  do.call(rbind, lapply(vars, function(v) {
    f <- reformulate(v)
    by_mean <- svyby(f, ~group, d, svymean, na.rm = TRUE)
    by_var <- svyby(f, ~group, d, svyvar, na.rm = TRUE)
    test <- svyttest(reformulate(outcome, v), d)
    data.frame(outcome = outcome, variable = v, group = c("all", as.character(by_mean$group)),
               mean = c(coef(svymean(f, d, na.rm = TRUE)), by_mean[[v]]),
               sd = sqrt(c(coef(svyvar(f, d, na.rm = TRUE)), by_var[[v]])),
               p_value = as.numeric(test$p.value))
  }))
}

# Weighted distribution of residence across wealth tertiles (Rao-Scott test)
residence_by_wealth <- function(design) {
  tab <- svytable(~residence + wealth, design)
  out <- as.data.frame(100 * prop.table(tab, 1))
  names(out)[3] <- "row_pct"
  out$p_value <- svychisq(~residence + wealth, design, statistic = "F")$p.value
  out
}

# ---- Table 1 by household wealth tertile ------------------------------------
# Distribution of the outcome sample (released scale) across wealth tertiles:
# unweighted n and survey-weighted column % per tertile, with the Rao-Scott
# corrected F test (categorical) or a design-based Wald F test (continuous).
# GDM on each scale is added as rows, each within its own outcome sample.

describe_by_wealth <- function(design, vars, cont_vars, outcomes, base = outcomes[[1]]) {
  d <- subset(design, !is.na(design$variables[[base]]))
  groups <- levels(droplevels(factor(d$variables$wealth)))
  cat_rows <- function(dd, v, label = v) {
    x <- dd$variables[[v]]; w <- dd$variables$wealth
    tab_w <- svytable(reformulate(c(v, "wealth")), dd)
    tab_n <- table(x, w)
    test <- suppressWarnings(svychisq(reformulate(c(v, "wealth")), dd, statistic = "F"))
    tot_w <- rowSums(tab_w)
    out <- data.frame(variable = label, level = rownames(tab_w), group = "Total",
                      n = as.numeric(rowSums(tab_n)), value = as.numeric(100 * tot_w / sum(tot_w)), sd = NA)
    for (g in groups) out <- rbind(out, data.frame(variable = label, level = rownames(tab_w), group = g,
      n = as.numeric(tab_n[, g]), value = as.numeric(100 * tab_w[, g] / sum(tab_w[, g])), sd = NA))
    out$p_value <- as.numeric(test$p.value)
    out
  }
  cat_part <- do.call(rbind, lapply(setdiff(vars, "wealth"), function(v) cat_rows(d, v)))
  cont_part <- do.call(rbind, lapply(cont_vars, function(v) {
    f <- reformulate(v)
    m_all <- svymean(f, d, na.rm = TRUE); v_all <- svyvar(f, d, na.rm = TRUE)
    by_m <- svyby(f, ~wealth, d, svymean, na.rm = TRUE); by_v <- svyby(f, ~wealth, d, svyvar, na.rm = TRUE)
    fit <- svyglm(reformulate("wealth", v), d)
    data.frame(variable = v, level = "Mean (SD)", group = c("Total", as.character(by_m$wealth)),
               n = NA, value = c(coef(m_all), by_m[[v]]), sd = sqrt(c(coef(v_all), by_v[[v]])),
               p_value = as.numeric(regTermTest(fit, ~wealth)$p))
  }))
  gdm_part <- do.call(rbind, lapply(names(outcomes), function(s) {
    o <- outcomes[[s]]
    dd <- subset(design, !is.na(design$variables[[o]]))
    dd$variables$gdm_status <- factor(ifelse(dd$variables[[o]] == 1, "GDM", "No GDM"), c("GDM", "No GDM"))
    r <- cat_rows(dd, "gdm_status", paste0("gdm_", s))
    r[r$level == "GDM", ]
  }))
  out <- rbind(cat_part, cont_part, gdm_part)
  rownames(out) <- NULL
  out
}

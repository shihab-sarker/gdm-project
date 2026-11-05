# src/common/regression.R
# Survey-weighted regression models of GDM and their diagnostics.
#
# Prevalence ratios (PR): modified Poisson regression (log link, robust
# variance; Zou 2004, Barros & Hirakata 2003). Because the outcome is common
# (about 60%), odds ratios would be far from the corresponding prevalence
# ratios, so the PR is the main effect measure. Odds ratios (OR) from
# logistic regression are reported for comparison with the anchor study.
# Variances are design-based (Taylor linearisation) and CIs use the design
# degrees of freedom (number of PSUs minus number of strata).

fit_pr <- function(design, covariates, outcome = "gdm") {
  f <- reformulate(covariates, outcome)
  # Starting values: log of the overall prevalence for the intercept, zero
  # for the other coefficients (helps the log-link model to converge).
  p0 <- weighted.mean(design$variables[[outcome]], weights(design), na.rm = TRUE)
  k <- ncol(model.matrix(f, model.frame(f, design$variables)))
  do.call(svyglm, list(formula = f, design = design, family = quasipoisson(link = "log"),
                       start = c(log(p0), rep(0, k - 1))))
}

fit_or <- function(design, covariates, outcome = "gdm") {
  svyglm(reformulate(covariates, outcome), design = design, family = quasibinomial())
}

# Exponentiated coefficients with 95% CIs on the design degrees of freedom
effect_table <- function(model, measure, label) {
  b <- coef(model)
  ci <- exp(confint(model, ddf = degf(model$survey.design)))
  out <- data.frame(model = label, measure = measure, term = names(b),
                    estimate = exp(b), ci_lower = ci[, 1], ci_upper = ci[, 2],
                    p_value = summary(model)$coefficients[, 4], n = length(model$y))
  out[out$term != "(Intercept)", ]
}

# Crude (one covariate at a time) estimates on the same women as the
# adjusted model, so crude and adjusted estimates are comparable.
crude_effects <- function(design, covariates, measure, outcome = "gdm") {
  keep <- complete.cases(design$variables[c(outcome, covariates)])
  d <- subset(design, keep)
  fit <- if (measure == "PR") fit_pr else fit_or
  do.call(rbind, lapply(covariates, function(v) effect_table(fit(d, v, outcome), measure, "crude")))
}

# Firth's penalised-likelihood logistic regression (Heinze & Schemper 2002)
# with the sampling weights. A check on small-sample bias of the logistic
# estimates; it does not use strata or PSUs.
firth_check <- function(sample, covariates, outcome = "gdm") {
  f <- reformulate(covariates, outcome)
  d <- sample[complete.cases(sample[all.vars(f)]), ]
  m <- logistf::logistf(f, data = d, weights = d$weight)
  data.frame(term = names(coef(m)), estimate = exp(coef(m)),
             ci_lower = exp(m$ci.lower), ci_upper = exp(m$ci.upper), p_value = m$prob)[-1, ]
}

# ---- Model diagnostics ----------------------------------------------------------

# Events per parameter (Peduzzi et al. 1996): the smaller outcome class
# divided by the number of coefficients other than the intercept.
events_per_parameter <- function(model) {
  y <- model$y
  data.frame(n = length(y), events = sum(y == 1), non_events = sum(y == 0),
             parameters = length(coef(model)) - 1,
             epv = min(sum(y == 1), sum(y == 0)) / (length(coef(model)) - 1))
}

# Collinearity. (1) Generalised VIF (Fox & Monette 1992) from the
# design-based covariance matrix of the logistic model; GVIF^(1/(2 df)) is
# comparable across terms with different df. (2) Survey-weighted VIF of
# Liao & Valliant (2012), which is the ordinary VIF times two design factors
# (zeta, varrho). The R-squared of each column on the others shows whether a
# large value comes from collinearity or from the design factors.
collinearity <- function(model_or, sample, cfg) {
  g <- car::vif(model_or)
  gvif <- data.frame(term = rownames(g), gvif = g[, 1], df = g[, 2], gvif_adjusted = g[, 3])
  m <- fit_or(svy_restricted(sample), cfg$covariates)
  X <- model.matrix(m)[, -1]
  lv <- svydiags::svyvif(mobj = m, X = X, w = weights(m$survey.design), stvar = "v022", clvar = "v021")
  lv <- lv[["Intercept adjusted"]]
  svy <- data.frame(column = rownames(lv), survey_vif = lv[, "svy.vif.m"],
                    ordinary_vif = lv[, "reg.vif.m"], zeta = lv[, "zeta"],
                    varrho = lv[, "varrho.m"], r_squared = lv[, "Rsq.m"])
  list(gvif = gvif, survey_vif = svy)
}

# Groups of predicted risk with (approximately) equal sums of weights
risk_groups <- function(p, w, g = 10) {
  o <- order(p)
  share <- numeric(length(p))
  share[o] <- cumsum(w[o]) / sum(w)
  pmin(ceiling(share * g), g)
}

# Archer-Lemeshow goodness-of-fit test for a survey logistic model
# (Archer & Lemeshow 2006): the residuals y - p are regressed on indicators
# of deciles of predicted risk and the decile coefficients are tested
# jointly with a design-based Wald F test.
archer_lemeshow <- function(sample, cfg, g = 10, outcome = "gdm") {
  keep <- complete.cases(sample[c(outcome, cfg$covariates)])
  d <- sample[keep, ]
  m <- fit_or(svy_restricted(d), cfg$covariates, outcome)
  d$p <- fitted(m)
  d$resid <- d[[outcome]] - d$p
  d$risk_group <- factor(risk_groups(d$p, d$weight, g))
  m_r <- svyglm(resid ~ 0 + risk_group, design = svy_restricted(d))
  t <- regTermTest(m_r, ~risk_group)
  data.frame(groups = g, F = as.numeric(t$Ftest), df1 = t$df, df2 = t$ddf, p_value = as.numeric(t$p))
}

# Weighted observed vs mean predicted risk by decile of predicted risk from
# the survey logistic model (data for the calibration plot)
calibration_groups <- function(sample, cfg, g = 10, outcome = "gdm") {
  keep <- complete.cases(sample[c(outcome, cfg$covariates)])
  d <- sample[keep, ]
  m <- fit_or(svy_restricted(d), cfg$covariates, outcome)
  d$p <- fitted(m)
  d$risk_group <- risk_groups(d$p, d$weight, g)
  des <- svy_restricted(d)
  do.call(rbind, lapply(sort(unique(d$risk_group)), function(k) {
    dk <- subset(des, risk_group == k)
    obs <- weighted_prevalence(dk, outcome)
    data.frame(group = k, n = obs$n, predicted = weighted.mean(d$p[d$risk_group == k], d$weight[d$risk_group == k]),
               observed = obs$prevalence, ci_lower = obs$ci_lower, ci_upper = obs$ci_upper)
  }))
}

# Area under the ROC curve from ranks (Mann-Whitney form)
auc <- function(y, p) {
  r <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Calibration intercept (calibration-in-the-large, slope fixed at 1) and
# calibration slope of predictions p for outcomes y (Steyerberg 2019).
calibration_stats <- function(y, p) {
  lp <- qlogis(p)
  one <- rep(1, length(y))
  c(intercept = unname(glm.fit(cbind(one), y, offset = lp, family = binomial())$coefficients[1]),
    slope = unname(glm.fit(cbind(one, lp), y, family = binomial())$coefficients[2]))
}

# Bootstrap internal validation (Harrell et al. 1996; Steyerberg 2019).
# The logistic model is refitted in B bootstrap samples; the optimism of each
# performance measure is its value in the bootstrap sample minus its value
# when the bootstrap model is applied to the original data. The model is
# unweighted: this assesses performance in the sample, not in the population.
# glm.fit() on rows of one model matrix gives the fits of glm() without building
# a model frame in each replicate.
bootstrap_validation <- function(sample, cfg, outcome = "gdm") {
  f <- reformulate(cfg$covariates, outcome)
  d <- sample[complete.cases(sample[all.vars(f)]), ]
  X <- model.matrix(f, d)
  y <- d[[outcome]]
  logit <- binomial()
  perf <- function(y, p) c(auc = auc(y, p), calibration_stats(y, p), brier = mean((y - p)^2))
  apparent <- perf(y, glm.fit(X, y, family = logit)$fitted.values)
  set.seed(cfg$seed)
  optimism <- replicate(cfg$bootstrap_B, {
    i <- sample(nrow(d), replace = TRUE)
    fb <- glm.fit(X[i, ], y[i], family = logit)
    b <- fb$coefficients
    b[is.na(b)] <- 0   # a covariate level absent from this bootstrap sample
    perf(y[i], fb$fitted.values) - perf(y, logit$linkinv(drop(X %*% b)))
  })
  data.frame(measure = names(apparent), apparent = apparent,
             optimism = rowMeans(optimism), corrected = apparent - rowMeans(optimism),
             B = cfg$bootstrap_B, row.names = NULL)
}

# Linearity of the continuous covariates on the logit scale, in the anchor
# specification (duration in weeks): Box-Tidwell terms x*log(x) for BMI and
# gestational weeks, and trimester indicators added to the linear duration
# term. Each added term is tested with a design-based Wald F test.
linearity_tests <- function(sample, women, cfg, outcome = "gdm") {
  base <- cfg$covariates_anchor
  d <- sample
  d$bmi_log_bmi <- d$bmi * log(d$bmi)
  d$weeks_log_weeks <- d$gest_weeks * log(d$gest_weeks)
  des <- svy_domain(d, women)
  test <- function(label, added) {
    m <- fit_or(des, c(base, added), outcome)
    t <- regTermTest(m, reformulate(added))
    data.frame(test = label, df = t$df, F = as.numeric(t$Ftest), p_value = as.numeric(t$p))
  }
  rbind(test("BMI: Box-Tidwell term (BMI x log BMI)", "bmi_log_bmi"),
        test("Gestational weeks: Box-Tidwell term (weeks x log weeks)", "weeks_log_weeks"),
        test("Trimester indicators added to linear duration", "trimester"))
}

# Influential observations under the survey design (Li & Valliant 2011):
# DFBETAS (flag |DFBETAS| above the svydiags cut-off) and modified Cook's
# distance (flag values above 3) for the main logistic model. The model is
# refitted without the flagged women and the adjusted PRs are compared.
influence_checks <- function(sample, women, cfg, outcome = "gdm") {
  keep <- complete.cases(sample[c(outcome, cfg$covariates)])
  d <- sample[keep, ]
  m <- fit_or(svy_restricted(d), cfg$covariates, outcome)
  dfb <- svydiags::svydfbetas(m, stvar = "v022", clvar = "v021")
  cook <- svydiags::svyCooksD(m, stvar = "v022", clvar = "v021")
  flag_dfb <- apply(abs(dfb$Dfbetas[-1, , drop = FALSE]), 2, max) > dfb$cutoff
  flag_cook <- as.numeric(cook) > 3
  terms <- c("wealthHighest", "residenceUrban", "trimester3rd")
  refit <- function(label, drop) {
    e <- effect_table(fit_pr(svy_domain(d[!drop, ], women), cfg$covariates, outcome), "PR", label)
    e <- e[e$term %in% terms, c("term", "estimate", "ci_lower", "ci_upper")]
    cbind(scenario = label, n_dropped = sum(drop), e)
  }
  rbind(refit("all women", rep(FALSE, nrow(d))),
        refit("without women flagged by DFBETAS", flag_dfb),
        refit("without women flagged by Cook's distance", flag_cook),
        refit("without women flagged by either", flag_dfb | flag_cook))
}

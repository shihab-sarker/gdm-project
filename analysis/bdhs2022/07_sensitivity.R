# 07_sensitivity.R
# Sensitivity analyses of prevalence, adjusted prevalence ratios and (where
# the sample changes) the Erreygers concentration index, each run on both
# glucose scales (released and unconverted).
#   S0  main analysis
#   S1  outcome includes the diabetes range (FPG 5.1 mmol/L or above)
#   S2  restricted to 6 or more completed months of gestation
#   S3  excluding completed months 1-2 (early recognition of pregnancy)
#   S4  age cut-off at 30 years
#   S5  continuous age instead of age group
#   S6  additional adjustment for education
#   S7  additional adjustment for hypertension
#   S8  measured BMI (no weight-gain correction) instead of estimated
#       pre-pregnancy BMI; trimester stays in the model
#   S9  excluding women whose test day was not a fasted-appointment day
#   S10 unweighted prevalence with an exact binomial CI
#   S11 excluding women who drank anything other than water within 8 hours
#   S12 national wealth quintiles (v190) in three groups instead of tertiles

source(here::here("src", "common", "setup.R"))
women  <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))
des  <- svy_domain(sample, women)
boot <- bootstrap_weights(sample, women, cfg)
x    <- cfg$covariates
all_rows <- rep(TRUE, nrow(sample))

scenarios <- list(
  S0 = list(label = "Main analysis", outcome = "gdm", covariates = x, rows = all_rows, ci = TRUE),
  S1 = list(label = "Outcome includes FPG 7.0 mmol/L or above", outcome = "gdm_incl_dm",
            covariates = x, rows = all_rows, ci = TRUE),
  S2 = list(label = "Restricted to 6 or more months of gestation", outcome = "gdm",
            covariates = setdiff(x, "trimester"), rows = sample$v214 >= 6, ci = TRUE),
  S3 = list(label = "Excluding months 1-2 of gestation", outcome = "gdm",
            covariates = x, rows = sample$v214 >= 3, ci = TRUE),
  S4 = list(label = "Age cut-off at 30 years", outcome = "gdm",
            covariates = sub("age_group", "age_group_alt", x), rows = all_rows, ci = FALSE),
  S5 = list(label = "Continuous age", outcome = "gdm",
            covariates = sub("age_group", "age_c", x), rows = all_rows, ci = FALSE),
  S6 = list(label = "Additionally adjusted for education", outcome = "gdm",
            covariates = c(x, "education"), rows = all_rows, ci = FALSE),
  S7 = list(label = "Additionally adjusted for hypertension", outcome = "gdm",
            covariates = c(x, "hypertension"), rows = all_rows, ci = FALSE),
  S8 = list(label = "Measured BMI instead of estimated pre-pregnancy BMI", outcome = "gdm",
            covariates = sub("^bmi$", "bmi_measured", x), rows = all_rows, ci = FALSE),
  S9 = list(label = "Excluding women not tested on a fasted-appointment day", outcome = "gdm",
            covariates = x, rows = sample$tested_on_fasted_day, ci = TRUE),
  S11 = list(label = "Excluding women who drank other than water within 8 hours", outcome = "gdm",
             covariates = x, rows = !sample$drank_8h, ci = TRUE),
  S12 = list(label = "National wealth quintiles (poorest two, middle, richest two)", outcome = "gdm",
             covariates = sub("^wealth$", "wealth_v190", x), rows = all_rows, ci = FALSE))

# Outcome variable of a scenario on a given scale
on_scale <- function(outcome, s) if (s == "unconverted") paste0(outcome, "_alt") else outcome

prev <- list(); eff <- list(); ineq <- list()
for (s in names(scales)) for (id in names(scenarios)) {
  sc <- scenarios[[id]]
  o <- on_scale(sc$outcome, s)
  key <- paste(s, id)
  d_s <- subset(des, sc$rows)
  prev[[key]] <- cbind(scale = s, id = id, scenario = sc$label, weighted_prevalence(d_s, o))
  eff[[key]] <- cbind(scale = s, id = id, scenario = sc$label,
                      effect_table(fit_pr(d_s, sc$covariates, o), "PR", "adjusted"))
  if (sc$ci) {
    keep <- sc$rows & !is.na(sample[[o]])
    d <- sample[keep, ]
    est <- 4 * weighted.mean(d[[o]], d$weight) * concentration_index(d[[o]], d$wealth_score, d$weight)
    reps <- bootstrap_erreygers(d[[o]], d$wealth_score, boot$weights[keep, ])
    ineq[[key]] <- cbind(scale = s, id = id, scenario = sc$label, bootstrap_summary(est, reps, boot))
  }
}

# S10: unweighted prevalence with an exact binomial CI
for (s in names(scales)) {
  y <- sample[[scales[[s]]]]
  y <- y[!is.na(y)]
  bt <- binom.test(sum(y), length(y))
  prev[[paste(s, "S10")]] <- data.frame(scale = s, id = "S10", scenario = "Unweighted (exact binomial CI)",
                                        n = length(y), events = sum(y), prevalence = mean(y), se = NA,
                                        ci_lower = bt$conf.int[1], ci_upper = bt$conf.int[2])
}

save_table(do.call(rbind, prev), "sensitivity_prevalence")
eff <- do.call(rbind, eff)
eff$term <- sub("^wealth_v190", "wealth", eff$term)   # S12: same rows of the table as the tertiles
save_table(eff, "sensitivity_effects")
save_table(do.call(rbind, ineq), "sensitivity_inequality")
print(do.call(rbind, prev))

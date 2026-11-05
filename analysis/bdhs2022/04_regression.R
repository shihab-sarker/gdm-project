# 04_regression.R
# Prevalence ratios (main effect measure) and odds ratios on both glucose
# scales:
#   crude     one covariate at a time
#   minimal   the exposure (wealth or residence) with age group and trimester
#             only, the confounders in the conceptual framework (total
#             association; residence and BMI may lie on the pathway)
#   adjusted  all covariates of the anchor study (Mazumder et al. 2022)
# Also: the anchor-study specification (weeks instead of trimester), a model
# additionally adjusted for hypertension, and Firth logistic regression.

source(here::here("src", "common", "setup.R"))
women  <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))
des <- svy_domain(sample, women)
x <- cfg$covariates
confounders <- c("age_group", "trimester")

minimal_effects <- function(des, outcome) {
  do.call(rbind, lapply(c("wealth", "residence"), function(v) {
    e <- effect_table(fit_pr(des, c(v, confounders), outcome), "PR", "minimal")
    e[startsWith(e$term, v), ]
  }))
}

effects <- list(); anchor <- list(); htn <- list(); firth <- list()
for (s in names(scales)) {
  o <- scales[[s]]
  model_pr <- fit_pr(des, x, o)
  model_or <- fit_or(des, x, o)
  effects[[s]] <- cbind(scale = s, rbind(
    crude_effects(des, x, "PR", o), minimal_effects(des, o), effect_table(model_pr, "PR", "adjusted"),
    crude_effects(des, x, "OR", o), effect_table(model_or, "OR", "adjusted")))
  anchor[[s]] <- cbind(scale = s, rbind(
    effect_table(fit_pr(des, cfg$covariates_anchor, o), "PR", "adjusted"),
    effect_table(fit_or(des, cfg$covariates_anchor, o), "OR", "adjusted")))
  x_htn <- c(x, "hypertension")
  htn[[s]] <- cbind(scale = s, rbind(
    effect_table(fit_pr(des, x_htn, o), "PR", "adjusted + hypertension"),
    effect_table(fit_or(des, x_htn, o), "OR", "adjusted + hypertension")))
  firth[[s]] <- cbind(scale = s, rbind(cbind(model = "main", firth_check(sample, x, o)),
                                       cbind(model = "+ hypertension", firth_check(sample, x_htn, o))))
  saveRDS(model_pr, file.path(dir_models, paste0("model_pr_", s, ".rds")))
  saveRDS(model_or, file.path(dir_models, paste0("model_or_", s, ".rds")))
}
effects <- do.call(rbind, effects)
print(effects[effects$model == "adjusted" & effects$measure == "PR", ])

# The released-scale model sample must equal the last row of the participant flow
stopifnot(unique(effects$n[effects$scale == "released" & effects$model == "adjusted"]) ==
            tail(read_table("sample_flow")$n_remaining, 1))

save_table(effects, "effects_main")
save_table(do.call(rbind, anchor), "effects_anchor")
save_table(do.call(rbind, htn), "effects_hypertension")
save_table(do.call(rbind, firth), "firth")

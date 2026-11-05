# 05_model_diagnostics.R
# Diagnostics of the adjusted models on both glucose scales: events per
# parameter, collinearity (GVIF and survey-weighted VIF), Archer-Lemeshow
# goodness of fit, bootstrap internal validation (discrimination and
# calibration), linearity of the continuous covariates, and influential
# observations.

source(here::here("src", "common", "setup.R"))
women  <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))

by_scale <- function(f) {
  do.call(rbind, lapply(names(scales), function(s) cbind(scale = s, f(scales[[s]], s))))
}
model_or <- function(s) readRDS(file.path(dir_models, paste0("model_or_", s, ".rds")))

save_table(by_scale(function(o, s) events_per_parameter(model_or(s))), "epv")

# Collinearity concerns the covariates only, which are the same on both scales
vif <- collinearity(model_or("released"), sample, cfg)
print(vif)
save_table(vif$gvif, "vif_generalised")
save_table(vif$survey_vif, "vif_survey")

save_table(by_scale(function(o, s) archer_lemeshow(sample, cfg, outcome = o)), "archer_lemeshow")
save_table(by_scale(function(o, s) bootstrap_validation(sample, cfg, o)), "bootstrap_validation")
save_table(by_scale(function(o, s) calibration_groups(sample, cfg, outcome = o)), "calibration_groups")
save_table(by_scale(function(o, s) linearity_tests(sample, women, cfg, o)), "linearity_tests")
save_table(by_scale(function(o, s) influence_checks(sample, women, cfg, o)), "influence_refits")

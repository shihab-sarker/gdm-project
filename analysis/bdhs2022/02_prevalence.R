# 02_prevalence.R
# Survey-weighted prevalence of GDM and of the other glucose outcomes on both
# glucose scales, prevalence by trimester, and checks of the survey design.

source(here::here("src", "common", "setup.R"))
women  <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))
des <- svy_domain(sample, women)

outcomes <- c(gdm              = "GDM, FPG 5.1-6.9 mmol/L, released scale",
              fpg_dm_range     = "FPG 7.0 mmol/L or above, released scale",
              gdm_incl_dm      = "FPG 5.1 mmol/L or above, released scale",
              gdm_alt          = "GDM, FPG 5.1-6.9 mmol/L, unconverted scale",
              fpg_dm_range_alt = "FPG 7.0 mmol/L or above, unconverted scale",
              gdm_incl_dm_alt  = "FPG 5.1 mmol/L or above, unconverted scale")
prevalence <- do.call(rbind, lapply(names(outcomes), function(v) {
  cbind(outcome = v, label = outcomes[[v]], weighted_prevalence(des, v))
}))
print(prevalence)
save_table(prevalence, "prevalence")
save_table(do.call(rbind, lapply(names(scales), function(s) {
  cbind(scale = s, prevalence_by(des, scales[[s]], "trimester"))
})), "prevalence_by_trimester")

save_table(do.call(rbind, lapply(names(scales), function(s) {
  cbind(scale = s, design_effect(des, sample, scales[[s]]))
})), "design_effect")
save_table(psu_structure(sample, women), "psu_structure")
save_table(design_comparison(sample, women, cfg), "design_comparison")

# 11_bdhs2017.R
# BDHS 2017-18 with the definitions and design-based methods of the 2022
# analysis, for the comparison in Discussion 4.1. Settings: "bdhs2017_18" in
# config/config.yml. Results go to output/bdhs2017_18/ and data/derived/bdhs2017_18/, so no
# 2022 output is touched. Skipped if the 2017-18 files are not available.

if (Sys.getenv("R_CONFIG_ACTIVE") %in% c("", "bdhs2022") && file.exists(here::here("data", "raw", "bdhs2017_18", "BDIR7RFL.DTA"))) {
  Sys.setenv(R_CONFIG_ACTIVE = "bdhs2017_18")
  source(here::here("src", "common", "setup.R"))
  Sys.unsetenv("R_CONFIG_ACTIVE")

  files <- read_and_merge_2017(cfg)
  women <- derive_variables(files$women, cfg)
  sample <- women[women$analytic, ]
  saveRDS(women, data_file("women"))
  saveRDS(sample, data_file("analytic_sample"))
  flow <- participant_flow(women, cfg)
  print(flow)
  save_table(flow, "sample_flow")

  # Prevalence
  des <- svy_domain(sample, women)
  prev <- rbind(cbind(group = "all", level = "all", weighted_prevalence(des, "gdm")),
                do.call(rbind, lapply(c("trimester", "residence", "wealth"), function(v) prevalence_by(des, "gdm", v))))
  print(prev)
  save_table(prev, "prevalence")
  save_table(cbind(outcome = "fpg_dm_range", weighted_prevalence(des, "fpg_dm_range")), "diabetes_range_pregnant")

  # Prevalence ratios
  x <- cfg$covariates
  adjusted <- fit_pr(des, x, "gdm")
  minimal <- do.call(rbind, lapply(c("wealth", "residence"), function(v) {
    e <- effect_table(fit_pr(des, c(v, "age_group", "trimester"), "gdm"), "PR", "minimal")
    e[startsWith(e$term, v), ]
  }))
  effects <- rbind(effect_table(adjusted, "PR", "adjusted"), minimal)
  effects$design_df <- degf(adjusted$survey.design)
  print(effects)
  save_table(effects, "effects")

  # Inequality
  boot <- bootstrap_weights(sample, women, cfg)
  keep <- !is.na(sample$gdm)
  d <- sample[keep, ]
  est <- inequality_indices(d$gdm, d$wealth_score, d$weight)
  reps <- bootstrap_indices(d$gdm, d$wealth_score, boot$weights[keep, ])
  ineq <- cbind(index = c("relative", "wagstaff", "erreygers"), rbind(
    bootstrap_summary(est$C_relative, reps[1, ], boot),
    bootstrap_summary(est$C_wagstaff, reps[2, ], boot),
    bootstrap_summary(est$C_erreygers, reps[3, ], boot)))
  print(ineq)
  save_table(ineq, "inequality")

  # Glucose in pregnant and non-pregnant women, and the report's diabetes figure
  dist <- fpg_distribution(women, "fpg")
  print(dist)
  save_table(dist, "fpg_distribution")
  dm <- diabetes_2017(files$pr)
  print(dm)
  save_table(dm, "diabetes_report_replication")

  # Rounding signature, with the 2022 variables as controls of known factor
  p22 <- zap_labels(read_dta(here::here("data", "raw", "bdhs2022", "BDPR81FL.DTA"), col_select = dplyr::all_of(c("sb267", "sb267g"))))
  ok22 <- !p22$sb267 %in% c(NA, 994:996) & !p22$sb267g %in% c(NA, 994:996)
  tenths17 <- round(files$pr$sb267g[!is.na(files$pr$sb267g)] * 10)
  sig <- rbind(rounding_signature(round(p22$sb267g[ok22] * 10), "2022 sb267g (known factor 1.11)"),
               rounding_signature(round(p22$sb267[ok22] / 18 * 10), "2022 sb267/18 (known factor 1.00)"),
               rounding_signature(tenths17, "2017-18 sb335b, all adults tested"))
  print(sig)
  save_table(sig, "rounding_signature")
}

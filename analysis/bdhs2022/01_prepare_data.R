# 01_prepare_data.R
# Read the BDHS 2022 women's and household-member files, link each woman to
# her biomarker record, derive the analysis variables and define the
# analytic sample. Output: data/derived/bdhs2022/women.rds, data/derived/bdhs2022/analytic_sample.rds and the
# sample-accounting tables.

source(here::here("src", "common", "setup.R"))

files <- read_and_merge(cfg)
women <- derive_variables(files$women, cfg)

flow <- participant_flow(women, cfg)
print(flow)
save_table(flow, "sample_flow")
save_table(fasting_check(women), "fasting_check")
save_table(diabetes_history_counts(women), "diabetes_history")
save_table(glucose_scale_relation(files$pr), "glucose_scale_relation")
save_table(weight_variables(cfg), "weight_variables")
save_table(dhs_report_replication(files$pr), "dhs_report_replication")

keep <- c("v001", "v002", "v003", "v021", "v022", "weight", "v190", "v213", "v214",
          "shbpbg", "sb236", "fasted", "fpg", "fpg_alt", "pre_dm",
          "gdm", "gdm_alt", "gdm_incl_dm", "fpg_dm_range", "gdm_incl_dm_alt", "fpg_dm_range_alt",
          "tested_on_fasted_day", "hours_ate", "hours_drank", "drank_never", "drank_8h", "gest_month", "age_c", "bmi_measured", "v025",
          "age", "age_group", "age_group_alt", "education", "working", "residence",
          "birth_order", "gest_weeks", "trimester", "bmi", "hypertension",
          "wealth_score", "wealth", "wealth_v190", "analytic")
women <- women[keep]
saveRDS(women, data_file("women"))
saveRDS(women[women$analytic, ], data_file("analytic_sample"))
cat("Analytic sample:", sum(women$analytic), "women\n")

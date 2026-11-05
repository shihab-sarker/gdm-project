# 09_tables_figures.R
# Manuscript tables (formatted CSVs) and figures, from the saved results.
# Both glucose scales appear side by side (tables) or overlaid (figures).
# Table 1: characteristics by household wealth tertile and GDM on both scales.

source(here::here("src", "common", "setup.R"))
women <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))
fig <- function(name) file.path(dir_figures, paste0(name, ".png"))

effects <- read_table("effects_main")
pr <- effects[effects$measure == "PR", ]
ci <- read_table("inequality_bootstrap")

save_table(table1_wealth(read_table("table1_by_wealth"), read_table("table1_categorical"), read_table("table1_continuous")),
           "manuscript_table1")
save_table(effects_table(pr, with_p = FALSE), "manuscript_table2")
save_table(effects_table(effects[effects$measure == "OR", ]), "manuscript_table_or")
save_table(inequality_table(ci, read_table("decomposition_main")), "manuscript_table3")
save_table(inequality_table(ci, read_table("decomposition_main"), full = TRUE), "manuscript_table_decomp_full")
save_table(inequality_table(ci, read_table("decomposition_with_wealth"), full = TRUE),
           "manuscript_table_decomp_wealth")
save_table(sensitivity_table(read_table("sensitivity_prevalence"), read_table("sensitivity_effects"),
                             read_table("sensitivity_inequality")), "manuscript_table_sensitivity")
save_table(effects_table(read_table("effects_anchor")), "manuscript_table_anchor")
save_table(effects_table(read_table("effects_hypertension")), "manuscript_table_hypertension")

plot_flow(read_table("sample_flow"), fig("Figure_1"), n_alt = sum(!is.na(sample$gdm_alt)))
plot_inequality(read_table("concentration_curve"), read_table("decomposition_main"), fig("Figure_2"))
plot_calibration(read_table("calibration_groups"), fig("Figure_SA1"))
plot_fpg_ecdf(women, fig("Figure_SA3"))

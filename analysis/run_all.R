# analysis/run_all.R
# Runs the complete analysis of both surveys, from the raw files to the manuscript
# tables and figures. Each script reads its inputs from data/ or output/ and
# writes its own results, so any script can also be run on its own.
#
# From the project root:  Rscript analysis/run_all.R
# To keep a log:          Rscript analysis/run_all.R > output/bdhs2022/run_log.txt 2>&1

rm(list = ls())

source(here::here("analysis", "bdhs2022", "01_prepare_data.R"))
source(here::here("analysis", "bdhs2022", "02_prevalence.R"))
source(here::here("analysis", "bdhs2022", "03_characteristics.R"))
source(here::here("analysis", "bdhs2022", "04_regression.R"))
source(here::here("analysis", "bdhs2022", "05_model_diagnostics.R"))
source(here::here("analysis", "bdhs2022", "06_inequality.R"))
source(here::here("analysis", "bdhs2022", "07_sensitivity.R"))
source(here::here("analysis", "bdhs2022", "08_validity_checks.R"))
source(here::here("analysis", "bdhs2022", "09_tables_figures.R"))
source(here::here("analysis", "bdhs2022", "10_reviewer_checks.R"))
source(here::here("analysis", "bdhs2017_18", "11_bdhs2017.R"))

warnings()
sessionInfo()

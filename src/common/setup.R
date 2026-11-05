# src/common/setup.R
# Sourced at the top of every analysis script. Loads the packages, reads the
# settings in config/config.yml, defines where files are written, and sources
# the analysis functions in src/common, src/bdhs2022 and src/bdhs2017_18.

suppressPackageStartupMessages({
  library(haven)    # reading DHS Stata files
  library(survey)   # design-based estimation (Lumley 2004)
})

# One settings block per survey; without R_CONFIG_ACTIVE the 2022 block is used.
cfg_name <- Sys.getenv("R_CONFIG_ACTIVE")
if (cfg_name %in% c("", "default")) cfg_name <- "bdhs2022"
cfg <- config::get(config = cfg_name, file = here::here("config", "config.yml"))

# Strata containing a single PSU in the analytic domain: centre that PSU at
# the grand mean instead of failing (survey package option).
options(survey.lonely.psu = "adjust")

# Output folders
dir_tables  <- here::here(cfg$dir_output, "tables")
dir_figures <- here::here(cfg$dir_output, "figures")
dir_models  <- here::here(cfg$dir_output, "models")
for (d in c(here::here(cfg$dir_data), dir_tables, dir_figures, dir_models)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Every result table is a CSV in <dir_output>/tables, so scripts communicate only
# through files on disk.
save_table <- function(x, name) {
  write.csv(x, file.path(dir_tables, paste0(name, ".csv")), row.names = FALSE)
  invisible(x)
}
read_table <- function(name) {
  read.csv(file.path(dir_tables, paste0(name, ".csv")), check.names = FALSE)
}
data_file <- function(name) here::here(cfg$dir_data, paste0(name, ".rds"))

for (f in c("data_preparation.R", "survey_design.R", "regression.R", "inequality.R", "validity_checks.R")) {
  source(here::here("src", "common", f))
}
source(here::here("src", "bdhs2017_18", "bdhs2017.R"))
source(here::here("src", "bdhs2022", "tables_figures.R"))

# The two glucose scales analysed in parallel throughout:
#   released    : sb267g, the DHS-released mmol/L value (sb267 x 1.11 / 18)
#   unconverted : sb267 / 18, without the whole-blood-to-plasma factor 1.11
# Outcome variables for GDM (FPG 5.1-6.9 mmol/L) on each scale:
scales <- c(released = "gdm", unconverted = "gdm_alt")

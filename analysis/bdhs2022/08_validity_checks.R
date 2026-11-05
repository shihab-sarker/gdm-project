# 08_validity_checks.R
# Checks on the glucose outcome (distribution in pregnant and non-pregnant
# women on both scales, digit preference, prevalence by month of pregnancy),
# timing of the test by socioeconomic position, and selection into the
# analytic sample.

source(here::here("src", "common", "setup.R"))
women  <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))

fpg <- rbind(fpg_distribution(women, "fpg"), fpg_distribution(women, "fpg_alt"))
print(fpg)
save_table(fpg, "fpg_distribution")
save_table(glucose_heaping(women), "glucose_heaping")
save_table(do.call(rbind, lapply(names(scales), function(s) {
  cbind(scale = s, prevalence_by_month(sample, women, scales[[s]]))
})), "prevalence_by_month")
save_table(timing_by_ses(sample, women), "timing_by_ses")
save_table(selection_comparison(women, cfg), "selection_comparison")
paired <- paired_classification(sample, women, cfg)
print(paired)
save_table(paired, "paired_classification")
common <- common_sample(sample, women, cfg)
print(common)
save_table(common, "common_sample")

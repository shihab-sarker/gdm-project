# 06_inequality.R
# Wealth-related inequality in GDM on both glucose scales: relative, Wagstaff
# and Erreygers concentration indices with Rao-Wu bootstrap 95% CIs,
# concentration curves, and decomposition of the Erreygers index (main:
# wealth tertiles excluded; supplementary: included).

source(here::here("src", "common", "setup.R"))
women  <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))
boot <- bootstrap_weights(sample, women, cfg)

indices <- list(); curve <- list(); ci <- list(); dec <- list(); dec_w <- list()
for (s in names(scales)) {
  o <- scales[[s]]
  keep <- !is.na(sample[[o]])
  d <- sample[keep, ]
  est <- inequality_indices(d[[o]], d$wealth_score, d$weight)
  indices[[s]] <- cbind(scale = s, est)
  curve[[s]] <- cbind(scale = s, concentration_curve(d[[o]], d$wealth_score, d$weight))
  reps <- bootstrap_indices(d[[o]], d$wealth_score, boot$weights[keep, ])
  ci[[s]] <- cbind(scale = s, index = c("relative", "wagstaff", "erreygers"), rbind(
    bootstrap_summary(est$C_relative, reps[1, ], boot),
    bootstrap_summary(est$C_wagstaff, reps[2, ], boot),
    bootstrap_summary(est$C_erreygers, reps[3, ], boot)))
  dec[[s]] <- cbind(scale = s, decomposition_with_ci(sample, cfg$decomp_covariates, boot, o))
  dec_w[[s]] <- cbind(scale = s, decomposition_with_ci(sample, c(cfg$decomp_covariates, "wealth"), boot, o))
}
ci <- do.call(rbind, ci)
print(do.call(rbind, indices)); print(ci)
save_table(do.call(rbind, indices), "inequality_indices")
save_table(do.call(rbind, curve), "concentration_curve")
save_table(ci, "inequality_bootstrap")
save_table(ci[ci$index == "erreygers", ], "erreygers_bootstrap")
save_table(do.call(rbind, dec), "decomposition_main")
save_table(do.call(rbind, dec_w), "decomposition_with_wealth")

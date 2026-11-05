# tests/make_synthetic_data.R
# Creates small synthetic IR and PR Stata files with the same variable names
# and value codes as BDHS 2022, so the pipeline and unit tests can run without
# the restricted DHS microdata. Values are random; no real record is used.
#
# Usage (from the project root): Rscript tests/make_synthetic_data.R
# Output: tests/synthetic/BDIR_SYN.DTA, tests/synthetic/BDPR_SYN.DTA

make_synthetic_bdhs <- function(out_dir = here::here("tests", "synthetic"),
                                n_strata = 20, clusters_per_stratum = 4,
                                women_per_cluster = 15, seed = 2022) {
  set.seed(seed)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  n_clusters <- n_strata * clusters_per_stratum
  cl <- data.frame(
    v001 = seq_len(n_clusters),
    v022 = rep(seq_len(n_strata), each = clusters_per_stratum),
    v025 = sample(1:2, n_clusters, replace = TRUE, prob = c(0.35, 0.65)),
    v024 = sample(1:8, n_clusters, replace = TRUE),
    v005 = round(stats::runif(n_clusters, 0.4e6, 2.0e6)),
    shbpbg = sample(c(0, 1), n_clusters, replace = TRUE, prob = c(0.3, 0.7))
  )
  ir <- cl[rep(seq_len(n_clusters), each = women_per_cluster), ]
  n <- nrow(ir)
  ir$v002 <- rep(seq_len(women_per_cluster), n_clusters) # one woman per household
  ir$v003 <- 2L
  ir$v021 <- ir$v001

  urban <- ir$v025 == 1
  ir$v191 <- round(stats::rnorm(n, ifelse(urban, 60000, -20000), 80000))
  ir$v190 <- as.integer(cut(ir$v191, stats::quantile(ir$v191, 0:5 / 5),
                            include.lowest = TRUE, labels = FALSE))
  ir$v012 <- sample(15:49, n, replace = TRUE)
  ir$v013 <- pmin(7L, (ir$v012 - 15L) %/% 5L + 1L)
  ir$v106 <- sample(0:3, n, replace = TRUE, prob = c(0.15, 0.3, 0.4, 0.15))
  ir$v714 <- sample(0:1, n, replace = TRUE, prob = c(0.7, 0.3))
  ir$v220 <- sample(0:6, n, replace = TRUE)
  # Oversample pregnancy so the analytic sample is large enough for models
  ir$v213 <- sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5))
  ir$v214 <- ifelse(ir$v213 == 1, sample(1:9, n, replace = TRUE), NA)
  ir$v438 <- round(stats::rnorm(n, 1510, 55))
  ir$v437 <- round(stats::rnorm(n, 520, 90))
  ir$v437[sample(n, 5)] <- 9999 # missing code
  ir$v438[sample(n, 5)] <- 9999

  # PR file: the women above plus male household members (hv104 = 1)
  wealth_z <- as.numeric(scale(ir$v191))
  mgdl <- round(stats::rnorm(n, 82 + 4 * wealth_z + 3 * urban - 0.8 * ifelse(is.na(ir$v214), 0, ir$v214), 12))
  biom <- ir$shbpbg == 1 & ir$v012 >= 18
  sb267 <- ifelse(biom, mgdl, NA)
  sb267[which(biom)[1:4]] <- 995 # missing code
  sb267g <- ifelse(sb267 %in% 994:996, sb267, round(sb267 * 1.11 / 18, 1))

  fast1 <- ifelse(biom, sample(0:1, n, replace = TRUE, prob = c(0.05, 0.95)), NA)
  fast2 <- ifelse(biom & fast1 == 0, sample(0:1, n, replace = TRUE, prob = c(0.5, 0.5)), NA)
  day1 <- ifelse(biom, sample(1:28, n, replace = TRUE), NA)
  day2 <- ifelse(!is.na(fast2), pmin(day1 + 1, 31), NA)
  test_day <- ifelse(fast1 %in% 1, day1, ifelse(fast2 %in% 1, day2, NA))

  sbp <- ifelse(biom, round(stats::rnorm(n, 112, 14)), NA)
  dbp <- ifelse(biom, round(stats::rnorm(n, 74, 10)), NA)

  women_pr <- data.frame(
    hv001 = ir$v001, hv002 = ir$v002, hvidx = ir$v003, hv104 = 2L,
    hv005 = ir$v005, hv105 = ir$v012, hv025 = ir$v025,
    shbpbg = ir$shbpbg,
    sb267 = sb267, sb267g = sb267g,
    sb261a = 2100, sb261b = 2100, sb261c = NA_real_, sb262a = 9595, sb262b = 9595, sb262c = NA_real_,
    sb263a = fast1, sb263b = fast2, sb263c = NA_real_, sb266t = 700,
    sb236 = ifelse(biom, sample(c(0, 1, 8), n, replace = TRUE, prob = c(0.95, 0.04, 0.01)), NA),
    sb240 = NA_real_,
    sb258d1 = day1, sb258d2 = day2, sb258d3 = NA_real_, sb266d = test_day,
    wbp13 = sbp, wbp14 = dbp,
    wbp15 = ifelse(biom, sample(0:1, n, replace = TRUE), NA),
    wbp19 = ifelse(biom, sample(0:1, n, replace = TRUE, prob = c(0.97, 0.03)), NA),
    wbp24 = sbp + ifelse(biom, sample(-3:3, n, replace = TRUE), NA),
    wbp25 = dbp + ifelse(biom, sample(-3:3, n, replace = TRUE), NA)
  )
  men_pr <- women_pr[seq_len(200), ]
  men_pr[, setdiff(names(men_pr), c("hv001", "hv002", "shbpbg"))] <- NA
  men_pr$hvidx <- 1L
  men_pr$hv104 <- 1L
  pr <- rbind(women_pr, men_pr)

  ir$shbpbg <- NULL # biomarker selection is a PR (household) variable only
  haven::write_dta(ir, file.path(out_dir, "BDIR_SYN.DTA"))
  haven::write_dta(pr, file.path(out_dir, "BDPR_SYN.DTA"))
  invisible(list(ir = ir, pr = pr))
}

if (sys.nframe() == 0L) {
  d <- make_synthetic_bdhs()
  cat("Synthetic IR rows:", nrow(d$ir), " PR rows:", nrow(d$pr), "\n")
}

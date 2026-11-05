ci_data <- function(n = 400) {
  set.seed(1)
  d <- data.frame(wealth_score = rnorm(n), weight = runif(n, 0.5, 2))
  d$residence <- factor(ifelse(d$wealth_score + rnorm(n) > 0.5, "Urban", "Rural"), c("Rural", "Urban"))
  d$age_group <- factor(sample(c("<25", "25+"), n, TRUE), c("<25", "25+"))
  d$gdm <- rbinom(n, 1, plogis(0.3 + 0.6 * d$wealth_score))
  d
}

test_that("with equal weights the index equals the textbook formula 2/(n mu) sum(y R) - 1", {
  d <- ci_data()
  n <- nrow(d); R <- (rank(d$wealth_score) - 0.5) / n
  expect_equal(concentration_index(d$gdm, d$wealth_score, rep(1, n)),
               2 / (n * mean(d$gdm)) * sum(d$gdm * R) - 1)
})

test_that("Erreygers and Wagstaff indices are rescaled relative indices", {
  d <- ci_data()
  r <- inequality_indices(d$gdm, d$wealth_score, d$weight)
  mu <- weighted.mean(d$gdm, d$weight)
  expect_equal(r$C_erreygers, 4 * mu * r$C_relative)
  expect_equal(r$C_wagstaff, r$C_relative / (1 - mu))
  expect_gt(r$C_erreygers, 0)   # outcome generated to rise with wealth
})

test_that("decomposition contributions and residual add up to the Erreygers index", {
  d <- ci_data()
  dec <- decompose_index(d, c("age_group", "residence"))
  expect_equal(sum(dec$contribution), inequality_indices(d$gdm, d$wealth_score, d$weight)$C_erreygers)
  expect_equal(sum(dec$percent), 100)
})

test_that("modified Poisson regression recovers a known prevalence ratio", {
  set.seed(2)
  d <- data.frame(x = factor(sample(c("a", "b"), 4000, TRUE)), v021 = rep(1:200, each = 20),
                  v022 = rep(1:20, each = 200), weight = 1)
  d$gdm <- rbinom(4000, 1, ifelse(d$x == "b", 0.6, 0.4))
  e <- effect_table(fit_pr(svy_restricted(d), "x"), "PR", "test")
  expect_equal(e$estimate, mean(d$gdm[d$x == "b"]) / mean(d$gdm[d$x == "a"]), tolerance = 1e-6)
  expect_true(e$ci_lower < 1.5 && e$ci_upper > 1.5)
})

test_that("Rao-Wu replicate weights equal those of survey::as.svrepdesign exactly", {
  set.seed(3)
  women <- data.frame(v021 = rep(1:12, each = 5), v022 = rep(1:3, each = 20), weight = runif(60, 0.5, 2))
  sample <- women[women$v021 != 4, ]   # PSU 4 has no woman in the domain
  r <- bootstrap_weights(sample, women, list(seed = 1, bootstrap_B = 20))
  frame <- domain_frame(sample, women)
  set.seed(1)
  ref <- as.svrepdesign(svydesign(ids = ~v021, strata = ~v022, weights = ~weight, data = frame),
                        type = "subbootstrap", replicates = 20)
  expect_identical(unname(r$weights), unname(weights(ref, "analysis")[frame$in_domain, ]))
  expect_equal(r$scale, ref$scale)
  expect_equal(r$rscales, ref$rscales)
})

test_that("calibration intercept and slope equal those of glm()", {
  set.seed(4)
  p <- runif(300, 0.05, 0.6)
  y <- rbinom(300, 1, p)
  lp <- qlogis(p)
  r <- calibration_stats(y, p)
  expect_equal(r[["intercept"]], unname(coef(glm(y ~ offset(lp), family = binomial))[1]))
  expect_equal(r[["slope"]], unname(coef(glm(y ~ lp, family = binomial))[2]))
})

test_that("rank-based AUC equals the proportion of correctly ordered pairs", {
  y <- c(0, 0, 1, 1, 0, 1); p <- c(0.1, 0.4, 0.35, 0.8, 0.2, 0.9)
  pairs <- outer(p[y == 1], p[y == 0], ">")
  expect_equal(auc(y, p), mean(pairs))
})

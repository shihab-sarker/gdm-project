# src/common/inequality.R
# Wealth-related inequality in GDM: concentration index, concentration
# curve and regression-based decomposition (O'Donnell et al. 2008).
#
# Ranking variable: household wealth index score (v191). With fractional
# rank R (weighted), the relative concentration index of y is
#     C = 2 cov_w(y, R) / mean_w(y).
# For a binary outcome the bounds of C depend on the prevalence (Wagstaff
# 2005), so the Erreygers (2009) corrected index is the main measure:
#     E = 4 mean_w(y) C   (binary y, so max - min = 1).
# The Wagstaff (2005) normalised index C / (1 - mean_w(y)) is also reported.

# Weighted fractional rank: cumulative weight share up to the midpoint of
# each woman's own weight, after sorting by the wealth score.
fractional_rank <- function(x, w) {
  o <- order(x)
  w <- w / sum(w)
  r <- numeric(length(x))
  r[o] <- cumsum(w[o]) - w[o] / 2
  r
}

weighted_cov <- function(a, b, w) {
  w <- w / sum(w)
  sum(w * (a - sum(w * a)) * (b - sum(w * b)))
}

# Relative concentration index of y by wealth score x, with weights w
concentration_index <- function(y, x, w) {
  concentration_given_rank(y, fractional_rank(x, w), w)
}

# The same index from ranks r already computed, for several variables at one rank
concentration_given_rank <- function(y, r, w) {
  2 * weighted_cov(y, r, w) / weighted.mean(y, w)
}

# All three indices for a binary outcome
inequality_indices <- function(y, x, w) {
  mu <- weighted.mean(y, w)
  C <- concentration_index(y, x, w)
  data.frame(n = length(y), prevalence = mu, C_relative = C,
             C_wagstaff = C / (1 - mu), C_erreygers = 4 * mu * C,
             bound_erreygers = 4 * mu * (1 - mu))
}

# Concentration curve: cumulative share of women ranked by wealth (x axis)
# against cumulative share of GDM cases (y axis), both weighted.
concentration_curve <- function(y, x, w) {
  o <- order(x)
  data.frame(cum_women = c(0, cumsum(w[o]) / sum(w)),
             cum_gdm = c(0, cumsum(w[o] * y[o]) / sum(w * y)))
}

# Erreygers index of every bootstrap replicate (replicate weights in columns)
bootstrap_erreygers <- function(y, x, W) {
  apply(W, 2, function(w) {
    keep <- w > 0
    4 * weighted.mean(y[keep], w[keep]) * concentration_index(y[keep], x[keep], w[keep])
  })
}

# Bootstrap 95% CIs for all three indices (relative, Wagstaff, Erreygers)
# (rows: relative, Wagstaff, Erreygers; one column per replicate)
bootstrap_indices <- function(y, x, W) {
  reps <- apply(W, 2, function(w) {
    keep <- w > 0
    r <- inequality_indices(y[keep], x[keep], w[keep])
    c(r$C_relative, r$C_wagstaff, r$C_erreygers)
  })
  reps
}

# Decomposition (Wagstaff, van Doorslaer & Watanabe 2003). In a weighted
# linear probability model y = a + sum_k b_k x_k + e, the concentration
# index splits into covariate contributions plus a residual:
#     C = sum_k (b_k mean(x_k) / mean(y)) C_k + C_e,
# where b_k mean(x_k) / mean(y) is the elasticity of y with respect to x_k
# and C_k the concentration index of x_k. Contributions are multiplied by
# 4 mean(y) so that they add up to the Erreygers index.
decompose_index <- function(d, covariates, w = d$weight, outcome = "gdm") {
  X <- decomposition_matrix(d, covariates)
  r <- decompose_terms(X, d[[outcome]], d$wealth_score, w)
  data.frame(term = c(colnames(X), "residual"),
             coefficient = c(r$b, NA), mean = c(r$mean_x, NA),
             elasticity = c(r$elasticity, NA), concentration_index = c(r$C_x, NA),
             contribution = r$contribution, percent = r$percent,
             row.names = NULL)
}

decomposition_matrix <- function(d, covariates) {
  model.matrix(reformulate(covariates), d)[, -1, drop = FALSE]
}

# Terms of the decomposition for covariate matrix X (no intercept column);
# contribution and percent end with the residual.
decompose_terms <- function(X, y, wealth_score, w) {
  b <- lm.wfit(cbind(1, X), y, w)$coefficients[-1]
  mu <- weighted.mean(y, w)
  r <- fractional_rank(wealth_score, w)
  mean_x <- apply(X, 2, weighted.mean, w = w)
  C_x <- apply(X, 2, concentration_given_rank, r = r, w = w)
  elasticity <- b * mean_x / mu
  contribution <- 4 * mu * elasticity * C_x
  total <- 4 * mu * concentration_given_rank(y, r, w)
  contribution <- c(contribution, total - sum(contribution))
  list(b = b, mean_x = mean_x, C_x = C_x, elasticity = elasticity,
       contribution = contribution, percent = 100 * contribution / total)
}

# Decomposition with Rao-Wu bootstrap 95% CIs for each contribution
decomposition_with_ci <- function(sample, covariates, boot, outcome = "gdm") {
  keep <- complete.cases(sample[c(outcome, covariates, "wealth_score")])
  d <- sample[keep, ]
  W <- boot$weights[keep, , drop = FALSE]
  est <- decompose_index(d, covariates, outcome = outcome)
  X <- decomposition_matrix(d, covariates)
  y <- d[[outcome]]
  reps <- apply(W, 2, function(w) {
    k <- w > 0
    r <- decompose_terms(X[k, , drop = FALSE], y[k], d$wealth_score[k], w[k])
    c(r$contribution, r$percent)
  })
  m <- nrow(est)
  est$contribution_ci_lower <- apply(reps[1:m, ], 1, quantile, 0.025, na.rm = TRUE)
  est$contribution_ci_upper <- apply(reps[1:m, ], 1, quantile, 0.975, na.rm = TRUE)
  est$percent_ci_lower <- apply(reps[m + 1:m, ], 1, quantile, 0.025, na.rm = TRUE)
  est$percent_ci_upper <- apply(reps[m + 1:m, ], 1, quantile, 0.975, na.rm = TRUE)
  cbind(n = nrow(d), est)
}

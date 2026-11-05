# src/bdhs2017_18/bdhs2017.R
# BDHS 2017-18 comparison. The recode names its biomarker variables differently,
# so they are renamed to their 2022 names on reading and the 2022 derivations are
# reused unchanged. The recode has one glucose variable, sb335b (mmol/L x 10, read
# as "released"); there is no mg/dL variable, so the two scales coincide. Blood
# pressure is not needed for the comparison and stays missing.

read_and_merge_2017 <- function(cfg) {
  ir <- zap_labels(read_dta(here::here(cfg$path_ir), col_select = dplyr::all_of(
    c("v001", "v002", "v003", "v005", "v021", "v022", "v012", "v025", "v106", "v190", "v191",
      "v213", "v214", "v220", "v437", "v438", "v714"))))
  pr <- zap_labels(read_dta(here::here(cfg$path_pr), col_select = dplyr::all_of(
    c("hv001", "hv002", "hvidx", "hv104", "hv005", "hv105", "hv025",
      "shelig", "sb326b", "sb327a", "sb334h", "sb335b"))))
  stopifnot(!anyDuplicated(ir[c("v001", "v002", "v003")]), !anyDuplicated(pr[c("hv001", "hv002", "hvidx")]))
  glucose <- replace(pr$sb335b, pr$sb335b %in% 994:996, NA) / 10
  pr <- dplyr::mutate(pr, shbpbg = shelig, sb236 = sb326b, sb240 = sb327a, sb263a = sb334h,
                      sb263b = NA_real_, sb263c = NA_real_, sb267g = glucose, sb267 = glucose * 18,
                      wbp19 = NA_real_, wbp24 = NA_real_, wbp25 = NA_real_)
  women <- dplyr::left_join(ir, pr, by = c(v001 = "hv001", v002 = "hv002", v003 = "hvidx"))
  stopifnot(nrow(women) == nrow(ir), all(women$hv104 %in% c(2, NA)))
  list(women = as.data.frame(women), pr = as.data.frame(pr))
}

# Diabetes among women aged 18+ (FPG >= 7.0 mmol/L or medication), household
# weight, the definition of the survey reports; the 2017-18 report gives 10%.
diabetes_2017 <- function(pr) {
  p <- pr[pr$hv104 %in% 2 & pr$hv105 >= 18 & pr$hv105 < 95 & pr$shbpbg %in% 1 & !is.na(pr$sb267g), ]
  data.frame(n = nrow(p), diabetes = 100 * weighted.mean(p$sb267g >= 7 | p$sb240 %in% 1, p$hv005 / 1e6),
             report_2017_18 = 10)
}

# Rounding signature of a glucose variable held in tenths of mmol/L. Integer
# mg/dL readings converted with a factor f reach some tenths from two integers
# and others from one, which leaves a repeating pattern in the counts. Counts per
# tenth are fitted by a Poisson model with a smooth density and, as offset, the
# number of integers reaching each tenth under f (or no offset: no pattern).
# The best f is the implied factor; lower deviance is a better fit.
rounding_signature <- function(tenths, label, lo = 35, hi = 90, grid = seq(0.90, 1.25, by = 0.005)) {
  t <- lo:hi
  count <- as.numeric(table(factor(tenths, t)))
  deviance_of <- function(offset) deviance(glm(count ~ splines::ns(t, 8), family = poisson, offset = offset))
  reach <- function(f) log(pmax(as.numeric(table(factor(round((1:400) * f / 18 * 10), t))), 0.5))
  fits <- vapply(grid, function(f) deviance_of(reach(f)), 0)
  none <- deviance_of(rep(0, length(t)))
  # no factor beats "no pattern" (f near 0.9 reaches every tenth equally, so it ties)
  best <- if (min(fits) < none - 1e-6) grid[which.min(fits)] else NA
  # x 1.11 applied to values already in tenths cannot reach some tenths
  unreachable <- setdiff(t, round((1:1000) * 1.11))
  data.frame(variable = label, n = sum(count), deviance_no_pattern = none,
             deviance_factor_1.00 = deviance_of(reach(1)), deviance_factor_1.11 = deviance_of(reach(1.11)),
             best_factor = best, df = length(t) - 9,
             n_unreachable_tenths = length(unreachable), n_at_unreachable_tenths = sum(count[t %in% unreachable]))
}

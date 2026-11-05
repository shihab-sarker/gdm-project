# Loaded by testthat before the tests. Uses the "test" settings, so any
# file written goes to tests/synthetic/out/ and never to the survey folders in output/.
library(testthat)
Sys.setenv(R_CONFIG_ACTIVE = "test")
source(here::here("src", "common", "setup.R"))

# One row per woman with the raw DHS columns used by derive_variables()
toy_women <- function(n = 6, ...) {
  d <- data.frame(v001 = 1, v002 = seq_len(n), v003 = 2, v005 = 1e6, v021 = 1, v022 = 1,
                  v012 = 25, v025 = 2, v106 = 2, v190 = 3, v191 = seq_len(n) * 1000,
                  v213 = 1, v214 = 5, v220 = 1, v437 = 600, v438 = 1500, v714 = 0,
                  shbpbg = 1, sb236 = 0, sb263a = 1, sb263b = NA, sb263c = NA,
                  sb267 = 90, sb267g = 5.55, wbp19 = 0, wbp24 = 110, wbp25 = 70)
  extra <- list(...)
  for (k in names(extra)) d[[k]] <- extra[[k]]
  d
}

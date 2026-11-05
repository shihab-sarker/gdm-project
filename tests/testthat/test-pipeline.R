# Runs analysis/run_all.R on synthetic data (settings "test", B = 50).
# Takes about a minute, so it runs only when RUN_PIPELINE_TEST=true.

test_that("the full analysis runs on synthetic data", {
  skip_if_not(Sys.getenv("RUN_PIPELINE_TEST") == "true", "set RUN_PIPELINE_TEST=true to run")
  source(here::here("tests", "make_synthetic_data.R"), local = TRUE)
  make_synthetic_bdhs()
  status <- withr::with_dir(here::here(), withr::with_envvar(c(R_CONFIG_ACTIVE = "test"), system2(
    file.path(R.home("bin"), "Rscript"), shQuote(file.path("analysis", "run_all.R")),
    stdout = FALSE, stderr = FALSE)))
  expect_equal(status, 0)
  flow <- read_table("sample_flow")
  eff <- read_table("effects_main")
  # both glucose scales are analysed; the released-scale model uses the final sample
  expect_setequal(unique(eff$scale), c("released", "unconverted"))
  expect_equal(unique(eff$n[eff$scale == "released" & eff$model == "adjusted"]), tail(flow$n_remaining, 1))
  expect_true(all(read_table("inequality_bootstrap")$failed == 0))
  expect_true(all(file.exists(file.path(dir_figures, paste0("Figure_", c("1", "2", "SA1", "SA2", "SA3"), ".pdf")))))
})

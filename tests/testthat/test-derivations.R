test_that("glucose: DHS missing codes become NA; alternative scale is sb267 / 18", {
  w <- derive_variables(toy_women(3, sb267 = c(90, 995, 120), sb267g = c(5.55, 995, 7.4)), cfg)
  expect_equal(w$fpg, c(5.55, NA, 7.4))
  expect_equal(w$fpg_alt, c(90, NA, 120) / 18)
})

test_that("outcomes follow the 5.1 and 7.0 mmol/L thresholds", {
  g <- c(5.0, 5.1, 6.9, 7.0, 6.0, 6.0)
  w <- derive_variables(toy_women(6, sb267g = g, sb236 = c(0, 0, 0, 0, 1, 8)), cfg)
  expect_equal(w$pre_dm, c(0, 0, 0, 0, 1, NA))
  expect_equal(w$gdm, c(0, 1, 1, NA, NA, NA))
  expect_equal(w$gdm_incl_dm, c(0, 1, 1, 1, NA, NA))
  expect_equal(w$fpg_dm_range, c(0, 0, 0, 1, NA, NA))
})

test_that("a missing diabetes answer counts as no pre-existing diabetes", {
  expect_equal(derive_variables(toy_women(3, sb236 = c(NA, 0, 1)), cfg)$pre_dm, c(0, 0, 1))
})

test_that("fasting at any appointment counts; no record fails", {
  w <- derive_variables(toy_women(4, sb263a = c(1, 0, 0, NA), sb263b = c(NA, 1, 0, NA)), cfg)
  expect_equal(w$fasted, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(w$analytic, c(TRUE, TRUE, FALSE, FALSE))
})

test_that("trimesters are months 1-3, 4-6 and 7+", {
  w <- derive_variables(toy_women(6, v214 = c(1, 3, 4, 6, 7, 9)), cfg)
  expect_equal(as.character(w$trimester), c("1st", "1st", "2nd", "2nd", "3rd", "3rd"))
})

test_that("pre-pregnancy BMI subtracts 3 kg plus 0.5 kg per week after week 20", {
  w <- derive_variables(toy_women(3, v214 = c(4, 7, 5), v437 = c(600, 600, 9999)), cfg)
  expect_equal(w$bmi[1], (60 - 3) / 1.5^2)
  expect_equal(w$bmi[2], (60 - 3 - 0.5 * (7.5 * 4.345 - 20)) / 1.5^2)
  expect_true(is.na(w$bmi[3]))
})

test_that("hypertension: final BP >= 140/90 or medication; missing BP is NA", {
  w <- derive_variables(toy_women(5, wbp24 = c(139, 140, 120, 995, 120),
                                  wbp25 = c(89, 80, 90, 95, 70), wbp19 = c(0, 0, 0, 0, 1)), cfg)
  expect_equal(as.character(w$hypertension), c("No", "Yes", "Yes", NA, "Yes"))
})

test_that("wealth tertiles are formed within the analytic sample only", {
  w <- derive_variables(toy_women(6, v213 = c(1, 1, 1, 1, 1, 0)), cfg)
  expect_true(is.na(w$wealth[6]))
  expect_equal(as.character(w$wealth[1:5]), c("Lowest", "Lowest", "Middle", "Highest", "Highest"))
})

test_that("participant flow counts are consistent", {
  w <- toy_women(10, v213 = c(0, rep(1, 9)), shbpbg = c(1, 0, rep(1, 8)),
                 v012 = c(25, 25, 17, rep(25, 7)),
                 sb267g = c(5, 5, 5, 995, 5, 5, 5.5, 7.2, 6, 6),
                 sb263a = c(1, 1, 1, 1, 0, 1, 1, 1, 1, 1),
                 sb236 = c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0), v437 = c(rep(600, 9), 9999))
  flow <- participant_flow(derive_variables(w, cfg), cfg)
  expect_equal(flow$n_remaining, c(10, 9, 8, 7, 6, 5, 4, 3, 2))
})

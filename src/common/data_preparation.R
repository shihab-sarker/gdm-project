# src/common/data_preparation.R
# Reading the BDHS 2022 files, linking women to their biomarker records,
# deriving the analysis variables and defining the analytic sample.
# DHS value codes follow the BDHS 2022 recode manual; the code meaning is
# written next to each use.

# ---- 1. Read and merge -------------------------------------------------------

# Women's file (IR): identifiers, design variables and covariates.
ir_vars <- c("v001", "v002", "v003",          # cluster, household, line number
             "v005", "v021", "v022",          # weight (6 decimals), PSU, stratum
             "v012", "v025", "v106", "v190", "v191", "v213", "v214",
             "v220", "v437", "v438", "v714")

# Household member file (PR): biomarker module for women aged 18+.
pr_vars <- c("hv001", "hv002", "hvidx", "hv104", # keys, sex of member
             "hv005", "hv105", "hv025",         # household weight, age of member, residence
             "sb240",                           # taking medication for diabetes
             "shbpbg",                          # household selected for biomarkers
             "sb236",                           # ever told she has diabetes
             "sb261a", "sb261b", "sb261c",      # clock time last ate (hhmm), appointment 1/2/3
             "sb262a", "sb262b", "sb262c",      # clock time last drank other than plain water (9595 = never)
             "sb263a", "sb263b", "sb263c",      # fasted >= 8 h at appointment 1/2/3
             "sb258d1", "sb258d2", "sb258d3",   # day of appointment 1/2/3
             "sb266d", "sb266t",                # day and clock time of the glucose test
             "sb267", "sb267g",                 # FPG in mg/dL; FPG in mmol/L (DHS)
             "wbp13", "wbp14",                  # 2nd BP reading, systolic/diastolic
             "wbp19",                           # taking BP medication
             "wbp24", "wbp25")                  # final (DHS) systolic/diastolic BP

read_and_merge <- function(cfg) {
  ir <- read_dta(here::here(cfg$path_ir), col_select = dplyr::all_of(ir_vars))
  pr <- read_dta(here::here(cfg$path_pr), col_select = dplyr::all_of(pr_vars))
  ir <- zap_labels(ir)
  pr <- zap_labels(pr)

  # Each woman in IR is the household member with the same cluster,
  # household and line number in PR.
  stopifnot(!anyDuplicated(ir[c("v001", "v002", "v003")]),
            !anyDuplicated(pr[c("hv001", "hv002", "hvidx")]))
  women <- dplyr::left_join(ir, pr, by = c(v001 = "hv001", v002 = "hv002", v003 = "hvidx"))

  # Checks: no woman lost or duplicated; every matched member is female (hv104 = 2)
  stopifnot(nrow(women) == nrow(ir), all(women$hv104 %in% c(2, NA)))
  cat("Women in IR:", nrow(ir), "| matched to PR:", sum(!is.na(women$hv104)), "\n")
  list(women = as.data.frame(women), pr = as.data.frame(pr))
}

# ---- 2. Derived variables ----------------------------------------------------

derive_variables <- function(w, cfg) {
  # DHS special codes become missing: 994-996 for glucose and blood pressure
  # (not present/refused/other), 9994-9999 for weight and height.
  na_codes <- function(x, codes) replace(x, x %in% codes, NA)
  glu_mmol <- na_codes(w$sb267g, 994:996)
  glu_mgdl <- na_codes(w$sb267, 994:996)

  # Fasting plasma glucose (mmol/L) on the released scale and the unconverted scale
  fpg_sb267g <- glu_mmol
  fpg_div18  <- glu_mgdl / 18
  released <- cfg$glucose_scale == "sb267g"
  w$fpg     <- if (released) fpg_sb267g else fpg_div18
  w$fpg_alt <- if (released) fpg_div18 else fpg_sb267g

  # Fasted >= 8 h (code 1) at any of the up to three appointments
  w$fasted <- (w$sb263a %in% 1) | (w$sb263b %in% 1) | (w$sb263c %in% 1)

  # Pre-existing diabetes (sb236): 1 = yes, 0 = no, 8/9 = don't know/refused
  # (missing). A missing answer is taken as "no" (counted in script 01).
  w$pre_dm <- ifelse(is.na(w$sb236) | w$sb236 == 0, 0, ifelse(w$sb236 == 1, 1, NA))

  # Outcomes. Women with pre-existing diabetes are excluded from all of them.
  #   gdm          : 5.1 <= FPG < 7.0 mmol/L; FPG >= 7.0 set to missing
  #   gdm_alt      : the same definition on the alternative glucose scale
  #   gdm_incl_dm  : FPG >= 5.1 including the diabetes range (sensitivity S1)
  #   fpg_dm_range : FPG >= 7.0 (diabetes range), reported separately
  gdm_def <- function(fpg) {
    ifelse(is.na(fpg) | w$pre_dm %in% c(1, NA) | fpg >= cfg$fpg_diabetes, NA,
           as.numeric(fpg >= cfg$fpg_gdm))
  }
  no_dm <- !is.na(w$fpg) & w$pre_dm %in% 0
  w$gdm          <- gdm_def(w$fpg)
  w$gdm_alt      <- gdm_def(w$fpg_alt)
  w$gdm_incl_dm  <- ifelse(no_dm, as.numeric(w$fpg >= cfg$fpg_gdm), NA)
  w$fpg_dm_range <- ifelse(no_dm, as.numeric(w$fpg >= cfg$fpg_diabetes), NA)
  # The same two outcomes on the alternative (unconverted) glucose scale
  no_dm_alt <- !is.na(w$fpg_alt) & w$pre_dm %in% 0
  w$gdm_incl_dm_alt  <- ifelse(no_dm_alt, as.numeric(w$fpg_alt >= cfg$fpg_gdm), NA)
  w$fpg_dm_range_alt <- ifelse(no_dm_alt, as.numeric(w$fpg_alt >= cfg$fpg_diabetes), NA)

  # Glucose test on the day of an appointment at which the woman reported
  # fasting >= 8 h (FALSE for the few women tested on another day).
  fasted_k <- cbind(w$sb263a %in% 1, w$sb263b %in% 1, w$sb263c %in% 1)
  col_or_na <- function(v) if (is.null(w[[v]])) rep(NA_real_, nrow(w)) else w[[v]]
  day_k    <- do.call(cbind, lapply(c("sb258d1", "sb258d2", "sb258d3"), col_or_na))
  same_day <- fasted_k & !is.na(day_k) & day_k == col_or_na("sb266d")
  same_day[is.na(same_day)] <- FALSE
  # without appointment and test dates (e.g. synthetic data) nobody is excluded
  w$tested_on_fasted_day <- if (is.null(w$sb266d)) rep(TRUE, nrow(w)) else rowSums(same_day) > 0

  # Hours from the last food and the last drink other than plain water to the
  # glucose test, at the appointment of the test (the first fasted appointment
  # for the few women tested on another day). Clock times are hhmm, so a time
  # later in the day than the test belongs to the previous day.
  first_true <- function(m) ifelse(rowSums(m) > 0, max.col(m, ties.method = "first"), NA)
  k <- ifelse(rowSums(same_day) > 0, first_true(same_day), first_true(fasted_k))
  at_k <- function(v) do.call(cbind, lapply(v, col_or_na))[cbind(seq_len(nrow(w)), k)]
  # any value that is not a clock time (e.g. 9595 never, 9898 don't know) is missing
  clock <- function(x) ifelse(x > 2359 | x %% 100 > 59, NA, x %/% 100 + (x %% 100) / 60)
  ate <- at_k(paste0("sb261", c("a", "b", "c")))
  drank <- at_k(paste0("sb262", c("a", "b", "c")))
  w$drank_never <- drank %in% 9595
  w$hours_ate   <- (clock(col_or_na("sb266t")) - clock(ate)) %% 24
  w$hours_drank <- (clock(col_or_na("sb266t")) - clock(drank)) %% 24
  w$drank_8h    <- !is.na(w$hours_drank) & w$hours_drank < 8

  # Duration of pregnancy: v214 is recorded in completed months, so a woman
  # reporting m months is, on average, m + 0.5 months pregnant.
  # Trimesters: completed months 1-3, 4-6, 7 and above (because months are
  # completed, the first group spans roughly weeks 4-17).
  # Weeks = (completed months + 0.5) x 4.345 (midpoint of the month).
  w$gest_weeks <- (w$v214 + 0.5) * 4.345
  w$gest_month <- w$v214
  w$trimester <- cut(w$v214, c(-Inf, 3, 6, Inf), labels = c("1st", "2nd", "3rd"))

  # Age (v012) in completed years
  w$age <- w$v012
  w$age_c <- w$v012 - 25                       # continuous age, centred (sensitivity)
  w$age_group     <- factor(ifelse(w$age >= 25, "25+", "<25"), levels = c("<25", "25+"))
  w$age_group_alt <- factor(ifelse(w$age >= 30, "30+", "<30"), levels = c("<30", "30+"))

  # Education (v106): 0 none, 1 primary, 2 secondary, 3 higher
  w$education <- factor(c("Primary or below", "Primary or below", "Secondary", "Higher")[w$v106 + 1],
                        levels = c("Primary or below", "Secondary", "Higher"))

  # Currently working (v714: 1 yes, 0 no)
  w$working <- factor(ifelse(w$v714 == 1, "Yes", "No"), levels = c("No", "Yes"))

  # Place of residence (v025: 1 urban, 2 rural); rural is the reference
  w$residence <- factor(ifelse(w$v025 == 1, "Urban", "Rural"), levels = c("Rural", "Urban"))

  # Birth order of the current pregnancy from the number of births so far (v220)
  w$birth_order <- factor(ifelse(w$v220 <= 1, "First", ifelse(w$v220 == 2, "Second", "Third or higher")),
                          levels = c("First", "Second", "Third or higher"))

  # Pre-pregnancy BMI, as in the anchor study: measured weight minus 3 kg,
  # and minus a further 0.5 kg for every week beyond week 20.
  weight_kg <- na_codes(w$v437, 9994:9999) / 10
  height_m  <- na_codes(w$v438, 9994:9999) / 1000
  gain <- 3 + 0.5 * pmax(w$gest_weeks - 20, 0)
  w$bmi <- (weight_kg - gain) / height_m^2
  w$bmi[w$bmi < 12] <- NA                       # implausible after correction
  # Measured (uncorrected) BMI during pregnancy, for a sensitivity analysis
  # with adjustment for trimester instead of the fixed weight-gain correction.
  w$bmi_measured <- weight_kg / height_m^2
  w$bmi_measured[w$bmi_measured < 12] <- NA

  # Hypertension: DHS final blood pressure (mean of readings 2 and 3)
  # >= 140/90 mmHg, or taking blood-pressure medication (wbp19 = 1).
  sbp <- na_codes(w$wbp24, 994:996)
  dbp <- na_codes(w$wbp25, 994:996)
  htn <- ifelse(is.na(sbp) | is.na(dbp), NA, ifelse(sbp >= 140 | dbp >= 90, "Yes", "No"))
  htn[w$wbp19 %in% 1] <- "Yes"
  w$hypertension <- factor(htn, levels = c("No", "Yes"))

  # Survey weight (v005 has six implied decimals) and wealth index score (v191)
  w$weight <- w$v005 / 1e6
  w$wealth_score <- w$v191

  # Analytic sample: currently pregnant (v213 = 1), in a household selected
  # for biomarkers (shbpbg = 1), aged 18+, with a valid glucose value after
  # fasting >= 8 h.
  w$analytic <- w$v213 %in% 1 & w$shbpbg %in% 1 & w$age >= 18 &
    !is.na(w$fpg) & w$fasted

  # Wealth tertiles: cut-points are the 1/3 and 2/3 quantiles of the wealth
  # score among women in the analytic sample (as in the anchor study).
  cuts <- quantile(w$wealth_score[w$analytic], c(1/3, 2/3))
  w$wealth <- cut(w$wealth_score, c(-Inf, cuts, Inf), labels = c("Lowest", "Middle", "Highest"))
  w$wealth[!w$analytic] <- NA

  # National wealth quintiles (v190) in three groups: poorest two, middle,
  # richest two. Level names match the tertiles so that terms line up in tables.
  w$wealth_v190 <- cut(w$v190, c(0, 2, 3, 5), labels = c("Lowest", "Middle", "Highest"))

  w
}

# ---- 3. Participant flow and the analytic sample -----------------------------

participant_flow <- function(w, cfg) {
  pregnant <- w$v213 %in% 1
  selected <- pregnant & w$shbpbg %in% 1
  adult    <- selected & w$age >= 18
  glucose  <- adult & !is.na(w$fpg)
  fasted   <- glucose & w$fasted
  no_dm    <- fasted & w$pre_dm %in% 0
  outcome  <- fasted & !is.na(w$gdm)
  complete <- outcome & complete.cases(w[cfg$covariates])
  stopifnot(identical(fasted, w$analytic))

  n <- c(nrow(w), sum(pregnant), sum(selected), sum(adult), sum(glucose),
         sum(fasted), sum(no_dm), sum(outcome), sum(complete))
  data.frame(
    step = c("Women interviewed",
             "Currently pregnant",
             "In household selected for biomarker testing",
             "Aged 18 years or older",
             "Valid fasting plasma glucose value",
             "Fasted at least 8 hours (analytic sample)",
             "No self-reported pre-existing diabetes",
             "FPG below 7.0 mmol/L (outcome sample)",
             "Complete covariate data (model sample)"),
    n_remaining = n,
    n_excluded = c(NA, -diff(n))
  )
}

# Agreement between the glucose test day and the day of an appointment at
# which the woman reported fasting >= 8 h (analytic sample).
fasting_check <- function(w) {
  a <- w[w$analytic, ]
  fasted_k <- cbind(a$sb263a %in% 1, a$sb263b %in% 1, a$sb263c %in% 1)
  day_k    <- cbind(a$sb258d1, a$sb258d2, a$sb258d3)
  same_day <- fasted_k & !is.na(day_k) & day_k == a$sb266d
  same_day[is.na(same_day)] <- FALSE
  first_fasted <- c("1st", "2nd", "3rd")[max.col(fasted_k, ties.method = "first")]
  tested_on    <- ifelse(rowSums(same_day) > 0,
                         paste0("yes (", c("1st", "2nd", "3rd")[max.col(same_day, ties.method = "first")], ")"),
                         "no")
  out <- as.data.frame(table(first_fasted_appointment = first_fasted,
                             test_on_day_of_fasted_appointment = tested_on))
  out[out$Freq > 0, ]
}

# Answers to sb236 (ever told she has diabetes) in the analytic sample
diabetes_history_counts <- function(w) {
  dx <- w$sb236[w$analytic]
  data.frame(answer = c("yes (excluded)", "no", "don't know/refused (excluded)", "missing (treated as no)"),
             n = c(sum(dx %in% 1), sum(dx %in% 0), sum(dx %in% c(8, 9)), sum(is.na(dx))))
}

# Relation between the two glucose variables among all household members
# with both values (confirms sb267g = sb267 x 1.11 / 18).
glucose_scale_relation <- function(pr) {
  ok <- !is.na(pr$sb267) & !is.na(pr$sb267g) & !pr$sb267 %in% 994:996 & !pr$sb267g %in% 994:996
  factor <- pr$sb267g[ok] / (pr$sb267[ok] / 18)
  data.frame(n_pairs = sum(ok), n_women = sum(pr$hv104[ok] %in% 2), n_men = sum(pr$hv104[ok] %in% 1),
             implied_factor_median = median(factor),
             implied_factor_p01 = unname(quantile(factor, 0.01)),
             implied_factor_p99 = unname(quantile(factor, 0.99)),
             share_equal_to_rounded_sb267x1.11div18 =
               mean(abs(pr$sb267g[ok] - round(pr$sb267[ok] * 1.11 / 18, 1)) <= 0.051))
}

# Sampling weights available in the IR and PR files (variable labels only):
# BDHS 2022 releases no weight specific to the biomarker subsample.
weight_variables <- function(cfg) {
  labs <- function(path, file) {
    d <- read_dta(here::here(path), n_max = 0)
    l <- vapply(d, function(x) if (is.null(attr(x, "label"))) "" else attr(x, "label"), "")
    data.frame(file = file, variable = names(l), label = unname(l))
  }
  all <- rbind(labs(cfg$path_ir, "IR"), labs(cfg$path_pr, "PR"))
  all[grepl("weight", all$label, ignore.case = TRUE) & grepl("decimal", all$label, ignore.case = TRUE), ]
}

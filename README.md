# Glucose calibration, test timing and socioeconomic patterning of fasting-glucose-defined gestational diabetes, BDHS 2022

R analysis code for the manuscript *Glucose calibration, test timing and socioeconomic patterning of fasting-glucose-defined gestational diabetes: A cross-sectional analysis of the 2022 Bangladesh Demographic and Health Survey*, prepared for *SSM – Population Health*.

The analysis studies pregnant women aged 18 years or older in the BDHS 2022 biomarker subsample. GDM is defined as a single fasting plasma glucose (FPG) of 5.1–6.9 mmol/L (FPG-defined GDM). Every result is reported on two glucose scales: the released variable and the unconverted scale. The analysis estimates:

- survey-weighted prevalence, overall, by trimester and by month
- prevalence ratios for household wealth and place of residence (crude, minimal and adjusted models)
- the Erreygers corrected concentration index, with the Wagstaff and relative indices, and its decomposition
- model diagnostics: collinearity, fit, influence, Firth regression and calibration
- sensitivity analyses S1–S12
- checks that do not depend on the glucose scale or on one threshold, and a design-based rerun of the 2017–18 survey

## Layout

```
config/config.yml          shared settings, then one block per survey (bdhs2022, bdhs2017_18) and one for the tests
src/common/                functions used by both surveys (setup.R sources every function file)
src/bdhs2022/              manuscript tables and figures
src/bdhs2017_18/           reading and renaming the 2017–18 recode, rounding-pattern test
analysis/bdhs2022/         steps 01-10; each reads from data/ or output/ and writes its own results
analysis/bdhs2017_18/      step 11, the 2017–18 rerun
analysis/run_all.R         runs 01-11 in order, for both surveys
tests/                     synthetic-data generator and testthat tests
data/raw/bdhs2022/         BDHS 2022 recode files (not tracked)
data/raw/bdhs2017_18/      BDHS 2017–18 recode files (not tracked)
data/derived/<survey>/     derived data sets (not tracked)
output/bdhs2022/tables     results as CSV
output/bdhs2022/figures    figures: vector PDF with embedded Times New Roman, and PNG at 1000 dpi
output/bdhs2017_18/tables  tables of the 2017–18 rerun
docs/manuscript/           current manuscript, supplement, title page, highlights; extra/ working documents; grammer/ grammar parts
docs/journal/              journal guide for authors and a sample paper
docs/STROBE/               STROBE checklist
docs/anchor-paper/         anchor paper
docs/data-description/     BDHS 2022 final report and DDI documentation
references/                EndNote exports of the reference list
archive/                   earlier code, documents and outputs (not tracked)
```

| Script | What it does |
|---|---|
| `src/common/data_preparation.R` | reads the IR and PR files, links them, derives variables and the participant flow |
| `src/common/survey_design.R` | domain design on all PSUs, Rao–Wu bootstrap replicates, prevalence, Table 1 tests |
| `src/common/regression.R` | survey-weighted Poisson PR, logistic OR, Firth regression, GVIF and survey VIF, Archer–Lemeshow test, bootstrap calibration, influence (DFBETAS, modified Cook's distance) |
| `src/common/inequality.R` | fractional rank, concentration indices, decomposition, bootstrap percentile CIs |
| `src/common/validity_checks.R` | FPG distributions, prevalence by month, test timing by residence and wealth, selection comparison, replication of the survey report |
| `src/bdhs2022/tables_figures.R` | number formatting, manuscript tables, figures |
| `src/bdhs2017_18/bdhs2017.R` | reading the 2017–18 recode files, diabetes in 2017–18, rounding-pattern test |

| Analysis step | Output |
|---|---|
| 01 prepare data | analytic sample, sample flow, fasting and diabetes-history checks, glucose-scale relation, replication of the BDHS 2022 report diabetes estimate |
| 02 prevalence | prevalence on both scales, by trimester and month; design effect; PSU structure |
| 03 characteristics | Table 1 on both scales |
| 04 regression | per scale: crude, minimal and adjusted PR; OR; specification with gestational weeks; model with hypertension; Firth |
| 05 model diagnostics | per scale: events per parameter, GVIF, survey VIF, Archer–Lemeshow, calibration, influence with Poisson refits |
| 06 inequality | per scale: relative, Wagstaff and Erreygers indices with bootstrap CIs, concentration curves, decomposition |
| 07 sensitivity | S1–S12, each on both scales |
| 08 validity checks | FPG distributions, prevalence by month, timing of the test by residence and wealth, selection comparison |
| 09 tables and figures | manuscript tables, `Figure_1` (sample flow), `Figure_2` (inequality), `Figure_SA1` (calibration), `Figure_SA3` (FPG distribution) |
| 10 reviewer checks | log-FPG regression, threshold sweep 4.4–5.4 mmol/L with `Figure_SA2`, FPG by month and trimester, fasting duration, slope per month, relative concentration index of continuous FPG |
| 11 BDHS 2017–18 | the same analysis on the 2017–18 survey, if its files are in `data/raw/bdhs2017_18/`; results in `output/bdhs2017_18/tables` |

## Final documents

The .docx files in `docs/manuscript/` are the source of truth for the text; no generator builds them, so edit them directly. The manuscript has 59 references. The supplement ends with its own methods note ("Supplementary methods: model checks") and reference list ("Supplementary references", 7 works); other works it cites are in the main list.

| File (GDM_BDHS2022_…) | Content | Pages |
|---|---|---|
| `manuscript.docx`, `manuscript.pdf` | manuscript with Tables 1–2 and figure legends | 33 |
| `supplementary.docx` | Tables SA1–SA16, Figs SA1–SA3, methods note and references | 24 |
| `title_page.docx` | title page and declarations | 2 |
| `highlights.docx` | highlights | 1 |
| `submission_todo.docx` | to-do before submission | 5 |
| `change_log.docx` | change log | 28 |
| `claim_verification.docx` | claim verification | 20 |
| `defence_guide_merged.docx` | reviewer questions and answers | 11 |
| `grammar_part1_title_abstract_introduction.docx` | text for the grammar check | 6 |
| `grammar_part2_methods.docx` | text for the grammar check | 5 |
| `grammar_part3_results.docx` | text for the grammar check | 3 |
| `grammar_part4_discussion.docx` | text for the grammar check | 5 |
| `grammar_part5_conclusions_statements_captions.docx` | text for the grammar check | 6 |

## Running

1. Request the BDHS 2022 data from the DHS Program (https://dhsprogram.com) and place `BDIR81FL.DTA` and `BDPR81FL.DTA` in `data/raw/bdhs2022/`. DHS microdata must never be committed; `data/` is git-ignored.
2. Restore packages with `Rscript -e "renv::restore()"`. The analysis used R 4.6.1 with survey 4.5, car 3.1-5, svydiags 0.7 and logistf 1.26.1 (versions are locked in `renv.lock`).
3. Run the analysis with `Rscript analysis/run_all.R` (about 40 seconds). To keep a log: `Rscript analysis/run_all.R > output/bdhs2022/run_log.txt 2>&1`.
4. Run the tests with `Rscript -e "testthat::test_dir('tests/testthat')"`. For the end-to-end run on synthetic data (B = 50), set `RUN_PIPELINE_TEST=true` first.
5. Figures are written in the upload format: one vector PDF with embedded fonts and one 1000-dpi PNG per figure.

## Main analytic decisions

- **Sample.** Ever-married women aged 18–49 years who were currently pregnant, lived in a household selected for biomarker testing, had a valid glucose value and had not eaten for at least 8 hours at one of up to three appointments. Women ever told that they had diabetes are excluded; FPG ≥ 7.0 mmol/L is reported separately as the diabetes range. Model samples: 264 (released) and 267 (unconverted).
- **Glucose scale.** Released-scale GDM (sb267g 5.1–6.9) equals an unconverted FPG of about 4.6–6.3 mmol/L, so the two scales differ only in where the threshold falls; their estimates bound the plausible value.
  - `sb267g` (mmol/L, released) equals `sb267 × 1.11 / 18`, the whole-blood-to-plasma factor, although the survey report states that the HemoCue system already gave plasma-equivalent values.
  - `sb267g` reproduces the report's diabetes estimate for women aged 18+ (16.8%); `sb267 / 18` gives 9.5%, close to the 2017–18 level (10%).
  - Neither scale is primary: the released scale is analysed because the DHS Program distributes and reports it, and the unconverted scale (`sb267 / 18`) in parallel throughout. The two are the bounds of the plausible prevalence; the 2017–18 match and the rounding-pattern test favour the lower bound.
- **Survey design.** Weight `v005`/1,000,000; PSU `v021`; the 16 strata defined for sampling errors `v022`. Estimation within the domain of analytic women, keeping all PSUs; CIs use the design degrees of freedom.
- **Exposures and covariates.** Wealth tertiles from unweighted cut-points of the wealth score in the analytic sample; the continuous score ranks women for the concentration index. Covariates follow Mazumder et al. (2022): age group, residence, wealth tertile, trimester (in place of gestational weeks) and estimated pre-pregnancy BMI. No variable selection.
- **Diagnostics.** Run on the logistic version of the adjusted model; the Poisson model is refitted without women flagged by DFBETAS or modified Cook's distance.
- **Decomposition.** Linear probability model; wealth tertiles are excluded because wealth is the ranking variable.
- **Uncertainty for the indices and decomposition.** Rao–Wu rescaled bootstrap within strata, 1,000 replicates, seed 123, percentile intervals.

## Data and code availability

The BDHS 2022 data are available from the DHS Program to registered users on application and may not be redistributed. This repository contains the analysis code, configuration and synthetic test data.

## Citation

If you use this code, please cite the article (details to be added on publication) and the dataset: National Institute of Population Research and Training (NIPORT) & ICF (2024). *Bangladesh Demographic and Health Survey 2022* (BDIR81FL, BDPR81FL). The DHS Program, ICF.

## Supplementary numbering

Supplementary tables and figures are numbered SA1, SA2, … in the order of their first citation in the main text. The CSV files in `output/bdhs2022/tables` keep their descriptive names, and sensitivity scenarios keep the labels S0–S12.

| Supplementary item | Output |
|---|---|
| Table SA1 | fasting_check, diabetes_history |
| Table SA2 | none (STROBE checklist) |
| Table SA3 | glucose_scale_relation, dhs_report_replication, fpg_distribution, glucose_heaping, output/bdhs2017_18/tables |
| Table SA4 | selection_comparison |
| Table SA5 | manuscript_table_hypertension, firth |
| Table SA6 | epv, vif_generalised, vif_survey, archer_lemeshow, linearity_tests, bootstrap_validation, influence_refits |
| Table SA7 | calibration_groups |
| Table SA8 | manuscript_table_sensitivity |
| Table SA9 | paired_classification, common_sample |
| Table SA10 | design_effect, psu_structure, design_comparison |
| Table SA11 | prevalence_by_month |
| Table SA12 | timing_by_ses |
| Table SA13 | manuscript_table_or |
| Table SA14 | manuscript_table_anchor |
| Table SA15 | reviewer_log_fpg, reviewer_threshold_sweep |
| Table SA16 | manuscript_table_decomp_full, manuscript_table_decomp_wealth |
| Fig. SA1 | output/bdhs2022/figures/Figure_SA1 (calibration) |
| Fig. SA2 | output/bdhs2022/figures/Figure_SA2 (threshold sweep) |
| Fig. SA3 | output/bdhs2022/figures/Figure_SA3 (FPG distribution) |

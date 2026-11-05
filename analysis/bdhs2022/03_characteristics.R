# 03_characteristics.R
# Table 1: characteristics of pregnant women and prevalence of GDM by
# characteristic, on both glucose scales.

source(here::here("src", "common", "setup.R"))
women  <- readRDS(data_file("women"))
sample <- readRDS(data_file("analytic_sample"))
des <- svy_domain(sample, women)

vars <- c("age_group", "education", "working", "residence", "wealth", "birth_order",
          "trimester", "hypertension")
t1_cat  <- do.call(rbind, lapply(scales, function(o) describe_categorical(des, vars, o)))
t1_cont <- do.call(rbind, lapply(scales, function(o) describe_continuous(des, c("age", "gest_weeks", "bmi"), o)))
print(t1_cat)
save_table(t1_cat, "table1_categorical")
save_table(t1_cont, "table1_continuous")
save_table(residence_by_wealth(des), "residence_by_wealth")

# Table 1 of the manuscript: characteristics by household wealth tertile
t1_wealth <- describe_by_wealth(des, vars, c("age", "gest_weeks", "bmi"), as.list(scales))
print(t1_wealth)
save_table(t1_wealth, "table1_by_wealth")

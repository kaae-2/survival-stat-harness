### Survival Stats - Report
# Purpose: Write a statistical report
# perform a statistical analysis of a dataset and present the results
# (in a separate document)
# assume no prior knowledge of R or the dataset
# write in academic language
library(haven)
source("helpers.R")
# Plan:
# 0.0 define the hypotheses
# 0.1 load the dataset
readRenviron(".Renviron")

df <- read_sav(get_env_required("SAV_DATA_PATH"))
date_df <- read_sav(get_env_required("SAV_DATE_PATH"))

# 0.2 clean the datasets, keeping only relevant columns
df <- rename_from_env(df, c(
  id = "COL_ID",
  included="COL_INCLUDED",
  type="COL_TYPE",
  sex="COL_SEX",
  diagnosis="COL_DATE_DIAGNOSIS",
  age="COL_AGE",
  risk_grp="COL_RISK_GROUP",
  event="COL_EVENT",
  date_event="COL_EVENT_DATE",
  date_death="COL_DEATH_DATE",
  lymf_count="COL_LYMF_COUNT"
))

df <- df[!is.na(df$included) & df$included == 1,]

date_df <- rename_from_env(date_df, c(
  id = "COL_ID",
  landmark_date="COL_LANDMARK_DATE",
  trsct_date="COL_TRSCT_DATE"
))


# 0.3 plan the report generations
# 1.0 plan out the models

# --- initial ideas for the report ---
# 1.1 Competing risks
# 1.2 landmark analysis
# 1.3 stratification based on ...
# 1.4 Goodness of Fit analysis
# 1.5 Restricted Mean Survival Time
# 1.6 testing of significance

# --- expected outputs ---
# 2.1 stratified hazard curves
# 2.2 Goodness of Fit analysis
# 2.3 Survival curves for restricted mean
# 2.4 Table of the different outcomes (hazard rates, confidence intervals), significance etc.

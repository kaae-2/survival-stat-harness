### Survival Stats - Report
# Purpose: Write a statistical report
# perform a statistical analysis of a dataset and present the results
# (in a separate document)
# assume no prior knowledge of R or the dataset
# write in academic language

# Plan:
# 0.0 define the hypotheses
# 0.1 load the dataset
# 0.2 clean the datasets, keeping only relevant columns
pipeline_env <- new.env(parent = globalenv())
sys.source("load.R", envir=pipeline_env)
df <- pipeline_env$load_dataset()

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

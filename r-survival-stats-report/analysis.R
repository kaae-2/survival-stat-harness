### Survival Stats - Report
# Purpose: Write a statistical report
# perform a statistical analysis of a dataset and present the results
# (in a separate document)
# assume no prior knowledge of R or the dataset
# write in academic language
#
# Analysis purpose: fit the planned survival analyses in one sequential script
# and keep the resulting objects available for inspection at the end.
#
# Plan:
# 0.0 define the hypotheses
# 0.1 load the dataset
# 0.2 clean the datasets, keeping only relevant columns
# 0.3 plan the report generations
# 0.4 make the data variables for survival stat
# 1.0 plan out the models
#
# --- initial ideas for the report ---
# 1.1 Competing risks
# 1.2 landmark analysis
# 1.3 stratification based on risk group and baseline covariates
# 1.4 Goodness of Fit analysis
# 1.5 Restricted Mean Survival Time
# 1.6 testing of significance
#
# --- expected outputs ---
# 2.1 stratified hazard curves
# 2.2 Goodness of Fit analysis
# 2.3 Survival curves for restricted mean
# 2.4 Table of the different outcomes (hazard rates, confidence intervals, significance etc.)

# --- setup / report scaffolding ---
# Use explicit namespace calls so rerunning the script is idempotent and does
# not depend on the current attached package state.
required_packages <- c("survival", "prodlim", "cmprsk", "mets")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

out_dir <- "out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

report_times_days <- c(365.25, 3 * 365.25, 5 * 365.25)
analysis_helpers <- new.env(parent = globalenv())
sys.source("analysis_helpers.R", envir = analysis_helpers, chdir = TRUE)

# --- 0.1 load the dataset ---
# The dataset is loaded a single time and then reused throughout the script.
pipeline_env <- new.env(parent = globalenv())
sys.source("load.R", envir = pipeline_env)
df <- pipeline_env$load_dataset()

# --- 0.2 and 0.4 clean and transform the analysis dataset ---
# Keep the original loaded object untouched and prepare a separate analysis data
# frame with labels and convenience variables used in several downstream blocks.
analysis_df <- transform(
  df,
  sex = factor(sex),
  risk_grp = factor(risk_grp),
  years = days / 365.25,
  totalyears = totaldays / 365.25,
  compriskyears = compriskdays / 365.25,
  log_lymf_count = ifelse(is.na(lymf_count), NA_real_, log1p(pmax(lymf_count, 0))),
  sct_by_landmark = ifelse(
    is.na(landmark_date),
    NA_character_,
    ifelse(!is.na(sct_date) & sct_date <= landmark_date, "SCT by landmark", "No SCT by landmark")
  )
)
analysis_df$sct_by_landmark <- factor(
  analysis_df$sct_by_landmark,
  levels = c("No SCT by landmark", "SCT by landmark")
)

# --- 0.3 cohort summaries for the report ---
analysis_cohort_summary <- data.frame(
  n_patients = nrow(analysis_df),
  n_events = sum(analysis_df$event == 1, na.rm = TRUE),
  n_censored = sum(analysis_df$event == 0, na.rm = TRUE),
  n_relapse = sum(analysis_df$comprisk == 1, na.rm = TRUE),
  n_death = sum(analysis_df$comprisk == 2, na.rm = TRUE),
  n_smn = sum(analysis_df$comprisk == 3, na.rm = TRUE)
)

analysis_risk_levels <- levels(analysis_df$risk_grp)
analysis_baseline_by_risk_group <- data.frame(
  risk_grp = analysis_risk_levels,
  n = as.integer(table(analysis_df$risk_grp)[analysis_risk_levels]),
  n_events = as.numeric(tapply(analysis_df$event == 1, analysis_df$risk_grp, sum, na.rm = TRUE)[analysis_risk_levels]),
  n_relapse = as.numeric(tapply(analysis_df$comprisk == 1, analysis_df$risk_grp, sum, na.rm = TRUE)[analysis_risk_levels]),
  n_death = as.numeric(tapply(analysis_df$comprisk == 2, analysis_df$risk_grp, sum, na.rm = TRUE)[analysis_risk_levels]),
  n_smn = as.numeric(tapply(analysis_df$comprisk == 3, analysis_df$risk_grp, sum, na.rm = TRUE)[analysis_risk_levels]),
  median_age = as.numeric(tapply(analysis_df$age, analysis_df$risk_grp, stats::median, na.rm = TRUE)[analysis_risk_levels]),
  q1_age = as.numeric(tapply(analysis_df$age, analysis_df$risk_grp, stats::quantile, probs = 0.25, na.rm = TRUE)[analysis_risk_levels]),
  q3_age = as.numeric(tapply(analysis_df$age, analysis_df$risk_grp, stats::quantile, probs = 0.75, na.rm = TRUE)[analysis_risk_levels]),
  median_lymf_count = as.numeric(tapply(analysis_df$lymf_count, analysis_df$risk_grp, stats::median, na.rm = TRUE)[analysis_risk_levels]),
  median_followup_days = as.numeric(tapply(analysis_df$totaldays, analysis_df$risk_grp, stats::median, na.rm = TRUE)[analysis_risk_levels])
)
analysis_baseline_dtable <- mets::dtable(analysis_df, ~ risk_grp + sex, level = 1)

analysis_helpers$save_csv(analysis_cohort_summary, "cohort_summary.csv")
analysis_helpers$save_csv(analysis_baseline_by_risk_group, "baseline_by_risk_group.csv")
analysis_helpers$save_text(capture.output(analysis_baseline_dtable), "baseline_dtable.txt")

# --- 1.3 and 1.6 overall event-free survival and significance testing ---
# This follows the standard course workflow: non-parametric curves first, then
# log-rank testing, then an adjusted Cox model.
overall_km_df <- subset(
  analysis_df,
  !is.na(days) & !is.na(event) & !is.na(risk_grp)
)

overall_model_df <- subset(
  overall_km_df,
  !is.na(sex) & !is.na(age) & !is.na(log_lymf_count)
)

overall_report_times_days <- report_times_days[
  report_times_days <= max(overall_km_df$days, na.rm = TRUE)
]

overall_km_fit <- prodlim::prodlim(
  survival::Surv(days, event) ~ risk_grp,
  data = overall_km_df
)

overall_cumhaz_fit <- survival::survfit(
  survival::Surv(days, event) ~ risk_grp,
  data = overall_km_df
)

overall_logrank_test <- survival::survdiff(
  survival::Surv(days, event) ~ risk_grp,
  data = overall_km_df
)

overall_adjusted_cox_fit <- survival::coxph(
  survival::Surv(days, event) ~ risk_grp + sex + age + log_lymf_count,
  data = overall_model_df,
  x = TRUE
)

overall_event_free_summary <- analysis_helpers$extract_prodlim_summary(
  overall_km_fit,
  overall_report_times_days
)

overall_adjusted_cox_table <- analysis_helpers$extract_cox_table(overall_adjusted_cox_fit)

analysis_helpers$save_csv(overall_event_free_summary, "overall_event_free_summary.csv")
analysis_helpers$save_csv(overall_adjusted_cox_table, "overall_adjusted_cox_table.csv")
analysis_helpers$save_text(capture.output(overall_logrank_test), "overall_logrank_test.txt")
analysis_helpers$save_text(capture.output(summary(overall_adjusted_cox_fit)), "overall_adjusted_cox_summary.txt")

# --- 2.1 stratified hazard curves ---
analysis_helpers$save_plot_png("overall_km_by_risk_group.png", function() {
  group_levels <- levels(overall_km_df$risk_grp)
  group_cols <- seq_along(group_levels)

  plot(
    overall_km_fit,
    col = group_cols,
    lty = 1,
    xlab = "Days since diagnosis",
    ylab = "Event-free survival probability",
    main = "Event-free survival by risk group"
  )
  graphics::legend(
    "bottomleft",
    legend = group_levels,
    col = group_cols,
    lty = 1,
    bty = "n"
  )
})

analysis_helpers$save_plot_png("overall_cumulative_hazard_by_risk_group.png", function() {
  group_levels <- levels(overall_km_df$risk_grp)
  group_cols <- seq_along(group_levels)

  graphics::plot(
    overall_cumhaz_fit,
    fun = "cumhaz",
    col = group_cols,
    lty = 1,
    xlab = "Days since diagnosis",
    ylab = "Cumulative hazard",
    main = "Cumulative hazard by risk group"
  )
  graphics::legend(
    "topleft",
    legend = group_levels,
    col = group_cols,
    lty = 1,
    bty = "n"
  )
})

# --- 1.1 competing risks ---
# This mirrors the course examples: non-parametric cumulative incidence curves,
# Gray's test, and then cause-specific Cox models for the main event types.
competing_risk_df <- subset(
  analysis_df,
  !is.na(compriskdays) & !is.na(comprisk) & !is.na(risk_grp)
)

competing_risk_report_times_days <- report_times_days[
  report_times_days <= max(competing_risk_df$compriskdays, na.rm = TRUE)
]

competing_risk_cif_fit <- prodlim::prodlim(
  prodlim::Hist(compriskdays, comprisk) ~ risk_grp,
  data = competing_risk_df
)

competing_risk_gray_fit <- cmprsk::cuminc(
  ftime = competing_risk_df$compriskdays,
  fstatus = competing_risk_df$comprisk,
  group = competing_risk_df$risk_grp,
  cencode = 0
)

competing_risk_gray_tests <- as.data.frame(competing_risk_gray_fit$Tests)
competing_risk_gray_tests$comparison <- rownames(competing_risk_gray_tests)
rownames(competing_risk_gray_tests) <- NULL
competing_risk_gray_timepoints <- analysis_helpers$extract_cuminc_timepoints(
  competing_risk_gray_fit,
  competing_risk_report_times_days
)

competing_risk_relapse_summary <- analysis_helpers$extract_prodlim_summary(competing_risk_cif_fit, competing_risk_report_times_days, cause = 1)
competing_risk_death_summary <- analysis_helpers$extract_prodlim_summary(competing_risk_cif_fit, competing_risk_report_times_days, cause = 2)
competing_risk_smn_summary <- analysis_helpers$extract_prodlim_summary(competing_risk_cif_fit, competing_risk_report_times_days, cause = 3)

competing_risk_relapse_cox_fit <- analysis_helpers$fit_cause_specific_cox(analysis_df, 1)
competing_risk_death_cox_fit <- analysis_helpers$fit_cause_specific_cox(analysis_df, 2)
competing_risk_smn_cox_fit <- analysis_helpers$fit_cause_specific_cox(analysis_df, 3)

competing_risk_relapse_cox_table <- analysis_helpers$extract_cox_table(competing_risk_relapse_cox_fit)
competing_risk_death_cox_table <- analysis_helpers$extract_cox_table(competing_risk_death_cox_fit)
competing_risk_smn_cox_table <- analysis_helpers$extract_cox_table(competing_risk_smn_cox_fit)

analysis_helpers$save_csv(competing_risk_gray_tests, "competing_risk_gray_tests.csv")
analysis_helpers$save_csv(competing_risk_gray_timepoints, "competing_risk_gray_timepoints.csv")
analysis_helpers$save_csv(competing_risk_relapse_cox_table, "competing_risk_relapse_cox_table.csv")
analysis_helpers$save_text(capture.output(summary(competing_risk_relapse_cox_fit)), "competing_risk_relapse_cox_summary.txt")
analysis_helpers$save_csv(competing_risk_death_cox_table, "competing_risk_death_cox_table.csv")
analysis_helpers$save_text(capture.output(summary(competing_risk_death_cox_fit)), "competing_risk_death_cox_summary.txt")
analysis_helpers$save_csv(competing_risk_smn_cox_table, "competing_risk_smn_cox_table.csv")
analysis_helpers$save_text(capture.output(summary(competing_risk_smn_cox_fit)), "competing_risk_smn_cox_summary.txt")

analysis_helpers$save_plot_png("competing_risk_cif_by_risk_group.png", function() {
  graphics::par(mfrow = c(1, 3))

  plot(
    competing_risk_cif_fit,
    cause = 1,
    xlab = "Days since diagnosis",
    ylab = "Cumulative incidence",
    main = "Relapse"
  )

  plot(
    competing_risk_cif_fit,
    cause = 2,
    xlab = "Days since diagnosis",
    ylab = "Cumulative incidence",
    main = "Death"
  )

  plot(
    competing_risk_cif_fit,
    cause = 3,
    xlab = "Days since diagnosis",
    ylab = "Cumulative incidence",
    main = "SMN"
  )
}, width = 2400, height = 900)

# --- 1.2 landmark analysis for SCT ---
# Patients are classified by SCT status at the landmark date and then followed
# from the landmark onward. This avoids using future SCT information at baseline.
landmark_source_df <- subset(
  analysis_df,
  !is.na(landmark_date) & !is.na(date) & !is.na(event) & !is.na(sct_by_landmark)
)

landmark_analysis_df <- transform(
  subset(landmark_source_df, date > landmark_date),
  time_from_landmark = as.numeric(date - landmark_date),
  event_from_landmark = event
)

landmark_group_counts <- data.frame(
  total_with_landmark_information = nrow(landmark_source_df),
  included_after_landmark = nrow(landmark_analysis_df),
  excluded_event_by_landmark = sum(
    landmark_source_df$date <= landmark_source_df$landmark_date & landmark_source_df$event == 1,
    na.rm = TRUE
  ),
  excluded_censored_by_landmark = sum(
    landmark_source_df$date <= landmark_source_df$landmark_date & landmark_source_df$event == 0,
    na.rm = TRUE
  )
)

landmark_sct_counts <- as.data.frame(table(landmark_analysis_df$sct_by_landmark, useNA = "ifany"))
names(landmark_sct_counts) <- c("sct_by_landmark", "n")

landmark_km_fit <- prodlim::prodlim(
  survival::Surv(time_from_landmark, event_from_landmark) ~ sct_by_landmark,
  data = landmark_analysis_df
)

landmark_model_df <- subset(
  landmark_analysis_df,
  !is.na(risk_grp) & !is.na(sex) & !is.na(age) & !is.na(log_lymf_count)
)

landmark_adjusted_cox_fit <- survival::coxph(
  survival::Surv(time_from_landmark, event_from_landmark) ~ sct_by_landmark + risk_grp + sex + age + log_lymf_count,
  data = landmark_model_df,
  x = TRUE
)

landmark_report_times_days <- report_times_days[
  report_times_days <= max(landmark_analysis_df$time_from_landmark, na.rm = TRUE)
]

landmark_event_free_summary <- analysis_helpers$extract_prodlim_summary(landmark_km_fit, landmark_report_times_days)

landmark_adjusted_cox_table <- analysis_helpers$extract_cox_table(landmark_adjusted_cox_fit)

analysis_helpers$save_csv(landmark_group_counts, "landmark_group_counts.csv")
analysis_helpers$save_csv(landmark_sct_counts, "landmark_sct_counts.csv")
analysis_helpers$save_csv(landmark_event_free_summary, "landmark_event_free_summary.csv")
analysis_helpers$save_csv(landmark_adjusted_cox_table, "landmark_adjusted_cox_table.csv")
analysis_helpers$save_text(capture.output(summary(landmark_adjusted_cox_fit)), "landmark_adjusted_cox_summary.txt")

analysis_helpers$save_plot_png("landmark_km_by_sct_status.png", function() {

  sct_levels <- levels(landmark_analysis_df$sct_by_landmark)
  sct_cols <- seq_along(sct_levels)

  plot(
    landmark_km_fit,
    col = sct_cols,
    lty = 1,
    xlab = "Days since landmark",
    ylab = "Post-landmark event-free survival",
    main = "Landmark analysis by SCT status"
  )
  graphics::legend(
    "bottomleft",
    legend = sct_levels,
    col = sct_cols,
    lty = 1,
    bty = "n"
  )
})

# --- 1.5 restricted mean survival time ---
# RMST provides an absolute time-scale summary that complements the hazard ratio.
rmst_tau_days <- min(5 * 365.25, max(overall_km_df$days, na.rm = TRUE))

overall_rmst_fit <- survival::survfit(
  survival::Surv(days, event) ~ risk_grp,
  data = overall_km_df
)

overall_rmst_table <- analysis_helpers$extract_rmst_table(overall_rmst_fit, rmst_tau_days)

competing_risk_time_lost_table <- do.call(
  rbind,
  list(
    analysis_helpers$extract_time_lost_by_cause(competing_risk_df, rmst_tau_days, 1, "Relapse"),
    analysis_helpers$extract_time_lost_by_cause(competing_risk_df, rmst_tau_days, 2, "Death"),
    analysis_helpers$extract_time_lost_by_cause(competing_risk_df, rmst_tau_days, 3, "SMN")
  )
)

analysis_helpers$save_csv(overall_rmst_table, "overall_rmst_table.csv")
analysis_helpers$save_csv(competing_risk_time_lost_table, "competing_risk_time_lost_table.csv")

analysis_helpers$save_plot_png("overall_rmst_curves_by_risk_group.png", function() {
  group_levels <- levels(overall_km_df$risk_grp)
  group_cols <- seq_along(group_levels)

  graphics::plot(
    overall_rmst_fit,
    col = group_cols,
    lty = 1,
    xlab = "Days since diagnosis",
    ylab = "Event-free survival probability",
    xlim = c(0, rmst_tau_days),
    main = "Event-free survival up to the RMST horizon"
  )
  graphics::abline(v = rmst_tau_days, lty = 2)
  graphics::legend(
    "bottomleft",
    legend = group_levels,
    col = group_cols,
    lty = 1,
    bty = "n"
  )
})

# --- 1.4 goodness-of-fit and model diagnostics ---
# The full Cox model is checked for proportional hazards, and the continuous
# covariates are checked for linearity on the log-hazard scale.
overall_phreg_fit <- mets::phreg(
  survival::Surv(days, event) ~ risk_grp + sex + age + log_lymf_count,
  data = overall_model_df
)

overall_ph_gof <- mets::gof(overall_phreg_fit)
overall_age_linearity_gof <- mets::gofZ.phreg(survival::Surv(days, event) ~ age, data = overall_model_df)
overall_log_lymf_linearity_gof <- mets::gofZ.phreg(
  survival::Surv(days, event) ~ log_lymf_count,
  data = overall_model_df
)

analysis_helpers$save_text(capture.output(summary(overall_ph_gof)), "overall_ph_gof_summary.txt")
analysis_helpers$save_text(capture.output(summary(overall_age_linearity_gof)), "overall_age_linearity_gof_summary.txt")
analysis_helpers$save_text(capture.output(summary(overall_log_lymf_linearity_gof)), "overall_log_lymf_linearity_gof_summary.txt")

analysis_helpers$save_plot_png("overall_ph_gof.png", function() {
  graphics::par(mfrow = grDevices::n2mfrow(length(stats::coef(overall_phreg_fit))))
  plot(overall_ph_gof)
}, width = 2200, height = 1800)

analysis_helpers$save_plot_png("overall_age_linearity_gof.png", function() {
  plot(overall_age_linearity_gof, type = "z")
}, width = 1400, height = 1000)

analysis_helpers$save_plot_png("overall_log_lymf_linearity_gof.png", function() {
  plot(overall_log_lymf_linearity_gof, type = "z")
}, width = 1400, height = 1000)

# --- 2.2 to 2.4 collected outputs for the final report ---
# Keep a single object that points to the main analysis results so the analysis
# can be inspected at the end without searching through every intermediate name.
analysis_outputs <- list(
  cohort_summary = analysis_cohort_summary,
  baseline_by_risk_group = analysis_baseline_by_risk_group,
  overall = list(
    km_fit = overall_km_fit,
    logrank_test = overall_logrank_test,
    adjusted_cox_fit = overall_adjusted_cox_fit,
    event_free_summary = overall_event_free_summary,
    adjusted_cox_table = overall_adjusted_cox_table
  ),
  competing_risks = list(
    cif_fit = competing_risk_cif_fit,
    gray_fit = competing_risk_gray_fit,
    gray_tests = competing_risk_gray_tests,
    gray_timepoints = competing_risk_gray_timepoints,
    relapse_summary = competing_risk_relapse_summary,
    death_summary = competing_risk_death_summary,
    smn_summary = competing_risk_smn_summary,
    relapse_cox_fit = competing_risk_relapse_cox_fit,
    relapse_cox_table = competing_risk_relapse_cox_table,
    death_cox_fit = competing_risk_death_cox_fit,
    death_cox_table = competing_risk_death_cox_table,
    smn_cox_fit = competing_risk_smn_cox_fit,
    smn_cox_table = competing_risk_smn_cox_table
  ),
  landmark = list(
    counts = landmark_group_counts,
    sct_counts = landmark_sct_counts,
    km_fit = landmark_km_fit,
    adjusted_cox_fit = landmark_adjusted_cox_fit,
    event_free_summary = landmark_event_free_summary,
    adjusted_cox_table = landmark_adjusted_cox_table
  ),
  rmst = list(
    tau_days = rmst_tau_days,
    overall_rmst_table = overall_rmst_table,
    competing_risk_time_lost_table = competing_risk_time_lost_table
  ),
  diagnostics = list(
    phreg_fit = overall_phreg_fit,
    proportional_hazards_gof = overall_ph_gof,
    age_linearity_gof = overall_age_linearity_gof,
    log_lymf_linearity_gof = overall_log_lymf_linearity_gof
  )
)

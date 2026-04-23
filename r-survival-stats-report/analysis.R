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

# mets helper functions may evaluate formulas containing cluster() via the search
# path, so expose the survival version without attaching the whole package.
cluster <- survival::cluster

# --- 0.1 load the dataset ---
# The dataset is loaded a single time and then reused throughout the script.
pipeline_env <- new.env(parent = globalenv())
sys.source("load.R", envir = pipeline_env)
df <- pipeline_env$load_dataset()

# --- 0.2 and 0.4 clean and transform the analysis dataset ---
# Keep only the analysis-specific convenience variables that are used later.
analysis_df <- transform(
  df,
  sex = factor(sex),
  risk_grp = factor(risk_grp),
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
{
  analysis_cohort_summary <- data.frame(
    n_patients = nrow(analysis_df),
    n_events = sum(analysis_df$event == 1, na.rm = TRUE),
    n_censored = sum(analysis_df$event == 0, na.rm = TRUE),
    n_relapse = sum(analysis_df$comprisk == 1, na.rm = TRUE),
    n_death = sum(analysis_df$comprisk == 2, na.rm = TRUE)
  )

  analysis_risk_levels <- levels(analysis_df$risk_grp)
  analysis_baseline_by_risk_group <- data.frame(
    risk_grp = analysis_risk_levels,
    n = as.integer(table(analysis_df$risk_grp)[analysis_risk_levels]),
    n_events = as.numeric(tapply(analysis_df$event == 1, analysis_df$risk_grp, sum, na.rm = TRUE)[analysis_risk_levels]),
    n_relapse = as.numeric(tapply(analysis_df$comprisk == 1, analysis_df$risk_grp, sum, na.rm = TRUE)[analysis_risk_levels]),
    n_death = as.numeric(tapply(analysis_df$comprisk == 2, analysis_df$risk_grp, sum, na.rm = TRUE)[analysis_risk_levels])
  )
  analysis_baseline_dtable <- mets::dtable(analysis_df, ~ risk_grp + sex, level = 1)

  analysis_helpers$save_csv(analysis_cohort_summary, "cohort_summary.csv")
  analysis_helpers$save_csv(analysis_baseline_by_risk_group, "baseline_by_risk_group.csv")
  analysis_helpers$save_text(capture.output(analysis_baseline_dtable), "baseline_dtable.txt")
}

# --- 1.3 and 1.6 overall event-free survival and significance testing ---
# This follows the standard course workflow: non-parametric curves first, then
# log-rank testing, then an adjusted Cox model.
{
  overall_surv_df <- subset(
    analysis_df,
    !is.na(days) & !is.na(event) & !is.na(risk_grp)
  )

  overall_model_df <- subset(
    overall_surv_df,
    !is.na(sex) & !is.na(age) & !is.na(log_lymf_count)
  )

  overall_report_times_days <- report_times_days[
    report_times_days <= max(overall_surv_df$days, na.rm = TRUE)
  ]

  overall_surv_fit <- survival::survfit(
    survival::Surv(days, event) ~ risk_grp,
    data = overall_surv_df
  )

  overall_logrank_test <- survival::survdiff(
    survival::Surv(days, event) ~ risk_grp,
    data = overall_surv_df
  )

  overall_adjusted_cox_fit <- survival::coxph(
    survival::Surv(days, event) ~ risk_grp + sex + age + log_lymf_count,
    data = overall_model_df,
    x = TRUE
  )

  overall_event_free_summary <- analysis_helpers$extract_surv_summary(
    overall_surv_fit,
    overall_report_times_days
  )
  overall_adjusted_cox_table <- analysis_helpers$extract_cox_table(overall_adjusted_cox_fit)

  analysis_helpers$save_csv(overall_event_free_summary, "overall_event_free_summary.csv")
  analysis_helpers$save_csv(overall_adjusted_cox_table, "overall_adjusted_cox_table.csv")
  analysis_helpers$save_text(capture.output(overall_logrank_test), "overall_logrank_test.txt")
  analysis_helpers$save_text(capture.output(summary(overall_adjusted_cox_fit)), "overall_adjusted_cox_summary.txt")
}

# --- 2.1 stratified hazard curves ---
{
  overall_surv_df <- subset(
    analysis_df,
    !is.na(days) & !is.na(event) & !is.na(risk_grp)
  )
  overall_surv_fit <- survival::survfit(
    survival::Surv(days, event) ~ risk_grp,
    data = overall_surv_df
  )
  group_levels <- levels(overall_surv_df$risk_grp)
  group_cols <- seq_along(group_levels)

  analysis_helpers$save_plot_with_side_legend_png("figure_01_overall_event_free_survival_by_risk_group.png", function() {
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

    graphics::plot(
      overall_surv_fit,
      col = group_cols,
      lty = 1,
      xlab = "Days since diagnosis",
      ylab = "Event-free survival probability",
      main = "Event-Free Survival by Risk Group After Diagnosis"
    )
  }, function() {
    graphics::legend(
      "center",
      legend = group_levels,
      col = group_cols,
      lty = 1,
      bty = "n",
      title = "Risk group"
    )
  })

  analysis_helpers$save_plot_with_side_legend_png("figure_02_overall_cumulative_hazard_by_risk_group.png", function() {
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

    graphics::plot(
      overall_surv_fit,
      fun = "cumhaz",
      col = group_cols,
      lty = 1,
      xlab = "Days since diagnosis",
      ylab = "Cumulative hazard",
      main = "Cumulative Event Hazard by Risk Group After Diagnosis"
    )
  }, function() {
    graphics::legend(
      "center",
      legend = group_levels,
      col = group_cols,
      lty = 1,
      bty = "n",
      title = "Risk group"
    )
  })
}

# --- 1.1 competing risks ---
# This mirrors the course examples: non-parametric cumulative incidence curves,
# Gray's test, and then cause-specific Cox models for the main event types.
# SMN is grouped with relapse in cause 1.
{
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

  competing_risk_relapse_summary <- analysis_helpers$extract_prodlim_summary(
    competing_risk_cif_fit,
    competing_risk_report_times_days,
    cause = 1
  )
  competing_risk_death_summary <- analysis_helpers$extract_prodlim_summary(
    competing_risk_cif_fit,
    competing_risk_report_times_days,
    cause = 2
  )

  competing_risk_relapse_cox_fit <- analysis_helpers$fit_cause_specific_cox(analysis_df, 1)
  competing_risk_death_cox_fit <- analysis_helpers$fit_cause_specific_cox(analysis_df, 2)

  competing_risk_relapse_cox_table <- analysis_helpers$extract_cox_table(competing_risk_relapse_cox_fit)
  competing_risk_death_cox_table <- analysis_helpers$extract_cox_table(competing_risk_death_cox_fit)

  analysis_helpers$save_csv(competing_risk_gray_tests, "competing_risk_gray_tests.csv")
  analysis_helpers$save_csv(competing_risk_gray_timepoints, "competing_risk_gray_timepoints.csv")
  analysis_helpers$save_csv(competing_risk_relapse_cox_table, "competing_risk_relapse_cox_table.csv")
  analysis_helpers$save_text(capture.output(summary(competing_risk_relapse_cox_fit)), "competing_risk_relapse_cox_summary.txt")
  analysis_helpers$save_csv(competing_risk_death_cox_table, "competing_risk_death_cox_table.csv")
  analysis_helpers$save_text(capture.output(summary(competing_risk_death_cox_fit)), "competing_risk_death_cox_summary.txt")
}

{
  competing_risk_df <- subset(
    analysis_df,
    !is.na(compriskdays) & !is.na(comprisk) & !is.na(risk_grp)
  )
  competing_risk_cif_fit <- prodlim::prodlim(
    prodlim::Hist(compriskdays, comprisk) ~ risk_grp,
    data = competing_risk_df
  )
  group_levels <- levels(competing_risk_df$risk_grp)
  group_cols <- seq_along(group_levels)

  analysis_helpers$save_plot_with_side_legend_png("figure_03_competing_risk_cif_relapse_or_smn_by_risk_group.png", function() {
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
    plot(
      competing_risk_cif_fit,
      cause = 1,
      legend = FALSE,
      xlab = "Days since diagnosis",
      ylab = "Cumulative incidence"
    )
    graphics::title(main = "Cumulative Incidence of Relapse or SMN by Risk Group")
  }, function() {
    graphics::legend(
      "center",
      legend = group_levels,
      lty = 1,
      col = group_cols,
      bty = "n",
      title = "Risk group"
    )
  })

  analysis_helpers$save_plot_with_side_legend_png("figure_04_competing_risk_cif_death_by_risk_group.png", function() {
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))
    plot(
      competing_risk_cif_fit,
      cause = 2,
      legend = FALSE,
      xlab = "Days since diagnosis",
      ylab = "Cumulative incidence"
    )
    graphics::title(main = "Cumulative Incidence of Death by Risk Group")
  }, function() {
    graphics::legend(
      "center",
      legend = group_levels,
      lty = 1,
      col = group_cols,
      bty = "n",
      title = "Risk group"
    )
  })
}

# --- 1.2 landmark analysis for SCT ---
# Patients are classified by SCT status at the landmark date and then followed
# from the landmark onward. This avoids using future SCT information at baseline.
{
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

  landmark_km_fit <- survival::survfit(
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

  landmark_event_free_summary <- analysis_helpers$extract_surv_summary(
    landmark_km_fit,
    landmark_report_times_days
  )
  landmark_adjusted_cox_table <- analysis_helpers$extract_cox_table(landmark_adjusted_cox_fit)

  analysis_helpers$save_csv(landmark_group_counts, "landmark_group_counts.csv")
  analysis_helpers$save_csv(landmark_sct_counts, "landmark_sct_counts.csv")
  analysis_helpers$save_csv(landmark_event_free_summary, "landmark_event_free_summary.csv")
  analysis_helpers$save_csv(landmark_adjusted_cox_table, "landmark_adjusted_cox_table.csv")
  analysis_helpers$save_text(capture.output(summary(landmark_adjusted_cox_fit)), "landmark_adjusted_cox_summary.txt")
}

{
  landmark_source_df <- subset(
    analysis_df,
    !is.na(landmark_date) & !is.na(date) & !is.na(event) & !is.na(sct_by_landmark)
  )
  landmark_analysis_df <- transform(
    subset(landmark_source_df, date > landmark_date),
    time_from_landmark = as.numeric(date - landmark_date),
    event_from_landmark = event
  )
  landmark_km_fit <- survival::survfit(
    survival::Surv(time_from_landmark, event_from_landmark) ~ sct_by_landmark,
    data = landmark_analysis_df
  )
  sct_levels <- levels(landmark_analysis_df$sct_by_landmark)
  sct_cols <- seq_along(sct_levels)

  analysis_helpers$save_plot_with_side_legend_png("figure_05_landmark_event_free_survival_by_sct_status.png", function() {
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

    graphics::plot(
      landmark_km_fit,
      col = sct_cols,
      lty = 1,
      xlab = "Days since landmark",
      ylab = "Post-landmark event-free survival",
      main = "Post-Landmark Event-Free Survival by SCT Status"
    )
  }, function() {
    graphics::legend(
      "center",
      legend = sct_levels,
      col = sct_cols,
      lty = 1,
      bty = "n",
      title = "SCT status"
    )
  })
}

# --- 1.5 restricted mean survival time ---
# RMST provides an absolute time-scale summary that complements the hazard ratio.
{
  rmst_overall_surv_df <- subset(
    analysis_df,
    !is.na(days) & !is.na(event) & !is.na(risk_grp)
  )
  rmst_overall_surv_fit <- survival::survfit(
    survival::Surv(days, event) ~ risk_grp,
    data = rmst_overall_surv_df
  )
  rmst_competing_risk_df <- subset(
    analysis_df,
    !is.na(compriskdays) & !is.na(comprisk) & !is.na(risk_grp)
  )

  rmst_tau_days <- min(5 * 365.25, max(rmst_overall_surv_df$days, na.rm = TRUE))

  overall_rmst_table <- analysis_helpers$extract_rmst_table(rmst_overall_surv_fit, rmst_tau_days)
  overall_rmtl_table <- transform(
    overall_rmst_table,
    restricted_mean_time_lost = days_lost
  )
  competing_risk_time_lost_relapse <- analysis_helpers$extract_time_lost_by_cause(
    rmst_competing_risk_df,
    rmst_tau_days,
    1,
    "Relapse or SMN"
  )
  competing_risk_time_lost_death <- analysis_helpers$extract_time_lost_by_cause(
    rmst_competing_risk_df,
    rmst_tau_days,
    2,
    "Death"
  )
  competing_risk_time_lost_table <- rbind(
    competing_risk_time_lost_relapse,
    competing_risk_time_lost_death
  )

  analysis_helpers$save_csv(overall_rmst_table, "overall_rmst_table.csv")
  analysis_helpers$save_csv(overall_rmtl_table, "overall_rmtl_table.csv")
  analysis_helpers$save_csv(competing_risk_time_lost_table, "competing_risk_time_lost_table.csv")
}

{
  rmst_overall_surv_df <- subset(
    analysis_df,
    !is.na(days) & !is.na(event) & !is.na(risk_grp)
  )
  rmst_overall_surv_fit <- survival::survfit(
    survival::Surv(days, event) ~ risk_grp,
    data = rmst_overall_surv_df
  )
  rmst_tau_days <- min(5 * 365.25, max(rmst_overall_surv_df$days, na.rm = TRUE))
  group_levels <- levels(rmst_overall_surv_df$risk_grp)
  group_cols <- seq_along(group_levels)

  analysis_helpers$save_plot_with_side_legend_png("figure_06_overall_rmst_survival_by_risk_group.png", function() {
    graphics::par(mar = c(5.1, 4.1, 4.1, 2.1))

    graphics::plot(
      rmst_overall_surv_fit,
      col = group_cols,
      lty = 1,
      xlab = "Days since diagnosis",
      ylab = "Event-free survival probability",
      xlim = c(0, rmst_tau_days),
      main = "Event-Free Survival by Risk Group Up to the RMST Horizon"
    )
    graphics::abline(v = rmst_tau_days, lty = 2)
  }, function() {
    graphics::legend(
      "center",
      legend = group_levels,
      col = group_cols,
      lty = 1,
      bty = "n",
      title = "Risk group"
    )
  })

  overall_rmst_table <- analysis_helpers$extract_rmst_table(rmst_overall_surv_fit, rmst_tau_days)
  overall_rmtl_values <- overall_rmst_table$days_lost
  names(overall_rmtl_values) <- overall_rmst_table$group

  analysis_helpers$save_plot_png("figure_07_overall_restricted_mean_time_lost_by_risk_group.png", function() {
    graphics::par(mar = c(8.1, 4.1, 4.1, 2.1))
    bar_cols <- seq_along(overall_rmtl_values)
    bar_positions <- graphics::barplot(
      overall_rmtl_values,
      col = bar_cols,
      las = 2,
      ylab = "Days lost up to tau",
      main = "Restricted Mean Time Lost by Risk Group"
    )
    graphics::abline(h = 0, lty = 1)
    graphics::text(
      x = bar_positions,
      y = overall_rmtl_values,
      labels = round(overall_rmtl_values, 1),
      pos = 3,
      cex = 0.8
    )
  }, width = 1800, height = 1200)

  rmst_competing_risk_df <- subset(
    analysis_df,
    !is.na(compriskdays) & !is.na(comprisk) & !is.na(risk_grp)
  )
  competing_risk_time_lost_relapse <- analysis_helpers$extract_time_lost_by_cause(
    rmst_competing_risk_df,
    rmst_tau_days,
    1,
    "Relapse or SMN"
  )
  competing_risk_time_lost_death <- analysis_helpers$extract_time_lost_by_cause(
    rmst_competing_risk_df,
    rmst_tau_days,
    2,
    "Death"
  )
  competing_risk_time_lost_table <- rbind(
    competing_risk_time_lost_relapse,
    competing_risk_time_lost_death
  )
  competing_risk_rmtl_matrix <- xtabs(
    restricted_mean_time_lost ~ cause + group,
    data = competing_risk_time_lost_table
  )
  bar_cols <- c("steelblue", "firebrick")

  analysis_helpers$save_plot_with_side_legend_png("figure_08_competing_risk_restricted_mean_time_lost_by_group_and_cause.png", function() {
    graphics::par(mar = c(8.1, 4.1, 4.1, 2.1))
    graphics::barplot(
      competing_risk_rmtl_matrix,
      beside = TRUE,
      col = bar_cols,
      las = 2,
      ylab = "Days lost up to tau",
      main = "Restricted Mean Time Lost by Risk Group and Cause"
    )
  }, function() {
    graphics::legend(
      "center",
      legend = rownames(competing_risk_rmtl_matrix),
      fill = bar_cols,
      bty = "n",
      title = "Cause"
    )
  }, width = 2200, height = 1200)
}

# --- 1.4 goodness-of-fit and model diagnostics ---
# The full Cox model is checked for proportional hazards, and the continuous
# covariates are checked for linearity on the log-hazard scale.
{
  diagnostics_model_df <- subset(
    analysis_df,
    !is.na(days) & !is.na(event) & !is.na(risk_grp) & !is.na(sex) & !is.na(age) & !is.na(log_lymf_count)
  )

  overall_phreg_fit <- mets::phreg(
    survival::Surv(days, event) ~ risk_grp + sex + age + log_lymf_count,
    data = diagnostics_model_df
  )

  overall_ph_gof <- mets::gof(overall_phreg_fit)
  overall_age_linearity_gof <- mets::gofZ.phreg(survival::Surv(days, event) ~ age, data = diagnostics_model_df)
  overall_log_lymf_linearity_gof <- mets::gofZ.phreg(
    survival::Surv(days, event) ~ log_lymf_count,
    data = diagnostics_model_df
  )

  analysis_helpers$save_text(capture.output(summary(overall_ph_gof)), "overall_ph_gof_summary.txt")
  analysis_helpers$save_text(capture.output(summary(overall_age_linearity_gof)), "overall_age_linearity_gof_summary.txt")
  analysis_helpers$save_text(capture.output(summary(overall_log_lymf_linearity_gof)), "overall_log_lymf_linearity_gof_summary.txt")
}

{
  diagnostics_model_df <- subset(
    analysis_df,
    !is.na(days) & !is.na(event) & !is.na(risk_grp) & !is.na(sex) & !is.na(age) & !is.na(log_lymf_count)
  )
  overall_phreg_fit <- mets::phreg(
    survival::Surv(days, event) ~ risk_grp + sex + age + log_lymf_count,
    data = diagnostics_model_df
  )
  overall_ph_gof <- mets::gof(overall_phreg_fit)
  overall_age_linearity_gof <- mets::gofZ.phreg(survival::Surv(days, event) ~ age, data = diagnostics_model_df)
  overall_log_lymf_linearity_gof <- mets::gofZ.phreg(
    survival::Surv(days, event) ~ log_lymf_count,
    data = diagnostics_model_df
  )

  analysis_helpers$save_plot_png("figure_09_overall_ph_gof.png", function() {
    graphics::par(mfrow = grDevices::n2mfrow(length(stats::coef(overall_phreg_fit))))
    plot(overall_ph_gof)
  }, width = 2200, height = 1800)

  analysis_helpers$save_plot_png("figure_10_overall_age_linearity_gof.png", function() {
    plot(overall_age_linearity_gof, type = "z")
  }, width = 1400, height = 1000)

  analysis_helpers$save_plot_png("figure_11_overall_log_lymf_linearity_gof.png", function() {
    plot(overall_log_lymf_linearity_gof, type = "z")
  }, width = 1400, height = 1000)
}

# --- 2.2 to 2.4 collected outputs for the final report ---
# Keep a single object that points to the main analysis results so the analysis
# can be inspected at the end without searching through every intermediate name.
analysis_outputs <- list(
  cohort_summary = analysis_cohort_summary,
  baseline_by_risk_group = analysis_baseline_by_risk_group,
  baseline_dtable = analysis_baseline_dtable,
  overall = list(
    surv_fit = overall_surv_fit,
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
    relapse_cox_fit = competing_risk_relapse_cox_fit,
    relapse_cox_table = competing_risk_relapse_cox_table,
    death_cox_fit = competing_risk_death_cox_fit,
    death_cox_table = competing_risk_death_cox_table
  ),
  landmark = list(
    counts = landmark_group_counts,
    sct_counts = landmark_sct_counts,
    surv_fit = landmark_km_fit,
    adjusted_cox_fit = landmark_adjusted_cox_fit,
    event_free_summary = landmark_event_free_summary,
    adjusted_cox_table = landmark_adjusted_cox_table
  ),
  rmst = list(
    tau_days = rmst_tau_days,
    overall_rmst_table = overall_rmst_table,
    overall_rmtl_table = overall_rmtl_table,
    competing_risk_time_lost_relapse = competing_risk_time_lost_relapse,
    competing_risk_time_lost_death = competing_risk_time_lost_death,
    competing_risk_time_lost_table = competing_risk_time_lost_table
  ),
  diagnostics = list(
    phreg_fit = overall_phreg_fit,
    proportional_hazards_gof = overall_ph_gof,
    age_linearity_gof = overall_age_linearity_gof,
    log_lymf_linearity_gof = overall_log_lymf_linearity_gof
  )
)


save_plot_png <- function(file_name, plot_code, width = 1600, height = 1200, res = 180) {
  grDevices::png(
    filename = file.path(out_dir, file_name),
    width = width,
    height = height,
    res = res
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  plot_code()
  invisible(file.path(out_dir, file_name))
}

save_csv <- function(x, file_name) {
  utils::write.csv(x, file.path(out_dir, file_name), row.names = FALSE)
}

save_text <- function(lines, file_name) {
  writeLines(lines, file.path(out_dir, file_name))
}

extract_prodlim_summary <- function(fit, times, cause = NULL) {
  fit_summary <- if (is.null(cause)) {
    summary(fit, times = times)
  } else {
    summary(fit, times = times, cause = cause)
  }

  as.data.frame(fit_summary)
}

extract_cuminc_timepoints <- function(fit, times) {
  fit_timepoints <- cmprsk::timepoints(fit, times = times)
  est_matrix <- fit_timepoints$est
  var_matrix <- fit_timepoints$var
  curves <- rownames(est_matrix)
  time_values <- as.numeric(colnames(est_matrix))

  data.frame(
    curve = rep(curves, each = length(time_values)),
    time = rep(time_values, times = length(curves)),
    estimate = as.vector(t(est_matrix)),
    variance = as.vector(t(var_matrix)),
    std_err = sqrt(as.vector(t(var_matrix)))
  )
}

extract_cox_table <- function(fit) {
  fit_summary <- summary(fit)
  coef_table <- fit_summary$coefficients
  conf_table <- fit_summary$conf.int

  data.frame(
    term = rownames(coef_table),
    hazard_ratio = conf_table[, "exp(coef)"],
    conf_low = conf_table[, "lower .95"],
    conf_high = conf_table[, "upper .95"],
    p_value = coef_table[, "Pr(>|z|)"],
    row.names = NULL
  )
}

extract_rmst_table <- function(fit, tau_days) {
  rmst_summary <- summary(fit, rmean = tau_days)$table

  if (is.null(dim(rmst_summary))) {
    rmst_summary <- t(rmst_summary)
  }

  rmst_df <- as.data.frame(rmst_summary)
  rmst_df$group <- rownames(rmst_df)
  rownames(rmst_df) <- NULL

  data.frame(
    group = rmst_df$group,
    rmst_days = rmst_df$rmean,
    rmst_std_err = rmst_df$`se(rmean)`,
    days_lost = tau_days - rmst_df$rmean,
    tau_days = tau_days
  )
}

fit_cause_specific_cox <- function(data, cause_code) {
  cause_df <- subset(
    data,
    !is.na(compriskdays) &
      !is.na(comprisk) &
      !is.na(risk_grp) &
      !is.na(sex) &
      !is.na(age) &
      !is.na(log_lymf_count)
  )

  survival::coxph(
    survival::Surv(compriskdays, comprisk == cause_code) ~ risk_grp + sex + age + log_lymf_count,
    data = cause_df,
    x = TRUE
  )
}

extract_time_lost_by_cause <- function(data, tau_days, cause_code, cause_label) {
  fit <- mets::resmeanIPCW(
    mets::Event(compriskdays, comprisk) ~ -1 + risk_grp,
    data = data,
    time = tau_days,
    cause = cause_code,
    cens.model = ~ risk_grp,
    model = "lin"
  )

  fit_estimate <- mets::estimate(fit)

  data.frame(
    group = names(stats::coef(fit_estimate)),
    cause = cause_label,
    restricted_mean_time_lost = as.numeric(stats::coef(fit_estimate)),
    tau_days = tau_days,
    row.names = NULL
  )
}

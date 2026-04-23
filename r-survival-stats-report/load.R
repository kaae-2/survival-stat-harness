load_dataset <- function() {
  # avoiding poluting working env
  helpers = new.env(parent = baseenv())
  sys.source("load_helpers.R", envir = helpers, chdir = TRUE)
  administrative_censor_date <- as.Date("2026-03-01")

  # read and join the 2 data sets
  readRenviron(".Renviron")

  data_df <- haven::read_sav(helpers$get_env_required("SAV_DATA_PATH"))
  date_df <- haven::read_sav(helpers$get_env_required("SAV_DATE_PATH"))

  helpers$write_source_metadata(data_df, "data_df")
  helpers$write_source_metadata(date_df, "date_df")

  header_reference <- data.frame(
    dataset = c(
      rep("data_df", length(names(data_df))),
      rep("date_df", length(names(date_df)))
    ),
    header = c(names(data_df), names(date_df))
  )

  if (tolower(Sys.getenv("WRITE_HEADER_REFERENCE", unset = "0")) %in% c("1", "true", "yes")) {
    dir.create("out", showWarnings = FALSE, recursive = TRUE)

    # Persist raw imported headers for later column mapping reference.
    utils::write.csv(
      header_reference,
      file.path("out", "data_headers_reference.csv"),
      row.names = FALSE
    )
  }

  data_df <- helpers$rename_from_env(
    data_df,
    c(
      id = "COL_ID",
      included = "COL_INCLUDED",
      sex = "COL_SEX",
      diagnosis_date = "COL_DATE_DIAGNOSIS",
      age = "COL_AGE",
      risk_grp = "COL_RISK_GROUP",
      lymf_count = "COL_LYMF_COUNT",
      delayed_exclusion_date = "COL_DELAYED_EXCLUSION_DATE"
    ),
    dataset_name = "data_df",
    header_reference = header_reference
  )

  data_df <- data_df[
    !is.na(data_df$included) &
      data_df$included == 1 &
      (is.na(data_df$risk_grp) | data_df$risk_grp != 0),
  ]

  date_df <- helpers$rename_from_env(
    date_df,
    c(
      id = "COL_ID",
      primary_event_status = "COL_PRIMARY_EVENT_STATUS",
      primary_event_date = "COL_PRIMARY_EVENT_DATE",
      followup_active_date = "COL_FOLLOWUP_ACTIVE_DATE",
      last_followup_date = "COL_LAST_FOLLOWUP_DATE",
      lost_to_followup = "COL_LOST_TO_FOLLOWUP",
      lost_to_followup_date = "COL_LOST_TO_FOLLOWUP_DATE",
      landmark_date = "COL_LANDMARK_DATE",
      trsct_date = "COL_TRSCT_DATE",
      bmt1_date = "COL_BMT1_DATE"
    ),
    dataset_name = "date_df",
    header_reference = header_reference
  )

  df <- dplyr::left_join(data_df, date_df, by = "id")

  # Build only the final survival-analysis fields used downstream.
  df <- dplyr::mutate(
    df,
    raw_sct_date = dplyr::coalesce(bmt1_date, trsct_date),
    raw_followdate = dplyr::case_when(
      lost_to_followup == 1 ~ lost_to_followup_date,
      !is.na(delayed_exclusion_date) ~ delayed_exclusion_date,
      followup_active_date > last_followup_date ~ followup_active_date,
      !is.na(last_followup_date) ~ last_followup_date,
      !is.na(followup_active_date) ~ followup_active_date,
      TRUE ~ administrative_censor_date
    ),
    followdate = dplyr::case_when(
      is.na(raw_followdate) ~ administrative_censor_date,
      raw_followdate > administrative_censor_date ~ administrative_censor_date,
      TRUE ~ raw_followdate
    ),
    primary_event_days = dplyr::if_else(
      is.na(primary_event_date),
      NA_real_,
      as.numeric(primary_event_date - diagnosis_date)
    ),
    event = dplyr::case_when(
      !is.na(primary_event_status) &
        primary_event_status > 0 &
        !is.na(primary_event_date) &
        primary_event_date <= administrative_censor_date ~ 1,
      !is.na(followdate) ~ 0,
      TRUE ~ NA_real_
    ),
    relapse = dplyr::case_when(
      !is.na(primary_event_status) &
        primary_event_status %in% c(3, 5) &
        !is.na(primary_event_date) &
        primary_event_date <= administrative_censor_date ~ 1,
      TRUE ~ 0
    ),
    comprisk = dplyr::case_when(
      !is.na(primary_event_status) &
        primary_event_status %in% c(3, 5) &
        !is.na(primary_event_date) &
        primary_event_date <= administrative_censor_date ~ 1,
      !is.na(primary_event_status) &
        primary_event_status %in% c(1, 2, 4) &
        !is.na(primary_event_date) &
        primary_event_date <= administrative_censor_date ~ 2,
      !is.na(followdate) ~ 0,
      TRUE ~ NA_real_
    ),
    sct_date = dplyr::if_else(
      !is.na(raw_sct_date) &
        primary_event_status %in% c(3, 5) &
        !is.na(primary_event_date) &
        primary_event_date <= administrative_censor_date &
        primary_event_days > 0 &
        as.numeric(raw_sct_date - diagnosis_date) > primary_event_days,
      as.Date(NA),
      raw_sct_date,
      missing = raw_sct_date
    ),
    timetosct = dplyr::if_else(
      is.na(sct_date),
      NA_real_,
      as.numeric(sct_date - diagnosis_date)
    ),
    days = dplyr::if_else(
      event == 1,
      primary_event_days,
      as.numeric(followdate - diagnosis_date),
      missing = NA_real_
    ),
    # totaldays currently follows the same time scale as days.
    totaldays = days,
    date = dplyr::if_else(
      event == 1,
      primary_event_date,
      followdate,
      missing = as.Date(NA)
    ),
    compriskdays = dplyr::case_when(
      comprisk == 0 ~ as.numeric(followdate - diagnosis_date),
      comprisk %in% c(1, 2) ~ primary_event_days
    )
  )

  dplyr::select(
    df,
    id,
    event,
    comprisk,
    sex,
    age,
    risk_grp,
    lymf_count,
    days,
    totaldays,
    compriskdays,
    date,
    relapse,
    timetosct,
    landmark_date,
    sct_date
  )
}

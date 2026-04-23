load_dataset <- function() {
  # avoiding poluting working env
  helpers = new.env(parent = baseenv())
  sys.source("load_helpers.R", envir=helpers, chdir=TRUE)

  #read and join the 2 data sets:
  readRenviron(".Renviron")
  
  data_df <- haven::read_sav(helpers$get_env_required("SAV_DATA_PATH"))
  date_df <- haven::read_sav(helpers$get_env_required("SAV_DATE_PATH"))

  # Persist raw imported headers for later column mapping reference.
  utils::write.csv(
    data.frame(
      dataset = c(
        rep("data_df", length(names(data_df))),
        rep("date_df", length(names(date_df)))
      ),
      header = c(names(data_df), names(date_df))
    ),
    "data_headers_reference.csv",
    row.names = FALSE
  )
  header_reference_path <- "data_headers_reference.csv"
  
  data_df <- helpers$rename_from_env(data_df, c(
    id = "COL_ID",
    included="COL_INCLUDED",
    type="COL_TYPE",
    sex="COL_SEX",
    diagnosis_date="COL_DATE_DIAGNOSIS",
    age="COL_AGE",
    risk_grp="COL_RISK_GROUP",
    date_death="COL_DEATH_DATE",
    lymf_count="COL_LYMF_COUNT",
    delayed_exclusion_date="COL_DELAYED_EXCLUSION_DATE"
  ), dataset_name = "data_df", reference_path = header_reference_path)
  
  data_df <- data_df[!is.na(data_df$included) & data_df$included == 1,]
  
  date_df <- helpers$rename_from_env(date_df, c(
    id = "COL_ID",
    primary_event_status="COL_PRIMARY_EVENT_STATUS",
    primary_event_date="COL_PRIMARY_EVENT_DATE",
    followup_active_date="COL_FOLLOWUP_ACTIVE_DATE",
    last_followup_date="COL_LAST_FOLLOWUP_DATE",
    lost_to_followup="COL_LOST_TO_FOLLOWUP",
    lost_to_followup_date="COL_LOST_TO_FOLLOWUP_DATE",
    landmark_date="COL_LANDMARK_DATE",
    trsct_date="COL_TRSCT_DATE",
    bmt1_date="COL_BMT1_DATE"
  ), dataset_name = "date_df", reference_path = header_reference_path)
  
  df <- dplyr::left_join(data_df, date_df, by = "id")
  
  
  #transforming data to a surv compatible object:
  df <- df %>%
    dplyr::mutate(
      sct_date = dplyr::coalesce(bmt1_date, trsct_date),
      timetosct = dplyr::if_else(
        is.na(sct_date),
        NA_real_,
        as.numeric(sct_date - diagnosis_date)
      ),
      event = dplyr::case_when(
        is.na(primary_event_status) ~ NA_real_,
        primary_event_status > 0 ~ 1,
        primary_event_status < 1 ~ 0
      ),
      relapse = dplyr::case_when(
        is.na(primary_event_status) ~ NA_real_,
        primary_event_status == 3 ~ 1,
        TRUE ~ 0
      ),
      followdate = dplyr::case_when(
        lost_to_followup == 1 ~ lost_to_followup_date,
        !is.na(delayed_exclusion_date) ~ delayed_exclusion_date,
        is.na(followup_active_date) ~ last_followup_date,
        is.na(last_followup_date) ~ followup_active_date,
        followup_active_date > last_followup_date ~ followup_active_date,
        TRUE ~ last_followup_date
      ),
      followtotaldate = followdate,
      followtotaldays = as.numeric(followtotaldate - diagnosis_date),
      followdays = as.numeric(followdate - diagnosis_date),
      timetoevent = dplyr::if_else(
        event == 1,
        as.numeric(primary_event_date - diagnosis_date),
        NA_real_
      ),
      totaldays = dplyr::if_else(event == 1, timetoevent, followtotaldays, missing = NA_real_),
      days = dplyr::if_else(event == 1, timetoevent, followdays, missing = NA_real_),
      date = dplyr::if_else(event == 1, primary_event_date, followdate, missing = as.Date(NA)),
      comprisk = dplyr::case_when(
        is.na(primary_event_status) ~ NA_real_,
        primary_event_status < 1 ~ 0,
        primary_event_status == 3 ~ 1,
        primary_event_status %in% c(1, 2, 4) ~ 2,
        primary_event_status == 5 ~ 3
      ),
      timetorelapse = dplyr::if_else(
        comprisk == 1,
        as.numeric(primary_event_date - diagnosis_date),
        NA_real_
      ),
      timetodeath = dplyr::if_else(
        comprisk == 2,
        as.numeric(primary_event_date - diagnosis_date),
        NA_real_
      ),
      timetosmn = dplyr::if_else(
        comprisk == 3,
        as.numeric(primary_event_date - diagnosis_date),
        NA_real_
      ),
      compriskdays = dplyr::case_when(
        comprisk == 0 ~ followdays,
        comprisk == 1 ~ timetorelapse,
        comprisk == 2 ~ timetodeath,
        comprisk == 3 ~ timetosmn
      ),
      compriskdate = dplyr::if_else(
        comprisk > 0,
        primary_event_date,
        followdate,
        missing = as.Date(NA)
      ),
      sct_date = dplyr::if_else(
        !is.na(timetorelapse) & timetosct > timetorelapse & timetorelapse > 0,
        as.Date(NA),
        sct_date,
        missing = sct_date
      )
    )
  
  # Keep only the derived survival-analysis fields used downstream.
  df <- df %>%
    dplyr::select(
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
  
  df
}

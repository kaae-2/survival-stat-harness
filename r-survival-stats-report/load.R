load_dataset <- function() {
  helpers = new.env(parent = baseenv())
  sys.source("helpers.R", envir=helpers, chdir=TRUE)
  readRenviron(".Renviron")
  
  df <- haven::read_sav(helpers$get_env_required("SAV_DATA_PATH"))
  date_df <- haven::read_sav(helpers$get_env_required("SAV_DATE_PATH"))
  
  df <- helpers$rename_from_env(df, c(
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
  
  date_df <- helpers$rename_from_env(date_df, c(
    id = "COL_ID",
    landmark_date="COL_LANDMARK_DATE",
    trsct_date="COL_TRSCT_DATE"
  ))
  
  df <- dplyr::left_join(df, date_df, by = "id")
  
  df
}
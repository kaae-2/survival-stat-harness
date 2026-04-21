# renaming columns from env files to avoid leaking information of the dataset

get_env_required <- function(x) {
  value <- Sys.getenv(x)
  if (value == "") stop(sprintf("Missing env var: %s", x))
  value
}

rename_from_env <- function(df, mapping) {
  raw_names <- vapply(mapping, get_env_required, character(1))
  
  missing_cols <- raw_names[!raw_names %in% names(df)]
  if (length(missing_cols) > 0) {
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
  }
  
  names(df)[match(raw_names, names(df))] <- names(mapping)
  df <- df[names(mapping)]
  df
}
# renaming columns from env files to avoid leaking information of the dataset

get_env_required <- function(x) {
  value <- Sys.getenv(x)
  if (value == "") stop(sprintf("Missing env var: %s", x))
  value
}

rename_from_env <- function(df, mapping, dataset_name = NULL, reference_path = NULL) {
  raw_names <- vapply(mapping, get_env_required, character(1))

  if (!is.null(dataset_name) && !is.null(reference_path) && file.exists(reference_path)) {
    header_reference <- utils::read.csv(reference_path, stringsAsFactors = FALSE)
    allowed_names <- header_reference$header[header_reference$dataset == dataset_name]
    wrong_dataset_cols <- raw_names[!raw_names %in% allowed_names]

    if (length(wrong_dataset_cols) > 0) {
      stop(
        "Columns mapped to ", dataset_name, " are not listed there in ",
        reference_path, ": ",
        paste(wrong_dataset_cols, collapse = ", ")
      )
    }
  }
  
  missing_cols <- raw_names[!raw_names %in% names(df)]
  if (length(missing_cols) > 0) {
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
  }
  
  names(df)[match(raw_names, names(df))] <- names(mapping)
  df <- df[names(mapping)]
  df
}

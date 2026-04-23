# renaming columns from env files to avoid leaking information of the dataset

get_env_required <- function(x) {
  value <- Sys.getenv(x)
  if (value == "") stop(sprintf("Missing env var: %s", x))
  value
}

rename_from_env <- function(df, mapping, dataset_name = NULL, header_reference = NULL) {
  raw_names <- vapply(mapping, get_env_required, character(1))

  if (!is.null(dataset_name) && !is.null(header_reference)) {
    allowed_names <- header_reference$header[header_reference$dataset == dataset_name]
    wrong_dataset_cols <- raw_names[!raw_names %in% allowed_names]

    if (length(wrong_dataset_cols) > 0) {
      stop(
        "Columns mapped to ", dataset_name, " are not listed there in the header reference: ",
        paste(wrong_dataset_cols, collapse = ", ")
      )
    }
  }

  missing_cols <- raw_names[!raw_names %in% names(df)]
  if (length(missing_cols) > 0) {
    stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
  }

  names(df)[match(raw_names, names(df))] <- names(mapping)
  df[names(mapping)]
}

extract_variable_metadata <- function(df, dataset_name) {
  data.frame(
    dataset = dataset_name,
    variable = names(df),
    class = vapply(df, function(x) paste(class(x), collapse = "|"), character(1)),
    label = vapply(df, function(x) {
      value <- attr(x, "label", exact = TRUE)
      if (is.null(value)) "" else as.character(value)
    }, character(1)),
    format_spss = vapply(df, function(x) {
      value <- attr(x, "format.spss", exact = TRUE)
      if (is.null(value)) "" else as.character(value)
    }, character(1)),
    display_width = vapply(df, function(x) {
      value <- attr(x, "display_width", exact = TRUE)
      if (is.null(value)) NA_real_ else as.numeric(value)
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
}

extract_value_labels <- function(df, dataset_name) {
  label_rows <- lapply(names(df), function(variable_name) {
    value_labels <- attr(df[[variable_name]], "labels", exact = TRUE)

    if (is.null(value_labels) || length(value_labels) == 0) {
      return(NULL)
    }

    data.frame(
      dataset = dataset_name,
      variable = variable_name,
      raw_value = unname(value_labels),
      value_label = names(value_labels),
      stringsAsFactors = FALSE
    )
  })

  label_rows <- Filter(Negate(is.null), label_rows)

  if (length(label_rows) == 0) {
    return(data.frame(
      dataset = character(0),
      variable = character(0),
      raw_value = numeric(0),
      value_label = character(0),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, label_rows)
}

write_source_metadata <- function(df, dataset_name, out_dir = "out") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  utils::write.csv(
    extract_variable_metadata(df, dataset_name),
    file.path(out_dir, paste0(dataset_name, "_variable_metadata.csv")),
    row.names = FALSE
  )

  utils::write.csv(
    extract_value_labels(df, dataset_name),
    file.path(out_dir, paste0(dataset_name, "_value_labels.csv")),
    row.names = FALSE
  )
}

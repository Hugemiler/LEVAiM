#' Construct a longitudinal multi-assay dataset
#'
#' Stores one or more feature tables together with sample metadata.
#'
#' Each assay must be a matrix-like object with rows as samples and columns as
#' features. Assays can represent taxa, pathways, UniRefs, KOs, ECs,
#' metabolites, biomarkers, or other longitudinally measured features.
#'
#' @param assays A named list of assay tables. Each assay must have rows as
#'   samples and columns as features. A single data.frame or matrix is accepted
#'   and wrapped as `list(features = assays)`.
#' @param metadata A data.frame with rows as samples.
#' @param time_col Name of the metadata column encoding time.
#' @param subject_col Name of the metadata column encoding subject ID.
#'
#' @return A `LongitudinalDataset` object.
#' @export
LongitudinalDataset <- function(
    assays,
    metadata,
    time_col = "time",
    subject_col = "subject"
) {
  if (is.data.frame(assays) || is.matrix(assays)) {
    assays <- list(features = assays)
  }

  if (!is.list(assays) || length(assays) == 0L) {
    cli::cli_abort("`assays` must be a non-empty named list of feature tables.")
  }

  if (is.null(names(assays)) || any(names(assays) == "")) {
    cli::cli_abort("`assays` must be a named list.")
  }

  assays <- lapply(assays, function(x) {
    as.data.frame(x, check.names = FALSE)
  })

  metadata <- as.data.frame(metadata, check.names = FALSE)

  validate_longitudinal_dataset(
    assays = assays,
    metadata = metadata,
    time_col = time_col,
    subject_col = subject_col
  )

  structure(
    list(
      assays = assays,
      metadata = metadata,
      time_col = time_col,
      subject_col = subject_col
    ),
    class = "LongitudinalDataset"
  )
}


#' Validate a longitudinal multi-assay dataset
#'
#' @param assays Named list of assay tables.
#' @param metadata Sample metadata table.
#' @param time_col Metadata column encoding time.
#' @param subject_col Metadata column encoding subject ID.
#'
#' @return Invisibly returns `TRUE`.
#' @keywords internal
validate_longitudinal_dataset <- function(
    assays,
    metadata,
    time_col,
    subject_col
) {
  if (is.null(rownames(metadata)) || any(rownames(metadata) == "")) {
    cli::cli_abort("`metadata` must have sample IDs as row names.")
  }

  missing_cols <- setdiff(c(time_col, subject_col), colnames(metadata))

  if (length(missing_cols) > 0L) {
    cli::cli_abort(
      "Missing required metadata column(s): {paste(missing_cols, collapse = ', ')}"
    )
  }

  for (assay_name in names(assays)) {
    assay <- assays[[assay_name]]

    if (is.null(rownames(assay)) || any(rownames(assay) == "")) {
      cli::cli_abort(
        "Assay {.val {assay_name}} must have sample IDs as row names."
      )
    }

    if (!identical(rownames(assay), rownames(metadata))) {
      cli::cli_abort(
        "Assay {.val {assay_name}} and `metadata` must have identical sample order."
      )
    }

    non_numeric <- !vapply(assay, is.numeric, logical(1))

    if (any(non_numeric)) {
      cli::cli_abort(
        "Assay {.val {assay_name}} has non-numeric feature column(s): {paste(names(assay)[non_numeric], collapse = ', ')}"
      )
    }
  }

  invisible(TRUE)
}


#' Print a LongitudinalDataset
#'
#' @param x A `LongitudinalDataset`.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.LongitudinalDataset <- function(x, ...) {
  cli::cli_h1("LongitudinalDataset")

  cli::cli_text("Samples: {nrow(x$metadata)}")
  cli::cli_text("Assays: {length(x$assays)}")
  cli::cli_text("Time column: {.field {x$time_col}}")
  cli::cli_text("Subject column: {.field {x$subject_col}}")

  cli::cli_ul(names(x$assays))

  invisible(x)
}

#' Define an assay feature response
#'
#' @param assay Name of the assay in a `LongitudinalDataset`.
#' @param feature Name of the feature within the assay.
#'
#' @return A LEVAiM response object.
#' @export
assay_feature <- function(assay, feature) {
  structure(
    list(
      source = "assay",
      assay = assay,
      feature = feature
    ),
    class = c("levaim_assay_feature", "levaim_response")
  )
}


#' Define a metadata feature response
#'
#' @param column Name of the metadata column.
#'
#' @return A LEVAiM response object.
#' @export
metadata_feature <- function(column) {
  structure(
    list(
      source = "metadata",
      column = column
    ),
    class = c("levaim_metadata_feature", "levaim_response")
  )
}


#' Define the subject identifier
#'
#' @param column Metadata column identifying subjects.
#'
#' @return A LEVAiM subject variable.
#' @export
subject_var <- function(column) {
  structure(
    list(column = column),
    class = c("levaim_subject_var", "levaim_var")
  )
}


#' Define the time variable
#'
#' @param column Metadata column encoding chronology.
#'
#' @return A LEVAiM time variable.
#' @export
time_var <- function(column) {
  structure(
    list(column = column),
    class = c("levaim_time_var", "levaim_var")
  )
}


#' Define a categorical model effect
#'
#' @param column Metadata column encoding a categorical/grouping effect.
#' @param reference Optional reference level.
#' @param role Covariate role. `"auto"` infers the role from the model frame;
#'   `"subject_static"` requires one value per subject; `"time_varying"` allows
#'   values to change within subject; `"sample_level"` marks sample-specific
#'   annotations; `"derived_exposure"` marks columns derived from a longitudinal
#'   exposure history.
#'
#' @return A LEVAiM covariate object.
#' @export
group_effect <- function(
    column,
    reference = NULL,
    role = c("auto", "subject_static", "time_varying", "sample_level", "derived_exposure")
) {
  role <- match.arg(role)

  structure(
    list(
      column = column,
      reference = reference,
      type = "group",
      role = role
    ),
    class = c("levaim_group_effect", "levaim_covariate")
  )
}


#' Define a continuous model effect
#'
#' @param column Metadata column encoding a continuous covariate.
#' @param scale Logical. Whether to scale the covariate before modeling.
#'
#' @return A LEVAiM covariate object.
#' @export
continuous_effect <- function(column, scale = FALSE) {
  structure(
    list(
      column = column,
      scale = scale,
      type = "continuous"
    ),
    class = c("levaim_continuous_effect", "levaim_covariate")
  )
}


#' Define a random effect
#'
#' @param column Metadata column encoding a random-effect grouping variable.
#'
#' @return A LEVAiM covariate object.
#' @export
random_effect <- function(column) {
  structure(
    list(
      column = column,
      type = "random"
    ),
    class = c("levaim_random_effect", "levaim_covariate")
  )
}


#' Define a longitudinal model specification
#'
#' Declares the response, subject identifier, chronology variable, and model
#' covariates. This object does not define whether trajectories are shared or
#' vary by covariates; that is handled later by `trajectory()`.
#'
#' @param response Response definition, usually from `assay_feature()` or
#'   `metadata_feature()`.
#' @param subject Subject identifier, from `subject_var()`.
#' @param time Chronology variable, from `time_var()`.
#' @param covariates List of covariate definitions.
#'
#' @return A `levaim_model_spec` object.
#' @export
model_spec <- function(
    response,
    subject,
    time,
    covariates = list()
) {
  validate_model_spec_inputs(
    response = response,
    subject = subject,
    time = time,
    covariates = covariates
  )

  structure(
    list(
      response = response,
      subject = subject,
      time = time,
      covariates = covariates
    ),
    class = "levaim_model_spec"
  )
}


#' Validate model specification inputs
#'
#' @keywords internal
validate_model_spec_inputs <- function(
    response,
    subject,
    time,
    covariates
) {
  if (!inherits(response, "levaim_response")) {
    cli::cli_abort("`response` must be created with `assay_feature()` or `metadata_feature()`.")
  }

  if (!inherits(subject, "levaim_subject_var")) {
    cli::cli_abort("`subject` must be created with `subject_var()`.")
  }

  if (!inherits(time, "levaim_time_var")) {
    cli::cli_abort("`time` must be created with `time_var()`.")
  }

  if (!is.list(covariates)) {
    cli::cli_abort("`covariates` must be a list of covariate objects.")
  }

  bad <- !vapply(covariates, inherits, logical(1), what = "levaim_covariate")

  if (any(bad)) {
    cli::cli_abort("All entries in `covariates` must be created with covariate constructors.")
  }

  invisible(TRUE)
}


#' Print a LEVAiM model specification
#'
#' @param x A `levaim_model_spec`.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.levaim_model_spec <- function(x, ...) {
  cli::cli_h1("LEVAiM model specification")

  cli::cli_h2("Response")

  if (x$response$source == "assay") {
    cli::cli_text("Assay: {.field {x$response$assay}}")
    cli::cli_text("Feature: {.field {x$response$feature}}")
  } else {
    cli::cli_text("Metadata column: {.field {x$response$column}}")
  }

  cli::cli_h2("Longitudinal structure")
  cli::cli_text("Subject: {.field {x$subject$column}}")
  cli::cli_text("Time: {.field {x$time$column}}")

  cli::cli_h2("Covariates")

  if (length(x$covariates) == 0L) {
    cli::cli_text("No covariates declared.")
  } else {
    for (cov in x$covariates) {
      role <- cov$role %||% "auto"
      cli::cli_text("- {.field {cov$column}} [{cov$type}; role={role}]")
    }
  }

  invisible(x)
}

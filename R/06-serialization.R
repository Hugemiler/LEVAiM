#####
# LEVAiM serialization
#####

#' Save a LEVAiM object
#'
#' Serializes a LEVAiM fit, result summary, model frame, or cross-validation
#' object with R's native RDS format.
#'
#' @param x A LEVAiM object.
#' @param file Path to the output `.rds` file.
#' @param ... Additional arguments passed to [saveRDS()].
#'
#' @return Invisibly returns `file`.
#' @export
save_levaim_object <- function(x, file, ...) {
  if (!is_levaim_serializable(x)) {
    cli::cli_abort(
      "`x` must be a LEVAiM fit, results, model frame, cross-validation, dataset, specification, trajectory, or control object."
    )
  }

  saveRDS(x, file = file, ...)

  invisible(file)
}


#' Load a LEVAiM object
#'
#' Reads a LEVAiM object saved with [save_levaim_object()] and validates that
#' the recovered object still has a recognized LEVAiM class.
#'
#' @param file Path to an `.rds` file produced by [save_levaim_object()].
#' @param ... Additional arguments passed to [readRDS()].
#'
#' @return The deserialized LEVAiM object.
#' @export
load_levaim_object <- function(file, ...) {
  x <- readRDS(file = file, ...)

  if (!is_levaim_serializable(x)) {
    cli::cli_abort(
      "The file does not contain a recognized LEVAiM object."
    )
  }

  x
}


#' Check whether an object uses a serializable LEVAiM class
#'
#' @param x An object.
#'
#' @return `TRUE` when `x` has a recognized LEVAiM class.
#' @keywords internal
is_levaim_serializable <- function(x) {
  inherits(
    x,
    c(
      "LongitudinalDataset",
      "levaim_model_spec",
      "levaim_trajectory",
      "levaim_trajectory_control",
      "levaim_model_frame",
      "levaim_fit",
      "levaim_results",
      "levaim_cv"
    )
  )
}

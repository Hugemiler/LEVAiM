#' Define a naive (non-varying) trajectory
#'
#' A naive trajectory uses a single shared time function across all samples.
#'
#' @return A LEVAiM trajectory design object.
#' @export
naive_trajectory <- function() {
  structure(
    list(type = "shared", by = NULL, interaction = FALSE),
    class = "levaim_trajectory_design"
  )
}


#' Define a covariate-varying trajectory
#'
#' A varying trajectory allows the time function to differ across levels of
#' one or more declared covariates.
#'
#' @param by Character vector of covariate column names.
#' @param interaction Logical. If `TRUE`, multiple `by` variables are combined
#'   into a single interaction factor.
#'
#' @return A LEVAiM trajectory design object.
#' @export
varying_trajectory <- function(by, interaction = length(by) > 1L) {
  if (!is.character(by) || length(by) == 0L) {
    cli::cli_abort("`by` must be a non-empty character vector.")
  }

  structure(
    list(
      type = "varying",
      by = by,
      interaction = interaction
    ),
    class = "levaim_trajectory_design"
  )
}


#' Create a trajectory from a model specification
#'
#' Combines a model specification with a trajectory design and validates that
#' requested trajectory modifiers are declared covariates.
#'
#' @param spec A `levaim_model_spec` object.
#' @param design A trajectory design from `naive_trajectory()` or
#'   `varying_trajectory()`.
#'
#' @return A `levaim_trajectory` object.
#' @export
trajectory <- function(spec, design = naive_trajectory()) {
  if (!inherits(spec, "levaim_model_spec")) {
    cli::cli_abort("`spec` must be created with `model_spec()`.")
  }

  if (!inherits(design, "levaim_trajectory_design")) {
    cli::cli_abort("`design` must be created with `naive_trajectory()` or `varying_trajectory()`.")
  }

  validate_trajectory_design(spec, design)

  structure(
    list(
      spec = spec,
      design = design
    ),
    class = "levaim_trajectory"
  )
}


#' Validate a trajectory design against a model specification
#'
#' @keywords internal
validate_trajectory_design <- function(spec, design) {
  declared_covariates <- vapply(
    spec$covariates,
    function(x) x$column,
    character(1)
  )

  if (design$type == "varying") {
    missing <- setdiff(design$by, declared_covariates)

    if (length(missing) > 0L) {
      cli::cli_abort(
        "Trajectory modifier(s) not declared as covariates: {paste(missing, collapse = ', ')}"
      )
    }

    cov_types <- vapply(
      spec$covariates,
      function(x) x$type,
      character(1)
    )
    names(cov_types) <- declared_covariates

    bad <- design$by[cov_types[design$by] == "random"]

    if (length(bad) > 0L) {
      cli::cli_abort(
        "Random effect(s) cannot define varying trajectories: {paste(bad, collapse = ', ')}"
      )
    }
  }

  invisible(TRUE)
}


#' Print a LEVAiM trajectory
#'
#' @param x A `levaim_trajectory`.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.levaim_trajectory <- function(x, ...) {
  cli::cli_h1("LEVAiM trajectory")

  cli::cli_text("Response:")
  if (x$spec$response$source == "assay") {
    cli::cli_text("  {x$spec$response$assay} / {x$spec$response$feature}")
  } else {
    cli::cli_text("  metadata / {x$spec$response$column}")
  }

  cli::cli_text("Time: {.field {x$spec$time$column}}")
  cli::cli_text("Subject: {.field {x$spec$subject$column}}")

  cli::cli_text("Trajectory design: {.field {x$design$type}}")

  if (x$design$type == "varying") {
    cli::cli_text("Varying by: {.field {paste(x$design$by, collapse = ' + ')}}")
    cli::cli_text("Interaction: {x$design$interaction}")
  }

  invisible(x)
}

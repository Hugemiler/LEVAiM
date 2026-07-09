#####
# Model frame construction
#####

#' Build a model frame for a LEVAiM trajectory
#'
#' Creates the response column `.y` and joins metadata. This function does not
#' precompute formula interactions; those are left to the modeling backend.
#'
#' @param ds A `LongitudinalDataset`.
#' @param traj A `levaim_trajectory`.
#' @param transform Response transformation. One of `"identity"`, `"log1p"`,
#'   or `"log10p"`.
#'
#' @return A data.frame ready for model fitting.
#' @keywords internal
make_model_frame <- function(
    ds,
    traj,
    transform = c("identity", "log1p", "log10p")
) {
  transform <- match.arg(transform)

  if (!inherits(ds, "LongitudinalDataset")) {
    cli::cli_abort("`ds` must be a `LongitudinalDataset` object.")
  }

  if (!inherits(traj, "levaim_trajectory")) {
    cli::cli_abort("`traj` must be a `levaim_trajectory` object.")
  }

  spec <- traj$spec

  y <- extract_response(ds, spec$response)
  y <- transform_response(y, transform)

  meta <- ds$metadata

  if (!identical(names(y), rownames(meta))) {
    cli::cli_abort("Response vector and metadata must have identical sample order.")
  }

  mf <- as.data.frame(meta, check.names = FALSE)
  mf$.y <- as.numeric(y)

  validate_model_frame_columns(mf, traj)

  mf
}

#' Build a LEVAiM model frame
#'
#' Materializes the response, metadata, trajectory, formula, and modeling
#' controls needed for fitting. This object is backend-aware but not yet fitted.
#'
#' @param ds A `LongitudinalDataset`.
#' @param traj A `levaim_trajectory`.
#' @param control A `levaim_trajectory_control`.
#'
#' @return A `levaim_model_frame` object.
#' @export
model_frame <- function(
    ds,
    traj,
    control = trajectory_control()
) {
  if (!inherits(control, "levaim_trajectory_control")) {
    cli::cli_abort("`control` must be created with `trajectory_control()`.")
  }

  data <- make_model_frame(
    ds = ds,
    traj = traj,
    transform = control$transform
  )

  spline_k <- NULL
  if (control$engine == "spline") {
    time_col <- traj$spec$time$column
    spline_k <- min(control$spline$k, length(unique(data[[time_col]])))
  }

  formula <- compile_formula(
    traj = traj,
    engine = control$engine,
    spline_k = spline_k
  )

  structure(
    list(
      data = data,
      trajectory = traj,
      control = control,
      formula = formula
    ),
    class = "levaim_model_frame"
  )
}

#' Extract response vector from a LongitudinalDataset
#'
#' @keywords internal
extract_response <- function(ds, response) {
  if (inherits(response, "levaim_assay_feature")) {
    assay_name <- response$assay
    feature_name <- response$feature

    if (!assay_name %in% names(ds$assays)) {
      cli::cli_abort("Assay {.val {assay_name}} not found in `ds$assays`.")
    }

    assay <- ds$assays[[assay_name]]

    if (!feature_name %in% colnames(assay)) {
      cli::cli_abort(
        "Feature {.val {feature_name}} not found in assay {.val {assay_name}}."
      )
    }

    y <- assay[[feature_name]]
    names(y) <- rownames(assay)

    return(y)
  }

  if (inherits(response, "levaim_metadata_feature")) {
    column <- response$column

    if (!column %in% colnames(ds$metadata)) {
      cli::cli_abort("Metadata response column {.val {column}} not found.")
    }

    y <- ds$metadata[[column]]
    names(y) <- rownames(ds$metadata)

    return(y)
  }

  cli::cli_abort("Unknown response type.")
}


#' Transform a response vector
#'
#' @keywords internal
transform_response <- function(y, transform) {
  if (!is.numeric(y)) {
    cli::cli_abort("The response must be numeric for trajectory modeling.")
  }

  switch(
    transform,
    identity = y,
    log1p = log1p(y),
    log10p = log10(y + 1)
  )
}


#' Validate model frame columns required by trajectory
#'
#' @keywords internal
validate_model_frame_columns <- function(mf, traj) {
  spec <- traj$spec
  design <- traj$design

  required <- c(
    spec$subject$column,
    spec$time$column,
    vapply(spec$covariates, function(x) x$column, character(1))
  )

  if (design$type == "varying") {
    required <- c(required, design$by)
  }

  missing <- setdiff(unique(required), colnames(mf))

  if (length(missing) > 0L) {
    cli::cli_abort(
      "Missing required model-frame column(s): {paste(missing, collapse = ', ')}"
    )
  }

  invisible(TRUE)
}

#' Print a LEVAiM model frame
#'
#' @param x A `levaim_model_frame`.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.levaim_model_frame <- function(x, ...) {
  cli::cli_h1("LEVAiM model frame")

  cli::cli_text("Engine: {.field {x$control$engine}}")
  cli::cli_text("Family: {.field {x$control$family}}")
  cli::cli_text("Transform: {.field {x$control$transform}}")
  cli::cli_text("Samples: {nrow(x$data)}")
  cli::cli_text("Formula: {.code {deparse(x$formula)}}")

  invisible(x)
}

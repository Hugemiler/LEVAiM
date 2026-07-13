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

  preprocessed <- preprocess_model_frame_annotations(
    data = data,
    traj = traj,
    control = control
  )

  data <- preprocessed$data
  traj <- preprocessed$trajectory

  validate_response_family(data$.y, control$family)

  spline_k <- NULL
  if (control$engine == "spline") {
    time_col <- traj$spec$time$column
    spline_k <- min(control$spline$k, length(unique(data[[time_col]])))
  }

  formula <- compile_formula(
    traj = traj,
    engine = control$engine,
    spline_k = spline_k,
    spline_basis = control$spline$basis,
    gp_kernel = control$gp$kernel,
    gp_basis_k = control$gp$basis_k
  )

  structure(
    list(
      data = data,
      trajectory = traj,
      control = control,
      formula = formula,
      preprocessing = preprocessed$preprocessing
    ),
    class = "levaim_model_frame"
  )
}


#' Validate the response domain for a model family
#'
#' @keywords internal
validate_response_family <- function(y, family) {
  observed <- y[!is.na(y)]

  if (length(observed) == 0L || any(!is.finite(observed))) {
    cli::cli_abort("The modeled response must contain finite observed values.")
  }

  if (family == "binomial" && any(!observed %in% c(0, 1))) {
    cli::cli_abort("`family = 'binomial'` requires a response containing only 0 and 1.")
  }

  if (family == "beta" && any(observed <= 0 | observed >= 1)) {
    cli::cli_abort(c(
      "`family = 'beta'` requires every response value to lie strictly between 0 and 1.",
      "i" = "Zeros and ones require an explicit zero/one-inflated model; LEVAiM does not transform them implicitly."
    ))
  }

  invisible(TRUE)
}


#' Preprocess sparse and missing model-frame annotations
#'
#' @keywords internal
preprocess_model_frame_annotations <- function(data, traj, control) {
  spec <- traj$spec
  design <- traj$design

  preprocessing <- list(
    na_action = control$na_action,
    dropped_rows = integer(),
    missing_levels = list(),
    collapsed_levels = list(),
    dropped_covariates = character(),
    covariate_roles = data.frame()
  )

  always_complete <- unique(c(".y", spec$subject$column, spec$time$column))
  continuous_cols <- vapply(
    spec$covariates,
    function(cov) {
      if (inherits(cov, "levaim_continuous_effect")) {
        return(cov$column)
      }
      NA_character_
    },
    character(1)
  )
  continuous_cols <- continuous_cols[!is.na(continuous_cols)]

  complete_cols <- unique(c(always_complete, continuous_cols))
  if (control$na_action == "complete") {
    complete_cols <- unique(c(
      complete_cols,
      vapply(spec$covariates, function(cov) cov$column, character(1)),
      design$by
    ))
  }

  keep <- stats::complete.cases(data[, complete_cols, drop = FALSE])
  if (!all(keep)) {
    preprocessing$dropped_rows <- which(!keep)
    data <- data[keep, , drop = FALSE]
  }

  role_diagnostics <- covariate_role_diagnostics(data, traj)
  validate_covariate_roles(role_diagnostics, traj)
  preprocessing$covariate_roles <- role_diagnostics

  data[[spec$subject$column]] <- as_sparse_annotation_factor(
    data[[spec$subject$column]],
    column = spec$subject$column,
    control = control,
    preprocessing = preprocessing,
    collapse_sparse = FALSE
  )$values

  if (nlevels(data[[spec$subject$column]]) < 2L) {
    cli::cli_abort(
      "Subject column {.field {spec$subject$column}} must contain at least two levels after preprocessing."
    )
  }

  categorical_cols <- unique(c(
    vapply(
      spec$covariates,
      function(cov) {
        if (inherits(cov, c("levaim_group_effect", "levaim_random_effect"))) {
          return(cov$column)
        }
        NA_character_
      },
      character(1)
    ),
    design$by
  ))
  categorical_cols <- categorical_cols[!is.na(categorical_cols)]

  for (column in categorical_cols) {
    prepared <- as_sparse_annotation_factor(
      data[[column]],
      column = column,
      control = control,
      preprocessing = preprocessing
    )

    data[[column]] <- prepared$values
    preprocessing <- prepared$preprocessing
  }

  scaled_cols <- character()
  for (cov in spec$covariates) {
    if (inherits(cov, "levaim_continuous_effect") && isTRUE(cov$scale)) {
      column <- cov$column
      data[[column]] <- as.numeric(scale(data[[column]]))
      scaled_cols <- c(scaled_cols, column)
    }
  }
  preprocessing$scaled_continuous <- scaled_cols

  invariant_cols <- unique(categorical_cols)[
    vapply(
      unique(categorical_cols),
      function(column) nlevels(data[[column]]) < 2L,
      logical(1)
    )
  ]

  if (length(invariant_cols) > 0L) {
    varying_invariant <- intersect(invariant_cols, design$by)

    if (length(varying_invariant) > 0L) {
      cli::cli_abort(
        "Varying trajectory column(s) have fewer than two levels after preprocessing: {paste(varying_invariant, collapse = ', ')}"
      )
    }

    if (isTRUE(control$annotations$drop_invariant_covariates)) {
      spec$covariates <- Filter(
        function(cov) {
          !cov$column %in% invariant_cols
        },
        spec$covariates
      )
      preprocessing$dropped_covariates <- invariant_cols
    } else {
      cli::cli_abort(
        "Covariate column(s) have fewer than two levels after preprocessing: {paste(invariant_cols, collapse = ', ')}"
      )
    }
  }

  traj$spec <- spec

  list(
    data = data,
    trajectory = traj,
    preprocessing = preprocessing
  )
}


#' Infer and validate covariate roles
#'
#' @keywords internal
covariate_role_diagnostics <- function(data, traj) {
  spec <- traj$spec
  subject_col <- spec$subject$column

  covariate_columns <- vapply(spec$covariates, function(cov) cov$column, character(1))
  declared_roles <- vapply(
    spec$covariates,
    function(cov) cov$role %||% "auto",
    character(1)
  )
  names(declared_roles) <- covariate_columns

  design_columns <- traj$design$by
  columns <- unique(c(covariate_columns, design_columns))

  if (length(columns) == 0L) {
    return(data.frame(
      column = character(),
      declared_role = character(),
      inferred_role = character(),
      n_subjects = integer(),
      subjects_with_changes = integer(),
      n_levels = integer(),
      used_as_trajectory_modifier = logical(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(columns, function(column) {
    subject_change <- values_change_within_subject(
      values = data[[column]],
      subjects = data[[subject_col]]
    )
    n_unique <- length(unique(as.character(data[[column]])[!is.na(data[[column]])]))

    inferred_role <- if (subject_change$subjects_with_changes == 0L) {
      "subject_static"
    } else if (n_unique > max(20L, floor(nrow(data) / 2))) {
      "sample_level"
    } else {
      "time_varying"
    }

    declared_role <- if (column %in% names(declared_roles)) {
      unname(declared_roles[[column]])
    } else {
      "auto"
    }

    data.frame(
      column = column,
      declared_role = declared_role,
      inferred_role = inferred_role,
      n_subjects = length(unique(data[[subject_col]])),
      subjects_with_changes = subject_change$subjects_with_changes,
      n_levels = n_unique,
      used_as_trajectory_modifier = column %in% design_columns,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' @keywords internal
values_change_within_subject <- function(values, subjects) {
  pieces <- split(as.character(values), subjects, drop = TRUE)
  changes <- vapply(pieces, function(x) {
    x <- x[!is.na(x)]
    length(unique(x)) > 1L
  }, logical(1))

  list(subjects_with_changes = sum(changes))
}


#' @keywords internal
validate_covariate_roles <- function(role_diagnostics, traj) {
  if (nrow(role_diagnostics) == 0L) {
    return(invisible(TRUE))
  }

  subject_static_mismatch <- role_diagnostics[
    role_diagnostics$declared_role == "subject_static" &
      role_diagnostics$subjects_with_changes > 0L,
    ,
    drop = FALSE
  ]

  if (nrow(subject_static_mismatch) > 0L) {
    cli::cli_abort(
      "Covariate(s) declared as `role = 'subject_static'` change within subject: {paste(subject_static_mismatch$column, collapse = ', ')}"
    )
  }

  direct_time_varying <- role_diagnostics[
    role_diagnostics$declared_role == "time_varying" &
      role_diagnostics$used_as_trajectory_modifier,
    ,
    drop = FALSE
  ]

  if (nrow(direct_time_varying) > 0L) {
    cli::cli_warn(
      "Time-varying covariate(s) used directly as trajectory modifiers: {paste(direct_time_varying$column, collapse = ', ')}. This models current observed state, not exposure history. Use `derive_time_varying_features()` for history-based questions."
    )
  }

  sample_modifier <- role_diagnostics[
    role_diagnostics$declared_role == "sample_level" &
      role_diagnostics$used_as_trajectory_modifier,
    ,
    drop = FALSE
  ]

  if (nrow(sample_modifier) > 0L) {
    cli::cli_warn(
      "Sample-level covariate(s) used as trajectory modifiers: {paste(sample_modifier$column, collapse = ', ')}. Confirm that this is intended."
    )
  }

  invisible(TRUE)
}


#' Convert a categorical annotation to a sparse-aware factor
#'
#' @keywords internal
as_sparse_annotation_factor <- function(
    x,
    column,
    control,
    preprocessing,
    collapse_sparse = TRUE
) {
  values <- as.character(x)
  values[!nzchar(values)] <- NA_character_

  missing_idx <- is.na(values)
  if (any(missing_idx)) {
    values[missing_idx] <- control$annotations$missing_level
    preprocessing$missing_levels[[column]] <- control$annotations$missing_level
  }

  if (isTRUE(collapse_sparse)) {
    counts <- table(values)
    sparse_levels <- names(counts)[counts < control$annotations$sparse_min_n]

    if (length(sparse_levels) > 0L) {
      values[values %in% sparse_levels] <- control$annotations$sparse_other_level
      preprocessing$collapsed_levels[[column]] <- sparse_levels
    }
  }

  list(
    values = factor(values),
    preprocessing = preprocessing
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

  if (!is.null(x$preprocessing)) {
    cli::cli_text("Dropped rows: {length(x$preprocessing$dropped_rows)}")

    if (length(x$preprocessing$dropped_covariates) > 0L) {
      cli::cli_text("Dropped covariates: {.field {paste(x$preprocessing$dropped_covariates, collapse = ', ')}}")
    }
  }

  invisible(x)
}

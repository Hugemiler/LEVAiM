#####
# Trajectory comparison summaries
#####

#' Compare fitted trajectories across features
#'
#' Builds a table-first comparison layer for fitted LEVAiM trajectories. The
#' default view compares feature trajectories within the same metadata level,
#' which supports heatmaps such as feature-by-feature distances for each group.
#'
#' @param fits A list of `levaim_fit` objects or a `levaim_feature_tests` object
#'   created with `keep_fits = TRUE`.
#' @param n Number of time grid points when `time_values` is not supplied.
#' @param time_values Optional numeric vector of time values.
#' @param windows Optional data.frame defining time windows. Accepted columns are
#'   `window_label`, `window_start`, and `window_end`, or `label`, `start`, and
#'   `end`.
#' @param view Comparison view. `"same_group_level"` compares features within
#'   the same trajectory group, `"all_trajectory_units"` compares every fitted
#'   trajectory unit, and `"overall"` currently supports shared trajectories.
#' @param metrics Metrics to compute.
#' @param standardize Logical. Whether to z-score each trajectory within each
#'   comparison window before computing correlations and distances.
#' @param inference Inference mode. `"descriptive"` reports correlation
#'   p-values from [stats::cor.test()] with an explicit fitted-grid warning.
#'   `"none"` suppresses p-values. `"subject_permutation"` is reserved for a
#'   future subject-aware resampling test.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#'
#' @return A `levaim_trajectory_comparison` data.frame. Each row is one
#'   feature-pair, trajectory-group, time-window, and metric comparison.
#' @export
compare_trajectories <- function(
    fits,
    n = 100L,
    time_values = NULL,
    windows = NULL,
    view = c("same_group_level", "overall", "all_trajectory_units"),
    metrics = c("fitted_correlation", "fitted_distance", "derivative_correlation"),
    standardize = TRUE,
    inference = c("descriptive", "none", "subject_permutation"),
    p_adjust_method = "BH"
) {
  view <- match.arg(view)
  inference <- match.arg(inference)

  if (inference == "subject_permutation") {
    cli::cli_abort(
      "`inference = 'subject_permutation'` is not implemented yet. Use {.val descriptive} or {.val none}."
    )
  }

  metrics <- match.arg(
    metrics,
    c("fitted_correlation", "fitted_distance", "derivative_correlation"),
    several.ok = TRUE
  )

  fits <- normalize_comparison_fits(fits)
  if (length(fits) < 2L) {
    cli::cli_abort("Trajectory comparison requires at least two fitted trajectories.")
  }

  units <- extract_comparison_units(
    fits = fits,
    n = n,
    time_values = time_values,
    include_derivatives = "derivative_correlation" %in% metrics
  )

  if (view == "overall" && any(units$key$trajectory_type != "shared")) {
    cli::cli_abort(
      "`view = 'overall'` currently supports shared trajectories only. Use {.val same_group_level} or {.val all_trajectory_units} for varying trajectories."
    )
  }

  windows <- normalize_comparison_windows(
    windows = windows,
    observed_time = units$aligned_effects$.time
  )

  table <- build_trajectory_comparison_table(
    key = units$key,
    effects = units$aligned_effects,
    derivatives = units$aligned_derivatives,
    windows = windows,
    view = view,
    metrics = metrics,
    standardize = standardize,
    inference = inference
  )

  table$q_value <- NA_real_
  valid_p <- is.finite(table$p_value)
  if (any(valid_p)) {
    table$q_value[valid_p] <- stats::p.adjust(
      table$p_value[valid_p],
      method = p_adjust_method
    )
  }

  table <- table[, comparison_table_columns(), drop = FALSE]

  attr(table, "settings") <- list(
    view = view,
    metrics = metrics,
    standardize = standardize,
    inference = inference,
    p_adjust_method = p_adjust_method,
    n = n,
    time_values = time_values,
    windows = windows
  )

  class(table) <- c("levaim_trajectory_comparison", class(table))
  table
}


#' Extract a trajectory distance matrix
#'
#' @param comparison A `levaim_trajectory_comparison` data.frame.
#' @param metric Metric used to select matrices. Currently defaults to
#'   `"fitted_distance"`.
#' @param matrix Optional matrix name. If omitted, the first matching matrix is
#'   returned.
#'
#' @return A numeric square matrix.
#' @export
trajectory_distance_matrix <- function(
    comparison,
    metric = "fitted_distance",
    matrix = NULL
) {
  if (!inherits(comparison, "levaim_trajectory_comparison")) {
    cli::cli_abort("`comparison` must be created with `compare_trajectories()`.")
  }

  matrices <- build_distance_matrices(comparison)
  if (length(matrices) == 0L) {
    cli::cli_abort("No trajectory distance matrices are available.")
  }

  keep <- vapply(
    matrices,
    function(x) identical(attr(x, "metric"), metric),
    logical(1)
  )
  matrices <- matrices[keep]

  if (length(matrices) == 0L) {
    cli::cli_abort("No distance matrix is available for metric {.val {metric}}.")
  }

  if (is.null(matrix)) {
    return(matrices[[1L]])
  }

  if (!matrix %in% names(matrices)) {
    cli::cli_abort(
      "Matrix {.val {matrix}} was not found. Available matrices: {paste(names(matrices), collapse = ', ')}"
    )
  }

  matrices[[matrix]]
}


#' @export
print.levaim_trajectory_comparison <- function(x, ...) {
  cli::cli_h1("LEVAiM trajectory comparison")
  cli::cli_text("Rows: {nrow(x)}")
  print(utils::head(as.data.frame(x), 10L), row.names = FALSE)
  invisible(x)
}


#' @keywords internal
normalize_comparison_fits <- function(fits) {
  if (inherits(fits, "levaim_fit")) {
    fits <- list(fit_1 = fits)
  }

  if (inherits(fits, "levaim_feature_tests")) {
    if (is.null(fits$fits)) {
      cli::cli_abort(
        "`fits` was created by `test_assay_features()` without retained fits. Re-run with `keep_fits = TRUE`."
      )
    }
    fits <- fits$fits
  }

  if (!is.list(fits) || length(fits) == 0L) {
    cli::cli_abort("`fits` must be a non-empty list of `levaim_fit` objects.")
  }

  keep <- vapply(fits, inherits, logical(1), what = "levaim_fit")
  fits <- fits[keep]

  if (length(fits) == 0L) {
    cli::cli_abort("`fits` must contain at least one `levaim_fit` object.")
  }

  if (is.null(names(fits))) {
    names(fits) <- rep("", length(fits))
  }

  empty_names <- names(fits) == ""
  if (any(empty_names)) {
    names(fits)[empty_names] <- paste0("fit_", which(empty_names))
  }

  fits
}


#' @keywords internal
extract_comparison_units <- function(
    fits,
    n,
    time_values,
    include_derivatives
) {
  keys <- list()
  effect_rows <- list()
  derivative_rows <- list()
  unit_i <- 0L

  for (fit_name in names(fits)) {
    fit <- fits[[fit_name]]
    mf <- fit$model_frame
    spec <- mf$trajectory$spec
    design <- mf$trajectory$design
    time_col <- spec$time$column
    response <- response_components(spec$response)
    effects <- trajectory_effects(fit, n = n, time_values = time_values)

    derivatives <- NULL
    if (isTRUE(include_derivatives) && fit$engine %in% c("spline", "gp")) {
      derivatives <- trajectory_derivatives(fit, n = n, time_values = time_values)
    }

    by_cols <- if (design$type == "varying") design$by else character()
    split_key <- comparison_split_key(effects, by_cols)
    effect_groups <- split(effects, split_key, drop = TRUE)
    derivative_groups <- if (!is.null(derivatives)) {
      split(derivatives, comparison_split_key(derivatives, by_cols), drop = TRUE)
    } else {
      list()
    }

    for (group_name in names(effect_groups)) {
      unit_i <- unit_i + 1L
      unit_id <- paste0("trajectory_", unit_i)
      effect_piece <- effect_groups[[group_name]]
      group_info <- comparison_group_info(effect_piece, by_cols)

      keys[[unit_i]] <- data.frame(
        trajectory_id = unit_id,
        source_name = fit_name,
        assay = response$assay,
        feature = response$feature,
        response_type = response$type,
        engine = fit$engine,
        trajectory_type = design$type,
        group_variable = group_info$variable,
        group_level = group_info$level,
        trajectory_group = group_info$label,
        stringsAsFactors = FALSE
      )

      effect_rows[[unit_i]] <- data.frame(
        trajectory_id = unit_id,
        .time = as.numeric(effect_piece[[time_col]]),
        .fitted = as.numeric(effect_piece$.fitted),
        stringsAsFactors = FALSE
      )

      derivative_piece <- derivative_groups[[group_name]]
      if (!is.null(derivative_piece)) {
        derivative_rows[[unit_i]] <- data.frame(
          trajectory_id = unit_id,
          .time = as.numeric(derivative_piece[[time_col]]),
          .derivative = as.numeric(derivative_piece$.derivative),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  key <- do.call(rbind, keys)
  rownames(key) <- NULL

  aligned_effects <- do.call(rbind, effect_rows)
  rownames(aligned_effects) <- NULL

  derivative_rows <- derivative_rows[!vapply(derivative_rows, is.null, logical(1))]
  aligned_derivatives <- if (length(derivative_rows) > 0L) {
    out <- do.call(rbind, derivative_rows)
    rownames(out) <- NULL
    out
  } else {
    data.frame(
      trajectory_id = character(),
      .time = numeric(),
      .derivative = numeric(),
      stringsAsFactors = FALSE
    )
  }

  list(
    key = key,
    aligned_effects = aligned_effects,
    aligned_derivatives = aligned_derivatives
  )
}


#' @keywords internal
response_components <- function(response) {
  if (inherits(response, "levaim_assay_feature")) {
    return(list(
      type = "assay",
      assay = response$assay,
      feature = response$feature
    ))
  }

  if (inherits(response, "levaim_metadata_feature")) {
    return(list(
      type = "metadata",
      assay = "metadata",
      feature = response$column
    ))
  }

  list(type = "unknown", assay = NA_character_, feature = "unknown")
}


#' @keywords internal
comparison_split_key <- function(data, by_cols) {
  if (length(by_cols) == 0L) {
    return(factor(rep("shared", nrow(data))))
  }

  interaction(data[by_cols], drop = TRUE)
}


#' @keywords internal
comparison_group_info <- function(data, by_cols) {
  if (length(by_cols) == 0L) {
    return(list(
      variable = NA_character_,
      level = "shared",
      label = "shared"
    ))
  }

  values <- vapply(by_cols, function(column) {
    as.character(data[[column]][[1L]])
  }, character(1))

  pieces <- paste0(by_cols, "=", values)

  list(
    variable = paste(by_cols, collapse = "+"),
    level = paste(values, collapse = "+"),
    label = paste(pieces, collapse = ";")
  )
}


#' @keywords internal
normalize_comparison_windows <- function(windows, observed_time) {
  if (is.null(windows)) {
    return(data.frame(
      window_label = "full",
      window_start = min(observed_time, na.rm = TRUE),
      window_end = max(observed_time, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
  }

  if (!is.data.frame(windows) || nrow(windows) == 0L) {
    cli::cli_abort("`windows` must be a non-empty data.frame.")
  }

  if (!"window_label" %in% names(windows) && "label" %in% names(windows)) {
    windows$window_label <- windows$label
  }
  if (!"window_start" %in% names(windows) && "start" %in% names(windows)) {
    windows$window_start <- windows$start
  }
  if (!"window_end" %in% names(windows) && "end" %in% names(windows)) {
    windows$window_end <- windows$end
  }

  required <- c("window_start", "window_end")
  missing <- setdiff(required, names(windows))
  if (length(missing) > 0L) {
    cli::cli_abort("`windows` must contain {.field {missing}}.")
  }

  if (!"window_label" %in% names(windows)) {
    windows$window_label <- paste0("window_", seq_len(nrow(windows)))
  }

  windows <- windows[c("window_label", "window_start", "window_end")]
  windows$window_label <- as.character(windows$window_label)
  windows$window_start <- as.numeric(windows$window_start)
  windows$window_end <- as.numeric(windows$window_end)

  if (any(!is.finite(windows$window_start)) || any(!is.finite(windows$window_end))) {
    cli::cli_abort("Window starts and ends must be finite numbers.")
  }
  if (any(windows$window_end <= windows$window_start)) {
    cli::cli_abort("Each window end must be greater than its start.")
  }

  windows
}


#' @keywords internal
build_trajectory_comparison_table <- function(
    key,
    effects,
    derivatives,
    windows,
    view,
    metrics,
    standardize,
    inference
) {
  pairs <- comparison_pairs(key, view)
  rows <- list()
  row_i <- 0L

  for (pair_i in seq_len(nrow(pairs))) {
    unit_1 <- pairs$trajectory_id_1[[pair_i]]
    unit_2 <- pairs$trajectory_id_2[[pair_i]]
    key_1 <- key[key$trajectory_id == unit_1, , drop = FALSE]
    key_2 <- key[key$trajectory_id == unit_2, , drop = FALSE]

    for (window_i in seq_len(nrow(windows))) {
      window <- windows[window_i, , drop = FALSE]

      fitted_pair <- align_unit_values(
        data = effects,
        unit_1 = unit_1,
        unit_2 = unit_2,
        value_col = ".fitted",
        window = window
      )

      derivative_pair <- NULL
      if ("derivative_correlation" %in% metrics) {
        derivative_pair <- align_unit_values(
          data = derivatives,
          unit_1 = unit_1,
          unit_2 = unit_2,
          value_col = ".derivative",
          window = window
        )
      }

      for (metric in metrics) {
        row_i <- row_i + 1L

        pair_data <- if (metric == "derivative_correlation") {
          derivative_pair
        } else {
          fitted_pair
        }

        rows[[row_i]] <- summarize_trajectory_pair(
          key_1 = key_1,
          key_2 = key_2,
          pair_data = pair_data,
          metric = metric,
          window = window,
          standardize = standardize,
          inference = inference
        )
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


#' @keywords internal
comparison_pairs <- function(key, view) {
  if (view == "same_group_level") {
    split_key <- paste(key$group_variable, key$group_level, sep = "\r")
    pieces <- split(key, split_key, drop = TRUE)
    pair_rows <- lapply(pieces, function(piece) {
      if (nrow(piece) < 2L) {
        return(NULL)
      }
      combinations <- utils::combn(piece$trajectory_id, 2L)
      out <- data.frame(
        trajectory_id_1 = combinations[1L, ],
        trajectory_id_2 = combinations[2L, ],
        stringsAsFactors = FALSE
      )
      feature_1 <- piece$feature[match(out$trajectory_id_1, piece$trajectory_id)]
      feature_2 <- piece$feature[match(out$trajectory_id_2, piece$trajectory_id)]
      out[feature_1 != feature_2, , drop = FALSE]
    })
    pair_rows <- pair_rows[!vapply(pair_rows, is.null, logical(1))]
    if (length(pair_rows) == 0L) {
      cli::cli_abort("No comparable feature pairs were found for the requested view.")
    }
    return(do.call(rbind, pair_rows))
  }

  combinations <- utils::combn(key$trajectory_id, 2L)
  data.frame(
    trajectory_id_1 = combinations[1L, ],
    trajectory_id_2 = combinations[2L, ],
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
align_unit_values <- function(data, unit_1, unit_2, value_col, window) {
  if (nrow(data) == 0L || !value_col %in% names(data)) {
    return(data.frame(.time = numeric(), x = numeric(), y = numeric()))
  }

  data <- data[
    data$.time >= window$window_start[[1L]] &
      data$.time <= window$window_end[[1L]],
    ,
    drop = FALSE
  ]

  x <- data[data$trajectory_id == unit_1, c(".time", value_col), drop = FALSE]
  y <- data[data$trajectory_id == unit_2, c(".time", value_col), drop = FALSE]

  names(x)[names(x) == value_col] <- "x"
  names(y)[names(y) == value_col] <- "y"

  merge(x, y, by = ".time")
}


#' @keywords internal
summarize_trajectory_pair <- function(
    key_1,
    key_2,
    pair_data,
    metric,
    window,
    standardize,
    inference
) {
  flags <- character()

  if (nrow(pair_data) < 3L) {
    flags <- c(flags, "insufficient_time_points")
  }

  x <- pair_data$x
  y <- pair_data$y

  if (isTRUE(standardize)) {
    standardized <- standardize_pair_values(x, y)
    x <- standardized$x
    y <- standardized$y
    flags <- c(flags, standardized$flags)
  }

  if (metric == "derivative_correlation" && nrow(pair_data) == 0L) {
    flags <- c(flags, "derivatives_not_available")
  }

  estimate <- NA_real_
  p_value <- NA_real_
  test_method <- "none"
  null_hypothesis <- NA_character_
  warning <- NA_character_

  valid <- is.finite(x) & is.finite(y)
  x_valid <- x[valid]
  y_valid <- y[valid]

  if (metric %in% c("fitted_correlation", "derivative_correlation")) {
    null_hypothesis <- "The compared fitted trajectories have zero Pearson correlation over the selected grid."

    if (length(x_valid) >= 3L &&
        stats::sd(x_valid) > 0 &&
        stats::sd(y_valid) > 0) {
      estimate <- stats::cor(x_valid, y_valid)

      if (inference == "descriptive") {
        test <- suppressWarnings(stats::cor.test(x_valid, y_valid))
        p_value <- test$p.value
        test_method <- "curve_cor_test_descriptive"
        warning <- "Fitted grid points are not independent; p-value is descriptive."
        flags <- c(flags, "fitted_grid_points_not_independent")
      }
    }
  } else if (metric == "fitted_distance") {
    if (length(x_valid) >= 2L) {
      estimate <- sqrt(mean((x_valid - y_valid)^2))
    }
    null_hypothesis <- "No inferential null hypothesis is attached to fitted trajectory distance."
  }

  if (!is.finite(estimate)) {
    flags <- c(flags, "estimate_not_available")
  }
  if (identical(key_1$feature[[1L]], key_2$feature[[1L]])) {
    flags <- c(flags, "within_feature_group_difference")
  }
  if (!identical(key_1$group_variable[[1L]], key_2$group_variable[[1L]])) {
    flags <- c(flags, "different_grouping_variables")
  }

  flags <- unique(flags)

  data.frame(
    trajectory_id_1 = key_1$trajectory_id,
    trajectory_id_2 = key_2$trajectory_id,
    feature_1 = key_1$feature,
    feature_2 = key_2$feature,
    assay_1 = key_1$assay,
    assay_2 = key_2$assay,
    trajectory_group = comparison_pair_group(key_1, key_2),
    group_variable = comparison_pair_value(key_1$group_variable, key_2$group_variable),
    group_level = comparison_pair_value(key_1$group_level, key_2$group_level),
    window_label = window$window_label,
    window_start = window$window_start,
    window_end = window$window_end,
    metric = metric,
    estimate = estimate,
    p_value = p_value,
    q_value = NA_real_,
    test_method = test_method,
    null_hypothesis = null_hypothesis,
    inference_warning = warning,
    n_time_points = length(x_valid),
    standardized = isTRUE(standardize),
    support_flags = paste(flags, collapse = ";"),
    statistical_summary = comparison_statistical_summary(metric, estimate),
    biological_hypothesis = comparison_biological_hypothesis(metric, estimate),
    limitations = comparison_limitations(metric, inference),
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
standardize_pair_values <- function(x, y) {
  flags <- character()

  out <- list(x = x, y = y, flags = flags)

  if (length(x) == 0L || length(y) == 0L) {
    return(out)
  }

  x_sd <- stats::sd(x, na.rm = TRUE)
  y_sd <- stats::sd(y, na.rm = TRUE)

  if (is.finite(x_sd) && x_sd > 0) {
    out$x <- as.numeric(scale(x))
  } else {
    out$x <- rep(0, length(x))
    out$flags <- c(out$flags, "zero_variance_trajectory_1")
  }

  if (is.finite(y_sd) && y_sd > 0) {
    out$y <- as.numeric(scale(y))
  } else {
    out$y <- rep(0, length(y))
    out$flags <- c(out$flags, "zero_variance_trajectory_2")
  }

  out
}


#' @keywords internal
comparison_pair_group <- function(key_1, key_2) {
  if (identical(key_1$trajectory_group[[1L]], key_2$trajectory_group[[1L]])) {
    return(key_1$trajectory_group[[1L]])
  }

  paste(key_1$trajectory_group[[1L]], key_2$trajectory_group[[1L]], sep = " vs ")
}


#' @keywords internal
comparison_pair_value <- function(x, y) {
  x <- x[[1L]]
  y <- y[[1L]]

  if (identical(x, y)) {
    return(x)
  }

  paste(x, y, sep = " vs ")
}


#' @keywords internal
comparison_statistical_summary <- function(metric, estimate) {
  if (!is.finite(estimate)) {
    return("The comparison estimate is not available for this pair and window.")
  }

  if (metric == "fitted_distance") {
    return(paste0(
      "The standardized fitted trajectory distance is ",
      format(round(estimate, 3), nsmall = 3),
      "."
    ))
  }

  strength <- if (abs(estimate) >= 0.7) {
    "strong"
  } else if (abs(estimate) >= 0.3) {
    "moderate"
  } else {
    "weak"
  }

  direction <- if (estimate >= 0) "positive" else "negative"

  paste0(
    "The fitted trajectories show a ",
    strength,
    " ",
    direction,
    " correlation over this window."
  )
}


#' @keywords internal
comparison_biological_hypothesis <- function(metric, estimate) {
  if (!is.finite(estimate)) {
    return("No biological trajectory hypothesis should be drawn from this unavailable estimate.")
  }

  if (metric == "fitted_distance") {
    return("Smaller distances are consistent with more similar temporal response shapes.")
  }

  if (estimate >= 0.7) {
    return("The features may share coordinated temporal behavior under this model specification.")
  }
  if (estimate <= -0.7) {
    return("The features may show opposing temporal behavior under this model specification.")
  }

  "The fitted trajectories do not provide strong evidence of coordinated temporal behavior."
}


#' @keywords internal
comparison_limitations <- function(metric, inference) {
  if (metric == "fitted_distance") {
    return("Distance is descriptive and depends on the selected model, grid, window, and standardization.")
  }

  if (inference == "descriptive") {
    return("Correlation p-values summarize fitted grid curves, not independent biological observations.")
  }

  "No inferential p-value was requested for this comparison."
}


#' @keywords internal
build_distance_matrices <- function(table) {
  distance_rows <- table[
    table$metric == "fitted_distance" & is.finite(table$estimate),
    ,
    drop = FALSE
  ]

  if (nrow(distance_rows) == 0L) {
    return(list())
  }

  matrix_group <- paste(
    distance_rows$window_label,
    distance_rows$trajectory_group,
    sep = " | "
  )
  pieces <- split(distance_rows, matrix_group, drop = TRUE)

  matrices <- lapply(pieces, function(piece) {
    ids <- sort(unique(c(piece$trajectory_id_1, piece$trajectory_id_2)))
    mat <- matrix(NA_real_, nrow = length(ids), ncol = length(ids))
    rownames(mat) <- ids
    colnames(mat) <- ids
    diag(mat) <- 0

    for (i in seq_len(nrow(piece))) {
      id_1 <- piece$trajectory_id_1[[i]]
      id_2 <- piece$trajectory_id_2[[i]]
      mat[id_1, id_2] <- piece$estimate[[i]]
      mat[id_2, id_1] <- piece$estimate[[i]]
    }

    attr(mat, "metric") <- "fitted_distance"
    attr(mat, "trajectory_group") <- piece$trajectory_group[[1L]]
    attr(mat, "window_label") <- piece$window_label[[1L]]
    mat
  })

  names(matrices) <- make.names(names(matrices), unique = TRUE)
  matrices
}


#' @keywords internal
comparison_table_columns <- function() {
  c(
    "trajectory_id_1",
    "trajectory_id_2",
    "feature_1",
    "feature_2",
    "assay_1",
    "assay_2",
    "trajectory_group",
    "group_variable",
    "group_level",
    "window_label",
    "window_start",
    "window_end",
    "metric",
    "estimate",
    "p_value",
    "q_value",
    "test_method",
    "null_hypothesis",
    "inference_warning",
    "n_time_points",
    "standardized",
    "support_flags",
    "statistical_summary",
    "biological_hypothesis",
    "limitations"
  )
}

#####
# Effect summaries and plots
#####

#' Extract trajectory effect predictions
#'
#' Creates a prediction grid over the trajectory time variable and returns
#' fitted values with optional standard errors.
#'
#' @param fit A `levaim_fit`.
#' @param n Number of time grid points when `time_values` is not supplied.
#' @param time_values Optional numeric vector of time values.
#' @param se_fit Logical. Whether to compute standard errors.
#' @param type Prediction type passed to the backend.
#' @param level Confidence level for `.lower` and `.upper`.
#' @param ... Additional arguments passed to backend prediction.
#'
#' @return A data.frame of trajectory effects.
#' @export
trajectory_effects <- function(
    fit,
    n = 100L,
    time_values = NULL,
    se_fit = TRUE,
    type = "response",
    level = 0.95,
    ...
) {
  if (!inherits(fit, "levaim_fit")) {
    cli::cli_abort("`fit` must be created with `fit_trajectory()`.")
  }

  if (!fit$engine %in% c("spline", "gp")) {
    cli::cli_abort("Trajectory effects are currently implemented for spline and GP fits.")
  }

  mf <- fit$model_frame
  spec <- mf$trajectory$spec
  design <- mf$trajectory$design
  time_col <- spec$time$column

  if (is.null(time_values)) {
    observed_time <- mf$data[[time_col]]
    time_values <- seq(
      min(observed_time, na.rm = TRUE),
      max(observed_time, na.rm = TRUE),
      length.out = as.integer(n)
    )
  }

  grid <- effect_reference_grid(mf, time_values)

  if (design$type == "varying") {
    grid[design$by] <- NULL

    by_values <- lapply(design$by, function(column) {
      values <- mf$data[[column]]
      if (is.factor(values)) {
        levels(values)
      } else {
        sort(unique(values))
      }
    })
    names(by_values) <- design$by

    by_grid <- expand.grid(by_values, KEEP.OUT.ATTRS = FALSE)
    grid <- merge(grid, by_grid, all = TRUE)

    for (column in design$by) {
      if (is.factor(mf$data[[column]])) {
        grid[[column]] <- factor(grid[[column]], levels = levels(mf$data[[column]]))
      }
    }
  }

  if (fit$engine == "gp") {
    pred <- stats::fitted(
      fit$fit,
      newdata = grid,
      ...
    )
  } else {
    pred <- stats::predict(
      fit$fit,
      newdata = grid,
      type = type,
      se.fit = se_fit,
      exclude = spline_random_effect_terms(mf),
      ...
    )
  }

  if (se_fit) {
    if (fit$engine == "gp") {
      grid$.fitted <- as.numeric(pred[, "Estimate"])
      grid$.se <- as.numeric(pred[, "Est.Error"])
      grid$.lower <- as.numeric(pred[, "Q2.5"])
      grid$.upper <- as.numeric(pred[, "Q97.5"])
    } else {
      grid$.fitted <- as.numeric(pred$fit)
      grid$.se <- as.numeric(pred$se.fit)
      z <- stats::qnorm((1 + level) / 2)
      grid$.lower <- grid$.fitted - z * grid$.se
      grid$.upper <- grid$.fitted + z * grid$.se
    }
  } else {
    if (fit$engine == "gp") {
      grid$.fitted <- as.numeric(pred[, "Estimate"])
    } else {
      grid$.fitted <- as.numeric(pred)
    }
  }

  grid
}


#' @keywords internal
effect_reference_grid <- function(mf, time_values) {
  data <- mf$data
  spec <- mf$trajectory$spec
  time_col <- spec$time$column

  reference <- lapply(data, reference_value)
  reference$.y <- NULL
  reference[[time_col]] <- time_values

  grid <- as.data.frame(reference, stringsAsFactors = FALSE)

  for (column in names(data)) {
    if (column %in% names(grid) && is.factor(data[[column]])) {
      grid[[column]] <- factor(grid[[column]], levels = levels(data[[column]]))
    }
  }

  grid
}


#' @keywords internal
reference_value <- function(x) {
  if (is.factor(x)) {
    return(factor(levels(x)[[1]], levels = levels(x)))
  }

  if (is.numeric(x)) {
    return(stats::median(x, na.rm = TRUE))
  }

  if (is.logical(x)) {
    return(stats::median(as.numeric(x), na.rm = TRUE) > 0)
  }

  unique(x)[[1]]
}


#' Plot trajectory effects
#'
#' @param fit A `levaim_fit`.
#' @param effects Optional output from `trajectory_effects()`.
#' @param ... Additional arguments passed to `trajectory_effects()`.
#'
#' @return Invisibly returns the plotted effects data.
#' @importFrom graphics legend lines
#' @export
plot_trajectory_effects <- function(fit, effects = NULL, ...) {
  if (is.null(effects)) {
    effects <- trajectory_effects(fit, ...)
  }

  mf <- fit$model_frame
  time_col <- mf$trajectory$spec$time$column
  by_cols <- if (mf$trajectory$design$type == "varying") {
    mf$trajectory$design$by
  } else {
    character()
  }

  if (length(by_cols) == 0L) {
    plot(
      effects[[time_col]],
      effects$.fitted,
      type = "l",
      xlab = time_col,
      ylab = "Fitted response"
    )
  } else {
    group <- interaction(effects[by_cols], drop = TRUE)
    groups <- levels(group)
    colors <- seq_along(groups)

    plot(
      effects[[time_col]],
      effects$.fitted,
      type = "n",
      xlab = time_col,
      ylab = "Fitted response"
    )

    for (i in seq_along(groups)) {
      keep <- group == groups[[i]]
      lines(effects[[time_col]][keep], effects$.fitted[keep], col = colors[[i]])
    }

    legend("topright", legend = groups, col = colors, lty = 1, bty = "n")
  }

  invisible(effects)
}


#' Extract pairwise trajectory contrasts over time
#'
#' Computes covariance-aware pairwise differences between fitted trajectory
#' groups on a common time grid. Pointwise p-values localize separation but do
#' not replace an omnibus trajectory test.
#'
#' @param fit A varying-trajectory spline `levaim_fit`.
#' @param n Number of grid points when `time_values` is not supplied.
#' @param time_values Optional numeric time grid.
#' @param levels Optional character vector of trajectory-group levels to compare.
#' @param level Confidence level for pointwise intervals.
#'
#' @return A data.frame with one row per time and group-pair contrast.
#' @export
trajectory_contrasts <- function(
    fit,
    n = 100L,
    time_values = NULL,
    levels = NULL,
    level = 0.95
) {
  if (!inherits(fit, "levaim_fit") || fit$engine != "spline") {
    cli::cli_abort("`fit` must be a fitted spline trajectory.")
  }

  mf <- fit$model_frame
  design <- mf$trajectory$design
  if (design$type != "varying" || length(design$by) != 1L) {
    cli::cli_abort("Trajectory contrasts currently require one varying-trajectory modifier.")
  }

  time_col <- mf$trajectory$spec$time$column
  group_col <- design$by[[1L]]
  observed_levels <- if (is.factor(mf$data[[group_col]])) {
    levels(mf$data[[group_col]])
  } else {
    sort(unique(as.character(mf$data[[group_col]])))
  }
  levels <- levels %||% observed_levels

  if (!is.character(levels) || length(levels) < 2L) {
    cli::cli_abort("`levels` must select at least two trajectory-group levels.")
  }
  missing_levels <- setdiff(levels, observed_levels)
  if (length(missing_levels) > 0L) {
    cli::cli_abort(
      "Unknown trajectory level(s): {paste(missing_levels, collapse = ', ')}"
    )
  }

  if (is.null(time_values)) {
    observed_time <- mf$data[[time_col]]
    time_values <- seq(
      min(observed_time, na.rm = TRUE),
      max(observed_time, na.rm = TRUE),
      length.out = as.integer(n)
    )
  }
  time_values <- sort(unique(as.numeric(time_values)))

  grid <- effect_reference_grid(mf, time_values)
  grid[[group_col]] <- NULL
  group_grid <- data.frame(.group = levels, stringsAsFactors = FALSE)
  names(group_grid) <- group_col
  grid <- merge(grid, group_grid, all = TRUE)
  if (is.factor(mf$data[[group_col]])) {
    grid[[group_col]] <- factor(grid[[group_col]], levels = observed_levels)
  }
  grid <- grid[order(grid[[time_col]], grid[[group_col]]), , drop = FALSE]

  exclude <- spline_random_effect_terms(mf)
  design_matrix <- stats::predict(
    fit$fit,
    newdata = grid,
    type = "lpmatrix",
    exclude = exclude
  )
  coefficients <- stats::coef(fit$fit)
  covariance <- stats::vcov(fit$fit)
  z_multiplier <- stats::qnorm((1 + level) / 2)

  pairs <- utils::combn(levels, 2L, simplify = FALSE)
  rows <- vector("list", length(pairs) * length(time_values))
  contrast_rows <- vector("list", length(rows))
  row_id <- 1L

  for (time_value in time_values) {
    for (pair in pairs) {
      first <- which(grid[[time_col]] == time_value & as.character(grid[[group_col]]) == pair[[1L]])
      second <- which(grid[[time_col]] == time_value & as.character(grid[[group_col]]) == pair[[2L]])
      contrast <- design_matrix[first, ] - design_matrix[second, ]
      estimate <- as.numeric(contrast %*% coefficients)
      se <- sqrt(max(0, as.numeric(contrast %*% covariance %*% contrast)))
      statistic <- if (is.finite(se) && se > 0) estimate / se else NA_real_

      rows[[row_id]] <- data.frame(
        time = time_value,
        group_1 = pair[[1L]],
        group_2 = pair[[2L]],
        estimate = estimate,
        std_error = se,
        lower = estimate - z_multiplier * se,
        upper = estimate + z_multiplier * se,
        z_value = statistic,
        p_value = if (is.finite(statistic)) {
          2 * stats::pnorm(-abs(statistic))
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
      contrast_rows[[row_id]] <- as.numeric(contrast)
      row_id <- row_id + 1L
    }
  }

  out <- do.call(rbind, rows)
  names(out)[names(out) == "time"] <- time_col
  attr(out, "contrast_matrix") <- do.call(rbind, contrast_rows)
  attr(out, "coefficient_covariance") <- covariance
  attr(out, "group_variable") <- group_col
  attr(out, "time_variable") <- time_col
  attr(out, "tested_levels") <- levels
  out
}


#' Estimate trajectory derivatives
#'
#' Estimates local rates of change using finite differences. Spline derivatives
#' use fitted effects. GP derivatives are calculated for each posterior fitted
#' trajectory and summarized with pointwise credible intervals and posterior
#' sign probabilities. This function reports derivative estimates; it does not
#' infer peaks, valleys, or biological shape classes.
#'
#' @param fit A `levaim_fit`.
#' @param n Number of time grid points when `time_values` is not supplied.
#' @param time_values Optional numeric vector of time values.
#' @param method Derivative method. Currently only `"finite_difference"`.
#' @param level Confidence level passed to [trajectory_effects()].
#' @param type Prediction type passed to the backend.
#' @param draws Maximum number of posterior draws used for GP derivatives.
#' @param seed Sampling seed when GP posterior draws are thinned.
#' @param ... Additional arguments passed to [trajectory_effects()].
#'
#' @return A data.frame with fitted trajectory effects and derivative columns.
#' @export
trajectory_derivatives <- function(
    fit,
    n = 100L,
    time_values = NULL,
    method = c("finite_difference"),
    level = 0.95,
    type = "response",
    draws = 500L,
    seed = 1L,
    ...
) {
  method <- match.arg(method)

  if (!inherits(fit, "levaim_fit")) {
    cli::cli_abort("`fit` must be created with `fit_trajectory()`.")
  }

  if (!fit$engine %in% c("spline", "gp")) {
    cli::cli_abort("Trajectory derivatives require a spline or GP fit.")
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    cli::cli_abort("`level` must lie strictly between 0 and 1.")
  }
  draws <- as.integer(draws)
  if (length(draws) != 1L || is.na(draws) || draws < 20L) {
    cli::cli_abort("`draws` must be an integer >= 20.")
  }

  effects <- trajectory_effects(
    fit,
    n = n,
    time_values = time_values,
    se_fit = TRUE,
    type = type,
    level = level,
    ...
  )

  mf <- fit$model_frame
  time_col <- mf$trajectory$spec$time$column
  by_cols <- if (mf$trajectory$design$type == "varying") {
    mf$trajectory$design$by
  } else {
    character()
  }

  if (fit$engine == "gp") {
    return(posterior_derivative_effects(
      fit = fit,
      effects = effects,
      time_col = time_col,
      by_cols = by_cols,
      draws = draws,
      level = level,
      seed = seed
    ))
  }

  out <- finite_difference_effects(
    effects = effects,
    time_col = time_col,
    by_cols = by_cols
  )
  out$.probability_positive <- NA_real_
  out$.derivative_method <- "fitted_finite_difference"
  out$.draws <- NA_integer_
  out
}


#' @keywords internal
posterior_derivative_effects <- function(
    fit,
    effects,
    time_col,
    by_cols = character(),
    draws = 500L,
    level = 0.95,
    seed = 1L
) {
  effects$.trajectory_group <- if (length(by_cols) == 0L) {
    factor("shared")
  } else {
    interaction(effects[by_cols], drop = TRUE)
  }
  fitted_draws <- trajectory_fitted_draws(fit, effects, draws, seed)
  probabilities <- c((1 - level) / 2, 1 - (1 - level) / 2)

  pieces <- lapply(split(seq_len(nrow(effects)), effects$.trajectory_group), function(indices) {
    ordering <- order(effects[[time_col]][indices])
    indices <- indices[ordering]
    piece <- effects[indices, , drop = FALSE]
    time <- piece[[time_col]]
    group_draws <- fitted_draws[, indices, drop = FALSE]
    derivative_draws <- apply(group_draws, 1L, function(values) {
      finite_difference_vector(time, values)
    })
    if (is.null(dim(derivative_draws))) {
      derivative_draws <- matrix(derivative_draws, ncol = 1L)
    }

    piece$.derivative <- rowMeans(derivative_draws)
    piece$.lower_derivative <- apply(
      derivative_draws,
      1L,
      stats::quantile,
      probs = probabilities[[1L]],
      na.rm = TRUE
    )
    piece$.upper_derivative <- apply(
      derivative_draws,
      1L,
      stats::quantile,
      probs = probabilities[[2L]],
      na.rm = TRUE
    )
    piece$.probability_positive <- rowMeans(derivative_draws > 0, na.rm = TRUE)
    piece$.sign <- classify_derivative_sign(
      piece$.derivative,
      lower = piece$.lower_derivative,
      upper = piece$.upper_derivative
    )
    piece$.derivative_method <- "posterior_finite_difference"
    piece$.draws <- ncol(derivative_draws)
    piece
  })

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out$.trajectory_group <- NULL
  out
}


#' @keywords internal
finite_difference_effects <- function(effects, time_col, by_cols = character()) {
  if (!time_col %in% names(effects)) {
    cli::cli_abort("Time column {.field {time_col}} was not found in effects data.")
  }

  if (!".fitted" %in% names(effects)) {
    cli::cli_abort("Effects data must contain a `.fitted` column.")
  }

  if (length(by_cols) == 0L) {
    effects$.trajectory_group <- factor("shared")
  } else {
    effects$.trajectory_group <- interaction(effects[by_cols], drop = TRUE)
  }

  pieces <- lapply(
    split(effects, effects$.trajectory_group),
    function(piece) {
      piece <- piece[order(piece[[time_col]]), , drop = FALSE]
      piece$.derivative <- finite_difference_vector(piece[[time_col]], piece$.fitted)

      if (all(c(".lower", ".upper") %in% names(piece))) {
        lower_derivative <- finite_difference_vector(piece[[time_col]], piece$.lower)
        upper_derivative <- finite_difference_vector(piece[[time_col]], piece$.upper)
        piece$.lower_derivative <- pmin(lower_derivative, upper_derivative)
        piece$.upper_derivative <- pmax(lower_derivative, upper_derivative)
        piece$.sign <- classify_derivative_sign(
          piece$.derivative,
          lower = piece$.lower_derivative,
          upper = piece$.upper_derivative
        )
      } else {
        piece$.sign <- classify_derivative_sign(piece$.derivative)
      }

      piece
    }
  )

  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out$.trajectory_group <- NULL
  out
}


#' @keywords internal
finite_difference_vector <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)

  if (length(x) != length(y)) {
    cli::cli_abort("`x` and `y` must have the same length.")
  }

  if (length(x) < 2L) {
    cli::cli_abort("At least two time values are required to estimate derivatives.")
  }

  if (any(diff(x) <= 0)) {
    cli::cli_abort("Time values must be strictly increasing within each trajectory.")
  }

  derivative <- numeric(length(x))
  derivative[[1]] <- (y[[2]] - y[[1]]) / (x[[2]] - x[[1]])

  last <- length(x)
  derivative[[last]] <- (y[[last]] - y[[last - 1L]]) / (x[[last]] - x[[last - 1L]])

  if (length(x) > 2L) {
    for (i in 2:(length(x) - 1L)) {
      derivative[[i]] <- (y[[i + 1L]] - y[[i - 1L]]) / (x[[i + 1L]] - x[[i - 1L]])
    }
  }

  derivative
}


#' @keywords internal
classify_derivative_sign <- function(derivative, lower = NULL, upper = NULL) {
  if (!is.null(lower) && !is.null(upper)) {
    sign <- ifelse(
      lower > 0,
      "increasing",
      ifelse(upper < 0, "decreasing", "uncertain")
    )
  } else {
    sign <- ifelse(
      derivative > 0,
      "increasing",
      ifelse(derivative < 0, "decreasing", "flat")
    )
  }

  sign[!is.finite(derivative)] <- "uncertain"
  sign
}

#####
# Trajectory moment summaries
#####

#' Summarize fitted trajectories with functional moments
#'
#' Produces a table-first summary of fitted trajectory level, timing, variation,
#' and rates of change over declared time windows. Spline uncertainty is drawn
#' from the fitted coefficient covariance; GP uncertainty uses posterior fitted
#' draws.
#'
#' @param fit A `levaim_fit`, a list of fits, or `levaim_feature_tests` created
#'   with `keep_fits = TRUE`.
#' @param n Number of grid points when `time_values` is not supplied.
#' @param time_values Optional numeric time grid.
#' @param windows Optional data.frame with `window_label`, `window_start`, and
#'   `window_end`, or `label`, `start`, and `end`.
#' @param levels Optional trajectory-group levels to summarize. Other fitted
#'   levels remain part of the model but are omitted from the report.
#' @param draws Number of uncertainty draws.
#' @param level Interval level.
#' @param seed Simulation seed.
#'
#' @return A `levaim_trajectory_moments` data.frame with one row per feature,
#'   trajectory group, window, and moment.
#' @export
trajectory_moments <- function(
    fit,
    n = 200L,
    time_values = NULL,
    windows = NULL,
    levels = NULL,
    draws = 500L,
    level = 0.95,
    seed = 1L
) {
  fits <- normalize_comparison_fits(fit)
  draws <- as.integer(draws)
  if (length(draws) != 1L || is.na(draws) || draws < 20L) {
    cli::cli_abort("`draws` must be an integer >= 20.")
  }
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    cli::cli_abort("`level` must lie strictly between 0 and 1.")
  }

  pieces <- lapply(seq_along(fits), function(i) {
    trajectory_moments_one(
      fit = fits[[i]],
      feature_name = names(fits)[[i]],
      n = n,
      time_values = time_values,
      windows = windows,
      levels = levels,
      draws = draws,
      level = level,
      seed = seed + i - 1L
    )
  })

  table <- do.call(rbind, lapply(pieces, `[[`, "table"))
  moment_draws <- do.call(rbind, lapply(pieces, `[[`, "draws"))
  rownames(table) <- NULL
  rownames(moment_draws) <- NULL
  attr(table, "moment_draws") <- moment_draws
  attr(table, "settings") <- list(
    n = n,
    time_values = time_values,
    windows = windows,
    levels = levels,
    draws = draws,
    level = level,
    seed = seed
  )
  class(table) <- c("levaim_trajectory_moments", class(table))
  table
}


#' Compare functional trajectory moments between groups
#'
#' @param moments Output from [trajectory_moments()].
#' @param reference Optional reference trajectory-group level. When supplied,
#'   only contrasts against that level are returned.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#'
#' @return A `levaim_moment_comparisons` data.frame.
#' @export
compare_trajectory_moments <- function(
    moments,
    reference = NULL,
    p_adjust_method = "BH"
) {
  if (!inherits(moments, "levaim_trajectory_moments")) {
    cli::cli_abort("`moments` must be created with `trajectory_moments()`.")
  }
  draw_matrix <- attr(moments, "moment_draws")
  if (is.null(draw_matrix) || nrow(draw_matrix) != nrow(moments)) {
    cli::cli_abort("Moment uncertainty draws are unavailable.")
  }

  split_columns <- c("feature", "window_label", "moment")
  split_key <- interaction(moments[split_columns], drop = TRUE, lex.order = TRUE)
  pieces <- split(seq_len(nrow(moments)), split_key)
  rows <- list()
  row_id <- 1L
  interval_level <- attr(moments, "settings")$level
  probabilities <- c((1 - interval_level) / 2, 1 - (1 - interval_level) / 2)

  for (indices in pieces) {
    group_levels <- moments$group_level[indices]
    if (length(unique(group_levels)) < 2L) {
      next
    }
    pairs <- utils::combn(indices, 2L, simplify = FALSE)
    for (pair in pairs) {
      first <- pair[[1L]]
      second <- pair[[2L]]
      if (!is.null(reference)) {
        if (!reference %in% moments$group_level[pair]) {
          next
        }
        if (moments$group_level[[first]] == reference) {
          pair <- rev(pair)
          first <- pair[[1L]]
          second <- pair[[2L]]
        }
      }

      difference_draws <- draw_matrix[first, ] - draw_matrix[second, ]
      finite <- is.finite(difference_draws)
      difference_draws <- difference_draws[finite]
      if (length(difference_draws) == 0L) {
        lower <- upper <- p_value <- probability_positive <- NA_real_
        status <- "not_estimable"
      } else {
        interval <- stats::quantile(difference_draws, probabilities, na.rm = TRUE)
        lower <- interval[[1L]]
        upper <- interval[[2L]]
        probability_positive <- mean(difference_draws > 0)
        p_value <- min(1, 2 * min(probability_positive, 1 - probability_positive))
        status <- "estimated"
      }

      rows[[row_id]] <- data.frame(
        feature = moments$feature[[first]],
        response = moments$response[[first]],
        group_variable = moments$group_variable[[first]],
        group_1 = moments$group_level[[first]],
        group_2 = moments$group_level[[second]],
        window_label = moments$window_label[[first]],
        window_start = moments$window_start[[first]],
        window_end = moments$window_end[[first]],
        moment = moments$moment[[first]],
        difference = moments$estimate[[first]] - moments$estimate[[second]],
        lower = lower,
        upper = upper,
        probability_positive = probability_positive,
        p_value = p_value,
        q_value = NA_real_,
        test_method = if (moments$engine[[first]] == "gp") {
          "posterior_tail_probability"
        } else {
          "coefficient_simulation_tail_probability"
        },
        inference_warning = if (moments$engine[[first]] == "gp") {
          "Posterior tail probability is not a frequentist sampling p-value."
        } else {
          "Simulation-tail p-value uses the fitted spline coefficient covariance approximation."
        },
        status = status,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
    }
  }

  if (length(rows) == 0L) {
    cli::cli_abort("No between-group moment comparisons are available.")
  }
  out <- do.call(rbind, rows)
  valid <- is.finite(out$p_value)
  out$q_value[valid] <- stats::p.adjust(out$p_value[valid], method = p_adjust_method)
  rownames(out) <- NULL
  class(out) <- c("levaim_moment_comparisons", class(out))
  out
}


#' @keywords internal
trajectory_moments_one <- function(
    fit,
    feature_name,
    n,
    time_values,
    windows,
    levels,
    draws,
    level,
    seed
) {
  effects <- trajectory_effects(
    fit,
    n = n,
    time_values = time_values,
    se_fit = TRUE
  )
  mf <- fit$model_frame
  time_col <- mf$trajectory$spec$time$column
  by_cols <- if (mf$trajectory$design$type == "varying") {
    mf$trajectory$design$by
  } else {
    character()
  }
  effects$.moment_group <- if (length(by_cols) == 0L) {
    "shared"
  } else {
    as.character(interaction(effects[by_cols], drop = TRUE, sep = ":"))
  }

  if (!is.null(levels)) {
    missing_levels <- setdiff(levels, unique(effects$.moment_group))
    if (length(missing_levels) > 0L) {
      cli::cli_abort("Unknown moment-report level(s): {paste(missing_levels, collapse = ', ')}")
    }
    effects <- effects[effects$.moment_group %in% levels, , drop = FALSE]
  }

  observed_time <- range(effects[[time_col]], na.rm = TRUE)
  windows <- normalize_moment_windows(windows, observed_time)
  fitted_draws <- trajectory_fitted_draws(fit, effects, draws, seed)
  response <- format_response_label(mf$trajectory$spec$response)
  feature <- if (mf$trajectory$spec$response$source == "assay") {
    mf$trajectory$spec$response$feature
  } else {
    feature_name
  }

  rows <- list()
  draw_rows <- list()
  row_id <- 1L
  probabilities <- c((1 - level) / 2, 1 - (1 - level) / 2)

  for (group_level in unique(effects$.moment_group)) {
    group_indices <- which(effects$.moment_group == group_level)
    group_effects <- effects[group_indices, , drop = FALSE]
    group_draws <- fitted_draws[, group_indices, drop = FALSE]

    for (window_id in seq_len(nrow(windows))) {
      window <- windows[window_id, , drop = FALSE]
      keep <- group_effects[[time_col]] >= window$window_start &
        group_effects[[time_col]] <= window$window_end
      time <- group_effects[[time_col]][keep]
      values <- group_effects$.fitted[keep]
      window_draws <- group_draws[, keep, drop = FALSE]

      if (length(time) < 3L) {
        cli::cli_abort(
          "Moment window {.val {window$window_label}} contains fewer than three grid points."
        )
      }
      ordering <- order(time)
      time <- time[ordering]
      values <- values[ordering]
      window_draws <- window_draws[, ordering, drop = FALSE]

      point <- calculate_trajectory_moments(time, values)
      simulated <- t(apply(window_draws, 1L, function(x) {
        calculate_trajectory_moments(time, x)$values
      }))
      colnames(simulated) <- names(point$values)

      for (moment_name in names(point$values)) {
        moment_draw <- simulated[, moment_name]
        interval <- if (all(!is.finite(moment_draw))) {
          c(NA_real_, NA_real_)
        } else {
          stats::quantile(moment_draw, probabilities, na.rm = TRUE)
        }
        rows[[row_id]] <- data.frame(
          feature = feature,
          response = response,
          engine = fit$engine,
          group_variable = if (length(by_cols) == 0L) "shared" else paste(by_cols, collapse = ":"),
          group_level = group_level,
          window_label = window$window_label,
          window_start = window$window_start,
          window_end = window$window_end,
          moment = moment_name,
          estimate = unname(point$values[[moment_name]]),
          lower = interval[[1L]],
          upper = interval[[2L]],
          status = point$status[[moment_name]],
          n_grid = length(time),
          stringsAsFactors = FALSE
        )
        draw_rows[[row_id]] <- moment_draw
        row_id <- row_id + 1L
      }
    }
  }

  list(table = do.call(rbind, rows), draws = do.call(rbind, draw_rows))
}


#' @keywords internal
trajectory_fitted_draws <- function(fit, effects, draws, seed) {
  summary_columns <- c(
    ".fitted", ".se", ".lower", ".upper", ".moment_group",
    ".trajectory_group"
  )
  grid <- effects[, setdiff(names(effects), summary_columns), drop = FALSE]

  if (fit$engine == "gp") {
    posterior <- stats::fitted(
      fit$fit,
      newdata = grid,
      summary = FALSE,
      re_formula = NA
    )
    if (nrow(posterior) > draws) {
      indices <- moment_sample_indices(nrow(posterior), draws, seed)
      posterior <- posterior[indices, , drop = FALSE]
    }
    return(posterior)
  }

  design_matrix <- stats::predict(
    fit$fit,
    newdata = grid,
    type = "lpmatrix",
    exclude = spline_random_effect_terms(fit$model_frame)
  )
  coefficients <- stats::coef(fit$fit)
  covariance <- stats::vcov(fit$fit)
  coefficient_draws <- multivariate_normal_draws(
    coefficients,
    covariance,
    draws,
    seed
  )
  t(design_matrix %*% coefficient_draws)
}


#' @keywords internal
multivariate_normal_draws <- function(mean, covariance, draws, seed) {
  decomposition <- eigen(covariance, symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  old_seed <- preserve_random_seed()
  on.exit(restore_random_seed(old_seed), add = TRUE)
  set.seed(seed)
  standard <- matrix(stats::rnorm(length(mean) * draws), nrow = length(mean))
  sweep(
    decomposition$vectors %*% (sqrt(values) * standard),
    1L,
    mean,
    "+"
  )
}


#' @keywords internal
moment_sample_indices <- function(n, size, seed) {
  old_seed <- preserve_random_seed()
  on.exit(restore_random_seed(old_seed), add = TRUE)
  set.seed(seed)
  sample.int(n, size)
}


#' @keywords internal
preserve_random_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}


#' @keywords internal
restore_random_seed <- function(seed) {
  if (is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", seed, envir = .GlobalEnv)
  }
}


#' @keywords internal
normalize_moment_windows <- function(windows, observed_time) {
  if (is.null(windows)) {
    return(data.frame(
      window_label = "overall",
      window_start = observed_time[[1L]],
      window_end = observed_time[[2L]],
      stringsAsFactors = FALSE
    ))
  }
  windows <- as.data.frame(windows, stringsAsFactors = FALSE)
  aliases <- c(label = "window_label", start = "window_start", end = "window_end")
  for (alias in names(aliases)) {
    if (alias %in% names(windows) && !aliases[[alias]] %in% names(windows)) {
      names(windows)[names(windows) == alias] <- aliases[[alias]]
    }
  }
  required <- c("window_label", "window_start", "window_end")
  if (!all(required %in% names(windows))) {
    cli::cli_abort("`windows` must contain label/start/end columns.")
  }
  if (any(!is.finite(windows$window_start)) || any(!is.finite(windows$window_end)) ||
      any(windows$window_start >= windows$window_end)) {
    cli::cli_abort("Every moment window must have finite `start < end` bounds.")
  }
  windows[, required, drop = FALSE]
}


#' @keywords internal
calculate_trajectory_moments <- function(time, values) {
  duration <- max(time) - min(time)
  auc <- trapezoid_integral(time, values)
  mean_level <- auc / duration
  centered_variance <- trapezoid_integral(time, (values - mean_level)^2) / duration
  derivative <- finite_difference_vector(time, values)
  amplitude <- max(values) - min(values)
  peak <- which.max(values)
  valley <- which.min(values)
  center_denominator <- trapezoid_integral(time, pmax(values, 0))
  center_of_mass <- if (is.finite(center_denominator) && center_denominator > 0) {
    trapezoid_integral(time, time * pmax(values, 0)) / center_denominator
  } else {
    NA_real_
  }
  trend <- stats::coef(stats::lm(values ~ time))[[2L]]
  flat <- !is.finite(amplitude) || amplitude <= sqrt(.Machine$double.eps)
  peak_status <- turning_point_status(peak, length(values), flat)
  valley_status <- turning_point_status(valley, length(values), flat)

  moment_values <- c(
    mean_level = mean_level,
    auc = auc,
    temporal_variance = centered_variance,
    amplitude = amplitude,
    center_of_mass = center_of_mass,
    linear_trend = trend,
    mean_velocity = (values[[length(values)]] - values[[1L]]) / duration,
    velocity_variance = stats::var(derivative),
    max_growth_rate = max(derivative),
    max_growth_time = time[[which.max(derivative)]],
    max_decline_rate = min(derivative),
    max_decline_time = time[[which.min(derivative)]],
    start_level = values[[1L]],
    end_level = values[[length(values)]],
    net_change = values[[length(values)]] - values[[1L]],
    peak_value = values[[peak]],
    peak_time = time[[peak]],
    valley_value = values[[valley]],
    valley_time = time[[valley]]
  )
  status <- rep("estimated", length(moment_values))
  names(status) <- names(moment_values)
  status[c("peak_value", "peak_time")] <- peak_status
  status[c("valley_value", "valley_time")] <- valley_status
  if (!is.finite(center_of_mass)) {
    status[["center_of_mass"]] <- "not_estimable"
  }

  list(values = moment_values, status = status)
}


#' @keywords internal
trapezoid_integral <- function(x, y) {
  last <- length(y)
  sum(diff(x) * (y[-last] + y[-1L]) / 2)
}


#' @keywords internal
turning_point_status <- function(index, n, flat) {
  if (flat) {
    return("flat_or_uncertain")
  }
  if (index %in% c(1L, n)) {
    return("boundary_only")
  }
  "interior_supported"
}

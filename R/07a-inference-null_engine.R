#####
# Design-aware trajectory null engine
#####

#' Configure trajectory null estimation
#'
#' Defines how centered-shape and integrated-level trajectory tests generate null
#' distribution. The robust Gaussian default uses subject-cluster wild
#' bootstrap residuals from a hypothesis-specific reduced longitudinal model.
#' Non-Gaussian responses use a conditional parametric bootstrap. The Gaussian
#' coefficient option is a fast approximation and is labeled as such in result
#' tables.
#'
#' @param method Null-generation method. `"auto"` selects
#'   `"cluster_wild_bootstrap"` for Gaussian responses and
#'   `"parametric_bootstrap"` otherwise. `"subject_label_permutation"`
#'   permutes a subject-static focal label while preserving complete subject
#'   trajectories and is available for integrated-level hypotheses.
#' @param simulations Number of null replicates.
#' @param seed Random seed used to construct the shared resampling plan.
#' @param cores Number of worker processes. Values greater than one use forked
#'   workers on non-Windows systems and fall back to serial execution on
#'   Windows.
#' @param store_null Logical. Whether to retain the complete null statistic and
#'   effect distributions.
#'
#' @return A `levaim_trajectory_null_control` object.
#' @export
trajectory_null_control <- function(
    method = c(
      "auto", "cluster_wild_bootstrap", "parametric_bootstrap",
      "subject_label_permutation", "gaussian_approximation"
    ),
    simulations = 999L,
    seed = 1L,
    cores = 1L,
    store_null = FALSE
) {
  method <- match.arg(method)
  simulations <- as.integer(simulations)
  seed <- as.integer(seed)
  cores <- as.integer(cores)

  if (length(simulations) != 1L || is.na(simulations) || simulations < 19L) {
    cli::cli_abort("`simulations` must be an integer >= 19.")
  }
  if (simulations < 199L) {
    cli::cli_warn(
      "Fewer than 199 null replicates provide coarse p-value resolution and should be used only for tests or development."
    )
  }
  if (length(seed) != 1L || is.na(seed)) {
    cli::cli_abort("`seed` must be one non-missing integer.")
  }
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    cli::cli_abort("`cores` must be a positive integer.")
  }
  if (!is.logical(store_null) || length(store_null) != 1L || is.na(store_null)) {
    cli::cli_abort("`store_null` must be `TRUE` or `FALSE`.")
  }

  structure(
    list(
      method = method,
      simulations = simulations,
      seed = seed,
      cores = cores,
      store_null = store_null
    ),
    class = "levaim_trajectory_null_control"
  )
}


#' @keywords internal
resolve_trajectory_null_control <- function(null_control, family, seed) {
  if (is.null(null_control)) {
    null_control <- trajectory_null_control(
      method = "auto",
      seed = seed
    )
  }
  if (!inherits(null_control, "levaim_trajectory_null_control")) {
    cli::cli_abort("`null_control` must be created with `trajectory_null_control()`.")
  }

  if (null_control$method == "auto") {
    null_control$method <- if (family == "gaussian") {
      "cluster_wild_bootstrap"
    } else {
      "parametric_bootstrap"
    }
  }
  if (null_control$method == "cluster_wild_bootstrap" && family != "gaussian") {
    cli::cli_abort(
      "`cluster_wild_bootstrap` currently supports Gaussian responses only. Use `method = 'parametric_bootstrap'`."
    )
  }

  null_control
}


#' @keywords internal
shape_contrast_components <- function(
    fit,
    levels = NULL,
    time_values = NULL,
    center_time_values = NULL
) {
  raw <- trajectory_contrasts(
    fit,
    n = 100L,
    time_values = time_values,
    levels = levels
  )
  contrast_matrix <- attr(raw, "contrast_matrix")
  centering <- trajectory_contrasts(
    fit,
    n = 100L,
    time_values = center_time_values,
    levels = levels
  )
  centering_matrix <- attr(centering, "contrast_matrix")
  time_col <- attr(raw, "time_variable")
  pair_key <- paste(raw$group_1, raw$group_2, sep = "\r")
  centering_pair_key <- paste(centering$group_1, centering$group_2, sep = "\r")
  centered_matrix <- contrast_matrix

  for (pair in unique(pair_key)) {
    rows <- which(pair_key == pair)
    center_rows <- which(centering_pair_key == pair)
    if (length(center_rows) == 0L) {
      cli::cli_abort("A requested trajectory pair is absent from the centering grid.")
    }
    weights <- trapezoid_window_weights(centering[[time_col]][center_rows])
    window_mean_contrast <- as.numeric(
      crossprod(weights, centering_matrix[center_rows, , drop = FALSE])
    )
    centered_matrix[rows, ] <- sweep(
      contrast_matrix[rows, , drop = FALSE],
      2L,
      window_mean_contrast,
      FUN = "-"
    )
  }

  coefficients <- stats::coef(fit$fit)
  covariance <- stats::vcov(fit$fit)
  estimates <- as.numeric(centered_matrix %*% coefficients)
  contrast_covariance <- centered_matrix %*% covariance %*% t(centered_matrix)
  standard_errors <- sqrt(pmax(0, diag(contrast_covariance)))
  usable <- is.finite(estimates) & is.finite(standard_errors) & standard_errors > 0
  z_value <- rep(NA_real_, length(estimates))
  z_value[usable] <- estimates[usable] / standard_errors[usable]

  out <- raw
  out$raw_estimate <- out$estimate
  out$estimate <- estimates
  out$std_error <- standard_errors
  out$z_value <- z_value
  out$p_value <- 2 * stats::pnorm(-abs(z_value))
  z_multiplier <- stats::qnorm(0.975)
  out$lower <- estimates - z_multiplier * standard_errors
  out$upper <- estimates + z_multiplier * standard_errors
  attr(out, "contrast_matrix") <- centered_matrix
  attr(out, "coefficient_covariance") <- covariance
  attr(out, "contrast_scale") <- "centered_shape"
  attr(out, "centering_time_values") <- centering[[time_col]]

  integrated <- vapply(unique(pair_key), function(pair) {
    rows <- which(pair_key == pair)
    weights <- trapezoid_window_weights(out[[time_col]][rows])
    sum(weights * out$estimate[rows]^2)
  }, numeric(1))

  list(
    contrasts = out,
    contrast_matrix = centered_matrix,
    contrast_covariance = contrast_covariance,
    statistic = if (any(usable)) max(abs(z_value[usable])) else NA_real_,
    integrated_shape_difference = if (length(integrated)) max(integrated) else NA_real_
  )
}


#' @keywords internal
trapezoid_window_weights <- function(time_values) {
  time_values <- as.numeric(time_values)
  n <- length(time_values)
  if (n == 1L) {
    return(1)
  }
  gaps <- diff(time_values)
  if (any(!is.finite(gaps)) || any(gaps <= 0)) {
    cli::cli_abort("Shape-test time values must be finite and strictly increasing.")
  }
  weights <- c(gaps[[1L]] / 2, (gaps[-length(gaps)] + gaps[-1L]) / 2, gaps[[length(gaps)]] / 2)
  weights / sum(weights)
}


#' @keywords internal
fit_shape_null_model <- function(mf) {
  null_trajectory <- mf$trajectory
  null_trajectory$design <- naive_trajectory()
  time_col <- null_trajectory$spec$time$column
  spline_k <- min(
    mf$control$spline$k,
    length(unique(mf$data[[time_col]]))
  )
  formula <- compile_formula(
    null_trajectory,
    engine = "spline",
    spline_k = spline_k,
    spline_basis = mf$control$spline$basis
  )
  fit <- mgcv::gam(
    formula = formula,
    data = mf$data,
    family = spline_family(mf$control$family),
    method = mf$control$spline$method,
    select = TRUE
  )
  list(fit = fit, formula = formula)
}


#' @keywords internal
estimate_shape_null <- function(
    fit,
    levels,
    time_values,
    center_time_values,
    null_control
) {
  observed <- shape_contrast_components(
    fit,
    levels,
    time_values,
    center_time_values
  )
  if (!is.finite(observed$statistic)) {
    cli::cli_abort("The observed trajectory shape statistic is not estimable.")
  }
  method <- null_control$method

  if (method == "gaussian_approximation") {
    calibrated <- simulate_gaussian_shape_null(
      observed = observed,
      simulations = null_control$simulations,
      seed = null_control$seed
    )
    return(list(
      observed = observed,
      p_value = calibrated$p_value,
      null_statistics = calibrated$null_statistics,
      null_effects = calibrated$null_effects,
      null_model_formula = NA,
      successful_simulations = length(calibrated$null_statistics),
      failed_simulations = 0L,
      method = method
    ))
  }

  mf <- fit$model_frame
  null_model <- fit_shape_null_model(mf)
  resamples <- shape_null_resample_plan(
    null_fit = null_model$fit,
    mf = mf,
    null_control = null_control
  )

  run_one <- function(i) {
    data <- mf$data
    data$.y <- resamples[, i]
    boot_fit <- tryCatch(
      mgcv::gam(
        formula = mf$formula,
        data = data,
        family = spline_family(mf$control$family),
        method = mf$control$spline$method,
        select = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) {
      return(c(statistic = NA_real_, effect = NA_real_))
    }
    components <- tryCatch(
      shape_null_statistics(boot_fit, observed),
      error = function(e) NULL
    )
    if (is.null(components)) {
      return(c(statistic = NA_real_, effect = NA_real_))
    }
    c(
      statistic = components$statistic,
      effect = components$integrated_shape_difference
    )
  }

  indices <- seq_len(null_control$simulations)
  if (null_control$cores > 1L && .Platform$OS.type != "windows") {
    simulated <- parallel::mclapply(indices, run_one, mc.cores = null_control$cores)
  } else {
    if (null_control$cores > 1L && .Platform$OS.type == "windows") {
      cli::cli_warn("Parallel null fitting falls back to serial execution on Windows.")
    }
    simulated <- lapply(indices, run_one)
  }
  simulated <- do.call(rbind, simulated)
  valid <- is.finite(simulated[, "statistic"])
  failed <- sum(!valid)
  if (sum(valid) < ceiling(0.9 * null_control$simulations)) {
    cli::cli_abort("More than 10% of null-model refits failed; the null distribution is not reliable.")
  }
  null_statistics <- simulated[valid, "statistic"]
  null_effects <- simulated[valid, "effect"]

  list(
    observed = observed,
    p_value = (1 + sum(null_statistics >= observed$statistic)) /
      (length(null_statistics) + 1),
    null_statistics = null_statistics,
    null_effects = null_effects,
    null_model_formula = null_model$formula,
    successful_simulations = length(null_statistics),
    failed_simulations = failed,
    method = method
  )
}


#' @keywords internal
fit_integrated_level_null_model <- function(mf) {
  null_trajectory <- mf$trajectory
  focal_columns <- null_trajectory$design$by
  keep <- vapply(
    null_trajectory$spec$covariates,
    function(covariate) !covariate$column %in% focal_columns,
    logical(1)
  )
  null_trajectory$spec$covariates <- null_trajectory$spec$covariates[keep]
  null_trajectory$design <- naive_trajectory()
  time_col <- null_trajectory$spec$time$column
  spline_k <- min(mf$control$spline$k, length(unique(mf$data[[time_col]])))
  formula <- compile_formula(
    null_trajectory,
    engine = "spline",
    spline_k = spline_k,
    spline_basis = mf$control$spline$basis
  )
  fit <- mgcv::gam(
    formula = formula,
    data = mf$data,
    family = spline_family(mf$control$family),
    method = mf$control$spline$method,
    select = TRUE
  )
  list(fit = fit, formula = formula)
}


#' @keywords internal
estimate_integrated_level_null <- function(
    fit,
    integrated_matrix,
    integrated_covariance,
    estimates,
    null_control,
    permutation_levels = NULL
) {
  standard_errors <- sqrt(pmax(0, diag(integrated_covariance)))
  usable <- is.finite(estimates) & is.finite(standard_errors) & standard_errors > 0
  if (!any(usable)) {
    cli::cli_abort("The integrated fitted-trajectory contrast is not estimable.")
  }
  observed_statistic <- max(abs(estimates[usable] / standard_errors[usable]))

  if (null_control$method == "gaussian_approximation") {
    calibrated <- calibrate_max_contrast(
      estimates = estimates,
      covariance = integrated_covariance,
      simulations = null_control$simulations,
      seed = null_control$seed,
      return_null = TRUE,
      minimum_simulations = 19L
    )
    return(list(
      observed_statistic = observed_statistic,
      p_value = calibrated$p_value,
      null_statistics = calibrated$null_statistics,
      null_model_formula = NA,
      successful_simulations = length(calibrated$null_statistics),
      failed_simulations = 0L,
      method = null_control$method
    ))
  }

  mf <- fit$model_frame
  if (null_control$method == "subject_label_permutation") {
    return(estimate_subject_label_permutation_null(
      fit = fit,
      integrated_matrix = integrated_matrix,
      observed_statistic = observed_statistic,
      null_control = null_control,
      permutation_levels = permutation_levels
    ))
  }
  null_model <- fit_integrated_level_null_model(mf)
  resamples <- shape_null_resample_plan(null_model$fit, mf, null_control)
  run_one <- function(i) {
    data <- mf$data
    data$.y <- resamples[, i]
    boot_fit <- tryCatch(
      mgcv::gam(
        formula = mf$formula,
        data = data,
        family = spline_family(mf$control$family),
        method = mf$control$spline$method,
        select = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) return(NA_real_)
    boot_estimates <- as.numeric(integrated_matrix %*% stats::coef(boot_fit))
    projected <- integrated_matrix %*% stats::vcov(boot_fit)
    boot_se <- sqrt(pmax(0, rowSums(projected * integrated_matrix)))
    boot_usable <- is.finite(boot_estimates) & is.finite(boot_se) & boot_se > 0
    if (!any(boot_usable)) return(NA_real_)
    max(abs(boot_estimates[boot_usable] / boot_se[boot_usable]))
  }

  indices <- seq_len(null_control$simulations)
  if (null_control$cores > 1L && .Platform$OS.type != "windows") {
    null_statistics <- unlist(parallel::mclapply(
      indices, run_one, mc.cores = null_control$cores
    ))
  } else {
    if (null_control$cores > 1L && .Platform$OS.type == "windows") {
      cli::cli_warn("Parallel null fitting falls back to serial execution on Windows.")
    }
    null_statistics <- unlist(lapply(indices, run_one))
  }
  valid <- is.finite(null_statistics)
  failed <- sum(!valid)
  if (sum(valid) < ceiling(0.9 * null_control$simulations)) {
    cli::cli_abort(
      "More than 10% of integrated-level null refits failed; inference is not reliable."
    )
  }
  null_statistics <- null_statistics[valid]
  list(
    observed_statistic = observed_statistic,
    p_value = (1 + sum(null_statistics >= observed_statistic)) /
      (length(null_statistics) + 1),
    null_statistics = null_statistics,
    null_model_formula = null_model$formula,
    successful_simulations = length(null_statistics),
    failed_simulations = failed,
    method = null_control$method
  )
}


#' @keywords internal
estimate_subject_label_permutation_null <- function(
    fit,
    integrated_matrix,
    observed_statistic,
    null_control,
    permutation_levels
) {
  mf <- fit$model_frame
  group_columns <- mf$trajectory$design$by
  if (length(group_columns) != 1L) {
    cli::cli_abort("Subject-label permutation requires one trajectory modifier.")
  }
  group_col <- group_columns[[1L]]
  subject_col <- mf$trajectory$spec$subject$column
  subjects <- as.character(mf$data[[subject_col]])
  groups <- as.character(mf$data[[group_col]])
  subject_groups <- split(groups, subjects)
  if (any(vapply(subject_groups, function(x) length(unique(x)) != 1L, logical(1)))) {
    cli::cli_abort(
      "Subject-label permutation requires a trajectory modifier that is constant within subject."
    )
  }
  subject_group <- vapply(subject_groups, `[[`, character(1), 1L)
  permutation_levels <- permutation_levels %||% unique(subject_group)
  eligible_subjects <- names(subject_group)[subject_group %in% permutation_levels]
  if (length(eligible_subjects) < 2L) {
    cli::cli_abort("Too few subjects are eligible for label permutation.")
  }

  old_seed <- preserve_random_seed()
  on.exit(restore_random_seed(old_seed), add = TRUE)
  set.seed(null_control$seed)
  permutations <- replicate(
    null_control$simulations,
    sample(subject_group[eligible_subjects], replace = FALSE),
    simplify = "matrix"
  )
  rownames(permutations) <- eligible_subjects

  run_one <- function(i) {
    data <- mf$data
    permuted_subject_group <- subject_group
    permuted_subject_group[eligible_subjects] <- permutations[, i]
    replacement <- permuted_subject_group[subjects]
    if (is.factor(data[[group_col]])) {
      replacement <- factor(replacement, levels = levels(data[[group_col]]))
    }
    data[[group_col]] <- replacement
    boot_fit <- tryCatch(
      mgcv::gam(
        formula = mf$formula,
        data = data,
        family = spline_family(mf$control$family),
        method = mf$control$spline$method,
        select = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) return(NA_real_)
    boot_estimates <- as.numeric(integrated_matrix %*% stats::coef(boot_fit))
    projected <- integrated_matrix %*% stats::vcov(boot_fit)
    boot_se <- sqrt(pmax(0, rowSums(projected * integrated_matrix)))
    usable <- is.finite(boot_estimates) & is.finite(boot_se) & boot_se > 0
    if (!any(usable)) return(NA_real_)
    max(abs(boot_estimates[usable] / boot_se[usable]))
  }

  indices <- seq_len(null_control$simulations)
  if (null_control$cores > 1L && .Platform$OS.type != "windows") {
    null_statistics <- unlist(parallel::mclapply(
      indices, run_one, mc.cores = null_control$cores
    ))
  } else {
    if (null_control$cores > 1L && .Platform$OS.type == "windows") {
      cli::cli_warn("Parallel null fitting falls back to serial execution on Windows.")
    }
    null_statistics <- unlist(lapply(indices, run_one))
  }
  valid <- is.finite(null_statistics)
  failed <- sum(!valid)
  if (sum(valid) < ceiling(0.9 * null_control$simulations)) {
    cli::cli_abort(
      "More than 10% of subject-label permutation refits failed; inference is not reliable."
    )
  }
  null_statistics <- null_statistics[valid]
  list(
    observed_statistic = observed_statistic,
    p_value = (1 + sum(null_statistics >= observed_statistic)) /
      (length(null_statistics) + 1),
    null_statistics = null_statistics,
    null_model_formula = NA,
    successful_simulations = length(null_statistics),
    failed_simulations = failed,
    method = null_control$method
  )
}


#' @keywords internal
simulate_gaussian_shape_null <- function(observed, simulations, seed) {
  covariance <- observed$contrast_covariance
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- max(decomposition$values, 0) * .Machine$double.eps^0.5
  keep <- decomposition$values > tolerance
  if (!any(keep)) {
    cli::cli_abort("The Gaussian shape-null covariance has no estimable directions.")
  }

  old_seed <- preserve_random_seed()
  on.exit(restore_random_seed(old_seed), add = TRUE)
  set.seed(seed)
  normal <- matrix(
    stats::rnorm(sum(keep) * simulations),
    nrow = sum(keep),
    ncol = simulations
  )
  draws <- decomposition$vectors[, keep, drop = FALSE] %*%
    (sqrt(decomposition$values[keep]) * normal)
  standard_errors <- sqrt(pmax(0, diag(covariance)))
  usable <- is.finite(standard_errors) & standard_errors > 0
  null_statistics <- apply(
    abs(sweep(draws[usable, , drop = FALSE], 1L, standard_errors[usable], "/")),
    2L,
    max
  )

  contrasts <- observed$contrasts
  time_col <- attr(contrasts, "time_variable")
  pair_key <- paste(contrasts$group_1, contrasts$group_2, sep = "\r")
  pair_effects <- vapply(unique(pair_key), function(pair) {
    rows <- which(pair_key == pair)
    weights <- trapezoid_window_weights(contrasts[[time_col]][rows])
    as.numeric(crossprod(weights, draws[rows, , drop = FALSE]^2))
  }, numeric(simulations))
  if (is.null(dim(pair_effects))) {
    null_effects <- pair_effects
  } else {
    null_effects <- apply(pair_effects, 1L, max)
  }

  list(
    p_value = (1 + sum(null_statistics >= observed$statistic)) /
      (simulations + 1),
    null_statistics = null_statistics,
    null_effects = null_effects
  )
}


#' @keywords internal
shape_null_statistics <- function(boot_fit, observed) {
  contrast_matrix <- observed$contrast_matrix
  estimates <- as.numeric(contrast_matrix %*% stats::coef(boot_fit))
  projected <- contrast_matrix %*% stats::vcov(boot_fit)
  standard_errors <- sqrt(pmax(0, rowSums(projected * contrast_matrix)))
  usable <- is.finite(estimates) & is.finite(standard_errors) & standard_errors > 0
  statistic <- if (any(usable)) {
    max(abs(estimates[usable] / standard_errors[usable]))
  } else {
    NA_real_
  }

  contrasts <- observed$contrasts
  time_col <- attr(contrasts, "time_variable")
  pair_key <- paste(contrasts$group_1, contrasts$group_2, sep = "\r")
  integrated <- vapply(unique(pair_key), function(pair) {
    rows <- which(pair_key == pair)
    weights <- trapezoid_window_weights(contrasts[[time_col]][rows])
    sum(weights * estimates[rows]^2)
  }, numeric(1))

  list(
    statistic = statistic,
    integrated_shape_difference = if (length(integrated)) {
      max(integrated)
    } else {
      NA_real_
    }
  )
}


#' @keywords internal
shape_null_resample_plan <- function(null_fit, mf, null_control) {
  simulations <- null_control$simulations
  old_seed <- preserve_random_seed()
  on.exit(restore_random_seed(old_seed), add = TRUE)
  set.seed(null_control$seed)

  fitted <- as.numeric(stats::predict(null_fit, type = "response"))
  family <- mf$control$family
  n <- nrow(mf$data)

  if (null_control$method == "cluster_wild_bootstrap") {
    subject_col <- mf$trajectory$spec$subject$column
    subjects <- as.character(mf$data[[subject_col]])
    subject_levels <- unique(subjects)
    weights <- matrix(
      sample(c(-1, 1), length(subject_levels) * simulations, replace = TRUE),
      nrow = length(subject_levels),
      ncol = simulations
    )
    residuals <- stats::residuals(null_fit, type = "response")
    return(vapply(seq_len(simulations), function(i) {
      fitted + residuals * weights[match(subjects, subject_levels), i]
    }, numeric(n)))
  }

  if (family == "gaussian") {
    sigma <- sqrt(summary(null_fit)$scale)
    return(matrix(
      stats::rnorm(n * simulations, mean = rep(fitted, simulations), sd = sigma),
      nrow = n,
      ncol = simulations
    ))
  }
  if (family == "binomial") {
    return(matrix(
      stats::rbinom(n * simulations, size = 1L, prob = rep(fitted, simulations)),
      nrow = n,
      ncol = simulations
    ))
  }
  if (family == "beta") {
    precision <- null_fit$family$getTheta(trans = TRUE)
    mu <- pmin(pmax(rep(fitted, simulations), .Machine$double.eps), 1 - .Machine$double.eps)
    return(matrix(
      stats::rbeta(n * simulations, mu * precision, (1 - mu) * precision),
      nrow = n,
      ncol = simulations
    ))
  }

  cli::cli_abort("Unsupported null-bootstrap family: {.val {family}}.")
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
  invisible(NULL)
}

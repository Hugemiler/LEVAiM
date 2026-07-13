#####
# High-level feature testing workflows
#####

#' Test assay features one model at a time
#'
#' Fits the same LEVAiM trajectory model to each selected feature in an assay and
#' collates feature-level evidence with multiple-testing correction.
#'
#' @param ds A `LongitudinalDataset`.
#' @param assay Name of the assay in `ds`.
#' @param features Character vector of feature names. Defaults to all features
#'   in the assay.
#' @param subject Subject variable, either a column name or `subject_var()`.
#' @param time Time variable, either a column name or `time_var()`.
#' @param covariates List of covariates created with `group_effect()`,
#'   `continuous_effect()`, or `random_effect()`.
#' @param design A trajectory design created with `naive_trajectory()` or
#'   `varying_trajectory()`.
#' @param control A `levaim_trajectory_control`.
#' @param hypothesis Hypothesis tested for every feature. One of
#'   `"overall_trajectory"`, `"trajectory_difference"`, or
#'   `"covariate_effect"`.
#' @param test_covariate Covariate tested by `"covariate_effect"`. For
#'   `"trajectory_difference"`, this may name the focal trajectory modifier.
#' @param test_levels Optional biological trajectory levels included in a
#'   `"trajectory_difference"` test. Omitted levels remain in the fitted model
#'   but do not contribute to the omnibus contrast.
#' @param test_time_values Optional numeric grid for trajectory-difference
#'   inference. Defaults to 100 points over the observed time domain.
#' @param test_simulations Number of Gaussian coefficient simulations used to
#'   calibrate the maximum standardized trajectory-contrast statistic.
#' @param test_seed Seed used for contrast-test simulation.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param keep_fits Logical. Whether to retain fitted objects and result objects
#'   in the returned object.
#' @param error_action One of `"continue"` or `"stop"`. `"continue"` records
#'   failed feature fits in the result table.
#'
#' @return A `levaim_feature_tests` object.
#' @export
test_assay_features <- function(
    ds,
    assay,
    features = NULL,
    subject,
    time,
    covariates = list(),
    design = naive_trajectory(),
    control = trajectory_control(),
    hypothesis = c("overall_trajectory", "trajectory_difference", "covariate_effect"),
    test_covariate = NULL,
    test_levels = NULL,
    test_time_values = NULL,
    test_simulations = 2000L,
    test_seed = 1L,
    p_adjust_method = "BH",
    keep_fits = FALSE,
    error_action = c("continue", "stop")
) {
  error_action <- match.arg(error_action)
  hypothesis <- match.arg(hypothesis)

  if (!inherits(ds, "LongitudinalDataset")) {
    cli::cli_abort("`ds` must be a `LongitudinalDataset` object.")
  }

  if (!assay %in% names(ds$assays)) {
    cli::cli_abort("Assay {.val {assay}} not found in `ds$assays`.")
  }

  assay_table <- ds$assays[[assay]]
  features <- features %||% colnames(assay_table)

  if (!is.character(features) || length(features) == 0L) {
    cli::cli_abort("`features` must be a non-empty character vector.")
  }

  missing_features <- setdiff(features, colnames(assay_table))
  if (length(missing_features) > 0L) {
    cli::cli_abort(
      "Feature(s) not found in assay {.val {assay}}: {paste(missing_features, collapse = ', ')}"
    )
  }

  subject <- normalize_subject_var(subject)
  time <- normalize_time_var(time)

  if (!is.list(covariates)) {
    cli::cli_abort("`covariates` must be a list of covariate objects.")
  }

  if (!inherits(design, "levaim_trajectory_design")) {
    cli::cli_abort("`design` must be created with `naive_trajectory()` or `varying_trajectory()`.")
  }

  if (!inherits(control, "levaim_trajectory_control")) {
    cli::cli_abort("`control` must be created with `trajectory_control()`.")
  }

  validate_feature_test_hypothesis(
    hypothesis,
    test_covariate,
    covariates,
    design,
    control
  )

  rows <- vector("list", length(features))
  fits <- vector("list", length(features))
  results <- vector("list", length(features))
  contrasts <- vector("list", length(features))
  names(fits) <- features
  names(results) <- features
  names(contrasts) <- features

  for (i in seq_along(features)) {
    feature <- features[[i]]

    feature_out <- tryCatch(
      fit_one_assay_feature(
        ds = ds,
        assay = assay,
        feature = feature,
        subject = subject,
        time = time,
        covariates = covariates,
        design = design,
        control = control,
        hypothesis = hypothesis,
        test_covariate = test_covariate,
        test_levels = test_levels,
        test_time_values = test_time_values,
        test_simulations = test_simulations,
        test_seed = test_seed,
        retain_fit = keep_fits
      ),
      error = function(e) {
        if (error_action == "stop") {
          stop(e)
        }

        list(error = conditionMessage(e))
      }
    )

    if (!is.null(feature_out$error)) {
      rows[[i]] <- data.frame(
        assay = assay,
        feature = feature,
        status = "error",
        p_value = NA_real_,
        test_statistic = NA_real_,
        test_term = NA_character_,
        hypothesis = hypothesis,
        tested_covariate = test_covariate %||% NA_character_,
        null_hypothesis = feature_test_null(hypothesis, test_covariate),
        test_method = if (hypothesis == "trajectory_difference") {
          "global_trajectory_max_contrast"
        } else {
          "nested_gam_likelihood_ratio_approximate"
        },
        inference_warning = if (hypothesis == "trajectory_difference") {
          "Joint Wald trajectory-contrast test could not be estimated for this feature."
        } else {
          feature_test_warning()
        },
        full_formula = NA_character_,
        reduced_formula = NA_character_,
        tested_levels = paste(test_levels %||% character(), collapse = ";"),
        max_abs_z = NA_real_,
        max_difference_time = NA_real_,
        integrated_squared_difference = NA_real_,
        n = NA_integer_,
        r_sq = NA_real_,
        deviance_explained = NA_real_,
        error = feature_out$error,
        stringsAsFactors = FALSE
      )

      next
    }

    rows[[i]] <- feature_out$row
    fits[[feature]] <- feature_out$fit
    results[[feature]] <- feature_out$results
    contrasts[[feature]] <- feature_out$contrasts
  }

  table <- do.call(rbind, rows)
  table$p_adjust <- NA_real_

  valid <- is.finite(table$p_value)
  if (any(valid)) {
    table$p_adjust[valid] <- stats::p.adjust(
      table$p_value[valid],
      method = p_adjust_method
    )
  }

  out <- list(
    table = table,
    assay = assay,
    features = features,
    p_adjust_method = p_adjust_method,
    hypothesis = hypothesis,
    test_covariate = test_covariate,
    test_levels = test_levels,
    test_time_values = test_time_values,
    test_simulations = test_simulations,
    test_seed = test_seed,
    control = control
  )

  if (isTRUE(keep_fits)) {
    out$fits <- fits
    out$results <- results
    out$contrasts <- contrasts
  }

  structure(out, class = "levaim_feature_tests")
}


#' @keywords internal
fit_one_assay_feature <- function(
    ds,
    assay,
    feature,
    subject,
    time,
    covariates,
    design,
    control,
    hypothesis,
    test_covariate,
    test_levels,
    test_time_values,
    test_simulations,
    test_seed,
    retain_fit
) {
  spec <- model_spec(
    response = assay_feature(assay, feature),
    subject = subject,
    time = time,
    covariates = covariates
  )

  traj <- trajectory(spec, design = design)
  mf <- model_frame(ds, traj, control = control)
  display_fit <- if (hypothesis == "trajectory_difference" || isTRUE(retain_fit)) {
    fit_trajectory(mf)
  } else {
    NULL
  }
  test <- if (hypothesis == "trajectory_difference") {
    test_trajectory_contrast_hypothesis(
      display_fit,
      levels = test_levels,
      time_values = test_time_values,
      simulations = test_simulations,
      seed = test_seed
    )
  } else {
    test_feature_hypothesis(mf, hypothesis, test_covariate)
  }
  fit <- if (!is.null(display_fit)) {
    display_fit
  } else {
    structure(
      list(fit = test$full_fit, model_frame = mf, engine = "spline"),
      class = "levaim_fit"
    )
  }
  results <- trajectory_results(fit)

  row <- data.frame(
    assay = assay,
    feature = feature,
    status = "ok",
    p_value = test$p_value,
    test_statistic = test$test_statistic,
    test_term = test$test_term,
    hypothesis = hypothesis,
    tested_covariate = test_covariate %||% NA_character_,
    null_hypothesis = feature_test_null(hypothesis, test_covariate),
    test_method = test$test_method,
    inference_warning = test$inference_warning,
    full_formula = paste(deparse(test$full_formula), collapse = " "),
    reduced_formula = paste(deparse(test$reduced_formula), collapse = " "),
    tested_levels = paste(test$tested_levels %||% character(), collapse = ";"),
    max_abs_z = test$max_abs_z %||% NA_real_,
    max_difference_time = test$max_difference_time %||% NA_real_,
    integrated_squared_difference = test$integrated_squared_difference %||% NA_real_,
    n = results$n,
    r_sq = results$r_sq,
    deviance_explained = results$deviance_explained,
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  list(
    row = row,
    fit = fit,
    results = results,
    contrasts = test$contrasts %||% NULL
  )
}


#' @keywords internal
test_trajectory_contrast_hypothesis <- function(
    fit,
    levels = NULL,
    time_values = NULL,
    simulations = 2000L,
    seed = 1L
) {
  contrasts <- trajectory_contrasts(
    fit,
    n = 100L,
    time_values = time_values,
    levels = levels
  )
  contrast_matrix <- attr(contrasts, "contrast_matrix")
  covariance <- attr(contrasts, "coefficient_covariance")
  contrast_estimates <- contrast_matrix %*% stats::coef(fit$fit)
  contrast_covariance <- contrast_matrix %*% covariance %*% t(contrast_matrix)
  calibrated <- calibrate_max_contrast(
    estimates = as.numeric(contrast_estimates),
    covariance = contrast_covariance,
    simulations = simulations,
    seed = seed
  )

  max_row <- if (all(!is.finite(contrasts$z_value))) {
    NA_integer_
  } else {
    which.max(abs(contrasts$z_value))
  }
  time_col <- fit$model_frame$trajectory$spec$time$column

  list(
    p_value = calibrated$p_value,
    test_statistic = calibrated$statistic,
    test_term = "trajectory_difference",
    test_method = "global_trajectory_max_contrast",
    inference_warning = paste(
      "Maximum standardized fitted-contrast test calibrated from the model",
      "coefficient covariance; pointwise p-values only localize separation."
    ),
    full_formula = fit$model_frame$formula,
    reduced_formula = NA,
    full_fit = fit$fit,
    tested_levels = attr(contrasts, "tested_levels"),
    max_abs_z = if (is.na(max_row)) NA_real_ else abs(contrasts$z_value[[max_row]]),
    max_difference_time = if (is.na(max_row)) NA_real_ else contrasts[[time_col]][[max_row]],
    integrated_squared_difference = mean(contrasts$estimate^2, na.rm = TRUE),
    contrasts = contrasts
  )
}


#' @keywords internal
calibrate_max_contrast <- function(estimates, covariance, simulations, seed) {
  simulations <- as.integer(simulations)
  if (length(simulations) != 1L || is.na(simulations) || simulations < 199L) {
    cli::cli_abort("`test_simulations` must be an integer >= 199.")
  }

  standard_errors <- sqrt(pmax(0, diag(covariance)))
  usable <- is.finite(estimates) & is.finite(standard_errors) & standard_errors > 0
  if (!any(usable)) {
    return(list(statistic = NA_real_, p_value = NA_real_))
  }

  correlation <- stats::cov2cor(covariance[usable, usable, drop = FALSE])
  decomposition <- eigen(correlation, symmetric = TRUE)
  tolerance <- max(decomposition$values, 0) * .Machine$double.eps^0.5
  keep <- decomposition$values > tolerance
  observed <- max(abs(estimates[usable] / standard_errors[usable]))

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  normal <- matrix(
    stats::rnorm(sum(keep) * simulations),
    nrow = sum(keep),
    ncol = simulations
  )
  simulated <- decomposition$vectors[, keep, drop = FALSE] %*%
    (sqrt(decomposition$values[keep]) * normal)
  simulated_max <- apply(abs(simulated), 2L, max)

  list(
    statistic = observed,
    p_value = (1 + sum(simulated_max >= observed)) / (simulations + 1)
  )
}


#' @keywords internal
validate_feature_test_hypothesis <- function(
    hypothesis,
    test_covariate,
    covariates,
    design,
    control
) {
  if (control$engine != "spline") {
    cli::cli_abort(
      "Explicit assay-wide hypothesis tests are currently implemented for `engine = 'spline'` only."
    )
  }

  covariate_columns <- vapply(covariates, `[[`, character(1), "column")

  if (hypothesis == "trajectory_difference") {
    if (design$type != "varying") {
      cli::cli_abort(
        "`hypothesis = 'trajectory_difference'` requires `varying_trajectory()`."
      )
    }
    if (!is.null(test_covariate) && !test_covariate %in% design$by) {
      cli::cli_abort("`test_covariate` must be one of the trajectory modifiers in `design$by`.")
    }
  }

  if (hypothesis == "covariate_effect") {
    if (is.null(test_covariate) || length(test_covariate) != 1L) {
      cli::cli_abort(
        "`hypothesis = 'covariate_effect'` requires one `test_covariate`."
      )
    }
    if (!test_covariate %in% covariate_columns) {
      cli::cli_abort("`test_covariate` must name a declared covariate.")
    }
  }

  invisible(TRUE)
}


#' @keywords internal
test_feature_hypothesis <- function(mf, hypothesis, test_covariate = NULL) {
  formulas <- feature_test_formulas(mf, hypothesis, test_covariate)
  test_data <- mf$data
  if (!is.null(formulas$test_group)) {
    test_data$.levaim_test_group <- formulas$test_group
  }

  reduced <- mgcv::gam(
    formula = formulas$reduced,
    data = test_data,
    family = spline_family(mf$control$family),
    method = "ML",
    select = TRUE
  )
  full <- mgcv::gam(
    formula = formulas$full,
    data = test_data,
    family = spline_family(mf$control$family),
    method = "ML",
    select = TRUE
  )

  reduced_log_lik <- stats::logLik(reduced)
  full_log_lik <- stats::logLik(full)
  statistic <- max(
    0,
    2 * (as.numeric(full_log_lik) - as.numeric(reduced_log_lik))
  )
  degrees_freedom <- full$rank - reduced$rank
  p_value <- if (is.finite(degrees_freedom) && degrees_freedom > 0) {
    stats::pchisq(statistic, df = degrees_freedom, lower.tail = FALSE)
  } else {
    NA_real_
  }

  list(
    p_value = p_value,
    test_statistic = statistic,
    test_term = hypothesis,
    test_method = "nested_gam_likelihood_ratio_approximate",
    inference_warning = feature_test_warning(),
    full_formula = formulas$full,
    reduced_formula = formulas$reduced,
    full_fit = full,
    tested_levels = NULL,
    contrasts = NULL
  )
}


#' @keywords internal
feature_test_warning <- function() {
  paste(
    "Approximate likelihood-ratio test using nominal model-rank degrees of freedom",
    "from penalized GAM fits; confirm priority findings with sensitivity analyses."
  )
}


#' @keywords internal
feature_test_formulas <- function(mf, hypothesis, test_covariate = NULL) {
  traj <- mf$trajectory
  control <- mf$control

  compile_one <- function(x) {
    compile_formula(
      x,
      engine = "spline",
      spline_k = control$spline$k,
      spline_basis = control$spline$basis
    )
  }

  if (hypothesis == "trajectory_difference") {
    reduced_traj <- traj
    reduced_traj$design <- naive_trajectory()
    reduced_formula <- compile_one(reduced_traj)
    test_group <- interaction(
      mf$data[traj$design$by],
      drop = TRUE,
      sep = ":"
    )
    test_group <- ordered(test_group, levels = levels(test_group))
    difference_term <- compile_difference_smooth(
      time_col = traj$spec$time$column,
      k = control$spline$k,
      basis = control$spline$basis
    )
    full_formula <- stats::update.formula(
      reduced_formula,
      paste(". ~ . +", difference_term)
    )
    return(list(
      full = full_formula,
      reduced = reduced_formula,
      test_group = test_group
    ))
  }

  if (hypothesis == "covariate_effect") {
    full_traj <- traj
    full_traj$design <- naive_trajectory()
    reduced_traj <- full_traj
    keep <- vapply(
      reduced_traj$spec$covariates,
      function(x) !identical(x$column, test_covariate),
      logical(1)
    )
    reduced_traj$spec$covariates <- reduced_traj$spec$covariates[keep]
    return(list(
      full = compile_one(full_traj),
      reduced = compile_one(reduced_traj),
      test_group = NULL
    ))
  }

  list(
    full = mf$formula,
    reduced = compile_spline_null_formula(traj),
    test_group = NULL
  )
}


#' @keywords internal
compile_difference_smooth <- function(time_col, k, basis) {
  if (basis == "gam") {
    return(paste0(
      "s(", time_col, ", by = .levaim_test_group, k = ", as.integer(k), ")"
    ))
  }

  spline <- if (basis == "natural") {
    paste0("splines::ns(", time_col, ", df = ", as.integer(k), ")")
  } else {
    paste0(
      "splines::bs(", time_col, ", df = ", as.integer(k), ", degree = 3)"
    )
  }

  paste0(".levaim_test_group:", spline)
}


#' @keywords internal
compile_spline_null_formula <- function(traj) {
  spec <- traj$spec
  rhs_terms <- c(
    compile_covariate_terms_spline(spec$covariates),
    paste0("s(", spec$subject$column, ", bs = 're')")
  )
  rhs_terms <- unique(rhs_terms[nzchar(rhs_terms)])
  stats::as.formula(paste(".y ~", paste(rhs_terms, collapse = " + ")))
}


#' @keywords internal
feature_test_null <- function(hypothesis, test_covariate = NULL) {
  switch(
    hypothesis,
    overall_trajectory = "The response has no longitudinal trajectory after adjustment for declared covariates.",
    trajectory_difference = "The response follows one shared longitudinal trajectory rather than group-specific trajectories.",
    covariate_effect = paste0(
      "The adjusted response does not differ by `",
      test_covariate,
      "` under a shared longitudinal trajectory."
    )
  )
}


#' @keywords internal
extract_feature_test <- function(results, time_col) {
  smooth <- results$smooth_terms

  if (!is.data.frame(smooth) || nrow(smooth) == 0L) {
    return(empty_feature_test())
  }

  terms <- rownames(smooth)
  p_col <- find_p_value_column(smooth)

  if (is.na(p_col)) {
    return(empty_feature_test())
  }

  trajectory_rows <- grepl(time_col, terms, fixed = TRUE)
  if (!any(trajectory_rows)) {
    trajectory_rows <- rep(TRUE, nrow(smooth))
  }

  p_values <- as.numeric(smooth[[p_col]])
  p_values[!trajectory_rows] <- NA_real_

  if (all(is.na(p_values))) {
    return(empty_feature_test())
  }

  best <- which.min(p_values)
  stat_col <- find_test_statistic_column(smooth)
  statistic <- if (is.na(stat_col)) NA_real_ else as.numeric(smooth[[stat_col]][[best]])

  list(
    p_value = p_values[[best]],
    test_statistic = statistic,
    test_term = terms[[best]]
  )
}


#' @keywords internal
empty_feature_test <- function() {
  list(
    p_value = NA_real_,
    test_statistic = NA_real_,
    test_term = NA_character_
  )
}


#' @keywords internal
find_p_value_column <- function(x) {
  candidates <- c("p-value", "p.value", "Pr(>Chi)", "Pr(>Chisq)", "Pr(>F)", "Pr(>|t|)", "Pr(>|z|)")
  matched <- intersect(candidates, colnames(x))

  if (length(matched) == 0L) {
    return(NA_character_)
  }

  matched[[1]]
}


#' @keywords internal
find_test_statistic_column <- function(x) {
  candidates <- c("F", "Chi.sq", "t value", "z value")
  matched <- intersect(candidates, colnames(x))

  if (length(matched) == 0L) {
    return(NA_character_)
  }

  matched[[1]]
}


#' @keywords internal
normalize_subject_var <- function(x) {
  if (inherits(x, "levaim_subject_var")) {
    return(x)
  }

  if (is.character(x) && length(x) == 1L) {
    return(subject_var(x))
  }

  cli::cli_abort("`subject` must be a column name or an object created with `subject_var()`.")
}


#' @keywords internal
normalize_time_var <- function(x) {
  if (inherits(x, "levaim_time_var")) {
    return(x)
  }

  if (is.character(x) && length(x) == 1L) {
    return(time_var(x))
  }

  cli::cli_abort("`time` must be a column name or an object created with `time_var()`.")
}


#' Print LEVAiM feature tests
#'
#' @param x A `levaim_feature_tests` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.levaim_feature_tests <- function(x, ...) {
  cli::cli_h1("LEVAiM feature tests")
  cli::cli_text("Assay: {.field {x$assay}}")
  cli::cli_text("Features: {length(x$features)}")
  cli::cli_text("P-value adjustment: {.field {x$p_adjust_method}}")

  ok <- sum(x$table$status == "ok", na.rm = TRUE)
  cli::cli_text("Successful fits: {ok}")

  print(x$table)

  invisible(x)
}


#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

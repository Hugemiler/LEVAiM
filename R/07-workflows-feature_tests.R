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
#'   `"overall_trajectory"`, `"trajectory_difference"`,
#'   `"trajectory_shape_difference"`, `"critical_period_level"`,
#'   `"occurrence_difference"`, `"positive_abundance_difference"`,
#'   `"persistent_level_difference"`, `"level_or_shape"`, or
#'   `"covariate_effect"`.
#' @param test_covariate Covariate tested by `"covariate_effect"`. For
#'   `"trajectory_difference"`, this may name the focal trajectory modifier.
#' @param test_levels Optional biological trajectory levels included in a
#'   group-trajectory test. Omitted levels remain in the fitted model but do
#'   not contribute to the omnibus contrast.
#' @param test_time_values Optional numeric grid for trajectory inference.
#'   Defaults to 100 points over the observed domain, or a grid inside
#'   `test_window` for critical-period and component hypotheses.
#' @param test_center_time_values Optional numeric grid defining the domain over
#'   which shape contrasts are centered. `NULL` uses the full observed time
#'   domain, even when `test_time_values` localizes inference to a narrower
#'   biological window.
#' @param test_simulations Number of Gaussian coefficient simulations used to
#'   calibrate the maximum standardized trajectory-contrast statistic.
#' @param test_seed Seed used for contrast-test simulation.
#' @param null_control Null-distribution controls created with
#'   [trajectory_null_control()]. Used by centered-shape and integrated-level
#'   trajectory hypotheses.
#' @param test_window Optional two-number closed time window. It is required by
#'   `"critical_period_level"` and may localize occurrence or
#'   positive-abundance tests.
#' @param detection_threshold Numeric abundance above which a feature is
#'   considered detected for occurrence and positive-abundance analyses.
#' @param composition_policy Whether assay values are used `"as_is"` or each
#'   sample is `"renormalize"`d to total 100 before analysis.
#' @param p_adjust_method Method passed to [stats::p.adjust()].
#' @param fdr_level FDR threshold used to diagnose whether the empirical
#'   p-value resolution is adequate for assay-wide testing.
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
    hypothesis = c(
      "overall_trajectory", "trajectory_difference",
      "trajectory_shape_difference", "critical_period_level",
      "occurrence_difference", "positive_abundance_difference",
      "persistent_level_difference", "level_or_shape", "covariate_effect"
    ),
    test_covariate = NULL,
    test_levels = NULL,
    test_time_values = NULL,
    test_center_time_values = NULL,
    test_simulations = 2000L,
    test_seed = 1L,
    null_control = NULL,
    test_window = NULL,
    detection_threshold = 0,
    composition_policy = c("as_is", "renormalize"),
    p_adjust_method = "BH",
    fdr_level = 0.05,
    keep_fits = FALSE,
    error_action = c("continue", "stop")
) {
  error_action <- match.arg(error_action)
  hypothesis <- match.arg(hypothesis)
  composition_policy <- match.arg(composition_policy)

  if (!inherits(ds, "LongitudinalDataset")) {
    cli::cli_abort("`ds` must be a `LongitudinalDataset` object.")
  }

  if (!assay %in% names(ds$assays)) {
    cli::cli_abort("Assay {.val {assay}} not found in `ds$assays`.")
  }

  validate_feature_test_window(test_window, hypothesis)
  if (!is.null(test_window) && !is.null(test_time_values) &&
      any(test_time_values < test_window[[1]] | test_time_values > test_window[[2]])) {
    cli::cli_abort("Every `test_time_values` entry must lie inside `test_window`.")
  }
  if (!is.numeric(detection_threshold) || length(detection_threshold) != 1L ||
      !is.finite(detection_threshold) || detection_threshold < 0) {
    cli::cli_abort("`detection_threshold` must be one finite non-negative number.")
  }

  assay_table <- apply_composition_policy(ds$assays[[assay]], composition_policy)
  ds$assays[[assay]] <- assay_table
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
  if (!is.null(test_window)) {
    observed_time <- range(ds$metadata[[time$column]], na.rm = TRUE)
    if (test_window[[1]] < observed_time[[1]] || test_window[[2]] > observed_time[[2]]) {
      cli::cli_abort("`test_window` must lie inside the observed time domain.")
    }
  }

  if (!is.list(covariates)) {
    cli::cli_abort("`covariates` must be a list of covariate objects.")
  }

  if (!inherits(design, "levaim_trajectory_design")) {
    cli::cli_abort("`design` must be created with `naive_trajectory()` or `varying_trajectory()`.")
  }

  if (!inherits(control, "levaim_trajectory_control")) {
    cli::cli_abort("`control` must be created with `trajectory_control()`.")
  }

  if (!is.numeric(fdr_level) || length(fdr_level) != 1L ||
      !is.finite(fdr_level) || fdr_level <= 0 || fdr_level >= 1) {
    cli::cli_abort("`fdr_level` must be one number strictly between 0 and 1.")
  }

  resolved_null_control <- NULL
  robust_null_hypotheses <- c(
    "trajectory_shape_difference", "critical_period_level",
    "occurrence_difference", "positive_abundance_difference",
    "persistent_level_difference"
  )
  if (hypothesis %in% robust_null_hypotheses) {
    null_family <- if (hypothesis == "occurrence_difference") {
      "binomial"
    } else {
      control$family
    }
    resolved_null_control <- resolve_trajectory_null_control(
      null_control = null_control,
      family = null_family,
      seed = test_seed
    )
    if (hypothesis == "trajectory_shape_difference" &&
        resolved_null_control$method == "subject_label_permutation") {
      cli::cli_abort(
        "`subject_label_permutation` is not valid for the centered-shape null, which permits group offsets."
      )
    }
  } else if (!is.null(null_control)) {
    cli::cli_abort(
      "`null_control` is used only with a shape or integrated-level trajectory hypothesis."
    )
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
  null_distributions <- vector("list", length(features))
  names(fits) <- features
  names(results) <- features
  names(contrasts) <- features
  names(null_distributions) <- features

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
        test_center_time_values = test_center_time_values,
        test_simulations = test_simulations,
        test_seed = test_seed,
        null_control = resolved_null_control,
        test_window = test_window,
        detection_threshold = detection_threshold,
        composition_policy = composition_policy,
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
        test_method = feature_test_method(hypothesis, resolved_null_control),
        inference_warning = feature_test_inference_warning(
          hypothesis,
          resolved_null_control,
          failed = TRUE
        ),
        full_formula = NA_character_,
        reduced_formula = NA_character_,
        tested_levels = paste(test_levels %||% character(), collapse = ";"),
        max_abs_z = NA_real_,
        max_difference_time = NA_real_,
        integrated_squared_difference = NA_real_,
        integrated_level_difference = NA_real_,
        integrated_squared_shape_difference = NA_real_,
        null_centered_shape_difference = NA_real_,
        null_mean_shape_difference = NA_real_,
        null_mean_statistic = NA_real_,
        null_method = resolved_null_control$method %||% NA_character_,
        null_simulations = NA_integer_,
        null_failed_simulations = NA_integer_,
        p_value_resolution = NA_real_,
        response_component = feature_test_response_component(hypothesis),
        window_start = if (is.null(test_window)) NA_real_ else test_window[[1]],
        window_end = if (is.null(test_window)) NA_real_ else test_window[[2]],
        effect_level_1 = NA_character_,
        effect_level_2 = NA_character_,
        model_effect = NA_real_,
        model_std_error = NA_real_,
        model_conf_low = NA_real_,
        model_conf_high = NA_real_,
        effect_scale = NA_character_,
        descriptive_difference = NA_real_,
        descriptive_summary = NA_character_,
        n_detected = NA_integer_,
        detection_rate = NA_real_,
        min_group_n = NA_integer_,
        min_group_detection_rate = NA_real_,
        composition_policy = composition_policy,
        assay_total_min = min(rowSums(assay_table, na.rm = TRUE)),
        assay_total_max = max(rowSums(assay_table, na.rm = TRUE)),
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
    null_distributions[[feature]] <- feature_out$null_distribution
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

  fdr_resolution <- feature_test_fdr_resolution(table, fdr_level)
  if (hypothesis %in% robust_null_hypotheses &&
      is.finite(fdr_resolution$single_hit_min_q) &&
      fdr_resolution$single_hit_min_q > fdr_level) {
    cli::cli_warn(c(
      "The simulated null is too coarse to support an isolated FDR discovery at q <= {fdr_level}.",
      "i" = "With {fdr_resolution$features_tested} tested features, the best attainable single-hit adjusted p-value is {format(fdr_resolution$single_hit_min_q, digits = 3)}.",
      "i" = "Increase `null_control$simulations`; this diagnostic does not invalidate larger groups of concordant small p-values."
    ))
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
    test_center_time_values = test_center_time_values,
    test_simulations = test_simulations,
    test_seed = test_seed,
    null_control = resolved_null_control,
    test_window = test_window,
    detection_threshold = detection_threshold,
    composition_policy = composition_policy,
    fdr_level = fdr_level,
    fdr_resolution = fdr_resolution,
    control = control
  )

  if (isTRUE(keep_fits)) {
    out$fits <- fits
    out$results <- results
    out$contrasts <- contrasts
  }
  if (!is.null(resolved_null_control) && isTRUE(resolved_null_control$store_null)) {
    out$null_distributions <- null_distributions
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
    test_center_time_values,
    test_simulations,
    test_seed,
    null_control,
    test_window,
    detection_threshold,
    composition_policy,
    retain_fit
) {
  prepared <- prepare_feature_response(
    ds, assay, feature, hypothesis, detection_threshold, control
  )
  ds <- prepared$ds
  control <- prepared$control

  spec <- model_spec(
    response = assay_feature(assay, feature),
    subject = subject,
    time = time,
    covariates = covariates
  )

  traj <- trajectory(spec, design = design)
  mf <- model_frame(ds, traj, control = control)
  is_contrast_test <- hypothesis %in% c(
    "trajectory_difference", "trajectory_shape_difference",
    "persistent_level_difference", "level_or_shape",
    "critical_period_level", "occurrence_difference",
    "positive_abundance_difference"
  )
  display_fit <- if (is_contrast_test || isTRUE(retain_fit)) {
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
  } else if (hypothesis == "trajectory_shape_difference") {
    test_trajectory_shape_hypothesis(
      display_fit,
      levels = test_levels,
      time_values = test_time_values,
      center_time_values = test_center_time_values,
      null_control = null_control
    )
  } else if (hypothesis == "persistent_level_difference") {
    test_persistent_level_hypothesis(
      display_fit,
      levels = test_levels,
      time_values = test_time_values,
      simulations = test_simulations,
      seed = test_seed,
      hypothesis = hypothesis,
      null_control = null_control
    )
  } else if (hypothesis == "level_or_shape") {
    test_trajectory_contrast_hypothesis(
      display_fit,
      levels = test_levels,
      time_values = test_time_values,
      simulations = test_simulations,
      seed = test_seed,
      hypothesis = "level_or_shape"
    )
  } else if (hypothesis %in% c(
    "critical_period_level", "occurrence_difference",
    "positive_abundance_difference"
  )) {
    test_persistent_level_hypothesis(
      display_fit,
      levels = test_levels,
      time_values = resolve_feature_test_time_values(
        display_fit, test_time_values, test_window
      ),
      simulations = test_simulations,
      seed = test_seed,
      hypothesis = hypothesis,
      null_control = null_control
    )
  } else {
    test_feature_hypothesis(mf, hypothesis, test_covariate)
  }
  fit_model_frame <- test$model_frame %||% mf
  fit <- if (!is.null(display_fit)) {
    display_fit
  } else {
    structure(
      list(fit = test$full_fit, model_frame = fit_model_frame, engine = "spline"),
      class = "levaim_fit"
    )
  }
  results <- trajectory_results(fit)
  audit <- feature_descriptive_audit(
    prepared$audit_ds, assay, feature, test_covariate, test_levels, test_window,
    detection_threshold
  )

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
    integrated_level_difference = test$integrated_level_difference %||% NA_real_,
    integrated_squared_shape_difference = test$integrated_squared_shape_difference %||% NA_real_,
    null_centered_shape_difference = test$null_centered_shape_difference %||% NA_real_,
    null_mean_shape_difference = test$null_mean_shape_difference %||% NA_real_,
    null_mean_statistic = test$null_mean_statistic %||% NA_real_,
    null_method = test$null_method %||% NA_character_,
    null_simulations = test$null_simulations %||% NA_integer_,
    null_failed_simulations = test$null_failed_simulations %||% NA_integer_,
    p_value_resolution = test$p_value_resolution %||% NA_real_,
    response_component = prepared$response_component,
    window_start = if (is.null(test_window)) NA_real_ else test_window[[1]],
    window_end = if (is.null(test_window)) NA_real_ else test_window[[2]],
    effect_level_1 = test$effect_level_1 %||% NA_character_,
    effect_level_2 = test$effect_level_2 %||% NA_character_,
    model_effect = test$model_effect %||% NA_real_,
    model_std_error = test$model_std_error %||% NA_real_,
    model_conf_low = test$model_conf_low %||% NA_real_,
    model_conf_high = test$model_conf_high %||% NA_real_,
    effect_scale = test$effect_scale %||% NA_character_,
    descriptive_difference = audit$descriptive_difference,
    descriptive_summary = audit$descriptive_summary,
    n_detected = audit$n_detected,
    detection_rate = audit$detection_rate,
    min_group_n = audit$min_group_n,
    min_group_detection_rate = audit$min_group_detection_rate,
    composition_policy = composition_policy,
    assay_total_min = min(rowSums(prepared$audit_ds$assays[[assay]], na.rm = TRUE)),
    assay_total_max = max(rowSums(prepared$audit_ds$assays[[assay]], na.rm = TRUE)),
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
    contrasts = test$contrasts %||% NULL,
    null_distribution = test$null_distribution %||% NULL
  )
}


#' @keywords internal
test_trajectory_shape_hypothesis <- function(
    fit,
    levels = NULL,
    time_values = NULL,
    center_time_values = NULL,
    null_control
) {
  estimated <- estimate_shape_null(
    fit,
    levels,
    time_values,
    center_time_values,
    null_control
  )
  observed <- estimated$observed
  contrasts <- observed$contrasts
  max_row <- if (all(!is.finite(contrasts$z_value))) {
    NA_integer_
  } else {
    which.max(abs(contrasts$z_value))
  }
  time_col <- fit$model_frame$trajectory$spec$time$column
  null_mean_effect <- if (any(is.finite(estimated$null_effects))) {
    mean(estimated$null_effects, na.rm = TRUE)
  } else {
    NA_real_
  }

  list(
    p_value = estimated$p_value,
    test_statistic = observed$statistic,
    test_term = "trajectory_shape_difference",
    test_method = feature_test_method("trajectory_shape_difference", null_control),
    inference_warning = feature_test_inference_warning(
      "trajectory_shape_difference",
      null_control
    ),
    full_formula = fit$model_frame$formula,
    reduced_formula = estimated$null_model_formula,
    full_fit = fit$fit,
    tested_levels = attr(contrasts, "tested_levels"),
    max_abs_z = if (is.na(max_row)) NA_real_ else abs(contrasts$z_value[[max_row]]),
    max_difference_time = if (is.na(max_row)) NA_real_ else contrasts[[time_col]][[max_row]],
    integrated_squared_difference = observed$integrated_shape_difference,
    integrated_squared_shape_difference = observed$integrated_shape_difference,
    null_centered_shape_difference = observed$integrated_shape_difference - null_mean_effect,
    null_mean_shape_difference = null_mean_effect,
    null_mean_statistic = mean(estimated$null_statistics, na.rm = TRUE),
    null_method = estimated$method,
    null_simulations = estimated$successful_simulations,
    null_failed_simulations = estimated$failed_simulations,
    p_value_resolution = 1 / (estimated$successful_simulations + 1),
    contrasts = contrasts,
    null_distribution = if (isTRUE(null_control$store_null)) {
      data.frame(
        statistic = estimated$null_statistics,
        integrated_squared_shape_difference = estimated$null_effects
      )
    } else {
      NULL
    }
  )
}


#' @keywords internal
test_trajectory_contrast_hypothesis <- function(
    fit,
    levels = NULL,
    time_values = NULL,
    simulations = 2000L,
    seed = 1L,
    hypothesis = "trajectory_difference"
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
    test_term = hypothesis,
    test_method = feature_test_method(hypothesis),
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
    effect_level_1 = if (is.na(max_row)) NA_character_ else contrasts$group_1[[max_row]],
    effect_level_2 = if (is.na(max_row)) NA_character_ else contrasts$group_2[[max_row]],
    model_effect = if (is.na(max_row)) NA_real_ else contrasts$estimate[[max_row]],
    model_std_error = if (is.na(max_row)) NA_real_ else contrasts$std_error[[max_row]],
    model_conf_low = if (is.na(max_row)) NA_real_ else contrasts$lower[[max_row]],
    model_conf_high = if (is.na(max_row)) NA_real_ else contrasts$upper[[max_row]],
    effect_scale = if (fit$model_frame$control$family == "binomial") {
      "log_odds"
    } else {
      paste0(fit$model_frame$control$transform, "_abundance")
    },
    contrasts = contrasts
  )
}


#' @keywords internal
test_persistent_level_hypothesis <- function(
    fit,
    levels = NULL,
    time_values = NULL,
    simulations = 2000L,
    seed = 1L,
    hypothesis = "persistent_level_difference",
    null_control = NULL
) {
  contrasts <- trajectory_contrasts(
    fit, n = 100L, time_values = time_values, levels = levels
  )
  contrast_matrix <- attr(contrasts, "contrast_matrix")
  covariance <- attr(contrasts, "coefficient_covariance")
  pair_id <- interaction(
    contrasts$group_1, contrasts$group_2, drop = TRUE, sep = " - "
  )
  pair_rows <- split(seq_len(nrow(contrasts)), pair_id)
  integrated_matrix <- do.call(rbind, lapply(
    pair_rows,
    function(rows) colMeans(contrast_matrix[rows, , drop = FALSE])
  ))
  estimates <- as.numeric(integrated_matrix %*% stats::coef(fit$fit))
  integrated_covariance <- integrated_matrix %*% covariance %*% t(integrated_matrix)
  estimated_null <- if (!is.null(null_control)) {
    estimate_integrated_level_null(
      fit, integrated_matrix, integrated_covariance, estimates, null_control,
      permutation_levels = attr(contrasts, "tested_levels")
    )
  } else {
    NULL
  }
  calibrated <- if (is.null(estimated_null)) {
    calibrate_max_contrast(estimates, integrated_covariance, simulations, seed)
  } else {
    list(statistic = estimated_null$observed_statistic, p_value = estimated_null$p_value)
  }
  standard_errors <- sqrt(pmax(0, diag(integrated_covariance)))
  z_values <- estimates / standard_errors
  strongest <- if (all(!is.finite(z_values))) NA_integer_ else which.max(abs(z_values))
  pair_names <- strsplit(names(pair_rows), " - ", fixed = TRUE)
  z_multiplier <- stats::qnorm(0.975)

  list(
    p_value = calibrated$p_value,
    test_statistic = calibrated$statistic,
    test_term = hypothesis,
    test_method = feature_test_method(hypothesis, null_control),
    inference_warning = feature_test_inference_warning(hypothesis, null_control),
    full_formula = fit$model_frame$formula,
    reduced_formula = if (is.null(estimated_null)) NA else
      estimated_null$null_model_formula,
    full_fit = fit$fit,
    tested_levels = attr(contrasts, "tested_levels"),
    max_abs_z = if (is.na(strongest)) NA_real_ else abs(z_values[[strongest]]),
    integrated_squared_difference = mean(estimates^2, na.rm = TRUE),
    integrated_level_difference = if (is.na(strongest)) NA_real_ else estimates[[strongest]],
    effect_level_1 = if (is.na(strongest)) NA_character_ else pair_names[[strongest]][[1]],
    effect_level_2 = if (is.na(strongest)) NA_character_ else pair_names[[strongest]][[2]],
    model_effect = if (is.na(strongest)) NA_real_ else estimates[[strongest]],
    model_std_error = if (is.na(strongest)) NA_real_ else standard_errors[[strongest]],
    model_conf_low = if (is.na(strongest)) NA_real_ else
      estimates[[strongest]] - z_multiplier * standard_errors[[strongest]],
    model_conf_high = if (is.na(strongest)) NA_real_ else
      estimates[[strongest]] + z_multiplier * standard_errors[[strongest]],
    effect_scale = if (fit$model_frame$control$family == "binomial") {
      "log_odds"
    } else {
      paste0(fit$model_frame$control$transform, "_abundance")
    },
    null_mean_statistic = if (is.null(estimated_null)) NA_real_ else
      mean(estimated_null$null_statistics, na.rm = TRUE),
    null_method = if (is.null(estimated_null)) "gaussian_approximation" else
      estimated_null$method,
    null_simulations = if (is.null(estimated_null)) simulations else
      estimated_null$successful_simulations,
    null_failed_simulations = if (is.null(estimated_null)) 0L else
      estimated_null$failed_simulations,
    p_value_resolution = if (is.null(estimated_null)) 1 / (simulations + 1) else
      1 / (estimated_null$successful_simulations + 1),
    null_distribution = if (!is.null(estimated_null) && isTRUE(null_control$store_null)) {
      data.frame(statistic = estimated_null$null_statistics)
    } else {
      NULL
    },
    contrasts = contrasts
  )
}


#' @keywords internal
calibrate_max_contrast <- function(
    estimates,
    covariance,
    simulations,
    seed,
    return_null = FALSE,
    minimum_simulations = 199L
) {
  simulations <- as.integer(simulations)
  if (length(simulations) != 1L || is.na(simulations) ||
      simulations < minimum_simulations) {
    cli::cli_abort(
      "`test_simulations` must be an integer >= {minimum_simulations}."
    )
  }

  standard_errors <- sqrt(pmax(0, diag(covariance)))
  usable <- is.finite(estimates) & is.finite(standard_errors) & standard_errors > 0
  if (!any(usable)) {
    return(list(
      statistic = NA_real_,
      p_value = NA_real_,
      null_statistics = if (isTRUE(return_null)) numeric() else NULL
    ))
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

  out <- list(
    statistic = observed,
    p_value = (1 + sum(simulated_max >= observed)) / (simulations + 1)
  )
  if (isTRUE(return_null)) {
    out$null_statistics <- simulated_max
  }
  out
}


#' @keywords internal
apply_composition_policy <- function(assay, policy) {
  assay <- as.data.frame(assay, check.names = FALSE)
  values <- as.matrix(assay)
  if (!is.numeric(values)) {
    cli::cli_abort("Assay composition policies require numeric values.")
  }
  if (policy == "renormalize") {
    if (any(values < 0, na.rm = TRUE)) {
      cli::cli_abort("`composition_policy = 'renormalize'` requires non-negative values.")
    }
    totals <- rowSums(values, na.rm = TRUE)
    usable <- is.finite(totals) & totals > 0
    values[usable, ] <- values[usable, , drop = FALSE] / totals[usable] * 100
  }
  as.data.frame(values, check.names = FALSE)
}


#' @keywords internal
validate_feature_test_window <- function(window, hypothesis) {
  if (is.null(window)) {
    if (hypothesis == "critical_period_level") {
      cli::cli_abort(
        "`hypothesis = 'critical_period_level'` requires `test_window = c(start, end)`."
      )
    }
    return(invisible(TRUE))
  }
  if (!is.numeric(window) || length(window) != 2L || any(!is.finite(window)) ||
      window[[1]] >= window[[2]]) {
    cli::cli_abort("`test_window` must contain two finite increasing time values.")
  }
  if (!hypothesis %in% c(
    "critical_period_level", "occurrence_difference",
    "positive_abundance_difference"
  )) {
    cli::cli_abort(
      "`test_window` is used only by critical-period or component trajectory hypotheses."
    )
  }
  invisible(TRUE)
}


#' @keywords internal
feature_test_response_component <- function(hypothesis) {
  switch(
    hypothesis,
    occurrence_difference = "occurrence",
    positive_abundance_difference = "positive_abundance",
    "abundance"
  )
}


#' @keywords internal
prepare_feature_response <- function(
    ds, assay, feature, hypothesis, detection_threshold, control
) {
  audit_ds <- ds
  response_component <- feature_test_response_component(hypothesis)
  values <- ds$assays[[assay]][[feature]]

  if (hypothesis == "occurrence_difference") {
    ds$assays[[assay]][[feature]] <- ifelse(
      is.na(values), NA_real_, as.numeric(values > detection_threshold)
    )
    control$family <- "binomial"
    control$transform <- "identity"
  } else if (hypothesis == "positive_abundance_difference") {
    if (control$family != "gaussian") {
      cli::cli_abort(
        "`positive_abundance_difference` currently requires `family = 'gaussian'`."
      )
    }
    values[values <= detection_threshold] <- NA_real_
    ds$assays[[assay]][[feature]] <- values
  }

  list(
    ds = ds,
    audit_ds = audit_ds,
    control = control,
    response_component = response_component
  )
}


#' @keywords internal
resolve_feature_test_time_values <- function(fit, time_values, window) {
  if (!is.null(time_values)) return(time_values)
  if (!is.null(window)) return(seq(window[[1]], window[[2]], length.out = 50L))
  NULL
}


#' @keywords internal
feature_descriptive_audit <- function(
    ds, assay, feature, group, levels, window, detection_threshold
) {
  values <- ds$assays[[assay]][[feature]]
  metadata <- ds$metadata
  if (!is.null(window)) {
    time_col <- ds$time_col
    keep <- metadata[[time_col]] >= window[[1]] & metadata[[time_col]] <= window[[2]]
    values <- values[keep]
    metadata <- metadata[keep, , drop = FALSE]
  }
  detected <- is.finite(values) & values > detection_threshold
  overall_n <- sum(is.finite(values))
  out <- list(
    n_detected = sum(detected),
    detection_rate = if (overall_n) sum(detected) / overall_n else NA_real_,
    min_group_n = NA_integer_,
    min_group_detection_rate = NA_real_,
    descriptive_difference = NA_real_,
    descriptive_summary = NA_character_
  )
  if (is.null(group) || !group %in% names(metadata)) return(out)

  groups <- as.character(metadata[[group]])
  levels <- levels %||% unique(groups[!is.na(groups)])
  summaries <- lapply(levels, function(level) {
    use <- groups == level & is.finite(values)
    positive <- use & values > detection_threshold
    n <- sum(use, na.rm = TRUE)
    data.frame(
      level = level,
      n = n,
      detected = sum(positive, na.rm = TRUE),
      detection_rate = if (n) sum(positive, na.rm = TRUE) / n else NA_real_,
      mean = if (n) mean(values[use], na.rm = TRUE) else NA_real_,
      positive_median = if (any(positive, na.rm = TRUE))
        stats::median(values[positive], na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  summaries <- do.call(rbind, summaries)
  out$min_group_n <- if (nrow(summaries)) min(summaries$n) else NA_integer_
  out$min_group_detection_rate <- if (any(is.finite(summaries$detection_rate)))
    min(summaries$detection_rate, na.rm = TRUE) else NA_real_
  out$descriptive_difference <- if (sum(is.finite(summaries$mean)) >= 2L)
    diff(range(summaries$mean, na.rm = TRUE)) else NA_real_
  out$descriptive_summary <- paste(sprintf(
    "%s:n=%d,det=%d,prev=%.3f,mean=%.4g,pos_median=%.4g",
    summaries$level, summaries$n, summaries$detected,
    summaries$detection_rate, summaries$mean, summaries$positive_median
  ), collapse = "; ")
  out
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

  if (hypothesis %in% c(
    "trajectory_difference", "trajectory_shape_difference",
    "persistent_level_difference", "level_or_shape",
    "critical_period_level", "occurrence_difference",
    "positive_abundance_difference"
  )) {
    if (design$type != "varying") {
      cli::cli_abort(
        "`hypothesis = '{hypothesis}'` requires `varying_trajectory()`."
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

  if (hypothesis %in% c(
    "critical_period_level", "occurrence_difference",
    "positive_abundance_difference"
  )) {
    if (is.null(test_covariate) || length(test_covariate) != 1L) {
      cli::cli_abort(
        "`hypothesis = '{hypothesis}'` requires one `test_covariate`."
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
    trajectory_shape_difference = paste(
      "Group trajectories have the same shape after centering over the declared",
      "centering domain; constant group offsets are allowed under the null."
    ),
    critical_period_level = paste0(
      "Adjusted abundance does not differ by `", test_covariate,
      "` inside the declared critical-period window."
    ),
    occurrence_difference = paste0(
      "Adjusted detection probability does not differ by `", test_covariate,
      "` in the declared analysis domain."
    ),
    positive_abundance_difference = paste0(
      "Among detected observations, adjusted abundance does not differ by `",
      test_covariate, "` in the declared analysis domain."
    ),
    persistent_level_difference = paste(
      "The uncentered time-integrated contrast between fitted group trajectories is zero."
    ),
    level_or_shape = paste(
      "The fitted group trajectories do not differ in level or shape over the declared domain."
    ),
    covariate_effect = paste0(
      "The adjusted response does not differ by `",
      test_covariate,
      "` under a shared longitudinal trajectory."
    )
  )
}


#' @keywords internal
feature_test_method <- function(hypothesis, null_control = NULL) {
  if (hypothesis == "trajectory_difference") {
    return("global_trajectory_max_contrast")
  }
  if (hypothesis == "level_or_shape") {
    return("global_trajectory_max_contrast")
  }
  if (hypothesis == "persistent_level_difference") {
    return(paste0(
      "global_integrated_level_",
      null_control$method %||% "gaussian_approximation"
    ))
  }
  if (hypothesis == "trajectory_shape_difference") {
    return(paste0("trajectory_shape_", null_control$method))
  }
  if (hypothesis %in% c(
    "critical_period_level", "occurrence_difference",
    "positive_abundance_difference"
  )) {
    return(paste0(
      "smooth_window_integrated_",
      null_control$method %||% "gaussian_approximation"
    ))
  }
  "nested_gam_likelihood_ratio_approximate"
}


#' @keywords internal
feature_test_inference_warning <- function(
    hypothesis,
    null_control = NULL,
    failed = FALSE
) {
  if (hypothesis %in% c("trajectory_difference", "level_or_shape")) {
    return(if (failed) {
      "Joint Wald trajectory-contrast test could not be estimated for this feature."
    } else {
      paste(
        "Maximum standardized fitted-contrast test calibrated from the model",
        "coefficient covariance; pointwise p-values only localize separation."
      )
    })
  }
  if (hypothesis == "persistent_level_difference") {
    return(if (failed) {
      "The integrated fitted-trajectory level contrast could not be estimated."
    } else {
      integrated_level_inference_warning(null_control)
    })
  }
  if (hypothesis == "trajectory_shape_difference") {
    if (failed) {
      return("Design-aware trajectory shape null could not be estimated for this feature.")
    }
    if (null_control$method == "gaussian_approximation") {
      return(paste(
        "Fast conditional Gaussian coefficient approximation; it is not an",
        "empirical repeated-measures null and should be confirmed by bootstrap."
      ))
    }
    return(paste(
      "Empirical reduced-model null preserving subject-level dependence;",
      "inference is conditional on the declared model and observed sampling design."
    ))
  }
  if (hypothesis %in% c(
    "critical_period_level", "occurrence_difference",
    "positive_abundance_difference"
  )) {
    return(if (failed) {
      "The declared window-specific component model could not be estimated."
    } else {
      paste(
        "The trajectory is fitted using all informative longitudinal observations;",
        "the uncentered contrast is integrated only over the declared window.",
        integrated_level_inference_warning(null_control)
      )
    })
  }
  feature_test_warning()
}


#' @keywords internal
integrated_level_inference_warning <- function(null_control) {
  if (is.null(null_control) || null_control$method == "gaussian_approximation") {
    return(paste(
      "Calibration uses a conditional Gaussian coefficient approximation, not",
      "an empirical repeated-measures null."
    ))
  }
  if (null_control$method == "subject_label_permutation") {
    return(paste(
      "Calibration permutes the focal label among eligible subjects while",
      "preserving each subject's complete trajectory and observed group counts."
    ))
  }
  paste(
    "Calibration uses a reduced model with the focal group removed and",
    "resampling that preserves the observed subject-level sampling structure."
  )
}


#' @keywords internal
feature_test_fdr_resolution <- function(table, fdr_level) {
  valid <- is.finite(table$p_value)
  resolutions <- table$p_value_resolution[
    valid & is.finite(table$p_value_resolution)
  ]
  minimum_p <- if (length(resolutions)) min(resolutions) else NA_real_
  n_tested <- sum(valid)
  single_hit_min_q <- minimum_p * n_tested

  data.frame(
    features_tested = n_tested,
    fdr_level = fdr_level,
    minimum_attainable_p = minimum_p,
    single_hit_min_q = single_hit_min_q,
    single_hit_resolution_adequate = is.finite(single_hit_min_q) &&
      single_hit_min_q <= fdr_level,
    stringsAsFactors = FALSE
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
  if (!is.null(x$null_control)) {
    cli::cli_text("Null method: {.field {x$null_control$method}}")
    cli::cli_text("Successful null replicates are reported per feature.")
  }

  ok <- sum(x$table$status == "ok", na.rm = TRUE)
  cli::cli_text("Successful fits: {ok}")

  if (!is.null(x$fdr_resolution) &&
      is.finite(x$fdr_resolution$single_hit_min_q)) {
    cli::cli_text(
      "Best attainable single-hit q-value: {format(x$fdr_resolution$single_hit_min_q, digits = 3)}"
    )
  }

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

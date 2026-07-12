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
    p_adjust_method = "BH",
    keep_fits = FALSE,
    error_action = c("continue", "stop")
) {
  error_action <- match.arg(error_action)

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

  rows <- vector("list", length(features))
  fits <- vector("list", length(features))
  results <- vector("list", length(features))
  names(fits) <- features
  names(results) <- features

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
        control = control
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
    control = control
  )

  if (isTRUE(keep_fits)) {
    out$fits <- fits
    out$results <- results
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
    control
) {
  spec <- model_spec(
    response = assay_feature(assay, feature),
    subject = subject,
    time = time,
    covariates = covariates
  )

  traj <- trajectory(spec, design = design)
  mf <- model_frame(ds, traj, control = control)
  fit <- fit_trajectory(mf)
  results <- trajectory_results(fit)
  test <- extract_feature_test(results, time$column)

  row <- data.frame(
    assay = assay,
    feature = feature,
    status = "ok",
    p_value = test$p_value,
    test_statistic = test$test_statistic,
    test_term = test$test_term,
    n = results$n,
    r_sq = results$r_sq,
    deviance_explained = results$deviance_explained,
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  list(
    row = row,
    fit = fit,
    results = results
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
  candidates <- c("p-value", "p.value", "Pr(>F)", "Pr(>|t|)", "Pr(>|z|)")
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

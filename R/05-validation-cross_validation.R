#####
# Validation helpers
#####

#' Create grouped cross-validation folds
#'
#' Splits rows into folds without splitting groups across train/test sets.
#'
#' @param data A data.frame or vector.
#' @param group Grouping variable. For data.frames, either a column name or a
#'   vector with one value per row.
#' @param v Number of folds.
#' @param seed Optional random seed.
#'
#' @return A list of fold definitions.
#' @export
grouped_cv_folds <- function(data, group, v = 5L, seed = NULL) {
  if (is.data.frame(data)) {
    if (is.character(group) && length(group) == 1L) {
      if (!group %in% colnames(data)) {
        cli::cli_abort("Grouping column {.field {group}} was not found in `data`.")
      }
      group_values <- data[[group]]
    } else {
      group_values <- group
    }
    n <- nrow(data)
  } else {
    group_values <- data
    n <- length(group_values)
  }

  if (length(group_values) != n) {
    cli::cli_abort("`group` must have one value per row.")
  }

  groups <- unique(group_values)
  groups <- groups[!is.na(groups)]

  if (length(groups) < 2L) {
    cli::cli_abort("Grouped CV requires at least two non-missing groups.")
  }

  v <- as.integer(v)
  if (length(v) != 1L || is.na(v) || v < 2L) {
    cli::cli_abort("`v` must be an integer >= 2.")
  }
  v <- min(v, length(groups))

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  groups <- sample(groups)
  fold_id <- rep(seq_len(v), length.out = length(groups))

  lapply(seq_len(v), function(i) {
    test_groups <- groups[fold_id == i]
    test <- which(group_values %in% test_groups)
    train <- setdiff(seq_len(n), test)

    list(
      fold = i,
      train = train,
      test = test,
      test_groups = test_groups
    )
  })
}


#' Cross-validate a fitted trajectory specification
#'
#' Runs grouped cross-validation on a `levaim_model_frame`. By default, folds are
#' grouped by the declared subject identifier.
#'
#' @param mf A `levaim_model_frame`.
#' @param v Number of folds.
#' @param group Optional grouping column. Defaults to the model subject column.
#' @param seed Optional random seed.
#' @param metrics Metrics to compute. Currently `"rmse"` and `"mae"`.
#' @param repeats Number of repeated grouped CV partitions.
#' @param keep_predictions Logical. Whether to retain held-out predictions.
#' @param ... Additional arguments passed to `predict()`.
#'
#' @return A `levaim_cv` object.
#' @export
cross_validate_trajectory <- function(
    mf,
    v = 5L,
    group = NULL,
    seed = NULL,
    metrics = c("rmse", "mae"),
    repeats = 1L,
    keep_predictions = TRUE,
    ...
) {
  if (!inherits(mf, "levaim_model_frame")) {
    cli::cli_abort("`mf` must be created with `model_frame()`.")
  }

  if (is.null(group)) {
    group <- mf$trajectory$spec$subject$column
  }

  metrics <- match.arg(metrics, c("rmse", "mae"), several.ok = TRUE)
  repeats <- as.integer(repeats)

  if (length(repeats) != 1L || is.na(repeats) || repeats < 1L) {
    cli::cli_abort("`repeats` must be an integer >= 1.")
  }

  if (!is.logical(keep_predictions) || length(keep_predictions) != 1L || is.na(keep_predictions)) {
    cli::cli_abort("`keep_predictions` must be `TRUE` or `FALSE`.")
  }

  fold_sets <- lapply(seq_len(repeats), function(repeat_id) {
    fold_seed <- if (is.null(seed)) {
      NULL
    } else {
      seed + repeat_id - 1L
    }

    folds <- grouped_cv_folds(mf$data, group = group, v = v, seed = fold_seed)
    lapply(folds, function(fold) {
      fold[["repeat"]] <- repeat_id
      fold
    })
  })

  folds <- unlist(fold_sets, recursive = FALSE)

  fold_outputs <- lapply(folds, function(fold) {
    evaluate_cv_fold(
      mf = mf,
      fold = fold,
      group = group,
      metrics = metrics,
      keep_predictions = keep_predictions,
      ...
    )
  })

  fold_metrics <- do.call(rbind, lapply(fold_outputs, `[[`, "metrics"))
  prediction_rows <- lapply(fold_outputs, `[[`, "predictions")
  prediction_rows <- prediction_rows[!vapply(prediction_rows, is.null, logical(1))]

  predictions <- if (length(prediction_rows) > 0L) {
    do.call(rbind, prediction_rows)
  } else {
    NULL
  }

  metric_summary <- summarize_cv_metrics(fold_metrics, metrics)

  structure(
    list(
      folds = folds,
      metrics = fold_metrics,
      summary = metric_summary,
      predictions = predictions,
      group = group,
      repeats = repeats,
      model_frame = mf
    ),
    class = "levaim_cv"
  )
}


#' @keywords internal
evaluate_cv_fold <- function(
    mf,
    fold,
    group,
    metrics,
    keep_predictions,
    ...
) {
  train_mf <- mf
  train_mf$data <- mf$data[fold$train, , drop = FALSE]

  test_data <- mf$data[fold$test, , drop = FALSE]
  fit <- fit_trajectory(train_mf)

  pred <- predict_for_validation(fit, test_data, ...)
  observed <- test_data$.y

  out <- data.frame(
    `repeat` = fold[["repeat"]],
    fold = fold$fold,
    n_train = length(fold$train),
    n_test = length(fold$test),
    check.names = FALSE
  )

  if ("rmse" %in% metrics) {
    out$rmse <- sqrt(mean((observed - pred)^2, na.rm = TRUE))
  }

  if ("mae" %in% metrics) {
    out$mae <- mean(abs(observed - pred), na.rm = TRUE)
  }

  predictions <- NULL
  if (isTRUE(keep_predictions)) {
    sample_id <- rownames(test_data)
    if (is.null(sample_id)) {
      sample_id <- as.character(fold$test)
    }

    predictions <- data.frame(
      `repeat` = fold[["repeat"]],
      fold = fold$fold,
      row = fold$test,
      sample_id = sample_id,
      group = as.character(test_data[[group]]),
      observed = observed,
      predicted = pred,
      residual = observed - pred,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  list(
    metrics = out,
    predictions = predictions
  )
}


#' @keywords internal
summarize_cv_metrics <- function(fold_metrics, metrics) {
  rows <- lapply(metrics, function(metric) {
    values <- fold_metrics[[metric]]

    data.frame(
      metric = metric,
      mean = mean(values, na.rm = TRUE),
      sd = stats::sd(values, na.rm = TRUE),
      min = min(values, na.rm = TRUE),
      max = max(values, na.rm = TRUE),
      n = sum(is.finite(values)),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


#' Extract CV predictions
#'
#' @param x A `levaim_cv` object.
#'
#' @return A data.frame of held-out predictions, or `NULL` if predictions were
#'   not retained.
#' @export
cv_predictions <- function(x) {
  if (!inherits(x, "levaim_cv")) {
    cli::cli_abort("`x` must be created with `cross_validate_trajectory()`.")
  }

  x$predictions
}


#' Summarize CV metrics
#'
#' @param x A `levaim_cv` object.
#'
#' @return A data.frame of aggregate CV metrics.
#' @export
cv_summary <- function(x) {
  if (!inherits(x, "levaim_cv")) {
    cli::cli_abort("`x` must be created with `cross_validate_trajectory()`.")
  }

  x$summary
}


#' @keywords internal
predict_for_validation <- function(fit, newdata, ...) {
  if (fit$engine == "spline") {
    prediction_data <- neutralize_random_effect_columns(fit$model_frame, newdata)

    return(as.numeric(stats::predict(
      fit$fit,
      newdata = prediction_data,
      exclude = spline_random_effect_terms(fit$model_frame),
      ...
    )))
  }

  as.numeric(predict(fit, newdata = newdata, ...))
}


#' @keywords internal
neutralize_random_effect_columns <- function(mf, newdata) {
  random_columns <- mf$trajectory$spec$subject$column
  random_columns <- c(
    random_columns,
    vapply(
      mf$trajectory$spec$covariates,
      function(cov) {
        if (inherits(cov, "levaim_random_effect")) {
          return(cov$column)
        }
        NA_character_
      },
      character(1)
    )
  )
  random_columns <- random_columns[!is.na(random_columns)]

  for (column in random_columns) {
    if (!column %in% colnames(newdata)) {
      next
    }

    reference <- mf$data[[column]]

    if (is.factor(reference)) {
      newdata[[column]] <- factor(as.character(reference[[1]]), levels = levels(reference))
    } else {
      newdata[[column]] <- reference[[1]]
    }
  }

  newdata
}


#' @keywords internal
spline_random_effect_terms <- function(mf) {
  covariates <- mf$trajectory$spec$covariates
  random_terms <- vapply(
    covariates,
    function(cov) {
      if (inherits(cov, "levaim_random_effect")) {
        return(paste0("s(", cov$column, ")"))
      }
      NA_character_
    },
    character(1)
  )
  random_terms <- random_terms[!is.na(random_terms)]

  c(
    paste0("s(", mf$trajectory$spec$subject$column, ")"),
    random_terms
  )
}


#' Print grouped cross-validation results
#'
#' @param x A `levaim_cv` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.levaim_cv <- function(x, ...) {
  cli::cli_h1("LEVAiM grouped cross-validation")
  cli::cli_text("Group: {.field {x$group}}")
  cli::cli_text("Repeats: {x$repeats}")
  cli::cli_text("Folds: {length(x$folds)}")
  print(x$summary)
  invisible(x)
}

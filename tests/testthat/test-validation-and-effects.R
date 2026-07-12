make_validation_example <- function(spline_basis = "gam", varying = TRUE) {
  set.seed(2)

  metadata <- expand.grid(
    subject = sprintf("S%02d", 1:10),
    week = c(0, 1, 2, 4, 8),
    KEEP.OUT.ATTRS = FALSE
  )
  metadata$sample_id <- sprintf("%s_W%s", metadata$subject, metadata$week)
  metadata$treatment <- rep(
    rep(c("control", "fiber"), each = 5),
    each = length(unique(metadata$week))
  )
  metadata$subject <- factor(metadata$subject)
  metadata$treatment <- factor(metadata$treatment)
  rownames(metadata) <- metadata$sample_id

  assay <- data.frame(
    feature_a = 1 +
      0.18 * metadata$week -
      0.01 * metadata$week^2 +
      ifelse(metadata$treatment == "fiber", 0.2 * log1p(metadata$week), 0) +
      rnorm(nrow(metadata), sd = 0.04)
  )
  rownames(assay) <- rownames(metadata)

  ds <- LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )

  spec <- model_spec(
    response = assay_feature("taxa", "feature_a"),
    subject = subject_var("subject"),
    time = time_var("week"),
    covariates = list(group_effect("treatment"))
  )

  traj <- if (varying) {
    trajectory(spec, varying_trajectory("treatment"))
  } else {
    trajectory(spec)
  }

  model_frame(
    ds,
    traj,
    trajectory_control(
      engine = "spline",
      spline_basis = spline_basis,
      spline_k = 4
    )
  )
}

test_that("grouped CV keeps subjects in one test fold", {
  mf <- make_validation_example()

  folds <- grouped_cv_folds(mf$data, group = "subject", v = 4, seed = 1)

  test_subjects <- unlist(lapply(folds, function(fold) {
    as.character(unique(mf$data$subject[fold$test]))
  }))

  expect_length(unique(test_subjects), nlevels(mf$data$subject))
  expect_equal(as.integer(sort(table(test_subjects))), rep(1L, nlevels(mf$data$subject)))
})

test_that("grouped CV fits and scores spline trajectories", {
  mf <- make_validation_example(varying = FALSE)

  cv <- cross_validate_trajectory(mf, v = 3, seed = 1)

  expect_s3_class(cv, "levaim_cv")
  expect_equal(nrow(cv$metrics), 3)
  expect_true(all(is.finite(cv$metrics$rmse)))
  expect_true(all(is.finite(cv$metrics$mae)))
  expect_equal(nrow(cv$summary), 2)
  expect_true(all(c("metric", "mean", "sd", "min", "max", "n") %in% colnames(cv$summary)))
  expect_equal(nrow(cv$predictions), nrow(mf$data))
  expect_true(all(c("observed", "predicted", "residual") %in% colnames(cv$predictions)))
  expect_s3_class(cv_summary(cv), "data.frame")
  expect_s3_class(cv_predictions(cv), "data.frame")
})

test_that("grouped CV supports repeated folds and optional prediction retention", {
  mf <- make_validation_example(varying = FALSE)

  cv <- cross_validate_trajectory(mf, v = 3, repeats = 2, seed = 1)

  expect_s3_class(cv, "levaim_cv")
  expect_equal(cv$repeats, 2)
  expect_equal(nrow(cv$metrics), 6)
  expect_equal(nrow(cv$predictions), nrow(mf$data) * 2)
  expect_equal(sort(unique(cv$metrics[["repeat"]])), c(1L, 2L))
  expect_equal(sort(unique(cv$predictions[["repeat"]])), c(1L, 2L))

  for (fold in cv$folds) {
    train_subjects <- unique(as.character(mf$data$subject[fold$train]))
    test_subjects <- unique(as.character(mf$data$subject[fold$test]))
    expect_length(intersect(train_subjects, test_subjects), 0)
  }

  cv_without_predictions <- cross_validate_trajectory(
    mf,
    v = 2,
    seed = 1,
    keep_predictions = FALSE
  )

  expect_null(cv_predictions(cv_without_predictions))
})

test_that("natural and cubic spline basis controls compile and fit", {
  natural_mf <- make_validation_example(spline_basis = "natural", varying = FALSE)
  cubic_mf <- make_validation_example(spline_basis = "cubic", varying = FALSE)

  expect_match(paste(deparse(natural_mf$formula), collapse = " "), "splines::ns", fixed = TRUE)
  expect_match(paste(deparse(cubic_mf$formula), collapse = " "), "splines::bs", fixed = TRUE)

  expect_s3_class(fit_trajectory(natural_mf), "levaim_fit")
  expect_s3_class(fit_trajectory(cubic_mf), "levaim_fit")
})

test_that("single-response models support multiple covariate features", {
  set.seed(3)

  metadata <- expand.grid(
    subject = sprintf("S%02d", 1:12),
    week = c(0, 1, 2, 4, 8),
    KEEP.OUT.ATTRS = FALSE
  )
  metadata$sample_id <- sprintf("%s_W%s", metadata$subject, metadata$week)
  metadata$treatment <- rep(
    rep(c("control", "fiber"), each = 6),
    each = length(unique(metadata$week))
  )
  metadata$delivery <- rep(
    rep(c("vaginal", "cesarean"), times = 6),
    each = length(unique(metadata$week))
  )
  metadata$baseline_diversity <- rep(runif(12, min = 0.2, max = 1), each = 5)
  metadata$batch <- factor(rep(c("B1", "B2", "B3"), length.out = nrow(metadata)))
  metadata$subject <- factor(metadata$subject)
  metadata$treatment <- factor(metadata$treatment)
  metadata$delivery <- factor(metadata$delivery)
  rownames(metadata) <- metadata$sample_id

  assay <- data.frame(
    feature_a = 1 +
      0.12 * metadata$week +
      ifelse(metadata$treatment == "fiber", 0.25, 0) +
      ifelse(metadata$delivery == "cesarean", -0.15, 0) +
      0.3 * metadata$baseline_diversity +
      rnorm(nrow(metadata), sd = 0.05)
  )
  rownames(assay) <- rownames(metadata)

  ds <- LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )

  spec <- model_spec(
    response = assay_feature("taxa", "feature_a"),
    subject = subject_var("subject"),
    time = time_var("week"),
    covariates = list(
      group_effect("treatment"),
      group_effect("delivery"),
      continuous_effect("baseline_diversity"),
      random_effect("batch")
    )
  )

  mf <- model_frame(
    ds,
    trajectory(
      spec,
      varying_trajectory(c("treatment", "delivery"))
    ),
    trajectory_control(engine = "spline", spline_k = 5)
  )

  formula_text <- paste(deparse(mf$formula), collapse = " ")
  fit <- fit_trajectory(mf)
  results <- trajectory_results(fit)

  expect_match(formula_text, "treatment", fixed = TRUE)
  expect_match(formula_text, "delivery", fixed = TRUE)
  expect_match(formula_text, "baseline_diversity", fixed = TRUE)
  expect_match(formula_text, "s(batch", fixed = TRUE)
  expect_match(formula_text, "interaction(treatment, delivery)", fixed = TRUE)

  expect_s3_class(fit, "levaim_fit")
  expect_s3_class(results, "levaim_results")
  expect_equal(results$response, "taxa / feature_a")
})

test_that("trajectory effects return prediction data for shared and varying fits", {
  shared_fit <- fit_trajectory(make_validation_example(varying = FALSE))
  varying_fit <- fit_trajectory(make_validation_example(varying = TRUE))

  direct_predictions <- predict(
    shared_fit,
    newdata = head(shared_fit$model_frame$data, 5)
  )

  shared_effects <- trajectory_effects(shared_fit, n = 12)
  varying_effects <- trajectory_effects(varying_fit, n = 12)

  expect_length(direct_predictions, 5)
  expect_true(all(is.finite(direct_predictions)))

  expect_equal(nrow(shared_effects), 12)
  expect_true(all(c(".fitted", ".se", ".lower", ".upper") %in% colnames(shared_effects)))
  expect_true(all(is.finite(shared_effects$.fitted)))

  expect_equal(nrow(varying_effects), 24)
  expect_true("treatment" %in% colnames(varying_effects))
  expect_equal(sort(unique(as.character(varying_effects$treatment))), c("control", "fiber"))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- plot_trajectory_effects(shared_fit, effects = shared_effects)
  expect_equal(nrow(plotted), nrow(shared_effects))
})

test_that("trajectory derivatives estimate local rates for shared and varying spline fits", {
  shared_fit <- fit_trajectory(make_validation_example(varying = FALSE))
  varying_fit <- fit_trajectory(make_validation_example(varying = TRUE))

  shared_derivatives <- trajectory_derivatives(shared_fit, n = 12)
  varying_derivatives <- trajectory_derivatives(varying_fit, n = 12)

  expect_equal(nrow(shared_derivatives), 12)
  expect_true(all(c(".derivative", ".lower_derivative", ".upper_derivative", ".sign") %in% colnames(shared_derivatives)))
  expect_true(all(is.finite(shared_derivatives$.derivative)))
  expect_true(all(shared_derivatives$.sign %in% c("increasing", "decreasing", "uncertain", "flat")))

  expect_equal(nrow(varying_derivatives), 24)
  expect_true("treatment" %in% colnames(varying_derivatives))
  expect_equal(sort(unique(as.character(varying_derivatives$treatment))), c("control", "fiber"))
  expect_true(all(is.finite(varying_derivatives$.derivative)))

  custom_time <- c(0, 2, 4, 8)
  custom_derivatives <- trajectory_derivatives(shared_fit, time_values = custom_time)
  expect_equal(custom_derivatives$week, custom_time)
})

test_that("trajectory derivatives are explicit about unsupported engines", {
  gp_fit <- structure(
    list(
      fit = NULL,
      model_frame = make_validation_example(varying = FALSE),
      engine = "gp"
    ),
    class = "levaim_fit"
  )

  expect_error(
    trajectory_derivatives(gp_fit, n = 5),
    "spline fits only"
  )
})

test_that("LEVAiM fits and results serialize with RDS round trips", {
  fit <- fit_trajectory(make_validation_example(varying = FALSE))
  results <- trajectory_results(fit)

  fit_path <- tempfile(fileext = ".rds")
  results_path <- tempfile(fileext = ".rds")

  expect_identical(save_levaim_object(fit, fit_path), fit_path)
  expect_identical(save_levaim_object(results, results_path), results_path)

  restored_fit <- load_levaim_object(fit_path)
  restored_results <- load_levaim_object(results_path)

  expect_s3_class(restored_fit, "levaim_fit")
  expect_s3_class(restored_results, "levaim_results")
  expect_equal(restored_fit$engine, fit$engine)
  expect_equal(restored_results$response, results$response)

  restored_predictions <- predict(
    restored_fit,
    newdata = head(restored_fit$model_frame$data, 4)
  )

  expect_length(restored_predictions, 4)
  expect_true(all(is.finite(restored_predictions)))

  expect_error(
    save_levaim_object(list(not = "levaim"), tempfile(fileext = ".rds")),
    "must be a LEVAiM"
  )

  invalid_path <- tempfile(fileext = ".rds")
  saveRDS(list(not = "levaim"), invalid_path)

  expect_error(
    load_levaim_object(invalid_path),
    "recognized LEVAiM object"
  )
})

test_that("model frames retain missing categorical annotations as explicit levels", {
  metadata <- data.frame(
    subject = factor(rep(sprintf("S%02d", 1:8), each = 4)),
    week = rep(c(0, 1, 2, 4), times = 8),
    feeding = c(rep("breast", 12), rep(NA_character_, 8), rep("formula", 11), "rare"),
    delivery = "vaginal",
    baseline = rep(seq(0.2, 0.9, length.out = 8), each = 4)
  )
  rownames(metadata) <- paste0(metadata$subject, "_", metadata$week)

  assay <- data.frame(
    feature_a = 1 + 0.1 * metadata$week + rnorm(nrow(metadata), sd = 0.02)
  )
  rownames(assay) <- rownames(metadata)

  ds <- LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )

  spec <- model_spec(
    response = assay_feature("taxa", "feature_a"),
    subject = subject_var("subject"),
    time = time_var("week"),
    covariates = list(
      group_effect("feeding"),
      group_effect("delivery"),
      continuous_effect("baseline", scale = TRUE)
    )
  )

  mf <- model_frame(
    ds,
    trajectory(spec, varying_trajectory("feeding")),
    trajectory_control(
      engine = "spline",
      sparse_min_n = 2,
      missing_level = "Unknown",
      sparse_other_level = "Other"
    )
  )

  formula_text <- paste(deparse(mf$formula), collapse = " ")

  expect_equal(nrow(mf$data), nrow(metadata))
  expect_true("Unknown" %in% levels(mf$data$feeding))
  expect_true("Other" %in% levels(mf$data$feeding))
  expect_true("feeding" %in% mf$trajectory$design$by)
  expect_true("delivery" %in% mf$preprocessing$dropped_covariates)
  expect_false(grepl("delivery", formula_text, fixed = TRUE))
  expect_s3_class(fit_trajectory(mf), "levaim_fit")
})

test_that("complete-case annotation handling remains available", {
  metadata <- data.frame(
    subject = factor(rep(sprintf("S%02d", 1:4), each = 3)),
    week = rep(c(0, 1, 2), times = 4),
    feeding = c("breast", NA, "formula", "breast", "formula", NA, rep("breast", 6))
  )
  rownames(metadata) <- paste0(metadata$subject, "_", metadata$week)

  assay <- data.frame(
    feature_a = 1 + 0.1 * metadata$week + rnorm(nrow(metadata), sd = 0.01)
  )
  rownames(assay) <- rownames(metadata)

  ds <- LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )

  spec <- model_spec(
    response = assay_feature("taxa", "feature_a"),
    subject = subject_var("subject"),
    time = time_var("week"),
    covariates = list(group_effect("feeding"))
  )

  mf <- model_frame(
    ds,
    trajectory(spec),
    trajectory_control(engine = "spline", na_action = "complete")
  )

  expect_equal(nrow(mf$data), nrow(metadata) - 2)
  expect_length(mf$preprocessing$dropped_rows, 2)
  expect_false("Unknown" %in% levels(mf$data$feeding))
})

test_that("assay feature testing fits one model per feature and adjusts p-values", {
  set.seed(7)

  metadata <- expand.grid(
    subject = sprintf("S%02d", 1:12),
    week = c(0, 1, 2, 4, 8),
    KEEP.OUT.ATTRS = FALSE
  )
  metadata$sample_id <- sprintf("%s_W%s", metadata$subject, metadata$week)
  metadata$treatment <- rep(
    rep(c("control", "fiber"), each = 6),
    each = length(unique(metadata$week))
  )
  metadata$subject <- factor(metadata$subject)
  metadata$treatment <- factor(metadata$treatment)
  rownames(metadata) <- metadata$sample_id

  assay <- data.frame(
    feature_time = 1 + 0.25 * metadata$week + rnorm(nrow(metadata), sd = 0.04),
    feature_group = 1 + ifelse(metadata$treatment == "fiber", 0.8, 0) + rnorm(nrow(metadata), sd = 0.04),
    feature_flat = rnorm(nrow(metadata), sd = 0.04),
    check.names = FALSE
  )
  rownames(assay) <- rownames(metadata)

  ds <- LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )

  batch <- test_assay_features(
    ds,
    assay = "taxa",
    features = c("feature_time", "feature_group", "feature_flat"),
    subject = "subject",
    time = "week",
    covariates = list(group_effect("treatment")),
    design = naive_trajectory(),
    control = trajectory_control(engine = "spline", spline_k = 4),
    keep_fits = TRUE
  )

  expect_s3_class(batch, "levaim_feature_tests")
  expect_equal(nrow(batch$table), 3)
  expect_true(all(batch$table$status == "ok"))
  expect_true(all(is.finite(batch$table$p_value)))
  expect_true(all(is.finite(batch$table$p_adjust)))
  expect_equal(batch$table$p_adjust, stats::p.adjust(batch$table$p_value, method = "BH"))
  expect_true(all(c("assay", "feature", "p_value", "p_adjust", "test_term") %in% colnames(batch$table)))
  expect_s3_class(batch$fits$feature_time, "levaim_fit")
  expect_s3_class(batch$results$feature_time, "levaim_results")
})

test_that("assay feature testing records feature-level errors when requested", {
  mf <- make_validation_example(varying = FALSE)
  ds <- LongitudinalDataset(
    assays = list(taxa = mf$data[".y"]),
    metadata = mf$data[, setdiff(names(mf$data), ".y"), drop = FALSE],
    time_col = "week",
    subject_col = "subject"
  )
  names(ds$assays$taxa) <- "feature_a"

  batch <- test_assay_features(
    ds,
    assay = "taxa",
    features = "feature_a",
    subject = "subject",
    time = "week",
    covariates = list(group_effect("treatment")),
    design = varying_trajectory("treatment"),
    control = trajectory_control(engine = "spline", drop_invariant_covariates = FALSE),
    error_action = "continue"
  )

  expect_equal(batch$table$status, "ok")

  expect_error(
    test_assay_features(
      ds,
      assay = "taxa",
      features = "missing_feature",
      subject = "subject",
      time = "week"
    ),
    "not found"
  )
})

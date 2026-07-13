stable_api_formals <- list(
  LongitudinalDataset = c("assays", "metadata", "time_col", "subject_col"),
  assay_feature = c("assay", "feature"),
  metadata_feature = c("column"),
  subject_var = c("column"),
  time_var = c("column"),
  group_effect = c("column", "reference", "role"),
  continuous_effect = c("column", "scale"),
  random_effect = c("column"),
  model_spec = c("response", "subject", "time", "covariates"),
  naive_trajectory = character(),
  varying_trajectory = c("by", "interaction"),
  trajectory = c("spec", "design"),
  trajectory_control = c(
    "engine", "family", "transform", "na_action", "spline_k",
    "spline_basis", "spline_method", "gp_kernel", "gp_basis_k", "chains", "iter",
    "cores", "adapt_delta", "max_treedepth", "seed", "sparse_min_n", "sparse_other_level", "missing_level",
    "drop_invariant_covariates"
  ),
  model_frame = c("ds", "traj", "control"),
  fit_trajectory = c("mf"),
  trajectory_results = c("fit"),
  grouped_cv_folds = c("data", "group", "v", "seed"),
  cross_validate_trajectory = c(
    "mf", "v", "group", "seed", "metrics", "repeats", "keep_predictions", "..."
  ),
  cv_predictions = c("x"),
  cv_summary = c("x"),
  trajectory_effects = c("fit", "n", "time_values", "se_fit", "type", "level", "..."),
  trajectory_contrasts = c("fit", "n", "time_values", "levels", "level"),
  trajectory_moments = c("fit", "n", "time_values", "windows", "levels", "draws", "level", "seed"),
  compare_trajectory_moments = c("moments", "reference", "p_adjust_method"),
  plot_trajectory_report = c("fit", "n", "time_values", "levels", "window", "points", "point_alpha", "palette"),
  plot_trajectory_contrasts = c("contrasts", "show_points", "palette"),
  plot_trajectory_moments = c("moments", "moment", "window", "feature", "palette"),
  plot_trajectory_moment_comparisons = c("comparisons", "moment", "window", "feature", "q_threshold"),
  trajectory_derivatives = c("fit", "n", "time_values", "method", "level", "type", "draws", "seed", "..."),
  plot_trajectory_effects = c("fit", "effects", "..."),
  test_assay_features = c(
    "ds", "assay", "features", "subject", "time", "covariates", "design",
    "control", "hypothesis", "test_covariate", "test_levels", "test_time_values",
    "test_simulations", "test_seed", "p_adjust_method", "keep_fits", "error_action"
  ),
  compare_trajectories = c(
    "fits", "n", "time_values", "windows", "view", "metrics", "standardize",
    "inference", "p_adjust_method"
  ),
  trajectory_distance_matrix = c("comparison", "metric", "matrix"),
  derive_time_varying_features = c("ds", "exposure", "summaries", "subject", "time"),
  exposure_window = c("start", "end", "method", "value"),
  exposure_ever = c("value"),
  exposure_last_observed = character(),
  exposure_transition = character(),
  exposure_derivations = c("ds")
)

test_that("stable exported APIs keep their core signatures", {
  namespace_exports <- getNamespaceExports("LEVAiM")

  for (fn in names(stable_api_formals)) {
    expect_true(fn %in% namespace_exports, info = fn)
    actual_formals <- names(formals(get(fn, envir = asNamespace("LEVAiM"))))
    if (is.null(actual_formals)) {
      actual_formals <- character()
    }

    expect_equal(actual_formals, stable_api_formals[[fn]], info = fn)
  }
})

test_that("stable S3 methods stay registered", {
  expect_true(is.function(getS3method("predict", "levaim_fit", optional = TRUE)))
  expect_true(is.function(getS3method("summary", "levaim_fit", optional = TRUE)))
  expect_true(is.function(getS3method("plot", "levaim_fit", optional = TRUE)))
  expect_true(is.function(getS3method("print", "levaim_cv", optional = TRUE)))
  expect_true(is.function(getS3method("print", "levaim_feature_tests", optional = TRUE)))
  expect_true(is.function(getS3method("print", "levaim_trajectory_comparison", optional = TRUE)))
})

test_that("stable spline API returns expected object classes", {
  metadata <- data.frame(
    subject = factor(rep(c("S1", "S2", "S3", "S4"), each = 4)),
    week = rep(c(0, 1, 2, 4), times = 4)
  )
  rownames(metadata) <- paste0(metadata$subject, "_", metadata$week)

  assay <- data.frame(
    feature_a = 1 + 0.2 * metadata$week + rnorm(nrow(metadata), sd = 0.01)
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
    time = time_var("week")
  )

  mf <- model_frame(
    ds,
    trajectory(spec),
    trajectory_control(engine = "spline", spline_k = 4)
  )

  fit <- fit_trajectory(mf)
  results <- trajectory_results(fit)
  effects <- trajectory_effects(fit, n = 5)
  derivatives <- trajectory_derivatives(fit, n = 5)
  cv <- cross_validate_trajectory(mf, v = 2, seed = 1)

  expect_s3_class(ds, "LongitudinalDataset")
  expect_s3_class(spec, "levaim_model_spec")
  expect_s3_class(mf, "levaim_model_frame")
  expect_s3_class(fit, "levaim_fit")
  expect_s3_class(results, "levaim_results")
  expect_s3_class(cv, "levaim_cv")
  expect_s3_class(effects, "data.frame")
  expect_s3_class(derivatives, "data.frame")
  expect_s3_class(cv_summary(cv), "data.frame")
  expect_s3_class(cv_predictions(cv), "data.frame")
})

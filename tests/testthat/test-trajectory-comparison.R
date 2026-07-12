make_comparison_dataset <- function() {
  metadata <- data.frame(
    subject = factor(rep(paste0("S", 1:8), each = 5)),
    week = rep(c(0, 1, 2, 3, 4), times = 8),
    feeding = factor(rep(c("breast", "formula"), each = 20))
  )
  rownames(metadata) <- paste0(metadata$subject, "_", metadata$week)

  offset <- ifelse(metadata$feeding == "formula", 0.3, 0)
  sample_index <- seq_len(nrow(metadata))
  assay <- data.frame(
    bug_a = 1 + 0.4 * metadata$week + offset + 0.04 * sin(sample_index),
    bug_b = 1.2 + 0.42 * metadata$week + offset * 0.9 + 0.04 * cos(sample_index),
    bug_c = 2.8 - 0.35 * metadata$week + offset * 0.2 + 0.04 * sin(sample_index * 0.7)
  )
  rownames(assay) <- rownames(metadata)

  LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )
}


test_that("compare_trajectories returns table-first summaries for retained batch fits", {
  ds <- make_comparison_dataset()

  batch <- test_assay_features(
    ds,
    assay = "taxa",
    features = c("bug_a", "bug_b", "bug_c"),
    subject = "subject",
    time = "week",
    design = naive_trajectory(),
    control = trajectory_control(engine = "spline", spline_k = 4),
    keep_fits = TRUE
  )

  comparison <- compare_trajectories(
    batch,
    n = 8,
    metrics = c("fitted_correlation", "fitted_distance", "derivative_correlation")
  )

  expect_s3_class(comparison, "levaim_trajectory_comparison")
  expect_s3_class(comparison, "data.frame")
  expect_named(
    comparison,
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
  )
  expect_true(any(comparison$metric == "fitted_distance"))
  expect_true(any(comparison$metric == "derivative_correlation"))
  expect_true(any(comparison$test_method == "curve_cor_test_descriptive"))
  expect_true(any(is.finite(comparison$q_value)))
  expect_true(any(grepl("descriptive", comparison$inference_warning)))

  distance <- trajectory_distance_matrix(comparison)
  expect_type(distance, "double")
  expect_equal(nrow(distance), 3)
  expect_equal(ncol(distance), 3)
  expect_equal(unname(diag(distance)), rep(0, 3))
})


test_that("compare_trajectories builds separate matrices within metadata levels", {
  ds <- make_comparison_dataset()

  batch <- test_assay_features(
    ds,
    assay = "taxa",
    features = c("bug_a", "bug_b", "bug_c"),
    subject = "subject",
    time = "week",
    covariates = list(group_effect("feeding")),
    design = varying_trajectory(by = "feeding"),
    control = trajectory_control(engine = "spline", spline_k = 4),
    keep_fits = TRUE
  )

  comparison <- compare_trajectories(
    batch,
    n = 8,
    view = "same_group_level",
    metrics = c("fitted_correlation", "fitted_distance")
  )

  expect_setequal(
    unique(comparison$group_level),
    c("breast", "formula")
  )
  distance <- trajectory_distance_matrix(comparison)
  expect_type(distance, "double")
  expect_true(all(comparison$group_variable == "feeding"))
  expect_false(any(comparison$feature_1 == comparison$feature_2))
})


test_that("compare_trajectories supports windows and non-inferential mode", {
  ds <- make_comparison_dataset()

  batch <- test_assay_features(
    ds,
    assay = "taxa",
    features = c("bug_a", "bug_b"),
    subject = "subject",
    time = "week",
    design = naive_trajectory(),
    control = trajectory_control(engine = "spline", spline_k = 4),
    keep_fits = TRUE
  )

  comparison <- compare_trajectories(
    batch,
    n = 9,
    windows = data.frame(
      label = c("early", "late"),
      start = c(0, 2),
      end = c(2, 4)
    ),
    metrics = "fitted_correlation",
    inference = "none"
  )

  expect_setequal(comparison$window_label, c("early", "late"))
  expect_true(all(is.na(comparison$p_value)))
  expect_true(all(comparison$test_method == "none"))
})


test_that("compare_trajectories rejects unavailable inputs clearly", {
  ds <- make_comparison_dataset()

  batch <- test_assay_features(
    ds,
    assay = "taxa",
    features = c("bug_a", "bug_b"),
    subject = "subject",
    time = "week",
    design = naive_trajectory(),
    control = trajectory_control(engine = "spline", spline_k = 4),
    keep_fits = FALSE
  )

  expect_error(
    compare_trajectories(batch),
    "keep_fits = TRUE"
  )

  batch_with_fits <- test_assay_features(
    ds,
    assay = "taxa",
    features = c("bug_a", "bug_b"),
    subject = "subject",
    time = "week",
    covariates = list(group_effect("feeding")),
    design = varying_trajectory(by = "feeding"),
    control = trajectory_control(engine = "spline", spline_k = 4),
    keep_fits = TRUE
  )

  expect_error(
    compare_trajectories(batch_with_fits, view = "overall"),
    "shared trajectories only"
  )

  expect_error(
    compare_trajectories(batch_with_fits, inference = "subject_permutation"),
    "not implemented yet"
  )
})

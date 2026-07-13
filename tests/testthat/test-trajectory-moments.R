make_moment_example <- function() {
  set.seed(22)
  metadata <- expand.grid(
    subject = sprintf("S%02d", 1:16),
    week = c(0, 1, 2, 4, 8),
    KEEP.OUT.ATTRS = FALSE
  )
  metadata$group <- ifelse(
    as.integer(sub("S", "", metadata$subject)) <= 8,
    "control",
    "fiber"
  )
  metadata$subject <- factor(metadata$subject)
  metadata$group <- factor(metadata$group)
  rownames(metadata) <- paste0(metadata$subject, "_", metadata$week)
  assay <- data.frame(
    feature_a = 0.4 + 0.15 * metadata$week - 0.015 * metadata$week^2 +
      ifelse(metadata$group == "fiber", 0.35 * exp(-((metadata$week - 4) / 2)^2), 0) +
      rnorm(nrow(metadata), sd = 0.05),
    row.names = rownames(metadata)
  )
  ds <- LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )
  spec <- model_spec(
    assay_feature("taxa", "feature_a"),
    subject_var("subject"),
    time_var("week"),
    list(group_effect("group"))
  )
  model_frame(
    ds,
    trajectory(spec, varying_trajectory("group")),
    trajectory_control(engine = "spline", spline_k = 5)
  )
}

test_that("trajectory moments summarize groups and windows in a long table", {
  fit <- fit_trajectory(make_moment_example())
  windows <- data.frame(
    label = c("early", "later"),
    start = c(0, 2),
    end = c(2, 8)
  )

  moments <- trajectory_moments(
    fit,
    n = 41,
    windows = windows,
    draws = 100,
    seed = 2
  )

  expect_s3_class(moments, "levaim_trajectory_moments")
  expect_equal(sort(unique(moments$group_level)), c("control", "fiber"))
  expect_equal(sort(unique(moments$window_label)), c("early", "later"))
  expect_true(all(c(
    "mean_level", "auc", "temporal_variance", "amplitude",
    "center_of_mass", "linear_trend", "mean_velocity",
    "velocity_variance", "max_growth_rate", "max_decline_rate",
    "start_level", "end_level", "net_change", "peak_value", "peak_time",
    "valley_value", "valley_time"
  ) %in% moments$moment))
  expect_true(all(c("estimate", "lower", "upper", "status") %in% names(moments)))
  expect_true(all(moments$status %in% c(
    "estimated", "interior_supported", "boundary_only",
    "flat_or_uncertain", "not_estimable"
  )))
  expect_equal(nrow(attr(moments, "moment_draws")), nrow(moments))
})

test_that("trajectory moment comparisons return draw-based differences and q-values", {
  fit <- fit_trajectory(make_moment_example())
  moments <- trajectory_moments(
    fit,
    time_values = seq(0, 8, length.out = 31),
    draws = 100,
    seed = 3
  )
  comparisons <- compare_trajectory_moments(moments, reference = "control")
  moment_plot <- plot_trajectory_moments(
    moments,
    moment = c("mean_level", "auc")
  )
  comparison_plot <- plot_trajectory_moment_comparisons(
    comparisons,
    moment = c("mean_level", "auc")
  )

  expect_s3_class(comparisons, "levaim_moment_comparisons")
  expect_true(all(comparisons$group_1 == "fiber"))
  expect_true(all(comparisons$group_2 == "control"))
  expect_true(all(c(
    "difference", "lower", "upper", "probability_positive", "p_value",
    "q_value", "test_method", "inference_warning", "status"
  ) %in% names(comparisons)))
  expect_true(all(comparisons$probability_positive >= 0 &
    comparisons$probability_positive <= 1))
  expect_true(all(comparisons$p_value >= 0 & comparisons$p_value <= 1))
  expect_true(all(comparisons$q_value >= 0 & comparisons$q_value <= 1))
  expect_s3_class(moment_plot, "ggplot")
  expect_s3_class(comparison_plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(moment_plot))
  expect_silent(ggplot2::ggplot_build(comparison_plot))
})

test_that("moment turning points distinguish interior and boundary extrema", {
  interior <- LEVAiM:::calculate_trajectory_moments(
    time = 0:4,
    values = c(0, 1, 3, 1, 0)
  )
  boundary <- LEVAiM:::calculate_trajectory_moments(
    time = 0:4,
    values = 0:4
  )
  flat <- LEVAiM:::calculate_trajectory_moments(
    time = 0:4,
    values = rep(1, 5)
  )

  expect_equal(interior$status[["peak_time"]], "interior_supported")
  expect_equal(boundary$status[["peak_time"]], "boundary_only")
  expect_equal(flat$status[["peak_time"]], "flat_or_uncertain")
  expect_equal(interior$values[["peak_time"]], 2)
})

test_that("moment reports round-trip through LEVAiM serialization", {
  fit <- fit_trajectory(make_moment_example())
  moments <- trajectory_moments(fit, n = 21, draws = 50, seed = 4)
  comparisons <- compare_trajectory_moments(moments)
  moment_path <- tempfile(fileext = ".rds")
  comparison_path <- tempfile(fileext = ".rds")

  save_levaim_object(moments, moment_path)
  save_levaim_object(comparisons, comparison_path)

  expect_s3_class(load_levaim_object(moment_path), "levaim_trajectory_moments")
  expect_s3_class(load_levaim_object(comparison_path), "levaim_moment_comparisons")
})

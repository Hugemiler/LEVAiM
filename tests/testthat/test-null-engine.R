make_shape_null_example <- function(kind = c("localized", "offset"), features = 1L) {
  kind <- match.arg(kind)
  metadata <- expand.grid(
    subject = sprintf("S%02d", seq_len(20)),
    week = 0:6,
    KEEP.OUT.ATTRS = FALSE
  )
  metadata$group <- factor(ifelse(
    as.integer(sub("S", "", metadata$subject)) <= 10,
    "A",
    "B"
  ))
  metadata$subject <- factor(metadata$subject)
  rownames(metadata) <- paste0(metadata$subject, "_", metadata$week)

  signal <- if (kind == "localized") {
    ifelse(
      metadata$group == "A",
      exp(-0.5 * ((metadata$week - 3) / 0.8)^2),
      0
    )
  } else {
    sin(metadata$week / 2) + ifelse(metadata$group == "A", 2, 0)
  }
  assay <- replicate(
    features,
    signal + stats::rnorm(nrow(metadata), sd = 0.12),
    simplify = "matrix"
  )
  colnames(assay) <- paste0("feature_", seq_len(features))
  rownames(assay) <- rownames(metadata)

  LongitudinalDataset(
    assays = list(taxa = as.data.frame(assay)),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )
}


run_shape_null_example <- function(ds, null_control, ...) {
  suppressWarnings(test_assay_features(
    ds,
    assay = "taxa",
    subject = "subject",
    time = "week",
    covariates = list(group_effect("group")),
    design = varying_trajectory("group"),
    control = trajectory_control(engine = "spline", spline_k = 6),
    hypothesis = "trajectory_shape_difference",
    test_covariate = "group",
    null_control = null_control,
    keep_fits = TRUE,
    ...
  ))
}


test_that("trajectory null controls validate robust and approximate engines", {
  control <- suppressWarnings(trajectory_null_control(
    method = "cluster_wild_bootstrap",
    simulations = 19,
    seed = 4,
    store_null = TRUE
  ))

  expect_s3_class(control, "levaim_trajectory_null_control")
  expect_equal(control$simulations, 19L)
  expect_true(control$store_null)
  expect_error(trajectory_null_control(simulations = 18), ">= 19")

  gaussian_default <- LEVAiM:::resolve_trajectory_null_control(NULL, "gaussian", 3)
  binomial_default <- LEVAiM:::resolve_trajectory_null_control(NULL, "binomial", 3)
  expect_equal(gaussian_default$method, "cluster_wild_bootstrap")
  expect_equal(binomial_default$method, "parametric_bootstrap")
  expect_equal(gaussian_default$simulations, 999L)
})


test_that("shape inference ignores constant offsets but detects localized change", {
  set.seed(19)
  offset <- run_shape_null_example(
    make_shape_null_example("offset"),
    suppressWarnings(trajectory_null_control(simulations = 19, seed = 4))
  )

  set.seed(18)
  localized <- run_shape_null_example(
    make_shape_null_example("localized"),
    suppressWarnings(trajectory_null_control(
      simulations = 19,
      seed = 4,
      store_null = TRUE
    )),
    test_time_values = seq(2.5, 3.5, length.out = 12)
  )

  expect_gt(offset$table$p_value, 0.1)
  expect_lt(abs(offset$table$null_centered_shape_difference), 0.01)
  expect_lte(localized$table$p_value, 0.1)
  expect_gt(localized$table$null_centered_shape_difference, 0)
  expect_equal(localized$table$null_method, "cluster_wild_bootstrap")
  expect_equal(localized$table$p_value_resolution, 1 / 20)
  expect_equal(nrow(localized$null_distributions$feature_1), 19)
  expect_match(localized$table$null_hypothesis, "constant group offsets")
  expect_equal(
    range(attr(localized$contrasts$feature_1, "centering_time_values")),
    c(0, 6)
  )
})


test_that("cluster wild bootstrap applies one multiplier per subject", {
  set.seed(5)
  ds <- make_shape_null_example("offset")
  spec <- model_spec(
    assay_feature("taxa", "feature_1"),
    subject_var("subject"),
    time_var("week"),
    list(group_effect("group"))
  )
  mf <- model_frame(
    ds,
    trajectory(spec, varying_trajectory("group")),
    trajectory_control(engine = "spline", spline_k = 6)
  )
  null_fit <- LEVAiM:::fit_shape_null_model(mf)$fit
  control <- suppressWarnings(trajectory_null_control(
    method = "cluster_wild_bootstrap",
    simulations = 19,
    seed = 8
  ))
  plan <- LEVAiM:::shape_null_resample_plan(null_fit, mf, control)
  fitted <- as.numeric(stats::predict(null_fit, type = "response"))
  residuals <- stats::residuals(null_fit, type = "response")
  usable <- abs(residuals) > 1e-10
  multipliers <- (plan[usable, 1] - fitted[usable]) / residuals[usable]
  by_subject <- split(multipliers, mf$data$subject[usable], drop = TRUE)

  expect_true(all(vapply(by_subject, function(x) length(unique(round(x, 8))) == 1L, logical(1))))
  expect_true(all(round(unlist(by_subject), 8) %in% c(-1, 1)))

  full_fit <- fit_trajectory(mf)
  full_components <- LEVAiM:::shape_contrast_components(full_fit)
  fast_components <- LEVAiM:::shape_null_statistics(
    full_fit$fit,
    full_components
  )
  expect_equal(fast_components$statistic, full_components$statistic)
  expect_equal(
    fast_components$integrated_shape_difference,
    full_components$integrated_shape_difference
  )
})


test_that("assay-wide null output reports Monte Carlo FDR resolution", {
  set.seed(21)
  ds <- make_shape_null_example("localized", features = 2)
  control <- suppressWarnings(trajectory_null_control(
    method = "gaussian_approximation",
    simulations = 19,
    seed = 3
  ))

  expect_warning(
    result <- test_assay_features(
      ds,
      assay = "taxa",
      subject = "subject",
      time = "week",
      covariates = list(group_effect("group")),
      design = varying_trajectory("group"),
      control = trajectory_control(engine = "spline", spline_k = 6),
      hypothesis = "trajectory_shape_difference",
      null_control = control
    ),
    "too coarse"
  )

  expect_equal(result$fdr_resolution$features_tested, 2)
  expect_equal(result$fdr_resolution$minimum_attainable_p, 0.05)
  expect_equal(result$fdr_resolution$single_hit_min_q, 0.1)
  expect_false(result$fdr_resolution$single_hit_resolution_adequate)
  expect_true(all(result$table$test_method ==
    "trajectory_shape_gaussian_approximation"))
  expect_true(all(grepl("not an empirical", result$table$inference_warning)))
  expect_true(all(is.finite(result$table$null_centered_shape_difference)))
})

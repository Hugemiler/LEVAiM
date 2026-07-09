test_that("spline trajectory pipeline fits and summarizes a feature", {
  set.seed(1)

  metadata <- expand.grid(
    subject = sprintf("S%02d", 1:8),
    week = c(0, 1, 2, 4, 8),
    KEEP.OUT.ATTRS = FALSE
  )

  metadata$sample_id <- sprintf("%s_W%s", metadata$subject, metadata$week)
  metadata$treatment <- rep(
    rep(c("control", "fiber"), each = 4),
    each = length(unique(metadata$week))
  )
  metadata$subject <- factor(metadata$subject)
  metadata$treatment <- factor(metadata$treatment)
  rownames(metadata) <- metadata$sample_id

  assay <- data.frame(
    feature_a = 1 +
      0.2 * metadata$week +
      ifelse(metadata$treatment == "fiber", 0.25 * log1p(metadata$week), 0) +
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
    covariates = list(group_effect("treatment"))
  )

  traj <- trajectory(
    spec,
    design = varying_trajectory("treatment")
  )

  mf <- model_frame(
    ds,
    traj,
    control = trajectory_control(engine = "spline")
  )

  fit <- fit_trajectory(mf)
  results <- trajectory_results(fit)

  expect_s3_class(mf, "levaim_model_frame")
  expect_s3_class(fit, "levaim_fit")
  expect_s3_class(results, "levaim_results")
  expect_equal(results$response, "taxa / feature_a")
  expect_equal(results$trajectory_type, "varying")
  expect_equal(results$n, nrow(metadata))
  expect_true(is.data.frame(results$smooth_terms))
})

test_that("default trajectory compiles to a shared smooth", {
  spec <- model_spec(
    response = metadata_feature("score"),
    subject = subject_var("subject"),
    time = time_var("week")
  )

  traj <- trajectory(spec)
  formula <- compile_formula(traj, engine = "spline")

  expect_equal(traj$design$type, "shared")
  expect_match(deparse(formula), "s\\(week\\)", fixed = FALSE)
  expect_match(deparse(formula), "s\\(subject, bs = \"re\"\\)", fixed = FALSE)
})

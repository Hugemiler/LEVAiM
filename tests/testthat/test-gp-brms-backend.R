make_gp_example <- function(kernel = "matern32") {
  set.seed(4)

  metadata <- expand.grid(
    subject = sprintf("S%02d", 1:4),
    week = c(0, 1, 2, 4),
    KEEP.OUT.ATTRS = FALSE
  )
  metadata$sample_id <- sprintf("%s_W%s", metadata$subject, metadata$week)
  metadata$subject <- factor(metadata$subject)
  rownames(metadata) <- metadata$sample_id

  assay <- data.frame(
    feature_a = sin(metadata$week / 2) + rnorm(nrow(metadata), sd = 0.05)
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

  model_frame(
    ds,
    trajectory(spec),
    trajectory_control(
      engine = "gp",
      gp_kernel = kernel,
      chains = 1,
      iter = 200,
      cores = 1,
      seed = 1
    )
  )
}

test_that("GP formulas compile requested brms covariance kernels", {
  rbf_mf <- make_gp_example(kernel = "rbf")
  matern_mf <- make_gp_example(kernel = "matern32")

  rbf_formula <- paste(deparse(rbf_mf$formula), collapse = " ")
  matern_formula <- paste(deparse(matern_mf$formula), collapse = " ")

  expect_match(rbf_formula, "gp\\(week, cov = \"exp_quad\"\\)", fixed = FALSE)
  expect_match(matern_formula, "gp\\(week, cov = \"matern32\"\\)", fixed = FALSE)
})

test_that("GP backend reports missing brms dependency clearly", {
  skip_if(requireNamespace("brms", quietly = TRUE))

  mf <- make_gp_example()

  expect_error(
    fit_trajectory(mf),
    "requires the .*brms.* package"
  )
})

test_that("brms GP backend can fit, predict, summarize, and produce effects", {
  skip_if_not_installed("brms")

  mf <- make_gp_example(kernel = "matern32")
  fit <- suppressWarnings(fit_trajectory(mf))

  pred <- predict(fit, newdata = head(mf$data, 3))
  results <- suppressWarnings(trajectory_results(fit))
  effects <- trajectory_effects(fit, n = 5)

  expect_s3_class(fit, "levaim_fit")
  expect_length(pred, 3)
  expect_s3_class(results, "levaim_results")
  expect_equal(results$engine, "gp")
  expect_equal(results$smooth_terms$kernel, "matern32")
  expect_equal(nrow(effects), 5)
  expect_true(all(c(".fitted", ".se", ".lower", ".upper") %in% colnames(effects)))
})

test_that("brms GP backend runs grouped cross-validation by subject", {
  skip_if_not_installed("brms")

  mf <- make_gp_example(kernel = "matern32")

  cv <- suppressWarnings(cross_validate_trajectory(
    mf,
    v = 2,
    seed = 1,
    allow_new_levels = TRUE
  ))

  expect_s3_class(cv, "levaim_cv")
  expect_equal(cv$group, "subject")
  expect_equal(nrow(cv$metrics), 2)
  expect_true(all(is.finite(cv$metrics$rmse)))
  expect_true(all(is.finite(cv$metrics$mae)))

  for (fold in cv$folds) {
    train_subjects <- unique(as.character(mf$data$subject[fold$train]))
    test_subjects <- unique(as.character(mf$data$subject[fold$test]))
    expect_length(intersect(train_subjects, test_subjects), 0)
  }
})

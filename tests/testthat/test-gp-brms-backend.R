make_gp_example <- function(
    kernel = "matern32",
    multiple_covariates = FALSE,
    family = "gaussian",
    basis_k = NULL
) {
  set.seed(4)

  metadata <- expand.grid(
    subject = sprintf("S%02d", 1:4),
    week = c(0, 1, 2, 4),
    KEEP.OUT.ATTRS = FALSE
  )
  metadata$sample_id <- sprintf("%s_W%s", metadata$subject, metadata$week)
  metadata$subject <- factor(metadata$subject)
  metadata$delivery <- factor(ifelse(
    metadata$subject %in% c("S01", "S02"),
    "vaginal",
    "cesarean"
  ))
  metadata$baseline <- seq(0.2, 0.8, length.out = 4)[as.integer(metadata$subject)]
  rownames(metadata) <- metadata$sample_id

  response <- switch(
    family,
    gaussian = sin(metadata$week / 2) + rnorm(nrow(metadata), sd = 0.05),
    binomial = stats::rbinom(
      nrow(metadata),
      1,
      stats::plogis(-0.5 + 0.4 * metadata$week)
    ),
    beta = pmin(
      pmax(stats::plogis(-1 + sin(metadata$week / 2)) + rnorm(nrow(metadata), sd = 0.02), 0.01),
      0.99
    )
  )
  assay <- data.frame(feature_a = response)
  rownames(assay) <- rownames(metadata)

  ds <- LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )

  covariates <- if (multiple_covariates) {
    list(group_effect("delivery"), continuous_effect("baseline"))
  } else {
    list()
  }
  spec <- model_spec(
    response = assay_feature("taxa", "feature_a"),
    subject = subject_var("subject"),
    time = time_var("week"),
    covariates = covariates
  )

  model_frame(
    ds,
    trajectory(spec),
    trajectory_control(
      engine = "gp",
      family = family,
      gp_kernel = kernel,
      gp_basis_k = basis_k,
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
  matern52_mf <- make_gp_example(kernel = "matern52")
  exponential_mf <- make_gp_example(kernel = "exponential", basis_k = 8)

  rbf_formula <- paste(deparse(rbf_mf$formula), collapse = " ")
  matern_formula <- paste(deparse(matern_mf$formula), collapse = " ")
  matern52_formula <- paste(deparse(matern52_mf$formula), collapse = " ")
  exponential_formula <- paste(deparse(exponential_mf$formula), collapse = " ")

  expect_match(rbf_formula, "gp\\(week, cov = \"exp_quad\"\\)", fixed = FALSE)
  expect_match(matern_formula, "gp\\(week, cov = \"matern32\"\\)", fixed = FALSE)
  expect_match(matern52_formula, "gp\\(week, cov = \"matern52\"\\)", fixed = FALSE)
  expect_match(
    exponential_formula,
    "gp\\(week, cov = \"exponential\", k = 8\\)",
    fixed = FALSE
  )
  expect_identical(exponential_mf$control$gp$approximation, "hsgp")
  expect_identical(exponential_mf$control$gp$basis_k, 8L)
  expect_equal(matern_mf$control$gp$adapt_delta, 0.95)
  expect_identical(matern_mf$control$gp$max_treedepth, 12L)
})

test_that("approximate and additional-kernel GP fits execute", {
  skip_if_not_installed("brms")

  mf <- make_gp_example(kernel = "exponential", basis_k = 8)
  fit <- suppressWarnings(fit_trajectory(mf))
  results <- suppressWarnings(trajectory_results(fit))

  expect_s3_class(fit, "levaim_fit")
  expect_equal(results$smooth_terms$kernel, "exponential")
  expect_equal(results$smooth_terms$approximation, "hsgp")
  expect_identical(results$smooth_terms$basis_k, 8L)
  expect_identical(fit$model_frame$control$gp$approximation, "hsgp")
  expect_length(predict(fit, newdata = head(mf$data, 2)), 2)
})

test_that("GP backend fits Bernoulli and beta response families", {
  skip_if_not_installed("brms")

  binomial_mf <- make_gp_example(family = "binomial", basis_k = 8)
  beta_mf <- make_gp_example(family = "beta", basis_k = 8)
  binomial_fit <- suppressWarnings(fit_trajectory(binomial_mf))
  beta_fit <- suppressWarnings(fit_trajectory(beta_mf))
  binomial_prediction <- predict(
    binomial_fit,
    newdata = head(binomial_mf$data, 3),
    type = "response"
  )
  beta_prediction <- predict(
    beta_fit,
    newdata = head(beta_mf$data, 3),
    type = "response"
  )

  expect_equal(binomial_fit$model_frame$control$family, "binomial")
  expect_equal(beta_fit$model_frame$control$family, "beta")
  expect_true(all(binomial_prediction >= 0 & binomial_prediction <= 1))
  expect_true(all(beta_prediction > 0 & beta_prediction < 1))
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
  derivatives <- trajectory_derivatives(fit, n = 7, draws = 20, seed = 2)
  moments <- trajectory_moments(fit, n = 7, draws = 20, seed = 2)

  expect_s3_class(fit, "levaim_fit")
  expect_length(pred, 3)
  expect_s3_class(results, "levaim_results")
  expect_equal(results$engine, "gp")
  expect_equal(results$smooth_terms$kernel, "matern32")
  expect_true(all(c(
    "max_rhat", "min_bulk_ess", "min_tail_ess", "divergences",
    "convergence_ok"
  ) %in% colnames(results$diagnostics)))
  expect_equal(nrow(effects), 5)
  expect_true(all(c(".fitted", ".se", ".lower", ".upper") %in% colnames(effects)))
  expect_equal(nrow(derivatives), 7)
  expect_true(all(c(
    ".derivative", ".lower_derivative", ".upper_derivative",
    ".probability_positive", ".sign", ".derivative_method", ".draws"
  ) %in% colnames(derivatives)))
  expect_true(all(is.finite(derivatives$.derivative)))
  expect_true(all(derivatives$.probability_positive >= 0 & derivatives$.probability_positive <= 1))
  expect_identical(unique(derivatives$.derivative_method), "posterior_finite_difference")
  expect_identical(unique(derivatives$.draws), 20L)
  expect_s3_class(moments, "levaim_trajectory_moments")
  expect_true(all(c("mean_level", "auc", "peak_time") %in% moments$moment))
})

test_that("fitted RBF GP supports multiple covariates, plotting, and serialization", {
  skip_if_not_installed("brms")

  mf <- make_gp_example(kernel = "rbf", multiple_covariates = TRUE)
  fit <- suppressWarnings(fit_trajectory(mf))
  results <- suppressWarnings(trajectory_results(fit))
  path <- tempfile(fileext = ".rds")

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  plotted <- plot_trajectory_effects(fit, n = 5)
  save_levaim_object(fit, path)
  restored <- load_levaim_object(path)
  restored_prediction <- predict(restored, newdata = head(mf$data, 2))

  expect_match(paste(deparse(mf$formula), collapse = " "), "delivery", fixed = TRUE)
  expect_match(paste(deparse(mf$formula), collapse = " "), "baseline", fixed = TRUE)
  expect_equal(results$smooth_terms$kernel, "rbf")
  expect_equal(nrow(plotted), 5)
  expect_s3_class(restored, "levaim_fit")
  expect_length(restored_prediction, 2)
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

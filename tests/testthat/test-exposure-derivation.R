make_exposure_dataset <- function() {
  metadata <- data.frame(
    sample = paste0("sample_", seq_len(8)),
    subject = rep(c("S1", "S2"), each = 4),
    week = rep(c(0, 2, 4, 6), times = 2),
    feeding = c(
      "breast", "breast", "mixed", "formula",
      "formula", "formula", "formula", "formula"
    ),
    delivery = rep(c("vaginal", "csection"), each = 4),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample

  assay <- data.frame(
    feature_a = seq_len(8),
    feature_b = rev(seq_len(8))
  )
  rownames(assay) <- rownames(metadata)

  LongitudinalDataset(
    assays = list(taxa = assay),
    metadata = metadata,
    time_col = "week",
    subject_col = "subject"
  )
}


test_that("derive_time_varying_features adds audited subject-level summaries", {
  ds <- make_exposure_dataset()

  derived <- derive_time_varying_features(
    ds,
    exposure = "feeding",
    summaries = list(
      feeding_3mo = exposure_window(1, 3, method = "modal"),
      ever_breast = exposure_ever("breast"),
      last_feeding = exposure_last_observed(),
      feeding_transition = exposure_transition()
    )
  )

  expect_s3_class(derived, "LongitudinalDataset")
  expect_true(all(c(
    "feeding_3mo",
    "ever_breast",
    "last_feeding",
    "feeding_transition"
  ) %in% names(derived$metadata)))

  expect_equal(
    unique(derived$metadata$feeding_3mo[derived$metadata$subject == "S1"]),
    "breast"
  )
  expect_equal(
    unique(derived$metadata$ever_breast[derived$metadata$subject == "S1"]),
    "TRUE"
  )
  expect_equal(
    unique(derived$metadata$ever_breast[derived$metadata$subject == "S2"]),
    "FALSE"
  )
  expect_equal(
    unique(derived$metadata$last_feeding[derived$metadata$subject == "S1"]),
    "formula"
  )
  expect_equal(
    unique(derived$metadata$feeding_transition[derived$metadata$subject == "S1"]),
    "breast->mixed->formula"
  )

  audit <- exposure_derivations(derived)
  expect_s3_class(audit, "data.frame")
  expect_setequal(audit$derived_column, c(
    "feeding_3mo",
    "ever_breast",
    "last_feeding",
    "feeding_transition"
  ))
})


test_that("model frames record and validate covariate roles", {
  ds <- make_exposure_dataset()
  ds <- derive_time_varying_features(
    ds,
    exposure = "feeding",
    summaries = list(ever_breast = exposure_ever("breast"))
  )

  spec <- model_spec(
    response = assay_feature("taxa", "feature_a"),
    subject = subject_var("subject"),
    time = time_var("week"),
    covariates = list(
      group_effect("delivery", role = "subject_static"),
      group_effect("ever_breast", role = "derived_exposure")
    )
  )

  mf <- model_frame(
    ds,
    trajectory(spec, design = varying_trajectory("ever_breast")),
    trajectory_control(engine = "spline", spline_k = 4)
  )

  roles <- mf$preprocessing$covariate_roles
  expect_s3_class(roles, "data.frame")
  expect_true("delivery" %in% roles$column)
  expect_true("ever_breast" %in% roles$column)
  expect_equal(
    roles$declared_role[roles$column == "delivery"],
    "subject_static"
  )
  expect_equal(
    roles$declared_role[roles$column == "ever_breast"],
    "derived_exposure"
  )
})


test_that("subject-static role rejects covariates that vary within subject", {
  ds <- make_exposure_dataset()

  spec <- model_spec(
    response = assay_feature("taxa", "feature_a"),
    subject = subject_var("subject"),
    time = time_var("week"),
    covariates = list(group_effect("feeding", role = "subject_static"))
  )

  expect_error(
    model_frame(
      ds,
      trajectory(spec),
      trajectory_control(engine = "spline", spline_k = 4)
    ),
    "change within subject"
  )
})


test_that("time-varying role warns when used as current-state trajectory modifier", {
  ds <- make_exposure_dataset()

  spec <- model_spec(
    response = assay_feature("taxa", "feature_a"),
    subject = subject_var("subject"),
    time = time_var("week"),
    covariates = list(group_effect("feeding", role = "time_varying"))
  )

  expect_warning(
    model_frame(
      ds,
      trajectory(spec, design = varying_trajectory("feeding")),
      trajectory_control(engine = "spline", spline_k = 4)
    ),
    "current observed state"
  )
})

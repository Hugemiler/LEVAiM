test_that("BackhedF example files can seed the first trajectory pipeline", {
  example_file <- function(filename) {
    installed <- system.file(
      "extdata", filename,
      package = "LEVAiM"
    )

    if (nzchar(installed)) {
      return(installed)
    }

    test_path("..", "..", "examples", "input_files", filename)
  }

  metadata_path <- example_file("2026-07-08-BackhedF_Metadata.csv")
  taxa_path <- example_file(
    "BackhedF_2015_mpa_vOct22_CHOCOPhlAnSGB_202403_lon_subsp_RELAB_merged_profile.tsv"
  )

  skip_if_not(file.exists(metadata_path))
  skip_if_not(file.exists(taxa_path))

  metadata <- read.csv(
    metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  taxa <- read_metaphlan(
    taxa_path,
    input_type = "merged",
    tax_level = "species",
    include_unclassified = FALSE
  )

  expect_true(all(metadata$sample %in% rownames(taxa)))

  metadata <- metadata[metadata$sample %in% rownames(taxa), , drop = FALSE]
  rownames(metadata) <- metadata$sample
  taxa <- taxa[rownames(metadata), , drop = FALSE]

  metadata$subject <- factor(metadata$subject)
  feature <- names(sort(colMeans(taxa), decreasing = TRUE))[1]

  ds <- LongitudinalDataset(
    assays = list(taxa = taxa),
    metadata = metadata,
    time_col = "ageMonths",
    subject_col = "subject"
  )

  spec <- model_spec(
    response = assay_feature("taxa", feature),
    subject = subject_var("subject"),
    time = time_var("ageMonths")
  )

  mf <- model_frame(
    ds,
    trajectory(spec),
    trajectory_control(
      engine = "spline",
      transform = "log1p",
      spline_k = 3
    )
  )

  fit <- fit_trajectory(mf)
  results <- trajectory_results(fit)

  expect_equal(nrow(mf$data), nrow(metadata))
  expect_equal(results$n, nrow(metadata))
  expect_equal(results$trajectory_type, "shared")
  expect_true(nrow(results$smooth_terms) >= 1)
})


test_that("BackhedF critical-period trajectories recover a feeding benchmark", {
  example_file <- function(filename) {
    installed <- system.file("extdata", filename, package = "LEVAiM")
    if (nzchar(installed)) return(installed)
    test_path("..", "..", "examples", "input_files", filename)
  }
  metadata_path <- example_file("2026-07-08-BackhedF_Metadata.csv")
  taxa_path <- example_file(
    "BackhedF_2015_mpa_vOct22_CHOCOPhlAnSGB_202403_lon_subsp_RELAB_merged_profile.tsv"
  )
  skip_if_not(file.exists(metadata_path))
  skip_if_not(file.exists(taxa_path))

  metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  taxa <- read_metaphlan(
    taxa_path, input_type = "merged", tax_level = "species",
    include_unclassified = FALSE
  )
  shared <- intersect(metadata$sample, rownames(taxa))
  metadata <- metadata[match(shared, metadata$sample), , drop = FALSE]
  rownames(metadata) <- metadata$sample
  taxa <- taxa[rownames(metadata), , drop = FALSE]
  metadata$subject <- factor(metadata$subject)
  ds <- LongitudinalDataset(
    assays = list(taxa = taxa), metadata = metadata,
    time_col = "ageMonths", subject_col = "subject"
  )
  ds <- derive_time_varying_features(
    ds, "feeding_state",
    list(feeding_4mo = exposure_window(3, 5, method = "modal"))
  )
  features <- colnames(taxa)[grepl(
    "s__Bifidobacterium_(breve|longum)$", colnames(taxa)
  )]

  result <- suppressWarnings(test_assay_features(
    ds, assay = "taxa", features = features,
    subject = "subject", time = "ageMonths",
    covariates = list(group_effect("feeding_4mo", role = "derived_exposure")),
    design = varying_trajectory("feeding_4mo"),
    control = trajectory_control(transform = "log1p", spline_k = 3,
      na_action = "explicit", sparse_min_n = 5),
    hypothesis = "critical_period_level",
    test_covariate = "feeding_4mo",
    test_levels = c("ExclBreastFed", "ExclFormulaFed", "Mixed"),
    test_window = c(3, 5),
    null_control = trajectory_null_control(
      method = "gaussian_approximation", simulations = 4999, seed = 11
    ),
    keep_fits = TRUE,
    error_action = "stop"
  ))
  breve <- result$table[grepl("Bifidobacterium_breve$", result$table$feature), ]

  expect_equal(nrow(result$table), 2L)
  expect_true(all(result$table$status == "ok"))
  expect_true(all(result$table$n == 278L))
  expect_equal(range(result$contrasts[[features[[1]]]]$ageMonths), c(3, 5))
  expect_lt(breve$p_value, 0.01)
  expect_gt(breve$model_effect, 0)
  expect_equal(breve$effect_level_1, "ExclBreastFed")
  expect_equal(breve$effect_level_2, "ExclFormulaFed")
})

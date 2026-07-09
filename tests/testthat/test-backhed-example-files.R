test_that("BackhedF example files can seed the first trajectory pipeline", {
  example_file <- function(filename) {
    installed <- system.file(
      "examples", "input_files", filename,
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
      spline_k = 6
    )
  )

  fit <- fit_trajectory(mf)
  results <- trajectory_results(fit)

  expect_equal(nrow(mf$data), nrow(metadata))
  expect_equal(results$n, nrow(metadata))
  expect_equal(results$trajectory_type, "shared")
  expect_true(nrow(results$smooth_terms) >= 1)
})

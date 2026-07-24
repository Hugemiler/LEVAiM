make_estimand_audit_example <- function(features = 1L, signal = TRUE) {
  subjects <- sprintf("S%03d", seq_len(60))
  metadata <- expand.grid(
    age = c(0.1, 4.1, 12.1),
    subject = subjects,
    KEEP.OUT.ATTRS = FALSE
  )
  subject_group <- setNames(
    c(rep("breast", 40), rep("mixed", 12), rep("formula", 8)),
    subjects
  )
  metadata$group <- factor(
    subject_group[metadata$subject],
    levels = c("breast", "mixed", "formula")
  )
  metadata$subject <- factor(metadata$subject)
  rownames(metadata) <- paste(metadata$subject, metadata$age, sep = "_")

  group_effect_size <- if (signal) {
    ifelse(metadata$group == "breast", 1.2,
      ifelse(metadata$group == "mixed", 0.5, 0)
    ) * exp(-0.5 * ((metadata$age - 4) / 1.8)^2)
  } else {
    0
  }
  subject_intercept <- stats::rnorm(length(subjects), sd = 0.2)
  baseline <- 0.35 + 0.15 * sin(metadata$age / 3) +
    subject_intercept[match(metadata$subject, subjects)]
  assay <- replicate(
    features,
    pmax(0, baseline + group_effect_size + stats::rnorm(nrow(metadata), sd = 0.25))
  )
  colnames(assay) <- paste0("feature_", seq_len(features))
  rownames(assay) <- rownames(metadata)

  LongitudinalDataset(
    assays = list(taxa = as.data.frame(assay)),
    metadata = metadata,
    time_col = "age",
    subject_col = "subject"
  )
}


estimand_audit_args <- function(ds, ...) {
  out <- list(
    ds = ds,
    assay = "taxa",
    subject = "subject",
    time = "age",
    covariates = list(group_effect("group")),
    design = varying_trajectory("group"),
    control = trajectory_control(transform = "log1p", spline_k = 3),
    test_covariate = "group",
    test_levels = c("breast", "mixed", "formula"),
    test_window = c(3, 5),
    null_control = suppressWarnings(trajectory_null_control(
      method = "gaussian_approximation",
      simulations = 199,
      seed = 4
    )),
    error_action = "stop"
  )
  replacements <- list(...)
  out[names(replacements)] <- replacements
  out
}


test_that("critical-period inference tests smooth trajectories without window subsetting", {
  set.seed(31)
  ds <- make_estimand_audit_example(signal = TRUE)
  result <- do.call(test_assay_features, estimand_audit_args(
    ds,
    hypothesis = "critical_period_level",
    keep_fits = TRUE
  ))

  expect_equal(result$table$n, nrow(ds$metadata))
  expect_equal(range(result$contrasts$feature_1$age), c(3, 5))
  expect_match(result$table$full_formula, "s\\(age")
  expect_equal(
    result$table$test_method,
    "smooth_window_integrated_gaussian_approximation"
  )
  expect_lt(result$table$p_value, 0.05)
  expect_gt(result$table$model_effect, 0)
  expect_true(all(c(
    "descriptive_summary", "model_std_error", "model_conf_low",
    "model_conf_high", "detection_rate", "composition_policy"
  ) %in% names(result$table)))
})


test_that("occurrence and positive abundance remain longitudinal component models", {
  set.seed(32)
  ds <- make_estimand_audit_example(signal = TRUE)
  occurrence <- do.call(test_assay_features, estimand_audit_args(
    ds,
    hypothesis = "occurrence_difference"
  ))
  positive <- do.call(test_assay_features, estimand_audit_args(
    ds,
    hypothesis = "positive_abundance_difference"
  ))

  expect_equal(occurrence$table$response_component, "occurrence")
  expect_equal(occurrence$table$effect_scale, "log_odds")
  expect_equal(occurrence$table$n, nrow(ds$metadata))
  expect_equal(positive$table$response_component, "positive_abundance")
  expect_lt(positive$table$n, nrow(ds$metadata))
  expect_equal(occurrence$table$window_start, 3)
  expect_equal(positive$table$window_end, 5)
})


test_that("integrated-level empirical null removes the focal group but keeps time", {
  set.seed(33)
  ds <- make_estimand_audit_example(signal = TRUE)
  args <- estimand_audit_args(
    ds,
    hypothesis = "critical_period_level",
    null_control = suppressWarnings(trajectory_null_control(
      method = "cluster_wild_bootstrap",
      simulations = 19,
      seed = 2,
      store_null = TRUE
    ))
  )
  result <- suppressWarnings(do.call(test_assay_features, args))

  expect_equal(result$table$null_method, "cluster_wild_bootstrap")
  expect_equal(result$table$null_simulations, 19L)
  expect_match(result$table$reduced_formula, "s\\(age")
  expect_false(grepl("group", result$table$reduced_formula, fixed = TRUE))
  expect_equal(nrow(result$null_distributions$feature_1), 19L)
})


test_that("subject-label permutation preserves trajectories for static exposures", {
  set.seed(35)
  ds <- make_estimand_audit_example(signal = TRUE)
  result <- suppressWarnings(do.call(test_assay_features, estimand_audit_args(
    ds,
    hypothesis = "critical_period_level",
    null_control = suppressWarnings(trajectory_null_control(
      method = "subject_label_permutation",
      simulations = 19,
      seed = 5,
      store_null = TRUE
    ))
  )))

  expect_equal(result$table$null_method, "subject_label_permutation")
  expect_equal(result$table$null_simulations, 19L)
  expect_lte(result$table$p_value, 0.10)
  expect_equal(nrow(result$null_distributions$feature_1), 19L)

  expect_error(
    do.call(test_assay_features, estimand_audit_args(
      ds,
      hypothesis = "trajectory_shape_difference",
      test_window = NULL,
      null_control = suppressWarnings(trajectory_null_control(
        method = "subject_label_permutation", simulations = 19
      ))
    )),
    "not valid for the centered-shape null"
  )
})


test_that("taxonomic aggregation and composition policies are explicit", {
  assay <- data.frame(
    "k__Bacteria|g__A|s__one" = c(2, 3),
    "k__Bacteria|g__A|s__two" = c(1, 1),
    "k__Bacteria|g__B|s__three" = c(1, 2),
    check.names = FALSE
  )
  rownames(assay) <- c("S1", "S2")
  genus <- aggregate_taxa(
    assay, tax_level = "genus", composition_policy = "renormalize"
  )

  expect_equal(names(genus), c("g__A", "g__B"))
  expect_equal(unname(rowSums(genus)), c(100, 100))
  expect_equal(attr(genus, "composition_policy"), "renormalize")
  expect_error(aggregate_taxa(data.frame(x = 1)), "No .* identifiers")

  community <- community_features(assay)
  expect_equal(
    names(community),
    c("richness", "shannon", "simpson", "dominance", "profiled_mass")
  )
  expect_equal(community$richness, c(3, 3))
  expect_equal(community$profiled_mass, c(4, 6))
  expect_true(all(community$shannon > 0))
  expect_true(all(community$simpson > 0 & community$simpson < 1))
  expect_true(all(community$dominance > 0 & community$dominance < 1))
})


test_that("null plasmode does not create widespread feeding discoveries", {
  set.seed(34)
  ds <- make_estimand_audit_example(features = 8L, signal = FALSE)
  result <- suppressWarnings(do.call(test_assay_features, estimand_audit_args(
    ds,
    hypothesis = "critical_period_level",
    null_control = suppressWarnings(trajectory_null_control(
      method = "cluster_wild_bootstrap",
      simulations = 19,
      seed = 8
    ))
  )))

  expect_true(all(result$table$status == "ok"))
  expect_gte(stats::median(result$table$p_value), 0.15)
  expect_lte(sum(result$table$p_value <= 0.10), 2L)
})

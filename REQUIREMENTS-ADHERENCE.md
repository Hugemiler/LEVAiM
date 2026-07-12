# LEVAiM Requirements Adherence Matrix

This document maps the D0.3 longitudinal microbiome modeling requirements to
current package evidence. It separates implementation progress from testing,
documentation, and approval so checkpoint claims remain auditable.

## Status Scale

- **Not started**: no user-facing implementation exists.
- **Prototype**: partial or internal implementation exists, but behavior is not
  stable enough for checkpoint credit.
- **Implemented**: user-facing implementation exists through exported functions
  or documented workflows.
- **Tested**: automated tests cover the expected behavior.
- **Documented**: README, vignette, or manual pages show intended use.
- **Approved**: reviewed and accepted by the designated executive-level reviewer.

## Current Evidence Inventory

- Package skeleton: `DESCRIPTION`, `NAMESPACE`, `R/`, `tests/`, `vignettes/`.
- Data parsing: `read_metaphlan()`, `read_humann()`, and
  `LongitudinalDataset()`.
- Preprocessing workflow: input parsing and alignment vignette in
  `vignettes/data-loading.Rmd`.
- Functional proof-of-concept vignettes use bundled real input files across
  `vignettes/model-specification.Rmd`, `vignettes/spline-fitting.Rmd`,
  `vignettes/validation.Rmd`, `vignettes/feature-testing.Rmd`,
  `vignettes/feeding-4mo-trajectories.Rmd`,
  `vignettes/trajectory-analysis.Rmd`, and
  `vignettes/gaussian-processes.Rmd`.
- Spline/GAM fitting: `trajectory_control(engine = "spline")`,
  `model_frame()`, `fit_trajectory()`, `trajectory_results()`.
- Sparse/missing annotation handling: `model_frame()` preprocesses categorical
  annotations consistently across spline and GP model frames, retaining missing
  categorical annotations as an explicit level by default and collapsing sparse
  levels to a configurable `"Other"` level.
- Time-varying exposure derivation: `derive_time_varying_features()` creates
  audited subject-level metadata from longitudinal exposures using
  `exposure_window()`, `exposure_ever()`, `exposure_last_observed()`, and
  `exposure_transition()`. `group_effect(role = ...)` lets `model_frame()`
  validate subject-static, time-varying, sample-level, and derived-exposure
  covariate usage.
- Per-feature workflow: `test_assay_features()` fits one model per assay
  feature, collates trajectory-level p-values, and applies FDR correction. The
  vignettes frame this as an all-feature screen before known or unexpected
  signals are interpreted.
- Trajectory comparison workflow: `compare_trajectories()` compares retained
  per-feature fits as fitted trajectory units and returns a long audit table
  with descriptive p-value labels and q-values. Heatmap-ready fitted distance
  matrices are derived from that table through `trajectory_distance_matrix()`.
- Spline basis controls: `trajectory_control(spline_basis = "gam")`,
  `trajectory_control(spline_basis = "natural")`, and
  `trajectory_control(spline_basis = "cubic")`.
- Grouped validation: `grouped_cv_folds()` and
  `cross_validate_trajectory()` support grouped folds, repeated CV,
  fold-level metrics, aggregate summaries, and optional held-out predictions.
- Prediction and plotting: `predict.levaim_fit()` and `plot.levaim_fit()`.
- Effect summaries and plots: `trajectory_effects()` and
  `plot_trajectory_effects()`.
- Derivative extraction: `trajectory_derivatives()` estimates finite-difference
  local rates of change for spline trajectory effects without making automatic
  peak/valley claims.
- Serialization: `save_levaim_object()` and `load_levaim_object()` round-trip
  recognized LEVAiM objects through RDS files.
- API lifecycle: `API-LIFECYCLE.md` plus signature and return-class tests in
  `tests/testthat/test-api-lifecycle.R`.
- Tests: `tests/testthat/test-pipeline.R` and
  `tests/testthat/test-backhed-example-files.R`,
  `tests/testthat/test-api-lifecycle.R`,
  `tests/testthat/test-validation-and-effects.R`, and
  `tests/testthat/test-gp-brms-backend.R`.
- Current verification: `devtools::test()` passes with 208 passing checks and
  one expected skip for the missing-`brms` branch when `brms` is installed;
  `R CMD check --no-manual` passes with one known NOTE for long example file
  names.

## D0.3.1 Gaussian Processes Mid-Development Checkpoint

| ID | Requirement | Current status | Evidence | Gap before checkpoint credit |
| --- | --- | --- | --- | --- |
| 0.3.1.1 | Package skeleton | Documented | Standard R package layout exists. | Approval signoff. |
| 0.3.1.2 | Data parsing and preprocessing pipeline | Tested + documented | MetaPhlAn/HUMAnN parsers, `LongitudinalDataset()`, sparse/missing annotation preprocessing, input-preparation vignette, model-specification vignette, Backhed example-file test. | Approval signoff. |
| 0.3.1.3 | GP regression with RBF/Matern32 kernels | Tested | `engine = "gp"` is wired to `brms::brm()` with `gp_kernel = "rbf"` mapped to `cov = "exp_quad"` and `gp_kernel = "matern32"` mapped to `cov = "matern32"`; formula behavior, missing-dependency behavior, and a fitted Matern32 `brms` smoke test are covered in `test-gp-brms-backend.R`. | Approval signoff; optionally add a fitted RBF smoke test for symmetric kernel coverage. |
| 0.3.1.4 | Grouped CV | Tested | `cross_validate_trajectory()` is engine-generic and supports grouped folds, repeated CV, aggregate summaries, and held-out prediction tables. Fitted GP grouped CV is covered in `test-gp-brms-backend.R` with held-out subjects and `allow_new_levels = TRUE`; grouped split behavior is also tested on spline model frames. | Approval signoff. |
| 0.3.1.5 | Prediction | Tested | GP `predict.levaim_fit()` dispatches through fitted `brms` objects; `test-gp-brms-backend.R` fits a GP and checks predictions. | Approval signoff. |
| 0.3.1.6 | Effect plots | Tested | `trajectory_effects()` and `plot_trajectory_effects()` support GP fits through fitted-value summaries; `test-gp-brms-backend.R` checks GP effect output columns and row counts after fitting. | Approval signoff; add explicit plot smoke test if required by reviewer. |
| 0.3.1.7 | Stable API | Prototype | Shared model specification/control API includes `engine = "gp"`; GP lifecycle remains prototype while broader GP coverage and reviewer feedback are pending. | Freeze GP-facing function signatures after fitted backend review. |

**Checkpoint summary:** 6/7 are tested or documented, and 1/7 remains prototype
pending GP API freeze. No items are marked approved.

## D0.3.2 Gaussian Processes End-of-Development Checkpoint

| ID | Requirement | Current status | Evidence | Gap before checkpoint credit |
| --- | --- | --- | --- | --- |
| 0.3.2.1 | Sparse GPs | Not started | No sparse GP implementation found. | Add sparse approximation strategy and tests. |
| 0.3.2.2 | Multi-feature fitting | Prototype | Interpreted as single-response modeling with multiple covariate features. GP formula compilation supports multiple covariates; fitted multi-covariate GP examples are not yet covered. | Add fitted multi-covariate GP examples/tests. |
| 0.3.2.3 | Classification + beta family | Not started | The brms GP backend currently gates execution to `family = "gaussian"`. | Add response-family support and examples/tests. |
| 0.3.2.4 | Additional kernels supported | Prototype | The brms GP backend supports the starting RBF/Matern32 set; additional kernels remain future work. | Add implemented kernels beyond the starting set. |
| 0.3.2.5 | Interpretability summaries | Prototype | `trajectory_results()` supports GP fits with kernel/covariance metadata and is covered by a fitted `brms` smoke test; `compare_trajectories()` can compare GP fitted effects but derivative summaries remain spline-only and GP-specific interpretation is not yet rich enough for end-development credit. | Add richer GP interpretability summaries beyond fitted-effect comparisons and kernel metadata. |
| 0.3.2.6 | High-level longitudinal microbiome workflows | Prototype | Functional real-file vignettes cover data loading, model specification, spline fitting, validation, feature testing, trajectory analysis, and GP setup. No fitted GP workflow vignette or GP-scale batch wrapper exists. | Add GP workflow wrappers/vignettes. |
| 0.3.2.7 | Serialization | Implemented | `save_levaim_object()` and `load_levaim_object()` support LEVAiM fits, results, model frames, datasets, specifications, trajectories, controls, and CV objects. Spline fit/result round trips are tested. | Add fitted GP serialization round-trip test before GP end-development credit. |

**Checkpoint summary:** 0/7 are ready for end-of-development credit. No items
are marked approved.

## D0.3.3 Spline Regression Mid-Development Checkpoint

| ID | Requirement | Current status | Evidence | Gap before checkpoint credit |
| --- | --- | --- | --- | --- |
| 0.3.3.1 | Package skeleton | Documented | Standard R package layout exists. | Approval signoff. |
| 0.3.3.2 | Data parsing and preprocessing pipeline | Tested + documented | MetaPhlAn/HUMAnN parsers, `LongitudinalDataset()`, sparse/missing annotation preprocessing, time-varying exposure derivation, `data-loading` and `model-specification` vignettes, Backhed example-file test. | Approval signoff. |
| 0.3.3.3 | Spline regression with natural/cubic starting set | Tested + documented | `spline_basis = "natural"` compiles through `splines::ns()` and `spline_basis = "cubic"` compiles through `splines::bs()`; both are covered by tests and generated manuals. | Approval signoff. |
| 0.3.3.4 | Grouped CV | Tested + documented | `grouped_cv_folds()` and `cross_validate_trajectory()` hold out whole subjects, support repeated CV, compute RMSE/MAE, retain held-out predictions, and expose `cv_summary()`/`cv_predictions()` accessors; tested in `test-validation-and-effects.R`. | Approval signoff. |
| 0.3.3.5 | Prediction | Tested + documented | `predict.levaim_fit()` is exported and directly tested on spline fits. | Approval signoff. |
| 0.3.3.6 | Effect plots | Tested + documented | `trajectory_effects()` returns inspectable prediction data and `plot_trajectory_effects()` plots it; both are tested. | Approval signoff. |
| 0.3.3.7 | Stable API | Tested + documented | `API-LIFECYCLE.md` classifies stable, experimental, and prototype surfaces; `test-api-lifecycle.R` guards stable exported signatures and return classes. | Approval signoff. |

**Checkpoint summary:** 7/7 now have tested or documented evidence. No items are
marked approved.

## D0.3.4 Spline Regression Later Checkpoint

| ID | Requirement | Current status | Evidence | Gap before checkpoint credit |
| --- | --- | --- | --- | --- |
| 0.3.4.1 | Additive models/GAM support | Tested + documented | `mgcv::gam()` backend, spline pipeline test, README, and the `spline-fitting` vignette. | Approval signoff. |
| 0.3.4.2 | Automatic smoothness selection | Tested + documented | Spline backend uses `method = "REML"` by default and `select = TRUE`; `spline_k` is capped by observed time points and spline controls are documented. | Approval signoff. |
| 0.3.4.3 | Multi-feature fitting | Tested + documented | Interpreted as single-response modeling with multiple covariate features. `model_spec(covariates = list(...))` supports categorical, continuous, and random-effect covariates; `varying_trajectory(by = c(...))` supports multiple declared trajectory modifiers; tests and the `model-specification` vignette cover this workflow. | Approval signoff. |
| 0.3.4.4 | Classification + beta family | Prototype | `family` control is accepted, but no classification/beta workflow or tests exist. | Validate supported families and add tests/examples. |
| 0.3.4.5 | Derivative extraction | Tested + documented | `trajectory_derivatives()` estimates finite-difference local rates of change for shared and varying spline trajectories; tests cover output shape, grouping columns, custom time grids, and unsupported GP engines. | Approval signoff; add GP posterior derivative support later if required. |
| 0.3.4.6 | Interpretability summaries | Tested + documented | `trajectory_results()` extracts fit metrics and smooth tables; `trajectory_derivatives()` summarizes local rates of change; `compare_trajectories()` returns feature-pair trajectory evidence as a long table with correlations, distances, derivative correlations, descriptive p-value labels, q-values, support flags, biological hypothesis text, and limitations. `trajectory_distance_matrix()` derives heatmap-ready matrices from that table. Tests cover table shape, metadata-level matrices, windows, non-inferential mode, and unavailable-input errors. | Approval signoff; add subject-aware permutation inference and plotting helpers if required by reviewer. |
| 0.3.4.7 | High-level microbiome workflows | Tested + documented | Functional real-file vignettes now cover data loading, assay-wide model specification, time-varying exposure derivation, all-feature per-feature testing, a feeding-around-4-months proof-of-concept, screen-derived spline follow-up, grouped validation, trajectory comparison, and GP setup; `test_assay_features()` provides an anpan_batch-style middle-layer workflow for all-feature assay testing with FDR correction; `compare_trajectories()` adds the downstream trajectory-comparison workflow for retained screen hits and metadata-level heatmap inputs. | Approval signoff; broaden workflow wrappers for additional common analysis patterns. |
| 0.3.4.8 | Serialization | Tested + documented | `save_levaim_object()` and `load_levaim_object()` serialize LEVAiM fits and results through RDS; spline fit/result round trips and post-load prediction are tested. | Approval signoff. |

**Checkpoint summary:** 7/8 are tested and documented, 1/8 is prototype. No
items are marked approved.

## Recommended Next Implementation Order

1. Review and sign off D0.3.3 spline mid-development evidence.
2. Add classification and beta-family validation/examples.
3. Add subject-aware permutation inference for trajectory comparisons.
4. Add fitted GP serialization round-trip coverage.
5. Add heatmap plotting helpers backed by `trajectory_distance_matrix()`.
6. Expand GP end-development work: sparse GP strategy, fitted multi-covariate
   GP examples, additional kernels, richer interpretability, and GP workflows.

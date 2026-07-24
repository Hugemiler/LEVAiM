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
  `vignettes/course-overview.Rmd`, `vignettes/model-specification.Rmd`,
  `vignettes/spline-fitting.Rmd`, `vignettes/validation.Rmd`,
  `vignettes/feeding-4mo-trajectories.Rmd`,
  `vignettes/response-families.Rmd`,
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
- Per-feature workflow: `test_assay_features()` fits one trajectory model per
  assay feature, tests explicit biological levels over a declared inference
  grid or critical-period window, localizes separation, and applies FDR
  correction. Critical-period level, occurrence, and positive-abundance tests
  fit full longitudinal smooths and integrate their uncentered contrasts over
  the declared window; they do not reduce the model to within-window
  proportions. Persistent-level and level-or-shape estimands remain available.
  Result tables combine fitted effects and uncertainty with raw sample counts,
  detection rates, descriptive group summaries, and explicit composition
  policy. `aggregate_taxa()` supplies genus and broader-rank sensitivity views;
  `community_features()` supplies prespecified richness, Shannon, Simpson,
  dominance, and profiled-mass endpoints. The centered shape hypothesis uses a
  reduced-model subject-cluster wild bootstrap
  for Gaussian responses, family-aware parametric bootstrap otherwise, and a
  separately labeled coefficient approximation for rapid screening. Result
  tables also support subject-label permutation for integrated-level hypotheses
  with subject-static focal exposures, preserving complete trajectories and
  observed group counts; that method is explicitly invalid for centered-shape
  nulls. Tables expose empirical p-value resolution, failed refits, integrated shape
  differences, and effects centered against their null expectation. The
  vignettes frame this as an outcome-blind, prevalence-defined feature screen
  before known or unexpected signals are interpreted; the executable feeding
  lecture tests the complete 66-species eligible universe and separates fast
  approximation from empirical-null confirmation.
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
  local rates for spline fitted effects and GP posterior fitted trajectories
  without making automatic peak/valley claims.
- Serialization: `save_levaim_object()` and `load_levaim_object()` round-trip
  recognized LEVAiM objects through RDS files.
- API lifecycle: `API-LIFECYCLE.md` plus signature and return-class tests in
  `tests/testthat/test-api-lifecycle.R`.
- Tests: `tests/testthat/test-pipeline.R` and
  `tests/testthat/test-backhed-example-files.R`,
  `tests/testthat/test-api-lifecycle.R`,
  `tests/testthat/test-validation-and-effects.R`,
  `tests/testthat/test-null-engine.R`, and
  `tests/testthat/test-gp-brms-backend.R`.
- Current verification: the full `devtools::test()` run passes 417 tests,
  including fitted RBF/Matern32 GP compile-and-sample paths, with one existing
  conditional GP skip and no failures or warnings. The revised feeding lecture
  rendered all assay-wide and 199-replicate community-null chunks end to end;
  all null refits succeeded. Source build succeeds, installed example data
  resolve through `extdata`, and code/documentation check stages pass. The
  check used `--no-vignettes --no-tests` after those stages were run separately,
  so its two warnings are the expected missing-`inst/doc` vignette warnings.

## D0.3.1 Gaussian Processes Mid-Development Checkpoint

| ID | Requirement | Current status | Evidence | Gap before checkpoint credit |
| --- | --- | --- | --- | --- |
| 0.3.1.1 | Package skeleton | Documented | Standard R package layout exists. | Approval signoff. |
| 0.3.1.2 | Data parsing and preprocessing pipeline | Tested + documented | MetaPhlAn/HUMAnN parsers, `LongitudinalDataset()`, sparse/missing annotation preprocessing, input-preparation vignette, model-specification vignette, Backhed example-file test. | Approval signoff. |
| 0.3.1.3 | GP regression with RBF/Matern32 kernels | Tested | `engine = "gp"` is wired to `brms::brm()` with fitted smoke tests for both RBF (`exp_quad`) and Matern32 covariances. | Approval signoff. |
| 0.3.1.4 | Grouped CV | Tested | `cross_validate_trajectory()` is engine-generic and supports grouped folds, repeated CV, aggregate summaries, and held-out prediction tables. Fitted GP grouped CV is covered in `test-gp-brms-backend.R` with held-out subjects and `allow_new_levels = TRUE`; grouped split behavior is also tested on spline model frames. | Approval signoff. |
| 0.3.1.5 | Prediction | Tested | GP `predict.levaim_fit()` dispatches through fitted `brms` objects; `test-gp-brms-backend.R` fits a GP and checks predictions. | Approval signoff. |
| 0.3.1.6 | Effect plots | Tested | `trajectory_effects()` and `plot_trajectory_effects()` are both exercised on fitted GP models. | Approval signoff. |
| 0.3.1.7 | Stable API | Tested + documented | `API-LIFECYCLE.md` freezes the 0.1.x spline and GP contracts; signature tests guard exported formals and return classes. The GP surface includes families, exact/approximate controls, kernels, diagnostics, prediction, interpretation, and serialization. | Approval signoff. |

**Checkpoint summary:** 7/7 are tested or documented and ready for acceptance
review. No items are marked approved.

## D0.3.2 Gaussian Processes End-of-Development Checkpoint

| ID | Requirement | Current status | Evidence | Gap before checkpoint credit |
| --- | --- | --- | --- | --- |
| 0.3.2.1 | Sparse GPs | Tested + documented | `gp_basis_k` requests the Hilbert-space approximate GP implemented by `brms::gp(k = ...)`; `NULL` retains the exact GP. Formula and fitted approximate-GP tests cover the path, and the GP vignette documents its approximation semantics and sensitivity obligation. | Approval signoff. |
| 0.3.2.2 | Multi-feature fitting | Tested + documented | Interpreted as single-response modeling with multiple covariate features. A fitted GP test includes categorical and continuous covariates, and the GP vignette documents the shared specification workflow. | Approval signoff. |
| 0.3.2.3 | Classification + beta family | Tested + documented | Fitted GP tests cover Bernoulli and beta likelihoods with response-scale predictions. The same explicit 0/1 and open-interval beta policy is enforced across engines and documented in `RESPONSE-FAMILIES.md` and the GP vignette. | Approval signoff. |
| 0.3.2.4 | Additional kernels supported | Tested + documented | Beyond the starting RBF/Matern32 set, Matern52 and exponential/Matern12 kernels compile through native `brms` covariance support; a fitted approximate exponential-GP test exercises an additional kernel end to end. | Approval signoff. |
| 0.3.2.5 | Interpretability summaries | Tested + documented | `trajectory_results()` reports kernel and sampler diagnostics; `trajectory_derivatives()` differentiates posterior fitted trajectories with credible intervals and sign probabilities; `trajectory_moments()` and `compare_trajectory_moments()` use posterior fitted draws. The executable `gaussian-processes` vignette fits the complete aligned Backhed cohort, enforces a diagnostic gate, and demonstrates prediction, plotting, derivatives, moments, and serialization. | Approval signoff. |
| 0.3.2.6 | High-level longitudinal microbiome workflows | Partial | The executable GP vignette fits the complete aligned Backhed cohort through parse, specification, compilation, diagnostics, prediction, posterior derivatives, moments, plotting, and serialization. The shared high-level API is complete for single-response follow-up; `test_assay_features()` remains spline-only for assay-scale batch inference. | Define and implement the computational and multiplicity contract for GP assay-wide batch inference. |
| 0.3.2.7 | Serialization | Tested + documented | Schema-1 envelopes record package version, class, and creation time; loaders reject unsupported future schemas and retain compatibility with legacy raw-object RDS files. Fitted spline and GP post-load prediction is tested. | Approval signoff. |

**Checkpoint summary:** 6/7 are tested and documented; 1/7 is partial because
GP assay-wide batch inference remains outside the frozen POC contract. No items
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
| 0.3.4.4 | Classification + beta family | Tested + documented | Spline binomial and beta models are fitted in tests. Model frames enforce 0/1 binomial responses and open-interval beta responses; zeros/ones are rejected rather than transformed implicitly. Grouped CV defaults to log loss/Brier score for binomial models. `RESPONSE-FAMILIES.md` fixes the boundary policy, and the `response-families` vignette demonstrates a real-data two-part occurrence/positive-abundance analysis. Explicit hurdle or zero/one-inflated likelihoods remain outside the current checkpoint scope. | Approval signoff. |
| 0.3.4.5 | Derivative extraction | Tested + documented | `trajectory_derivatives()` estimates finite-difference local rates of change for shared and varying spline trajectories; tests cover output shape, grouping columns, custom time grids, and unsupported GP engines. | Approval signoff; add GP posterior derivative support later if required. |
| 0.3.4.6 | Interpretability summaries | Tested + documented | `trajectory_results()` extracts fit metrics and smooth tables; `trajectory_derivatives()` summarizes local rates; `trajectory_contrasts()` localizes group differences; `trajectory_moments()` reports level, AUC, temporal variation, timing, trend, velocity, and conditional turning points with uncertainty; `compare_trajectory_moments()` returns long draw-based contrast tables; ggplot engines cover curves, localized contrasts, moments, and moment differences; feature-pair evidence and distance matrices remain available. | Approval signoff. |
| 0.3.4.7 | High-level microbiome workflows | Tested + documented | Functional real-file vignettes cover data loading, assay-wide model specification, time-varying exposure derivation, prevalence-defined per-feature testing, the feeding-around-4-months proof of concept, spline follow-up, grouped validation, trajectory comparison, and a fitted full-cohort GP follow-up; `test_assay_features()` provides an anpan_batch-style spline workflow with FDR correction; `compare_trajectories()` produces downstream trajectory-comparison and heatmap inputs. | Approval signoff; GP assay-wide batch inference remains the explicit D0.3.2.6 extension. |
| 0.3.4.8 | Serialization | Tested + documented | `save_levaim_object()` and `load_levaim_object()` serialize LEVAiM fits and results through RDS; spline fit/result round trips and post-load prediction are tested. | Approval signoff. |

**Checkpoint summary:** 8/8 are tested and documented. The checkpoint is ready
for acceptance review; no items are marked approved.

## Recommended Next Implementation Order

1. Review and sign off D0.3.1 and D0.3.3 mid-development evidence.
2. Review and sign off D0.3.4 spline later-checkpoint evidence.
3. Review the six completed D0.3.2 GP end-development rows.
4. Define whether GP assay-wide batch inference is required for D0.3.2.6 or is
   a post-POC scalability extension.
5. Add subject-aware inference for fitted-trajectory comparisons as analytical
   hardening after the acceptance-candidate release.

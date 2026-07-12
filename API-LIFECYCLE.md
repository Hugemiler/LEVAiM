# LEVAiM API Lifecycle

This document defines the stability level of exported LEVAiM functions for the
D0.3 longitudinal microbiome modeling checkpoints.

## Lifecycle Levels

- **Stable**: The function name, required arguments, return class, and primary
  behavior are intended to remain compatible within the 0.1.x series. Additive
  optional arguments are allowed when they do not change existing behavior.
- **Experimental**: The function is exported and tested, but names, defaults, or
  return details may change before a checkpoint is approved.
- **Prototype**: The function or argument exists to reserve an API direction, but
  the feature is not implemented enough for checkpoint credit.

Breaking changes to stable APIs require:

1. an entry in `NEWS.md`,
2. updated tests,
3. updated vignettes or documentation,
4. explicit project-owner approval.

## Stable API

These functions form the stable spline mid-development API:

| Function | Stability | Stable contract |
| --- | --- | --- |
| `LongitudinalDataset()` | Stable | Constructs a `LongitudinalDataset` from aligned assay tables and metadata. |
| `assay_feature()` | Stable | Declares an assay-backed response. |
| `metadata_feature()` | Stable | Declares a metadata-backed response. |
| `subject_var()` | Stable | Declares the subject identifier column. |
| `time_var()` | Stable | Declares the chronology column. |
| `group_effect()` | Stable | Declares a categorical covariate. |
| `continuous_effect()` | Stable | Declares a continuous covariate. |
| `random_effect()` | Stable | Declares a random-effect grouping covariate. |
| `model_spec()` | Stable | Combines response, subject, time, and covariates. |
| `naive_trajectory()` | Stable | Declares a shared time trajectory. |
| `varying_trajectory()` | Stable | Declares a covariate-varying time trajectory. |
| `trajectory()` | Stable | Combines a model specification with a trajectory design. |
| `trajectory_control()` | Stable | Declares engine, family, transformation, spline controls, GP controls, and sparse/missing annotation preprocessing controls. |
| `model_frame()` | Stable | Builds a backend-aware `levaim_model_frame`. |
| `fit_trajectory()` | Stable for `engine = "spline"` | Fits a spline trajectory and returns a `levaim_fit`. |
| `predict.levaim_fit()` | Stable for spline fits | Predicts from a fitted spline trajectory. |
| `summary.levaim_fit()` | Stable for spline fits | Returns the backend summary for spline fits. |
| `plot.levaim_fit()` | Stable for spline fits | Draws the backend spline plot. |
| `trajectory_results()` | Stable for spline fits | Returns a `levaim_results` summary. |
| `grouped_cv_folds()` | Stable | Creates grouped folds without splitting groups across train/test sets. |
| `cross_validate_trajectory()` | Stable for spline model frames | Runs grouped/repeated CV and returns `levaim_cv` with fold metrics, summaries, and optional held-out predictions. |
| `cv_summary()` | Stable for spline CV objects | Returns aggregate CV metric summaries. |
| `cv_predictions()` | Stable for spline CV objects | Returns held-out CV predictions when retained. |
| `trajectory_effects()` | Stable for spline fits | Returns an inspectable trajectory effect data frame. |
| `trajectory_derivatives()` | Stable for spline fits | Returns finite-difference local rate-of-change estimates without peak/valley interpretation. |
| `plot_trajectory_effects()` | Stable for spline fits | Plots trajectory effects and invisibly returns the effects data. |

## Experimental API

These functions are available but may change before later checkpoints:

| Function | Stability | Reason |
| --- | --- | --- |
| `read_metaphlan()` | Experimental | Parsing is functional and tested for current example files; additional MetaPhlAn variants may require argument/default changes. |
| `read_humann()` | Experimental | HUMAnN parsing is functional, but larger workflow coverage is still developing. |
| `save_levaim_object()` | Experimental | RDS-based serialization is tested for spline fits and results, but broader object-versioning policy is still developing. |
| `load_levaim_object()` | Experimental | RDS-based deserialization validates LEVAiM classes, but compatibility guarantees across object schema changes are still developing. |
| `test_assay_features()` | Experimental | Per-feature assay testing is tested for spline workflows, but additional test types and GP-scale performance patterns are still developing. |
| `compare_trajectories()` | Experimental | Trajectory comparison long tables are tested for spline workflows; subject-aware permutation inference and broader GP-scale examples are still developing. |
| `trajectory_distance_matrix()` | Experimental | Matrix extraction derives heatmap-ready matrices from current fitted-distance comparison tables; naming and matrix selection details may evolve with plotting helpers. |
| `derive_time_varying_features()` | Experimental | Exposure-history derivation is tested for initial window, ever, last-observed, and transition summaries; additional summary types and audit policies may evolve. |
| `exposure_window()` / `exposure_ever()` / `exposure_last_observed()` / `exposure_transition()` | Experimental | Constructors create lightweight exposure-summary specifications consumed by `derive_time_varying_features()`. |

## Prototype API

These controls are implemented enough for proof-of-concept use, but the GP
surface is not checkpoint-stable:

| Surface | Stability | Reason |
| --- | --- | --- |
| `trajectory_control(engine = "gp")` | Prototype | GP fitting is wired through `brms` and has fitted smoke-test coverage, but remains early and requires optional Stan tooling. |
| `trajectory_control(gp_kernel = ...)` | Prototype | Kernel choices are mapped to `brms` covariances for RBF and Matern32; broader kernel coverage and sparse approximations are future work. |
| `compile_formula(..., engine = "gp")` | Prototype | GP formula generation is tested for `brms` RBF/Matern32 terms, but the GP API is not yet stable. |
| `fit_trajectory(..., engine = "gp")` | Prototype | Fitted Matern32 GP smoke tests cover fit, prediction, results, effects, and grouped CV; additional families, sparse GP, and reviewer signoff remain open. |

## Proof-of-Concept Workflows

The package currently documents the workflow as real-file functional modules:

| Vignette | Role |
| --- | --- |
| `data-loading` | Input parsing, sample-name cleaning, sample alignment, and `LongitudinalDataset` construction. |
| `model-specification` | Declaring the assay-wide feeding trajectory screen: subject/time variables, covariate roles, time-varying exposure derivations, trajectory design, and preprocessing controls. |
| `feature-testing` | Testing every taxonomic feature with the same model, then labeling known early colonizers and unexpected signals after ranking. |
| `feeding-4mo-trajectories` | Worked proof-of-concept deriving feeding state from the 3-to-5 month window, screening taxonomic trajectories, interpreting hallmark and unexpected taxa, and generating trajectory-comparison inputs. |
| `spline-fitting` | Following up a screen-derived feature with results, effects, plots, and derivatives. |
| `validation` | Subject-grouped folds, repeated cross-validation, aggregate summaries, and held-out predictions for screen-derived follow-up models. |
| `trajectory-analysis` | Long trajectory-comparison tables and fitted-distance matrix generation for retained screen hits. |
| `gaussian-processes` | GP follow-up setup for screened trajectories, including RBF/Matern32 controls and brms/Stan fitting entry points. |

The proof-of-concept analysis layer follows a table-first rule. Per-feature
screening can retain fitted models with `keep_fits = TRUE`, and
`compare_trajectories()` turns those fits into a long comparison table with
feature pairs, metadata group labels, time windows, metrics, descriptive
p-values when requested, FDR-adjusted q-values, support flags, biological
hypothesis text, and limitations. Heatmap-ready distance matrices are derived
from that table through `trajectory_distance_matrix()` instead of replacing it.
Correlation p-values are explicitly labeled as descriptive fitted-curve tests
because fitted grid points are not independent observations.

## Sparse and Missing Annotations

Model frames apply the same annotation preprocessing before spline or GP
backends are called. Missing categorical annotations are retained as an explicit
level by default, sparse categorical levels are collapsed to a configurable
catch-all level, and invariant non-trajectory covariates can be dropped before
formula compilation. Essential fields such as response, time, subject, and
continuous covariates still require observed values unless future imputation
support is added.

## Signature Guardrails

The stable API signatures are guarded by automated tests in
`tests/testthat/test-api-lifecycle.R`. Those tests check exported function names,
core formal arguments, S3 method registration, and expected return classes for a
minimal spline pipeline.

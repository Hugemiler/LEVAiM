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

These functions form the frozen 0.1.x proof-of-concept API:

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
| `fit_trajectory()` | Stable for spline and GP engines | Fits a trajectory and returns a `levaim_fit`; GP supports exact or Hilbert-space approximate fitting and Gaussian, binomial, or beta families. |
| `predict.levaim_fit()` | Stable | Returns response-scale predictions from spline and GP fits. |
| `summary.levaim_fit()` | Stable | Returns the corresponding backend summary. |
| `plot.levaim_fit()` | Stable | Draws the corresponding backend plot. |
| `trajectory_results()` | Stable | Returns an engine-aware `levaim_results` summary. |
| `grouped_cv_folds()` | Stable | Creates grouped folds without splitting groups across train/test sets. |
| `cross_validate_trajectory()` | Stable for spline model frames | Runs grouped/repeated CV and returns `levaim_cv` with fold metrics, summaries, and optional held-out predictions. |
| `cv_summary()` | Stable for spline CV objects | Returns aggregate CV metric summaries. |
| `cv_predictions()` | Stable for spline CV objects | Returns held-out CV predictions when retained. |
| `trajectory_effects()` | Stable | Returns an inspectable trajectory effect data frame for spline and GP fits. |
| `trajectory_derivatives()` | Stable | Returns finite-difference local rates without peak/valley interpretation. GP fits differentiate posterior fitted trajectories and report pointwise credible intervals and sign probabilities. |
| `plot_trajectory_effects()` | Stable | Plots trajectory effects and invisibly returns the effects data. |
| `save_levaim_object()` / `load_levaim_object()` | Stable | Writes schema-versioned RDS envelopes, rejects unsupported future schemas, and reads legacy raw LEVAiM RDS objects. |

## Experimental API

These functions are available but may change before later checkpoints:

| Function | Stability | Reason |
| --- | --- | --- |
| `read_metaphlan()` | Experimental | Parsing is functional and tested for current example files; additional MetaPhlAn variants may require argument/default changes. |
| `aggregate_taxa()` | Experimental | Aggregates MetaPhlAn lineages to a declared rank with an explicit as-supplied or renormalized composition policy. |
| `community_features()` | Experimental | Derives prespecified richness, Shannon, Simpson, dominance, and profiled-mass endpoints for the shared longitudinal workflow. |
| `read_humann()` | Experimental | HUMAnN parsing is functional, but larger workflow coverage is still developing. |
| `trajectory_null_control()` | Experimental | Configures design-aware cluster wild bootstrap, subject-label permutation for subject-static integrated-level hypotheses, family-aware parametric bootstrap, or an explicitly labeled Gaussian coefficient approximation. Centered-shape tests reject label permutation because their null permits group offsets. |
| `test_assay_features()` | Experimental | Per-feature spline testing supports overall, any-separation, centered-shape, critical-period level, occurrence, positive-abundance, persistent-level, level-or-shape, and adjusted-covariate hypotheses. Critical-period windows localize contrasts of full longitudinal smooths rather than subsetting models to cross-sectional samples. Tables combine fitted effects, uncertainty, raw support, p/q-values, null diagnostics, and composition policy; GP-scale inference remains future work. |
| `trajectory_contrasts()` | Experimental | Returns covariance-aware pointwise group contrasts over time for localization after an omnibus trajectory test. |
| `trajectory_moments()` | Experimental | Returns engine-neutral functional summaries over declared windows using spline coefficient draws or GP posterior fitted draws. |
| `compare_trajectory_moments()` | Experimental | Returns draw-based between-group moment differences with intervals, tail probabilities, q-values, and explicit inferential warnings. |
| `plot_trajectory_report()` | Experimental | Returns a ggplot trajectory panel with raw observations, fitted curves, uncertainty, selected levels, and optional biological-window highlighting. |
| `plot_trajectory_contrasts()` | Experimental | Returns faceted time-localized fitted differences with pointwise uncertainty and zero-crossing markers. |
| `plot_trajectory_moments()` | Experimental | Returns faceted moment estimate and uncertainty panels from the long moments table. |
| `plot_trajectory_moment_comparisons()` | Experimental | Returns forest-style group-difference panels from the moment-comparison table. |
| `compare_trajectories()` | Experimental | Trajectory comparison long tables are tested for spline workflows; subject-aware permutation inference and broader GP-scale examples are still developing. |
| `trajectory_distance_matrix()` | Experimental | Matrix extraction derives heatmap-ready matrices from current fitted-distance comparison tables; naming and matrix selection details may evolve with plotting helpers. |
| `derive_time_varying_features()` | Experimental | Exposure-history derivation is tested for initial window, ever, last-observed, and transition summaries; additional summary types and audit policies may evolve. |
| `exposure_window()` / `exposure_ever()` / `exposure_last_observed()` / `exposure_transition()` | Experimental | Constructors create lightweight exposure-summary specifications consumed by `derive_time_varying_features()`. |

## GP Backend Boundary

The frozen GP surface supports exact and Hilbert-space approximate models,
Gaussian/binomial/beta families, and RBF, exponential, Matern-3/2, and
Matern-5/2 kernels through `brms`. GP-scale assay-wide batch inference remains
experimental because its computational and multiplicity contract is not yet
settled.

## Proof-of-Concept Workflows

The package currently documents the workflow as real-file functional modules:

| Vignette | Role |
| --- | --- |
| `course-overview` | Syllabus, recurring biological question, course order, feature-universe rule, and analytical contract. |
| `data-loading` | Input parsing, sample-name cleaning, sample alignment, and `LongitudinalDataset` construction. |
| `model-specification` | Declaring the assay-wide feeding trajectory screen: subject/time variables, covariate roles, time-varying exposure derivations, trajectory design, and preprocessing controls. |
| `feeding-4mo-trajectories` | All-eligible assay-wide critical-period screen using full smooth trajectories, two-part occurrence/positive-abundance sensitivity, composition and genus-rank checks, FDR correction, and post-ranking hallmark audit. |
| `spline-fitting` | Following up a screen-derived feature with results, effects, plots, and derivatives. |
| `validation` | Subject-grouped folds, repeated cross-validation, aggregate summaries, and held-out predictions for screen-derived follow-up models. |
| `trajectory-analysis` | Six-feature follow-up producing long trajectory-comparison tables and feeding-group-specific fitted-distance matrices. |
| `gp-crash-course` | Conceptual introduction to GP trajectory assumptions, kernels, likelihoods, approximation, diagnostics, derivatives, and biological interpretation. |
| `gaussian-processes` | Executable full-cohort Backhed GP follow-up with diagnostic gating, prediction, posterior derivatives, moments, plotting, and serialization. |
| `response-families` | Real-data spline workflow separating binomial occurrence from beta-distributed positive abundance under the explicit zero/one boundary policy. |

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

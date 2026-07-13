
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LEVAiM

<!-- badges: start -->

<!-- badges: end -->

LEVAiM is an R package for longitudinal microbiome modeling. It is being
built as a bioBakery-oriented extension for organizing microbial feature
tables, declaring time-aware models, and fitting feature trajectories
over repeated samples.

The package currently supports:

- MetaPhlAn and HUMAnN-style input parsing;
- longitudinal metadata and assay alignment;
- spline/GAM trajectory modeling;
- Gaussian-process trajectory modeling through `brms`;
- subject-grouped cross-validation;
- sparse and missing categorical annotation handling;
- derivation of subject-level features from time-varying exposures;
- assay-wide per-feature testing with FDR correction;
- prediction, trajectory effects, localized contrasts, functional moment
  reports, and compact result summaries.

## Vignettes

This README is the project orientation page. The operational
walkthroughs live in the package vignettes:

``` r
vignette("data-loading", package = "LEVAiM")
vignette("model-specification", package = "LEVAiM")
vignette("feature-testing", package = "LEVAiM")
vignette("feeding-4mo-trajectories", package = "LEVAiM")
vignette("spline-fitting", package = "LEVAiM")
vignette("validation", package = "LEVAiM")
vignette("trajectory-analysis", package = "LEVAiM")
vignette("gaussian-processes", package = "LEVAiM")
```

Use `data-loading` to move from files to a `LongitudinalDataset`. Then
use the modeling vignettes to declare an assay-wide feeding trajectory
screen, test all taxonomic features, and only then follow up known or
unexpected signals with fitting, validation, trajectory comparison, and
Gaussian-process setup. The `feeding-4mo-trajectories` vignette is the
worked proof-of-concept: it derives feeding state around 4 months,
screens taxonomic trajectories, and interprets the screen as a
trajectory-comparison problem.

If the vignettes have not been built yet, render them from the package
root:

``` r
rmarkdown::render("vignettes/data-loading.Rmd")
rmarkdown::render("vignettes/model-specification.Rmd")
rmarkdown::render("vignettes/feature-testing.Rmd")
rmarkdown::render("vignettes/feeding-4mo-trajectories.Rmd")
rmarkdown::render("vignettes/spline-fitting.Rmd")
rmarkdown::render("vignettes/validation.Rmd")
rmarkdown::render("vignettes/trajectory-analysis.Rmd")
rmarkdown::render("vignettes/gaussian-processes.Rmd")
```

## Installation

Install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("Hugemiler/LEVAiM")
```

LEVAiM depends on R packages used by the core spline workflow:

``` r
install.packages(c(
  "cli",
  "data.table",
  "mgcv",
  "testthat",
  "knitr",
  "rmarkdown"
))
```

Gaussian-process models require `brms` and a working Stan backend. For
example:

``` r
install.packages("brms")
install.packages(
  "cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
)

cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan(cores = 2)
```

Stan is only needed for `trajectory_control(engine = "gp")`. Spline
models do not require it.

## Basic Workflow

LEVAiM follows the same pipeline for spline and Gaussian-process models:

``` r
library(LEVAiM)

ds <- LongitudinalDataset(
  assays = list(taxa = taxa),
  metadata = metadata,
  time_col = "ageMonths",
  subject_col = "subject"
)

spec <- model_spec(
  response = assay_feature("taxa", "feature_name"),
  subject = subject_var("subject"),
  time = time_var("ageMonths"),
  covariates = list(group_effect("feeding_state"))
)

mf <- model_frame(
  ds,
  trajectory(spec, design = varying_trajectory("feeding_state")),
  control = trajectory_control(
    engine = "spline",
    na_action = "explicit",
    missing_level = "Unknown",
    sparse_min_n = 2
  )
)

fit <- fit_trajectory(mf)
trajectory_results(fit)
trajectory_effects(fit)
trajectory_derivatives(fit)
```

For assay-wide screening, fit the same trajectory model one feature at a
time and collate feature-level evidence:

``` r
feature_tests <- test_assay_features(
  ds,
  assay = "taxa",
  subject = "subject",
  time = "ageMonths",
  covariates = list(group_effect("feeding_state")),
  design = varying_trajectory("feeding_state"),
  control = trajectory_control(engine = "spline"),
  hypothesis = "trajectory_difference",
  test_covariate = "feeding_state",
  keep_fits = TRUE
)

feature_tests$table
```

The screen tests one declared null per feature and records the
hypothesis, tested levels, inference grid, maximum separation and its
time, p-value, q-value, and inferential limitation.
Trajectory-difference tests use a maximum standardized fitted-contrast
statistic calibrated from the model covariance; they do not select the
smallest p-value among group-specific smooth terms.

Those retained fits can then be compared as fitted trajectories. The
default comparison view produces feature-by-feature evidence within each
metadata level, plus distance matrices that can be used directly for
heatmaps.

``` r
trajectory_comparison <- compare_trajectories(
  feature_tests,
  metrics = c("fitted_correlation", "fitted_distance", "derivative_correlation")
)

trajectory_comparison
trajectory_distance_matrix(trajectory_comparison)
```

Functional moments summarize fitted level, AUC, temporal variation,
timing, trend, and velocity over declared biological windows. The same
long-table API uses coefficient-covariance draws for splines and
posterior fitted draws for GPs.

``` r
moments <- trajectory_moments(
  feature_tests,
  windows = data.frame(
    label = c("early", "feeding_window", "later"),
    start = c(0, 3, 5),
    end = c(3, 5, 12)
  ),
  draws = 500,
  seed = 1
)

moment_comparisons <- compare_trajectory_moments(
  moments,
  reference = "ExclFormulaFed"
)

plot_trajectory_moments(
  moments,
  moment = c("mean_level", "auc", "linear_trend"),
  window = "feeding_window"
)

plot_trajectory_moment_comparisons(
  moment_comparisons,
  moment = c("mean_level", "auc", "linear_trend"),
  window = "feeding_window"
)
```

Peak and valley fields carry explicit statuses such as
`interior_supported`, `boundary_only`, and `flat_or_uncertain`; boundary
maxima are not reported as biological peaks.

Correlation p-values in this table are descriptive summaries of fitted
grid curves; they are labeled as such because grid points are not
independent biological observations.

Grouped cross-validation keeps repeated measurements from the same
subject together and can retain held-out predictions for audit and
plotting:

``` r
cv <- cross_validate_trajectory(
  mf,
  v = 3,
  repeats = 2,
  seed = 1
)

cv_summary(cv)
cv_predictions(cv)
```

For Gaussian-process modeling, switch the engine and choose a starting
kernel:

``` r
gp_mf <- model_frame(
  ds,
  trajectory(spec),
  control = trajectory_control(
    engine = "gp",
    gp_kernel = "matern32", # also rbf, matern52, exponential
    gp_basis_k = NULL       # integer for Hilbert-space approximation
  )
)
```

GP posterior derivatives and functional moments use posterior fitted
draws. The full-cohort Backhed workflow in
`vignette("gaussian-processes", package = "LEVAiM")` includes an
explicit R-hat, effective-sample-size, and divergence gate before
interpretation. Both engines support Gaussian, binomial, and beta
responses under the explicit boundary policy in `RESPONSE-FAMILIES.md`.

## API Status

The spline mid-development API is documented in `API-LIFECYCLE.md`. In
short:

- spline modeling, prediction, effects, and grouped CV are the stable
  mid-development surface;
- MetaPhlAn/HUMAnN parsers are experimental while more input variants
  are added;
- Gaussian-process modeling is a `brms`-backed prototype pending broader
  fitted testing and review.

Requirement coverage is tracked in `REQUIREMENTS-ADHERENCE.md`.

Spline response families are specified explicitly. Binomial responses
must be exactly 0/1, while beta responses must lie strictly between 0
and 1. LEVAiM never adds pseudocounts or moves boundary values silently.
See `RESPONSE-FAMILIES.md` and the real-data
occurrence/positive-abundance workflow in
`vignette("response-families", package = "LEVAiM")`.

## Frequently Asked Questions

### I only want spline models. Do I need Stan?

No, Stan is only needed for `engine = "gp"`.

### How are missing annotations handled?

By default, missing categorical annotations are retained as an explicit
`"Unknown"` level and sparse categorical levels are collapsed to
`"Other"`. Rows are still removed when essential fields such as
response, time, subject, or continuous covariates are missing. Use
`na_action = "complete"` for strict complete-case annotation handling.

### Where should I start?

Start with:

``` r
vignette("data-loading", package = "LEVAiM")
```

Then move through the functional workflow vignettes:

``` r
vignette("model-specification", package = "LEVAiM")
vignette("feature-testing", package = "LEVAiM")
vignette("feeding-4mo-trajectories", package = "LEVAiM")
vignette("spline-fitting", package = "LEVAiM")
vignette("response-families", package = "LEVAiM")
vignette("gp-crash-course", package = "LEVAiM")
vignette("gaussian-processes", package = "LEVAiM")
vignette("validation", package = "LEVAiM")
vignette("trajectory-analysis", package = "LEVAiM")
```

### Why does package check mention long example filenames?

The included bioBakery-style example files preserve their original
descriptive names. R reports those paths as non-portable in source
tarballs, but the package tests and vignettes still run.

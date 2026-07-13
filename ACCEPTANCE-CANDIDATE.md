# LEVAiM Proof-of-Concept Acceptance Candidate

## Release Scope

This document records the technical acceptance candidate for the D0.3
longitudinal microbiome modeling proof of concept. Technical completion and
executive approval are reported separately.

The candidate supports the complete lifecycle:

`parse -> derive -> specify -> build -> compile -> fit -> validate -> predict -> interpret -> compare -> serialize`

## Checkpoint Position

| Checkpoint | Technical evidence | Remaining action |
| --- | --- | --- |
| D0.3.1 GP mid-development | 7/7 tested or documented | Executive review and signoff |
| D0.3.2 GP end-development | 6/7 complete, 1/7 partial | Decide and implement GP assay-wide batch inference if required |
| D0.3.3 spline mid-development | 7/7 tested or documented | Executive review and signoff |
| D0.3.4 spline later checkpoint | 8/8 tested and documented | Executive review and signoff |

## GP Acceptance Evidence

- Exact RBF, Matern-3/2, Matern-5/2, and exponential covariance paths compile.
- A fitted additional-kernel test exercises exponential covariance.
- `gp_basis_k` enables the `brms` Hilbert-space approximate GP path; exact GP
  remains the default.
- Gaussian, Bernoulli, and beta GP fits are tested.
- Single-response models accept multiple categorical and continuous covariates.
- Subject-grouped GP cross-validation holds out complete subjects.
- Prediction, effects, posterior derivatives, moments, plotting, diagnostics,
  and serialization are tested.
- The full aligned Backhed *B. longum* workflow passed its diagnostic gate with
  maximum R-hat 1.003, minimum bulk ESS 1692, minimum tail ESS 1880, and zero
  divergences in the acceptance run.
- The `gp-crash-course` vignette documents the statistical assumptions that
  connect kernels, likelihoods, irregular sampling, diagnostics, posterior
  derivatives, and biological claims.

## Spline Acceptance Evidence

- Natural, cubic, and GAM spline paths are tested.
- Automatic smoothness selection, multiple covariates, Gaussian/binomial/beta
  families, grouped CV, prediction, effects, derivatives, interpretation, and
  serialization are tested and documented.
- The Backhed proof of concept performs assay-wide feature testing, FDR
  correction, localized trajectory contrasts, functional moments, plotting,
  trajectory comparison, and distance-matrix generation.

## Frozen Contract

`API-LIFECYCLE.md` defines the frozen 0.1.x API. Stable signatures are guarded
by `tests/testthat/test-api-lifecycle.R`. Serialization schema 1 records object
class, package version, and creation time, rejects unsupported future schemas,
and reads legacy schema-0 raw objects.

## Explicit Limitation

`test_assay_features()` currently provides assay-wide batch inference for
spline models. GP follow-up is complete for one response at a time, but GP-scale
assay-wide fitting and multiplicity handling remain outside the frozen POC
contract. This limitation is reported as partial D0.3.2.6 evidence, not hidden
inside a completed claim.

## Approval Record

Technical evidence does not constitute executive approval. Approval status
remains open until the designated reviewer records acceptance against the
requirements matrix.

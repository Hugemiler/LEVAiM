# Response-Family and Boundary Policy

LEVAiM treats the response family as part of the model specification. The
package validates response domains before fitting and never repairs an
incompatible response silently.

## Gaussian

`family = "gaussian"` accepts finite numeric observations. The response may use
`transform = "identity"`, `"log1p"`, or `"log10p"`. The transformation and its
interpretive consequences must be reported with the model.

## Binomial

`family = "binomial"` requires every observed response to be exactly 0 or 1.
Predictions are event probabilities on the response scale. Grouped
cross-validation defaults to log loss and Brier score.

Counts, proportions with denominators, and fractional labels are not accepted
as binary outcomes. They require a separately specified sampling model.

## Beta

`family = "beta"` requires every observed response to satisfy `0 < y < 1`.
Exact zeros and ones are rejected before fitting. LEVAiM does not add
pseudocounts, apply continuity corrections, or squeeze values into the open
interval automatically.

The policy is deliberate:

- replacing zero with an arbitrary epsilon changes the estimand;
- microbiome zeros may represent absence, undersampling, or a detection limit;
- values at one may represent a genuine boundary event rather than an extreme
  draw from a continuous beta distribution;
- an implicit correction hides an analytical decision from the audit trail.

For zero-containing relative abundance, choose a model that matches the
question. A two-part analysis can model detection with a binomial trajectory and
positive abundance with a conditional beta trajectory. A joint analysis of the
boundary mass and continuous component requires an explicit hurdle or
zero/one-inflated family, which is not currently implemented in the LEVAiM
spline API.

## Transformations

Binomial and beta models require `transform = "identity"`. Response
transformations followed by a non-Gaussian likelihood are rejected. Missing
responses are removed by the model-frame preprocessing path; boundary values
are invalid data for beta regression and are not treated as missing.

## Current Backend Scope

Gaussian, binomial, and beta families are implemented for both spline and GP
engines. The GP backend maps binary outcomes to the Bernoulli likelihood and
uses `brms` beta regression for responses strictly inside `(0, 1)`. The same
zero/one boundary policy applies to both engines.

See `vignette("response-families", package = "LEVAiM")` for a real-data
occurrence and positive-abundance workflow.

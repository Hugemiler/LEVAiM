#####
# Trajectory modeling control
#####

#' Define trajectory modeling controls
#'
#' Stores engine, family, transformation, missing-data behavior, and
#' backend-specific fitting options.
#'
#' @param engine Modeling engine. One of `"spline"` or `"gp"`.
#' @param family Response family. Both engines support `"gaussian"`,
#'   `"binomial"`, and `"beta"`. Binomial responses must contain exactly 0 and 1. Beta
#'   responses must lie strictly inside `(0, 1)`; boundary values are rejected
#'   rather than transformed implicitly.
#' @param transform Response transformation. One of `"identity"`, `"log1p"`,
#'   or `"log10p"`.
#' @param na_action Missing-data handling. `"explicit"` keeps missing
#'   categorical annotations as `missing_level`; `"complete"` removes rows with
#'   missing model variables.
#' @param spline_k Basis dimension for spline smooths.
#' @param spline_basis Spline basis. `"gam"` uses `mgcv::s()`, `"natural"`
#'   uses `splines::ns()`, and `"cubic"` uses `splines::bs()`.
#' @param spline_method Fitting method for `mgcv::gam()`, usually `"REML"` or
#'   `"ML"`.
#' @param gp_kernel Gaussian-process covariance kernel. One of `"matern32"`,
#'   `"rbf"`, `"matern52"`, or `"exponential"`.
#' @param gp_basis_k Optional number of Hilbert-space basis functions. `NULL`
#'   fits an exact GP; an integer of at least 5 requests the scalable
#'   approximate GP implemented by `brms::gp(k = ...)`.
#' @param chains Number of MCMC chains for GP/brms backend.
#' @param iter Number of MCMC iterations for GP/brms backend.
#' @param cores Number of cores for GP/brms backend.
#' @param adapt_delta Target acceptance probability for GP/Stan sampling.
#' @param max_treedepth Maximum GP/Stan tree depth.
#' @param seed Optional random seed.
#' @param sparse_min_n Minimum count for a categorical annotation level. Levels
#'   with fewer samples are collapsed to `sparse_other_level`.
#' @param sparse_other_level Label used for collapsed sparse categorical levels.
#' @param missing_level Label used for missing categorical annotations when
#'   `na_action = "explicit"`.
#' @param drop_invariant_covariates Logical. Whether to remove non-trajectory
#'   covariates that have fewer than two levels after preprocessing.
#'
#' @return A `levaim_trajectory_control` object.
#' @export
trajectory_control <- function(
    engine = c("spline", "gp"),
    family = "gaussian",
    transform = c("identity", "log1p", "log10p"),
    na_action = c("explicit", "complete"),
    spline_k = 10L,
    spline_basis = c("gam", "natural", "cubic"),
    spline_method = c("REML", "ML", "GCV.Cp"),
    gp_kernel = c("matern32", "rbf", "matern52", "exponential"),
    gp_basis_k = NULL,
    chains = 4L,
    iter = 2000L,
    cores = 1L,
    adapt_delta = 0.95,
    max_treedepth = 12L,
    seed = NULL,
    sparse_min_n = 2L,
    sparse_other_level = "Other",
    missing_level = "Unknown",
    drop_invariant_covariates = TRUE
) {
  engine <- match.arg(engine)
  transform <- match.arg(transform)
  na_action <- match.arg(na_action)
  spline_basis <- match.arg(spline_basis)
  spline_method <- match.arg(spline_method)
  gp_kernel <- match.arg(gp_kernel)

  validate_trajectory_control(
    engine = engine,
    family = family,
    transform = transform,
    na_action = na_action,
    spline_k = spline_k,
    spline_basis = spline_basis,
    spline_method = spline_method,
    gp_kernel = gp_kernel,
    gp_basis_k = gp_basis_k,
    chains = chains,
    iter = iter,
    cores = cores,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    seed = seed,
    sparse_min_n = sparse_min_n,
    sparse_other_level = sparse_other_level,
    missing_level = missing_level,
    drop_invariant_covariates = drop_invariant_covariates
  )

  structure(
    list(
      engine = engine,
      family = family,
      transform = transform,
      na_action = na_action,
      spline = list(
        k = as.integer(spline_k),
        basis = spline_basis,
        method = spline_method
      ),
      gp = list(
        kernel = gp_kernel,
        basis_k = if (is.null(gp_basis_k)) NULL else as.integer(gp_basis_k),
        approximation = if (is.null(gp_basis_k)) "exact" else "hsgp",
        chains = as.integer(chains),
        iter = as.integer(iter),
        cores = as.integer(cores),
        adapt_delta = adapt_delta,
        max_treedepth = as.integer(max_treedepth)
      ),
      seed = seed,
      annotations = list(
        sparse_min_n = as.integer(sparse_min_n),
        sparse_other_level = sparse_other_level,
        missing_level = missing_level,
        drop_invariant_covariates = drop_invariant_covariates
      )
    ),
    class = "levaim_trajectory_control"
  )
}


#' Validate trajectory control arguments
#'
#' @keywords internal
validate_trajectory_control <- function(
    engine,
    family,
    transform,
    na_action,
    spline_k,
    spline_basis,
    spline_method,
    gp_kernel,
    gp_basis_k,
    chains,
    iter,
    cores,
    adapt_delta,
    max_treedepth,
    seed,
    sparse_min_n,
    sparse_other_level,
    missing_level,
    drop_invariant_covariates
) {
  if (!is.character(family) || length(family) != 1L ||
      !family %in% c("gaussian", "binomial", "beta")) {
    cli::cli_abort("`family` must be one of `gaussian`, `binomial`, or `beta`.")
  }

  if (!is.null(gp_basis_k) &&
      (!is.numeric(gp_basis_k) || length(gp_basis_k) != 1L ||
       is.na(gp_basis_k) || gp_basis_k < 5L)) {
    cli::cli_abort("`gp_basis_k` must be `NULL` or a numeric scalar >= 5.")
  }

  if (family != "gaussian" && transform != "identity") {
    cli::cli_abort("Non-Gaussian families currently require `transform = 'identity'`.")
  }

  if (!is.numeric(spline_k) || length(spline_k) != 1L || spline_k < 3L) {
    cli::cli_abort("`spline_k` must be a numeric scalar >= 3.")
  }

  if (!is.character(spline_basis) || length(spline_basis) != 1L) {
    cli::cli_abort("`spline_basis` must be a character scalar.")
  }

  if (!is.numeric(chains) || length(chains) != 1L || chains < 1L) {
    cli::cli_abort("`chains` must be a positive numeric scalar.")
  }

  if (!is.numeric(iter) || length(iter) != 1L || iter < 1L) {
    cli::cli_abort("`iter` must be a positive numeric scalar.")
  }

  if (!is.numeric(cores) || length(cores) != 1L || cores < 1L) {
    cli::cli_abort("`cores` must be a positive numeric scalar.")
  }

  if (!is.numeric(adapt_delta) || length(adapt_delta) != 1L ||
      adapt_delta <= 0 || adapt_delta >= 1) {
    cli::cli_abort("`adapt_delta` must lie strictly between 0 and 1.")
  }

  if (!is.numeric(max_treedepth) || length(max_treedepth) != 1L ||
      max_treedepth < 1L) {
    cli::cli_abort("`max_treedepth` must be a positive numeric scalar.")
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    cli::cli_abort("`seed` must be `NULL` or a numeric scalar.")
  }

  if (!is.numeric(sparse_min_n) || length(sparse_min_n) != 1L || sparse_min_n < 1L) {
    cli::cli_abort("`sparse_min_n` must be a numeric scalar >= 1.")
  }

  if (!is.character(sparse_other_level) || length(sparse_other_level) != 1L || !nzchar(sparse_other_level)) {
    cli::cli_abort("`sparse_other_level` must be a non-empty character scalar.")
  }

  if (!is.character(missing_level) || length(missing_level) != 1L || !nzchar(missing_level)) {
    cli::cli_abort("`missing_level` must be a non-empty character scalar.")
  }

  if (!is.logical(drop_invariant_covariates) || length(drop_invariant_covariates) != 1L || is.na(drop_invariant_covariates)) {
    cli::cli_abort("`drop_invariant_covariates` must be `TRUE` or `FALSE`.")
  }

  invisible(TRUE)
}


#' Print trajectory modeling controls
#'
#' @param x A `levaim_trajectory_control` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.levaim_trajectory_control <- function(x, ...) {
  cli::cli_h1("LEVAiM trajectory control")

  cli::cli_text("Engine: {.field {x$engine}}")
  cli::cli_text("Family: {.field {x$family}}")
  cli::cli_text("Transform: {.field {x$transform}}")
  cli::cli_text("NA action: {.field {x$na_action}}")

  if (x$engine == "spline") {
    cli::cli_h2("Spline controls")
    cli::cli_text("Basis dimension k: {x$spline$k}")
    cli::cli_text("Basis: {.field {x$spline$basis}}")
    cli::cli_text("Method: {.field {x$spline$method}}")
  }

  if (x$engine == "gp") {
    cli::cli_h2("GP controls")
    cli::cli_text("Kernel: {.field {x$gp$kernel}}")
    cli::cli_text("Approximation: {.field {x$gp$approximation}}")
    if (!is.null(x$gp$basis_k)) {
      cli::cli_text("Hilbert-space basis functions: {x$gp$basis_k}")
    }
    cli::cli_text("Chains: {x$gp$chains}")
    cli::cli_text("Iterations: {x$gp$iter}")
    cli::cli_text("Cores: {x$gp$cores}")
    cli::cli_text("Adapt delta: {x$gp$adapt_delta}")
    cli::cli_text("Maximum tree depth: {x$gp$max_treedepth}")
  }

  if (!is.null(x$seed)) {
    cli::cli_text("Seed: {x$seed}")
  }

  cli::cli_h2("Annotation preprocessing")
  cli::cli_text("Sparse level minimum n: {x$annotations$sparse_min_n}")
  cli::cli_text("Sparse level label: {.field {x$annotations$sparse_other_level}}")
  cli::cli_text("Missing level label: {.field {x$annotations$missing_level}}")
  cli::cli_text("Drop invariant covariates: {x$annotations$drop_invariant_covariates}")

  invisible(x)
}

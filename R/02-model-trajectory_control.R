#####
# Trajectory modeling control
#####

#' Define trajectory modeling controls
#'
#' Stores engine, family, transformation, missing-data behavior, and
#' backend-specific fitting options.
#'
#' @param engine Modeling engine. One of `"spline"` or `"gp"`.
#' @param family Model family. For now, usually `"gaussian"`.
#' @param transform Response transformation. One of `"identity"`, `"log1p"`,
#'   or `"log10p"`.
#' @param na_action Missing-data handling. Currently `"complete"`.
#' @param spline_k Basis dimension for spline smooths.
#' @param spline_method Fitting method for `mgcv::gam()`, usually `"REML"` or
#'   `"ML"`.
#' @param gp_kernel Gaussian process kernel. Placeholder for GP backend.
#' @param chains Number of MCMC chains for GP/brms backend.
#' @param iter Number of MCMC iterations for GP/brms backend.
#' @param cores Number of cores for GP/brms backend.
#' @param seed Optional random seed.
#'
#' @return A `levaim_trajectory_control` object.
#' @export
trajectory_control <- function(
    engine = c("spline", "gp"),
    family = "gaussian",
    transform = c("identity", "log1p", "log10p"),
    na_action = c("complete"),
    spline_k = 10L,
    spline_method = c("REML", "ML", "GCV.Cp"),
    gp_kernel = c("matern32", "rbf"),
    chains = 4L,
    iter = 2000L,
    cores = 1L,
    seed = NULL
) {
  engine <- match.arg(engine)
  transform <- match.arg(transform)
  na_action <- match.arg(na_action)
  spline_method <- match.arg(spline_method)
  gp_kernel <- match.arg(gp_kernel)

  validate_trajectory_control(
    engine = engine,
    family = family,
    transform = transform,
    na_action = na_action,
    spline_k = spline_k,
    spline_method = spline_method,
    gp_kernel = gp_kernel,
    chains = chains,
    iter = iter,
    cores = cores,
    seed = seed
  )

  structure(
    list(
      engine = engine,
      family = family,
      transform = transform,
      na_action = na_action,
      spline = list(
        k = as.integer(spline_k),
        method = spline_method
      ),
      gp = list(
        kernel = gp_kernel,
        chains = as.integer(chains),
        iter = as.integer(iter),
        cores = as.integer(cores)
      ),
      seed = seed
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
    spline_method,
    gp_kernel,
    chains,
    iter,
    cores,
    seed
) {
  if (!is.character(family) || length(family) != 1L) {
    cli::cli_abort("`family` must be a character scalar.")
  }

  if (!is.numeric(spline_k) || length(spline_k) != 1L || spline_k < 3L) {
    cli::cli_abort("`spline_k` must be a numeric scalar >= 3.")
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

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L)) {
    cli::cli_abort("`seed` must be `NULL` or a numeric scalar.")
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
    cli::cli_text("Method: {.field {x$spline$method}}")
  }

  if (x$engine == "gp") {
    cli::cli_h2("GP controls")
    cli::cli_text("Kernel: {.field {x$gp$kernel}}")
    cli::cli_text("Chains: {x$gp$chains}")
    cli::cli_text("Iterations: {x$gp$iter}")
    cli::cli_text("Cores: {x$gp$cores}")
  }

  if (!is.null(x$seed)) {
    cli::cli_text("Seed: {x$seed}")
  }

  invisible(x)
}

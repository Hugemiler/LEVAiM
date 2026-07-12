#####
# Model fitting
#####

#' Fit a LEVAiM trajectory
#'
#' Fits a trajectory model using the backend specified by the associated
#' trajectory control object.
#'
#' @param mf A `levaim_model_frame`.
#'
#' @return A `levaim_fit`.
#' @importFrom stats formula predict
#' @export
fit_trajectory <- function(mf) {

  if (!inherits(mf, "levaim_model_frame")) {
    cli::cli_abort(
      "`mf` must be created with `model_frame()`."
    )
  }

  switch(
    mf$control$engine,

    spline = fit_trajectory_spline(mf),

    gp = fit_trajectory_gp(mf)
  )
}

#' Fit a spline trajectory
#'
#' @param mf A `levaim_model_frame`.
#'
#' @return A `levaim_fit`.
#'
#' @keywords internal
fit_trajectory_spline <- function(mf) {

  fit <- mgcv::gam(

    formula = mf$formula,

    data = mf$data,

    family = mf$control$family,

    method = mf$control$spline$method,

    select = TRUE

  )

  structure(

    list(

      fit = fit,

      model_frame = mf,

      engine = "spline"

    ),

    class = "levaim_fit"

  )

}

#' @keywords internal
fit_trajectory_gp <- function(mf) {
  if (!requireNamespace("brms", quietly = TRUE)) {
    cli::cli_abort(c(
      "The GP backend requires the {.pkg brms} package.",
      "i" = "Install {.pkg brms} and a Stan backend, then retry with `engine = 'gp'`."
    ))
  }

  if (mf$control$family != "gaussian") {
    cli::cli_abort("The brms GP backend currently supports `family = 'gaussian'`.")
  }

  fit <- brms::brm(
    formula = mf$formula,
    data = mf$data,
    family = "gaussian",
    chains = mf$control$gp$chains,
    iter = mf$control$gp$iter,
    cores = mf$control$gp$cores,
    seed = mf$control$seed,
    refresh = 0
  )

  structure(
    list(
      fit = fit,
      model_frame = mf,
      engine = "gp"
    ),
    class = "levaim_fit"
  )
}


#' @keywords internal
brms_gp_covariance <- function(kernel) {
  switch(
    kernel,
    rbf = "exp_quad",
    matern32 = "matern32",
    cli::cli_abort("Unknown GP kernel: {.val {kernel}}.")
  )
}

#' @export
print.levaim_fit <- function(x, ...) {

  cli::cli_h1("LEVAiM fit")

  cli::cli_text(
    "Engine: {.field {x$engine}}"
  )

  cli::cli_text(
    "Backend: {.field {class(x$fit)[1]}}"
  )

  cli::cli_text(
    "Formula:"
  )

  if (x$engine == "gp") {
    print(x$model_frame$formula)
  } else {
    print(formula(x$fit))
  }

  invisible(x)

}

#' @export
summary.levaim_fit <- function(object, ...) {

  if (object$engine == "spline") {

    return(
      summary(object$fit)
    )

  }

  if (object$engine == "gp") {
    return(summary(object$fit))
  }

  NextMethod()

}

#' @export
predict.levaim_fit <- function(
    object,
    newdata = NULL,
    ...
) {

  if (object$engine == "gp") {
    fitted <- stats::fitted(object$fit, newdata = newdata, ...)
    return(as.numeric(fitted[, "Estimate"]))
  }

  predict(object$fit, newdata = newdata, ...)

}

#' @export
plot.levaim_fit <- function(
    x,
    ...
) {

  if (x$engine == "gp") {
    return(plot_trajectory_effects(x, ...))
  }

  plot(x$fit, pages = 1, ...)

}

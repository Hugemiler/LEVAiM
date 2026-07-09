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

  cli::cli_abort(
    "GP backend has not yet been implemented."
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

  print(formula(x$fit))

  invisible(x)

}

#' @export
summary.levaim_fit <- function(object, ...) {

  if (object$engine == "spline") {

    return(
      summary(object$fit)
    )

  }

  NextMethod()

}

#' @export
predict.levaim_fit <- function(
    object,
    newdata = NULL,
    ...
) {

  predict(
    object$fit,
    newdata = newdata,
    ...
  )

}

#' @export
plot.levaim_fit <- function(
    x,
    ...
) {

  plot(
    x$fit,
    pages = 1,
    ...
  )

}

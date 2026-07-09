#####
# LEVAiM analytical results
#####

#' Extract analytical results from a LEVAiM fit
#'
#' Produces a lightweight LEVAiM-native summary from a fitted trajectory model.
#'
#' @param fit A `levaim_fit`.
#'
#' @return A `levaim_results` object.
#' @export
trajectory_results <- function(fit) {
  if (!inherits(fit, "levaim_fit")) {
    cli::cli_abort("`fit` must be created with `fit_trajectory()`.")
  }

  switch(
    fit$engine,
    spline = trajectory_results_spline(fit),
    gp = cli::cli_abort("GP results are not implemented yet.")
  )
}


#' Extract spline/GAM results
#'
#' @keywords internal
trajectory_results_spline <- function(fit) {
  gam_summary <- summary(fit$fit)

  mf <- fit$model_frame
  traj <- mf$trajectory
  spec <- traj$spec
  design <- traj$design

  response_label <- format_response_label(spec$response)

  smooth_table <- as.data.frame(gam_summary$s.table)
  parametric_table <- as.data.frame(gam_summary$p.table)

  structure(
    list(
      engine = fit$engine,
      response = response_label,
      formula = mf$formula,
      trajectory_type = design$type,
      trajectory_by = if (design$type == "varying") design$by else NULL,
      family = mf$control$family,
      transform = mf$control$transform,
      n = nrow(mf$data),
      r_sq = gam_summary$r.sq,
      deviance_explained = gam_summary$dev.expl,
      parametric_terms = parametric_table,
      smooth_terms = smooth_table,
      backend_summary = gam_summary
    ),
    class = "levaim_results"
  )
}


#' Format response label
#'
#' @keywords internal
format_response_label <- function(response) {
  if (inherits(response, "levaim_assay_feature")) {
    return(paste0(response$assay, " / ", response$feature))
  }

  if (inherits(response, "levaim_metadata_feature")) {
    return(paste0("metadata / ", response$column))
  }

  "unknown response"
}


#' Print LEVAiM results
#'
#' @param x A `levaim_results` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.levaim_results <- function(x, ...) {
  cli::cli_h1("LEVAiM trajectory results")

  cli::cli_text("Response: {.field {x$response}}")
  cli::cli_text("Engine: {.field {x$engine}}")
  cli::cli_text("Family: {.field {x$family}}")
  cli::cli_text("Transform: {.field {x$transform}}")
  cli::cli_text("Samples: {x$n}")

  cli::cli_text("Formula:")
  cli::cli_text("{deparse(x$formula)}")

  cli::cli_h2("Model fit")
  cli::cli_text("R-squared: {round(x$r_sq, 4)}")
  cli::cli_text("Deviance explained: {round(100 * x$deviance_explained, 2)}%")

  cli::cli_h2("Trajectory design")
  cli::cli_text("Type: {.field {x$trajectory_type}}")

  if (!is.null(x$trajectory_by)) {
    cli::cli_text("Varying by: {.field {paste(x$trajectory_by, collapse = ' × ')}}")
  }

  cli::cli_h2("Smooth terms")
  print(x$smooth_terms)

  invisible(x)
}

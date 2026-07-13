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
    gp = trajectory_results_gp(fit)
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


#' Extract Gaussian process results
#'
#' @keywords internal
trajectory_results_gp <- function(fit) {
  mf <- fit$model_frame
  traj <- mf$trajectory
  spec <- traj$spec
  design <- traj$design

  response_label <- format_response_label(spec$response)
  fitted <- as.numeric(stats::fitted(fit$fit, newdata = mf$data)[, "Estimate"])
  observed <- mf$data$.y
  residual <- observed - fitted
  total <- observed - mean(observed)
  r_sq <- 1 - sum(residual^2) / sum(total^2)

  kernel_table <- data.frame(
    kernel = mf$control$gp$kernel,
    brms_covariance = brms_gp_covariance(mf$control$gp$kernel),
    approximation = mf$control$gp$approximation,
    basis_k = if (is.null(mf$control$gp$basis_k)) NA_integer_ else mf$control$gp$basis_k
  )
  backend_summary <- summary(fit$fit)
  diagnostics <- gp_fit_diagnostics(fit$fit, backend_summary)
  fixed <- backend_summary$fixed
  parametric_table <- if (is.null(fixed) || nrow(fixed) == 0L) {
    data.frame()
  } else {
    data.frame(
      term = rownames(fixed),
      estimate = as.numeric(fixed[, "Estimate"]),
      rhat = if ("Rhat" %in% colnames(fixed)) as.numeric(fixed[, "Rhat"]) else NA_real_,
      bulk_ess = if ("Bulk_ESS" %in% colnames(fixed)) as.numeric(fixed[, "Bulk_ESS"]) else NA_real_,
      tail_ess = if ("Tail_ESS" %in% colnames(fixed)) as.numeric(fixed[, "Tail_ESS"]) else NA_real_,
      row.names = NULL
    )
  }

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
      r_sq = r_sq,
      deviance_explained = r_sq,
      parametric_terms = parametric_table,
      smooth_terms = kernel_table,
      diagnostics = diagnostics,
      backend_summary = backend_summary
    ),
    class = "levaim_results"
  )
}


#' @keywords internal
gp_fit_diagnostics <- function(brms_fit, backend_summary = summary(brms_fit)) {
  matrices <- list(backend_summary$fixed, backend_summary$spec_pars)
  matrices <- matrices[vapply(
    matrices,
    function(x) is.matrix(x) || is.data.frame(x),
    logical(1)
  )]

  rhat <- unlist(lapply(matrices, function(x) {
    if ("Rhat" %in% colnames(x)) x[, "Rhat"] else numeric()
  }))
  bulk_ess <- unlist(lapply(matrices, function(x) {
    if ("Bulk_ESS" %in% colnames(x)) x[, "Bulk_ESS"] else numeric()
  }))
  tail_ess <- unlist(lapply(matrices, function(x) {
    if ("Tail_ESS" %in% colnames(x)) x[, "Tail_ESS"] else numeric()
  }))

  sampler <- tryCatch(brms::nuts_params(brms_fit), error = function(e) NULL)
  divergences <- if (is.null(sampler)) {
    NA_integer_
  } else {
    as.integer(sum(
      sampler$Parameter == "divergent__" & sampler$Value > 0,
      na.rm = TRUE
    ))
  }

  max_rhat <- diagnostic_max(rhat)
  min_bulk_ess <- diagnostic_min(bulk_ess)
  min_tail_ess <- diagnostic_min(tail_ess)
  rhat_ok <- is.finite(max_rhat) && max_rhat <= 1.01
  ess_ok <- is.finite(min_bulk_ess) && is.finite(min_tail_ess) &&
    min_bulk_ess >= 400 && min_tail_ess >= 400
  divergence_ok <- identical(divergences, 0L)

  data.frame(
    max_rhat = max_rhat,
    min_bulk_ess = min_bulk_ess,
    min_tail_ess = min_tail_ess,
    divergences = divergences,
    rhat_ok = rhat_ok,
    ess_ok = ess_ok,
    divergence_ok = divergence_ok,
    convergence_ok = rhat_ok && ess_ok && divergence_ok,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
diagnostic_max <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}


#' @keywords internal
diagnostic_min <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
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
  cat(deparse(x$formula), sep = "\n")

  cli::cli_h2("Model fit")
  cli::cli_text("R-squared: {round(x$r_sq, 4)}")
  cli::cli_text("Deviance explained: {round(100 * x$deviance_explained, 2)}%")

  cli::cli_h2("Trajectory design")
  cli::cli_text("Type: {.field {x$trajectory_type}}")

  if (!is.null(x$trajectory_by)) {
    cli::cli_text("Varying by: {.field {paste(x$trajectory_by, collapse = ' x ')}}")
  }

  cli::cli_h2("Smooth terms")
  print(x$smooth_terms)

  invisible(x)
}

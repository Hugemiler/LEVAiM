#####
# Formula compilation
#####

#' Compile a trajectory object into a model formula
#'
#' @param traj A `levaim_trajectory` object.
#' @param engine Modeling engine. One of `"spline"` or `"gp"`.
#' @param spline_k Optional basis dimension for spline smooths.
#'
#' @return A formula.
#' @keywords internal
compile_formula <- function(
    traj,
    engine = c("spline", "gp"),
    spline_k = NULL
) {
  engine <- match.arg(engine)

  if (!inherits(traj, "levaim_trajectory")) {
    cli::cli_abort("`traj` must be a `levaim_trajectory` object.")
  }

  switch(
    engine,
    spline = compile_formula_spline(traj, k = spline_k),
    gp = compile_formula_gp(traj)
  )
}


#' Compile a spline/GAM formula
#'
#' Compiles a LEVAiM trajectory into an `mgcv::gam()`-style formula.
#'
#' @param traj A `levaim_trajectory` object.
#' @param k Optional basis dimension for spline smooths.
#'
#' @return A formula.
#' @keywords internal
compile_formula_spline <- function(traj, k = NULL) {
  spec <- traj$spec
  design <- traj$design

  rhs_terms <- c(
    compile_covariate_terms_spline(spec$covariates),
    compile_time_terms_spline(spec$time$column, design, k = k),
    paste0("s(", spec$subject$column, ", bs = 're')")
  )

  rhs_terms <- unique(rhs_terms[nzchar(rhs_terms)])

  stats::as.formula(
    paste(".y ~", paste(rhs_terms, collapse = " + "))
  )
}


#' Compile a GP formula
#'
#' Compiles a LEVAiM trajectory into a `brms::brm()`-style formula.
#'
#' @param traj A `levaim_trajectory` object.
#'
#' @return A formula.
#' @keywords internal
compile_formula_gp <- function(traj) {
  spec <- traj$spec
  design <- traj$design

  rhs_terms <- c(
    compile_covariate_terms_gp(spec$covariates),
    compile_time_terms_gp(spec$time$column, design),
    paste0("(1 | ", spec$subject$column, ")")
  )

  rhs_terms <- unique(rhs_terms[nzchar(rhs_terms)])

  stats::as.formula(
    paste(".y ~", paste(rhs_terms, collapse = " + "))
  )
}


#' Compile covariate terms for spline/GAM models
#'
#' @keywords internal
compile_covariate_terms_spline <- function(covariates) {
  if (length(covariates) == 0L) {
    return(character())
  }

  vapply(
    covariates,
    function(cov) {
      if (inherits(cov, "levaim_random_effect")) {
        return(paste0("s(", cov$column, ", bs = 're')"))
      }

      if (inherits(cov, "levaim_group_effect")) {
        return(cov$column)
      }

      if (inherits(cov, "levaim_continuous_effect")) {
        return(cov$column)
      }

      cli::cli_abort("Unknown covariate type for column {.field {cov$column}}.")
    },
    character(1)
  )
}


#' Compile covariate terms for GP/brms models
#'
#' @keywords internal
compile_covariate_terms_gp <- function(covariates) {
  if (length(covariates) == 0L) {
    return(character())
  }

  vapply(
    covariates,
    function(cov) {
      if (inherits(cov, "levaim_random_effect")) {
        return(paste0("(1 | ", cov$column, ")"))
      }

      if (inherits(cov, "levaim_group_effect")) {
        return(cov$column)
      }

      if (inherits(cov, "levaim_continuous_effect")) {
        return(cov$column)
      }

      cli::cli_abort("Unknown covariate type for column {.field {cov$column}}.")
    },
    character(1)
  )
}


#' Compile time trajectory terms for spline/GAM models
#'
#' @keywords internal
compile_time_terms_spline <- function(time_col, design, k = NULL) {
  k_arg <- if (is.null(k)) "" else paste0(", k = ", as.integer(k))

  if (design$type %in% c("shared", "naive")) {
    return(paste0("s(", time_col, k_arg, ")"))
  }

  if (design$type == "varying") {
    by_expr <- trajectory_by_expr(design)

    return(paste0("s(", time_col, ", by = ", by_expr, k_arg, ")"))
  }

  cli::cli_abort("Unknown trajectory design type: {.val {design$type}}.")
}


#' Compile time trajectory terms for GP/brms models
#'
#' @keywords internal
compile_time_terms_gp <- function(time_col, design) {
  if (design$type %in% c("shared", "naive")) {
    return(paste0("gp(", time_col, ")"))
  }

  if (design$type == "varying") {
    by_expr <- trajectory_by_expr(design)

    return(paste0("gp(", time_col, ", by = ", by_expr, ")"))
  }

  cli::cli_abort("Unknown trajectory design type: {.val {design$type}}.")
}


#' Build a symbolic trajectory-by expression
#'
#' For one variable, returns that variable name. For interaction trajectories,
#' returns an inline symbolic interaction expression instead of materializing a
#' helper column.
#'
#' @keywords internal
trajectory_by_expr <- function(design) {
  if (design$type != "varying") {
    return(NULL)
  }

  if (length(design$by) == 1L && !isTRUE(design$interaction)) {
    return(design$by[[1]])
  }

  paste0("interaction(", paste(design$by, collapse = ", "), ")")
}

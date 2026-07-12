#####
# Formula compilation
#####

#' Compile a trajectory object into a model formula
#'
#' @param traj A `levaim_trajectory` object.
#' @param engine Modeling engine. One of `"spline"` or `"gp"`.
#' @param spline_k Optional basis dimension for spline smooths.
#' @param spline_basis Spline basis for spline models.
#' @param gp_kernel GP kernel for GP models.
#'
#' @return A formula.
#' @importFrom splines ns bs
#' @keywords internal
compile_formula <- function(
    traj,
    engine = c("spline", "gp"),
    spline_k = NULL,
    spline_basis = "gam",
    gp_kernel = "matern32"
) {
  engine <- match.arg(engine)

  if (!inherits(traj, "levaim_trajectory")) {
    cli::cli_abort("`traj` must be a `levaim_trajectory` object.")
  }

  switch(
    engine,
    spline = compile_formula_spline(traj, k = spline_k, basis = spline_basis),
    gp = compile_formula_gp(traj, kernel = gp_kernel)
  )
}


#' Compile a spline/GAM formula
#'
#' Compiles a LEVAiM trajectory into an `mgcv::gam()`-style formula.
#'
#' @param traj A `levaim_trajectory` object.
#' @param k Optional basis dimension for spline smooths.
#' @param basis Spline basis. One of `"gam"`, `"natural"`, or `"cubic"`.
#'
#' @return A formula.
#' @keywords internal
compile_formula_spline <- function(traj, k = NULL, basis = "gam") {
  basis <- match.arg(basis, c("gam", "natural", "cubic"))
  spec <- traj$spec
  design <- traj$design

  rhs_terms <- c(
    compile_covariate_terms_spline(spec$covariates),
    compile_time_terms_spline(spec$time$column, design, k = k, basis = basis),
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
#' @param kernel GP kernel. One of `"rbf"` or `"matern32"`.
#'
#' @return A formula.
#' @keywords internal
compile_formula_gp <- function(traj, kernel = "matern32") {
  spec <- traj$spec
  design <- traj$design

  rhs_terms <- c(
    compile_covariate_terms_gp(spec$covariates),
    compile_time_terms_gp(spec$time$column, design, kernel = kernel),
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
compile_time_terms_spline <- function(time_col, design, k = NULL, basis = "gam") {
  basis <- match.arg(basis, c("gam", "natural", "cubic"))
  k_arg <- if (is.null(k)) "" else paste0(", k = ", as.integer(k))

  if (basis == "natural") {
    df_arg <- if (is.null(k)) "" else paste0(", df = ", as.integer(k))
    time_term <- paste0("splines::ns(", time_col, df_arg, ")")
  } else if (basis == "cubic") {
    df_arg <- if (is.null(k)) "" else paste0(", df = ", as.integer(k))
    time_term <- paste0("splines::bs(", time_col, df_arg, ", degree = 3)")
  } else {
    time_term <- paste0("s(", time_col, k_arg, ")")
  }

  if (design$type %in% c("shared", "naive")) {
    return(time_term)
  }

  if (design$type == "varying") {
    by_expr <- trajectory_by_expr(design)

    if (basis == "gam") {
      return(paste0("s(", time_col, ", by = ", by_expr, k_arg, ")"))
    }

    return(paste0(by_expr, ":", time_term))
  }

  cli::cli_abort("Unknown trajectory design type: {.val {design$type}}.")
}


#' Compile time trajectory terms for GP/brms models
#'
#' @keywords internal
compile_time_terms_gp <- function(time_col, design, kernel = "matern32") {
  cov <- brms_gp_covariance(kernel)

  if (design$type %in% c("shared", "naive")) {
    return(paste0("gp(", time_col, ", cov = '", cov, "')"))
  }

  if (design$type == "varying") {
    by_expr <- trajectory_by_expr(design)

    return(paste0("gp(", time_col, ", by = ", by_expr, ", cov = '", cov, "')"))
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

#####
# Publication-oriented trajectory plots
#####

#' Plot fitted trajectories with observations and uncertainty
#'
#' @param fit A fitted `levaim_fit`.
#' @param n Number of fitted grid points.
#' @param time_values Optional numeric fitted time grid.
#' @param levels Optional trajectory-group levels to display.
#' @param window Optional length-two numeric vector highlighted on the time axis.
#' @param points Logical. Whether to overlay model-frame observations.
#' @param point_alpha Observation opacity.
#' @param palette Optional named or unnamed discrete color vector.
#'
#' @return A `ggplot` object.
#' @importFrom rlang .data
#' @export
plot_trajectory_report <- function(
    fit,
    n = 200L,
    time_values = NULL,
    levels = NULL,
    window = NULL,
    points = TRUE,
    point_alpha = 0.28,
    palette = NULL
) {
  if (!inherits(fit, "levaim_fit")) {
    cli::cli_abort("`fit` must be created with `fit_trajectory()`.")
  }
  effects <- trajectory_effects(fit, n = n, time_values = time_values)
  mf <- fit$model_frame
  time_col <- mf$trajectory$spec$time$column
  by_cols <- if (mf$trajectory$design$type == "varying") {
    mf$trajectory$design$by
  } else {
    character()
  }
  group_col <- ".trajectory_group"
  effects[[group_col]] <- if (length(by_cols) == 0L) {
    "shared"
  } else {
    as.character(interaction(effects[by_cols], drop = TRUE, sep = ":"))
  }
  observations <- mf$data
  observations[[group_col]] <- if (length(by_cols) == 0L) {
    "shared"
  } else {
    as.character(interaction(observations[by_cols], drop = TRUE, sep = ":"))
  }

  if (!is.null(levels)) {
    effects <- effects[effects[[group_col]] %in% levels, , drop = FALSE]
    observations <- observations[observations[[group_col]] %in% levels, , drop = FALSE]
  }
  if (nrow(effects) == 0L) {
    cli::cli_abort("No fitted trajectory levels remain to plot.")
  }

  colors <- trajectory_plot_palette(unique(effects[[group_col]]), palette)
  response_label <- short_response_label(mf$trajectory$spec$response)
  plot <- ggplot2::ggplot()
  if (!is.null(window)) {
    if (!is.numeric(window) || length(window) != 2L || window[[1L]] >= window[[2L]]) {
      cli::cli_abort("`window` must be a numeric vector with `start < end`.")
    }
    plot <- plot + ggplot2::annotate(
      "rect",
      xmin = window[[1L]], xmax = window[[2L]],
      ymin = -Inf, ymax = Inf,
      fill = "#D9E7EE", alpha = 0.45
    )
  }
  if (isTRUE(points)) {
    plot <- plot + ggplot2::geom_point(
      data = observations,
      ggplot2::aes(
        x = .data[[time_col]],
        y = .data$.y,
        color = .data[[group_col]]
      ),
      alpha = point_alpha,
      size = 1.25
    )
  }
  plot +
    ggplot2::geom_ribbon(
      data = effects,
      ggplot2::aes(
        x = .data[[time_col]],
        ymin = .data$.lower,
        ymax = .data$.upper,
        fill = .data[[group_col]],
        group = .data[[group_col]]
      ),
      alpha = 0.16,
      color = NA
    ) +
    ggplot2::geom_line(
      data = effects,
      ggplot2::aes(
        x = .data[[time_col]],
        y = .data$.fitted,
        color = .data[[group_col]],
        group = .data[[group_col]]
      ),
      linewidth = 0.9
    ) +
    ggplot2::scale_color_manual(values = colors, name = plot_group_label(by_cols)) +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::labs(
      title = response_label,
      subtitle = "Fitted longitudinal trajectories with pointwise uncertainty",
      x = time_col,
      y = paste0("Modeled response (", mf$control$transform, ")")
    ) +
    trajectory_plot_theme()
}


#' Plot pairwise fitted trajectory contrasts over time
#'
#' @param contrasts Output from [trajectory_contrasts()].
#' @param show_points Logical. Mark grid points whose pointwise interval excludes
#'   zero.
#' @param palette Optional discrete color vector for contrast pairs.
#'
#' @return A `ggplot` object.
#' @export
plot_trajectory_contrasts <- function(
    contrasts,
    show_points = TRUE,
    palette = NULL
) {
  required <- c("group_1", "group_2", "estimate", "lower", "upper", "p_value")
  if (!is.data.frame(contrasts) || !all(required %in% names(contrasts))) {
    cli::cli_abort("`contrasts` must be created with `trajectory_contrasts()`.")
  }
  time_col <- attr(contrasts, "time_variable")
  if (is.null(time_col) || !time_col %in% names(contrasts)) {
    candidates <- setdiff(names(contrasts), c(required, "std_error", "z_value"))
    numeric_candidates <- candidates[vapply(contrasts[candidates], is.numeric, logical(1))]
    time_col <- numeric_candidates[[1L]]
  }
  contrasts$.contrast <- paste(contrasts$group_1, contrasts$group_2, sep = " - ")
  contrasts$.separated <- contrasts$lower > 0 | contrasts$upper < 0
  colors <- trajectory_plot_palette(unique(contrasts$.contrast), palette)

  plot <- ggplot2::ggplot(
    contrasts,
    ggplot2::aes(
      x = .data[[time_col]],
      y = .data$estimate,
      color = .data$.contrast,
      fill = .data$.contrast
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "#555555", linewidth = 0.45) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      alpha = 0.16,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(stats::as.formula("~ .contrast"), scales = "free_y") +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::scale_fill_manual(values = colors, guide = "none") +
    ggplot2::labs(
      title = "Trajectory differences over time",
      subtitle = "Intervals are pointwise; the omnibus test remains the inferential anchor",
      x = time_col,
      y = "Fitted difference"
    ) +
    trajectory_plot_theme()

  if (isTRUE(show_points)) {
    plot <- plot + ggplot2::geom_point(
      data = contrasts[contrasts$.separated, , drop = FALSE],
      shape = 21,
      fill = "white",
      size = 1.8,
      stroke = 0.65
    )
  }
  plot
}


#' Plot trajectory moment estimates
#'
#' @param moments Output from [trajectory_moments()].
#' @param moment Character vector of moments to display.
#' @param window Optional window labels to display.
#' @param feature Optional feature names to display.
#' @param palette Optional group color vector.
#'
#' @return A faceted `ggplot` object.
#' @export
plot_trajectory_moments <- function(
    moments,
    moment = c("mean_level", "auc", "amplitude", "linear_trend"),
    window = NULL,
    feature = NULL,
    palette = NULL
) {
  if (!inherits(moments, "levaim_trajectory_moments")) {
    cli::cli_abort("`moments` must be created with `trajectory_moments()`.")
  }
  data <- moments[moments$moment %in% moment, , drop = FALSE]
  if (!is.null(window)) {
    data <- data[data$window_label %in% window, , drop = FALSE]
  }
  if (!is.null(feature)) {
    data <- data[data$feature %in% feature, , drop = FALSE]
  }
  if (nrow(data) == 0L) {
    cli::cli_abort("No requested moment rows are available.")
  }
  data$.feature <- short_taxon_name(data$feature)
  data$.panel <- paste(data$.feature, data$window_label, data$moment, sep = " | ")
  colors <- trajectory_plot_palette(unique(data$group_level), palette)

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$estimate,
      y = .data$group_level,
      color = .data$group_level
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "#D0D0D0", linewidth = 0.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      width = 0.18,
      linewidth = 0.65,
      orientation = "y"
    ) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::facet_wrap(stats::as.formula("~ .panel"), scales = "free_x") +
    ggplot2::scale_color_manual(values = colors, guide = "none") +
    ggplot2::labs(
      title = "Trajectory moment report",
      subtitle = "Points are fitted moments; intervals propagate model uncertainty",
      x = "Moment estimate",
      y = NULL
    ) +
    trajectory_plot_theme()
}


#' Plot between-group trajectory moment differences
#'
#' @param comparisons Output from [compare_trajectory_moments()].
#' @param moment Character vector of moments to display.
#' @param window Optional window labels to display.
#' @param feature Optional feature names to display.
#' @param q_threshold Optional q-value threshold used to emphasize estimates.
#'
#' @return A faceted `ggplot` object.
#' @export
plot_trajectory_moment_comparisons <- function(
    comparisons,
    moment = c("mean_level", "auc", "amplitude", "linear_trend"),
    window = NULL,
    feature = NULL,
    q_threshold = NULL
) {
  if (!inherits(comparisons, "levaim_moment_comparisons")) {
    cli::cli_abort("`comparisons` must be created with `compare_trajectory_moments()`.")
  }
  data <- comparisons[comparisons$moment %in% moment, , drop = FALSE]
  if (!is.null(window)) {
    data <- data[data$window_label %in% window, , drop = FALSE]
  }
  if (!is.null(feature)) {
    data <- data[data$feature %in% feature, , drop = FALSE]
  }
  if (nrow(data) == 0L) {
    cli::cli_abort("No requested moment-comparison rows are available.")
  }
  data$.comparison <- paste(data$group_1, data$group_2, sep = " - ")
  data$.feature <- short_taxon_name(data$feature)
  data$.panel <- paste(data$.feature, data$window_label, data$moment, sep = " | ")
  data$.supported <- if (is.null(q_threshold)) {
    data$lower > 0 | data$upper < 0
  } else {
    is.finite(data$q_value) & data$q_value <= q_threshold
  }

  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$difference,
      y = .data$.comparison,
      color = .data$.supported
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "#555555", linewidth = 0.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      width = 0.18,
      linewidth = 0.65,
      orientation = "y"
    ) +
    ggplot2::geom_point(size = 2.3) +
    ggplot2::facet_wrap(stats::as.formula("~ .panel"), scales = "free_x") +
    ggplot2::scale_color_manual(
      values = c(`TRUE` = "#0072B2", `FALSE` = "#7F7F7F"),
      name = "Supported"
    ) +
    ggplot2::labs(
      title = "Between-group moment differences",
      subtitle = "Intervals use fitted uncertainty draws; inferential method is recorded in the source table",
      x = "Group difference",
      y = NULL
    ) +
    trajectory_plot_theme()
}


#' @keywords internal
trajectory_plot_palette <- function(groups, palette = NULL) {
  defaults <- c("#0072B2", "#D55E00", "#009E73", "#7F7F7F", "#CC79A7", "#E69F00")
  if (is.null(palette)) {
    palette <- rep(defaults, length.out = length(groups))
  }
  if (length(palette) < length(groups)) {
    cli::cli_abort("`palette` must provide at least one color per displayed group.")
  }
  if (is.null(names(palette))) {
    names(palette) <- groups
  }
  palette[groups]
}


#' @keywords internal
trajectory_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}


#' @keywords internal
plot_group_label <- function(by_cols) {
  if (length(by_cols) == 0L) "Trajectory" else paste(by_cols, collapse = " : ")
}


#' @keywords internal
short_response_label <- function(response) {
  if (response$source == "assay") {
    return(short_taxon_name(response$feature))
  }
  response$column
}


#' @keywords internal
short_taxon_name <- function(x) {
  sub(".*[|]s__", "", x)
}

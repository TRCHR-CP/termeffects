#' Heatmap of a prediction surface
#'
#' Shows how the prediction changes over two continuous focal variables, with
#' all other predictors held constant (listed in the subtitle, or in the
#' caption of [term_figure()] when there are many). Contour lines make it
#' easier to read off values than colour alone.
#'
#' @param x A [term_predictions()] result with two continuous focal
#'   variables. A third, categorical focal variable is shown as panels.
#' @param contours Draw contour lines? `TRUE`, `FALSE`, or the number of
#'   contour levels.
#' @param cuts Optional named list marking where cuts are taken, e.g.
#'   `list(wt = c(2.5, 3.5))` draws dashed lines at those values. Use the same
#'   values in [term_predictions()] and [term_curves()] to show the cuts.
#' @param show_data Overlay the observed values of the two focal variables,
#'   to show where the surface is supported by data and where it extrapolates.
#' @param fill What the colour shows: the `"estimate"` (default), the
#'   `"band_width"` (width of the simultaneous band, on the outcome scale) or
#'   the `"std.error"` (on the scale the band is computed on). The last two
#'   show where the surface is estimated precisely.
#' @param show_support Fade the parts of the surface where there are few
#'   observations: outside the smallest region that holds `support_prob` of
#'   the data, estimated with a kernel density of the varying predictors.
#' @param support_prob Share of the data inside the region that is not faded.
#' @param layers Optional [term_layer()] result (or a list of them) whose
#'   curves are drawn on the heatmap.
#' @param show_direction Draw arrows along layer curves, pointing in the
#'   direction of their first variable.
#'
#' @return A ggplot object, which can be modified with the usual ggplot2
#'   functions.
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * splines::ns(wt, 3) + am, data = mtcars)
#' pred <- term_predictions(fit, c("hp", "wt"), n = 40)
#' term_heatmap(pred, cuts = list(wt = c(2.5, 3.5)), show_data = TRUE)
#' term_heatmap(pred, fill = "band_width", show_support = TRUE)
#' @export
term_heatmap <- function(x, contours = TRUE, cuts = NULL, show_data = FALSE,
                         fill = c("estimate", "band_width", "std.error"),
                         show_support = FALSE, support_prob = 0.975,
                         layers = NULL, show_direction = TRUE) {
  fill <- match.arg(fill)
  focal <- attr(x, "focal")
  if (length(focal) < 2 || !is.numeric(x[[focal[1]]]) || !is.numeric(x[[focal[2]]])) {
    stop("A heatmap needs two continuous focal variables.", call. = FALSE)
  }
  vx <- focal[1]
  vy <- focal[2]

  plot_data <- as.data.frame(x)
  plot_data$.fill <- switch(fill,
    estimate = plot_data$estimate,
    band_width = plot_data$band.high - plot_data$band.low,
    std.error = plot_data$std.error
  )
  fill_name <- switch(fill,
    estimate = attr(x, "outcome_label"),
    band_width = paste0("Width of ", 100 * attr(x, "ci_level"), "% simultaneous band"),
    std.error = "Standard error"
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[vx]], y = .data[[vy]])) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$.fill), interpolate = TRUE) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  p <- p + if (fill == "estimate") {
    ggplot2::scale_fill_viridis_c(name = .wrap(fill_name))
  } else {
    # uncertainty: one hue, light (precise) to dark (uncertain)
    ggplot2::scale_fill_gradient(name = .wrap(fill_name), low = "#F7FBFF", high = "#08306B")
  }

  if (!isFALSE(contours)) {
    bins <- if (is.numeric(contours)) contours else 8
    p <- p + ggplot2::geom_contour(ggplot2::aes(z = .data$.fill),
                                   bins = bins, colour = "white", linewidth = 0.3, alpha = 0.8)
  }

  if (show_support) {
    support <- .support(x, support_prob)
    if (!is.null(support)) {
      plot_data$.veil <- ifelse(support, 0, 0.6)
      p <- p + ggplot2::geom_raster(data = plot_data, ggplot2::aes(alpha = .data$.veil),
                                    fill = "white", interpolate = TRUE) +
        ggplot2::scale_alpha_identity()
    }
  }

  if (show_data && !is.null(attr(x, "data"))) {
    observed <- attr(x, "data")
    if (all(c(vx, vy) %in% names(observed))) {
      p <- p + ggplot2::geom_point(data = observed, ggplot2::aes(x = .data[[vx]], y = .data[[vy]]),
                                   inherit.aes = FALSE, size = 0.6, alpha = 0.4, colour = "black")
    }
  }

  for (name in names(cuts)) {
    if (name == vx) {
      p <- p + ggplot2::geom_vline(xintercept = cuts[[name]], linetype = "dashed", colour = "white")
    } else if (name == vy) {
      p <- p + ggplot2::geom_hline(yintercept = cuts[[name]], linetype = "dashed", colour = "white")
    } else {
      stop("`cuts` must name one of the two heatmap variables: ", vx, " or ", vy, ".", call. = FALSE)
    }
  }

  if (!is.null(layers)) {
    p <- .add_layers(p, layers, vx, vy, show_direction)
  }

  if (length(focal) >= 3) {
    p <- p + ggplot2::facet_wrap(focal[3], labeller = ggplot2::label_both)
  }

  p + ggplot2::coord_cartesian(xlim = range(plot_data[[vx]]), ylim = range(plot_data[[vy]]), expand = FALSE) +
    ggplot2::labs(x = vx, y = vy, subtitle = .held_subtitle(x),
                  caption = if (show_support) paste0("Faded: outside the region holding ",
                                                     100 * support_prob, "% of the observations")) +
    ggplot2::theme_minimal()
}


# curves of one or more layers drawn on a heatmap
.add_layers <- function(p, layers, vx, vy, show_direction) {
  if (inherits(layers, "term_predictions")) {
    layers <- list(layers)
  }
  paths <- lapply(seq_along(layers), function(i) {
    layer <- layers[[i]]
    if (!all(c(vx, vy) %in% names(layer))) {
      stop("A layer must follow the heatmap variables ", vx, " and ", vy, ".", call. = FALSE)
    }
    focal <- attr(layer, "focal")
    curve_vars <- intersect(focal[-1], names(layer))
    curve_vars <- curve_vars[vapply(curve_vars, function(v) !v %in% c(vx, vy), logical(1))]
    d <- as.data.frame(layer)
    label <- if (length(curve_vars)) {
      paste0(curve_vars[1], " = ", format(d[[curve_vars[1]]], digits = 3, trim = TRUE))
    } else {
      rep(focal[1], nrow(d))
    }
    d$.curve <- if (length(layers) > 1) paste0("Layer ", i, ": ", label) else label
    d <- d[order(d$.curve, d[[focal[1]]]), , drop = FALSE]
    d
  })
  paths <- do.call(rbind, lapply(paths, function(d) d[c(vx, vy, ".curve", setdiff(names(d), c(vx, vy, ".curve")))]))
  paths$.curve <- factor(paths$.curve, levels = unique(paths$.curve))
  arrow <- if (show_direction) {
    ggplot2::arrow(length = ggplot2::unit(0.12, "inches"), type = "closed")
  }
  p + ggplot2::geom_path(data = paths, ggplot2::aes(x = .data[[vx]], y = .data[[vy]],
                                                    group = .data$.curve, colour = .data$.curve),
                         inherit.aes = FALSE, linewidth = 0.9, arrow = arrow) +
    .curve_colours("Curve", nlevels(paths$.curve), aesthetics = "colour")
}


# which grid points lie in the region holding `prob` of the data, by a
# Gaussian kernel density of the continuous predictors that vary on the grid.
# Kernel and bandwidth follow Lukas Bolsinger's version in p-ggeffects; the
# threshold is the (1 - prob) quantile of the density at the observations
# (the usual sample estimate of a highest-density region).
.support <- function(x, prob) {
  observed <- attr(x, "data")
  if (is.null(observed)) {
    warning("The model data are not available, so the support cannot be shown.", call. = FALSE)
    return(NULL)
  }
  vars <- names(observed)[vapply(names(observed), function(v) {
    is.numeric(observed[[v]]) && v %in% names(x) && length(unique(x[[v]])) > 1
  }, logical(1))]
  if (!length(vars)) {
    return(NULL)
  }
  X <- as.matrix(observed[vars])
  X <- X[stats::complete.cases(X), , drop = FALSE]
  density <- .kernel_density(X)
  threshold <- stats::quantile(density(X), 1 - prob, names = FALSE)
  density(as.matrix(as.data.frame(x)[vars])) >= threshold
}

.kernel_density <- function(X) {
  n <- nrow(X)
  d <- ncol(X)
  H <- n^(-2 / (d + 4)) * stats::cov(X)
  H_inv <- solve(H)
  norm <- sqrt((2 * pi)^d * det(H))
  function(points) {
    vapply(seq_len(nrow(points)), function(i) {
      delta <- sweep(X, 2, points[i, ])
      mean(exp(-0.5 * rowSums((delta %*% H_inv) * delta))) / norm
    }, numeric(1))
  }
}


#' Curves through a prediction surface ("cuts") with confidence bands
#'
#' Plots the prediction against the first focal variable, with one curve per
#' value of the second focal variable. Each curve is a cut through the
#' surface shown by [term_heatmap()].
#'
#' The **simultaneous** band (default) covers all plotted curves at once with
#' the stated probability, so it can be used to judge their shapes, e.g.
#' whether a flat line fits inside it. The **pointwise** interval only holds
#' at each single value.
#'
#' @param x A [term_predictions()], [term_layer()] or
#'   [test_simultaneous_predictions()] result. The second focal variable, if
#'   any, should have a few values (use `values = ` in [term_predictions()]).
#'   For a test result, the reference function is drawn as a dashed line.
#' @param band Which uncertainty to show: `"simultaneous"`, `"pointwise"`,
#'   `"both"` (simultaneous shaded, pointwise dashed) or `"none"`.
#' @param show_direction Draw arrows pointing in the direction of the first
#'   variable (useful for layers, whose first variable is a position along a
#'   curve).
#'
#' @return A ggplot object.
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * wt + am, data = mtcars)
#' cuts <- term_predictions(fit, c("hp", "wt"), values = list(wt = c(2.5, 3.5)))
#' term_curves(cuts, band = "both")
#' @export
term_curves <- function(x, band = c("simultaneous", "pointwise", "both", "none"),
                        show_direction = FALSE) {
  band <- match.arg(band)
  focal <- attr(x, "focal")
  vx <- focal[1]
  if (!is.numeric(x[[vx]])) {
    stop("The first focal variable must be continuous.", call. = FALSE)
  }
  x <- as.data.frame(x)

  group <- if (length(focal) >= 2) focal[2] else NULL
  if (!is.null(group)) {
    if (length(unique(x[[group]])) > 8) {
      stop("`", group, "` has too many values for separate curves. Choose a few with ",
           "`values = list(", group, " = c(...))` in term_predictions(), ",
           "or use term_heatmap().", call. = FALSE)
    }
    x[[group]] <- factor(x[[group]], levels = unique(sort(x[[group]])))
  }
  colour_aes <- if (is.null(group)) NULL else ggplot2::aes(colour = .data[[group]], fill = .data[[group]])

  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[vx]], y = .data$estimate, group = if (is.null(group)) 1 else .data[[group]])) +
    colour_aes

  if (band %in% c("simultaneous", "both")) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$band.low, ymax = .data$band.high),
                                  alpha = 0.2, colour = NA)
  }
  if (band == "pointwise") {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$conf.low, ymax = .data$conf.high),
                                  alpha = 0.2, colour = NA)
  }
  if (band == "both") {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = .data$conf.low), linetype = "dashed", linewidth = 0.4) +
      ggplot2::geom_line(ggplot2::aes(y = .data$conf.high), linetype = "dashed", linewidth = 0.4)
  }
  if ("reference" %in% names(x)) {
    p <- p + ggplot2::geom_line(ggplot2::aes(y = .data$reference), linetype = "longdash",
                                linewidth = 0.6, colour = "grey20")
  }

  if (length(focal) >= 3) {
    p <- p + ggplot2::facet_wrap(focal[3], labeller = ggplot2::label_both)
  }

  arrow <- if (show_direction) ggplot2::arrow(length = ggplot2::unit(0.1, "inches"), type = "closed")
  caption <- .band_caption(x, band)
  if ("reference" %in% names(x)) {
    caption <- paste0(c(caption, "Long dashes: reference function"), collapse = "\n")
  }
  if (!is.null(group)) {
    p <- p + .curve_colours(group, nlevels(x[[group]]))
  } else {
    p <- p + ggplot2::aes(colour = NULL, fill = NULL)
  }
  p + ggplot2::geom_line(linewidth = 0.8, arrow = arrow) +
    ggplot2::labs(x = vx, y = attr(x, "outcome_label"), subtitle = .held_subtitle(x),
                  caption = caption) +
    ggplot2::theme_minimal()
}


#' Plot term predictions
#'
#' Chooses [term_heatmap()] when the first two focal variables are continuous
#' with many values, and [term_curves()] otherwise.
#'
#' @param x A [term_predictions()] result.
#' @param type `"auto"`, `"heatmap"` or `"curves"`.
#' @param ... Passed to [term_heatmap()] or [term_curves()].
#' @return A ggplot object.
#' @export
plot.term_predictions <- function(x, type = c("auto", "heatmap", "curves"), ...) {
  type <- match.arg(type)
  if (type == "auto") {
    type <- if (.is_surface(x)) "heatmap" else "curves"
  }
  if (type == "heatmap") term_heatmap(x, ...) else term_curves(x, ...)
}


#' Plot a simultaneous test
#'
#' @param x A [test_simultaneous_predictions()] result.
#' @param type `"curves"`: the predictions with their band and the reference
#'   function (dashed). `"difference"`: estimate minus reference, with the
#'   band shifted the same way; the reference fits inside the band wherever
#'   the band contains zero. For a surface, `"difference"` is a heatmap in
#'   which points where the reference lies outside the band are marked.
#'   `"auto"` uses `"difference"` for surfaces and `"curves"` otherwise.
#' @param ... Passed to [term_curves()].
#' @return A ggplot object.
#' @export
plot.term_test <- function(x, type = c("auto", "curves", "difference"), ...) {
  type <- match.arg(type)
  surface <- .is_surface(x)
  if (type == "auto") {
    type <- if (surface) "difference" else "curves"
  }
  if (type == "curves") {
    return(term_curves(x, ...))
  }
  subtitle <- sprintf("Simultaneous test: max |t| = %s, p = %s",
                      paste(sprintf("%.2f", attr(x, "statistic")), collapse = ", "),
                      paste(format.pval(attr(x, "p_value"), digits = 2, eps = 1e-4), collapse = ", "))
  focal <- attr(x, "focal")
  d <- as.data.frame(x)
  if (surface) {
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[focal[1]]], y = .data[[focal[2]]])) +
      ggplot2::geom_raster(ggplot2::aes(fill = .data$difference), interpolate = TRUE) +
      ggplot2::scale_fill_gradient2(name = "Estimate minus\nreference",
                                    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0) +
      ggplot2::geom_point(data = d[!d$inside, , drop = FALSE], shape = 4, size = 1, colour = "grey20") +
      ggplot2::coord_cartesian(expand = FALSE) +
      ggplot2::labs(caption = "Crosses: the reference lies outside the simultaneous band")
    if (length(focal) >= 3) {
      p <- p + ggplot2::facet_wrap(focal[3], labeller = ggplot2::label_both)
    }
    return(p + ggplot2::labs(subtitle = subtitle) + ggplot2::theme_minimal())
  }
  d$estimate <- d$difference
  d$band.low <- d$band.low - d$reference
  d$band.high <- d$band.high - d$reference
  d$conf.low <- d$conf.low - d$reference
  d$conf.high <- d$conf.high - d$reference
  d$reference <- NULL
  attributes(d)[c("focal", "held", "held_rule", "ci_level")] <- attributes(x)[c("focal", "held", "held_rule", "ci_level")]
  attr(d, "outcome_label") <- "Estimate minus reference"
  term_curves(d, ...) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(title = subtitle)
}


#' @export
plot.term_layer <- function(x, ...) {
  term_curves(x, ...)
}


.is_surface <- function(x) {
  focal <- attr(x, "focal")
  length(focal) >= 2 && is.numeric(x[[focal[1]]]) && is.numeric(x[[focal[2]]]) &&
    length(unique(x[[focal[2]]])) > 8
}


# Okabe-Ito colours in a fixed order (colour-blind safe); at most 8 curves
.curve_colours <- function(name, n, aesthetics = c("colour", "fill")) {
  palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
  if (n > length(palette)) {
    return(NULL)
  }
  ggplot2::scale_colour_manual(name = name, values = palette, aesthetics = aesthetics)
}


# "Held constant: age = 61, sex = Male"; with many held variables, a pointer
# to the full list instead
.held_subtitle <- function(x, max_shown = 4) {
  held <- .held_text(x, rules = FALSE)
  if (!length(held)) {
    return(NULL)
  }
  if (length(held) > max_shown) {
    return(.wrap(paste(length(held), "other covariates held at their median or most common",
                       "value (see the caption from term_figure())"), 70))
  }
  .wrap(paste0("Held constant: ", paste(held, collapse = ", ")), 90)
}


.band_caption <- function(x, band) {
  level <- paste0(100 * attr(x, "ci_level"), "%")
  switch(band,
    simultaneous = paste("Shaded:", level, "simultaneous confidence band"),
    pointwise = paste("Shaded:", level, "pointwise confidence intervals"),
    both = paste("Shaded:", level, "simultaneous band; dashed:", level, "pointwise intervals"),
    none = NULL
  )
}


.wrap <- function(text, width = 25) {
  paste(strwrap(text, width), collapse = "\n")
}

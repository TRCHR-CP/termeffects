#' Heatmap of a prediction surface
#'
#' Shows how the prediction changes over two continuous focal variables, with
#' all other predictors held constant (listed in the subtitle). Contour lines
#' make it easier to read off values than colour alone.
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
#'
#' @return A ggplot object, which can be modified with the usual ggplot2
#'   functions.
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * splines::ns(wt, 3) + am, data = mtcars)
#' pred <- term_predictions(fit, c("hp", "wt"), n = 40)
#' term_heatmap(pred, cuts = list(wt = c(2.5, 3.5)), show_data = TRUE)
#' @export
term_heatmap <- function(x, contours = TRUE, cuts = NULL, show_data = FALSE) {
  focal <- attr(x, "focal")
  if (length(focal) < 2 || !is.numeric(x[[focal[1]]]) || !is.numeric(x[[focal[2]]])) {
    stop("A heatmap needs two continuous focal variables.", call. = FALSE)
  }
  vx <- focal[1]
  vy <- focal[2]

  p <- ggplot2::ggplot(x, ggplot2::aes(x = .data[[vx]], y = .data[[vy]])) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$estimate), interpolate = TRUE) +
    ggplot2::scale_fill_viridis_c(name = .wrap(attr(x, "outcome_label"))) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))

  if (!isFALSE(contours)) {
    bins <- if (is.numeric(contours)) contours else 8
    p <- p + ggplot2::geom_contour(ggplot2::aes(z = .data$estimate),
                                   bins = bins, colour = "white", linewidth = 0.3, alpha = 0.8)
  }

  if (show_data && !is.null(attr(x, "data"))) {
    observed <- attr(x, "data")
    p <- p + ggplot2::geom_point(data = observed, ggplot2::aes(x = .data[[vx]], y = .data[[vy]]),
                                 inherit.aes = FALSE, size = 0.6, alpha = 0.4, colour = "black")
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

  if (length(focal) >= 3) {
    p <- p + ggplot2::facet_wrap(focal[3], labeller = ggplot2::label_both)
  }

  p + ggplot2::labs(x = vx, y = vy, subtitle = .held_subtitle(x)) +
    ggplot2::theme_minimal()
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
#' @param x A [term_predictions()] result. The second focal variable, if any,
#'   should have a few values (use `values = ` in [term_predictions()]).
#' @param band Which uncertainty to show: `"simultaneous"`, `"pointwise"`,
#'   `"both"` (simultaneous shaded, pointwise dashed) or `"none"`.
#'
#' @return A ggplot object.
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * wt + am, data = mtcars)
#' cuts <- term_predictions(fit, c("hp", "wt"), values = list(wt = c(2.5, 3.5)))
#' term_curves(cuts, band = "both")
#' @export
term_curves <- function(x, band = c("simultaneous", "pointwise", "both", "none")) {
  band <- match.arg(band)
  focal <- attr(x, "focal")
  vx <- focal[1]
  if (!is.numeric(x[[vx]])) {
    stop("The first focal variable must be continuous.", call. = FALSE)
  }

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

  if (length(focal) >= 3) {
    p <- p + ggplot2::facet_wrap(focal[3], labeller = ggplot2::label_both)
  }

  p + ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::labs(x = vx, y = attr(x, "outcome_label"), subtitle = .held_subtitle(x),
                  caption = .band_caption(x, band)) +
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
    focal <- attr(x, "focal")
    surface <- length(focal) >= 2 && is.numeric(x[[focal[1]]]) && is.numeric(x[[focal[2]]]) &&
      length(unique(x[[focal[2]]])) > 8
    type <- if (surface) "heatmap" else "curves"
  }
  if (type == "heatmap") term_heatmap(x, ...) else term_curves(x, ...)
}


# "Held constant: age = 61.2, sex = female"
.held_subtitle <- function(x) {
  held <- attr(x, "held")
  if (!length(held)) {
    return(NULL)
  }
  values <- vapply(held, function(v) {
    if (is.numeric(v)) format(signif(v, 3)) else as.character(v)
  }, character(1))
  .wrap(paste0("Held constant: ", paste0(names(held), " = ", values, collapse = ", ")), 90)
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

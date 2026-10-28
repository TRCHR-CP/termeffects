#' Predictions along curves through a heatmap ("layers")
#'
#' A heatmap shows a prediction surface over two variables, say age and
#' creatinine. A layer follows curves through that surface, for example
#' patients with the same eGFR, which combines age and creatinine. Each curve
#' is described by a parameterization: the two heatmap variables written as
#' functions of at most two new variables, e.g. a position `t` along the
#' curve and a value `r` that picks the curve.
#'
#' The other predictors are held at the same values as in the heatmap, so the
#' layer is a set of cuts through the same surface. Plot it on the heatmap
#' with `term_heatmap(x, layers = layer)` and as curves with
#' [term_curves()].
#'
#' The layer's band is simultaneous over the points of the layer (not over
#' the whole heatmap), so it is the right band for judging the curves.
#'
#' The concept and the restrictions are from Lukas Bolsinger's
#' `ggeffects_layer` in p-ggeffects.
#'
#' @param x A [term_predictions()] result with two continuous focal
#'   variables, as drawn by [term_heatmap()].
#' @param parameterization Named list giving one or both heatmap variables as
#'   functions (or strings such as `"t + 60"`) of at most two new variables.
#' @param terms The new variables: first the position along a curve, then
#'   (optionally) the variable that picks the curve. Defaults to the new
#'   variables in the order they appear.
#' @param values,n Values of the new variables, as in [term_predictions()].
#'   Give `values` for the variable that picks the curves.
#' @param ... Passed to [term_predictions()], e.g. `ci_method`. By default the
#'   layer uses the same confidence level, scale and method as `x`.
#'
#' @return A [term_predictions()] result of class `term_layer`, whose columns
#'   also contain the two heatmap variables along each curve.
#'
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * splines::ns(wt, 3), data = mtcars)
#' surface <- term_predictions(fit, c("hp", "wt"), n = 30)
#' # curves of constant power-to-weight ratio hp / wt = r
#' layer <- term_layer(surface, list(hp = function(t, r) r * t, wt = function(t) t),
#'                     terms = c("t", "r"), values = list(t = seq(2, 4, by = 0.1), r = c(30, 45, 60)))
#' term_heatmap(surface, layers = layer)
#' term_curves(layer)
#' @export
term_layer <- function(x, parameterization, terms = NULL, values = NULL, n = 30, ...) {
  base <- attr(x, "focal")
  if (length(base) < 2 || !is.numeric(x[[base[1]]]) || !is.numeric(x[[base[2]]])) {
    stop("`x` must have two continuous focal variables (a heatmap).", call. = FALSE)
  }
  base <- base[1:2]
  if (!length(parameterization) || !all(names(parameterization) %in% base)) {
    stop("`parameterization` may only give the heatmap variables ", toString(base), ".", call. = FALSE)
  }

  new_parameterizations <- as.functions(parameterization)
  # within one parameterization, an argument named like a parameterized
  # heatmap variable would be that variable's held value, not its value on
  # the curve, which is almost never what is meant
  clash <- intersect(.get_parameters(new_parameterizations), names(parameterization))
  clash <- clash[vapply(clash, function(v) {
    !identical(body(new_parameterizations[[v]]), quote((v)))
  }, logical(1))]
  if (length(clash)) {
    stop("The parameterization uses ", toString(clash), " as an input while also replacing it. ",
         "Write every heatmap variable in terms of the new variables, e.g. ",
         "list(age = function(t) t, creatinine = function(t, egfr) f(egfr, age = t)).", call. = FALSE)
  }
  pmodel <- parameterize(attr(x, "model"), parameterization)
  new_vars <- .get_parameters(new_parameterizations)
  new_vars <- setdiff(new_vars, c(base, names(attr(x, "held"))))
  if (length(new_vars) > 2) {
    stop("A layer can have at most two new variables (along a curve, and between curves); got ",
         toString(new_vars), ".", call. = FALSE)
  }
  if (is.null(terms)) {
    terms <- new_vars
  }
  # a heatmap variable that is not parameterized stays a free coordinate
  free <- setdiff(base, names(parameterization))
  if (length(free) && !any(free %in% terms)) {
    terms <- c(terms, free)
  }
  # panels of the heatmap (third focal variable) are kept
  terms <- c(terms, setdiff(attr(x, "focal")[-(1:2)], terms))

  args <- list(...)
  defaults <- list(ci_level = attr(x, "ci_level"), scale = attr(x, "scale"),
                   ci_method = attr(x, "ci_method"))
  args <- c(args, defaults[setdiff(names(defaults), names(args))])

  held <- attr(x, "held")
  held <- held[setdiff(names(held), terms)]
  out <- do.call(term_predictions, c(list(pmodel, terms, values = values, at = held, n = n), args))
  # position of each point in the heatmap's coordinates (needed when those
  # are themselves variables of a parameterization, and not model predictors)
  positions <- .parameterize_data_frame(as.data.frame(out), new_parameterizations)
  out[names(parameterization)] <- positions[names(parameterization)]
  attr(out, "base_focal") <- base
  attr(out, "held_rule") <- attr(x, "held_rule")[names(held)]
  class(out) <- c("term_layer", class(out))
  out
}

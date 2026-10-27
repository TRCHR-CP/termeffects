#' Predictions with pointwise and simultaneous confidence bands
#'
#' Computes model predictions over a grid of focal variables, holding all
#' other predictors constant (see [term_grid()]), together with
#'
#' * a **pointwise** confidence interval at each grid point, and
#' * a **simultaneous** confidence band that covers the true values at all
#'   grid points at once with probability `ci_level`.
#'
#' Use the pointwise interval to make a statement about one point. Use the
#' simultaneous band to judge the shape of a whole curve or surface, for
#' example whether a flat line or a reference curve fits inside it.
#'
#' Both intervals are computed on the link scale (e.g. log-odds) and then
#' transformed, so on the response scale they always respect its limits
#' (e.g. probabilities stay within 0 and 1).
#'
#' The simultaneous band is simultaneous over the grid points supplied, with
#' the other predictors held at the stated values. A finer grid gives a band
#' that is closer to one over the whole continuous range.
#'
#' @inheritParams term_grid
#' @param newdata Optional grid to predict on, e.g. from [term_grid()] or a
#'   hand-made data frame tracing a curve through a surface. If given,
#'   `values`, `at` and `n` are ignored. It must contain all predictors.
#' @param ci_level Confidence level for both intervals.
#' @param scale `"response"` (default) gives predictions on the outcome scale
#'   (e.g. probabilities for logistic models, hazard ratios for Cox models,
#'   the mean outcome for `rms::orm` models); `"link"` gives them on the
#'   scale of the linear predictor.
#' @param nsim Number of simulation draws for the simultaneous critical value.
#'
#' @return A data frame of class `term_predictions` with the grid columns and
#'   * `estimate`: the prediction,
#'   * `std.error`: its standard error **on the link scale**,
#'   * `conf.low`, `conf.high`: pointwise interval,
#'   * `band.low`, `band.high`: simultaneous band.
#'
#'   Attributes: `critical_value` (simultaneous), `pointwise_critical_value`,
#'   `ci_level`, `df`, `scale`, `focal`, `held`.
#'
#' @details Cox models: the estimate is the hazard ratio (or log hazard ratio
#'   with `scale = "link"`) relative to a subject with average values of the
#'   model's design-matrix columns, matching `predict(type = "lp")`.
#'
#'   `rms::orm` models: the link scale is rms's linear predictor (using the
#'   model's reference intercept) and the response scale is the mean outcome,
#'   from `Hmisc::Mean()`. Intervals are computed for the linear predictor and
#'   transformed, so they stay within the observed range of the outcome. As in
#'   `rms::Predict(fit, fun = Mean(fit))`, they include the uncertainty of the
#'   covariate effects and the reference intercept, but not that of the other
#'   intercepts.
#'
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * wt + am, data = mtcars)
#' set.seed(1)
#' pred <- term_predictions(fit, c("hp", "wt"), values = list(wt = c(2.5, 3.5)))
#' head(pred)
#' attr(pred, "critical_value")
#' @export
term_predictions <- function(model,
                             terms,
                             values = NULL,
                             at = NULL,
                             n = 30,
                             newdata = NULL,
                             data = NULL,
                             ci_level = 0.95,
                             scale = c("response", "link"),
                             nsim = 10000) {
  scale <- match.arg(scale)
  if (is.null(newdata)) {
    newdata <- term_grid(model, terms, values = values, at = at, n = n, data = data)
  } else if (!missing(terms) && is.null(attr(newdata, "focal"))) {
    attr(newdata, "focal") <- terms
  }

  design <- .design(model, newdata)
  eta <- drop(design$X %*% design$beta)
  se <- sqrt(pmax(rowSums((design$X %*% design$vcov) * design$X), 0))

  df <- .residual_df(model)
  alpha <- 1 - ci_level
  crit_pointwise <- if (is.finite(df)) stats::qt(1 - alpha / 2, df) else stats::qnorm(1 - alpha / 2)
  crit <- .critical_value(design$X, design$vcov, ci_level = ci_level, df = df, nsim = nsim)

  transform <- if (scale == "response") .inverse_link(model) else identity

  out <- as.data.frame(newdata)
  out$estimate <- transform(eta)
  out$std.error <- se
  out$conf.low <- transform(eta - crit_pointwise * se)
  out$conf.high <- transform(eta + crit_pointwise * se)
  out$band.low <- transform(eta - crit * se)
  out$band.high <- transform(eta + crit * se)

  attr(out, "critical_value") <- crit
  attr(out, "pointwise_critical_value") <- crit_pointwise
  attr(out, "ci_level") <- ci_level
  attr(out, "df") <- df
  attr(out, "scale") <- scale
  attr(out, "focal") <- attr(newdata, "focal")
  attr(out, "held") <- attr(newdata, "held")
  attr(out, "outcome_label") <- .outcome_label(model, scale)
  attr(out, "data") <- tryCatch(.model_data(model, data), error = function(e) NULL)
  class(out) <- c("term_predictions", "data.frame")
  out
}


# axis/legend label describing what `estimate` is
.outcome_label <- function(model, scale) {
  if (scale == "link") {
    return("Linear predictor")
  }
  if (inherits(model, "coxph")) {
    return("Hazard ratio (vs. average subject)")
  }
  response <- all.vars(stats::formula(model))[1]
  if (inherits(model, "orm")) {
    return(paste("Mean", response))
  }
  if (inherits(model, c("glm", "gam", "glmerMod")) &&
      stats::family(model)$family %in% c("binomial", "quasibinomial")) {
    return(paste("Probability of", response))
  }
  paste("Predicted", response)
}


#' @export
print.term_predictions <- function(x, n = 10, ...) {
  focal <- attr(x, "focal")
  held <- attr(x, "held")
  cat("Predictions over", toString(focal), "(", nrow(x), "grid points,", attr(x, "scale"), "scale )\n")
  if (length(held)) {
    held_text <- vapply(held, function(v) {
      if (is.numeric(v)) format(signif(v, 4)) else as.character(v)
    }, character(1))
    cat("Held constant:", paste0(names(held), " = ", held_text, collapse = ", "), "\n")
  }
  cat(sprintf(
    "Critical values for %g%% intervals: pointwise %.3f, simultaneous %.3f\n\n",
    100 * attr(x, "ci_level"), attr(x, "pointwise_critical_value"), attr(x, "critical_value")
  ))
  cols <- c(focal, "estimate", "conf.low", "conf.high", "band.low", "band.high")
  print(utils::head(as.data.frame(x)[intersect(cols, names(x))], n), digits = 4, row.names = FALSE)
  if (nrow(x) > n) cat("... with", nrow(x) - n, "more rows\n")
  invisible(x)
}

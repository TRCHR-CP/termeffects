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
#' example whether a flat line or a reference curve fits inside it (see
#' [test_simultaneous_predictions()]).
#'
#' The simultaneous band is simultaneous over the grid points supplied, with
#' the other predictors held at the stated values. A finer grid gives a band
#' that is closer to one over the whole continuous range.
#'
#' @section How the intervals are computed:
#' With `ci_method = "analytic"` (default), intervals come from the
#' coefficients' covariance matrix. They are computed on the link scale
#' (e.g. log-odds) and then transformed, so they respect the limits of the
#' outcome (e.g. probabilities stay within 0 and 1). The simultaneous critical
#' value is found by simulating from the joint distribution of the estimates.
#'
#' With `ci_method = "bootstrap"`, the model is refitted `nboot` times on
#' resampled data (for `rms` fits with [rms::bootcov()], which needs the
#' model fitted with `x = TRUE, y = TRUE`). Standard errors are the standard
#' deviations of the bootstrap predictions, and the simultaneous critical
#' value is the `ci_level` quantile of the largest standardized deviation of
#' a bootstrap curve from the estimate (a "sup-t" band). Not available for
#' mixed models.
#'
#' `rms::orm` models on the response scale (the mean outcome): the analytic
#' intervals are computed for rms's linear predictor (reference intercept
#' plus covariate effects) and transformed with `Hmisc::Mean()`, as
#' `rms::Predict(fit, fun = Mean(fit))` does. They ignore the uncertainty of
#' the other intercepts, but in simulations (`validation/08_orm_interval_methods.R`)
#' they had the most accurate coverage. `ci_method = "delta"` instead uses
#' the delta method on the mean scale with the covariance of *all*
#' coefficients, as `rms::Mean(fit)(X = , conf.int = )` does; its symmetric
#' intervals covered slightly less often than stated. Bootstrap intervals
#' for orm are also computed on the mean scale. Mean-scale intervals are
#' clipped to the observed outcome range.
#'
#' Cox models: the estimate is the hazard ratio (or log hazard ratio with
#' `scale = "link"`) relative to a subject with average values of the
#' model's design-matrix columns, matching `predict(type = "lp")`.
#'
#' @inheritParams term_grid
#' @param newdata Optional grid to predict on, e.g. from [term_grid()] or a
#'   hand-made data frame tracing a curve through a surface. If given,
#'   `values`, `at` and `n` are ignored. It must contain all predictors (or,
#'   for a parameterized model, the variables of its parameterizations).
#' @param ci_level Confidence level for both intervals.
#' @param scale `"response"` (default) gives predictions on the outcome scale
#'   (e.g. probabilities for logistic models, hazard ratios for Cox models,
#'   the mean outcome for `rms::orm` models); `"link"` gives them on the
#'   scale of the linear predictor.
#' @param ci_method `"analytic"` (default), `"bootstrap"` or, for `rms::orm`
#'   models only, `"delta"`; see the section on how the intervals are
#'   computed.
#' @param ci_scope Over which points the band is simultaneous: `"global"`
#'   (default) all grid points together; `"faceted"` separately within each
#'   combination of the third and later focal variables; `"curve"` separately
#'   for each curve, i.e. each combination of the second and later focal
#'   variables. Narrower scopes give narrower bands, which only hold within
#'   each facet or curve.
#' @param nsim Number of simulation draws for the analytic simultaneous
#'   critical value.
#' @param nboot Number of bootstrap refits for `ci_method = "bootstrap"`.
#'
#' @return A data frame of class `term_predictions` with the grid columns
#'   (for a parameterized model also the predictors computed from them) and
#'   * `estimate`: the prediction,
#'   * `std.error`: its standard error on the scale the intervals are
#'     computed on (the link scale, except for delta-method and bootstrap
#'     intervals of the `orm` mean, where it is the mean scale),
#'   * `conf.low`, `conf.high`: pointwise interval,
#'   * `band.low`, `band.high`: simultaneous band.
#'
#'   Attributes include `critical_value` (simultaneous; one per facet or curve
#'   if `ci_scope` is not `"global"`), `pointwise_critical_value`,
#'   `ci_level`, `ci_method`, `ci_scope`, `scale`, `focal`, `held`,
#'   `held_rule` and the `model`.
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
                             ci_method = c("analytic", "bootstrap", "delta"),
                             ci_scope = c("global", "faceted", "curve"),
                             nsim = 10000,
                             nboot = 1000) {
  scale <- match.arg(scale)
  ci_method <- match.arg(ci_method)
  ci_scope <- match.arg(ci_scope)
  input <- model
  model <- original_model(model)

  if (is.null(newdata)) {
    newdata <- term_grid(input, terms, values = values, at = at, n = n, data = data)
  } else if (!missing(terms) && is.null(attr(newdata, "focal"))) {
    attr(newdata, "focal") <- terms
  }
  focal <- attr(newdata, "focal")
  grid <- as.data.frame(newdata)
  if (is.parameterized_model(input)) {
    grid <- parameterize_data(input, grid)
  }

  inference <- .inference(model, grid, scale, ci_method, nsim, nboot, data)
  df <- inference$df
  alpha <- 1 - ci_level
  crit_pointwise <- if (is.finite(df)) stats::qt(1 - alpha / 2, df) else stats::qnorm(1 - alpha / 2)

  # simultaneous critical value, for all points together or per facet/curve
  scope_vars <- switch(ci_scope, global = character(0), faceted = focal[-(1:2)], curve = focal[-1])
  scope_vars <- intersect(scope_vars, names(grid))
  group <- if (length(scope_vars)) interaction(grid[scope_vars], drop = TRUE, lex.order = TRUE) else
    factor(rep("all", nrow(grid)))
  max_t <- lapply(split(seq_len(nrow(grid)), group), inference$max_t)
  crit_by_group <- vapply(max_t, function(m) {
    if (length(m)) as.numeric(stats::quantile(m, ci_level, names = FALSE)) else NA_real_
  }, numeric(1))
  crit <- unname(crit_by_group[as.character(group)])

  eta <- inference$estimate
  se <- inference$se
  transform <- inference$transform

  out <- grid
  out$estimate <- transform(eta)
  out$std.error <- se
  out$conf.low <- transform(eta - crit_pointwise * se)
  out$conf.high <- transform(eta + crit_pointwise * se)
  out$band.low <- transform(eta - crit * se)
  out$band.high <- transform(eta + crit * se)

  attr(out, "critical_value") <- if (ci_scope == "global") unname(crit_by_group) else crit_by_group
  attr(out, "pointwise_critical_value") <- crit_pointwise
  attr(out, "ci_level") <- ci_level
  attr(out, "ci_method") <- ci_method
  attr(out, "ci_scope") <- ci_scope
  attr(out, "scope_group") <- group
  attr(out, "max_t") <- max_t
  attr(out, "to_inference_scale") <- inference$to_inference_scale
  attr(out, "df") <- df
  attr(out, "nboot") <- inference$nboot
  attr(out, "vcov_parts") <- inference$vcov_parts
  attr(out, "scale") <- scale
  attr(out, "focal") <- focal
  attr(out, "held") <- attr(newdata, "held")
  attr(out, "held_rule") <- attr(newdata, "held_rule")
  attr(out, "outcome_label") <- .outcome_label(model, scale)
  attr(out, "model") <- input
  attr(out, "data") <- tryCatch(.model_data(model, data), error = function(e) NULL)
  class(out) <- c("term_predictions", "data.frame")
  out
}


# Estimates and standard errors on the scale the intervals are built on,
# a function mapping that scale to the reported one, and a function giving
# simulated or bootstrap maxima of |deviation| / SE for a subset of points.
.inference <- function(model, grid, scale, ci_method, nsim, nboot, data) {
  is_orm <- inherits(model, "orm")
  if (ci_method == "delta" && !is_orm) {
    stop("ci_method = \"delta\" is only for rms::orm models; for other models the ",
         "analytic intervals already are delta-method intervals.", call. = FALSE)
  }
  # the orm mean computed on its own scale (delta method or bootstrap), as
  # opposed to intervals for the linear predictor transformed by Mean()
  mean_scale <- is_orm && scale == "response" && ci_method != "analytic"
  df <- if (ci_method == "analytic") .residual_df(model) else Inf

  if (mean_scale) {
    X <- .orm_x(model, grid)
    y_range <- range(as.numeric(model$yunique))
    transform <- function(m) pmin(pmax(m, y_range[1]), y_range[2])
    to_inference_scale <- identity
  } else {
    design <- .design(model, grid)
    X <- design$X
    transform <- if (scale == "response") .inverse_link(model) else identity
    to_inference_scale <- if (scale == "link") identity else if (is_orm) .orm_mean_inverse(model) else .link(model)
  }

  if (ci_method != "bootstrap") {
    if (mean_scale) {
      delta <- .orm_mean(model, X, gradient = TRUE)
      estimate <- delta$mean
      G <- delta$gradient
      V <- .orm_vcov_all(model)
    } else {
      estimate <- drop(X %*% design$beta)
      G <- X
      V <- design$vcov
    }
    se <- sqrt(pmax(rowSums((G %*% V) * G), 0))
    max_t <- function(rows) .max_t_draws(G[rows, , drop = FALSE], V, df = df, nsim = nsim)
    nboot <- NULL
  } else {
    coefs <- .bootstrap_coefs(model, nboot, data)
    if (mean_scale) {
      estimate <- .orm_mean(model, X)
      draws <- t(apply(coefs, 1, function(b) .orm_mean(model, X, coefs = b)))
    } else {
      estimate <- drop(X %*% design$beta)
      draws <- coefs[, names(design$beta), drop = FALSE] %*% t(X)
    }
    draws <- matrix(draws, nrow = nrow(coefs))
    se <- apply(draws, 2, stats::sd)
    max_t <- function(rows) {
      .max_t_bootstrap(draws[, rows, drop = FALSE], estimate[rows], se[rows])
    }
    nboot <- nrow(coefs)
  }

  list(estimate = estimate, se = se, df = df, transform = transform,
       to_inference_scale = to_inference_scale, max_t = max_t, nboot = nboot,
       vcov_parts = if (ci_method != "bootstrap") list(G = G, V = V))
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
  cat("Predictions over", toString(focal), "(", nrow(x), "grid points,", attr(x, "scale"), "scale )\n")
  held <- .held_text(x)
  if (length(held)) {
    cat(strwrap(paste("Held constant:", toString(held)), exdent = 2), sep = "\n")
  }
  crit <- attr(x, "critical_value")
  method <- if (identical(attr(x, "ci_method"), "bootstrap")) {
    paste0("bootstrap, ", attr(x, "nboot"), " refits")
  } else {
    "analytic"
  }
  cat(sprintf(
    "Critical values for %g%% intervals (%s): pointwise %.3f, simultaneous %s%s\n\n",
    100 * attr(x, "ci_level"), method, attr(x, "pointwise_critical_value"),
    paste(sprintf("%.3f", crit), collapse = ", "),
    if (length(crit) > 1) paste0(" (one per ", sub("faceted", "facet", attr(x, "ci_scope")), ")") else ""
  ))
  cols <- c(focal, "estimate", "conf.low", "conf.high", "band.low", "band.high")
  print(utils::head(as.data.frame(x)[intersect(cols, names(x))], n), digits = 4, row.names = FALSE)
  if (nrow(x) > n) cat("... with", nrow(x) - n, "more rows\n")
  invisible(x)
}


# "age = 60 (median)", one string per held variable
.held_text <- function(x, rules = TRUE) {
  held <- attr(x, "held")
  if (!length(held)) {
    return(character(0))
  }
  rule <- attr(x, "held_rule")
  vapply(names(held), function(name) {
    v <- held[[name]]
    value <- if (is.numeric(v)) format(signif(v, 3)) else as.character(v)
    how <- if (rules && !is.null(rule[name]) && !is.na(rule[name])) paste0(" (", rule[[name]], ")") else ""
    paste0(name, " = ", value, how)
  }, character(1), USE.NAMES = FALSE)
}


#' Covariance matrix of predictions
#'
#' The joint covariance of the predictions in a [term_predictions()] result,
#' on the scale its intervals are computed on (the scale of `std.error`).
#' Ported from Lukas Bolsinger's `vcov()` for parameterized predictions in
#' p-ggeffects.
#'
#' @param object A [term_predictions()] result with analytic intervals.
#' @param ... Not used.
#' @return A square matrix with one row and column per grid point.
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3), data = mtcars)
#' pred <- term_predictions(fit, "hp", n = 5)
#' all.equal(sqrt(diag(vcov(pred))), pred$std.error)
#' @export
vcov.term_predictions <- function(object, ...) {
  parts <- attr(object, "vcov_parts")
  if (is.null(parts)) {
    stop("The covariance is only stored for analytic intervals, not bootstrap ones.", call. = FALSE)
  }
  unname(parts$G %*% parts$V %*% t(parts$G))
}

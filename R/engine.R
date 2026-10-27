# Prediction engine.
#
# For each supported model class, .design() returns everything needed to
# compute predictions and their joint uncertainty on the link scale:
#   X     design matrix for `newdata` (rows = grid points, cols = coefficients)
#   beta  coefficient estimates
#   vcov  covariance matrix of `beta`
# so that the linear predictor is X %*% beta with covariance X V X'.
#
# Only stats-level tools (model.matrix, coef, vcov) are used, so termeffects
# does not depend on any third-party prediction package.

.design <- function(model, newdata) {
  UseMethod(".design")
}

.design.default <- function(model, newdata) {
  .check_supported(model)
}

.check_supported <- function(model) {
  if (!inherits(model, c("lm", "merMod", "coxph", "gam", "orm"))) {
    stop("Models of class `", class(model)[1], "` are not supported yet. ",
         "Supported: lm, glm, lmerMod, glmerMod, coxph, gam (mgcv), orm (rms).", call. = FALSE)
  }
  invisible(TRUE)
}

.design.lm <- function(model, newdata) {
  .check_no_offset(model)
  X <- .model_matrix(model, newdata)
  .drop_aliased(X, stats::coef(model), stats::vcov(model))
}

.design.merMod <- function(model, newdata) {
  .check_no_offset(model)
  X <- .model_matrix(model, newdata)
  .drop_aliased(X, lme4::fixef(model), as.matrix(stats::vcov(model)))
}

.design.coxph <- function(model, newdata) {
  specials <- attr(stats::terms(model), "specials")
  used <- names(Filter(Negate(is.null), specials))
  if (length(used)) {
    stop("Cox models with ", toString(paste0(used, "()")),
         " terms are not supported yet.", call. = FALSE)
  }
  .check_no_offset(model)
  X <- .model_matrix(model, newdata)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  # centre at the sample means, as predict.coxph(type = "lp") does, so the
  # estimate is the log hazard ratio relative to an average subject
  X <- sweep(X, 2, model$means[colnames(X)])
  .drop_aliased(X, stats::coef(model), stats::vcov(model))
}

.design.orm <- function(model, newdata) {
  if (!is.null(stats::getCall(model)$offset) || !is.null(attr(stats::terms(model), "offset"))) {
    stop("Models with offsets are not supported yet.", call. = FALSE)
  }
  # rms's linear predictor uses the reference ("middle") intercept, and
  # vcov.orm() returns the covariance of that intercept plus the slopes
  X <- cbind(1, stats::predict(model, newdata = newdata, type = "x"))
  slopes <- stats::coef(model)[-seq_len(model$non.slopes)]
  beta <- c(stats::coef(model)[model$interceptRef], slopes)
  vcov <- stats::vcov(model)
  colnames(X) <- names(beta)
  dimnames(vcov) <- list(names(beta), names(beta))
  list(X = X, beta = beta, vcov = vcov)
}

.design.gam <- function(model, newdata) {
  X <- stats::predict(model, newdata = newdata, type = "lpmatrix")
  if (any(attr(X, "model.offset") != 0)) {
    stop("Models with offsets are not supported yet.", call. = FALSE)
  }
  list(X = unclass(X)[, , drop = FALSE], beta = stats::coef(model), vcov = stats::vcov(model))
}


# inverse link: maps the linear predictor to the response scale
.inverse_link <- function(model) {
  if (inherits(model, "coxph")) {
    return(exp)
  }
  if (inherits(model, "orm")) {
    # expected value of the outcome; increasing in the linear predictor
    return(Hmisc::Mean(model))
  }
  if (inherits(model, c("glm", "gam", "glmerMod"))) {
    return(stats::family(model)$linkinv)
  }
  identity
}


# residual degrees of freedom for t-based intervals; Inf means normal-based
.residual_df <- function(model) {
  if (inherits(model, "gam") || inherits(model, "glm")) {
    # the dispersion is estimated for these families, so intervals use t
    if (stats::family(model)$family %in% c("gaussian", "Gamma", "inverse.gaussian") ||
        grepl("^quasi", stats::family(model)$family)) {
      return(stats::df.residual(model))
    }
    return(Inf)
  }
  if (inherits(model, "lm")) {
    return(stats::df.residual(model))
  }
  Inf
}


# terms of the fixed (population-level) part of the model
.fixed_terms <- function(model) {
  if (inherits(model, "merMod")) {
    return(stats::terms(model, fixed.only = TRUE))
  }
  if (inherits(model, "gam")) {
    return(stats::terms(stats::formula(model)))
  }
  stats::terms(model)
}


# fixed-effects model matrix for new data, reusing the fitted model's spline
# knots (predvars), factor levels and contrasts
.model_matrix <- function(model, newdata) {
  tt <- stats::delete.response(.fixed_terms(model))
  original <- stats::model.matrix(model)
  # factor levels are keyed by model-frame column name, e.g. "factor(am)"
  xlev <- lapply(Filter(is.factor, stats::model.frame(model)), levels)
  frame_columns <- vapply(as.list(attr(tt, "variables"))[-1], deparse1, character(1))
  mf <- stats::model.frame(tt, newdata, xlev = xlev[names(xlev) %in% frame_columns],
                           na.action = stats::na.pass)
  stats::model.matrix(tt, mf, contrasts.arg = attr(original, "contrasts"))
}


# keep only estimable coefficients and line up X, beta and vcov by name
.drop_aliased <- function(X, beta, vcov) {
  keep <- names(beta)[!is.na(beta)]
  missing <- setdiff(keep, colnames(X))
  if (length(missing)) {
    stop("Could not build the design matrix for: ", toString(missing), call. = FALSE)
  }
  list(
    X = X[, keep, drop = FALSE],
    beta = beta[keep],
    vcov = as.matrix(vcov)[keep, keep, drop = FALSE]
  )
}


# variables stored as numbers but turned into factors in the formula, e.g.
# `am` in `y ~ factor(am)`; these are treated as categorical in the grid
.factor_in_formula <- function(model) {
  frame <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  columns <- names(Filter(is.factor, frame))
  vars <- lapply(columns, function(column) {
    v <- all.vars(str2lang(column))
    if (length(v) == 1 && v != column) v
  })
  intersect(unlist(vars), .predictor_names(model))
}


.check_no_offset <- function(model) {
  tt <- .fixed_terms(model)
  if (!is.null(attr(tt, "offset")) || !is.null(stats::getCall(model)$offset)) {
    stop("Models with offsets are not supported yet.", call. = FALSE)
  }
}

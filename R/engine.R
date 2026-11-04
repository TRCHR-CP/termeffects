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
  if (!inherits(model, c("lm", "merMod", "coxph", "gam", "orm", "rms", "lme", "logistf"))) {
    stop("Models of class `", class(model)[1], "` are not supported yet. ",
         "Supported: lm, glm (including brglm2 and geepack::geeglm), lmerMod, glmerMod, ",
         "nlme::lme, coxph, gam (mgcv), logistf, and rms's orm, lrm, cph and ols.", call. = FALSE)
  }
  invisible(TRUE)
}

# lm, glm, and models built on them (e.g. geepack::geeglm, brglm2 fits)
.design.lm <- function(model, newdata) {
  X <- .model_matrix(model, newdata)
  .with_offset(.drop_aliased(X, stats::coef(model), stats::vcov(model)), model, newdata)
}

.design.merMod <- function(model, newdata) {
  X <- .model_matrix(model, newdata)
  .with_offset(.drop_aliased(X, lme4::fixef(model), as.matrix(stats::vcov(model))), model, newdata)
}

# nlme::lme: fixed effects only (population-level prediction)
.design.lme <- function(model, newdata) {
  .check_no_offset(model)
  X <- .model_matrix(model, newdata)
  .drop_aliased(X, nlme::fixef(model), as.matrix(stats::vcov(model)))
}

# Firth's penalized logistic regression (logistf). The covariance is that of
# the penalized fit, so the intervals are Wald-type, not logistf's default
# profile penalized likelihood intervals.
.design.logistf <- function(model, newdata) {
  .check_no_offset(model)
  X <- .model_matrix(model, newdata)
  .drop_aliased(X, stats::coef(model), stats::vcov(model))
}

.design.coxph <- function(model, newdata) {
  tt <- stats::terms(model)
  specials <- attr(tt, "specials")
  # strata() only changes the baseline hazard, so hazard ratios are unaffected
  used <- setdiff(names(Filter(Negate(is.null), specials)), "strata")
  if (length(used)) {
    stop("Cox models with ", toString(paste0(used, "()")),
         " terms are not supported yet.", call. = FALSE)
  }
  if (any(grepl("strata(", attr(tt, "term.labels")[attr(tt, "order") > 1], fixed = TRUE))) {
    stop("Cox models with interactions with strata() are not supported yet.", call. = FALSE)
  }
  .check_no_offset(model)
  X <- .model_matrix(model, newdata)
  X <- X[, colnames(X) %in% names(stats::coef(model)), drop = FALSE]
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
  X <- cbind(1, .rms_x(model, newdata))
  slopes <- stats::coef(model)[-seq_len(model$non.slopes)]
  beta <- c(stats::coef(model)[model$interceptRef], slopes)
  vcov <- stats::vcov(model)
  colnames(X) <- names(beta)
  dimnames(vcov) <- list(names(beta), names(beta))
  list(X = X, beta = beta, vcov = vcov)
}

# rms's lrm (binary outcome), ols and cph; orm has its own method. rms builds
# the design matrix itself (rcs() and friends keep their knots there).
.design.rms <- function(model, newdata) {
  .check_no_offset(model)
  beta <- stats::coef(model)
  X <- .rms_x(model, newdata)
  if (inherits(model, "cph")) {
    # log hazard ratio relative to an average subject, as for coxph
    X <- sweep(X, 2, model$means)
  } else {
    if (isTRUE(model$non.slopes > 1)) {
      stop("lrm models with an ordinal outcome are not supported; fit the model with rms::orm().",
           call. = FALSE)
    }
    X <- cbind(1, X)
  }
  if (ncol(X) != length(beta)) {
    stop("Could not build the design matrix of this ", class(model)[1], " model.", call. = FALSE)
  }
  colnames(X) <- names(beta)
  list(X = X, beta = beta, vcov = as.matrix(stats::vcov(model))[names(beta), names(beta)])
}

# design matrix of an rms model for new data (without intercepts). For orm,
# rms inverts the information matrix here, which fails for badly scaled models.
.rms_x <- function(model, newdata) {
  tryCatch(
    stats::predict(model, newdata = newdata, type = "x"),
    error = function(e) {
      stop("rms could not invert the model's information matrix (", conditionMessage(e), "). ",
           "This usually means a predictor is on a very large scale (e.g. counts in the ",
           "thousands): rescale it or use its log, and refit.", call. = FALSE)
    }
  )
}


# Mean outcome of an orm model and its gradient with respect to all
# coefficients (every intercept and the slopes), for the delta method. This is
# the computation of rms's Mean(fit)(X = , conf.int = ), which also uses the
# full covariance matrix; it is redone here to get the joint covariance.
.orm_mean <- function(model, X, coefs = stats::coef(model), gradient = FALSE) {
  ns <- model$non.slopes
  alpha <- coefs[seq_len(ns)]
  xb <- drop(X %*% coefs[-seq_len(ns)])
  fam <- model$famfunctions
  if (is.null(fam)) {
    fam <- rms::probabilityFamilies$logistic
  }
  cumprob <- eval(fam[[1]])
  y <- as.numeric(model$yunique)
  eta <- outer(xb, alpha, "+")
  # P(Y >= y_(j+1)), so the mean is y_1 + sum_j P_j (y_(j+1) - y_j)
  P <- matrix(cumprob(eta), ncol = ns)
  mean <- y[1] + drop(P %*% diff(y))
  if (!gradient) {
    return(mean)
  }
  density <- eval(fam[[5]])
  d_alpha <- sweep(matrix(density(eta), ncol = ns), 2, diff(y), "*")
  list(mean = mean, gradient = cbind(d_alpha, rowSums(d_alpha) * X))
}


# inverse of Hmisc::Mean(model): the linear predictor giving a mean outcome
# (Mean is increasing in it); means outside the outcome range give -Inf/Inf
.orm_mean_inverse <- function(model) {
  mean_fun <- Hmisc::Mean(model)
  ref <- stats::coef(model)[model$interceptRef]
  alpha <- stats::coef(model)[seq_len(model$non.slopes)]
  lp <- seq(min(ref - alpha) - 20, max(ref - alpha) + 20, length.out = 4001)
  m <- mean_fun(lp)
  inverse <- stats::approxfun(m, lp, ties = mean)
  y_range <- range(as.numeric(model$yunique))
  function(mean) {
    out <- inverse(mean)
    out[mean <= y_range[1]] <- -Inf
    out[mean >= y_range[2]] <- Inf
    out
  }
}


# Covariance of all orm coefficients (every intercept and the slopes). Same as
# vcov(model, intercepts = "all"), but the information matrix is rescaled to
# unit diagonal before inverting, so it also works when predictors on very
# different scales (e.g. platelets per mL) make the raw matrix numerically
# singular.
.orm_vcov_all <- function(model) {
  info <- as.matrix(rms::infoMxop(model$info.matrix))
  d <- 1 / sqrt(diag(info))
  V <- solve(info * outer(d, d)) * outer(d, d)
  dimnames(V) <- list(names(stats::coef(model)), names(stats::coef(model)))
  V
}


# Coefficients refitted on bootstrap samples (rows = samples). rms fits use
# rms::bootcov(), which also handles orm intercepts missing from a resample;
# other models are refitted on resampled rows of their data.
.bootstrap_coefs <- function(model, nboot, data = NULL) {
  if (.is_clustered(model)) {
    stop("Bootstrap intervals are not available for models of clustered or repeated data ",
         "(mixed models, GEE, Cox models with (start, stop] data or cluster()), where ",
         "clusters rather than rows would have to be resampled. Use ci_method = \"analytic\".",
         call. = FALSE)
  }
  if (inherits(model, "rms")) {
    if (is.null(model$x) || is.null(model$y)) {
      stop("For bootstrap intervals, refit the rms model with `x = TRUE, y = TRUE`.", call. = FALSE)
    }
    utils::capture.output(boot <- suppressMessages(suppressWarnings(rms::bootcov(model, B = nboot))))
    coefs <- boot$boot.Coef
    return(coefs[stats::complete.cases(coefs), , drop = FALSE])
  }

  full <- .model_data(model, data, all_columns = TRUE)
  used <- intersect(all.vars(stats::formula(model)), names(full))
  full <- full[stats::complete.cases(full[used]), , drop = FALSE]
  call <- stats::getCall(model)
  env <- environment(stats::formula(model))
  original <- stats::coef(model)
  coefs <- lapply(seq_len(nboot), function(b) {
    call$data <- full[sample.int(nrow(full), replace = TRUE), , drop = FALSE]
    fit <- tryCatch(suppressWarnings(eval(call, env)), error = function(e) NULL)
    if (is.null(fit)) {
      return(NULL)
    }
    stats::coef(fit)[names(original)]
  })
  coefs <- do.call(rbind, coefs)
  if (is.null(coefs)) {
    stop("All bootstrap refits failed.", call. = FALSE)
  }
  # refits in which a coefficient could not be estimated are dropped
  coefs[stats::complete.cases(coefs[, !is.na(original), drop = FALSE]), , drop = FALSE]
}


# models whose rows are not independent, so rows cannot be bootstrapped
.is_clustered <- function(model) {
  inherits(model, c("merMod", "lme", "geeglm")) ||
    (inherits(model, "coxph") && (identical(ncol(model$y), 3L) || !is.null(model$naive.var)))
}


# family of a generalized linear model (also for lrm and logistf, which are
# logistic models without a family() method); NULL for other models
.family_of <- function(model) {
  if (inherits(model, c("lrm", "logistf"))) {
    return(stats::binomial())
  }
  if (inherits(model, c("glm", "gam", "glmerMod"))) {
    return(stats::family(model))
  }
  NULL
}


# link function: maps the response scale back to the linear predictor
.link <- function(model) {
  if (inherits(model, "coxph")) {
    return(log)
  }
  fam <- .family_of(model)
  if (!is.null(fam)) {
    return(fam$linkfun)
  }
  identity
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
  fam <- .family_of(model)
  if (!is.null(fam)) {
    return(fam$linkinv)
  }
  identity
}


# residual degrees of freedom for t-based intervals; Inf means normal-based
.residual_df <- function(model) {
  if (inherits(model, c("geeglm", "lrm", "logistf"))) {
    # robust (sandwich) or likelihood-based covariance: normal-based intervals
    return(Inf)
  }
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
  if (inherits(model, "lme")) {
    return(model$terms)
  }
  if (inherits(model, "coxph") && !inherits(model, "cph")) {
    # strata() terms only change the baseline hazard: leave them out, so they
    # are neither in the design matrix nor held covariates
    tt <- stats::terms(model)
    strata <- grep("strata(", attr(tt, "term.labels"), fixed = TRUE)
    if (length(strata)) {
      tt <- stats::drop.terms(tt, strata, keep.response = TRUE)
    }
    return(tt)
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
  frame <- .original_frame(model)
  # factor levels are keyed by model-frame column name, e.g. "factor(am)"
  xlev <- lapply(Filter(is.factor, frame), levels)
  frame_columns <- vapply(as.list(attr(tt, "variables"))[-1], deparse1, character(1))
  mf <- stats::model.frame(tt, newdata, xlev = xlev[names(xlev) %in% frame_columns],
                           na.action = stats::na.pass)
  X <- stats::model.matrix(tt, mf, contrasts.arg = .original_contrasts(model, frame))
  attr(X, "offset") <- stats::model.offset(mf)
  X
}

# model frame of the fitted model (nlme::lme keeps none, so it is rebuilt)
.original_frame <- function(model) {
  if (inherits(model, "lme")) {
    return(stats::model.frame(model$terms, nlme::getData(model)))
  }
  stats::model.frame(model)
}

.original_contrasts <- function(model, frame) {
  if (inherits(model, "lme")) {
    return(model$contrasts)
  }
  tryCatch(attr(stats::model.matrix(model), "contrasts"),
           error = function(e) attr(stats::model.matrix(.fixed_terms(model), frame), "contrasts"))
}

# offset for new data: offset() terms in the formula plus an `offset =`
# argument of the call, evaluated at the grid's values (offset variables are
# held like other covariates unless set with `at`)
.with_offset <- function(design, model, newdata) {
  offset <- attr(.model_matrix(model, newdata), "offset")
  call_offset <- stats::getCall(model)$offset
  if (!is.null(call_offset)) {
    extra <- eval(call_offset, newdata, environment(stats::formula(model)))
    offset <- if (is.null(offset)) extra else offset + extra
  }
  design$offset <- offset
  design
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

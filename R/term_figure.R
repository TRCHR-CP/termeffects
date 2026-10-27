#' A figure with its plots, covariate list and caption
#'
#' Bundles everything needed to present a non-linear interaction:
#'
#' * `$heatmap`: the prediction surface ([term_heatmap()]),
#' * `$curves`: cuts through it with their bands ([term_curves()]), if
#'   `cuts` are given,
#' * `$uncertainty`: a heatmap of the band width, showing where the surface
#'   is estimated precisely,
#' * `$covariates`: a data frame of every covariate held constant, its value,
#'   and why that value (median, most common, set by user),
#' * `$caption`: a figure legend in words: what is predicted, from which
#'   model, how the two focal variables enter the model, at what values the
#'   other covariates are held and which intervals are shown,
#' * `$model`: a one-line description of the model.
#'
#' The plots themselves only name the held covariates when there are a few;
#' with many, the caption lists them.
#'
#' @param x A [term_predictions()] result with two continuous focal
#'   variables (a surface).
#' @param cuts Optional named list of values of one focal variable at which to
#'   cut the surface, e.g. `list(sodium = c(130, 137, 142))`. The curves are
#'   computed with the same settings as `x`.
#' @param band Band shown in `$curves` (see [term_curves()]).
#' @param ... Passed to [term_heatmap()], e.g. `show_data = TRUE`.
#'
#' @return A list of class `term_figure`. Printing it draws the plots and
#'   prints the caption.
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * splines::ns(wt, 3) + am + qsec, data = mtcars)
#' surface <- term_predictions(fit, c("hp", "wt"), n = 30)
#' fig <- term_figure(surface, cuts = list(wt = c(2.5, 3.5)))
#' fig$heatmap
#' fig$covariates
#' cat(fig$caption)
#' @export
term_figure <- function(x, cuts = NULL, band = "simultaneous", ...) {
  if (!.is_surface(x) && !(length(attr(x, "focal")) >= 2 &&
                           is.numeric(x[[attr(x, "focal")[2]]]))) {
    stop("`x` must be a surface: two continuous focal variables.", call. = FALSE)
  }
  focal <- attr(x, "focal")
  curves <- NULL
  if (!is.null(cuts)) {
    if (length(cuts) != 1 || !names(cuts) %in% focal[1:2]) {
      stop("`cuts` must give values of one of ", toString(focal[1:2]), ".", call. = FALSE)
    }
    along <- setdiff(focal[1:2], names(cuts))
    cut_values <- list(sort(unique(as.data.frame(x)[[along]])))
    names(cut_values) <- along
    cut_values[names(cuts)] <- cuts
    curves <- .recompute(x, terms = c(along, names(cuts), focal[-(1:2)]), values = cut_values)
  }

  covariates <- .covariate_table(x)
  structure(list(
    heatmap = term_heatmap(x, cuts = cuts, ...),
    curves = if (!is.null(curves)) term_curves(curves, band = band),
    uncertainty = term_heatmap(x, fill = "band_width", contours = FALSE),
    covariates = covariates,
    caption = .figure_caption(x, curves, band),
    model = .model_description(attr(x, "model"))
  ), class = "term_figure")
}


#' @export
print.term_figure <- function(x, ...) {
  print(x$heatmap)
  if (!is.null(x$curves)) {
    print(x$curves)
  }
  cat(strwrap(x$caption), sep = "\n")
  invisible(x)
}


# the same predictions as `x` on another grid
.recompute <- function(x, terms, values) {
  held <- attr(x, "held")
  out <- term_predictions(attr(x, "model"), terms, values = values, at = held,
                          ci_level = attr(x, "ci_level"), scale = attr(x, "scale"),
                          ci_method = attr(x, "ci_method"), ci_scope = attr(x, "ci_scope"))
  attr(out, "held_rule") <- attr(x, "held_rule")
  out
}


.covariate_table <- function(x) {
  held <- attr(x, "held")
  if (!length(held)) {
    return(data.frame(covariate = character(0), value = character(0), chosen_as = character(0)))
  }
  rule <- attr(x, "held_rule")
  data.frame(
    covariate = names(held),
    value = vapply(held, function(v) if (is.numeric(v)) format(signif(v, 4)) else as.character(v),
                   character(1)),
    chosen_as = unname(ifelse(is.na(rule[names(held)]), "", rule[names(held)])),
    row.names = NULL
  )
}


.model_description <- function(model) {
  pm <- is.parameterized_model(model)
  model <- original_model(model)
  type <- if (inherits(model, "orm")) {
    "proportional odds model (rms::orm)"
  } else if (inherits(model, "coxph")) {
    "Cox proportional hazards model (survival::coxph)"
  } else if (inherits(model, "gam")) {
    paste0("generalized additive model (mgcv::gam, ", stats::family(model)$family, ")")
  } else if (inherits(model, "glmerMod")) {
    paste0("generalized linear mixed model (lme4::glmer, ", stats::family(model)$family, ")")
  } else if (inherits(model, "lmerMod")) {
    "linear mixed model (lme4::lmer)"
  } else if (inherits(model, "glm")) {
    fam <- stats::family(model)
    if (fam$family == "binomial" && fam$link == "logit") "logistic regression (glm)" else
      paste0("generalized linear model (glm, ", fam$family, " family, ", fam$link, " link)")
  } else {
    "linear regression"
  }
  n <- tryCatch(stats::nobs(model), error = function(e) NA)
  paste0(if (pm) "parameterized " else "", type,
         if (!is.na(n)) paste0(", n = ", n), ": ",
         paste(deparse(stats::formula(model), width.cutoff = 500), collapse = " "))
}


# how the focal variables enter the model, from the model's term labels
.focal_terms_text <- function(model, focal) {
  labels <- attr(stats::delete.response(.fixed_terms(model)), "term.labels")
  vars <- lapply(labels, function(l) all.vars(str2lang(l)))
  main <- vapply(focal, function(v) {
    own <- labels[vapply(vars, function(u) identical(u, v), logical(1))]
    if (!length(own)) return(paste(v, "(only in interactions)"))
    if (identical(own, v)) paste(v, "(linear)") else paste0(v, " as ", toString(own))
  }, character(1))
  both <- labels[vapply(vars, function(u) all(focal %in% u), logical(1))]
  interaction <- if (length(both)) {
    paste0(" They interact in the model (", toString(both), ").")
  } else {
    " The model has no interaction between them, so the surface is additive."
  }
  paste0("The focal variables enter the model as ", paste(main, collapse = " and "), ".", interaction)
}


.figure_caption <- function(x, curves, band) {
  focal <- attr(x, "focal")
  model <- attr(x, "model")
  level <- paste0(100 * attr(x, "ci_level"), "%")
  held <- .held_text(x)
  method <- if (identical(attr(x, "ci_method"), "bootstrap")) {
    paste0(" from ", attr(x, "nboot"), " bootstrap refits")
  } else {
    ""
  }
  crit <- attr(x, if (is.null(curves)) "critical_value" else "critical_value")
  text <- c(
    paste0(attr(x, "outcome_label"), " over ", focal[1], " and ", focal[2],
           if (length(focal) >= 3) paste0(", in panels by ", toString(focal[-(1:2)])),
           ", from a ", .model_description(model), "."),
    if (is.parameterized_model(model)) {
      p <- parameterizations(model)
      paste0("Predictors computed from the plotted variables: ",
             paste0(names(p), " = ", vapply(p, function(f) paste(deparse(.unwrap_body(body(f))), collapse = ""),
                                            character(1)), collapse = "; "), ".")
    },
    .focal_terms_text(original_model(model), intersect(focal[1:2], .predictor_names(original_model(model)))),
    if (length(held)) paste0("All other covariates are held constant: ", toString(held), ".")
    else "There are no other covariates.",
    if (!is.null(curves)) {
      paste0("Curves: cuts through the surface at ", attr(curves, "focal")[2], " = ",
             toString(format(unique(as.data.frame(curves)[[attr(curves, "focal")[2]]]), trim = TRUE)), ". ",
             .band_sentence(curves, band, level, method))
    } else {
      paste0("The ", level, " simultaneous band", method, " over this surface has critical value ",
             paste(sprintf("%.2f", crit), collapse = ", "), " (pointwise: ",
             sprintf("%.2f", attr(x, "pointwise_critical_value")), ").")
    }
  )
  paste(text, collapse = " ")
}

.band_sentence <- function(curves, band, level, method) {
  crit <- paste(sprintf("%.2f", attr(curves, "critical_value")), collapse = ", ")
  switch(band,
    simultaneous = paste0("Shaded: ", level, " simultaneous confidence band", method,
                          ", covering all curves at once (critical value ", crit, ")."),
    pointwise = paste0("Shaded: ", level, " pointwise confidence intervals", method, "."),
    both = paste0("Shaded: ", level, " simultaneous band", method, " (critical value ", crit,
                  "); dashed: ", level, " pointwise intervals."),
    none = ""
  )
}

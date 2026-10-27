#' Build a prediction grid for focal terms
#'
#' Creates a data frame in which the focal variables vary over a regular grid
#' and every other predictor in the model is held constant. This is the grid
#' a heatmap (two focal terms) or a set of curves (one continuous focal term
#' plus a few values of a second one) is drawn on.
#'
#' Held-constant values: numeric predictors are set to their mean, factor and
#' character predictors to their reference (first) level, and logical
#' predictors to `FALSE`. Use `at` to choose other values.
#'
#' @param model A fitted model (`lm`, `glm`, `lmerMod`, `glmerMod`, `coxph`,
#'   or `mgcv::gam`).
#' @param terms Character vector with the names of the focal variables.
#' @param values Optional named list giving the values to use for some or all
#'   focal variables, e.g. `list(age = c(50, 60, 70))`.
#' @param at Optional named list giving the values at which non-focal
#'   variables are held, e.g. `list(sex = "female")`.
#' @param n Number of grid points for each continuous focal variable without
#'   explicit `values`. Points are spread evenly over the observed range.
#' @param data The data the model was fitted on. Only needed if it cannot be
#'   recovered from the model call.
#'
#' @return A data frame with class `term_grid`. Attribute `focal` holds the
#'   focal variable names and attribute `held` the values of the other
#'   predictors.
#'
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * splines::ns(wt, 3) + am, data = mtcars)
#' grid <- term_grid(fit, c("hp", "wt"), n = 5)
#' head(grid)
#' attr(grid, "held")
#' @export
term_grid <- function(model, terms, values = NULL, at = NULL, n = 30, data = NULL) {
  .check_supported(model)
  data <- .model_data(model, data)
  predictors <- .predictor_names(model)

  categorical <- .factor_in_formula(model)
  data[categorical] <- lapply(data[categorical], function(x) factor(x, levels = sort(unique(x))))

  unknown <- setdiff(c(terms, names(values), names(at)), names(data))
  if (length(unknown)) {
    stop("Variable(s) not found in the model data: ", toString(unknown), call. = FALSE)
  }
  if (length(setdiff(names(values), terms))) {
    stop("`values` may only contain focal variables named in `terms`.", call. = FALSE)
  }
  if (length(intersect(names(at), terms))) {
    stop("Focal variables cannot also be given in `at`; use `values` instead.", call. = FALSE)
  }

  focal_values <- lapply(terms, function(term) {
    if (!is.null(values[[term]])) {
      return(.match_type(values[[term]], data[[term]]))
    }
    .focal_values(data[[term]], n)
  })
  names(focal_values) <- terms

  held_names <- setdiff(predictors, terms)
  held <- lapply(held_names, function(name) {
    if (!is.null(at[[name]])) {
      return(.match_type(at[[name]], data[[name]]))
    }
    .typical_value(data[[name]])
  })
  names(held) <- held_names

  grid <- expand.grid(c(focal_values, held), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  attr(grid, "focal") <- terms
  attr(grid, "held") <- held
  class(grid) <- c("term_grid", "data.frame")
  grid
}


# values over which a focal variable varies
.focal_values <- function(x, n) {
  if (is.numeric(x)) {
    rng <- range(x, na.rm = TRUE)
    return(seq(rng[1], rng[2], length.out = n))
  }
  if (is.logical(x)) {
    return(c(FALSE, TRUE))
  }
  x <- as.factor(x)
  factor(levels(x), levels = levels(x))
}


# value at which a non-focal variable is held constant
.typical_value <- function(x) {
  if (is.numeric(x)) {
    return(mean(x, na.rm = TRUE))
  }
  if (is.logical(x)) {
    return(FALSE)
  }
  x <- as.factor(x)
  factor(levels(x)[1], levels = levels(x))
}


# give user-supplied values the same type (and factor levels) as the data
.match_type <- function(value, x) {
  if (is.factor(x) || is.character(x)) {
    lev <- levels(as.factor(x))
    bad <- setdiff(as.character(value), lev)
    if (length(bad)) {
      stop("Value(s) ", toString(bad), " are not levels of the variable. ",
           "Levels are: ", toString(lev), call. = FALSE)
    }
    return(factor(value, levels = lev))
  }
  value
}


# the data the model was fitted on, restricted to complete predictor rows
.model_data <- function(model, data = NULL) {
  if (is.null(data)) {
    call_data <- stats::getCall(model)$data
    if (!is.null(call_data)) {
      data <- tryCatch(
        eval(call_data, environment(stats::formula(model))),
        error = function(e) NULL
      )
    }
  }
  vars <- .predictor_names(model)
  if (!is.data.frame(data)) {
    # some models (e.g. mgcv::gam) keep the raw variables in their model frame
    frame <- tryCatch(stats::model.frame(model), error = function(e) NULL)
    if (is.data.frame(frame) && all(vars %in% names(frame))) {
      data <- frame
    }
  }
  if (!is.data.frame(data)) {
    stop("Could not recover the model data. Please pass it via `data = `.", call. = FALSE)
  }
  vars <- intersect(vars, names(data))
  data <- as.data.frame(data)[vars]
  data[stats::complete.cases(data), , drop = FALSE]
}


# names of the variables on the right-hand side of the model formula
.predictor_names <- function(model) {
  tt <- .fixed_terms(model)
  all.vars(stats::delete.response(tt))
}

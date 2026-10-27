#' Build a prediction grid for focal terms
#'
#' Creates a data frame in which the focal variables vary over a regular grid
#' and every other predictor in the model is held constant. This is the grid
#' a heatmap (two focal terms) or a set of curves (one continuous focal term
#' plus a few values of a second one) is drawn on.
#'
#' Held-constant values: numeric predictors are set to their median, and
#' factor, character and logical predictors to their most common value (the
#' mode; the first level if tied). Use `at` to choose other values.
#'
#' For a [parameterized model][parameterize()], the terms can be its new
#' variables (e.g. `"egfr"`) as well as untransformed predictors. Without
#' `values`, a new variable spans a range centred at zero within which the
#' parameterized predictors stay inside their observed range (the search
#' written for p-ggeffects). Usually it is clearer to give the range with
#' `values`.
#'
#' @param model A fitted model (`lm`, `glm`, `rms::orm`, `lmerMod`,
#'   `glmerMod`, `coxph`, or `mgcv::gam`) or a [parameterized model][parameterize()].
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
#'   focal variable names, attribute `held` the values of the other
#'   variables and attribute `held_rule` how each was chosen (`"median"`,
#'   `"most common"` or `"set by user"`).
#'
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) * splines::ns(wt, 3) + am, data = mtcars)
#' grid <- term_grid(fit, c("hp", "wt"), n = 5)
#' head(grid)
#' attr(grid, "held")
#' @export
term_grid <- function(model, terms, values = NULL, at = NULL, n = 30, data = NULL) {
  pmodel <- NULL
  if (is.parameterized_model(model)) {
    pmodel <- model
    model <- original_model(model)
  }
  .check_supported(model)
  data <- .model_data(model, data)
  predictors <- .predictor_names(model)

  categorical <- .factor_in_formula(model)
  data[categorical] <- lapply(data[categorical], function(x) factor(x, levels = sort(unique(x))))

  # variables that make up the grid: predictors, except those computed by a
  # parameterization, plus the variables they are computed from
  replaced <- character(0)
  inputs <- character(0)
  if (!is.null(pmodel)) {
    replaced <- names(parameterizations(pmodel))
    inputs <- parameter_variables(pmodel)
  }
  parameters <- setdiff(inputs, names(data))
  variables <- union(setdiff(predictors, replaced), inputs)

  computed <- intersect(c(terms, names(values), names(at)), setdiff(replaced, inputs))
  if (length(computed)) {
    stop(toString(computed), " is computed by the parameterization. ",
         "Use the parameterization's variables instead: ", toString(inputs), ".",
         call. = FALSE)
  }
  unknown <- setdiff(c(terms, names(values), names(at)), c(names(data), parameters))
  if (length(unknown)) {
    stop("Variable(s) not found in the model data: ", toString(unknown), call. = FALSE)
  }
  if (length(setdiff(names(values), terms))) {
    stop("`values` may only contain focal variables named in `terms`.", call. = FALSE)
  }
  if (length(intersect(names(at), terms))) {
    stop("Focal variables cannot also be given in `at`; use `values` instead.", call. = FALSE)
  }

  held <- list()
  rule <- character(0)
  for (name in setdiff(variables, terms)) {
    if (!is.null(at[[name]])) {
      held[[name]] <- if (name %in% names(data)) .match_type(at[[name]], data[[name]]) else at[[name]]
      rule[[name]] <- "set by user"
    } else if (name %in% names(data)) {
      held[[name]] <- .typical_value(data[[name]])
      rule[[name]] <- .typical_rule(data[[name]])
    } else {
      stop("`", name, "` is a variable of the parameterization and not in the data, ",
           "so it has no typical value. Make it a focal term or give its value with ",
           "`at = list(", name, " = ...)`.", call. = FALSE)
    }
  }

  free <- setdiff(intersect(terms, parameters), names(values))
  ranges <- if (length(free)) .parameter_ranges(pmodel, free, held, data) else list()

  focal_values <- lapply(terms, function(term) {
    if (!is.null(values[[term]])) {
      return(if (term %in% names(data)) .match_type(values[[term]], data[[term]]) else values[[term]])
    }
    if (term %in% names(data)) {
      return(.focal_values(data[[term]], n))
    }
    if (is.null(ranges[[term]])) {
      stop("Could not find a range for `", term, "`. Give it with `values = `.", call. = FALSE)
    }
    seq(ranges[[term]][1], ranges[[term]][2], length.out = n)
  })
  names(focal_values) <- terms

  grid <- expand.grid(c(focal_values, held), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  attr(grid, "focal") <- terms
  attr(grid, "held") <- held
  attr(grid, "held_rule") <- rule
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
    return(stats::median(x, na.rm = TRUE))
  }
  # categories: the most common one (the first level if tied)
  counts <- table(x)
  mode <- names(counts)[which.max(counts)]
  if (is.logical(x)) {
    return(as.logical(mode))
  }
  x <- as.factor(x)
  factor(mode, levels = levels(x))
}

.typical_rule <- function(x) {
  if (is.numeric(x)) "median" else "most common"
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


# the data the model was fitted on; by default only complete rows of the
# predictors, with `all_columns = TRUE` the whole data frame (for refitting)
.model_data <- function(model, data = NULL, all_columns = FALSE) {
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
  data <- as.data.frame(data)
  if (all_columns) {
    return(data)
  }
  vars <- intersect(vars, names(data))
  data <- data[vars]
  data[stats::complete.cases(data), , drop = FALSE]
}


# names of the variables on the right-hand side of the model formula
.predictor_names <- function(model) {
  tt <- .fixed_terms(model)
  all.vars(stats::delete.response(tt))
}

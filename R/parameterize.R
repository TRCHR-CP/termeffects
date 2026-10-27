# Parameterized models: written by Lukas Bolsinger for p-ggeffects
# (R/parameterize.R). Ported with only these changes: insight helpers replaced
# by base R, misspelt parameterization names are an error instead of being
# dropped silently, and the examples use termeffects functions.

.abort <- function(...) {
  stop(paste0(...), call. = FALSE)
}

#' Parameterized model
#'
#' A `parameterized_model` represents a pair consisting of a model and a
#' collection of parameterizations of its continuous predictor variables.
#' Viewing the model as a function, a parameterized_model defines a new
#' model by pulling back the original model along the parameterization.
#'
#' @details
#' The parameterizations are functions that define new model variables in
#' terms of the original model variables. Only variables that occur as
#' predictors in the original model are parameterized. The variables on
#' which the parameterizations depend are retained as the parameters of the
#' parameterized_model.
#'
#' Predictions from a `parameterized_model` are obtained by first applying
#' its parameterizations to newdata and then passing the resulting data to
#' the prediction method of the original model. Consequently, newdata
#' should contain the variables required by the parameterizations rather than
#' the parameterized variables themselves.
#'
#' The model terms of a `parameterized_model` are those of the original
#' model pulled pack by the parameterization. In particular, the model
#' coefficients are retained. Parameterization changes how the corresponding
#' predictor variables are represented for prediction, but does not modify the
#' underlying model object or its terms.
#'
#' A `parameterized_model` retains a history of parameterizations applied to
#' it. When additional parameterizations are applied to an existing
#' parameterized_model, they are composed with the existing parameterizations
#' where appropriate, and the new parameterizations are added to the history.
#' This allows the sequence of parameterizations to be inspected or reconstructed.
#'
#' @seealso
#' [parameterize()], [parameterize_data()], [parameterizations()],
#' [parameter_variables()], [original_model()]
#'
#' @name parameterized_model-class
#' @docType class
NULL

#' Parameterize a model
#'
#' Adds parameterizations to the continuous predictor variables of a model.
#' The resulting object is a parameterized_model, which represents a model
#' together with a collection of parameterizations of its predictor
#' variables.
#'
#' Parameterizations are supplied as a named list. Each element defines a
#' parameterized predictor, with the name of the element giving the
#' corresponding predictor variable. Parameterizations may be supplied as
#' functions, expressions, or character strings and are converted to
#' functions internally.
#'
#' If model is already a parameterized_model, the new parameterizations
#' are composed with the existing parameterizations. The new
#' parameterizations are also appended to the parameterization history.
#'
#' @param model A model object, or an existing parameterized_model.
#' @param parameterizations A named list of parameterizations. The names
#' identify the continuous predictor variables to be parameterized.
#' Parameterizations may be functions, expressions, or character strings.
#'
#' @return
#' A parameterized_model object. If parameterizations is NULL, model
#' is returned unchanged, and does not inherit `parameterized_model`.
#'
#' @seealso
#' [parameterized_model()], [parameterize_data()],
#' [parameterizations()], [parameter_variables()], [original_model()]
#'
#' @examples
#' # Parameterize a predictor as a function of two underlying variables
#'
#' model <- lm(y ~ x,
#'             data = data.frame(y = 1:10, x = 10:1))
#'
#' pmodel <- parameterize(model,
#'                        parameterizations = list(x = \(a, b) a + b))
#'
#' ppmodel <- parameterize(pmodel,
#'                         parameterizations = list(a = "sin(t)"))
#'
#' predict(ppmodel,
#'         newdata = data.frame(t = -5:5, b = 0:10))
#'
#'
#' @export
parameterize <- function(model, parameterizations = NULL) {
  if(is.null(parameterizations)){
    return(model)
  }

  if (!is.list(parameterizations)) {
    .abort("`parameterizations` must be a named list.")
  }

  if (is.null(names(parameterizations))) {
    .abort("`parameterizations` must be a named list.")
  }

  if (anyNA(names(parameterizations)) || any(names(parameterizations) == "")) {
    .abort("Every parameterization must have a name.")
  }

  if (anyDuplicated(names(parameterizations))) {
    .abort("Parameterization names must be distinct.")
  }


  # a misspelt name would otherwise be dropped silently
  known <- .predictor_names(original_model(model))
  if (inherits(model, "parameterized_model")) {
    known <- c(known, model$parameters)
  }
  unknown <- setdiff(names(parameterizations), known)
  if (length(unknown)) {
    .abort("Cannot parameterize ", toString(unknown), ": not a predictor of the model",
           if (inherits(model, "parameterized_model")) " or a parameter of the parameterized model", ".")
  }

  new_parameterizations <- parameterizations <- as.functions(parameterizations)

  if (inherits(model, "parameterized_model")) {

    history <- model$parameterization_history

    parameterizations <- .compose_parameterizations(
      model$parameterization,
      new_parameterizations
    )

    model <- model$model
  }else{
    history <- list()
  }

  parameterized_model(
    model = model,
    parameterization = parameterizations,
    parameterization_history = c(history,list(new_parameterizations))
  )
}


#' Construct a parameterized model
#'
#' Constructs a `parameterized_model` from an original model and a collection
#' of parameterization functions.
#'
#' This is the low-level constructor used by [parameterize()].
#' Parameterized_model() expects the parameterizations to already be represented
#' as functions.
#'
#' Only parameterizations corresponding to variables occurring in the model
#' are retained. The variables required as arguments to the retained
#' parameterization functions are stored as the model's parameters.
#'
#' @param model A model or parameterized_model.
#' @param parameterization A named list of functions defining
#' parameterizations of the model's continuous predictor variables.
#' @param parameterization_history A list recording the successive
#' parameterizations applied to the model.
#'
#' @return
#' A parameterized_model object containing the original model, its current
#' parameterizations, the variables required by those parameterizations, and
#' its parameterization history.
#'
#' @seealso
#' [parameterize()]
#'
#' @export
parameterized_model <- function(model, parameterization = NULL, parameterization_history = list()) {
  if(is.null(parameterization)){
    return(model)
  }

  if (!is.list(parameterization) || is.null(names(parameterization))){
    .abort("`parameterization` must be a named list.")
  }

  if (anyDuplicated(names(parameterization))) {
    .abort("Parameterization names must be distinct.")
  }

  if (!all(vapply(parameterization, is.function, logical(1)))) {
    .abort("All parameterizations must be functions.")
  }


  # TODO: add checking later

  # model_data <- insight::get_data(model)
  # if(any(!names(parameterization %in% names(model_data)[sapply(model_data, is.numeric)]))){
  #   .abort("Non-numeric variables cannot be parameterized")
  # }


  parameterization <-
    parameterization[names(parameterization) %in% .predictor_names(model)]

  parameters <- .get_parameters(parameterization)

  structure(
    list(
      model = model,
      parameterization = parameterization,
      parameters = parameters,
      parameterization_history = parameterization_history
    ),
    class = c("parameterized_model", class(model))
  )
}


#' Predict from a parameterized model
#'
#' Applies the parameterizations to newdata and generates predictions from
#' the original model.
#'
#' @param object A parameterized_model.
#' @param newdata A data frame containing the parameter variables (including
#' unparameterized variables).
#' @param ... Arguments passed to the prediction method of the original model.
#'
#' @return The predictions from the original model.
#'
#'
#' @examples
#' data <- data.frame(
#'   y = 1:10,
#'   x = 1:10)
#'
#' model <- lm(y ~ x, data = data)
#' model <- parameterize(model, list(x = function(a, b) a + b))
#'
#' predict(
#'   model,
#'   newdata = data.frame(a = 1:5, b = 1:5))
#'
#' @exportS3Method stats::predict
predict.parameterized_model <- function(object, newdata, ...) {
  if (missing(newdata)) {
    .abort("`newdata` must be supplied for a parameterized model.")
  }

  model <- object$model
  newdata <- parameterize_data(object, newdata)

  newdata <- newdata[names(newdata) %in% .predictor_names(model)]

  predict(object$model, newdata = newdata, ...)
}


#' Get the original model
#'
#' Extracts the original model from a parameterized_model.
#'
#' @param x A parameterized_model.
#'
#' @return The original model.
#'
#' @export
original_model <- function(x) {
  if(is.parameterized_model(x))
    return(x$model)

  x
}


#' Get parameterizations
#'
#' Extracts the current parameterizations, or the compositions of the last
#' level parameterizations in the parameterization history.
#'
#' @param x A parameterized_model.
#' @param level The number of parameterization steps to return.
#'
#' @return A named list of parameterization functions.
#'
#' @examples
#' data <- data.frame(
#' y = 1:10,
#' x = 1:10
#' )
#'
#' model <- lm(y ~ x, data = data)
#' model <- parameterize(
#' model,
#' list(x = function(a, b) a + b)
#' )
#'
#' parameterizations(model)
#'
#' # Add another parameterization
#' model <- parameterize(
#' model,
#' list(x = function(c) c)
#' )
#'
#' # Get the parameterizations from the most recent step
#' parameterizations(model, level = 1)
#'
#' # Get the parameterizations from the last two steps, composed together
#' parameterizations(model, level = 2)
#'
#'
#' @export
parameterizations <- function(x, level = NULL) {
  .check_parameterized_model(x)

  if(is.null(level)){
    return(x$parameterization)
  }else{
    if(! is.numeric(level)) .abort("level must be numeric")

    params <- x$parameterization_history
    params <- x$parameterization_history[(length(params) - level + 1): length(params)]

    while (length(params) > 1) {
      params[[1]] <- .compose_parameterizations(
        params[[1]],
        params[[2]]
      )
      params <- params[-2]
    }
    return(params)
  }
}


#' Get parameter variables
#'
#' Extracts the variables required by the parameterizations of a
#' parameterized_model.
#'
#' @param x A parameterized_model.
#'
#' @return A character vector of parameter variables.
#'
#' @export
parameter_variables <- function(x) {
  .check_parameterized_model(x)

  x$parameters
}


#' Apply parameterizations to new data
#'
#' Applies the parameterizations of a parameterized_model to newdata.
#'
#' @param object A parameterized_model.
#' @param newdata A data frame containing the variables required by the
#' parameterizations.
#' @param level The number of parameterization steps to apply.
#'
#' @return A data frame with the parameterized variables added.
#'
#' @export
parameterize_data <- function(object, newdata, level = NULL) {
  if(!is.parameterized_model(object)) {
    return(newdata)
  }

  if (!is.data.frame(newdata)) {
    newdata <- as.data.frame(newdata)
  }

  .parameterize_data_frame(
    newdata,
    parameterizations(object, level = level)
  )
}

.parameterize_data_frame <- function(data, parameterization) {

  for (nm in names(parameterization)) {
    fun <- parameterization[[nm]]
    parameters <- names(formals(fun))

    for(params in parameters){
      data[params] <- as.numeric(as.character(data[,params]))
    }

    missing_parameters <- setdiff(
      parameters,
      names(data)
    )

    if (length(missing_parameters) > 0) {
      .abort(
        paste0(
          "The following parameters are missing from `data`: ",
          paste(missing_parameters, collapse = ", "),
          "."
        )
      )
    }

    data[[nm]] <- do.call(
      fun,
      as.list(data[parameters])
    )
  }

  data
}


#' Print a parameterized model
#'
#' Prints a summary of a parameterized_model, including the original model,
#' parameter variables, and current parameterizations.
#'
#' @param x A parameterized_model.
#' @param show_history Logical; whether to print the parameterization history.
#' @param ... Currently ignored.
#'
#' @return Invisibly returns x.
#'
#' @exportS3Method base::print
print.parameterized_model <- function(x, show_history = FALSE,  ...) {
  cat("Parameterized model\n\n")

  cat("Model:\n")
  print(x$model)


  cat("\nParameters:\n")
  print(x$parameters)

  cat("\nParameterizations:\n")

  for (nm in names(x$parameterization)) {
    fun <- x$parameterization[[nm]]

    cat(
      "$", nm, " <- ",
      paste(deparse(body(fun)), collapse = " "),
      "; (",
      paste(names(formals(fun)), collapse = ", "),
      ")\n",
      sep = ""
    )
  }

  if(show_history){
    cat("\n\n")
    .print_parameterization_history(x$parameterization_history)

  }



  invisible(x)
}



# ============================================================
# Helpers
# ============================================================


.as_parameterization_function <- function(x) {
  if (is.function(x)) {
    body(x) <- .parenthesize_body(body(x))
    return(x)
  }

  if (is.character(x)) {
    if (length(x) != 1L) {
      .abort(
        "Character parameterizations must have length one."
      )
    }

    x <- parse(text = x)[[1L]]
  }

  if (is.expression(x)) {
    if (length(x) != 1L) {
      .abort(
        "Expression parameterizations must have length one."
      )
    }

    x <- x[[1L]]
  }

  if (!is.call(x) && !is.name(x)) {
    .abort(
      "Parameterizations must be functions, expressions, or character strings."
    )
  }

  parameters <- all.vars(x)

  fun <- .make_function(
    parameters = parameters,
    body = x,
    envir = parent.frame()
  )

  fun
}

.parenthesize_body <- function(body) {
  if (is.call(body) &&
      identical(body[[1L]], as.name("(")) &&
      length(body) == 2L) {
    return(body)
  }

  as.call(list(as.name("("), body))
}


as.functions <- function(parameterizations) {
  result <- lapply(parameterizations, .as_parameterization_function)

  names(result) <- names(parameterizations)

  result
}


.get_parameters <- function(parameterization) {
  unique(
    unlist(
      lapply(parameterization, function(x) names(formals(x))),
      use.names = FALSE
    )
  )
}

.compose_parameterizations <- function(old, new) {

  old_names <- names(old)

  # Names appearing anywhere in old:
  old_used_names <- unique(c(
    old_names,
    unlist(lapply(old, function(f) names(formals(f))),
           use.names = FALSE)
  ))

  result <- lapply(old, function(f) {

    old_parameters <- names(formals(f))

    replacements <- new[
      intersect(old_parameters, names(new))
    ]

    # Nothing to replace
    if (!length(replacements)) {
      return(f)
    }

    body <- .unwrap_body(body(f))

    # Substitute matching formals
    for (nm in names(replacements)) {
      body <- .substitute_expression(
        body,
        nm,
        .unwrap_body(body(replacements[[nm]]))
      )
    }

    # Remaining old parameters + parameters introduced by replacements
    parameters <- c(
      setdiff(old_parameters, names(replacements)),
      .get_parameters(replacements)
    )

    .make_function(
      parameters = unique(parameters),
      body = body,
      envir = environment(f)
    )
  })

  names(result) <- old_names

  # Add new functions whose names appear nowhere in old
  additions <- new[
    setdiff(names(new), old_used_names)
  ]

  c(result, additions)
}


.unwrap_body <- function(x) {
  if (is.call(x) &&
      identical(x[[1L]], as.name("(")) &&
      length(x) == 2L) {
    return(x[[2L]])
  }

  x
}

.substitute_expression <- function(expr, name, replacement) {
  if (is.symbol(expr)) {
    if (identical(as.character(expr), name)) {
      return(replacement)
    }

    return(expr)
  }

  if (!is.call(expr)) {
    return(expr)
  }

  as.call(
    lapply(
      expr,
      function(x) {
        .substitute_expression(
          x,
          name,
          replacement
        )
      }
    )
  )
}

.make_function <- function(parameters, body, envir = parent.frame()) {
  body <- .parenthesize_body(
    .unwrap_body(body)
  )

  formals <- as.pairlist(
    setNames(
      vector("list", length(parameters)),
      parameters
    )
  )

  eval(
    call("function", formals, body),
    envir = envir
  )
}

.parameterize_expression <- function(expr, parameterization) {
  if (is.symbol(expr) && (name <- as.character(expr)) %in% names(parameterization)) {
    f <- parameterization[[name]]
    # arg <- names(formals(f))[1]
    return(do.call("substitute", list(body(f)[[2:length(body(f))]], list(expr))))
  }
  if (is.call(expr))
    expr[-1] <- lapply(expr[-1], .parameterize_expression, parameterization)
  expr
}

.parameterize_formula <- function(formula, parameterization) {
  if (!inherits(formula, "formula")) {
    .abort("`formula` must be a formula.")
  }

  if (!is.list(parameterization) || is.null(names(parameterization))) {
    .abort(
      "`parameterization` must be a named list of functions."
    )
  }

  if (!all(vapply(parameterization, is.function, logical(1)))) {
    .abort(
      "All elements of `parameterization` must be functions."
    )
  }

  result <- formula

  result[[length(result)]] <- .parameterize_formula_node(
    result[[length(result)]],
    lapply(parameterization, body)
  )

  result
}


.parameterize_variables <- function(variables, parameterization) {
  vapply(
    variables,
    function(variable) {
      if (variable %in% names(parameterization)) {
        paste(
          deparse(body(parameterization[[variable]])),
          collapse = ""
        )
      } else {
        variable
      }
    },
    character(1)
  )
}


.parameterize_formula_node <- function(expr, parameterization) {
  if (is.symbol(expr)) {
    name <- as.character(expr)

    if (name %in% names(parameterization)) {
      return(
        parameterization[[name]]
      )
    }

    return(expr)
  }

  if (!is.call(expr)) {
    return(expr)
  }

  # Do not replace the function name itself.
  #
  # For example:
  #
  #   log(x)
  #
  # should recurse into x, but should not try to
  # parameterize "log".
  result <- as.list(expr)

  if (length(result) > 1L) {
    result[-1L] <- lapply(
      result[-1L],
      .parameterize_formula_node,
      parameterization = parameterization
    )
  }

  as.call(result)
}

.variable_parameters <- function(variables, parameterization) {
  unique(unlist(
    lapply(
      variables,
      function(variable) {
        if (variable %in% names(parameterization)) {
          names(formals(parameterization[[variable]]))
        } else {
          variable
        }
      }
    ),
    use.names = FALSE
  ))
}

parameterize_variables <- function(variables, parameterization) {
  out <- unique(
    unlist(
      lapply(
        variables,
        function(variable) {
          if (!variable %in% names(parameterization)) {
            return(variable)
          }

          .parameterize_variables(
            names(formals(parameterization[[variable]])),
            parameterization
          )
        }
      ),
      use.names = FALSE
    )
  )

  out
}

#' Test for a parameterized model
#'
#' Tests whether an object inherits from the parameterized_model class.
#'
#' @param x An R object.
#'
#' @return A logical value indicating whether x is a parameterized_model.
#'
#' @export
is.parameterized_model <- function(x){
  inherits(x, "parameterized_model")
}

.check_parameterized_model <- function(x) {
  if (!is.parameterized_model(x)) {
    .abort("`x` must be a `parameterized_model`.")
  }
}

#' @export
`$.parameterized_model` <- function(x, name){
  if(name %in% names(x)){
    x[[name]]
  }else{
    x$model[[name]]
  }

}




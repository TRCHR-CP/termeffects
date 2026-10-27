#' Test whether a reference function fits inside the simultaneous band
#'
#' Compares the predictions with a reference function `f` at every grid point
#' at once. The null hypothesis is that the true prediction equals `f` at all
#' grid points (with the other predictors held as in `x`). It is rejected at
#' level `1 - ci_level` exactly when `f` leaves the simultaneous band
#' somewhere.
#'
#' The test statistic is the largest standardized distance between estimate
#' and reference, `max |estimate - f| / SE`, computed on the scale the band
#' is built on (e.g. log-odds for logistic models). Its p-value comes from the
#' same simulation (or bootstrap) as the band's critical value.
#'
#' This tests equality to one fully specified function, e.g. "no effect" or a
#' published curve. It is not a general test of linearity: a linear function
#' with estimated coefficients would need to be fitted, and its uncertainty
#' accounted for.
#'
#' The idea and the original implementation (comparison of `f` with the band)
#' are by Lukas Bolsinger for p-ggeffects; the p-value is new.
#'
#' @param x A [term_predictions()] result.
#' @param f A function whose argument names are variables of the original
#'   model (e.g. `function(age, bmi) 50 + 0 * age`), returning the reference
#'   values on the scale of `estimate`. Held and parameterized predictors can
#'   be used as well as focal ones.
#'
#' @return `x` with added columns `reference` (the values of `f`),
#'   `difference` (`estimate - reference`) and `inside` (whether the reference
#'   lies inside the band at that point), and class `term_test`. Attributes
#'   `statistic`, `p_value` (one per facet or curve if the band's `ci_scope`
#'   is not global) and `all_inside`.
#'
#' @examples
#' fit <- lm(mpg ~ splines::ns(hp, 3) + wt, data = mtcars)
#' set.seed(1)
#' pred <- term_predictions(fit, "hp")
#' # is mpg flat in hp (at the median weight)?
#' flat <- test_simultaneous_predictions(pred, function(hp) mean(pred$estimate) + 0 * hp)
#' flat
#' plot(flat)
#' @export
test_simultaneous_predictions <- function(x, f) {
  if (!inherits(x, "term_predictions")) {
    stop("`x` must be a term_predictions() result.", call. = FALSE)
  }
  args <- names(formals(f))
  missing_args <- setdiff(args, names(x))
  if (length(missing_args)) {
    stop("`f` uses ", toString(missing_args), ", which ", if (length(missing_args) > 1) "are" else "is",
         " not a variable of the predictions. Available: ",
         toString(setdiff(names(x), c("estimate", "std.error", "conf.low", "conf.high",
                                      "band.low", "band.high"))), ".", call. = FALSE)
  }

  inputs <- lapply(as.data.frame(x)[args], function(v) if (is.factor(v)) as.character(v) else v)
  reference <- do.call(f, inputs)
  if (length(reference) == 1) {
    reference <- rep(reference, nrow(x))
  }
  if (length(reference) != nrow(x)) {
    # fall back to calling f one row at a time
    reference <- vapply(seq_len(nrow(x)), function(i) {
      do.call(f, lapply(inputs, `[`, i))
    }, numeric(1))
  }

  to_scale <- attr(x, "to_inference_scale")
  t_values <- (to_scale(x$estimate) - suppressWarnings(to_scale(reference))) / x$std.error
  t_values[x$std.error == 0] <- 0

  group <- attr(x, "scope_group")
  max_t <- attr(x, "max_t")
  statistic <- c(tapply(abs(t_values), group, max))
  p_value <- vapply(names(max_t), function(g) {
    draws <- max_t[[g]]
    (1 + sum(draws >= statistic[[g]])) / (1 + length(draws))
  }, numeric(1))

  x$reference <- reference
  x$difference <- x$estimate - reference
  x$inside <- reference >= x$band.low & reference <= x$band.high

  attr(x, "statistic") <- if (length(statistic) == 1) unname(statistic[[1]]) else statistic[names(max_t)]
  attr(x, "p_value") <- if (length(p_value) == 1) unname(p_value) else p_value
  attr(x, "all_inside") <- all(x$inside)
  class(x) <- c("term_test", class(x))
  x
}


#' @export
print.term_test <- function(x, ...) {
  level <- paste0(100 * attr(x, "ci_level"), "%")
  cat("Simultaneous test: does the reference function fit inside the", level, "band?\n")
  held <- .held_text(x)
  if (length(held)) {
    cat(strwrap(paste("Held constant:", toString(held)), exdent = 2), sep = "\n")
  }
  stat <- attr(x, "statistic")
  p <- attr(x, "p_value")
  lines <- sprintf("max |t| = %.2f, p = %s", stat, format.pval(p, digits = 2, eps = 1e-4))
  if (length(stat) > 1) {
    lines <- paste0(names(stat), ": ", lines)
  }
  cat(lines, sep = "\n")
  cat(sprintf("The reference lies inside the band at %d of %d grid points.\n\n",
              sum(x$inside), nrow(x)))
  cols <- c(attr(x, "focal"), "estimate", "reference", "band.low", "band.high", "inside")
  print(utils::head(as.data.frame(x)[intersect(cols, names(x))], 10), digits = 4, row.names = FALSE)
  invisible(x)
}

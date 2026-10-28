# Default ranges for the new variables of a parameterized model.
#
# .get_excribed_box() was written by Lukas Bolsinger for p-ggeffects
# (R/utils_get_model_range.R) and is copied with only insight messages
# replaced. .parameter_ranges() is new: it calls it once per
# parameterization, holding arguments that are data variables at their
# grid values, and keeps the box that fits all parameterizations.

.parameter_ranges <- function(pmodel, free, fixed, data) {
  ranges <- list()
  for (target in names(parameterizations(pmodel))) {
    f <- parameterizations(pmodel)[[target]]
    args <- names(formals(f))
    vary <- intersect(args, free)
    if (!length(vary) || !is.numeric(data[[target]])) {
      next
    }
    held <- fixed[intersect(setdiff(args, vary), names(fixed))]
    # f as a function of the free variables only
    g <- function() NULL
    formals(g) <- stats::setNames(rep(list(quote(expr = )), length(vary)), vary)
    body(g) <- quote(do.call(f, c(mget(vary), held)))
    box <- .get_excribed_box(g, range(data[[target]], na.rm = TRUE), verbose = FALSE)
    for (p in names(box)) {
      # a parameter used by several parameterizations gets the box that keeps
      # all of their images inside the observed data range
      ranges[[p]] <- if (is.null(ranges[[p]])) box[[p]] else
        c(max(ranges[[p]][1], box[[p]][1]), min(ranges[[p]][2], box[[p]][2]))
    }
  }
  ranges
}

.get_excribed_box <- function(f, target,
                               r0 = 0.1,
                               growth = 1.2,
                               shrink = 0.9,
                               eps = 1e-6,
                               max_grow = 100,
                               max_shrink = 100,
                               n_random = 500,
                               n_starts = 10,
                               seed = NULL,
                               verbose = TRUE) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  if (length(target) != 2) {
    .abort("target must be have length(target) == 2")
  }

  target_lower <- min(target)
  target_upper <- max(target)

  # Get argument names and dimension
  arg_names <- names(formals(f))
  arg_names <- arg_names[arg_names != "..."]
  d <- length(arg_names)

  if (d == 0) {
    .abort("Could not determine the dimension of f")
  }

  # Turn f(x1, x2, ...) into a function of a vector
  f_vec <- function(x) {
    do.call(f, as.list(x))
  }

  # Initial radii
  if (length(r0) == 1) {
    r <- rep(r0, d)
  } else if (length(r0) == d) {
    r <- r0
  } else {
    .abort("r0 must have length 1 or length equal to dimension of f")
  }

  if (any(r <= 0)) {
    .abort("All radii must be positive")
  }

  # Generate random points in the box
  random_points <- function(r, n) {

    X <- matrix(
      runif(n * d),
      nrow = n,
      ncol = d
    )

    for (j in seq_len(d)) {
      X[, j] <- (2 * X[, j] - 1) * r[j]
    }

    return(X)
  }

  # Estimate minimum and maximum of f over the box
  image_bounds <- function(r) {

    X <- random_points(r, n_random)
    values <- apply(X, 1, f_vec)

    good <- is.finite(values)

    if (!any(good)) {
      .abort("f returned no finite values")
    }

    X <- X[good, , drop = FALSE]
    values <- values[good]

    n <- min(n_starts, length(values))

    min_indices <- order(values)[seq_len(n)]
    max_indices <- order(values, decreasing = TRUE)[seq_len(n)]

    min_values <- numeric(n)
    max_values <- numeric(n)

    # Find minimum
    for (i in seq_len(n)) {

      opt <- optim(
        par = X[min_indices[i], ],
        fn = f_vec,
        method = "L-BFGS-B",
        lower = -r,
        upper = r
      )

      min_values[i] <- opt$value
    }

    # Find maximum
    for (i in seq_len(n)) {

      opt <- optim(
        par = X[max_indices[i], ],
        fn = function(x) -f_vec(x),
        method = "L-BFGS-B",
        lower = -r,
        upper = r
      )

      max_values[i] <- -opt$value
    }

    return(c(
      lower = min(min_values),
      upper = max(max_values)
    ))
  }

  # Check whether the image is contained in the target range
  covers <- function(image) {

    image["lower"] >= target_lower - eps &&
      image["upper"] <= target_upper + eps
  }

  # Initial box
  image <- image_bounds(r)

  if (verbose) {
    cat(
      "initial: r =",
      paste(round(r, 5), collapse = ", "),
      " image = [",
      round(image["lower"], 5),
      ", ",
      round(image["upper"], 5),
      "]\n"
    )
  }

  # Find a feasible starting box if the initial box is too large
  shrink_iter <- 0

  while (!covers(image) && shrink_iter < max_shrink) {
    shrink_iter <- shrink_iter + 1
    r <- r * shrink
    image <- image_bounds(r)
  }

  if (!covers(image)) {
    warning("Could not find a rectangle whose image is contained in target", " Give the range with `values = `.", call. = FALSE)
    result <- lapply(r, function(radius) c(-radius, radius))
    names(result) <- arg_names
    return(result)
  }

  # Grow dimensions individually while keeping the image inside the target
  grow_iter <- 0

  for (iter in seq_len(max_grow)) {

    changed <- FALSE

    for (j in seq_len(d)) {

      r_new <- r
      r_new[j] <- r_new[j] * growth

      image_new <- image_bounds(r_new)

      if (covers(image_new)) {

        r <- r_new
        image <- image_new
        changed <- TRUE

        if (verbose) {
          cat("grow", iter, "dimension", arg_names[j],
              ": r =", paste(round(r, 5), collapse = ", "),
              " image = [", round(image["lower"], 5), ", ",
              round(image["upper"], 5), "]\n")
        }
      }
    }

    if (!changed) {
      break
    }
  }

  # Return named list of intervals
  result <- lapply(r, function(radius) {
    c(-radius, radius)
  })

  names(result) <- arg_names

  return(result)
}

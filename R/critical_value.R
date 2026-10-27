#' Critical value for a simultaneous confidence band
#'
#' Finds the multiplier `c` such that the band `estimate +/- c * SE` covers
#' the true values at *all* grid points at once with probability `ci_level`.
#' This is the `ci_level` quantile of `max_i |Z_i| / SE_i`, where `Z` is
#' normally distributed with the joint covariance of the estimates (or
#' multivariate t when `df` is finite). The quantile is found by simulation.
#'
#' The result depends on the random number generator. Call [set.seed()]
#' first for reproducible values; with the default `nsim` the simulation
#' error is about +/- 0.02.
#'
#' @param X Design matrix (grid points in rows).
#' @param vcov Covariance matrix of the coefficients.
#' @param ci_level Simultaneous coverage probability.
#' @param df Degrees of freedom; `Inf` for a normal reference distribution.
#' @param nsim Number of simulation draws.
#'
#' @return A single number.
#' @noRd
.critical_value <- function(X, vcov, ci_level = 0.95, df = Inf, nsim = 10000) {
  # V = L L', so the estimates X beta have covariance (X L)(X L)'
  e <- eigen((vcov + t(vcov)) / 2, symmetric = TRUE)
  L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), nrow = length(e$values))
  A <- X %*% L
  se <- sqrt(rowSums(A^2))

  # points with no uncertainty (e.g. the Cox reference point) cannot exceed
  # the band, so they are left out of the maximum
  A <- A[se > sqrt(.Machine$double.eps) * max(se), , drop = FALSE]
  A <- A / sqrt(rowSums(A^2))
  if (nrow(A) == 0) {
    return(NA_real_)
  }

  # simulate in chunks to keep memory use small on large grids
  chunk <- 2000
  max_t <- unlist(lapply(seq(1, nsim, by = chunk), function(start) {
    m <- min(chunk, nsim - start + 1)
    Z <- matrix(stats::rnorm(m * ncol(A)), nrow = m) %*% t(A)
    if (is.finite(df)) {
      Z <- Z / sqrt(stats::rchisq(m, df) / df)
    }
    apply(abs(Z), 1, max)
  }))
  as.numeric(stats::quantile(max_t, ci_level, names = FALSE))
}

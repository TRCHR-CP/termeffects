# Same comparison as 08_orm_interval_methods.R, for a model as complex as the
# heart-failure example in the how-to vignette: n = 300, a skewed predictor
# like creatinine, one like sodium, a 4 x 4 knot spline interaction and six
# adjustment covariates (about 24 coefficients), with a 10-level outcome.
# In the real data the bootstrap band was much wider than the default band;
# this checks which of them has the stated coverage in such a model.
#
# The true linear predictor lies in the model space (linear terms and their
# product are contained in rcs(x, 4) * rcs(z, 4)), so there is no bias.
pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(rms))
RNGkind("L'Ecuyer-CMRG")
set.seed(2027)
nrep <- 300
nboot <- 1000
cores <- max(1, parallel::detectCores() - 1)

n <- 300
values <- c(15, 20, 25, 30, 35, 40, 50, 60, 70, 80)
cuts <- qlogis(c(0.05, 0.12, 0.22, 0.35, 0.5, 0.63, 0.78, 0.9, 0.96))
lp_true <- function(x, z, a1 = 0, c1 = 0) {
  -0.8 * (x - 1.2) + 0.12 * (z - 137) + 0.1 * (x - 1.2) * (z - 137) + 0.4 * a1 + 0.3 * c1
}
true_mean <- function(lp) {
  p_upper <- sapply(cuts, function(cut) plogis(lp - cut))
  drop((cbind(1, p_upper) - cbind(p_upper, 0)) %*% values)
}

simulate <- function() {
  d <- data.frame(x = exp(rnorm(n, 0.2, 0.45)), z = rnorm(n, 137, 4.5),
                  a1 = rbinom(n, 1, 0.4), a2 = rbinom(n, 1, 0.3), a3 = rbinom(n, 1, 0.6),
                  c1 = rnorm(n), c2 = rnorm(n), c3 = rnorm(n))
  latent <- lp_true(d$x, d$z, d$a1, d$c1) + rlogis(n)
  d$y <- values[findInterval(latent, cuts) + 1]
  d
}

one_run <- function() {
  d <- simulate()
  fit <- tryCatch(orm(y ~ rcs(x, 4) * rcs(z, 4) + a1 + a2 + a3 + rcs(c1, 3) + rcs(c2, 3) + rcs(c3, 3),
                      data = d, x = TRUE, y = TRUE), error = function(e) NULL)
  if (is.null(fit) || isTRUE(fit$fail)) {
    return(NULL)
  }
  grid <- list(x = seq(0.6, 3.5, length.out = 21), z = c(130, 137, 142))
  held <- list(a1 = 0, a2 = 0, a3 = 0, c1 = 0, c2 = 0, c3 = 0)
  get <- function(method) {
    tryCatch(term_predictions(fit, c("x", "z"), values = grid, at = held, ci_method = method,
                              nsim = 2000, nboot = nboot), error = function(e) NULL)
  }
  res <- list(lp_wald = get("analytic"), delta = get("delta"), bootstrap = get("bootstrap"))
  truth <- true_mean(lp_true(res$lp_wald$x, res$lp_wald$z))
  mid <- which(res$lp_wald$x == grid$x[11] & res$lp_wald$z == 137)
  unlist(lapply(names(res), function(m) {
    p <- res[[m]]
    if (is.null(p)) {
      return(stats::setNames(rep(NA, 3), paste0(m, c("", "_point", "_width"))))
    }
    stats::setNames(c(all(truth >= p$band.low - 1e-8 & truth <= p$band.high + 1e-8),
                      truth[mid] >= p$conf.low[mid] & truth[mid] <= p$conf.high[mid],
                      mean(p$band.high - p$band.low)),
                    paste0(m, c("", "_point", "_width")))
  }))
}

runs <- parallel::mclapply(seq_len(nrep), function(i) one_run(), mc.cores = cores)
runs <- do.call(rbind, Filter(Negate(is.null), runs))
m <- colMeans(runs, na.rm = TRUE)
res <- data.frame(scenario = "complex_lvef_like", method = c("lp_wald", "delta", "bootstrap"),
                  simultaneous_coverage = m[c("lp_wald", "delta", "bootstrap")],
                  pointwise_coverage_at_centre = m[c("lp_wald_point", "delta_point", "bootstrap_point")],
                  mean_band_width = m[c("lp_wald_width", "delta_width", "bootstrap_width")],
                  runs = colSums(!is.na(runs[, c("lp_wald", "delta", "bootstrap")])), row.names = NULL)
print(res, digits = 3)
write.csv(res, "validation/output/10_orm_intervals_complex_model.csv", row.names = FALSE)

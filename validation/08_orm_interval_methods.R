# Which intervals for the mean of an rms::orm model are right?
#
# Data are generated from a true proportional odds model, so the true mean
# outcome at every grid point is known exactly. For each simulated data set
# we compute the 95% simultaneous band over three curves (21 x-values each)
# in three ways and record whether it contains all 63 true means:
#
#   lp_wald   termeffects' default (ci_method = "analytic"): band for rms's
#             linear predictor (reference intercept + slopes), transformed
#             with Mean(); what rms::Predict(fun = Mean(fit)) does. Ignores
#             the uncertainty of the other intercepts.
#   delta     ci_method = "delta": delta method on the mean scale with the
#             covariance of ALL intercepts and slopes (as rms's
#             Mean(fit)(X = , conf.int = )).
#   bootstrap ci_method = "bootstrap": rms::bootcov() refits, sup-t band on
#             the mean scale.
#
# Result of the first run (nboot = 200): lp_wald 0.973 / 0.947, delta
# 0.940 / 0.917, bootstrap 0.930 / 0.910 (LVEF-like / continuous), which is
# why the linear-predictor method is the default. This version uses
# nboot = 1000, the package default.
#
# Also recorded: pointwise coverage at one point and the average band width.
# Target for the simultaneous coverage: 0.95 (about 0.925-0.975 with 300 runs).
pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(rms))
RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
nrep <- 300
nboot <- 1000
cores <- max(1, parallel::detectCores() - 1)

lp_true <- function(x, z) 2 * x - 3 * x^2 + 2 * x * z
grid_x <- seq(0, 1, length.out = 21)
grid_z <- c(0.2, 0.5, 0.8)

scenarios <- list(
  # LVEF-like: 10 ordered values between 15 and 80, small sample
  lvef_like = list(n = 200, values = c(15, 20, 25, 30, 35, 40, 50, 60, 70, 80),
                   cuts = qlogis(c(0.05, 0.12, 0.22, 0.35, 0.5, 0.63, 0.78, 0.9, 0.96))),
  # continuous outcome: every observation its own value
  continuous = list(n = 300, values = NULL, cuts = NULL)
)

true_mean <- function(sc, lp) {
  if (is.null(sc$values)) {
    return(lp)  # Y = lp + logistic error, mean = lp
  }
  # P(Y = value_k) from the latent logistic variable lp + e cut at `cuts`
  p_upper <- sapply(sc$cuts, function(cut) plogis(lp - cut))  # P(Y > value_k)
  p <- cbind(1, p_upper) - cbind(p_upper, 0)
  drop(p %*% sc$values)
}

simulate <- function(sc) {
  d <- data.frame(x = runif(sc$n), z = runif(sc$n))
  latent <- lp_true(d$x, d$z) + rlogis(sc$n)
  d$y <- if (is.null(sc$values)) latent else sc$values[findInterval(latent, sc$cuts) + 1]
  d
}

one_run <- function(sc) {
  d <- simulate(sc)
  fit <- orm(y ~ pol(x, 2) * z, data = d, x = TRUE, y = TRUE)
  values <- list(x = grid_x, z = grid_z)
  delta <- term_predictions(fit, c("x", "z"), values = values, nsim = 2000, ci_method = "delta")
  lp <- term_predictions(fit, c("x", "z"), values = values, nsim = 2000)
  boot <- tryCatch(term_predictions(fit, c("x", "z"), values = values, ci_method = "bootstrap",
                                    nboot = nboot), error = function(e) NULL)
  truth <- true_mean(sc, lp_true(delta$x, delta$z))
  covers <- function(lo, hi) all(truth >= lo - 1e-8 & truth <= hi + 1e-8)
  mid <- which(delta$x == 0.5 & delta$z == 0.5)
  c(lp_wald = covers(lp$band.low, lp$band.high),
    delta = covers(delta$band.low, delta$band.high),
    bootstrap = if (is.null(boot)) NA else covers(boot$band.low, boot$band.high),
    lp_wald_point = truth[mid] >= lp$conf.low[mid] & truth[mid] <= lp$conf.high[mid],
    delta_point = truth[mid] >= delta$conf.low[mid] & truth[mid] <= delta$conf.high[mid],
    bootstrap_point = if (is.null(boot)) NA else truth[mid] >= boot$conf.low[mid] & truth[mid] <= boot$conf.high[mid],
    lp_wald_width = mean(lp$band.high - lp$band.low),
    delta_width = mean(delta$band.high - delta$band.low),
    bootstrap_width = if (is.null(boot)) NA else mean(boot$band.high - boot$band.low))
}

res <- do.call(rbind, lapply(names(scenarios), function(name) {
  runs <- parallel::mclapply(seq_len(nrep), function(i) one_run(scenarios[[name]]), mc.cores = cores)
  runs <- do.call(rbind, runs)
  m <- colMeans(runs, na.rm = TRUE)
  data.frame(scenario = name, method = c("lp_wald", "delta", "bootstrap"),
             simultaneous_coverage = m[c("lp_wald", "delta", "bootstrap")],
             pointwise_coverage_at_centre = m[c("lp_wald_point", "delta_point", "bootstrap_point")],
             mean_band_width = m[c("lp_wald_width", "delta_width", "bootstrap_width")],
             runs = nrow(runs), row.names = NULL)
}))
print(res, digits = 3)
write.csv(res, "validation/output/08_orm_interval_methods.csv", row.names = FALSE)

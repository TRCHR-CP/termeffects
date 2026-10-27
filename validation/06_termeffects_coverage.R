# Coverage of termeffects' simultaneous bands: how often does the band contain
# the true curve at ALL grid points? Target: about 0.95 (0.92-0.98 is fine
# with 300 runs). Pointwise intervals should cover the whole curve far less often.
# The true models lie inside the fitted model space, so there is no bias.
pkgload::load_all(".", quiet = TRUE)
set.seed(1)
nsim <- 300

coverage <- function(simulate, fit_model, truth, scale) {
  hits <- replicate(nsim, {
    d <- simulate()
    fit <- fit_model(d)
    p <- term_predictions(fit, c("x", "z"), values = list(z = c(0.2, 0.5, 0.8)),
                          n = 21, nsim = 2000, scale = scale)
    true <- truth(p$x, p$z)
    c(band = all(true >= p$band.low & true <= p$band.high),
      pointwise = all(true >= p$conf.low & true <= p$conf.high))
  })
  rowMeans(hits)
}

lm_cov <- coverage(
  simulate = function() {
    d <- data.frame(x = runif(300), z = runif(300))
    d$y <- 1 + 2 * d$x - 3 * d$x^2 + 2 * d$x * d$z + d$z^2 + rnorm(300)
    d
  },
  fit_model = function(d) lm(y ~ poly(x, 2, raw = TRUE) * poly(z, 2, raw = TRUE), data = d),
  truth = function(x, z) 1 + 2 * x - 3 * x^2 + 2 * x * z + z^2,
  scale = "response"
)

glm_cov <- coverage(
  simulate = function() {
    d <- data.frame(x = runif(600), z = runif(600))
    d$y <- rbinom(600, 1, plogis(-1 + d$x + 2 * d$x * d$z))
    d
  },
  fit_model = function(d) glm(y ~ x * z, family = binomial, data = d),
  truth = function(x, z) plogis(-1 + x + 2 * x * z),
  scale = "response"
)

# orm: continuous outcome with logistic errors, so the proportional odds model
# is exactly right and the true mean is the linear predictor
orm_cov <- coverage(
  simulate = function() {
    d <- data.frame(x = runif(300), z = runif(300))
    d$y <- 2 * d$x + 3 * d$x * d$z + rlogis(300, scale = 0.5)
    d
  },
  fit_model = function(d) rms::orm(y ~ x * z, data = d),
  truth = function(x, z) 2 * x + 3 * x * z,
  scale = "response"
)

res <- data.frame(model = c("lm, polynomial interaction", "logistic glm", "rms::orm (mean)"),
                  simultaneous_band = c(lm_cov["band"], glm_cov["band"], orm_cov["band"]),
                  pointwise = c(lm_cov["pointwise"], glm_cov["pointwise"], orm_cov["pointwise"]))
print(res)
write.csv(res, "validation/output/06_termeffects_coverage.csv", row.names = FALSE)

test_that("one grid point gives the ordinary pointwise critical value", {
  X <- matrix(1, 1, 1)
  set.seed(1)
  expect_equal(.critical_value(X, matrix(1), nsim = 50000), qnorm(0.975), tolerance = 0.02)
  set.seed(1)
  expect_equal(.critical_value(X, matrix(1), df = 10, nsim = 50000), qt(0.975, 10), tolerance = 0.03)
})

test_that("independent points match the exact Sidak value", {
  X <- diag(5)
  set.seed(2)
  sidak <- qnorm(1 - (1 - 0.95^(1 / 5)) / 2)
  expect_equal(.critical_value(X, diag(5), nsim = 50000), sidak, tolerance = 0.02)
})

test_that("agrees with mvtnorm::qmvt for correlated estimates", {
  skip_if_not_installed("mvtnorm")
  fit <- lm(mpg ~ splines::ns(hp, 3) * wt, data = mtcars)
  grid <- term_grid(fit, c("hp", "wt"), values = list(wt = c(2.5, 3.5)), n = 6)
  design <- .design(fit, grid)
  Sigma <- design$X %*% design$vcov %*% t(design$X)
  set.seed(3)
  reference <- mvtnorm::qmvt(0.95, tail = "both.tails", corr = cov2cor(Sigma),
                             df = df.residual(fit))$quantile
  set.seed(3)
  ours <- .critical_value(design$X, design$vcov, df = df.residual(fit), nsim = 50000)
  expect_equal(ours, reference, tolerance = 0.03)
})

test_that("a zero-variance grid point is ignored", {
  X <- rbind(c(0, 0), diag(2))
  set.seed(4)
  expect_true(is.finite(.critical_value(X, diag(2), nsim = 5000)))
})

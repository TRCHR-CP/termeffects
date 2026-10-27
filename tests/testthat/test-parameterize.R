# The first three tests and the range and layer tests are Lukas Bolsinger's
# tests from p-ggeffects (tests/testthat/test-parameterized_model.R); the
# range and layer tests are adapted to term_grid() and term_layer().

test_that("parameterized models predict from parameter variables", {
  dat <- data.frame(
    y = 1:10,
    x = 1:10
  )
  fit <- lm(y ~ x, data = dat)
  pfit <- parameterize(fit, list(x = function(a, b) a + b))
  newdata <- data.frame(a = 1:5, b = 1:5)

  expect_s3_class(pfit, "parameterized_model")
  expect_equal(
    predict(pfit, newdata),
    predict(fit, data.frame(x = newdata$a + newdata$b))
  )
})


test_that("parameterization helpers expose model metadata", {
  fit <- lm(y ~ x, data = data.frame(y = 1:10, x = 1:10))
  pfit <- parameterize(fit, list(x = function(a, b) a + b))

  expect_equal(parameter_variables(pfit), c("a", "b"))
  expect_equal(names(parameterizations(pfit)), "x")
  expect_identical(original_model(pfit), fit)
  expect_equal(
    parameterize_data(pfit, data.frame(a = 1:3, b = 2:4))$x,
    c(3, 5, 7)
  )
})


test_that("parameterizations compose and retain history", {
  fit <- lm(y ~ x, data = data.frame(y = 1:10, x = 1:10))
  pfit <- parameterize(fit, list(x = function(a, b) a + b))
  pfit2 <- parameterize(pfit, list(a = function(t) sin(t)))

  expect_equal(parameter_variables(pfit2), c("b", "t"))
  expect_length(pfit2$parameterization_history, 2)
  expect_equal(
    parameterize_data(pfit2, data.frame(t = 1, b = 2))$x,
    sin(1) + 2
  )
  expect_output(print(pfit2, show_history = TRUE), "Stage 2")
})


test_that("strings define parameterizations and misspelt names are errors", {
  fit <- lm(y ~ x, data = data.frame(y = 1:10, x = 1:10))
  pfit <- parameterize(fit, list(x = "2 * s"))
  expect_equal(parameterize_data(pfit, data.frame(s = 3))$x, 6)
  expect_error(parameterize(fit, list(xx = "2 * s")), "not a predictor")
})


test_that("parameterized ranges stay inside the observed bounding box", {
  set.seed(123)
  dat <- data.frame(
    y = rnorm(100),
    x = runif(100, -2, 2),
    z = runif(100, -2, 2)
  )
  fit <- lm(y ~ x + z, data = dat)
  pfit <- parameterize(
    fit,
    list(
      x = function(r, t) r + t^2,
      z = function(r, t) r^2 + t
    )
  )

  grid <- term_grid(pfit, c("r", "t"), n = 10)
  transformed_x <- with(grid, r + t^2)
  transformed_z <- with(grid, r^2 + t)

  expect_gte(min(transformed_x), min(dat$x))
  expect_lte(max(transformed_x), max(dat$x))
  expect_gte(min(transformed_z), min(dat$z))
  expect_lte(max(transformed_z), max(dat$z))
})


test_that("predictions of a parameterized model equal those on the transformed data", {
  set.seed(1)
  d <- data.frame(x = runif(80), z = runif(80))
  d$y <- d$x + d$z^2 + rnorm(80)
  fit <- lm(y ~ splines::ns(x, 3) * z, data = d)
  pfit <- parameterize(fit, list(x = function(t, c) c * t^2))
  pred <- term_predictions(pfit, c("t", "c"), values = list(t = seq(0.1, 0.9, 0.2), c = c(0.5, 1)))
  direct <- predict(fit, data.frame(x = pred$c * pred$t^2, z = median(d$z)), se.fit = TRUE)
  expect_equal(pred$estimate, unname(direct$fit))
  expect_equal(pred$std.error, unname(direct$se.fit))
  expect_equal(pred$x, pred$c * pred$t^2)
})


test_that("variables of a parameterization need values unless they are focal", {
  fit <- lm(mpg ~ hp + wt, data = mtcars)
  pfit <- parameterize(fit, list(hp = function(t, c) c * t))
  expect_error(term_grid(pfit, "t", values = list(t = 1:3)), "no typical value")
  expect_equal(nrow(term_grid(pfit, "t", values = list(t = 1:3), at = list(c = 50))), 3)
  expect_error(term_grid(pfit, "hp"), "computed by the parameterization")
})


test_that("layers follow parameterized curves through the heatmap", {
  set.seed(123)
  dat <- data.frame(
    y = rnorm(100),
    x = runif(100, -2, 2),
    z = runif(100, -2, 2)
  )
  fit <- lm(y ~ x + z, data = dat)
  out <- term_predictions(fit, c("x", "z"), n = 10)
  layer <- term_layer(out, list(x = "r + t^2", z = "r^2 + t"), terms = c("t", "r"),
                      values = list(t = seq(-1, 1, by = 0.25), r = c(-1, 0, 1)))

  expect_s3_class(layer, "term_layer")
  expect_equal(sort(unique(layer$r)), c(-1, 0, 1))
  expect_equal(layer$x, layer$r + layer$t^2)
  expect_true(all(is.finite(layer$estimate)))

  p <- term_heatmap(out, layers = layer)
  expect_true("GeomPath" %in% vapply(p$layers, function(l) class(l$geom)[1], character(1)))
  expect_silent(ggplot2::ggplot_build(p))
  expect_s3_class(term_curves(layer), "ggplot")

  expect_error(term_layer(out, list(y = "t")), "heatmap variables")
  expect_error(term_layer(out, list(x = function(t) t, z = function(t, x) x + t)), "as an input")
})

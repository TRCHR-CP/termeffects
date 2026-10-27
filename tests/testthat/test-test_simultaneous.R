fit <- lm(mpg ~ splines::ns(hp, 3) + wt, data = mtcars)

test_that("the test agrees with the band and gives a p-value", {
  set.seed(1)
  pred <- term_predictions(fit, "hp", n = 25)
  # the estimate itself: statistic 0, p = 1
  self <- test_simultaneous_predictions(pred, function(hp) pred$estimate)
  expect_equal(attr(self, "statistic"), 0)
  expect_equal(attr(self, "p_value"), 1)
  expect_true(attr(self, "all_inside"))

  flat <- test_simultaneous_predictions(pred, function(hp) mean(pred$estimate) + 0 * hp)
  expect_s3_class(flat, "term_test")
  # rejected exactly when the reference leaves the band somewhere
  expect_equal(attr(flat, "statistic") > attr(pred, "critical_value"), !attr(flat, "all_inside"))
  expect_equal(attr(flat, "p_value") < 0.05, !attr(flat, "all_inside"))
  expect_output(print(flat), "max \\|t\\|")
  expect_s3_class(plot(flat), "ggplot")
  expect_s3_class(plot(flat, "difference"), "ggplot")
})

test_that("the test works on the link scale of a GLM and per curve", {
  g <- glm(am ~ wt + hp, family = binomial, data = mtcars)
  set.seed(2)
  pred <- term_predictions(g, c("wt", "hp"), values = list(hp = c(100, 200)), n = 10,
                           ci_scope = "curve")
  expect_length(attr(pred, "critical_value"), 2)
  half <- test_simultaneous_predictions(pred, function(wt) 0.5 + 0 * wt)
  expect_length(attr(half, "p_value"), 2)
  t_manual <- abs(qlogis(pred$estimate) - qlogis(0.5)) / pred$std.error
  expect_equal(unname(attr(half, "statistic")), unname(c(tapply(t_manual, pred$hp, max))))
})

test_that("f must use variables of the predictions", {
  pred <- term_predictions(fit, "hp", n = 5)
  expect_error(test_simultaneous_predictions(pred, function(cyl) cyl), "cyl")
  expect_s3_class(test_simultaneous_predictions(pred, function(hp, wt) 20 + 0 * hp + 0 * wt), "term_test")
})

test_that("surface tests plot as a difference heatmap", {
  f2 <- lm(mpg ~ splines::ns(hp, 3) * wt, data = mtcars)
  surface <- term_predictions(f2, c("hp", "wt"), n = 12)
  res <- test_simultaneous_predictions(surface, function(hp, wt) 20 + 0 * hp)
  layers <- vapply(plot(res)$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomRaster" %in% layers)
})

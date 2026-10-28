fit <- lm(mpg ~ splines::ns(hp, 3) * wt + qsec, data = mtcars)

test_that("term_heatmap builds a surface with contours, cuts and data", {
  pred <- term_predictions(fit, c("hp", "wt"), n = 20)
  p <- term_heatmap(pred, cuts = list(wt = c(2.5, 3.5), hp = 150), show_data = TRUE)
  expect_s3_class(p, "ggplot")
  layers <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_true(all(c("GeomRaster", "GeomContour", "GeomPoint", "GeomHline", "GeomVline") %in% layers))
  expect_match(p$labels$subtitle, "qsec = ")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("term_heatmap rejects non-surfaces and unknown cut variables", {
  expect_error(term_heatmap(term_predictions(fit, "hp")), "two continuous")
  pred <- term_predictions(fit, c("hp", "wt"), n = 5)
  expect_error(term_heatmap(pred, cuts = list(qsec = 18)), "cuts")
})

test_that("term_curves draws one curve per cut with the chosen bands", {
  pred <- term_predictions(fit, c("hp", "wt"), values = list(wt = c(2.5, 3.5)), n = 20)
  p <- term_curves(pred, band = "both")
  layers <- vapply(p$layers, function(l) class(l$geom)[1], character(1))
  expect_equal(sum(layers == "GeomRibbon"), 1)
  expect_equal(sum(layers == "GeomLine"), 3)
  expect_match(p$labels$caption, "simultaneous")
  built <- ggplot2::ggplot_build(p)
  expect_equal(length(unique(built$data[[1]]$group)), 2)

  expect_null(term_curves(pred, band = "none")$labels$caption)
  expect_error(term_curves(term_predictions(fit, c("hp", "wt"), n = 20)), "too many values")
})

test_that("plot() picks heatmap or curves", {
  surface <- term_predictions(fit, c("hp", "wt"), n = 15)
  cuts <- term_predictions(fit, c("hp", "wt"), values = list(wt = c(2.5, 3.5)), n = 15)
  heat_layers <- vapply(plot(surface)$layers, function(l) class(l$geom)[1], character(1))
  curve_layers <- vapply(plot(cuts)$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomRaster" %in% heat_layers)
  expect_true("GeomRibbon" %in% curve_layers)
})

test_that("outcome labels describe the scale", {
  g <- glm(vs ~ wt, family = binomial, data = mtcars)
  expect_equal(attr(term_predictions(g, "wt", n = 3), "outcome_label"), "Probability of vs")
  expect_equal(attr(term_predictions(g, "wt", n = 3, scale = "link"), "outcome_label"), "Linear predictor")
  expect_equal(attr(term_predictions(fit, "hp", n = 3), "outcome_label"), "Predicted mpg")
})

test_that("layers with different variables can be drawn together", {
  surface <- term_predictions(fit, c("hp", "wt"), n = 15)
  one <- term_layer(surface, list(hp = function(w) 100 + 50 * w, wt = function(w) 3 - w),
                    terms = "w", values = list(w = seq(0, 1, by = 0.1)))
  family <- term_layer(surface, list(hp = function(w, s) s + 50 * w, wt = function(w) 3 - w),
                       terms = c("w", "s"), values = list(w = seq(0, 1, by = 0.1), s = c(100, 150)))
  p <- term_heatmap(surface, layers = list(one, family))
  expect_length(unique(ggplot2::layer_data(p, length(p$layers))$group), 3)
})

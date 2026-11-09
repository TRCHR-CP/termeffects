test_that("term_figure bundles plots, covariates and a caption", {
  d <- transform(mtcars, am = factor(am))
  fit <- lm(mpg ~ splines::ns(hp, 3) * splines::ns(wt, 3) + qsec + am + gear + carb + drat, data = d)
  surface <- term_predictions(fit, c("hp", "wt"), n = 15)
  fig <- term_figure(surface, cuts = list(wt = c(2.5, 3.5)), show_data = TRUE)

  expect_s3_class(fig, "term_figure")
  expect_s3_class(fig$heatmap, "ggplot")
  expect_s3_class(fig$curves, "ggplot")
  expect_s3_class(fig$uncertainty, "ggplot")
  expect_equal(fig$covariates$covariate, c("qsec", "am", "gear", "carb", "drat"))
  expect_equal(fig$covariates$chosen_as[fig$covariates$covariate == "am"], "most common")
  expect_equal(fig$covariates$value[fig$covariates$covariate == "am"], "0")
  expect_match(fig$caption, "Predicted mpg over hp and wt")
  expect_match(fig$caption, "interact in the model")
  expect_match(fig$caption, "qsec = 17.7 (median)", fixed = TRUE)
  expect_match(fig$caption, "cuts through the surface at wt = 2.5, 3.5")
  expect_match(fig$model, "^linear regression, n = 32")
  # all held covariates are listed in the subtitle, without pointing elsewhere
  expect_match(fig$heatmap$labels$subtitle, "Held constant: qsec = 17.7, am = 0, gear = 4", fixed = TRUE)
  expect_no_match(fig$heatmap$labels$subtitle, "term_figure", fixed = TRUE)
  expect_output(print(fig), "Predicted mpg")
})

test_that("captions describe additive surfaces and orm models", {
  fit <- lm(mpg ~ hp + wt, data = mtcars)
  fig <- term_figure(term_predictions(fit, c("hp", "wt"), n = 5))
  expect_match(fig$caption, "In the model, hp enters linearly and wt linearly.", fixed = TRUE)
  expect_match(fig$caption, "no interaction")
  expect_match(fig$caption, "There are no other covariates")

  skip_if_not_installed("rms")
  o <- rms::orm(ejection_fraction ~ rms::rcs(creatinine, 4) + age, data = heartfailure)
  fig <- term_figure(term_predictions(o, c("creatinine", "age"), n = 5))
  expect_match(fig$model, "proportional odds model (rms::orm), n = 299", fixed = TRUE)
  expect_match(fig$caption, "Mean ejection_fraction")
})

test_that("Cox model descriptions give patients and events, not nobs()", {
  skip_if_not_installed("survival")
  hf <- heartfailure
  hf$died <- as.numeric(hf$death_event == "Yes")
  fit <- survival::coxph(survival::Surv(time, died) ~ creatinine + sodium, data = hf)
  expect_match(termeffects:::.model_description(fit), "n = 299, 96 events", fixed = TRUE)
})

test_that("heatmap caption counts observations outside the plotted range", {
  fit <- lm(mpg ~ splines::ns(hp, 3) * splines::ns(wt, 3), data = mtcars)
  full <- term_predictions(fit, c("hp", "wt"), n = 10)
  expect_null(term_heatmap(full, show_data = TRUE)$labels$caption)
  part <- term_predictions(fit, c("hp", "wt"), values = list(hp = seq(100, 200, length.out = 10)), n = 10)
  outside <- sum(mtcars$hp < 100 | mtcars$hp > 200)
  expect_equal(term_heatmap(part, show_data = TRUE)$labels$caption,
               paste(outside, "of 32 observations lie outside the plotted range"))
  expect_null(term_heatmap(part)$labels$caption)
})

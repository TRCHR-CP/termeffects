# Each supported model class is checked against the model's own predict()
# method for estimates and standard errors.

test_that("lm: estimates and SEs match predict(), bands are wider than pointwise", {
  fit <- lm(mpg ~ splines::ns(hp, 3) * wt + factor(am), data = mtcars)
  set.seed(1)
  pred <- term_predictions(fit, c("hp", "wt"), values = list(wt = c(2.5, 3.5)), n = 15)
  ref <- predict(fit, newdata = pred, se.fit = TRUE)

  expect_equal(pred$estimate, unname(ref$fit))
  expect_equal(pred$std.error, unname(ref$se.fit))
  expect_equal(attr(pred, "pointwise_critical_value"), qt(0.975, df.residual(fit)))
  expect_equal(pred$conf.low, unname(ref$fit - qt(0.975, df.residual(fit)) * ref$se.fit))

  crit <- attr(pred, "critical_value")
  bonferroni <- qt(1 - 0.025 / nrow(pred), df.residual(fit))
  expect_gt(crit, attr(pred, "pointwise_critical_value"))
  expect_lt(crit, bonferroni)
  expect_true(all(pred$band.low <= pred$conf.low & pred$band.high >= pred$conf.high))
})

test_that("glm: bands are built on the link scale and stay inside [0, 1]", {
  set.seed(3)
  d <- data.frame(x = runif(600), z = runif(600))
  d$y <- rbinom(600, 1, plogis(-1 + d$x + 2 * d$x * d$z))
  fit <- glm(y ~ x * z, family = binomial, data = d)
  pred <- term_predictions(fit, c("x", "z"), values = list(z = c(0.2, 0.5, 0.8)), n = 21)
  eta <- predict(fit, newdata = pred, se.fit = TRUE)
  crit <- attr(pred, "critical_value")

  expect_equal(pred$estimate, plogis(unname(eta$fit)))
  expect_equal(pred$std.error, unname(eta$se.fit))
  expect_equal(pred$band.low, plogis(unname(eta$fit - crit * eta$se.fit)))
  expect_equal(pred$band.high, plogis(unname(eta$fit + crit * eta$se.fit)))
  expect_true(all(pred$band.low >= 0 & pred$band.high <= 1))
  expect_equal(attr(pred, "df"), Inf)
  expect_gt(crit, qnorm(0.975))

  link <- term_predictions(fit, c("x", "z"), values = list(z = 0.5), n = 5, scale = "link")
  expect_equal(link$estimate, unname(predict(fit, newdata = link)))
})

test_that("lmer: matches population-level predict() and its covariance", {
  skip_if_not_installed("lme4")
  fit <- lme4::lmer(Reaction ~ splines::ns(Days, 3) + (1 | Subject), data = lme4::sleepstudy)
  pred <- term_predictions(fit, "Days", n = 10)
  expect_equal(pred$estimate, unname(predict(fit, newdata = pred, re.form = NA)))

  X <- model.matrix(~ splines::ns(Days, 3), data = rbind(lme4::sleepstudy["Days"], pred["Days"]))
  X <- X[-seq_len(nrow(lme4::sleepstudy)), ]
  expect_equal(pred$std.error, unname(sqrt(diag(X %*% as.matrix(vcov(fit)) %*% t(X)))))
})

test_that("glmer: response scale uses the inverse link", {
  skip_if_not_installed("lme4")
  fit <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                     family = binomial, data = lme4::cbpp)
  pred <- term_predictions(fit, "period")
  expect_equal(pred$estimate, unname(predict(fit, newdata = pred, re.form = NA, type = "response")))
  expect_true(all(pred$band.low >= 0 & pred$band.high <= 1))
})

test_that("coxph: matches predict(type = 'lp') and gives hazard ratios", {
  skip_if_not_installed("survival")
  lung2 <- na.omit(survival::lung[c("time", "status", "age", "ph.karno", "sex")])
  fit <- survival::coxph(survival::Surv(time, status) ~ splines::ns(age, 3) * ph.karno + sex,
                         data = lung2)
  pred <- term_predictions(fit, c("age", "ph.karno"), values = list(ph.karno = c(60, 80, 100)),
                           n = 10, scale = "link")
  ref <- predict(fit, newdata = pred, type = "lp", se.fit = TRUE)

  expect_equal(pred$estimate, unname(ref$fit))
  expect_equal(pred$std.error, unname(ref$se.fit))

  hr <- term_predictions(fit, "age", n = 5)
  expect_equal(hr$estimate, unname(predict(fit, newdata = hr, type = "risk")))
})

test_that("coxph with strata() gives a clear error", {
  skip_if_not_installed("survival")
  fit <- survival::coxph(survival::Surv(time, status) ~ age + survival::strata(sex),
                         data = survival::lung)
  expect_error(term_predictions(fit, "age"), "not supported")
})

test_that("mgcv gam: matches predict(se.fit = TRUE)", {
  skip_if_not_installed("mgcv")
  set.seed(2)
  d <- data.frame(x = runif(300), z = runif(300))
  d$y <- sin(2 * pi * d$x) * d$z + rnorm(300, sd = 0.3)
  fit <- mgcv::gam(y ~ te(x, z), data = d)
  pred <- term_predictions(fit, c("x", "z"), n = 8)
  ref <- predict(fit, newdata = pred, se.fit = TRUE)

  expect_equal(pred$estimate, unname(as.numeric(ref$fit)))
  expect_equal(pred$std.error, unname(as.numeric(ref$se.fit)))
})

test_that("rms::orm: matches predict(), Mean() and Mean()'s delta-method limits", {
  skip_if_not_installed("rms")
  rcs <- rms::rcs
  set.seed(1)
  d <- data.frame(x = runif(300), z = runif(300), g = factor(sample(c("a", "b"), 300, TRUE)))
  d$y <- round(pmin(100, pmax(0, 50 + 30 * sin(3 * d$x) * d$z + rnorm(300, sd = 15))))
  fit <- rms::orm(y ~ rcs(x, 4) * z + g, data = d)

  link <- term_predictions(fit, c("x", "z"), values = list(z = c(0.2, 0.8)), n = 10, scale = "link")
  ref <- predict(fit, newdata = link, se.fit = TRUE)
  expect_equal(link$estimate, unname(ref$linear.predictors))
  expect_equal(link$std.error, unname(ref$se.fit))

  pred <- term_predictions(fit, c("x", "z"), values = list(z = c(0.2, 0.8)), n = 10)
  mean_fun <- Hmisc::Mean(fit)
  expect_equal(pred$estimate, unname(mean_fun(ref$linear.predictors)))
  expect_true(all(pred$band.low >= min(d$y) & pred$band.high <= max(d$y)))
  expect_true(all(pred$band.low <= pred$conf.low & pred$band.high >= pred$conf.high))

  # pointwise limits are rms's delta method with all intercepts:
  # Mean(fit)(X = , conf.int = 0.95)
  X <- predict(fit, newdata = pred, type = "x")
  limits <- attr(mean_fun(X = X, conf.int = 0.95), "limits")
  expect_equal(pred$conf.low, unname(limits$lower), tolerance = 1e-6)
  expect_equal(pred$conf.high, unname(limits$upper), tolerance = 1e-6)
})

test_that("bootstrap intervals agree with analytic ones for a linear model", {
  fit <- lm(mpg ~ splines::ns(hp, 3) + wt, data = mtcars)
  set.seed(3)
  boot <- term_predictions(fit, "hp", n = 10, ci_method = "bootstrap", nboot = 400)
  analytic <- term_predictions(fit, "hp", n = 10)
  expect_equal(boot$estimate, analytic$estimate)
  # bootstrap SEs are noisy, and larger at the edges of a spline in a small
  # sample (32 cars), but of the same size overall
  ratio <- boot$std.error / analytic$std.error
  expect_true(median(ratio) > 0.7 && median(ratio) < 1.5)
  expect_equal(attr(boot, "ci_method"), "bootstrap")
  expect_lte(attr(boot, "nboot"), 400)
  expect_gt(attr(boot, "critical_value"), attr(boot, "pointwise_critical_value"))
  expect_output(print(boot), "bootstrap")
})

test_that("bootstrap intervals for orm use bootcov and need x, y", {
  skip_if_not_installed("rms")
  set.seed(4)
  o <- rms::orm(ejection_fraction ~ rms::rcs(creatinine, 3) + age, data = heartfailure)
  expect_error(term_predictions(o, "creatinine", ci_method = "bootstrap", nboot = 10), "x = TRUE")
  o <- rms::orm(ejection_fraction ~ rms::rcs(creatinine, 3) + age, data = heartfailure, x = TRUE, y = TRUE)
  boot <- term_predictions(o, "creatinine", n = 8, ci_method = "bootstrap", nboot = 50)
  expect_true(all(boot$band.low <= boot$estimate & boot$estimate <= boot$band.high))
  expect_true(all(boot$band.low >= 14 & boot$band.high <= 80))
})

test_that("bootstrap is refused for mixed models", {
  skip_if_not_installed("lme4")
  m <- lme4::lmer(Reaction ~ Days + (1 | Subject), data = lme4::sleepstudy)
  expect_error(term_predictions(m, "Days", ci_method = "bootstrap"), "mixed models")
})

test_that("ci_scope gives one critical value per facet or curve", {
  fit <- lm(mpg ~ hp * wt * factor(am), data = mtcars)
  set.seed(5)
  values <- list(wt = c(2.5, 3.5))
  global <- term_predictions(fit, c("hp", "wt", "am"), values = values, n = 10)
  faceted <- term_predictions(fit, c("hp", "wt", "am"), values = values, n = 10, ci_scope = "faceted")
  curve <- term_predictions(fit, c("hp", "wt", "am"), values = values, n = 10, ci_scope = "curve")
  expect_length(attr(global, "critical_value"), 1)
  expect_length(attr(faceted, "critical_value"), 2)
  expect_length(attr(curve, "critical_value"), 4)
  # narrower scopes need smaller multipliers
  expect_true(all(attr(curve, "critical_value") < attr(global, "critical_value")))
})

test_that("newdata can trace an arbitrary curve through the surface", {
  fit <- lm(mpg ~ hp * wt, data = mtcars)
  t <- seq(0, 1, length.out = 20)
  curve <- data.frame(hp = 100 + 100 * t, wt = 2 + t^2)
  pred <- term_predictions(fit, newdata = curve)
  expect_equal(pred$estimate, unname(predict(fit, newdata = curve)))
  expect_equal(nrow(pred), 20)
})

test_that("unsupported models and offsets give clear errors", {
  fit <- lm(mpg ~ hp + offset(wt), data = mtcars)
  expect_error(term_predictions(fit, "hp"), "offsets")
  fit <- nls(mpg ~ a * exp(b * hp), data = mtcars, start = list(a = 40, b = -0.005))
  expect_error(term_predictions(fit, "hp", data = mtcars), "not supported")
})

test_that("print shows the critical values and held values", {
  fit <- lm(mpg ~ hp * wt + qsec, data = mtcars)
  set.seed(1)
  expect_output(print(term_predictions(fit, c("hp", "wt"), n = 3)), "simultaneous")
  expect_output(print(term_predictions(fit, c("hp", "wt"), n = 3)), "qsec = ")
})

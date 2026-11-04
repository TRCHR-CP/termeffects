# Each model type: termeffects' estimates and standard errors (on the link
# scale) must equal those of the model's own predict() method or of X V X'.

hf <- heartfailure
hf$died <- as.numeric(hf$death_event == "Yes")
set.seed(11)
hf$site <- factor(sample(1:8, nrow(hf), replace = TRUE))
values <- list(creatinine = c(0.8, 1.5, 3), sodium = c(130, 140))

link_predictions <- function(model, ...) {
  term_predictions(model, c("creatinine", "sodium"), values = values, scale = "link", ...)
}
grid_data <- function(pred) {
  d <- as.data.frame(pred)
  d[setdiff(names(d), c("estimate", "std.error", "conf.low", "conf.high", "band.low", "band.high"))]
}

test_that("Poisson GLM with an offset: offset variables are held, or set with `at`", {
  fit <- glm(died ~ splines::ns(creatinine, 3) + sodium + offset(log(time)), poisson, data = hf)
  pred <- link_predictions(fit)
  expect_equal(attr(pred, "held")$time, median(hf$time))
  own <- predict(fit, grid_data(pred), se.fit = TRUE)
  expect_equal(pred$estimate, unname(own$fit), ignore_attr = TRUE)
  expect_equal(pred$std.error, unname(own$se.fit), ignore_attr = TRUE)
  # rate per year
  per_year <- term_predictions(fit, c("creatinine", "sodium"), values = values, at = list(time = 365))
  expect_equal(per_year$estimate, exp(pred$estimate) * 365 / median(hf$time))
})

test_that("an `offset =` argument is used as well", {
  fit <- glm(died ~ splines::ns(creatinine, 3) + sodium, poisson, data = hf, offset = log(time))
  pred <- link_predictions(fit)
  expect_equal(pred$estimate, unname(predict(fit, grid_data(pred))), ignore_attr = TRUE)
})

test_that("stratified Cox model: hazard ratios as for an unstratified prediction", {
  skip_if_not_installed("survival")
  fit <- survival::coxph(survival::Surv(time, died) ~ splines::ns(creatinine, 3) + sodium +
                           survival::strata(sex), data = hf)
  pred <- link_predictions(fit)
  X <- model.matrix(fit)
  X_new <- termeffects:::.design(fit, grid_data(pred))$X
  expect_equal(colnames(X_new), names(coef(fit)))
  expect_equal(pred$estimate, drop(X_new %*% coef(fit)), ignore_attr = TRUE)
  expect_equal(pred$std.error, sqrt(rowSums((X_new %*% vcov(fit)) * X_new)), ignore_attr = TRUE)
  expect_equal(pred$estimate, unname(predict(fit, grid_data(pred), type = "lp", reference = "sample")), ignore_attr = TRUE)
})

test_that("Cox model with time-dependent covariates (tmerge data, as in survival's timedep vignette)", {
  skip_if_not_installed("survival")
  pbc <- survival::pbc
  base <- subset(pbc, id <= 312, select = c(id:sex, stage))
  pbc2 <- survival::tmerge(base, base, id = id, death = event(time, status))
  pbc2 <- survival::tmerge(pbc2, survival::pbcseq, id = id, bili = tdc(day, bili),
                           albumin = tdc(day, albumin))
  fit <- survival::coxph(survival::Surv(tstart, tstop, death == 2) ~ splines::ns(log(bili), 3) *
                           albumin + age, data = pbc2, cluster = id)
  pred <- term_predictions(fit, c("bili", "albumin"), values = list(bili = c(0.5, 2, 8), albumin = c(3, 4)),
                           scale = "link")
  own <- predict(fit, grid_data(pred), type = "lp", se.fit = TRUE, reference = "sample")
  expect_equal(pred$estimate, own$fit, ignore_attr = TRUE)
  expect_equal(pred$std.error, own$se.fit, ignore_attr = TRUE)
  # time variables and the event are part of the outcome, not covariates
  expect_setequal(names(attr(pred, "held")), "age")
  expect_error(term_predictions(fit, "bili", ci_method = "bootstrap"), "clustered")
})

test_that("Firth logistic regression (logistf and brglm2) agrees with the model", {
  skip_if_not_installed("logistf")
  fit <- logistf::logistf(died ~ splines::ns(creatinine, 3) + sodium + sex, data = hf)
  pred <- link_predictions(fit)
  X <- termeffects:::.design(fit, grid_data(pred))$X
  expect_equal(pred$estimate, drop(X %*% coef(fit)), ignore_attr = TRUE)
  expect_equal(pred$estimate, unname(predict(fit, grid_data(pred), type = "link")), ignore_attr = TRUE)
  expect_equal(pred$std.error, sqrt(rowSums((X %*% vcov(fit)) * X)), ignore_attr = TRUE)
  expect_equal(attr(term_predictions(fit, "creatinine", n = 3), "outcome_label"), "Probability of died")

  skip_if_not_installed("brglm2")
  br <- glm(died ~ splines::ns(creatinine, 3) + sodium + sex, binomial, data = hf, method = brglm2::brglmFit)
  own <- predict(br, grid_data(link_predictions(br)), se.fit = TRUE)
  expect_equal(link_predictions(br)$estimate, unname(own$fit), ignore_attr = TRUE)
  expect_equal(link_predictions(br)$std.error, unname(own$se.fit), ignore_attr = TRUE)
  # both are Firth's method, so the two packages agree
  expect_equal(link_predictions(br)$estimate, pred$estimate, tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("nlme::lme gives population-level predictions", {
  skip_if_not_installed("nlme")
  fit <- nlme::lme(ejection_fraction ~ splines::ns(creatinine, 3) + sodium + sex,
                   random = ~ 1 | site, data = hf)
  pred <- link_predictions(fit)
  expect_equal(pred$estimate, unname(predict(fit, grid_data(pred), level = 0)), ignore_attr = TRUE)
  X <- termeffects:::.design(fit, grid_data(pred))$X
  expect_equal(pred$std.error, sqrt(rowSums((X %*% vcov(fit)) * X)), ignore_attr = TRUE)
  expect_error(term_predictions(fit, "creatinine", ci_method = "bootstrap"), "clustered")
})

test_that("lme4::glmer gives population-level predictions", {
  skip_if_not_installed("lme4")
  fit <- suppressWarnings(suppressMessages(lme4::glmer(died ~ splines::ns(creatinine, 3) + sodium + (1 | site),
                                      data = hf, family = binomial)))
  pred <- link_predictions(fit)
  expect_equal(pred$estimate, unname(predict(fit, grid_data(pred), re.form = NA)), ignore_attr = TRUE)
})

test_that("GEE (geepack) uses the robust covariance and normal intervals", {
  skip_if_not_installed("geepack")
  d <- hf[order(hf$site), ]
  fit <- geepack::geeglm(died ~ splines::ns(creatinine, 3) + sodium, binomial, data = d,
                         id = site, corstr = "exchangeable")
  pred <- link_predictions(fit)
  X <- termeffects:::.design(fit, grid_data(pred))$X
  expect_equal(pred$estimate, unname(predict(fit, grid_data(pred))), ignore_attr = TRUE)
  expect_equal(pred$std.error, sqrt(rowSums((X %*% vcov(fit)) * X)), ignore_attr = TRUE)
  expect_equal(attr(pred, "pointwise_critical_value"), qnorm(0.975))
  expect_error(term_predictions(fit, "creatinine", ci_method = "bootstrap"), "clustered")
})

test_that("rms lrm, ols and cph agree with rms's predict()", {
  skip_if_not_installed("rms")
  fits <- list(
    lrm = rms::lrm(died ~ rms::rcs(creatinine, 4) + sodium + sex, data = hf, x = TRUE, y = TRUE),
    ols = rms::ols(ejection_fraction ~ rms::rcs(creatinine, 4) + sodium + sex, data = hf, x = TRUE, y = TRUE),
    cph = rms::cph(survival::Surv(time, died) ~ rms::rcs(creatinine, 4) + sodium + sex, data = hf,
                   x = TRUE, y = TRUE)
  )
  for (fit in fits) {
    pred <- link_predictions(fit)
    own <- predict(fit, grid_data(pred), type = "lp", se.fit = TRUE)
    expect_equal(pred$estimate, unname(own$linear.predictors), ignore_attr = TRUE)
    expect_equal(pred$std.error, unname(own$se.fit), ignore_attr = TRUE)
  }
  expect_equal(attr(term_predictions(fits$lrm, "creatinine", n = 3), "outcome_label"), "Probability of died")
  expect_true(attr(term_predictions(fits$cph, "creatinine", n = 3), "log_scale"))
  hf$ef3 <- cut(hf$ejection_fraction, c(0, 30, 45, 100))
  expect_error(term_predictions(rms::lrm(ef3 ~ creatinine, data = hf), "creatinine"), "orm")
})

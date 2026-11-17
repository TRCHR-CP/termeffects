# Independent check of the simultaneous bands: rms::Predict(conf.type =
# "simultaneous") computes them with multcomp (numerical integration of the
# multivariate t / normal distribution), termeffects by simulation. On the
# same model and grid, the estimates and standard errors must be identical and
# the critical values agree up to simulation error. Small grids only: rms's
# simultaneous bands fail for larger grids (see vignette("coverage_study")).

rms_vs_termeffects <- function(fit) {
  creatinine <- seq(0.6, 3, length.out = 10)
  sodium <- c(130, 137, 142)
  simultaneous <- rms::Predict(fit, creatinine = creatinine, sodium = sodium, conf.type = "simultaneous")
  pointwise <- rms::Predict(fit, creatinine = creatinine, sodium = sodium)
  set.seed(1)
  te <- term_predictions(fit, c("creatinine", "sodium"),
                         values = list(creatinine = creatinine, sodium = sodium),
                         scale = "link", nsim = 1e5)
  by <- c("creatinine", "sodium")
  m <- merge(merge(as.data.frame(simultaneous)[c(by, "yhat", "upper")],
                   as.data.frame(pointwise)[c(by, "upper")], by = by, suffixes = c("", "_pw")),
             as.data.frame(te), by = by)
  se_rms <- (m$upper_pw - m$yhat) / attr(te, "pointwise_critical_value")
  list(estimate = max(abs(m$yhat - m$estimate)),
       se = max(abs(se_rms / m$std.error - 1)),
       crit_rms = mean((m$upper - m$yhat) / se_rms),
       crit_te = attr(te, "critical_value"))
}

test_that("simultaneous bands agree with rms::Predict(conf.type = 'simultaneous')", {
  skip_on_cran()
  skip_if_not_installed("rms")
  skip_if_not_installed("multcomp")
  suppressMessages(library(rms))
  hf <- heartfailure
  hf$death <- as.numeric(hf$death_event == "Yes")
  dd <- rms::datadist(hf)
  old <- options(datadist = dd)
  on.exit(options(old))

  fits <- list(
    ols = ols(ejection_fraction ~ rcs(creatinine, 4) * rcs(sodium, 4) + rcs(age, 4) + sex + diabetes,
              data = hf, x = TRUE, y = TRUE),
    lrm = lrm(death ~ rcs(creatinine, 4) * rcs(sodium, 3) + rcs(age, 3) + sex,
              data = hf, x = TRUE, y = TRUE),
    cph = cph(Surv(time, death) ~ rcs(creatinine, 4) * rcs(sodium, 3) + rcs(age, 3) + sex,
              data = hf, x = TRUE, y = TRUE)
  )
  for (name in names(fits)) {
    r <- rms_vs_termeffects(fits[[name]])
    expect_lt(r$estimate, 1e-8, label = paste(name, "estimates"))
    expect_lt(r$se, 1e-8, label = paste(name, "standard errors"))
    expect_lt(abs(r$crit_te / r$crit_rms - 1), 0.01, label = paste(name, "critical value"))
  }
})

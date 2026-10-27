test_that("focal terms vary, other predictors are held at typical values", {
  d <- transform(mtcars, am = factor(am, labels = c("auto", "manual")))
  fit <- lm(mpg ~ hp * wt + am + qsec, data = d)
  grid <- term_grid(fit, c("hp", "wt"), n = 4)

  expect_equal(nrow(grid), 16)
  expect_equal(range(grid$hp), range(d$hp))
  expect_equal(unique(grid$qsec), mean(d$qsec))
  expect_equal(as.character(unique(grid$am)), "auto")
  expect_equal(levels(grid$am), c("auto", "manual"))
  expect_equal(attr(grid, "focal"), c("hp", "wt"))
  expect_named(attr(grid, "held"), c("am", "qsec"))
})

test_that("values and at override the defaults", {
  d <- transform(mtcars, am = factor(am, labels = c("auto", "manual")))
  fit <- lm(mpg ~ hp * wt + am, data = d)
  grid <- term_grid(fit, c("hp", "wt"), values = list(wt = c(2, 3)), at = list(am = "manual"), n = 5)

  expect_equal(sort(unique(grid$wt)), c(2, 3))
  expect_equal(as.character(unique(grid$am)), "manual")
  expect_equal(nrow(grid), 10)
})

test_that("factor focal terms use all levels", {
  d <- transform(mtcars, cyl = factor(cyl))
  fit <- lm(mpg ~ hp * cyl, data = d)
  grid <- term_grid(fit, c("hp", "cyl"), n = 3)
  expect_equal(levels(grid$cyl), c("4", "6", "8"))
  expect_equal(nrow(grid), 9)
})

test_that("helpful errors for misspelled variables and bad levels", {
  d <- transform(mtcars, am = factor(am))
  fit <- lm(mpg ~ hp + am, data = d)
  expect_error(term_grid(fit, "horsepower"), "not found")
  expect_error(term_grid(fit, "hp", at = list(am = "2")), "not levels")
  expect_error(term_grid(fit, "hp", at = list(hp = 100)), "Focal variables")
})

test_that("data can be supplied when it cannot be recovered from the call", {
  # with a spline, the raw `hp` is not in the model frame either
  fit <- local({
    dd <- mtcars
    lm(mpg ~ splines::ns(hp, 2), data = dd)
  })
  rm_env <- new.env()
  environment(fit$terms) <- rm_env
  expect_error(term_grid(fit, "hp", data = NULL), "Could not recover")
  expect_s3_class(term_grid(fit, "hp", data = mtcars), "term_grid")
})

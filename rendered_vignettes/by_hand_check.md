# Checking termeffects by hand

This document recomputes what `term_predictions()` returns, step by
step, with plain R and without any termeffects code, and compares every
number. It checks the method rather than the code: if every step agrees,
the package computes what it says it computes. Run it chunk by chunk in
RStudio (`validation/by_hand_check.Rmd`), or read the rendered version
in `rendered_vignettes/by_hand_check.md`.

The steps are the ones every model in termeffects goes through:

1.  the prediction grid: the focal variables vary, the others are held
    at their median or most common value;
2.  the design matrix $`X`$ of the grid, so that the predictions are
    $`X\hat\beta`$;
3.  the standard errors,
    $`\mathrm{SE} = \sqrt{\operatorname{diag}(X V X^\top)}`$ with $`V`$
    = `vcov(fit)`;
4.  the pointwise intervals, $`X\hat\beta \pm t_{0.975}\,\mathrm{SE}`$;
5.  the simultaneous critical value $`c`$: the 95% quantile of the
    largest $`|t|`$ over the grid, where
    $`t = X(\beta^* - \hat\beta)/\mathrm{SE}`$ for coefficients
    $`\beta^*`$ drawn from their estimated distribution;
6.  the band, $`X\hat\beta \pm c\,\mathrm{SE}`$, and for a logistic
    model its transformation to probabilities.

Part 1 does this for a linear model, part 2 for a logistic model.

## Part 1: a linear model

``` r
library(termeffects)
library(rms)        # for rcs()

hf <- heartfailure
fit <- lm(ejection_fraction ~ rcs(creatinine, 4) * sodium + rcs(age, 3) + sex, data = hf)
```

### What termeffects returns

Three curves: creatinine from 0.6 to 3 at sodium 130, 137 and 142, with
100,000 simulation draws for the critical value.

``` r
set.seed(1)
te <- term_predictions(fit, c("creatinine", "sodium"),
                       values = list(creatinine = seq(0.6, 3, length.out = 25),
                                     sodium = c(130, 137, 142)),
                       nsim = 1e5)
te <- as.data.frame(te)   # keeps the attributes we need below
head(te)
#>   creatinine sodium age  sex estimate std.error conf.low conf.high band.low
#> 1        0.6    130  60 Male 38.20247  6.158318 26.08145  50.32349 20.14104
#> 2        0.7    130  60 Male 37.93629  4.301968 29.46901  46.40358 25.31926
#> 3        0.8    130  60 Male 37.65227  2.795191 32.15067  43.15386 29.45439
#> 4        0.9    130  60 Male 37.26116  2.216619 32.89833  41.62398 30.76014
#> 5        1.0    130  60 Male 36.65588  2.358972 32.01287  41.29890 29.73737
#> 6        1.1    130  60 Male 35.78324  2.294177 31.26776  40.29872 29.05476
#>   band.high
#> 1  56.26390
#> 2  50.55332
#> 3  45.85014
#> 4  43.76217
#> 5  43.57440
#> 6  42.51172
```

### Step 1: the grid

``` r
most_common <- function(x) names(which.max(table(x)))
grid <- expand.grid(creatinine = seq(0.6, 3, length.out = 25), sodium = c(130, 137, 142))
grid$age <- median(hf$age)
grid$sex <- factor(most_common(hf$sex), levels = levels(hf$sex))
head(grid)
#>   creatinine sodium age  sex
#> 1        0.6    130  60 Male
#> 2        0.7    130  60 Male
#> 3        0.8    130  60 Male
#> 4        0.9    130  60 Male
#> 5        1.0    130  60 Male
#> 6        1.1    130  60 Male

# same points, same held values?
all.equal(grid$creatinine, te$creatinine)
#> [1] TRUE
all.equal(grid$sodium, te$sodium)
#> [1] TRUE
c(age = unique(te$age), sex = as.character(unique(te$sex)))
#>    age    sex 
#>   "60" "Male"
```

### Step 2: design matrix and predictions

`model.matrix()` builds the spline columns with the knots stored in the
fitted model, so the grid gets the same basis as the data.

``` r
X <- model.matrix(delete.response(terms(fit)), data = grid, xlev = fit$xlevels)
beta <- coef(fit)
estimate <- drop(X %*% beta)

# also the same as predict()
max(abs(estimate - predict(fit, newdata = grid)))
#> [1] 0
# difference from termeffects
max(abs(estimate - te$estimate))
#> [1] 0
```

### Step 3: standard errors

``` r
V <- vcov(fit)
se <- sqrt(diag(X %*% V %*% t(X)))

max(abs(se - te$std.error))
#> [1] 9.769963e-14
```

### Step 4: pointwise intervals

A linear model has residual degrees of freedom, so the quantile is from
the t distribution.

``` r
df <- df.residual(fit)
t_pointwise <- qt(0.975, df)
t_pointwise
#> [1] 1.968235

max(abs((estimate - t_pointwise * se) - te$conf.low))
#> [1] 1.918465e-13
max(abs((estimate + t_pointwise * se) - te$conf.high))
#> [1] 1.98952e-13
```

### Step 5: the simultaneous critical value

Draw coefficients $`\beta^*`$ from their estimated distribution, a
multivariate t with the model’s degrees of freedom:
$`\beta^* = \hat\beta + R^\top z /
\sqrt{w/\mathit{df}}`$, with $`R`$ the Cholesky factor of $`V`$, $`z`$
standard normal and $`w \sim \chi^2_{\mathit{df}}`$. For each draw,
compute the standardized deviation at every grid point and keep the
largest one. The critical value is the 95% quantile of these maxima.

``` r
set.seed(2)
nsim <- 1e5
R <- chol(V)
z <- matrix(rnorm(nsim * length(beta)), nsim)
w <- rchisq(nsim, df)
deviation <- (z %*% R) / sqrt(w / df)          # beta* - beta-hat, one row per draw
t_draws <- (deviation %*% t(X)) / rep(se, each = nsim)
max_t <- apply(abs(t_draws), 1, max)
c_by_hand <- quantile(max_t, 0.95, names = FALSE)

c(by_hand = c_by_hand, termeffects = attr(te, "critical_value"))
#>     by_hand termeffects 
#>    2.939793    2.932850
```

The two values come from different random draws, so they differ a
little. How much is acceptable? With 100,000 draws, the 95% confidence
interval for the true 95% quantile is given by two order statistics of
the simulated maxima:

``` r
k <- nsim * 0.95 + c(-1.96, 1.96) * sqrt(nsim * 0.95 * 0.05)
sort(max_t)[round(k)]
#> [1] 2.929771 2.949910
```

termeffects’ critical value should lie in, or very close to, this
interval (its own simulation has a similar error). For comparison, a
Bonferroni correction over the 75 points would use 3.44: the
simultaneous band is much narrower because neighbouring points are
strongly correlated.

### Step 6: the band

The band uses termeffects’ critical value, so this checks the
arithmetic:

``` r
c_te <- attr(te, "critical_value")
max(abs((estimate - c_te * se) - te$band.low))
#> [1] 2.877698e-13
max(abs((estimate + c_te * se) - te$band.high))
#> [1] 2.913225e-13
```

## Part 2: a logistic model

For a logistic model everything is computed on the log-odds scale, and
the estimate and the limits are then transformed to probabilities. There
are no residual degrees of freedom, so the normal distribution replaces
the t.

``` r
hf$death <- as.numeric(hf$death_event == "Yes")
fit2 <- glm(death ~ rcs(creatinine, 4) * sodium + rcs(age, 3) + sex, family = binomial, data = hf)

set.seed(3)
te2 <- as.data.frame(term_predictions(fit2, c("creatinine", "sodium"),
                                      values = list(creatinine = seq(0.6, 3, length.out = 25),
                                                    sodium = c(130, 137, 142)),
                                      nsim = 1e5))

X2 <- model.matrix(delete.response(terms(fit2)), data = grid, xlev = fit2$xlevels)
eta <- drop(X2 %*% coef(fit2))
se2 <- sqrt(diag(X2 %*% vcov(fit2) %*% t(X2)))

# estimate: probability = plogis(log-odds)
max(abs(plogis(eta) - te2$estimate))
#> [1] 1.110223e-16
# std.error is reported on the log-odds scale, where the band is built
max(abs(se2 - te2$std.error))
#> [1] 2.020606e-14
# pointwise limits: plogis(eta -/+ 1.96 SE)
max(abs(plogis(eta - qnorm(0.975) * se2) - te2$conf.low))
#> [1] 9.603429e-15
max(abs(plogis(eta + qnorm(0.975) * se2) - te2$conf.high))
#> [1] 6.661338e-15
```

The critical value, by hand with normal draws:

``` r
set.seed(4)
R2 <- chol(vcov(fit2))
z <- matrix(rnorm(nsim * length(coef(fit2))), nsim)
t_draws2 <- ((z %*% R2) %*% t(X2)) / rep(se2, each = nsim)
max_t2 <- apply(abs(t_draws2), 1, max)
c(by_hand = quantile(max_t2, 0.95, names = FALSE), termeffects = attr(te2, "critical_value"))
#>     by_hand termeffects 
#>    2.923809    2.929398
sort(max_t2)[round(k)]   # 95% interval for the by-hand value
#> [1] 2.914407 2.934756

c2 <- attr(te2, "critical_value")
max(abs(plogis(eta - c2 * se2) - te2$band.low))
#> [1] 1.271205e-14
max(abs(plogis(eta + c2 * se2) - te2$band.high))
#> [1] 7.549517e-15
```

## Result

| Step                                             | Difference |
|:-------------------------------------------------|:-----------|
| 1 grid                                           | 0.0e+00    |
| 2 estimates (linear)                             | 0.0e+00    |
| 3 standard errors (linear)                       | 9.8e-14    |
| 4 pointwise (linear)                             | 2.0e-13    |
| 5 critical value (linear), relative difference   | 2.4e-03    |
| 6 band (linear)                                  | 2.9e-13    |
| 2 estimates (logistic)                           | 1.1e-16    |
| 3 standard errors (logistic)                     | 2.0e-14    |
| 4 pointwise (logistic)                           | 9.6e-15    |
| 5 critical value (logistic), relative difference | 1.9e-03    |
| 6 band (logistic)                                | 1.3e-14    |

All differences except the critical values are rounding error (around
1e-12 or smaller). The critical values differ only by simulation error,
well below 1%, and termeffects’ value lies within the 95% interval of
the by-hand simulation (step 5 of each part). Other models (Cox, mixed,
orm) use the same steps with their own design matrix and covariance;
`vignette("coverage_study")` checks them against the truth in
simulations and against `rms::Predict()`.

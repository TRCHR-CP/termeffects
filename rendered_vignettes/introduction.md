# termeffects: introduction

> **About this vignette.** Written by Lukas Bolsinger for p-ggeffects
> (original file `introduction_p_ggeffects.Rmd`) and ported to
> termeffects by Brigitte Mueller. His text is kept; the code now uses
> the termeffects functions (`ggpredict()` became `term_predictions()`,
> `plot(..., "heatmap")` became `term_heatmap()`, layers became
> `term_layer()`), and held covariates are at their median or most
> common value. Also changed: the description of the methods and of the
> differences from ggeffects, and one example that called
> `test_predictions()` (which failed on parameterized models) is
> removed. Natural splines (`ns()`) are replaced by restricted cubic
> splines (`rcs()`), as used in our team. The model with `%ia%` is now
> fitted with `ols()` instead of `lm()`: in `lm()`, `%ia%` recomputes
> the spline knots from the new data at prediction time, so the
> predictions (and the earlier figures) were wrong. Every change is
> marked in `validation/output/vignette_track_changes.html` (struck
> through: removed; underlined: added), or run
> `git diff 2de18bc -- vignettes/introduction.Rmd`.

# Adjusted predictions

After fitting a model, it is useful to generate model-based estimates
(expected values, or adjusted predictions) of the response variable for
different cobinations of predictor values. Such estimates can be used to
make inferences about relationships between variables–adjusted
precitions tell you: what is the expected outcome for certain values or
levels of my predictors? For complex models, visualizations of adjusted
predictions are often far easier to understand and intuititively
demonstrate associations of predictors and outcomes.

## Methods

`termeffects` computes adjusted predictions at representative values:
the predictors of interest vary over a grid, and all other predictors
are held at their median (numeric) or most common value (categorical),
or at values chosen with `at`. The result is returned as a data frame
with consistent structure, with a view towards ggplot plotting.

We call predictors of interest (ie. the non-marginalized terms) focal
terms.

# From p-ggeffects to termeffects

`termeffects` began as `p-ggeffects`, a ‘generalization’ of the
`ggeffects` package with a focus on non-linear interactions, and
non-trivial marginalizations. Particularly, the deviations from
`ggeffects` are focused on the introduction of `parameterized_model`
objects, regular grid generation for the purpose of 2d-heatmap plotting,
and simultaneous inference. `termeffects` keeps these features but no
longer depends on `ggeffects`: it computes the predictions and their
joint uncertainty itself.

## Parameterized models

A `parameterized_model` is a data structure consisting of a tuple of a
fitted regression model, a map into the variable space of the
model–particularly, it is viewed as a continuous parameterization of the
continuous predictor variables of the models. Moreover, the
`parameterized_model` structure tracks various auxiliary data.

We view parameterized models as a ‘generalization’ of generalized linear
models. Consider the example $`y \sim f_\beta(x) + \varepsilon`$, where
$`f_\beta`$ is the function $`x \mapsto \beta^T x`$. Call the domain of
interest of $`f_\beta`$, $`\Omega'`$. For a mapping
$`\phi: \Omega \to \Omega'`$, by composition, we obtain a ‘model’
$`y \sim f_\beta \circ \phi (z) + \varepsilon`$. This is the classical
example of a generalized additive model–the ‘pullback’ of the model
$`f_\beta`$ by the basis functions $`\phi`$. Notice however, that a data
set $`(Y, X)`$ in $`\Omega'`$ does not in general induce an estimate of
$`\beta`$ (viewed as the coefficients of $`f\circ \phi`$). Indeed, the
points of $`X`$ may not lie in $`\phi(\Omega)`$, and even if so, there
may not be a unique preimage.

In this light, rather than as a regression model, we view parameterized
models instead as generalizations of marginalization by taking
$`x_i = c`$. Indeed, taking $`f_\beta: \mathbb R^n \to \mathbb R`$ as
above, the describe marginalization is equivalent to the
parameterization $`\phi: \mathbb R^{n-1} \to \mathbb R^n`$ given by
$`\phi(x_1, \ldots, x_{n-1}) = (x_1, \ldots, x_{i-1}, c, x_i, \ldots, x_{n-1})`$.
Following the philosophy of the `ggeffects` package, for a more complex
model $`f_\beta`$, we may visualize and interpret the effect of a
(non-primary) variable of interest by plotting its marginalization at a
finite (\>1) set of values.

### Basic methods

We demonstrate some of the basic methods of `parameterized_model`
objects. Consider some non-linear regression model
$`y \sim f_\beta(w, h, \ldots) + \varepsilon`$ dependent on weight and
height mesurements. Consider the parameterization
$`(c, t, \ldots) \mapsto (w = ct^2, h = t, \ldots)`$. Fixing $`c`$ (and
the other variables, $`\ldots`$), we obtain a 1 parameteric family of
observations with fixed BMI and varying height. Doing similarly with
$`t`$, we parameterize the observations with varying weight/BMI and
fixed height.

``` r
w <- rnorm(500, mean = 70, sd = 4)
h <- rnorm(500, mean = 170, sd = 10)

a <- rnorm(500)

y <- a + w/h^2 + h + w + w^2 * h + rnorm(500)


# model example
model <- ols(y ~ rcs(w, 4) + rcs(h, 4) + rcs(w, 4) %ia% rcs(h, 4) + a)

# parameterize model
pmodel1 <- parameterize(model, parameterizations = list(w = \(t,c) c*t^2, h = \(t) t))

# parameterized model attributes
print(pmodel1, show_history = TRUE)
#> Parameterized model
#> 
#> Model:
#> Linear Regression Model
#> 
#> ols(formula = y ~ rcs(w, 4) + rcs(h, 4) + rcs(w, 4) %ia% rcs(h, 
#>     4) + a)
#> 
#>                    Model Likelihood    Discrimination    
#>                          Ratio Test           Indexes    
#> Obs       500    LR chi2    5109.08    R2       1.000    
#> sigma643.6014    d.f.            12    R2 adj   1.000    
#> d.f.      487    Pr(> chi2)  0.0000    g   118790.029    
#> 
#> Residuals
#> 
#>      Min       1Q   Median       3Q      Max 
#> -1195.87  -422.48    50.51   354.88  7443.82 
#> 
#> 
#>           Coef       S.E.       t      Pr(>|t|)
#> Intercept 13307.3134 33589.2033   0.40 0.6921  
#> w           -70.6463   505.9989  -0.14 0.8890  
#> w'        -3747.8989  1579.5233  -2.37 0.0180  
#> w''       13158.8003  5860.1185   2.25 0.0252  
#> h         -4231.0573   205.4127 -20.60 <0.0001 
#> h'         2049.5490   487.8048   4.20 <0.0001 
#> h''       -9695.6103  2111.5310  -4.59 <0.0001 
#> w * h       129.3205     3.0807  41.98 <0.0001 
#> w * h'      -28.7097     7.0092  -4.10 <0.0001 
#> w * h''     136.2667    30.4018   4.48 <0.0001 
#> w' * h       41.0340     9.2103   4.46 <0.0001 
#> w'' * h    -101.5784    34.2189  -2.97 0.0031  
#> a           -23.1594    28.3581  -0.82 0.4145  
#> 
#> 
#> Parameters:
#> [1] "t" "c"
#> 
#> Parameterizations:
#> $w <- (c * t^2); (t, c)
#> $h <- (t); (t)
#> 
#> 
#> Parameterization history:
#> 
#> Stage 1        
#> ---------------
#> 
#> w [c * t^2] <── t, c
#> h [t] <── t

# underlying model
original_model(pmodel1)
#> Linear Regression Model
#> 
#> ols(formula = y ~ rcs(w, 4) + rcs(h, 4) + rcs(w, 4) %ia% rcs(h, 
#>     4) + a)
#> 
#>                    Model Likelihood    Discrimination    
#>                          Ratio Test           Indexes    
#> Obs       500    LR chi2    5109.08    R2       1.000    
#> sigma643.6014    d.f.            12    R2 adj   1.000    
#> d.f.      487    Pr(> chi2)  0.0000    g   118790.029    
#> 
#> Residuals
#> 
#>      Min       1Q   Median       3Q      Max 
#> -1195.87  -422.48    50.51   354.88  7443.82 
#> 
#> 
#>           Coef       S.E.       t      Pr(>|t|)
#> Intercept 13307.3134 33589.2033   0.40 0.6921  
#> w           -70.6463   505.9989  -0.14 0.8890  
#> w'        -3747.8989  1579.5233  -2.37 0.0180  
#> w''       13158.8003  5860.1185   2.25 0.0252  
#> h         -4231.0573   205.4127 -20.60 <0.0001 
#> h'         2049.5490   487.8048   4.20 <0.0001 
#> h''       -9695.6103  2111.5310  -4.59 <0.0001 
#> w * h       129.3205     3.0807  41.98 <0.0001 
#> w * h'      -28.7097     7.0092  -4.10 <0.0001 
#> w * h''     136.2667    30.4018   4.48 <0.0001 
#> w' * h       41.0340     9.2103   4.46 <0.0001 
#> w'' * h    -101.5784    34.2189  -2.97 0.0031  
#> a           -23.1594    28.3581  -0.82 0.4145

# model parameterizations
parameterizations(pmodel1)
#> $w
#> function (t, c) 
#> (c * t^2)
#> <environment: 0x61bf4b7cef48>
#> 
#> $h
#> function (t) 
#> (t)
#> <environment: 0x61bf4b7cef48>

# new model variables
parameter_variables(pmodel1)
#> [1] "t" "c"
```

The methods are as follows:

- **`parameterize()`** — *Parameterize one or more predictor variables
  of a fitted model.* Takes a model and a named list of functions,
  expressions, or character strings defining the parameterizations.
  These are converted to functions and associated with predictor
  variables of the original model. When applied to an existing
  `parameterized_model`, the new parameterizations are composed with the
  existing ones and recorded in the parameterization history.

- **`parameterized_model()`** — *Construct a `parameterized_model`
  directly.* This is the low-level constructor underlying
  `parameterize()`. It takes an original model and a named list of
  parameterization functions, retains those corresponding to model
  predictors, determines the variables required by those functions, and
  stores the resulting parameterizations and their history alongside the
  original model.

- **`predict()`** — *Generate predictions from a parameterized model.*
  The `predict.parameterized_model()` method first transforms `newdata`
  according to the model’s parameterizations and then delegates
  prediction to the original model’s `predict()` method. Consequently,
  `newdata` is specified in terms of the parameter variables rather than
  the parameterized predictors.

- **`original_model()`** — *Extract the underlying fitted model.*
  Returns the original model stored inside a `parameterized_model`. If
  the supplied object is not a `parameterized_model`, it is returned
  unchanged.

- **`parameterizations()`** — *Extract the parameterization functions
  associated with a model.* By default, returns the model’s current
  parameterizations. The `level` argument can instead be used to
  retrieve and compose a specified number of the most recent
  parameterization steps from the model’s parameterization history.

- **`parameter_variables()`** — *Return the variables on which the
  parameterized model depends.* Extracts the unique variables appearing
  as arguments of the model’s current parameterization functions. These
  are the variables that must be supplied when constructing prediction
  data for the parameterized predictors.

- **`parameterize_data()`** — *Apply a model’s parameterizations to a
  data frame.* Evaluates each parameterization using the corresponding
  parameter variables in `newdata` and adds the resulting parameterized
  predictors to the data. The optional `level` argument controls how
  many levels of the parameterization history are applied.

- **`print()`** — *Display a summary of a `parameterized_model`.* The
  `print.parameterized_model()` method prints the underlying model, its
  parameter variables, and the current parameterization functions.
  Setting `show_history = TRUE` additionally displays the sequence of
  parameterizations previously applied to the model.

- **`is.parameterized_model()`** — *Test whether an object is a
  `parameterized_model`.* Returns a logical value indicating whether the
  object inherits from the `"parameterized_model"` class.

## Predictions for parameterized models

`term_predictions()` works for parameterized model objects as for fitted
models, with some additional features:

- Regular grids (for heatmap visualization); the ranges of new variables
  are found by a search, or given with `values`
- No support for offset terms
- Computation of critical values for ‘simultaneous’ confidence bands by
  simulation from the coefficients’ distribution, or by non-parametric
  bootstrap.

Additionally, implements simultaneous testing for parameterized models

``` r
x <- rnorm(500)
z <- rnorm(500)
y <- x + z + rnorm(500)

model <- lm(y ~ x + z)
pmodel <- parameterize(model, parameterizations = list(x = \(t, s) t^3 + s, z = \(t) t^2))

ggout <- term_predictions(pmodel, terms = c("t", "s"), values = list(s = c(-1, 0, 1)))

ggout
#> Predictions over t, s ( 90 grid points, response scale )
#> Critical values for 95% intervals (analytic): pointwise 1.965, simultaneous 2.574
#> 
#>        t  s estimate conf.low conf.high band.low band.high
#>  -1.2839 -1  -1.4941   -1.804   -1.1847   -1.899   -1.0888
#>  -1.1954 -1  -1.3093   -1.581   -1.0372   -1.666   -0.9529
#>  -1.1068 -1  -1.1646   -1.404   -0.9249   -1.478   -0.8506
#>  -1.0183 -1  -1.0559   -1.268   -0.8438   -1.334   -0.7781
#>  -0.9297 -1  -0.9790   -1.168   -0.7902   -1.226   -0.7316
#>  -0.8412 -1  -0.9299   -1.100   -0.7601   -1.152   -0.7075
#>  -0.7526 -1  -0.9044   -1.059   -0.7498   -1.107   -0.7019
#>  -0.6641 -1  -0.8982   -1.041   -0.7555   -1.085   -0.7112
#>  -0.5755 -1  -0.9074   -1.041   -0.7734   -1.083   -0.7319
#>  -0.4870 -1  -0.9276   -1.055   -0.7999   -1.095   -0.7603
#> ... with 80 more rows

ggout |> as.data.frame() |> head()
#>            t  s         x         z   estimate  std.error  conf.low  conf.high  band.low  band.high
#> 1 -1.2839185 -1 -3.116471 1.6484466 -1.4941338 0.15747983 -1.803542 -1.1847255 -1.899440 -1.0888276
#> 2 -1.1953724 -1 -2.708086 1.4289151 -1.3092623 0.13845402 -1.581290 -1.0372350 -1.665602 -0.9529229
#> 3 -1.1068263 -1 -2.355933 1.2250644 -1.1645541 0.12198214 -1.404218 -0.9248898 -1.478500 -0.8506084
#> 4 -1.0182802 -1 -2.055849 1.0368945 -1.0558597 0.10791957 -1.267895 -0.8438249 -1.333613 -0.7781069
#> 5 -0.9297341 -1 -1.803667 0.8644054 -0.9790299 0.09611929 -1.167880 -0.7901797 -1.226412 -0.7316475
#> 6 -0.8411880 -1 -1.595222 0.7075972 -0.9299154 0.08642611 -1.099721 -0.7601098 -1.152350 -0.7074803


### prediction testing
# requires function in original variables
test_simultaneous_predictions(ggout, f = \(x, z) x - z) |> as.data.frame() |> head()
#>            t  s         x         z   estimate  std.error  conf.low  conf.high  band.low  band.high reference difference inside
#> 1 -1.2839185 -1 -3.116471 1.6484466 -1.4941338 0.15747983 -1.803542 -1.1847255 -1.899440 -1.0888276 -4.764918   3.270784  FALSE
#> 2 -1.1953724 -1 -2.708086 1.4289151 -1.3092623 0.13845402 -1.581290 -1.0372350 -1.665602 -0.9529229 -4.137001   2.827738  FALSE
#> 3 -1.1068263 -1 -2.355933 1.2250644 -1.1645541 0.12198214 -1.404218 -0.9248898 -1.478500 -0.8506084 -3.580998   2.416444  FALSE
#> 4 -1.0182802 -1 -2.055849 1.0368945 -1.0558597 0.10791957 -1.267895 -0.8438249 -1.333613 -0.7781069 -3.092744   2.036884  FALSE
#> 5 -0.9297341 -1 -1.803667 0.8644054 -0.9790299 0.09611929 -1.167880 -0.7901797 -1.226412 -0.7316475 -2.668073   1.689043  FALSE
#> 6 -0.8411880 -1 -1.595222 0.7075972 -0.9299154 0.08642611 -1.099721 -0.7601098 -1.152350 -0.7074803 -2.302819   1.372904  FALSE
test_simultaneous_predictions(ggout, f = \(x, z) x  +  z) -> test

attr(test, "all_inside")
#> [1] TRUE
test |> as.data.frame() |> head()
#>            t  s         x         z   estimate  std.error  conf.low  conf.high  band.low  band.high  reference  difference inside
#> 1 -1.2839185 -1 -3.116471 1.6484466 -1.4941338 0.15747983 -1.803542 -1.1847255 -1.899440 -1.0888276 -1.4680244 -0.02610933   TRUE
#> 2 -1.1953724 -1 -2.708086 1.4289151 -1.3092623 0.13845402 -1.581290 -1.0372350 -1.665602 -0.9529229 -1.2791705 -0.03009182   TRUE
#> 3 -1.1068263 -1 -2.355933 1.2250644 -1.1645541 0.12198214 -1.404218 -0.9248898 -1.478500 -0.8506084 -1.1308690 -0.03368503   TRUE
#> 4 -1.0182802 -1 -2.055849 1.0368945 -1.0558597 0.10791957 -1.267895 -0.8438249 -1.333613 -0.7781069 -1.0189546 -0.03690510   TRUE
#> 5 -0.9297341 -1 -1.803667 0.8644054 -0.9790299 0.09611929 -1.167880 -0.7901797 -1.226412 -0.7316475 -0.9392617 -0.03976816   TRUE
#> 6 -0.8411880 -1 -1.595222 0.7075972 -0.9299154 0.08642611 -1.099721 -0.7601098 -1.152350 -0.7074803 -0.8876250 -0.04229033   TRUE
```

Remark that the obtained critical value and respective (hyperbolic)
simultaneous confidence interval are simultaneous over only the
generated grid. In particular, non-focal terms are marginalized to
constants.

### Simultaneous confidence bands

#### Definition and Motivation:

For a regression model of the form $`y \sim f_\beta(x) + \varepsilon`$,
the standard confidence band for predictions, $`f_{\hat\beta}(x_0)`$ is
the set:
``` math
  C_{\beta}(x_0) = [\hat f_\beta(x_0) - t_{1-\alpha/2, n-p} SE(\hat f_\beta)(x_0), \hat f_\beta(x_0) + t_{1-\alpha/2, n-p} SE(\hat f_\beta)(x_0)]
```

The set-valued function $`C_\beta: \Omega \to 2^\mathbb{R}`$ is the
pointwise confidence interval. In particular, one only draws pointwise
inferences about $`f_\beta`$, effectively, by means of hypothesis tests
of the for $`H_0: f_\beta(x) = f_0`$, which we reject iff
$`f_\beta(x) \not\in C_\beta(x)`$. One may naively attempt to extend
this pointwise test to a global test by means of multiple joint
hypothesis tests, $`H_0`$ & $`H_1`$ & $`\cdots`$ & $`H_n`$–corrected
with multiple hypothesis corrections. In the case where we correct
uniformly accros the domain of interest, $`\Omega`$, we obtain the
standard hyperbolic simultaneous confidence band: the set-valued
function $`C_\beta: \Omega \to 2^{\mathbb R}`$ given by

``` math
  C_\beta(x_0) = [\hat f_\beta(x_0) - c_\alpha SE(\hat f_\beta)(x_0), \hat f_\beta(x_0) + c_\alpha SE(\hat f_\beta)(x_0)]
```
where $`c_\alpha`$ is chosen such that is the 95-percent quantile of
$`\sup_{x\in \Omega}\frac{\lVert f_\beta(x) - \hat f_\beta(x)\rVert}{SE \hat f_\beta(x)}`$.
By construction, particulary, the supremum, we observe that this
confidence band “considers” the covariance structure of the regression
model.

Observe that there is not an elementary closed form for the
distribution, in general. Hence, we compute it via simulation or
bootstrap based methods. For nice models (ie. generlized linear models)
we can sample from the distribution of the regression coefficients, for
less nice models, we must refit several models, then estimate. In
particular, in termeffects, we implement two methods (in addition to the
pointwise confidence band). For generalized linear models the band is
computed on the scale of the linear predictor and then transformed, so
it respects the range of the outcome (e.g. probabilities stay within 0
and 1).

``` r

x <- rnorm(500)
z <- rnorm(500)
y <- rnorm(500) + x^2 - 0.5 * z

model <- lm(y ~ x + I(x^2) + I(x * z) + z)

out2 <- term_predictions(model, terms = c("x", "z"), values = list(z = c(-1, 0, 1)), ci_level = 0.95, ci_method = "bootstrap", nboot = 500)
out3 <- term_predictions(model, terms = c("x", "z"), values = list(z = c(-1, 0, 1)), ci_level = 0.95, ci_method = "analytic")


head(out3)
#> Predictions over x, z ( 6 grid points, response scale )
#> Critical values for 95% intervals (analytic): pointwise 1.965, simultaneous 2.894
#> 
#>       x  z estimate conf.low conf.high band.low band.high
#>  -2.968 -1    9.838    9.110    10.565    8.766    10.909
#>  -2.761 -1    8.579    7.934     9.224    7.629     9.529
#>  -2.553 -1    7.411    6.842     7.980    6.573     8.249
#>  -2.345 -1    6.334    5.835     6.833    5.599     7.069
#>  -2.137 -1    5.347    4.912     5.782    4.707     5.988
#>  -1.930 -1    4.451    4.074     4.828    3.896     5.007

term_curves(out3, band = "pointwise")
```

![](introduction_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
term_curves(out3, band = "simultaneous")
```

![](introduction_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r


sim1 <- test_simultaneous_predictions(out3, f = \(x, z) x - 0.5 * z)
head(sim1)
#> Simultaneous test: does the reference function fit inside the 95% band?
#> max |t| = 42.69, p = <1e-04
#> The reference lies inside the band at 0 of 6 grid points.
#> 
#>       x  z estimate reference band.low band.high inside
#>  -2.968 -1    9.838    -2.468    8.766    10.909  FALSE
#>  -2.761 -1    8.579    -2.261    7.629     9.529  FALSE
#>  -2.553 -1    7.411    -2.053    6.573     8.249  FALSE
#>  -2.345 -1    6.334    -1.845    5.599     7.069  FALSE
#>  -2.137 -1    5.347    -1.637    4.707     5.988  FALSE
#>  -1.930 -1    4.451    -1.430    3.896     5.007  FALSE
plot(sim1)
```

![](introduction_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->
In addition, by specifying `ci_scope = c("global", "faceted", "curve")`,
one may impose further marginalization (with respect to the uncertainty
estimation). for the purpose of inference:

    out <- term_predictions(model, terms = c("x", "z", "a", "b"), values = list(a = c(-1, 0, 1), b = c(-1, 0, 1)), ci_scope = "faceted")

provides simultaneous confidence bands
$`C^\alpha_{a, b}: \Omega|_{a, b} \to 2^{\mathbb R}`$ for
$`(a,b) \in \{-1, 0, 1\}^2`$.

#### Utility

In the case of regression models with more than 4 variables, the
simultaneous confidence band based inference implemented in termeffects
allows far more restricted inferences. Indeed, the underlying
foundational object, the `term_predictions` object, is a marginalized
data grid, together with its respective prediction. The coefficient
based inference is determined from the observations, whereas (still
essentially coefficient based) simultaneous confidence testing requires
the additional specification of the domain (which can be lost).

## term_layer objects and plotting

Given some base `term_predictions` object (an adjusted prediction grid
over two continuous focal terms), a `term_layer` object holds
predictions along curves through it, where the underlying model is a
‘parameterization’ of the base object’s underlying model subject to some
restrictions: - the held-constant variables are all equal to those of
the base - only the two primary variables were parameterized - the total
number of new parameters is at most two In particular, each layer either
induces a 1-parametric curve or 2-parameteric curve (which induces a set
of 1-parametric curves) in the domain of the base model.

`term_layer` objects were introduced due to the additional difficulty of
interpreting heatmap visualizations. The main difficulties are two fold
(but both the result of colour scale usage): difficulty of
distinguishing color values, and diffiliculty contrasting values between
points. `term_layer` objects track data for plotting, allowing the
generatation of ‘anotated’ heatmaps–heatmaps with overlays of discrete
families of 1-parameteric curves, and plots of prediction data
variability along the curves.

Note: marginalizations of predictions are parameteric. For example, in
the case of standard marginalization $`x_i = c`$, this method
(`term_layer()`) requires parameterizations:
`list(x1 = "t1", ..., xn = "tn")`.

### Plotting

As described, `termeffects` provides *two* plot types for surfaces, and
various subfeatures, namely:

- Heatmaps (`term_heatmap()`)
  - standard error or band-width plotting (on seperate heatmap,
    `fill = "std.error"` or `"band_width"`)
  - kernel smoothing density estimation (`show_support = TRUE`)
- Heatmaps with overlays (`layers =`)
- Difference plotting for predicted value vs estimated function
  (`plot(test, "difference")`)
- directional plotting (`show_direction = TRUE`)

Ex:

``` r
x <- rnorm(500)
z <- rnorm(500)

y <- x - z + rnorm(500, sd = 2)


model <- lm(y ~ x + z)
pmodel <- parameterize(model, list(x = "s", z = "s - t"))

pout <- term_predictions(pmodel, terms = c("t", "s"))


term_heatmap(pout, show_support = TRUE)
```

![](introduction_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r

ppout1 <- term_layer(pout,
                     parameterization = list(t = "w^2",
                                             s = "w^3"),
                     terms = "w",
                     values = list(w = seq(-1, 1, by = 0.05)))

ppout2 <- term_layer(pout,
                     parameterization = list(t = "r + w^2",
                                             s = "r + w^3"),
                     terms = c("w", "r"),
                     values = list(w = seq(-1, 1, by = 0.05), r = c(-1, 0, 1)))


ppout1 |> head()
#> Predictions over w ( 6 grid points, response scale )
#> Critical values for 95% intervals (analytic): pointwise 1.965, simultaneous 2.443
#> 
#>      w estimate conf.low conf.high band.low band.high
#>  -1.00   0.8454   0.3862    1.3045   0.2745    1.4162
#>  -0.95   0.7662   0.3548    1.1777   0.2547    1.2778
#>  -0.90   0.6904   0.3214    1.0593   0.2316    1.1491
#>  -0.85   0.6179   0.2863    0.9494   0.2057    1.0301
#>  -0.80   0.5488   0.2497    0.8479   0.1770    0.9206
#>  -0.75   0.4833   0.2119    0.7547   0.1458    0.8207

ppout2 |> head()
#> Predictions over w, r ( 6 grid points, response scale )
#> Critical values for 95% intervals (analytic): pointwise 1.965, simultaneous 2.587
#> 
#>      w  r estimate conf.low conf.high band.low band.high
#>  -1.00 -1  -0.2420  -0.8107   0.32661  -0.9907  0.506656
#>  -0.95 -1  -0.3212  -0.8421   0.19970  -1.0070  0.364628
#>  -0.90 -1  -0.3970  -0.8752   0.08114  -1.0267  0.232557
#>  -0.85 -1  -0.4696  -0.9100  -0.02915  -1.0494  0.110300
#>  -0.80 -1  -0.5386  -0.9459  -0.13128  -1.0749 -0.002312
#>  -0.75 -1  -0.6041  -0.9828  -0.22540  -1.1027 -0.105487

term_heatmap(pout, layers = ppout1)
```

![](introduction_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
term_curves(ppout1, show_direction = TRUE)
```

![](introduction_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r

term_heatmap(pout, layers = ppout2, show_direction = FALSE)
```

![](introduction_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

# Parameterized models

> **About this vignette.** Written by Lukas Bolsinger for p-ggeffects
> (original file `introduction_parameterized_model.Rmd`) and ported to
> termeffects by Brigitte Mueller. His text is kept; the code now uses
> the termeffects functions (`ggpredict()` became `term_predictions()`,
> `plot(..., "heatmap")` became `term_heatmap()`, layers became
> `term_layer()`), and held covariates are at their median or most
> common value. Also corrected: a typo (`pmodel`) that stopped two
> chunks from running, and simulated weights with mean 7000 instead of
> 70. Natural splines (`ns()`) are replaced by restricted cubic splines
> (`rcs()`), as used in our team. The model with `%ia%` is now fitted
> with `ols()` instead of `lm()`: in `lm()`, `%ia%` recomputes the
> spline knots from the new data at prediction time, so the predictions
> (and the earlier figures) were wrong. Every change is marked in
> `validation/output/vignette_track_changes.html` (struck through:
> removed; underlined: added), or run
> `git diff 2de18bc -- vignettes/parameterized_models.Rmd`.

# parameterized_model objects

## Overview

A `parameterized_model` is a data structure consiting of a tuple of a
fitted regression model, and a collection of maps that parameterize the
continuous predictive variables of the model. Moreover, there exists
some auxiliary tracked information, ie. new parameters, and a
parameterization history.

The primary objective of parameterized models is the visualization of
higher order data in terms of lower order parameteric families via
marginalization. As an example, consider the fitted regression model:
``` math
y = f(x_1, \ldots, x_n; \beta)\quad\text{ with parameterizations}\quad x_i = h_i(z_1, \ldots, z_m)
```
Viewing the regression model as a function, notice that the data of the
parameterized model determines the data of a model by pulling back by
the $`h_i`$. We do note however, that the observed data for which
$`f_\beta = f(-; \beta)`$ was fitted may not be contained in the image
of $`h = (h_i)_i`$, and so the parameterized model $`f_\beta \circ h`$
is almost never determined by the initial observations. For this reason,
this parameterization is primarily a tool for developing intuitions and
visualizations.

Indeed, in sufficiently nice scenarios (for instance if $`h`$ is
smooth), the rank of the parameterized model (viewed as the composition
$`f_\beta \circ h`$) is less than or equal to the rank of $`f_\beta`$.
Particularly, by fixing arguments in $`z`$ (to obtain planes), we obtain
more varied marginalizations of $`f_\beta`$.

## Simple Example

Suppose we have some regression model
$`y \sim f_\beta(\operatorname{weight}, \operatorname{height}, \ldots)`$.
With parameterizations
$`(t, c, \ldots) \mapsto (ct^2, t, \ldots) = (\operatorname{weight}(t), \operatorname{height}(t), \ldots)`$
(where the other function are the identity), we obtain co-level sets in
c and $`\ldots`$ of subsets with fixed BMI, and co-level sets in $`t`$
and $`\ldots`$ of subsets with fixed height, but variable weight and
BMI.

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
#> <environment: 0x61bf597825a0>
#> 
#> $h
#> function (t) 
#> (t)
#> <environment: 0x61bf597825a0>

# new model variables
parameter_variables(pmodel1)
#> [1] "t" "c"
```

Furthermore, for a data grid
`data = data.frame(t = ..., c = ..., a = ...)`, with
`predict(pmodel, newdata = data)`, we obtain predictions from the
mechanism:

``` r
# predict.parameterized_model <- function(object, newdata, ...) {
#  ...
#
#  newdata <- parameterize_data(object, newdata)
#
#  predict(object$model, newdata = newdata, ...)
#}

# t is height in cm, so c is BMI in kg/cm^2 (BMI 18-22 kg/m^2 = 0.0018-0.0022)
df <- as.data.frame(
    expand.grid(t = 150:200, c = c(18, 20, 22) / 1e4, a = 0)
  )

head(parameterize_data(pmodel1, newdata = df))
#>     t      c a       w   h
#> 1 150 0.0018 0 40.5000 150
#> 2 151 0.0018 0 41.0418 151
#> 3 152 0.0018 0 41.5872 152
#> 4 153 0.0018 0 42.1362 153
#> 5 154 0.0018 0 42.6888 154
#> 6 155 0.0018 0 43.2450 155

head(predict(pmodel1, newdata = df))
#>        1        2        3        4        5        6 
#> 161409.8 172957.9 184716.6 196687.4 208871.6 221270.9
```

where `parameterize_data` evaluates the parameterization on the new
data.

## term_predictions with parameterized_models.

`term_predictions()` and `term_grid()` are implemented for
parameterized_model objects. As with fitted regression models, the user
specifies their primary terms of interest
`term_predictions(model, terms, ...)` together with auxiliary options.
These terms are variables of the parameterizations, or those
untransformed variables:

    term_predictions(pmodel, terms = c("t", "c"), values = list(c = c(1, 2, 3)), ...)

For standard regression models, for unconstrained terms–such as
“t”–bounds, and regularly spaced sample points are determined by the
data (the extrema of the observations). As discussed above, such a
method does not extend to parameterized_models. Rather, for all new
parameter variables, the largest box (centred at the origin in the
parameter space) whose image stays within the extrema of the
observations is computed. Where this is not possible (consider a
parameterized point or circle), a warning asks for the range to be given
with `values`.

## term_layer objects

Given some base `term_predictions` object (an adjusted prediction grid
over two continuous focal terms, as drawn by `term_heatmap()`), a
`term_layer` object holds predictions along curves through it, where the
underlying model is a ‘parameterization’ of the base object’s underlying
model subject to some restrictions: - the held-constant variables are
all equal to those of the base - only the two primary variables were
parameterized - the total number of new parameters is at most two In
particular, each layer either induces a 1-parametric curve or
2-parameteric curve (which induces a set of 1-parametric curves) in the
domain of the base model.

Ex:

``` r
x <- rnorm(500)
z <- rnorm(500)

y <- x - z + rnorm(500, sd = 2)


model <- lm(y ~ x + z)
pmodel <- parameterize(model, list(x = "t^3 + s", z = "s - t"))

pout <- term_predictions(pmodel, terms = c("t", "s"))

ppout1 <- term_layer(pout,
                     parameterization = list(t = "sin(w)",
                                             s = "cos(w)"),
                     terms = "w",
                     values = list(w = seq(0, 2 * pi, length.out = 60)))

ppout2 <- term_layer(pout,
                     parameterization = list(t = "r * sin(w)",
                                             s = "r * cos(w)"),
                     terms = c("w", "r"),
                     values = list(w = seq(0, 2 * pi, length.out = 60), r = c(0.5, 1, 1.5)))


ppout1
#> Predictions over w ( 60 grid points, response scale )
#> Critical values for 95% intervals (analytic): pointwise 1.965, simultaneous 2.621
#> 
#>       w  estimate conf.low conf.high band.low band.high
#>  0.0000 -0.098288 -0.39898    0.2024 -0.49936    0.3028
#>  0.1065  0.006792 -0.28239    0.2960 -0.37892    0.3925
#>  0.2130  0.117524 -0.15974    0.3948 -0.25230    0.4873
#>  0.3195  0.238946 -0.02715    0.5050 -0.11597    0.5939
#>  0.4260  0.374631  0.11783    0.6314  0.03211    0.7172
#>  0.5325  0.526143  0.27576    0.7765  0.19218    0.8601
#>  0.6390  0.692706  0.44517    0.9402  0.36254    1.0229
#>  0.7455  0.871120  0.62259    1.1197  0.53963    1.2026
#>  0.8520  1.055925  0.80287    1.3090  0.71840    1.3935
#>  0.9585  1.239817  0.97952    1.5001  0.89264    1.5870
#> ... with 50 more rows

ppout2
#> Predictions over w, r ( 180 grid points, response scale )
#> Critical values for 95% intervals (analytic): pointwise 1.965, simultaneous 2.774
#> 
#>       w   r   estimate conf.low conf.high   band.low band.high
#>  0.0000 0.5 -0.1053510 -0.31900    0.1083 -0.4070406    0.1963
#>  0.1065 0.5 -0.0532577 -0.26276    0.1562 -0.3490820    0.2426
#>  0.2130 0.5 -0.0009594 -0.20605    0.2041 -0.2905564    0.2886
#>  0.3195 0.5  0.0517380 -0.14892    0.2524 -0.2315975    0.3351
#>  0.4260 0.5  0.1048559 -0.09157    0.3013 -0.1725088    0.3822
#>  0.5325 0.5  0.1581803 -0.03444    0.3508 -0.1138033    0.4302
#>  0.6390 0.5  0.2112246  0.02182    0.4006 -0.0562186    0.4787
#>  0.7455 0.5  0.2632234  0.07631    0.4501 -0.0007066    0.5272
#>  0.8520 0.5  0.3131582  0.12793    0.4984  0.0516033    0.5747
#>  0.9585 0.5  0.3598115  0.17543    0.5442  0.0994605    0.6202
#> ... with 170 more rows

term_heatmap(pout, layers = ppout1, show_support = TRUE, support_prob = 0.975)
```

![](parameterized_models_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r

term_curves(ppout2)
```

![](parameterized_models_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

# Plotting

termeffects provides *two* plot types for surfaces: heatmaps, and
heatmaps with overlays of 1-parametric curves. Additionally, as a
substitute for plotting observed data, it estimates the density of the
initial observations using kernel smoothing, and fades the regions where
it is low.

Heatmaps are drawn on regular grids, which `term_predictions()`
generates.

``` r

w <- rnorm(500, mean = 70, sd = 4)
h <- rnorm(500, mean = 170, sd = 10)

a <- rnorm(500)

y <- a + w/h^2 + h + w + w^2 * h + rnorm(500)


# model example
model <- ols(y ~ rcs(w, 4) + rcs(h, 4) + rcs(w, 4) %ia% rcs(h, 4) + a)

# parameterize model
pmodel1 <- parameterize(model, parameterizations = list(w = \(t,c) c*t^2, h = \(t) t))

pout <- term_predictions(pmodel1, terms = c("c", "t"),
                         values = list(c = seq(0.0020, 0.0030, length.out = 30),
                                       t = seq(150, 190, length.out = 30)),
                         at = list(a = 1))

term_heatmap(pout, show_support = TRUE)
```

![](parameterized_models_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

# Simultaneous confidence bands

Every `term_predictions()` result also contains a simultaneous
confidence band (`band.low`, `band.high`), which covers the predictions
at all grid points at once. For parameterized models it is simultaneous
over the grid of the new variables. The vignette “termeffects:
introduction” explains how the band is computed and how to test a
reference function against it.

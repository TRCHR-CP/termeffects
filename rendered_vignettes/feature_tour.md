# A tour of every feature

This vignette uses every function in termeffects once, on the
`heartfailure` data. It is meant for trying the package out: open
`vignettes/feature_tour.Rmd` in RStudio and run it chunk by chunk. Each
section ends with **Check:**, saying what you should see. Features that
came from Lukas Bolsinger’s p-ggeffects are marked *(from p-ggeffects)*.

The other vignettes explain the ideas; this one only shows that
everything works.

## 0. Before you start: the automatic checks

These run from the package folder in RStudio:

``` r
devtools::load_all()   # load the current code without installing
devtools::test()       # about 160 automatic tests; all should pass
devtools::check()      # full R CMD check: 0 errors expected
devtools::build_vignettes()   # renders all vignettes, including this one
```

Then install the package (`devtools::install()`) and run this vignette.

## 1. One model to start with

A linear model of ejection fraction (LVEF) with a spline interaction
between creatinine and sodium, adjusted for age, sex and diabetes. It is
quick to fit, so it is used for most of the tour; section 7 repeats the
key steps for other model types.

``` r
library(termeffects)
library(rms)
library(splines)

hf <- heartfailure
fit <- lm(ejection_fraction ~ ns(creatinine, 3) * ns(sodium, 3) + ns(age, 3) + sex + diabetes,
          data = hf)
```

## 2. The prediction grid: `term_grid()`

``` r
grid <- term_grid(fit, c("creatinine", "sodium"), n = 5)
head(grid)
#>   creatinine sodium age  sex diabetes
#> 1      0.500 113.00  60 Male       No
#> 2      2.725 113.00  60 Male       No
#> 3      4.950 113.00  60 Male       No
#> 4      7.175 113.00  60 Male       No
#> 5      9.400 113.00  60 Male       No
#> 6      0.500 121.75  60 Male       No
attr(grid, "held")
#> $age
#> [1] 60
#> 
#> $sex
#> [1] Male
#> Levels: Female Male
#> 
#> $diabetes
#> [1] No
#> Levels: No Yes
attr(grid, "held_rule")
#>           age           sex      diabetes 
#>      "median" "most common" "most common"
```

**Check:** creatinine and sodium vary; age is held at its median (60),
sex and diabetes at their most common value (Male, No).

## 3. Predictions: `term_predictions()`

``` r
set.seed(1)
cr <- seq(0.5, 4, length.out = 40)     # where almost all patients are
na <- seq(125, 148, length.out = 40)
surface <- term_predictions(fit, c("creatinine", "sodium"), values = list(creatinine = cr, sodium = na))
surface
#> Predictions over creatinine, sodium ( 1600 grid points, response scale )
#> Held constant: age = 60 (median), sex = Male (most common), diabetes =
#>   No (most common)
#> Critical values for 95% intervals (analytic): pointwise 1.969, simultaneous 3.310
#> 
#>  creatinine sodium estimate conf.low conf.high band.low band.high
#>      0.5000    125  -39.340 -114.289     35.61 -165.353     86.67
#>      0.5897    125  -20.883  -79.076     37.31 -118.725     76.96
#>      0.6795    125   -3.267  -45.762     39.23  -74.715     68.18
#>      0.7692    125   12.668  -16.240     41.58  -35.935     61.27
#>      0.8590    125   26.080    7.197     44.96   -5.668     57.83
#>      0.9487    125   36.128   22.046     50.21   12.452     59.80
#>      1.0385    125   42.003   28.598     55.41   19.465     64.54
#>      1.1282    125   43.881   30.770     56.99   21.838     65.92
#>      1.2179    125   43.131   31.197     55.07   23.066     63.20
#>      1.3077    125   41.195   30.894     51.50   23.875     58.51
#> ... with 1590 more rows
head(as.data.frame(surface))
#>   creatinine sodium age  sex diabetes  estimate std.error    conf.low conf.high
#> 1  0.5000000    125  60 Male       No -39.33974 38.073489 -114.288696  35.60922
#> 2  0.5897436    125  60 Male       No -20.88285 29.561882  -79.076418  37.31072
#> 3  0.6794872    125  60 Male       No  -3.26678 21.587425  -45.762360  39.22880
#> 4  0.7692308    125  60 Male       No  12.66764 14.684822  -16.239931  41.57521
#> 5  0.8589744    125  60 Male       No  26.07959  9.592180    7.197055  44.96212
#> 6  0.9487179    125  60 Male       No  36.12824  7.153666   22.046003  50.21047
#>      band.low band.high
#> 1 -165.352771  86.67329
#> 2 -118.724745  76.95905
#> 3  -74.715365  68.18180
#> 4  -35.935180  61.27046
#> 5   -5.667954  57.82713
#> 6   12.451523  59.80495
```

**Check:** the header names the held values and gives two critical
values, 1.97 for pointwise intervals and a larger one for the
simultaneous band. The columns are `estimate`, `std.error`,
`conf.low`/`conf.high` (pointwise) and `band.low`/`band.high`
(simultaneous).

The covariance matrix of the predictions *(from p-ggeffects)*:

``` r
V <- vcov(surface)
dim(V)
#> [1] 1600 1600
all.equal(sqrt(diag(V)), surface$std.error)
#> [1] TRUE
```

**Check:** a 1600 × 1600 matrix whose diagonal gives the standard errors
(`TRUE`).

## 4. Heatmaps: `term_heatmap()`

``` r
term_heatmap(surface)
```

![](feature_tour_files/figure-gfm/heatmap-basic-1.png)<!-- -->

``` r
term_heatmap(surface, contours = FALSE)
```

![](feature_tour_files/figure-gfm/heatmap-basic-2.png)<!-- -->

Data points, the region with data (*from p-ggeffects*), and marks for
cuts:

``` r
term_heatmap(surface, show_data = TRUE, show_support = TRUE)
```

![](feature_tour_files/figure-gfm/heatmap-support-1.png)<!-- -->

``` r
term_heatmap(surface, show_support = TRUE, support_prob = 0.9)
```

![](feature_tour_files/figure-gfm/heatmap-support-2.png)<!-- -->

``` r
term_heatmap(surface, cuts = list(sodium = c(130, 137, 142)))
```

![](feature_tour_files/figure-gfm/heatmap-support-3.png)<!-- -->

**Check:** the faded area grows when `support_prob` drops from 0.975 to
0.9. The cuts are dashed horizontal lines.

Where is the surface estimated precisely? *(from p-ggeffects)*

``` r
term_heatmap(surface, fill = "band_width")
```

![](feature_tour_files/figure-gfm/heatmap-uncertainty-1.png)<!-- -->

``` r
term_heatmap(surface, fill = "std.error")
```

![](feature_tour_files/figure-gfm/heatmap-uncertainty-2.png)<!-- -->

**Check:** both are largest where there are few patients (high
creatinine, and the edges of the sodium range).

`plot()` picks the right plot for the object:

``` r
plot(surface)
```

![](feature_tour_files/figure-gfm/plot-auto-1.png)<!-- -->

## 5. Cuts: `term_curves()`

``` r
cuts <- term_predictions(fit, c("creatinine", "sodium"),
                         values = list(creatinine = cr, sodium = c(130, 137, 142)))
term_curves(cuts)                        # simultaneous band (default)
```

![](feature_tour_files/figure-gfm/curves-1.png)<!-- -->

``` r
term_curves(cuts, band = "pointwise")
```

![](feature_tour_files/figure-gfm/curves-2.png)<!-- -->

``` r
term_curves(cuts, band = "both")
```

![](feature_tour_files/figure-gfm/curves-3.png)<!-- -->

``` r
term_curves(cuts, band = "none")
```

![](feature_tour_files/figure-gfm/curves-4.png)<!-- -->

``` r
plot(cuts)
```

![](feature_tour_files/figure-gfm/curves-5.png)<!-- -->

**Check:** with `band = "both"` the simultaneous band is always the
wider one.

### Other held values: `at`

``` r
women <- term_predictions(fit, c("creatinine", "sodium"),
                          values = list(creatinine = cr, sodium = c(130, 137, 142)),
                          at = list(sex = "Female", diabetes = "Yes", age = 75))
women
#> Predictions over creatinine, sodium ( 120 grid points, response scale )
#> Held constant: age = 75 (set by user), sex = Female (set by user),
#>   diabetes = Yes (set by user)
#> Critical values for 95% intervals (analytic): pointwise 1.969, simultaneous 2.934
#> 
#>  creatinine sodium estimate conf.low conf.high band.low band.high
#>      0.5000    130    28.57    7.444     49.70   -2.923     60.06
#>      0.5897    130    31.76   15.358     48.17    7.308     56.22
#>      0.6795    130    34.79   22.716     46.87   16.790     52.80
#>      0.7692    130    37.49   28.970     46.02   24.788     50.20
#>      0.8590    130    39.70   33.455     45.95   30.390     49.01
#>      0.9487    130    41.25   35.753     46.75   33.054     49.45
#>      1.0385    130    41.99   36.432     47.55   33.705     50.27
#>      1.1282    130    41.95   36.416     47.48   33.702     50.19
#>      1.2179    130    41.40   36.136     46.66   33.555     49.24
#>      1.3077    130    40.63   35.705     45.55   33.290     47.96
#> ... with 110 more rows
term_curves(women)
```

![](feature_tour_files/figure-gfm/at-1.png)<!-- -->

**Check:** the header now says Female, Yes and 75, marked as set by you.

### Panels and the scope of the band *(from p-ggeffects)*

A third focal variable makes panels. `ci_scope` says over what the band
is simultaneous: everything (`"global"`), each panel (`"faceted"`) or
each curve (`"curve"`).

``` r
by_sex <- lapply(c(global = "global", faceted = "faceted", curve = "curve"), function(scope)
  term_predictions(fit, c("creatinine", "sodium", "sex"),
                   values = list(creatinine = cr, sodium = c(130, 137, 142)),
                   ci_scope = scope))
lapply(by_sex, attr, "critical_value")
#> $global
#> [1] 3.124853
#> 
#> $faceted
#>   Female     Male 
#> 2.981446 3.010862 
#> 
#> $curve
#> 130.Female   130.Male 137.Female   137.Male 142.Female   142.Male 
#>   2.642687   2.636335   2.563638   2.592961   2.614318   2.621560
term_curves(by_sex$global)
```

![](feature_tour_files/figure-gfm/scope-1.png)<!-- -->

**Check:** one critical value for `"global"`, one per sex for
`"faceted"`, one per curve for `"curve"`, and they get smaller in that
order.

### The link scale

``` r
glm_fit <- glm(death_event ~ ns(creatinine, 3) * ns(sodium, 3) + age + sex,
               data = hf, family = binomial)
risk <- term_predictions(glm_fit, c("creatinine", "sodium"),
                         values = list(creatinine = cr, sodium = c(130, 137, 142)))
term_curves(risk)
```

![](feature_tour_files/figure-gfm/link-1.png)<!-- -->

``` r
term_curves(term_predictions(glm_fit, c("creatinine", "sodium"), scale = "link",
                             values = list(creatinine = cr, sodium = c(130, 137, 142))))
```

![](feature_tour_files/figure-gfm/link-2.png)<!-- -->

**Check:** the first plot shows probabilities whose bands stay between 0
and 1; the second shows log-odds.

## 6. Interval methods

The bootstrap *(from p-ggeffects)*: the model is refitted on resampled
data.

``` r
set.seed(2)
cuts_boot <- term_predictions(fit, c("creatinine", "sodium"),
                              values = list(creatinine = cr, sodium = c(130, 137, 142)),
                              ci_method = "bootstrap", nboot = 200)
attr(cuts_boot, "critical_value")
#> [1] 3.547659
term_curves(cuts_boot)
```

![](feature_tour_files/figure-gfm/bootstrap-1.png)<!-- -->

**Check:** in the middle, where most patients are, the bands look like
the analytic ones in section 5. Above creatinine 3 (16 patients) and
below 0.7 (5 patients), the bootstrap bands are wider. The analytic band
assumes the model is right and the spread of LVEF is the same
everywhere; the bootstrap refits the model on resampled patients, and
where only a few patients determine the curve, the refits differ a lot
depending on which of them are drawn. Where the two disagree, trust the
curves less.

For `orm` models there are three methods (see the how-to vignette for
the simulation comparing them):

``` r
hf$platelets <- hf$platelets / 1000
ofit <- orm(ejection_fraction ~ rcs(creatinine, 4) * rcs(sodium, 4) + rcs(age, 4) +
              sex + diabetes + anaemia + hblood_pressure + smoking +
              rcs(platelets, 3) + rcs(log(cpk_enzyme), 3),
            data = hf, x = TRUE, y = TRUE)
orm_cuts <- lapply(c(analytic = "analytic", delta = "delta", bootstrap = "bootstrap"), function(m) {
  set.seed(3)
  term_predictions(ofit, c("creatinine", "sodium"),
                   values = list(creatinine = seq(0.5, 4, length.out = 25), sodium = c(130, 137, 142)),
                   ci_method = m, nboot = 200)
})
sapply(orm_cuts, function(p) mean(p$band.high - p$band.low))   # average band width
#>  analytic     delta bootstrap 
#>  18.67659  19.09745  31.55171
term_curves(orm_cuts$analytic)
```

![](feature_tour_files/figure-gfm/orm-1.png)<!-- -->

**Check:** the bootstrap bands are clearly wider than the other two.
(Here the delta bands are about as wide as the default ones; in the
simulations they were too narrow, which is why they are not the
default.)

## 7. Other model types

### Survival: a Cox model

A Cox model (`survival::coxph`) describes how the risk of death changes
with the predictors. termeffects shows the **hazard ratio** against an
average patient: 2 means twice the hazard of death at any time, 0.5 half
of it. Ratios are plotted on a log axis, so that 2 and 0.5 are equally
far from 1 (the dotted line, no difference).

``` r
library(survival)
hf$died <- as.numeric(hf$death_event == "Yes")
cox <- coxph(Surv(time, died) ~ ns(creatinine, 3) * ns(sodium, 3) + age + sex, data = hf)
term_curves(term_predictions(cox, c("creatinine", "sodium"),
                             values = list(creatinine = cr, sodium = c(130, 137, 142))))
```

![](feature_tour_files/figure-gfm/other-models-1.png)<!-- -->

**Check:** the y-axis is logarithmic (1, 10, 100, …). The bands are very
wide below creatinine 0.7, where there are only 5 patients.

### Smooth surfaces: a generalized additive model

A generalized additive model (`mgcv::gam`) does not need knots or
degrees of freedom chosen in advance. Each smooth term is a spline with
many basis functions, and a penalty for wiggliness is estimated from the
data, so the curve is only as flexible as the data support. `s(age)` is
a smooth curve in age. `te(creatinine, sodium)` is a smooth *surface*
over both variables together (a tensor-product spline), the GAM way of
writing a non-linear interaction. GAMs suit exploring a surface; spline
models with a fixed number of knots, as above, are easier to report.

``` r
gam_fit <- mgcv::gam(ejection_fraction ~ te(creatinine, sodium) + s(age) + sex, data = hf)
term_heatmap(term_predictions(gam_fit, c("creatinine", "sodium"), values = list(creatinine = cr, sodium = na)))
```

![](feature_tour_files/figure-gfm/gam-1.png)<!-- -->

**Check:** a smoother surface than the spline model in section 4. Mixed
models (`lme4`) work the same way; the `heartfailure` data have no
clusters, so they are covered by the automatic tests only.

## 8. A figure with its caption: `term_figure()`

``` r
ofit_surface <- term_predictions(ofit, c("creatinine", "sodium"),
                                 values = list(creatinine = cr, sodium = na))
fig <- term_figure(ofit_surface, cuts = list(sodium = c(130, 137, 142)), show_data = TRUE)
fig
```

![](feature_tour_files/figure-gfm/figure-1.png)<!-- -->![](feature_tour_files/figure-gfm/figure-2.png)<!-- -->

    #> Mean ejection_fraction over creatinine and sodium, from a proportional
    #> odds model (rms::orm), n = 299: ejection_fraction ~ rcs(creatinine, 4)
    #> * rcs(sodium, 4) + rcs(age, 4) + sex + diabetes + anaemia +
    #> hblood_pressure + smoking + rcs(platelets, 3) + rcs(log(cpk_enzyme),
    #> 3). In the model, creatinine enters as rcs(creatinine, 4) and sodium as
    #> rcs(sodium, 4). They interact in the model (rcs(creatinine,
    #> 4):rcs(sodium, 4)). All other covariates are held constant: age = 60
    #> (median), sex = Male (most common), diabetes = No (most common),
    #> anaemia = No (most common), hblood_pressure = No (most common), smoking
    #> = No (most common), platelets = 262 (median), cpk_enzyme = 250
    #> (median). Curves: cuts through the surface at sodium = 130, 137, 142.
    #> Shaded: 95% simultaneous confidence band, covering all curves at once
    #> (critical value 2.96).
    fig$heatmap

![](feature_tour_files/figure-gfm/figure-3.png)<!-- -->

``` r
fig$curves
```

![](feature_tour_files/figure-gfm/figure-4.png)<!-- -->

``` r
fig$uncertainty
```

![](feature_tour_files/figure-gfm/figure-5.png)<!-- -->

``` r
fig$covariates
#>         covariate value   chosen_as
#> 1             age    60      median
#> 2             sex  Male most common
#> 3        diabetes    No most common
#> 4         anaemia    No most common
#> 5 hblood_pressure    No most common
#> 6         smoking    No most common
#> 7       platelets   262      median
#> 8      cpk_enzyme   250      median
fig$model
#> [1] "proportional odds model (rms::orm), n = 299: ejection_fraction ~ rcs(creatinine, 4) * rcs(sodium, 4) + rcs(age, 4) + sex + diabetes + anaemia + hblood_pressure + smoking + rcs(platelets, 3) + rcs(log(cpk_enzyme), 3)"
```

``` r
cat(fig$caption)
```

Mean ejection_fraction over creatinine and sodium, from a proportional
odds model (rms::orm), n = 299: ejection_fraction ~ rcs(creatinine, 4)
\* rcs(sodium, 4) + rcs(age, 4) + sex + diabetes + anaemia +
hblood_pressure + smoking + rcs(platelets, 3) + rcs(log(cpk_enzyme), 3).
In the model, creatinine enters as rcs(creatinine, 4) and sodium as
rcs(sodium, 4). They interact in the model (rcs(creatinine,
4):rcs(sodium, 4)). All other covariates are held constant: age = 60
(median), sex = Male (most common), diabetes = No (most common), anaemia
= No (most common), hblood_pressure = No (most common), smoking = No
(most common), platelets = 262 (median), cpk_enzyme = 250 (median).
Curves: cuts through the surface at sodium = 130, 137, 142. Shaded: 95%
simultaneous confidence band, covering all curves at once (critical
value 2.96).

**Check:** with eight held covariates, the heatmap subtitle points to
the caption instead of listing them, and `fig$covariates` lists all
eight with how each value was chosen.

## 9. Parameterized models *(from p-ggeffects)*

A parameterized model gives the fitted model new inputs without
refitting it. Here creatinine is described by its logarithm:

``` r
pfit <- parameterize(fit, list(creatinine = function(log_cr) exp(log_cr)))
pfit
#> Parameterized model
#> 
#> Model:
#> 
#> Call:
#> lm(formula = ejection_fraction ~ ns(creatinine, 3) * ns(sodium, 
#>     3) + ns(age, 3) + sex + diabetes, data = hf)
#> 
#> Coefficients:
#>                       (Intercept)                 ns(creatinine, 3)1  
#>                         -270.9880                           -41.7106  
#>                ns(creatinine, 3)2                 ns(creatinine, 3)3  
#>                          928.3082                           908.3326  
#>                    ns(sodium, 3)1                     ns(sodium, 3)2  
#>                          168.9955                           594.0910  
#>                    ns(sodium, 3)3                        ns(age, 3)1  
#>                           88.1682                             1.5544  
#>                       ns(age, 3)2                        ns(age, 3)3  
#>                            5.2880                             7.5979  
#>                           sexMale                        diabetesYes  
#>                           -3.3220                            -0.1435  
#> ns(creatinine, 3)1:ns(sodium, 3)1  ns(creatinine, 3)2:ns(sodium, 3)1  
#>                           13.3967                          -439.6060  
#> ns(creatinine, 3)3:ns(sodium, 3)1  ns(creatinine, 3)1:ns(sodium, 3)2  
#>                         -368.2904                            60.1433  
#> ns(creatinine, 3)2:ns(sodium, 3)2  ns(creatinine, 3)3:ns(sodium, 3)2  
#>                        -1807.9032                         -1779.1195  
#> ns(creatinine, 3)1:ns(sodium, 3)3  ns(creatinine, 3)2:ns(sodium, 3)3  
#>                           35.1593                          -287.4718  
#> ns(creatinine, 3)3:ns(sodium, 3)3  
#>                         -330.1304  
#> 
#> 
#> Parameters:
#> [1] "log_cr"
#> 
#> Parameterizations:
#> $creatinine <- (exp(log_cr)); (log_cr)
is.parameterized_model(pfit)
#> [1] TRUE
parameter_variables(pfit)
#> [1] "log_cr"
parameterizations(pfit)
#> $creatinine
#> function (log_cr) 
#> (exp(log_cr))
#> <environment: 0x5a4db8ced018>
class(original_model(pfit))
#> [1] "lm"
```

Parameterizing again composes the two steps: here, log creatinine is
given as the number of doublings from 1 mg/dL.

``` r
pfit2 <- parameterize(pfit, list(log_cr = function(doublings) doublings * log(2)))
print(pfit2, show_history = TRUE)
#> Parameterized model
#> 
#> Model:
#> 
#> Call:
#> lm(formula = ejection_fraction ~ ns(creatinine, 3) * ns(sodium, 
#>     3) + ns(age, 3) + sex + diabetes, data = hf)
#> 
#> Coefficients:
#>                       (Intercept)                 ns(creatinine, 3)1  
#>                         -270.9880                           -41.7106  
#>                ns(creatinine, 3)2                 ns(creatinine, 3)3  
#>                          928.3082                           908.3326  
#>                    ns(sodium, 3)1                     ns(sodium, 3)2  
#>                          168.9955                           594.0910  
#>                    ns(sodium, 3)3                        ns(age, 3)1  
#>                           88.1682                             1.5544  
#>                       ns(age, 3)2                        ns(age, 3)3  
#>                            5.2880                             7.5979  
#>                           sexMale                        diabetesYes  
#>                           -3.3220                            -0.1435  
#> ns(creatinine, 3)1:ns(sodium, 3)1  ns(creatinine, 3)2:ns(sodium, 3)1  
#>                           13.3967                          -439.6060  
#> ns(creatinine, 3)3:ns(sodium, 3)1  ns(creatinine, 3)1:ns(sodium, 3)2  
#>                         -368.2904                            60.1433  
#> ns(creatinine, 3)2:ns(sodium, 3)2  ns(creatinine, 3)3:ns(sodium, 3)2  
#>                        -1807.9032                         -1779.1195  
#> ns(creatinine, 3)1:ns(sodium, 3)3  ns(creatinine, 3)2:ns(sodium, 3)3  
#>                           35.1593                          -287.4718  
#> ns(creatinine, 3)3:ns(sodium, 3)3  
#>                         -330.1304  
#> 
#> 
#> Parameters:
#> [1] "doublings"
#> 
#> Parameterizations:
#> $creatinine <- (exp(doublings * log(2))); (doublings)
#> 
#> 
#> Parameterization history:
#> 
#> Stage 1                        Stage 2                        
#> ----------------------------   -------------------------------
#> 
#> creatinine [exp(log_cr)] <─────┐
#>                                log_cr [doublings * log(2)] <── doublings
nd <- data.frame(doublings = c(-1, 0, 1, 2), sodium = 137, age = 60,
                 sex = "Male", diabetes = "No")
parameterize_data(pfit2, nd)
#>   doublings sodium age  sex diabetes creatinine
#> 1        -1    137  60 Male       No        0.5
#> 2         0    137  60 Male       No        1.0
#> 3         1    137  60 Male       No        2.0
#> 4         2    137  60 Male       No        4.0
predict(pfit2, nd)
#>        1        2        3        4 
#> 47.76279 37.04315 33.59173 35.69862
predict(fit, transform(nd, creatinine = 2^doublings))
#>        1        2        3        4 
#> 47.76279 37.04315 33.59173 35.69862
```

**Check:** `parameterize_data()` adds `creatinine` = 0.5, 1, 2, 4, and
the two `predict()` calls give the same numbers.

Predictions and plots work as for an ordinary model. Without `values`,
the range of a new variable is found by a search: a range centred at 0
within which creatinine stays inside its observed range.

``` r
p_surface <- term_predictions(pfit2, c("doublings", "sodium"))
range(p_surface$doublings)
#> [1] -0.89161  0.89161
term_heatmap(p_surface, show_support = TRUE)
```

![](feature_tour_files/figure-gfm/parameterized-predictions-1.png)<!-- -->

``` r
term_curves(term_predictions(pfit2, c("doublings", "sodium"),
                             values = list(doublings = seq(-1, 2, length.out = 30),
                                           sodium = c(130, 137, 142))))
```

![](feature_tour_files/figure-gfm/parameterized-predictions-2.png)<!-- -->

**Check:** `doublings` runs from about −0.9 to 0.9, i.e. creatinine from
0.54 to 1.86 mg/dL. Usually it is clearer to give the range with
`values`, as in the second plot.

## 10. Curves through a heatmap: `term_layer()` *(from p-ggeffects)*

A layer follows paths through the heatmap and predicts along them. One
path, a patient whose creatinine rises while sodium falls:

``` r
worsening <- term_layer(surface,
                        parameterization = list(creatinine = function(w) 1 + w,
                                                sodium = function(w) 140 - 5 * w),
                        terms = "w", values = list(w = seq(0, 1.5, by = 0.05)))
head(worsening)
#> Predictions over w ( 6 grid points, response scale )
#> Held constant: age = 60 (median), sex = Male (most common), diabetes =
#>   No (most common)
#> Critical values for 95% intervals (analytic): pointwise 1.969, simultaneous 2.547
#> 
#>     w estimate conf.low conf.high band.low band.high
#>  0.00    38.62    35.35     41.88    34.39     42.84
#>  0.05    38.40    35.06     41.74    34.07     42.72
#>  0.10    38.04    34.69     41.40    33.70     42.38
#>  0.15    37.62    34.31     40.94    33.33     41.91
#>  0.20    37.18    33.93     40.43    32.98     41.38
#>  0.25    36.74    33.56     39.92    32.63     40.85
term_heatmap(surface, layers = worsening)
```

![](feature_tour_files/figure-gfm/layer-one-1.png)<!-- -->

``` r
term_curves(worsening, show_direction = TRUE)
```

![](feature_tour_files/figure-gfm/layer-one-2.png)<!-- -->

``` r
plot(worsening)
```

![](feature_tour_files/figure-gfm/layer-one-3.png)<!-- -->

A family of paths, starting from different creatinine values:

``` r
paths <- term_layer(surface,
                    parameterization = list(creatinine = function(w, start) start + w,
                                            sodium = function(w) 140 - 5 * w),
                    terms = c("w", "start"),
                    values = list(w = seq(0, 1.5, by = 0.05), start = c(0.8, 1.2, 1.6)))
term_heatmap(surface, layers = paths)
```

![](feature_tour_files/figure-gfm/layer-family-1.png)<!-- -->

``` r
term_heatmap(surface, layers = list(worsening, paths), show_direction = FALSE)
```

![](feature_tour_files/figure-gfm/layer-family-2.png)<!-- -->

``` r
term_curves(paths)
```

![](feature_tour_files/figure-gfm/layer-family-3.png)<!-- -->

**Check:** the paths run from upper left to lower right, with arrows in
the first plot and without in the last; `term_curves()` shows LVEF along
each path.

## 11. Testing a reference curve: `test_simultaneous_predictions()` *(from p-ggeffects)*

Is a flat LVEF of 38% compatible with the curves, i.e. does it lie
inside all the simultaneous bands?

``` r
flat <- test_simultaneous_predictions(cuts, f = function(creatinine, sodium) 38 + 0 * creatinine)
flat
#> Simultaneous test: does the reference function fit inside the 95% band?
#> Held constant: age = 60 (median), sex = Male (most common), diabetes =
#>   No (most common)
#> max |t| = 3.24, p = 0.024
#> The reference lies inside the band at 100 of 120 grid points.
#> 
#>  creatinine sodium estimate reference band.low band.high inside
#>      0.5000    130    23.24        38   -8.622     55.11   TRUE
#>      0.5897    130    26.44        38    1.841     51.03   TRUE
#>      0.6795    130    29.47        38   11.540     47.39   TRUE
#>      0.7692    130    32.17        38   19.707     44.63   TRUE
#>      0.8590    130    34.38        38   25.347     43.40   TRUE
#>      0.9487    130    35.93        38   27.863     43.99   TRUE
#>      1.0385    130    36.66        38   28.358     44.97   TRUE
#>      1.1282    130    36.62        38   28.311     44.93   TRUE
#>      1.2179    130    36.07        38   28.207     43.93   TRUE
#>      1.3077    130    35.30        38   28.032     42.57   TRUE
attr(flat, "all_inside")
#> [1] FALSE
attr(flat, "p_value")
#> [1] 0.02409759
plot(flat)
```

![](feature_tour_files/figure-gfm/test-1.png)<!-- -->

``` r
plot(flat, type = "difference")
```

![](feature_tour_files/figure-gfm/test-2.png)<!-- -->

``` r
plot(test_simultaneous_predictions(surface, f = function(creatinine, sodium) 38 + 0 * creatinine))
```

![](feature_tour_files/figure-gfm/test-3.png)<!-- -->

**Check:** the p-value is below 0.05 exactly when `all_inside` is
`FALSE`. The difference plot of a surface is a heatmap in blue and red,
white where the prediction equals the reference.

## 12. Errors you should get

These calls are wrong on purpose; each should stop with a clear message.

``` r
term_heatmap(term_predictions(fit, "creatinine"))                        # one focal variable
#> Error: A heatmap needs two continuous focal variables.
term_predictions(fit, c("creatinine", "sodium"), ci_method = "delta")      # delta is for orm
#> Error: ci_method = "delta" is only for rms::orm models; for other models the analytic intervals already are delta-method intervals.
parameterize(fit, list(weight = function(t) t))                          # not in the model
#> Error: Cannot parameterize weight: not a predictor of the model.
term_layer(surface, list(age = function(t) t))                           # age is not on the heatmap
#> Error: `parameterization` may only give the heatmap variables creatinine, sodium.
test_simultaneous_predictions(cuts, f = function(cyl) cyl)               # unknown variable
#> Error: `f` uses cyl, which is not a variable of the predictions. Available: creatinine, sodium, age, sex, diabetes.
```

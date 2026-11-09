# termeffects feature demonstration: Cardiac systolic function and renal health

> **About this vignette.** Written by Lukas Bolsinger for p-ggeffects
> (original file `p_ggeffects_demonstration.Rmd`) and ported to
> termeffects by Brigitte Mueller. His text is kept; the code now uses
> the termeffects functions (`ggpredict()` became `term_predictions()`,
> `plot(..., "heatmap")` became `term_heatmap()`, layers became
> `term_layer()`), and held covariates are at their median or most
> common value. Also corrected: the eGFR formula is now the 2021 CKD-EPI
> equation (it said so before, but used a made-up linear formula);
> Figure 3 now plots eGFR (it plotted creatinine); the constant-eGFR
> curves in Figure 2 are now a `term_layer()` instead of hand-drawn
> contours; outdated limitations are updated. Natural splines (`ns()`)
> are replaced by restricted cubic splines (`rcs()`), as used in our
> team. Every change is marked in
> `validation/output/vignette_track_changes.html` (struck through:
> removed; underlined: added), or run
> `git diff 2de18bc -- vignettes/heart_failure_egfr.Rmd`.

# Introduction

This case study demonstrates a workflow for **termeffects** in which
fitting a regression model is deliberately separated from specifying
prediction domains and drawing figures. We analyze the `heartfailure`
dataset distributed with termeffects (originally from `dlookr`),
treating **left ventricular ejection fraction (LVEF)** as a continuous
outcome. We ask how its *adjusted association* with renal markers varies
across creatinine, serum sodium, and age.

The example illustrates six tools:

1.  **Spline interactions** allow a flexible two-variable response
    surface.
2.  **`term_predictions()`** produces adjusted predictions on regular
    grids for heatmaps.
3.  **`parameterize()`** transforms a *fitted* model into clinically
    meaningful prediction coordinates without refitting it.
4.  **Heatmaps and parameterized overlays** (`term_heatmap()`,
    `term_layer()`) display both the fitted surface and nonlinear
    constant-eGFR trajectories in the original predictor coordinates.
5.  **`term_curves(..., band = "simultaneous")`** demonstrates
    simultaneous uncertainty over a *specified prediction grid*.
6.  **`test_simultaneous_predictions(..., f = f0)`** demonstrates
    simultaneous comparison with a specified reference function,
    including a linear reference surface. This is not a general test of
    linearity.

**Scope.** Every plot has at most two continuous focal coordinates (or
one continuous coordinate and age strata). All other predictors are held
at stated reference values. In the model-fitting stage, no plotting
decisions are made. The resulting analysis is associative, not causal.
The data are observational and small (299 patients); apparent nonlinear
patterns may be unstable.

# Part I — Data and model (no plotting)

## 1. Load packages and inspect the dataset

``` r
knitr::opts_chunk$set(fig.width = 8, fig.height = 5, collapse = TRUE,
                      comment = "#>")

library(termeffects)
library(rms)
library(ggplot2)
```

``` r
data("heartfailure", package = "termeffects")
raw <- as.data.frame(heartfailure)
names(raw)
#>  [1] "age"               "anaemia"           "cpk_enzyme"       
#>  [4] "diabetes"          "ejection_fraction" "hblood_pressure"  
#>  [7] "platelets"         "creatinine"        "sodium"           
#> [10] "sex"               "smoking"           "time"             
#> [13] "death_event"
str(raw)
#> 'data.frame':    299 obs. of  13 variables:
#>  $ age              : int  75 55 65 50 65 90 75 60 65 80 ...
#>  $ anaemia          : Factor w/ 2 levels "No","Yes": 1 1 1 2 2 2 2 2 1 2 ...
#>  $ cpk_enzyme       : num  582 7861 146 111 160 ...
#>  $ diabetes         : Factor w/ 2 levels "No","Yes": 1 1 1 1 2 1 1 2 1 1 ...
#>  $ ejection_fraction: num  20 38 20 20 20 40 15 60 65 35 ...
#>  $ hblood_pressure  : Factor w/ 2 levels "No","Yes": 2 1 1 1 1 2 1 1 1 2 ...
#>  $ platelets        : num  265000 263358 162000 210000 327000 ...
#>  $ creatinine       : num  1.9 1.1 1.3 1.9 2.7 2.1 1.2 1.1 1.5 9.4 ...
#>  $ sodium           : num  130 136 129 137 116 132 137 131 138 133 ...
#>  $ sex              : Factor w/ 2 levels "Female","Male": 2 2 2 2 1 2 2 2 1 2 ...
#>  $ smoking          : Factor w/ 2 levels "No","Yes": 1 1 2 1 1 2 1 2 1 2 ...
#>  $ time             : int  4 6 7 7 8 8 10 10 10 10 ...
#>  $ death_event      : Factor w/ 2 levels "No","Yes": 2 2 2 2 2 2 2 2 2 2 ...
```

We use the conventional column names in the public heart-failure
clinical-records dataset.

``` r
required <- c("ejection_fraction", "creatinine", "sodium",
              "age", "sex", "diabetes")
stopifnot(all(required %in% names(raw)))
dat <- raw[stats::complete.cases(raw[, required]), required]
# Normalize the public dataset to the conventional names used below,
# and encode the categorical predictors as 0/1 for easy fixed-value
# conditioning in the parameterized predictions.
dat <- stats::setNames(dat, c("ejection_fraction", "serum_creatinine",
                              "serum_sodium", "age", "sex", "diabetes"))
dat$sex <- as.integer(dat$sex == "Male")
dat$diabetes <- as.integer(dat$diabetes == "Yes")
stopifnot(all(dat$sex %in% c(0, 1)))
stopifnot(all(dat$diabetes %in% c(0, 1)))
summary(dat)
#>  ejection_fraction serum_creatinine  serum_sodium        age       
#>  Min.   :14.00     Min.   :0.500    Min.   :113.0   Min.   :40.00  
#>  1st Qu.:30.00     1st Qu.:0.900    1st Qu.:134.0   1st Qu.:51.00  
#>  Median :38.00     Median :1.100    Median :137.0   Median :60.00  
#>  Mean   :38.08     Mean   :1.394    Mean   :136.6   Mean   :60.83  
#>  3rd Qu.:45.00     3rd Qu.:1.400    3rd Qu.:140.0   3rd Qu.:70.00  
#>  Max.   :80.00     Max.   :9.400    Max.   :148.0   Max.   :95.00  
#>       sex            diabetes     
#>  Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :1.0000   Median :0.0000  
#>  Mean   :0.6488   Mean   :0.4181  
#>  3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :1.0000   Max.   :1.0000
```

## 2. Fit one regression model

The outcome is LVEF (%). We allow a nonlinear creatinine–sodium
interaction and a nonlinear creatinine–age interaction, with modest
spline degrees of freedom. Sex and diabetes are adjustment variables.
This model is intentionally illustrative; its adequacy should be
assessed before substantive clinical interpretation.

``` r
model <- lm(
  ejection_fraction ~
    rcs(serum_creatinine, 4) +
    rcs(serum_sodium, 4) +
    rcs(age, 4) +
    rcs(serum_creatinine, 4) * rcs(serum_sodium, 4) +
    rcs(serum_creatinine, 4) * rcs(age, 4) +
    sex + diabetes,
  data = dat
)
summary(model)
#> 
#> Call:
#> lm(formula = ejection_fraction ~ rcs(serum_creatinine, 4) + rcs(serum_sodium, 
#>     4) + rcs(age, 4) + rcs(serum_creatinine, 4) * rcs(serum_sodium, 
#>     4) + rcs(serum_creatinine, 4) * rcs(age, 4) + sex + diabetes, 
#>     data = dat)
#> 
#> Residuals:
#>     Min      1Q  Median      3Q     Max 
#> -22.999  -7.121  -1.660   6.816  41.268 
#> 
#> Coefficients:
#>                                                                                 Estimate
#> (Intercept)                                                                   -1.294e+03
#> rcs(serum_creatinine, 4)serum_creatinine                                       1.399e+03
#> rcs(serum_creatinine, 4)serum_creatinine'                                     -1.549e+04
#> rcs(serum_creatinine, 4)serum_creatinine''                                     3.781e+04
#> rcs(serum_sodium, 4)serum_sodium                                               9.755e+00
#> rcs(serum_sodium, 4)serum_sodium'                                             -1.588e+01
#> rcs(serum_sodium, 4)serum_sodium''                                             2.357e+01
#> rcs(age, 4)age                                                                 8.325e-01
#> rcs(age, 4)age'                                                               -3.019e+00
#> rcs(age, 4)age''                                                               9.188e+00
#> sex                                                                           -3.435e+00
#> diabetes                                                                      -2.704e-01
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium     -1.024e+01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium     1.135e+02
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium   -2.766e+02
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium'     1.647e+01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium'   -1.432e+02
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium'   3.335e+02
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium''   -2.031e+01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium''  -1.403e+02
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium''  5.109e+02
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age                       -8.199e-01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age                       5.176e+00
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age                     -1.160e+01
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age'                       3.628e+00
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age'                     -4.812e+01
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age'                     1.243e+02
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age''                     -1.142e+01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age''                     1.878e+02
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age''                   -4.983e+02
#>                                                                               Std. Error
#> (Intercept)                                                                    1.175e+03
#> rcs(serum_creatinine, 4)serum_creatinine                                       1.336e+03
#> rcs(serum_creatinine, 4)serum_creatinine'                                      1.861e+04
#> rcs(serum_creatinine, 4)serum_creatinine''                                     4.796e+04
#> rcs(serum_sodium, 4)serum_sodium                                               8.954e+00
#> rcs(serum_sodium, 4)serum_sodium'                                              2.337e+01
#> rcs(serum_sodium, 4)serum_sodium''                                             1.721e+02
#> rcs(age, 4)age                                                                 3.877e+00
#> rcs(age, 4)age'                                                                1.441e+01
#> rcs(age, 4)age''                                                               4.579e+01
#> sex                                                                            1.483e+00
#> diabetes                                                                       1.429e+00
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium      1.015e+01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium     1.412e+02
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium    3.640e+02
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium'     2.637e+01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium'    3.706e+02
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium'   9.592e+02
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium''    1.915e+02
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium''   2.658e+03
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium''  6.867e+03
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age                        4.445e+00
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age                       7.083e+01
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age                      1.872e+02
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age'                       1.635e+01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age'                      2.482e+02
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age'                     6.514e+02
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age''                      5.154e+01
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age''                     7.588e+02
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age''                    1.982e+03
#>                                                                               t value
#> (Intercept)                                                                    -1.102
#> rcs(serum_creatinine, 4)serum_creatinine                                        1.047
#> rcs(serum_creatinine, 4)serum_creatinine'                                      -0.833
#> rcs(serum_creatinine, 4)serum_creatinine''                                      0.788
#> rcs(serum_sodium, 4)serum_sodium                                                1.089
#> rcs(serum_sodium, 4)serum_sodium'                                              -0.679
#> rcs(serum_sodium, 4)serum_sodium''                                              0.137
#> rcs(age, 4)age                                                                  0.215
#> rcs(age, 4)age'                                                                -0.209
#> rcs(age, 4)age''                                                                0.201
#> sex                                                                            -2.316
#> diabetes                                                                       -0.189
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium      -1.008
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium      0.804
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium    -0.760
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium'      0.625
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium'    -0.386
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium'    0.348
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium''    -0.106
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium''   -0.053
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium''   0.074
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age                        -0.184
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age                        0.073
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age                      -0.062
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age'                        0.222
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age'                      -0.194
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age'                      0.191
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age''                      -0.222
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age''                      0.248
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age''                    -0.251
#>                                                                               Pr(>|t|)
#> (Intercept)                                                                     0.2715
#> rcs(serum_creatinine, 4)serum_creatinine                                        0.2959
#> rcs(serum_creatinine, 4)serum_creatinine'                                       0.4058
#> rcs(serum_creatinine, 4)serum_creatinine''                                      0.4313
#> rcs(serum_sodium, 4)serum_sodium                                                0.2769
#> rcs(serum_sodium, 4)serum_sodium'                                               0.4974
#> rcs(serum_sodium, 4)serum_sodium''                                              0.8912
#> rcs(age, 4)age                                                                  0.8301
#> rcs(age, 4)age'                                                                 0.8342
#> rcs(age, 4)age''                                                                0.8411
#> sex                                                                             0.0213
#> diabetes                                                                        0.8501
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium       0.3143
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium      0.4221
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium     0.4479
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium'      0.5328
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium'     0.6995
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium'    0.7284
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium''     0.9156
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium''    0.9579
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium''   0.9407
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age                         0.8538
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age                        0.9418
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age                       0.9506
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age'                        0.8246
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age'                       0.8464
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age'                      0.8488
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age''                       0.8247
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age''                      0.8047
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age''                     0.8016
#>                                                                                
#> (Intercept)                                                                    
#> rcs(serum_creatinine, 4)serum_creatinine                                       
#> rcs(serum_creatinine, 4)serum_creatinine'                                      
#> rcs(serum_creatinine, 4)serum_creatinine''                                     
#> rcs(serum_sodium, 4)serum_sodium                                               
#> rcs(serum_sodium, 4)serum_sodium'                                              
#> rcs(serum_sodium, 4)serum_sodium''                                             
#> rcs(age, 4)age                                                                 
#> rcs(age, 4)age'                                                                
#> rcs(age, 4)age''                                                               
#> sex                                                                           *
#> diabetes                                                                       
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium      
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium     
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium    
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium'     
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium'    
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium'   
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium''    
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium''   
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium''  
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age                        
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age                       
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age                      
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age'                       
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age'                      
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age'                     
#> rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age''                      
#> rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age''                     
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age''                    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error: 11.48 on 269 degrees of freedom
#> Multiple R-squared:  0.1504, Adjusted R-squared:  0.0588 
#> F-statistic: 1.642 on 29 and 269 DF,  p-value: 0.0235
```

**Modeling caveat:** Two 3-by-3 interactions can still be ambitious for
299 observations. Inspect residuals, collinearity, and prediction
uncertainty; if necessary, remove one interaction or reduce basis
dimensions. Do not choose the model merely to produce visually
interesting heatmaps.

``` r
plot(fitted(model), resid(model),
     xlab = "Fitted LVEF (%)", ylab = "Residual",
     main = "Residuals versus fitted values")
abline(h = 0, lty = 2)
```

![](heart_failure_egfr_files/figure-gfm/model-diagnostics-1.png)<!-- -->

# Part II — Prediction definitions (no plotting)

## 3. Choose explicit adjustment values

All non-focal predictors are **fixed at a common reference profile**,
rather than averaged over the patient distribution. The reference
profile is descriptive, not necessarily a typical individual.

``` r
ref_age <- median(dat$age)
ref_sodium <- median(dat$serum_sodium)
ref_creatinine <- median(dat$serum_creatinine)
ref_sex <- 1                 # Male for eGFR illustrations
ref_diabetes <- 0
```

**Prediction API note:** `term_predictions()` holds non-focal predictors
at their median (numeric) or most common value (categorical) unless they
are given in `at`. Below, sex and diabetes are passed explicitly in
`at`, so every prediction uses the same reference profile; each result
prints the values it used.

## 4. Predictions on the original coordinates

``` r
# The output objects are data frames containing predicted values,
# pointwise confidence limits, simultaneous band limits, and the
# focal-variable grid coordinates.
reference_profile <- list(sex = ref_sex, diabetes = ref_diabetes)

pred_creat_sodium <- term_predictions(
  model,
  terms = c("serum_creatinine", "serum_sodium"),
  at = c(reference_profile, list(age = ref_age))
)

pred_age_creat <- term_predictions(
  model,
  terms = c("age", "serum_creatinine"),
  at = c(reference_profile, list(serum_sodium = ref_sodium))
)
```

## 5. Define the nonlinear eGFR transformation

For demonstration we use the **2021 CKD-EPI creatinine equation**, which
depends on serum creatinine (mg/dL), age (years), and sex. This is an
estimated kidney-function measure, not a direct measurement. Its
application to these historical records is retrospective.

``` math
\text{eGFR} = 142 \cdot \min(C/\kappa, 1)^{\alpha} \cdot \max(C/\kappa, 1)^{-1.200} \cdot 0.9938^{A} \cdot (1.012 \text{ if female}),
```
with $`\kappa = 0.7`$, $`\alpha = -0.241`$ for women and
$`\kappa = 0.9`$, $`\alpha = -0.302`$ for men.

``` r
# 2021 CKD-EPI creatinine equation (sex: 1 = male, 0 = female).
egfr_ckdepi <- function(creatinine, age, sex) {
  kappa <- ifelse(sex == 1, 0.9, 0.7)
  alpha <- ifelse(sex == 1, -0.302, -0.241)
  142 * pmin(creatinine / kappa, 1)^alpha * pmax(creatinine / kappa, 1)^-1.200 *
    0.9938^age * ifelse(sex == 1, 1, 1.012)
}

# Inverse transformation used to define the parameterized model: the
# creatinine at which a patient of this age and sex has the given eGFR.
creatinine_from_egfr <- function(egfr, age, sex) {
  kappa <- ifelse(sex == 1, 0.9, 0.7)
  alpha <- ifelse(sex == 1, -0.302, -0.241)
  ratio <- egfr / (142 * 0.9938^age * ifelse(sex == 1, 1, 1.012))
  kappa * ifelse(ratio >= 1, ratio^(1 / alpha), ratio^(1 / -1.200))
}

# check: the inverse undoes the equation
creatinine_from_egfr(egfr_ckdepi(c(0.6, 1, 2.5), age = 65, sex = 1), age = 65, sex = 1)
#> [1] 0.6 1.0 2.5
```

The sex-dependent transformation is **not** obtained by averaging sex
first. For clarity, all eGFR plots below use a fixed male reference
profile. The original model still includes sex as an adjustment
variable.

## 6. Parameterize the already-fitted model

This is the central termeffects step. We do **not** refit `model`: we
change only the inputs used for prediction. For a given eGFR and age,
the parameterization computes the corresponding serum creatinine.

``` r
pmodel_egfr <- parameterize(
  model,
  parameterizations = list(
    serum_creatinine = function(egfr, age) {
      creatinine_from_egfr(egfr, age, sex = ref_sex)
    },
    sex = function(egfr, age) {
      rep(ref_sex, length.out = max(length(egfr), length(age)))
    }
  )
)

print(pmodel_egfr, show_history = TRUE)
#> Parameterized model
#> 
#> Model:
#> 
#> Call:
#> lm(formula = ejection_fraction ~ rcs(serum_creatinine, 4) + rcs(serum_sodium, 
#>     4) + rcs(age, 4) + rcs(serum_creatinine, 4) * rcs(serum_sodium, 
#>     4) + rcs(serum_creatinine, 4) * rcs(age, 4) + sex + diabetes, 
#>     data = dat)
#> 
#> Coefficients:
#>                                                                   (Intercept)  
#>                                                                    -1.294e+03  
#>                                      rcs(serum_creatinine, 4)serum_creatinine  
#>                                                                     1.399e+03  
#>                                     rcs(serum_creatinine, 4)serum_creatinine'  
#>                                                                    -1.549e+04  
#>                                    rcs(serum_creatinine, 4)serum_creatinine''  
#>                                                                     3.781e+04  
#>                                              rcs(serum_sodium, 4)serum_sodium  
#>                                                                     9.755e+00  
#>                                             rcs(serum_sodium, 4)serum_sodium'  
#>                                                                    -1.588e+01  
#>                                            rcs(serum_sodium, 4)serum_sodium''  
#>                                                                     2.357e+01  
#>                                                                rcs(age, 4)age  
#>                                                                     8.325e-01  
#>                                                               rcs(age, 4)age'  
#>                                                                    -3.019e+00  
#>                                                              rcs(age, 4)age''  
#>                                                                     9.188e+00  
#>                                                                           sex  
#>                                                                    -3.435e+00  
#>                                                                      diabetes  
#>                                                                    -2.704e-01  
#>     rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium  
#>                                                                    -1.024e+01  
#>    rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium  
#>                                                                     1.135e+02  
#>   rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium  
#>                                                                    -2.766e+02  
#>    rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium'  
#>                                                                     1.647e+01  
#>   rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium'  
#>                                                                    -1.432e+02  
#>  rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium'  
#>                                                                     3.335e+02  
#>   rcs(serum_creatinine, 4)serum_creatinine:rcs(serum_sodium, 4)serum_sodium''  
#>                                                                    -2.031e+01  
#>  rcs(serum_creatinine, 4)serum_creatinine':rcs(serum_sodium, 4)serum_sodium''  
#>                                                                    -1.403e+02  
#> rcs(serum_creatinine, 4)serum_creatinine'':rcs(serum_sodium, 4)serum_sodium''  
#>                                                                     5.109e+02  
#>                       rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age  
#>                                                                    -8.199e-01  
#>                      rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age  
#>                                                                     5.176e+00  
#>                     rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age  
#>                                                                    -1.160e+01  
#>                      rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age'  
#>                                                                     3.628e+00  
#>                     rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age'  
#>                                                                    -4.812e+01  
#>                    rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age'  
#>                                                                     1.243e+02  
#>                     rcs(serum_creatinine, 4)serum_creatinine:rcs(age, 4)age''  
#>                                                                    -1.142e+01  
#>                    rcs(serum_creatinine, 4)serum_creatinine':rcs(age, 4)age''  
#>                                                                     1.878e+02  
#>                   rcs(serum_creatinine, 4)serum_creatinine'':rcs(age, 4)age''  
#>                                                                    -4.983e+02  
#> 
#> 
#> Parameters:
#> [1] "egfr" "age" 
#> 
#> Parameterizations:
#> $serum_creatinine <- ({     creatinine_from_egfr(egfr, age, sex = ref_sex) }); (egfr, age)
#> $sex <- ({     rep(ref_sex, length.out = max(length(egfr), length(age))) }); (egfr, age)
#> 
#> 
#> Parameterization history:
#> 
#> Stage 1                                                                    
#> ---------------------------------------------------------------------------
#> 
#> serum_creatinine [{    creatinine_from_egfr(egfr, age, sex = ref_sex)}] <── egfr, age
#> sex [{    rep(ref_sex, length.out = max(length(egfr), length(age)))}] <── egfr, age
parameter_variables(pmodel_egfr)
#> [1] "egfr" "age"
```

The transformed model depends on `egfr` and `age`, in addition to the
original untransformed predictors such as sodium and diabetes. Its
fitted coefficients remain those of `model`.

``` r
pred_egfr <- term_predictions(
  pmodel_egfr,
  terms = c("egfr", "age"),
  values = list(egfr = seq(20, 90, by = 2), age = c(50, 65, 80)),
  at = list(serum_sodium = ref_sodium, diabetes = ref_diabetes)
)
pred_egfr
#> Predictions over egfr, age ( 108 grid points, response scale )
#> Held constant: serum_sodium = 137 (set by user), diabetes = 0 (set by
#>   user)
#> Critical values for 95% intervals (analytic): pointwise 1.969, simultaneous 3.000
#> 
#>  egfr age estimate conf.low conf.high band.low band.high
#>    20  50    38.95    29.77     48.14    24.95     52.96
#>    22  50    37.51    28.97     46.05    24.50     50.52
#>    24  50    36.28    28.14     44.43    23.88     48.69
#>    26  50    35.24    27.31     43.17    23.16     47.32
#>    28  50    34.37    26.55     42.18    22.46     46.28
#>    30  50    33.67    25.93     41.41    21.87     45.47
#>    32  50    33.13    25.47     40.79    21.45     44.81
#>    34  50    32.72    25.17     40.28    21.21     44.24
#>    36  50    32.43    25.02     39.85    21.14     43.73
#>    38  50    32.24    25.01     39.48    21.22     43.27
#> ... with 98 more rows
```

# Part III — Visualization

## Figure 1 — Nonlinear creatinine × sodium heatmap

**Focal variables:** creatinine and sodium. **Other predictors:**
adjusted at reference values. **Question:** Does the adjusted
association of creatinine with LVEF vary by sodium concentration?

``` r
term_heatmap(pred_creat_sodium, show_support = TRUE, show_data = TRUE)
```

<figure>
<img
src="heart_failure_egfr_files/figure-gfm/plot-heatmap-1.png"
alt="Adjusted LVEF surface in original predictor coordinates. Faded: outside the region holding 97.5% of the patients." />
<figcaption aria-hidden="true">Adjusted LVEF surface in original
predictor coordinates. Faded: outside the region holding 97.5% of the
patients.</figcaption>
</figure>

The density overlay is important: a smooth model can extrapolate into
combinations of creatinine and sodium that are poorly represented in the
sample. Interpret those regions cautiously.

## Figure 2 — Creatinine × age heatmap with constant-eGFR curves

**Focal variables:** creatinine and age. **Overlay:** constant eGFR at
30, 45, 60, and 90 mL/min/1.73 m². **Other predictors:** adjusted; sex
fixed to male.

This is an annotated original-coordinate heatmap. The curves are **not**
extra axes or extra facets: each curve is a one-dimensional path in the
same two-dimensional predictor plane.

``` r
# A layer follows curves through the heatmap: age runs along each curve
# (t), and creatinine is the value giving that curve's eGFR at that age.
egfr_layer <- term_layer(
  pred_age_creat,
  parameterization = list(
    age = function(t) t,
    serum_creatinine = function(t, egfr) creatinine_from_egfr(egfr, age = t, sex = ref_sex)
  ),
  terms = c("t", "egfr"),
  values = list(t = seq(min(dat$age), max(dat$age), by = 1), egfr = c(30, 45, 60, 90))
)
```

``` r
term_heatmap(pred_age_creat, layers = egfr_layer, show_support = TRUE) +
  labs(x = "Age (years)", y = "Serum creatinine (mg/dL)")
```

<figure>
<img
src="heart_failure_egfr_files/figure-gfm/plot-egfr-overlays-1.png"
alt="Original-coordinate prediction heatmap with nonlinear constant-eGFR curves." />
<figcaption aria-hidden="true">Original-coordinate prediction heatmap
with nonlinear constant-eGFR curves.</figcaption>
</figure>

**How this connects to termeffects layers:** The overlaid curves are
level sets of a nonlinear parameterization of the original predictors,
computed by `term_layer()`. The layer also holds the predictions along
each curve, with a band that is simultaneous over the curves:

``` r
term_curves(egfr_layer) + labs(x = "Age (years)")
```

<figure>
<img
src="heart_failure_egfr_files/figure-gfm/plot-egfr-layer-curves-1.png"
alt="Predicted LVEF along the constant-eGFR curves, by age." />
<figcaption aria-hidden="true">Predicted LVEF along the constant-eGFR
curves, by age.</figcaption>
</figure>

## Figure 3 — LVEF versus eGFR with simultaneous confidence bands

**Focal variables:** eGFR and age (three representative ages). **Other
predictors:** adjusted; sex fixed to male. **Question:** How do adjusted
predictions vary over kidney-function values, and how does simultaneous
uncertainty differ from pointwise uncertainty?

``` r
term_curves(pred_egfr, band = "pointwise")
```

<figure>
<img
src="heart_failure_egfr_files/figure-gfm/plot-egfr-pointwise-1.png"
alt="Pointwise confidence intervals for age-stratified eGFR predictions." />
<figcaption aria-hidden="true">Pointwise confidence intervals for
age-stratified eGFR predictions.</figcaption>
</figure>

``` r
term_curves(pred_egfr, band = "simultaneous")
```

<figure>
<img
src="heart_failure_egfr_files/figure-gfm/plot-egfr-simultaneous-1.png"
alt="Simultaneous confidence bands over the specified eGFR-by-age prediction grid." />
<figcaption aria-hidden="true">Simultaneous confidence bands over the
specified eGFR-by-age prediction grid.</figcaption>
</figure>

The simultaneous band is calibrated for the **finite generated
prediction grid**, not automatically for every point in a continuous
domain. The `ci_scope` setting controls whether uncertainty is
simultaneous across the full set of age curves (`"global"`, the default)
or separately within each age group (`"curve"`). State that choice
explicitly in any published figure.

## Figure 4 — Simultaneous test against a specified linear surface

We test the **fixed, prespecified** linear reference surface
``` math
f_0(C,S) = 40 - 3(C-1) + 0.2(S-137),
```
with LVEF measured in percentage points. The test is simultaneous over
the **creatinine–sodium grid**, with other predictors adjusted as in
Figure 1. This is a test of equality to this *particular* linear
function, not a general test that the response is linear.

``` r
f0 <- function(serum_creatinine, serum_sodium) {
  40 - 3 * (serum_creatinine - 1) + 0.2 * (serum_sodium - 137)
}

sim_test <- test_simultaneous_predictions(
  pred_creat_sodium,
  f = f0
)
sim_test
#> Simultaneous test: does the reference function fit inside the 95% band?
#> Held constant: age = 60 (set by user), sex = 1 (set by user), diabetes
#>   = 0 (set by user)
#> max |t| = 3.41, p = 0.042
#> The reference lies inside the band at 899 of 900 grid points.
#> 
#>  serum_creatinine serum_sodium estimate reference band.low band.high inside
#>            0.5000          113  -52.738     36.70 -312.422    206.95   TRUE
#>            0.8069          113    9.227     35.78  -70.827     89.28   TRUE
#>            1.1138          113   40.598     34.86  -23.371    104.57   TRUE
#>            1.4207          113   37.929     33.94    3.974     71.88   TRUE
#>            1.7276          113   34.150     33.02    6.444     61.85   TRUE
#>            2.0345          113   29.838     32.10    0.817     58.86   TRUE
#>            2.3414          113   25.122     31.18   -4.917     55.16   TRUE
#>            2.6483          113   20.131     30.26  -15.914     56.18   TRUE
#>            2.9552          113   14.993     29.33  -35.050     65.04   TRUE
#>            3.2621          113    9.823     28.41  -58.722     78.37   TRUE
plot(sim_test)
```

![](heart_failure_egfr_files/figure-gfm/simultaneous-reference-test-1.png)<!-- -->

The reference function’s argument names match the **original
predictors**. The plotted deviations and simultaneous band are only for
the grid that was generated.

## 7. Limits of simultaneous inference after marginalization

The `test_simultaneous_predictions()` implementation works on a
`term_predictions` or `term_layer` object, **not on the entire original
regression function**. These objects store predictions on a chosen focal
grid after the remaining model predictors have been fixed at
representative values. This imposes important limitations:

- **Restricted domain.** The simultaneous test concerns only the
  generated grid of focal terms, not all combinations of predictors in
  the original model. A conclusion for the creatinine–sodium surface
  does not establish a corresponding conclusion for every age, sex, or
  diabetes status.
- **Marginalization is baked in.** The non-focal variables have already
  been fixed at the prediction object’s reference values. The test
  cannot recover their full variation, interactions, or the original
  higher-dimensional domain from that object.
- **Scope depends on the prediction object.** A `term_layer` can
  restrict the domain further to selected parameterized curves. A
  simultaneous band over such a layer does not provide simultaneous
  coverage over the entire background heatmap or over unplotted curves.
- **Finite-grid inference.** The critical value is computed for the
  generated prediction points; it is not automatically valid uniformly
  over a continuous region between those points.
- **A specified reference is not a general shape test.** Comparing
  against one linear `f0` only assesses that particular surface, at the
  chosen adjustment values. Rejecting it does not prove nonlinearity,
  and failing to reject it does not establish linearity. A fully general
  test of linearity would need to account for the family of possible
  linear functions and any uncertainty in estimated reference
  coefficients.

Thus, the demonstration above is intentionally a **simultaneous test
against a fixed reference function on a marginalized grid**, rather than
a global test of linearity of the fitted multivariable regression. The
marginalization choices and grid must be reported alongside the result.

# Part IV — Interpretation and reproducibility

## 8. Questions to answer after rendering

- **Figure 1:** Is there evidence of a nonlinear creatinine–sodium
  interaction? Which regions have adequate observational support?
- **Figure 2:** Along a constant-eGFR curve, does predicted LVEF change
  appreciably with age? If so, eGFR alone does not capture the fitted
  age–creatinine association.
- **Figure 3:** How much wider are simultaneous bands than pointwise
  intervals? Which portions of the prediction domain are too uncertain
  for strong conclusions?
- **Figure 4:** Is the specified linear reference surface compatible
  with the fitted predictions over the grid? What does the test
  establish—and what does it not establish—about linearity?

## 9. Important qualifications

1.  **No causal interpretation.** These are conditional associations in
    an observational dataset.
2.  **Small sample.** Nonlinear interactions may be weakly identified.
    Consider reduced splines, bootstrap stability checks, and
    out-of-sample validation.
3.  **Sparse support.** Limit grids and overlays to observed or
    well-supported regions; not all constant-eGFR curves are represented
    in the data.
4.  **LVEF is the response.** The example does not use death or
    follow-up time, even though those fields are available.
5.  **Consistent adjustment.** Verify the adjustment profile printed
    with each `term_predictions()` result, particularly for sex and
    diabetes. Fix sex before interpreting sex-specific eGFR contours.
6.  **Separate fitting from presentation.** To change the grid, eGFR
    transformation, or plotting style, update Parts II–III. The fitted
    `model` in Part I remains unchanged.

## 10. Optional extension

If a later dataset includes **height and weight**, the *same conceptual
workflow* can be extended to a height × weight heatmap with constant-BMI
and constant-body-surface-area curves. Do not add those variables to
this dataset: they are not recorded in `heartfailure`.

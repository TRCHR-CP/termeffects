# termeffects 0.0.0.9000

* New `term_grid()`: builds a prediction grid in which the focal variables
  vary and all other predictors are held at typical values (or values you
  choose with `at`).
* New `term_predictions()`: predictions over such a grid with pointwise
  confidence intervals and simultaneous confidence bands. Works for `lm`,
  `glm`, `rms::orm`, `lme4::lmer`/`glmer`, `survival::coxph` and `mgcv::gam`.
* New `term_heatmap()`, `term_curves()` and `plot()`: ggplot2 heatmaps of a
  prediction surface and curves ("cuts") through it with their bands.
* New vignette "Seeing non-linear interactions: heatmaps and cuts".

# termeffects 0.0.0.9000

* New `term_grid()`: builds a prediction grid in which the focal variables
  vary and all other predictors are held at their median (numeric) or most
  common value (categorical), or at values you choose with `at`.
* New `term_predictions()`: predictions over such a grid with pointwise
  confidence intervals and simultaneous confidence bands. Works for `lm`,
  `glm` (any family, with offsets; also brglm2 and `geepack::geeglm`),
  `logistf`, `rms::orm`/`lrm`/`cph`/`ols`, `lme4::lmer`/`glmer`,
  `nlme::lme`, `survival::coxph` (also stratified, and with time-dependent
  covariates in (start, stop] format) and `mgcv::gam`.
  Models using rms's `%ia%` outside an rms fit are refused, because their
  spline knots would be recomputed from the prediction grid.
  Intervals are analytic (default) or bootstrap (`ci_method = "bootstrap"`);
  for `orm` also the delta method with all intercepts (`ci_method = "delta"`).
  `ci_scope` makes the band simultaneous per facet or per curve.
* New `term_heatmap()`, `term_curves()` and `plot()`: ggplot2 heatmaps of a
  prediction surface and curves ("cuts") through it with their bands. Heatmaps
  can show the band width or standard error, fade regions with few data, and
  overlay curves.
  The subtitle lists all covariates held constant; with `show_data = TRUE`
  the heatmap caption says how many observations lie outside the plot.
* New `term_figure()`: the plots together with a table of the held covariates
  and a caption describing the model and the figure.
* From p-ggeffects (Lukas Bolsinger): `parameterize()` and helpers for
  parameterized models; `term_layer()` for predictions along curves through a
  heatmap; `test_simultaneous_predictions()` for comparing a reference
  function with the simultaneous band (now with a p-value); `vcov()` of
  predictions.
* New data set `heartfailure`.
* Vignettes: "Seeing non-linear interactions: heatmaps and cuts" (how-to),
  "A tour of every feature" (every function once, with what to check),
  "Supported model types, with examples", and
  the student's four vignettes ported to termeffects.

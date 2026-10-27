# termeffects development log

A high-level record of what was changed and why, newest first. NEWS.md is the
user-facing summary; this log also covers decisions, validation and tooling.

---

## 2026-10-27 — orm models, plots and the how-to vignette

### Added
- **`rms::orm` support** (Brigitte's main model for bounded continuous
  outcomes). The default scale is the mean outcome, via `Hmisc::Mean()`.
  Bands are built on rms's linear predictor (reference intercept plus
  covariate effects) and transformed, so they stay within the outcome's
  observed range. This is the same convention as
  `rms::Predict(fun = Mean(fit))`, and is documented in the help page.
- **Plots (ggplot2):**
  - `term_heatmap()`: the surface, with contour lines. Optionally adds
    dashed lines marking cuts and the observed data points (to show where
    the surface extrapolates). A third categorical focal variable becomes
    panels.
  - `term_curves()`: cuts through the surface, one curve per value of the
    second variable, with the simultaneous band, pointwise intervals, both
    or none.
  - `plot()`: picks between the two.
  - Every plot states in its subtitle which variables are held constant, and
    at what values, and in its caption which interval is shown. Axis and
    legend labels describe the scale (e.g. "Mean score", "Probability of …").
- **Vignette** "Seeing non-linear interactions: heatmaps and cuts": a first
  draft of the how-to for PIs and analysts, using an orm example (fit, heatmap,
  cuts, pointwise vs. simultaneous in plain language, other models, caveats).

### Decisions
- Model fitting is unchanged: users fit with rms, stats, survival, lme4 or
  mgcv as usual. termeffects only works with the fitted model.
- ggplot2 is used for all plots (Imports). The only other hard dependencies
  are base R's `stats` and `utils`. rms/Hmisc, lme4, survival and mgcv are
  optional (Suggests), needed only by people who use those models.

### Tests and validation
- 77 tests, all passing. orm predictions, standard errors and intervals
  match rms's own `predict()`/`Mean()` exactly. The plot tests check the
  layers, labels and errors.
- Coverage simulation now includes orm: 98.3% for a 95% band. This is
  conservative (slightly wide, never too narrow). Likely reason: the mean
  averages over all intercepts, so it is estimated more precisely than the
  single reference intercept whose uncertainty the band includes.
- `validation/07_example_plots.R`: renders an example orm heatmap and cut
  plot to `validation/output/plots/` for visual checking.
- `R CMD check`: 0 errors. The only warning (qpdf) and note (timestamp) come
  from the Docker container, not the package.

---

## 2026-10-27 — Core engine: predictions and simultaneous bands

### Decisions
- **No dependency on ggeffects.** ggeffects will not be maintained in future,
  so termeffects does not use it, not even its public functions. We also do
  not use marginaleffects, whose interface changes often.
- **Own small prediction engine.** For each supported model type, predictions
  and their joint uncertainty come directly from the model's design matrix
  and coefficient covariance (base R `model.matrix()`, `coef()`, `vcov()`).
  The package currently depends only on base R (`stats`, `utils`). This keeps
  it stable when other packages are updated.
- **Bands are built on the link scale, then transformed.** This fixes the
  main error found in the student's code (see Validation), and keeps
  probability bands inside 0–1.
- **Simultaneous over the chosen grid, others held constant.** A band covers
  the whole curve or surface over the grid points at once, with the other
  predictors fixed at typical values. This matches the package goal (heatmaps
  and cuts through surfaces "holding everything else constant").

### Added
- `term_grid()`: builds the grid. Focal variables span their observed range
  (or values you give). Other predictors are held at their mean (numbers) or
  reference level (categories), or at values you set with `at`.
- `term_predictions()`: estimate, pointwise interval and simultaneous band at
  each grid point, on the response scale (default) or the link scale. A
  hand-made grid (`newdata`) can trace any curve through a surface; this is
  the basis for "cuts" and, later, the student's parameterized curves.
- Supported models: `lm`, `glm`, `lme4::lmer`/`glmer` (population-level),
  `survival::coxph` (hazard ratio vs. an average subject), `mgcv::gam`.
  Clear errors for unsupported cases (other classes, offsets, Cox `strata()`).
- Variables turned into factors inside the formula (e.g. `factor(am)`) are
  treated as categories.

### Tests (53, all passing)
- Estimates and standard errors match each model's own `predict()` for lm,
  glm, lmer, glmer, coxph and gam.
- Critical value matches exact answers (single point, independent points)
  and `mvtnorm::qmvt` for correlated spline predictions.
- GLM bands equal the textbook link-scale construction and stay in 0–1.

### Validation
- `validation/06_termeffects_coverage.R`: across 300 simulated datasets the
  simultaneous band covered the whole true curve set 95.3% (linear model)
  and 95.7% (logistic model) of the time; target 95%. Pointwise intervals
  covered it only 49% and 72% of the time, showing why the band is needed.
- `R CMD check`: 0 errors, 0 warnings.

### Package housekeeping
- DESCRIPTION: removed `ggeffects`, `insight`, `ggplot2` (ggplot2 returns with
  the plotting step); test-only packages listed under Suggests.
- `.Rbuildignore`: this log, `Claude_memory.md`, `how_i_started.R`.

---

## 2026-10-27 — Review of the student's fork (no package code yet)

### Validation scripts (`validation/`)
- `00_render_student_vignettes.R` (new): renders the student's four
  vignettes against his code into `validation/output/vignettes/`, showing
  errors inline. Renders from a temporary copy so the archive is untouched.
- `01_run_student_tests.R` (rewritten): runs the full ggeffects test suite on
  clean copies of the fork and of the upstream version it started from, and
  lists failures the fork introduced. The old version wrote files into the
  archive and stopped after 10 failures.
- `04_check_glm.R`: added the bootstrap critical-value check.

### Findings about the student's code
- Linear models: simultaneous critical values are correct.
- Logistic/GLM: bands mix scales (response-scale prediction ± link-scale SE);
  the bootstrap critical value is ~0.7 instead of ~2.7.
- Cox models: `ggpredict()` fails on any `coxph` (regression vs. upstream).
- `ggaverage()`, `ggemmeans()`, `ggeffect()` are not registered and fail.
- Upstream fails 1 of 220 test blocks; the fork fails 121 of 228.
- Archive: removed two stray snapshot files left by an earlier test run, so
  it matches GitHub again.

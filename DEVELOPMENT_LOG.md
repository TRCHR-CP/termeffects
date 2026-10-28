# termeffects development log

A high-level record of what was changed and why, newest first. NEWS.md is the
user-facing summary; this log also covers decisions, validation and tooling.

---

## 2026-10-27 — Feature tour vignette for testing

New vignette `vignettes/feature_tour.Rmd` ("A tour of every feature"): uses
every exported function once on the `heartfailure` data. Each section says what
the reader should see ("Check:"), so it can serve as a manual test. It also
lists the automatic checks (`devtools::test()`, `check()`, `build_vignettes()`)
and ends with calls that must fail with a clear message. Features from
p-ggeffects are marked.

Problems found by writing it, now fixed:

* `term_heatmap(layers = list(a, b))` failed when the layers had different
  variables (e.g. a single path and a family of paths): rows were combined
  with all their columns. Only the drawn columns are combined now; test added.
* The `term_figure()` caption read "enter the model as creatinine as
  rcs(creatinine, 4)". Now: "In the model, creatinine enters as
  rcs(creatinine, 4) and sodium as rcs(sodium, 4)."
* The difference heatmap of `plot(test)` did not say which values the other
  covariates were held at. It now has the test result as title and the held
  values as subtitle, like the curves version.

Noted, not changed: a layer parameterization whose input is named like a held
covariate (e.g. `function(age)` while age is held at 60) uses the held value.
This is consistent with how held values work elsewhere.

---

## 2026-10-27 — The student's work ported; orm intervals checked; figures

Git history on branch `port-student-work` (newest last):
`374612e` engine → `2de18bc` student's four vignettes, unchanged →
`980ec9a` his features ported → `f0e71af` orm intervals and fixes →
`91785a1` vignettes ported → `f4efe6d` vcov and check fixes.

### Where each part of the student's work is now

| p-ggeffects (Lukas Bolsinger) | termeffects | Status |
|:--|:--|:--|
| `parameterize()`, `parameterized_model()`, `predict()`, `original_model()`, `parameterizations()`, `parameter_variables()`, `parameterize_data()`, `is.parameterized_model()`, composition, `$` access | same functions (`R/parameterize.R`) | Copied; only insight helpers replaced. A misspelt name is now an error (he listed the silent drop as a limitation). |
| Parameterization history printer | `print(pmodel, show_history = TRUE)` (`R/print_parameterization_history.R`) | Copied unchanged (box characters as `\u` escapes for CRAN). |
| Range search for new variables (`.get_excribed_box`) | `term_grid()` / `term_predictions()` on a parameterized model | Copied; called per parameterization, with data variables held at their grid values. |
| `ggpredict()` on parameterized models | `term_predictions(pmodel, ...)` | Rebuilt on the new engine. |
| `ggeffects_layer` (curves through a heatmap) | `term_layer()`, `term_heatmap(layers = )`, `term_curves(layer)` | Rebuilt. His plot showed only the first layer and failed with several; both work now. Catches a parameterization that uses a replaced variable as input. |
| Simultaneous bands, `ci_method = "coefficient"` | `ci_method = "analytic"` | Rebuilt; correct now for GLMs and Cox models (link scale). |
| `ci_method = "bootstrap"` | `ci_method = "bootstrap"` | Rebuilt; refits by `rms::bootcov()` or resampling rows, sup-t band on the correct scale. His version mixed scales. |
| `ci_scope = "faceted"` | `ci_scope = "faceted"`, plus `"curve"` | Implemented (his was forced to global). |
| `test_simultaneous_predictions()` | same name | Rebuilt; adds a p-value and works on the link scale. |
| Difference plot of a test (`use_diff`) | `plot(test, "difference")` | Rebuilt; difference heatmap for surfaces. |
| Heatmap plot type | `term_heatmap()` | Rebuilt (viridis instead of red-white-blue). |
| Standard-error heatmap (`show_ci`) | `term_heatmap(fill = "std.error")` or `"band_width"`; `term_figure()$uncertainty` | Separate plots instead of a patchwork stack. |
| Kernel-density veil (`show_kernel_density`, `threshold_prob`) | `term_heatmap(show_support = , support_prob = )` | His kernel and bandwidth; density of the varying predictors only, and the usual highest-density-region threshold. |
| Direction arrows (`show_direction`) | `term_heatmap(show_direction = )`, `term_curves(show_direction = )` | Rebuilt. |
| `vcov()` of predictions | `vcov(term_predictions(...))` | Rebuilt (analytic intervals). |
| His tests for parameterized models | `tests/testthat/test-parameterize.R` | Three copied, three adapted. |
| His four vignettes | same vignettes, ported | See below and the track-changes page. |

Not carried over: the upstream ggeffects machinery (`ggemmeans()`, `ggaverage()`,
`ggeffect()`, `test_predictions()` pairwise contrasts, `pool_predictions()`,
`print_md()`/`print_html()`). These were ggeffects features rather than his,
and his layer versions of them were stubs or failed. Also not carried over:
`inherit_scb` (a layer reusing the heatmap's critical value). A layer gets its
own band, simultaneous over its curves, which is the right band for judging
them.

### Vignettes, with track changes
- Commit `2de18bc` holds the student's four vignettes exactly as he wrote
  them (renamed only). The commits after it change them, so
  `git diff 2de18bc -- vignettes/` shows every change.
- `validation/output/vignette_track_changes.html` shows the same changes as
  a redline: deletions struck through, insertions underlined. Regenerate it
  with `python validation/09_vignette_track_changes.py`.
- Besides the API changes, these corrections:
  - **eGFR vignette:** it said it used the 2021 CKD-EPI equation but used a
    made-up linear formula. It now uses the real equation and its inverse.
  - **eGFR vignette, Figure 3:** it said "LVEF versus eGFR" but plotted
    creatinine. It now plots eGFR.
  - **eGFR vignette, Figure 2:** the constant-eGFR overlay used hand-drawn
    contours because his layer API "was not confirmed". It is now a
    `term_layer()`, with its own curves plot.
  - **Parameterized-models vignette:** the `pmodel` typo (2 errors) is fixed,
    and weight now has mean 70 instead of 7000, an apparent typo given the
    first example.
  - **Limitations sections:** items that no longer hold (forced global scope,
    silently ignored names, bootstrap as a hidden fallback) are updated.

### orm intervals: is the bootstrap needed?
Brigitte's question. Simulations with a known true mean
(`validation/08_orm_interval_methods.R`, `validation/10_orm_intervals_complex_model.R`),
95% simultaneous bands, 300 data sets each:

| Setting | Default (linear predictor) | Bootstrap, 1000 refits | Delta, all intercepts |
|:--|:--|:--|:--|
| 10-level outcome, n = 200, small model | 0.943 | 0.927 | 0.877 |
| Continuous outcome, n = 300, small model | 0.950 | 0.943 | 0.930 |
| 10-level outcome, n = 300, 4x4-knot spline interaction + 6 covariates (about 24 coefficients) | 0.937 (width 24.7) | 0.933 (width 32.1) | 0.633 (width 22.6) |

Monte Carlo error is about +/- 0.025. An earlier run with only 200 bootstrap
refits gave bootstrap coverage of 0.93 / 0.91.

- **Default (`ci_method = "analytic"`):** intervals for the linear predictor,
  transformed by `Mean()`. This is what `rms::Predict(fun = Mean(fit))` does.
- **`ci_method = "delta"`:** the delta method on the mean scale with all
  intercepts, as `rms::Mean(fit)(X = , conf.int = )`. Available for comparison.
- **`ci_method = "bootstrap"`:** `rms::bootcov()` refits.

**Conclusion:** the bootstrap is not needed. The default reaches the stated
coverage, also in a model as complex as the LVEF example. The bootstrap
reaches it too, but with bands up to 30% wider. The delta method with all
intercepts is too narrow and is documented as not recommended. (For a
while today the delta method was the default; the simulation reversed that.)

### Figures with their covariates and caption
- `term_figure()` returns `$heatmap`, `$curves`, `$uncertainty`, `$covariates`
  (each held covariate, its value and why), `$caption` (a figure legend in
  words: what is predicted, the model, how the focal variables enter it, the
  held values, the intervals) and `$model`.
- With more than four held covariates, plot subtitles point to the caption
  instead of listing them.

### Other changes
- **Held values:** numeric covariates are held at their **median** (was the
  mean), categorical ones at their **most common** value (was the first
  level), as Brigitte asked.
- **Data:** `heartfailure` is bundled (CC BY 4.0, with attribution), so the
  vignettes do not need dlookr.
- **Model data:** models fitted without `data =` now work.
- **orm numerics:** when rms cannot invert a badly scaled model, the error
  says to rescale.
- **How-to vignette:** rewritten with an orm model of ejection fraction (ten
  covariates), `term_figure()`, and the interval comparison.

### Tests
- About 160 tests, all passing.
- `R CMD check`: 0 errors. The remaining warning (qpdf) and note (timestamp)
  come from the Docker container.

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

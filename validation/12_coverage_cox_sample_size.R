# Follow-up to 11_coverage_study.R: the analytic band for Cox models covered
# less than 95% at n = 300. Is that the small-sample behaviour of the Cox
# estimates (bias that shrinks with n), or an error in termeffects (which
# would not shrink)? Same data-generating mechanism, surface estimand,
# n = 300, 1000, 4000; 1000 runs each. Results: vignettes/coverage_study_cox.csv.
# Run with OPENBLAS_NUM_THREADS=1; about 15 minutes on 7 cores.
src <- readLines("validation/11_coverage_study.R")
eval(parse(text = src[1:(grep("^# ---- scenarios", src) - 1)]))
set.seed(20261106)
one <- function(n) {
  d <- outcomes$cox(population(n))
  fit <- fit_model("cox", d)
  p <- term_predictions(fit, c("x1", "x2"), values = grids$surface, scale = "link", nsim = 5000)
  z <- (p$estimate - true_values("cox", fit, p)) / p$std.error
  c(n = n, events = sum(d$status), mean_z = mean(z), band = all(abs(z) <= attr(p, "critical_value")),
    pointwise_mean = mean(abs(z) <= attr(p, "pointwise_critical_value")))
}
res <- do.call(rbind, lapply(c(300, 1000, 4000), function(n) {
  r <- do.call(rbind, mclapply(seq_len(runs_analytic), function(i) one(n), mc.cores = cores))
  data.frame(n = n, runs = nrow(r), events = mean(r[, "events"]),
             mean_standardized_error = mean(r[, "mean_z"]),
             mcse_mean_standardized_error = sd(r[, "mean_z"]) / sqrt(nrow(r)),
             coverage_band = mean(r[, "band"]), mcse_band = sqrt(mean(r[, "band"]) * (1 - mean(r[, "band"])) / nrow(r)),
             coverage_pointwise_mean = mean(r[, "pointwise_mean"]))
}))
print(res, digits = 3)
if (!quick) write.csv(res, "vignettes/coverage_study_cox.csv", row.names = FALSE)

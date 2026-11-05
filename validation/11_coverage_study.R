# Coverage of termeffects' simultaneous confidence bands: a simulation study.
#
# Described for publication in vignettes/coverage_study.Rmd, which reads the
# results written here (vignettes/coverage_study_results.csv). Run from the
# package folder with one BLAS thread per process (OPENBLAS_NUM_THREADS=1),
# otherwise the parallel runs slow each other down; about 2 hours on 7 cores:
#   Rscript validation/11_coverage_study.R            (full study)
#   Rscript validation/11_coverage_study.R quick      (5 runs per scenario, to test)
#
# Structure follows Morris, White & Crowther (2019), "Using simulation
# studies to evaluate statistical methods", Stat Med 38:2074-2102 (ADEMP).

pkgload::load_all(".", quiet = TRUE)
suppressMessages(library(rms))
library(survival)
library(parallel)

quick <- "quick" %in% commandArgs(trailingOnly = TRUE)

runs_analytic <- if (quick) 5 else 1000
runs_bootstrap <- if (quick) 5 else 500
nboot <- 500
cores <- 7
RNGkind("L'Ecuyer-CMRG")
set.seed(20261105)

# ---- data-generating mechanism ---------------------------------------------
# x1, x2: the two focal variables (correlated, standard normal margins)
# w: an adjustment covariate, correlated with x2
population <- function(n) {
  x1 <- rnorm(n)
  x2 <- 0.4 * x1 + sqrt(1 - 0.4^2) * rnorm(n)
  w <- 0.3 * x2 + sqrt(1 - 0.3^2) * rnorm(n)
  data.frame(x1 = x1, x2 = x2, w = w)
}

# restricted cubic splines with 4 knots at the 5, 35, 65, 95% quantiles of a
# standard normal (Harrell's default positions), fixed in advance
knots <- qnorm(c(0.05, 0.35, 0.65, 0.95))
k <- paste(deparse(knots, control = "digits17"), collapse = "")
rhs <- sprintf("rcs(x1, %s) * rcs(x2, %s)", k, k)

# The true surface f(x1, x2) is the best approximation, within this spline
# model, of a smooth non-linear interaction g. The fitted model is therefore
# correctly specified, and coverage measures the inference, not model error.
g <- function(x1, x2) 0.8 * tanh(1.5 * x1) + 0.4 * x2 + 0.5 * tanh(x1) * x2 - 0.2 * x2^2
pop <- population(2e5)
pop$g <- g(pop$x1, pop$x2)
truth_fit <- lm(as.formula(paste("g ~", rhs)), data = pop)
f_true <- function(d) unname(predict(truth_fit, d))
beta_w <- 0.5

outcomes <- list(
  linear = function(d) { d$y <- f_true(d) + beta_w * d$w + rnorm(nrow(d)); d },
  logistic = function(d) { d$y <- rbinom(nrow(d), 1, plogis(-1 + f_true(d) + beta_w * d$w)); d },
  cox = function(d) {
    t_event <- rexp(nrow(d), 0.1 * exp(f_true(d) + beta_w * d$w))
    t_cens <- runif(nrow(d), 0, 15)
    d$time <- pmin(t_event, t_cens)
    d$status <- as.numeric(t_event <= t_cens)
    d
  },
  orm = function(d) { d$y <- f_true(d) + beta_w * d$w + rlogis(nrow(d), scale = 0.5); d },
  mixed = function(d) {
    d$cluster <- rep(seq_len(nrow(d) / 10), each = 10)
    u <- rnorm(nrow(d) / 10, sd = 0.8)
    d$y <- f_true(d) + beta_w * d$w + u[d$cluster] + rnorm(nrow(d))
    d
  }
)

fit_model <- function(outcome, d) {
  form <- function(lhs) as.formula(paste(lhs, "~", rhs, "+ w"))
  switch(outcome,
    linear = lm(form("y"), data = d),
    logistic = glm(form("y"), family = binomial, data = d),
    cox = coxph(form("Surv(time, status)"), data = d),
    orm = orm(form("y"), data = d, x = TRUE, y = TRUE),
    mixed = suppressMessages(lme4::lmer(as.formula(paste("y ~", rhs, "+ w + (1 | cluster)")), data = d))
  )
}

# ---- estimands ---------------------------------------------------------------
# A: the surface over x1, x2 in their central 90% (20 x 20 grid)
# B: three curves along x1, at x2 = -1, 0, 1 (30 points each)
q <- qnorm(c(0.05, 0.95))
grids <- list(
  surface = list(x1 = seq(q[1], q[2], length.out = 20), x2 = seq(q[1], q[2], length.out = 20)),
  curves = list(x1 = seq(q[1], q[2], length.out = 30), x2 = c(-1, 0, 1))
)

# true values on the scale of the estimates (link scale; for orm the mean).
# Monotone transformations (inverse link, exp) do not change coverage.
true_values <- function(outcome, fit, p) {
  grid <- as.data.frame(p)
  if (outcome == "cox") {
    # hazard ratio vs. the sample's average covariate vector (as coxph centres)
    X <- termeffects:::.design.coxph(fit, grid)$X
    beta <- c(coef(truth_fit)[-1], w = beta_w)
    return(drop(X %*% beta[colnames(X)]))
  }
  f_true(grid) + beta_w * grid$w + if (outcome == "logistic") -1 else 0
}

one_run <- function(outcome, n, method) {
  d <- outcomes[[outcome]](population(n))
  fit <- fit_model(outcome, d)
  scale <- if (outcome == "orm") "response" else "link"
  out <- lapply(names(grids), function(estimand) {
    p <- term_predictions(fit, c("x1", "x2"), values = grids[[estimand]], scale = scale,
                          ci_method = method, nsim = 5000, nboot = nboot)
    truth <- true_values(outcome, fit, p)
    inside_band <- truth >= p$band.low & truth <= p$band.high
    inside_pw <- truth >= p$conf.low & truth <= p$conf.high
    data.frame(estimand = estimand, band = all(inside_band), pointwise_all = all(inside_pw),
               pointwise_mean = mean(inside_pw), critical_value = attr(p, "critical_value"),
               pointwise_critical_value = attr(p, "pointwise_critical_value"),
               events = if (outcome %in% c("logistic", "cox")) sum(d[[if (outcome == "cox") "status" else "y"]]) else NA)
  })
  do.call(rbind, out)
}

# ---- scenarios -----------------------------------------------------------------
scenarios <- rbind(
  expand.grid(outcome = names(outcomes), n = c(300, 1000), method = "analytic", stringsAsFactors = FALSE),
  expand.grid(outcome = c("linear", "logistic", "cox", "orm"), n = c(300, 1000), method = "bootstrap",
              stringsAsFactors = FALSE)
)

results <- list()
failures <- list()
for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i, ]
  runs <- if (s$method == "bootstrap") runs_bootstrap else runs_analytic
  started <- Sys.time()
  res <- mclapply(seq_len(runs), function(r) {
    tryCatch(cbind(run = r, one_run(s$outcome, s$n, s$method)), error = function(e) conditionMessage(e))
  }, mc.cores = cores)
  failed <- vapply(res, is.character, logical(1))
  failures[[i]] <- cbind(s, failed = sum(failed), first_error = if (any(failed)) res[[which(failed)[1]]] else NA)
  results[[i]] <- cbind(s, do.call(rbind, res[!failed]), row.names = NULL)
  message(sprintf("%-9s n = %4d %-9s %4d runs, %d failed, %.1f min", s$outcome, s$n, s$method, runs,
                  sum(failed), as.numeric(difftime(Sys.time(), started, units = "mins"))))
}
results <- do.call(rbind, results)
failures <- do.call(rbind, failures)

# per-run results (for checking) and the summary used by the vignette
dir.create("validation/output", showWarnings = FALSE)
write.csv(results, "validation/output/11_coverage_runs.csv", row.names = FALSE)
write.csv(failures, "validation/output/11_coverage_failures.csv", row.names = FALSE)

summary_table <- do.call(rbind, lapply(split(results, list(results$outcome, results$n, results$method, results$estimand),
                                             drop = TRUE), function(r) {
  m <- nrow(r)
  f <- failures[failures$outcome == r$outcome[1] & failures$n == r$n[1] & failures$method == r$method[1], ]
  data.frame(outcome = r$outcome[1], n = r$n[1], method = r$method[1], estimand = r$estimand[1], runs = m,
             failed = f$failed,
             coverage_band = mean(r$band), mcse_band = sqrt(mean(r$band) * (1 - mean(r$band)) / m),
             coverage_pointwise_all = mean(r$pointwise_all),
             coverage_pointwise_mean = mean(r$pointwise_mean),
             critical_value = mean(r$critical_value),
             pointwise_critical_value = mean(r$pointwise_critical_value),
             events = mean(r$events))
}))
summary_table <- summary_table[order(summary_table$method, summary_table$outcome, summary_table$estimand, summary_table$n), ]
if (!quick) {
  write.csv(summary_table, "vignettes/coverage_study_results.csv", row.names = FALSE)
}
print(summary_table, digits = 3, row.names = FALSE)
writeLines(capture.output(sessionInfo()), "validation/output/11_session_info.txt")

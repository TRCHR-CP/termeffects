# Compares the student's simultaneous critical value with independent calculations.
pkgload::load_all("/home/rstudio/data/Staff/Brigitte Mueller/R-Packages/p-ggeffects-archive", quiet = TRUE)  # the fork, loaded as 'ggeffects'
library(splines)

set.seed(2026)
n   <- 400
dat <- data.frame(x = runif(n), z = runif(n))
dat$y <- 1 + sin(2 * pi * dat$x) * (1 + dat$z) + rnorm(n, sd = 0.5)
fit <- lm(y ~ ns(x, df = 4) * ns(z, df = 3), data = dat)

x_vals <- seq(0, 1, by = 0.05)
z_vals <- c(0.2, 0.5, 0.8)

# 1. Student's implementation (run twice: it sets no seed, so expect small differences)
pred <- ggpredict(fit, terms = c("x [0:1 by=0.05]", "z [0.2, 0.5, 0.8]"),
                  ci_method = "coefficient")
crit_student   <- unname(attr(pred, "critical_values"))
crit_student_2 <- unname(attr(ggpredict(fit, terms = c("x [0:1 by=0.05]", "z [0.2, 0.5, 0.8]"),
                                        ci_method = "coefficient"), "critical_values"))

crit_student_boot <- unname(attr(ggpredict(fit, terms = c("x [0:1 by=0.05]", "z [0.2, 0.5, 0.8]"),
                                           ci_method = "bootstrap"), "critical_values"))

# Rebuild the same grid independently and check predictions and SEs agree
grid  <- data.frame(x = pred$x, z = as.numeric(as.character(pred$group)))
X     <- model.matrix(delete.response(terms(fit)), grid)
stopifnot(isTRUE(all.equal(drop(X %*% coef(fit)), pred$predicted,
                           check.attributes = FALSE, tolerance = 1e-6)))
Sigma <- X %*% vcov(fit) %*% t(X)
se    <- sqrt(diag(Sigma))
stopifnot(isTRUE(all.equal(se, pred$std.error, check.attributes = FALSE, tolerance = 1e-6)))
R <- cov2cor(Sigma)

# 2. Hand-coded max-|t| simulation (normal and t reference)
sim_maxt <- function(R, nsim = 50000, df = Inf) {
  e    <- eigen(R, symmetric = TRUE)
  keep <- e$values > 1e-10
  L    <- e$vectors[, keep, drop = FALSE] %*% diag(sqrt(e$values[keep]), sum(keep))
  Z    <- matrix(rnorm(nsim * sum(keep)), nsim) %*% t(L)
  if (is.finite(df)) Z <- Z / sqrt(rchisq(nsim, df) / df)
  apply(abs(Z), 1, max)
}
crit_sim_normal <- unname(quantile(sim_maxt(R), 0.95))
crit_sim_t      <- unname(quantile(sim_maxt(R, df = df.residual(fit)), 0.95))

# 3. mvtnorm (the method modelbased's "sup-t" uses) — may fail on a singular matrix
crit_qmvt <- tryCatch(
  mvtnorm::qmvt(0.95, tail = "both.tails", corr = R, df = df.residual(fit))$quantile,
  error = function(e) { message("qmvt: ", conditionMessage(e)); NA })

# 4. modelbased p_adjust = "sup-t"           (<- syntax may need adjusting)
by_spec <- c(paste0("x = c(", paste(x_vals, collapse = ", "), ")"),
             paste0("z = c(", paste(z_vals, collapse = ", "), ")"))
crit_modelbased <- tryCatch({
  mb <- modelbased::estimate_means(fit, by = by_spec, p_adjust = "sup-t")
  median((mb$CI_high - mb$CI_low) / (2 * mb$SE))
}, error = function(e) { message("modelbased: ", conditionMessage(e)); NA })

# 5. marginaleffects + multcomp single-step  (<- may need adjusting)
crit_marginaleffects <- tryCatch({
  p <- marginaleffects::predictions(fit, newdata = marginaleffects::datagrid(x = x_vals, z = z_vals, model = fit))
  h <- marginaleffects::hypotheses(p, multcomp = "single-step")
  median((h$conf.high - h$conf.low) / (2 * h$std.error))
}, error = function(e) { message("marginaleffects: ", conditionMessage(e)); NA })

out <- data.frame(
  method = c("pointwise (reference)", "student (run 1)", "student (run 2)","student bootstrap",
             "hand simulation, normal", "hand simulation, t",
             "mvtnorm::qmvt", "modelbased sup-t", "marginaleffects single-step"),
  critical_value = round(c(qnorm(0.975), crit_student, crit_student_2,crit_student_boot,
                           crit_sim_normal, crit_sim_t, crit_qmvt,
                           crit_modelbased, crit_marginaleffects), 3)
)
print(out)
write.csv(out, "validation/output/02_critical_values.csv", row.names = FALSE)


#The student’s value should be within about ±0.03 of “hand simulation, normal”. The t-based methods (t simulation, qmvt, modelbased) come out slightly larger, which is expected.
#If modelbased shows about 1.96, it ignored p_adjust, meaning sup-t isn’t supported for predictions. That’s worth knowing too.
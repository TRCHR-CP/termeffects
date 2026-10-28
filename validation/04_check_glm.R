pkgload::load_all("/home/rstudio/data/Staff/Brigitte Mueller/R-Packages/p-ggeffects-archive", quiet = TRUE)
set.seed(3)
dat   <- data.frame(x = runif(600), z = runif(600))
dat$y <- rbinom(600, 1, plogis(-1 + dat$x + 2 * dat$x * dat$z))
fit   <- glm(y ~ x * z, family = binomial, data = dat)

pred <- ggpredict(fit, terms = c("x [0:1 by=0.05]", "z [0.2, 0.5, 0.8]"), ci_method = "coefficient")
crit <- unname(attr(pred, "critical_values"))

grid <- data.frame(x = pred$x, z = as.numeric(as.character(pred$group)))
eta  <- predict(fit, grid, se.fit = TRUE)                 # link scale
correct_low  <- plogis(eta$fit - crit * eta$se.fit)
correct_high <- plogis(eta$fit + crit * eta$se.fit)

cmp <- data.frame(predicted = pred$predicted,
                  student_low = pred$simconf.low, correct_low,
                  student_high = pred$simconf.high, correct_high)
print(head(cmp, 10))
cat("Student matches correct band:", isTRUE(all.equal(pred$simconf.low, correct_low, check.attributes = FALSE)), "\n")
cat("Student band inside [0,1]:", all(pred$simconf.low >= 0 & pred$simconf.high <= 1), "\n")
cat("Simultaneous band contains pointwise band:",
    all(pred$simconf.low <= pred$conf.low & pred$simconf.high >= pred$conf.high), "\n")
write.csv(cmp, "validation/output/04_glm.csv", row.names = FALSE)

# Bootstrap path: a 95% simultaneous critical value must be larger than the
# pointwise 1.96, and should be close to the coefficient-based value above.
pred_boot <- ggpredict(fit, terms = c("x [0:1 by=0.05]", "z [0.2, 0.5, 0.8]"),
                       ci_method = "bootstrap", iterations = 200, verbose = FALSE)
crit_boot <- unname(attr(pred_boot, "critical_values"))
cat("Critical value, coefficient:", round(crit, 3),
    " bootstrap:", round(crit_boot, 3), " (must be > 1.96)\n")

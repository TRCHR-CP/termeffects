pkgload::load_all("/home/rstudio/data/Staff/Brigitte Mueller/R-Packages/p-ggeffects-archive", quiet = TRUE)
set.seed(1)

truth_fun <- function(x, z) 1 + 2 * x - 3 * x^2 + 2 * x * z + z^2   # inside the model space, so no bias
nsim <- 200
covered_scb <- covered_pw <- logical(nsim)

for (i in seq_len(nsim)) {
  dat   <- data.frame(x = runif(300), z = runif(300))
  dat$y <- truth_fun(dat$x, dat$z) + rnorm(300)
  fit   <- lm(y ~ poly(x, 2, raw = TRUE) * poly(z, 2, raw = TRUE), data = dat)
  pred  <- suppressMessages(ggpredict(fit, terms = c("x [0:1 by=0.05]", "z [0.2, 0.5, 0.8]"),
                                      ci_method = "coefficient", verbose = FALSE))
  truth <- truth_fun(pred$x, as.numeric(as.character(pred$group)))
  covered_scb[i] <- all(truth >= pred$simconf.low & truth <= pred$simconf.high)
  covered_pw[i]  <- all(truth >= pred$conf.low    & truth <= pred$conf.high)
  if (i %% 20 == 0) cat(i, "")
}

res <- data.frame(band = c("simultaneous", "pointwise"),
                  coverage_of_whole_surface = c(mean(covered_scb), mean(covered_pw)))
print(res)
write.csv(res, "validation/output/03_coverage.csv", row.names = FALSE)

#multaneous about 0.95 (0.92–0.98 is fine with 200 runs); pointwise clearly lower. If simultaneous is well below 0.92, the bands are too narrow.
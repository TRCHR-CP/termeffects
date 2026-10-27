pkgload::load_all("/home/rstudio/data/Staff/Brigitte Mueller/R-Packages/p-ggeffects-archive", quiet = TRUE)
library(splines); library(survival)
lung2 <- na.omit(lung[c("time", "status", "age", "ph.karno")])
fit_cox <- coxph(Surv(time, status) ~ ns(age, 3) * ns(ph.karno, 2), data = lung2)
r <- try(ggpredict(fit_cox, terms = c("age", "ph.karno [60, 80, 100]"), ci_method = "coefficient"))
cat("Cox:", if (inherits(r, "try-error")) "FAILS" else "works", "\n")
# add: lme4::lmer(...), glmmTMB, mgcv::gam(...), rms::cph(...) — whatever TRCP uses
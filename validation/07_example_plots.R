# Example plots for visual checking: an orm model for a 0-100 outcome.
pkgload::load_all(".", quiet = TRUE)
library(rms)
set.seed(7)
n <- 400
d <- data.frame(age = runif(n, 30, 85), bmi = runif(n, 18, 40),
                sex = factor(sample(c("female", "male"), n, TRUE)))
d$score <- round(pmin(100, pmax(0, 60 + 25 * sin((d$age - 30) / 18) * (d$bmi - 29) / 11 +
                                 3 * (d$sex == "male") + rnorm(n, sd = 12))))
fit <- orm(score ~ rcs(age, 4) * rcs(bmi, 4) + sex, data = d)

dir.create("validation/output/plots", showWarnings = FALSE)
surface <- term_predictions(fit, c("age", "bmi"), n = 60)
cut_values <- list(bmi = c(22, 30, 38))
cuts <- term_predictions(fit, c("age", "bmi"), values = cut_values, n = 60)
print(cuts, n = 3)

ggplot2::ggsave("validation/output/plots/orm_heatmap.png",
                term_heatmap(surface, cuts = cut_values, show_data = TRUE), width = 7, height = 5, dpi = 110)
ggplot2::ggsave("validation/output/plots/orm_cuts.png",
                term_curves(cuts, band = "both"), width = 7, height = 5, dpi = 110)

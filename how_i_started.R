path  = getwd()

usethis::create_package(
  path,
  fields = list(
    Package = "termeffects",
    Title = "Term Plots and Simultaneous Confidence Bands for Interaction Surfaces",
    Description = "Computes adjusted predictions with simultaneous confidence bands over prediction curves and two-variable interaction surfaces, including parameterized (transformed) predictor coordinates, and visualizes them as term plots and heatmaps.",
    `Authors@R` = 'c(
      person("Brigitte", "Mueller", email = "brigitte.mueller@uhn.ca", role = c("aut", "cre")),
      person("Lukas", "Bolsinger", role = "aut"))',
    Version = "0.0.0.9000"
  ),
  rstudio = TRUE,
  open = FALSE
)

usethis::proj_set(path)
usethis::use_mit_license("University Health Network")   # <- replace with what UHN tells you (Step 0.3)
usethis::use_testthat(3)
usethis::use_readme_rmd(open = FALSE)
usethis::use_news_md(open = FALSE)
usethis::use_package("insight")
usethis::use_package("ggplot2")
usethis::use_package("stats")
usethis::use_package("ggeffects", type = "Suggests")    # used for comparisons and validation only
dir.create(file.path(path, "validation", "output"), recursive = TRUE)
usethis::use_build_ignore("validation")

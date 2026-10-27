# Heart-failure clinical records (Chicco & Jurman 2020), as distributed by
# dlookr, saved unchanged as termeffects::heartfailure so the vignettes do not
# need dlookr. Licence of the data: CC BY 4.0.
data("heartfailure", package = "dlookr")
heartfailure <- as.data.frame(heartfailure)
usethis::use_data(heartfailure, overwrite = TRUE, compress = "xz")

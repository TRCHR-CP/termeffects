#' Heart-failure clinical records
#'
#' Clinical records of 299 patients with heart failure (Faisalabad, Pakistan,
#' 2015), as distributed by the dlookr package. Used in the vignettes, with
#' the ejection fraction as the outcome.
#'
#' @format A data frame with 299 rows and 13 variables:
#' \describe{
#'   \item{age}{age (years)}
#'   \item{anaemia}{decrease of red blood cells or haemoglobin (No, Yes)}
#'   \item{cpk_enzyme}{creatine phosphokinase in the blood (mcg/L)}
#'   \item{diabetes}{diabetes (No, Yes)}
#'   \item{ejection_fraction}{left ventricular ejection fraction (percent)}
#'   \item{hblood_pressure}{hypertension (No, Yes)}
#'   \item{platelets}{platelets in the blood (per microlitre)}
#'   \item{creatinine}{serum creatinine (mg/dL)}
#'   \item{sodium}{serum sodium (mEq/L)}
#'   \item{sex}{Female, Male}
#'   \item{smoking}{smoker (No, Yes)}
#'   \item{time}{follow-up (days)}
#'   \item{death_event}{died during follow-up (No, Yes)}
#' }
#' @source Chicco D, Jurman G (2020). Machine learning can predict survival of
#'   patients with heart failure from serum creatinine and ejection fraction
#'   alone. *BMC Medical Informatics and Decision Making* 20, 16.
#'   \doi{10.1186/s12911-020-1023-5}. Data licensed under CC BY 4.0, via
#'   Kaggle (<https://www.kaggle.com/andrewmvd/heart-failure-clinical-data>)
#'   and the dlookr package.
"heartfailure"

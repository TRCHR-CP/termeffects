# Renders the student's vignettes against the fork's source code, without
# writing anything into the archive folder.
#
# - The fork is loaded with pkgload::load_all(), so the vignettes use the
#   student's code and not the CRAN version of ggeffects.
# - Each vignette is copied to a temporary folder and rendered there in a
#   fresh R session (callr), so one vignette cannot affect another.
# - The finished HTML files go to validation/output/vignettes/.
#
# Run from the termeffects folder:  source("validation/00_render_student_vignettes.R")

archive <- "/home/rstudio/data/Staff/Brigitte Mueller/R-Packages/p-ggeffects-archive"
out_dir <- normalizePath(file.path("validation", "output", "vignettes"), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

vignettes <- c(
  "introduction_p_ggeffects.Rmd",
  "introduction_parameterized_model.Rmd",
  "p_ggeffects_heart_failure.Rmd",
  "p_ggeffects_demonstration.Rmd"
)

render_one <- function(archive, vignette, out_dir) {
  work_dir <- tempfile("vignette_")
  dir.create(work_dir)
  file.copy(file.path(archive, "vignettes", vignette), work_dir)

  pkgload::load_all(archive, quiet = TRUE)
  set.seed(2026)

  # error = TRUE: a failing chunk shows its error in the HTML and rendering
  # continues, so broken chunks are visible instead of stopping the vignette
  format <- rmarkdown::html_document(
    toc = TRUE, toc_float = TRUE, self_contained = TRUE
  )
  format$knitr$opts_chunk$error <- TRUE

  rmarkdown::render(
    file.path(work_dir, vignette),
    output_format = format,
    output_dir = out_dir,
    intermediates_dir = work_dir,
    envir = globalenv(),
    quiet = TRUE
  )
}

status <- vapply(vignettes, function(v) {
  cat("Rendering", v, "... ")
  res <- tryCatch(
    {
      callr::r(render_one, args = list(archive, v, out_dir), show = FALSE)
      "ok"
    },
    error = function(e) paste("FAILED:", conditionMessage(e))
  )
  cat(res, "\n")
  res
}, character(1))

print(data.frame(vignette = names(status), status = unname(status)), right = FALSE)

# Render every vignette (and validation/feature_tour.Rmd) to rendered_vignettes/ as GitHub Markdown with PNG
# figures, readable on GitHub, in VS Code's Markdown preview and in RStudio's
# visual editor. Run from the package folder:
#   Rscript validation/render_vignettes.R              (all vignettes)
#   Rscript validation/render_vignettes.R model_types  (only some)

devtools::load_all(quiet = TRUE)
out <- "rendered_vignettes"
dir.create(out, showWarnings = FALSE)

# the package vignettes, and the feature tour (kept in validation/, not shipped)
rmds <- c(list.files("vignettes", pattern = "\\.Rmd$", full.names = TRUE),
          "validation/feature_tour.Rmd")
only <- commandArgs(trailingOnly = TRUE)
if (length(only)) {
  rmds <- rmds[tools::file_path_sans_ext(basename(rmds)) %in% only]
}

for (rmd in rmds) {
  name <- tools::file_path_sans_ext(basename(rmd))
  title <- rmarkdown::yaml_front_matter(rmd)$title
  unlink(file.path(out, paste0(name, "_files")), recursive = TRUE)  # no stale figures

  md <- rmarkdown::render(rmd, output_dir = out, quiet = TRUE, envir = new.env(),
                          output_format = rmarkdown::md_document(variant = "gfm"))
  # Markdown output drops the title; put it back as a heading. Figure links
  # (Markdown and <img src=...>) come out as absolute paths; make them
  # relative so they work anywhere.
  lines <- gsub("(\\]\\(|src=\")[^)\"]*/rendered_vignettes/", "\\1", readLines(md))
  writeLines(c(paste("#", title), "", lines), md)
  message("rendered ", name)
}

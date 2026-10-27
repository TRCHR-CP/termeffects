# Runs the ggeffects test suite on the student's fork AND on the upstream
# ggeffects commit his work started from, then lists the tests that fail only
# in the fork (i.e. failures the fork introduced).
#
# Both versions are exported with `git archive` into a temporary folder, so the
# archive folder is never written to (devtools::test() writes snapshot files
# into the package it tests).

archive  <- "/home/rstudio/data/Staff/Brigitte Mueller/R-Packages/p-ggeffects-archive"
upstream <- "392198e1^"  # parent of the student's first commit
work     <- tempfile("student_tests_")
dir.create(work)

export_version <- function(rev, dest) {
  dir.create(dest)
  cmd <- sprintf("cd %s && git archive %s | tar -x -C %s",
                 shQuote(archive), shQuote(rev), shQuote(dest))
  stopifnot(system(cmd) == 0)
  dest
}

run_tests <- function(path) {
  callr::r(function(path) {
    testthat::set_max_fails(Inf)
    res <- devtools::test(path, stop_on_failure = FALSE,
                          reporter = testthat::SilentReporter$new())
    df <- as.data.frame(res)
    data.frame(file = df$file, test = df$test, n = df$nb, failed = df$failed,
               skipped = df$skipped, error = df$error, warnings = df$warning)
  }, args = list(path))
}

fork_res     <- run_tests(export_version("HEAD", file.path(work, "fork")))
upstream_res <- run_tests(export_version(upstream, file.path(work, "upstream")))

is_bad <- function(d) d$failed > 0 | d$error
key    <- function(d) paste(d$file, d$test, sep = " :: ")
new_failures <- setdiff(key(fork_res)[is_bad(fork_res)], key(upstream_res)[is_bad(upstream_res)])

write.csv(fork_res, "validation/output/01_tests_fork.csv", row.names = FALSE)
write.csv(upstream_res, "validation/output/01_tests_upstream.csv", row.names = FALSE)
write.csv(data.frame(failing_in_fork_only = new_failures),
          "validation/output/01_tests_new_failures.csv", row.names = FALSE)

cat("Fork:     ", nrow(fork_res), "test blocks,", sum(is_bad(fork_res)), "failing\n")
cat("Upstream: ", nrow(upstream_res), "test blocks,", sum(is_bad(upstream_res)), "failing\n")
cat("Failing in fork only:", length(new_failures), "\n")
print(sort(table(sub(" :: .*", "", new_failures)), decreasing = TRUE))

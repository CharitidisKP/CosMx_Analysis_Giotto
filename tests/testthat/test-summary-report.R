# Phase 13 — Render_Merged_Summary.R + Rmd template. Tests the skip-and-log
# fallback when rmarkdown / pandoc isn't available, and that the helper
# parses without errors.

test_that("Render_Merged_Summary.R parses and exposes render_merged_summary", {
  env <- .source_helper(file.path("Helper_Scripts", "Render_Merged_Summary.R"))
  expect_true(exists("render_merged_summary", envir = env, inherits = FALSE))
})

test_that("render_merged_summary skips when rmarkdown missing", {
  skip_if(requireNamespace("rmarkdown", quietly = TRUE),
          "rmarkdown installed; this test only verifies the skip path")
  env <- .source_helper(file.path("Helper_Scripts", "Render_Merged_Summary.R"))
  res <- expect_invisible(env$render_merged_summary(
    merged_output_dir = tempdir()))
  expect_null(res)
})

test_that("Rmd template front matter parses as YAML", {
  template <- file.path(.repo_dir(), "Helper_Scripts", "templates",
                        "merged_run_summary.Rmd")
  skip_if_not(file.exists(template), "Rmd template not found")
  lines <- readLines(template, warn = FALSE)
  delims <- which(lines == "---")
  expect_gte(length(delims), 2L)
  yaml_text <- paste(lines[(delims[1] + 1L):(delims[2] - 1L)], collapse = "\n")
  parsed <- yaml::yaml.load(yaml_text)
  expect_true("title" %in% names(parsed))
  expect_true("output" %in% names(parsed))
  expect_true("params" %in% names(parsed))
  expect_true("merged_output_dir" %in% names(parsed$params))
  expect_true("run_label" %in% names(parsed$params))
})

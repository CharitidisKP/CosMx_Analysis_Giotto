#!/usr/bin/env Rscript
# Render Helper_Scripts/extract_results_summary.Rmd to HTML.
# Usage: Rscript Helper_Scripts/Render_Results_Summary.R [project_dir]

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || identical(a, "")) b else a
}

render_results_summary <- function(project_dir = NULL,
                                   rmd        = NULL,
                                   output_dir = NULL,
                                   output_file = "extract_results_summary.html",
                                   extra_params = list()) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    cat("Warning: rmarkdown not installed; skipping results summary render.\n")
    return(invisible(NULL))
  }
  if (!isTRUE(rmarkdown::pandoc_available())) {
    cat("Warning: pandoc not available; cannot render results summary.\n")
    return(invisible(NULL))
  }

  if (is.null(project_dir) || !nzchar(project_dir)) {
    env_pd <- Sys.getenv("COSMX_PROJECT_DIR", unset = "")
    project_dir <- if (nzchar(env_pd)) env_pd
                   else normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."),
                                      mustWork = FALSE)
  }
  project_dir <- normalizePath(project_dir, mustWork = TRUE)

  if (is.null(rmd) || !nzchar(rmd)) {
    rmd <- file.path(project_dir, "Helper_Scripts", "extract_results_summary.Rmd")
  }
  if (!file.exists(rmd)) {
    stop("Rmd template not found: ", rmd)
  }

  if (is.null(output_dir) || !nzchar(output_dir)) {
    cfg_path <- file.path(project_dir, "Parameters", "config.yaml")
    cfg <- if (file.exists(cfg_path) && requireNamespace("yaml", quietly = TRUE)) {
      yaml::read_yaml(cfg_path)
    } else list()
    out_root <- if (!is.null(cfg$paths$output_dir) && nzchar(cfg$paths$output_dir)) {
      normalizePath(file.path(dirname(cfg_path), cfg$paths$output_dir), mustWork = FALSE)
    } else file.path(project_dir, "Output")
    output_dir <- file.path(out_root, "results_summary")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  params <- modifyList(list(project_dir = project_dir), extra_params)

  cat("Rendering results summary...\n")
  cat("  Rmd     : ", rmd, "\n", sep = "")
  cat("  Output  : ", file.path(output_dir, output_file), "\n", sep = "")

  out_path <- tryCatch({
    rmarkdown::render(
      input         = rmd,
      output_file   = output_file,
      output_dir    = output_dir,
      params        = params,
      envir         = new.env(),
      quiet         = TRUE
    )
    file.path(output_dir, output_file)
  }, error = function(e) {
    cat("Warning: results summary render failed (non-fatal): ", conditionMessage(e), "\n", sep = "")
    NULL
  })

  csvs <- list.files(output_dir, pattern = "\\.csv$", full.names = FALSE)
  cat("CSVs in ", output_dir, ": ", length(csvs), "\n", sep = "")
  invisible(out_path)
}

if (!interactive() && identical(sys.nframe(), 0L)) {
  args <- commandArgs(trailingOnly = TRUE)
  pd <- if (length(args) >= 1L) args[[1L]] else NULL
  render_results_summary(project_dir = pd)
}

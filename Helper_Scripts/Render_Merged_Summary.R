#!/usr/bin/env Rscript
# ==============================================================================
# Render_Merged_Summary.R
# Renders the merged-run summary HTML report from an R Markdown template.
#
# Stage 1, Phase 13: knits diagnostics, composition, niche frequencies,
# top DE, paired LRs, paired pathways, and the smiDE blocking audit into
# `<merged_output_dir>/merged_<label>_summary.html`.
#
# Uses R Markdown (not Quarto) so the existing R image works without
# adding `quarto-cli`. pandoc must be available — Install_Stage1_Dependencies.R
# asserts this at install time.
# ==============================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}

#' Render the merged-run summary HTML.
#'
#' @param merged_output_dir Path to the merged-pipeline output directory
#'   (the parent of 10_Merged/, 14_BANKSY/, 15_Composition/, etc.).
#' @param run_label Human-readable run label (used in the report header
#'   and the output filename). Defaults to basename(merged_output_dir).
#' @param template Optional override path to an Rmd template. Defaults to
#'   Helper_Scripts/templates/merged_run_summary.Rmd.
render_merged_summary <- function(merged_output_dir,
                                  run_label = NULL,
                                  template  = NULL) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    cat("⚠ rmarkdown not installed — skipping merged-run summary report.\n")
    return(invisible(NULL))
  }
  if (!isTRUE(rmarkdown::pandoc_available())) {
    cat("⚠ pandoc not available — cannot render summary report.\n")
    return(invisible(NULL))
  }
  helpers_dir <- dirname(normalizePath(sys.frame(1)$ofile %||%
                                        getOption("script.dir") %||% ".",
                                      mustWork = FALSE))
  if (is.null(template) || !nzchar(template)) {
    template <- file.path(helpers_dir, "templates", "merged_run_summary.Rmd")
  }
  if (!file.exists(template)) {
    cat("⚠ Summary template not found at: ", template, "\n", sep = "")
    return(invisible(NULL))
  }
  if (is.null(run_label) || !nzchar(run_label)) {
    run_label <- basename(merged_output_dir)
  }
  out_html <- file.path(merged_output_dir,
                        paste0(run_label, "_summary.html"))

  cat("Rendering merged-run summary...\n")
  cat("  Template: ", template, "\n", sep = "")
  cat("  Output  : ", out_html, "\n", sep = "")

  tryCatch({
    rmarkdown::render(
      input         = template,
      output_file   = basename(out_html),
      output_dir    = dirname(out_html),
      params        = list(
        merged_output_dir = merged_output_dir,
        run_label         = run_label
      ),
      envir         = new.env(),
      quiet         = TRUE
    )
    cat("✓ Summary report written: ", out_html, "\n\n", sep = "")
  }, error = function(e) {
    cat("⚠ Summary rendering failed (non-fatal): ", conditionMessage(e),
        "\n\n", sep = "")
  })
  invisible(out_html)
}

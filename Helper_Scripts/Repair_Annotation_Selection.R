#!/usr/bin/env Rscript
# Repair_Annotation_Selection.R
#
# Backfills missing `07_Annotation/annotation_selection.json` files by
# re-running `select_best_annotation()` against each sample's step-07
# checkpoint. Written to recover samples produced by older pipeline
# versions that didn't yet write the selection JSON, and as a safety net
# for any future silent-skip bugs.
#
# Designed for non-interactive use (Rscript). Because `interactive()`
# returns FALSE under Rscript, the patched select_best_annotation()
# (07_Annotation.R:437) will auto-select the top-scoring annotation
# without prompting.
#
# Usage:
#   Rscript Helper_Scripts/Repair_Annotation_Selection.R \
#     --output-root /mnt/home/koncha/P_lab/CosMx_analysis/Output
#
# Optional flags:
#   --samples id1,id2,id3   restrict to specific sample IDs
#   --force                 overwrite existing JSONs (default: skip)
#   --conf-threshold N      override the default 0.8

# ------------------------------------------------------------------------------
# Locate Scripts directory (same resolution logic as Manual_Single_Sample_Test.R)
# ------------------------------------------------------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 ||
                             (length(x) == 1 && is.atomic(x) && is.na(x))) y else x

this_script <- (function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) >= 1) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]),
                         winslash = "/", mustWork = FALSE))
  }
  NULL
})()

helper_dir  <- if (!is.null(this_script)) dirname(this_script) else getwd()
scripts_dir <- dirname(helper_dir)
if (!file.exists(file.path(scripts_dir, "07_Annotation.R"))) {
  scripts_dir <- Sys.getenv("COSMX_PROJECT_DIR", unset = getwd())
}
if (!file.exists(file.path(scripts_dir, "07_Annotation.R"))) {
  stop("Could not locate 07_Annotation.R. Set COSMX_PROJECT_DIR or run this ",
       "script from its checked-in location (Helper_Scripts/).")
}

source(file.path(scripts_dir, "Helper_Scripts", "Pipeline_Utils.R"))

# ------------------------------------------------------------------------------
# CLI arg parsing (lightweight)
# ------------------------------------------------------------------------------
parse_cli <- function(args) {
  out <- list(
    output_root    = NULL,
    samples        = NULL,
    force          = FALSE,
    conf_threshold = 0.8
  )
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--output-root") {
      out$output_root <- args[[i + 1L]]; i <- i + 2L
    } else if (a == "--samples") {
      out$samples <- strsplit(args[[i + 1L]], ",", fixed = TRUE)[[1]]
      out$samples <- trimws(out$samples); i <- i + 2L
    } else if (a == "--force") {
      out$force <- TRUE; i <- i + 1L
    } else if (a == "--conf-threshold") {
      out$conf_threshold <- as.numeric(args[[i + 1L]]); i <- i + 2L
    } else if (a %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript Repair_Annotation_Selection.R --output-root <path>\n",
        "       [--samples id1,id2] [--force] [--conf-threshold 0.8]\n",
        sep = ""
      )
      quit(status = 0L)
    } else {
      warning("Unknown argument: ", a, call. = FALSE)
      i <- i + 1L
    }
  }
  if (is.null(out$output_root) || !nzchar(out$output_root)) {
    stop("--output-root is required. Run with --help for usage.")
  }
  out
}

opts <- parse_cli(commandArgs(trailingOnly = TRUE))
output_root <- normalizePath(opts$output_root, winslash = "/", mustWork = TRUE)

# ------------------------------------------------------------------------------
# Load config for score_weights (optional)
# ------------------------------------------------------------------------------
config_path <- file.path(dirname(output_root), "Parameters", "config.yaml")
if (!file.exists(config_path)) {
  config_path <- file.path(scripts_dir, "Parameters", "config.yaml")
}
score_weights <- NULL
if (file.exists(config_path) && requireNamespace("yaml", quietly = TRUE)) {
  cfg <- tryCatch(yaml::read_yaml(config_path), error = function(e) NULL)
  score_weights <- cfg$parameters$annotation$score_weights
}

# ------------------------------------------------------------------------------
# Discover samples
# ------------------------------------------------------------------------------
sample_dirs <- list.dirs(output_root, recursive = FALSE, full.names = TRUE)
sample_dirs <- sample_dirs[grepl("/Sample_", sample_dirs)]
if (!is.null(opts$samples)) {
  sample_dirs <- sample_dirs[sub(".*/Sample_", "", sample_dirs) %in% opts$samples]
}
if (length(sample_dirs) == 0L) {
  stop("No Sample_* directories found under ", output_root,
       if (!is.null(opts$samples)) paste0(" matching --samples filter.") else "")
}

# ------------------------------------------------------------------------------
# Pull in 07_Annotation.R AFTER we've located the Scripts dir.
# Sourcing 07 triggers Giotto loading, so we do it once, upfront.
# ------------------------------------------------------------------------------
options(cosmx.disable_cli = TRUE)   # don't let 07 try to parse argv as its own
suppressPackageStartupMessages({
  source(file.path(scripts_dir, "07_Annotation.R"))
})
stopifnot(exists("select_best_annotation", mode = "function"))

# ------------------------------------------------------------------------------
# Per-sample repair
# ------------------------------------------------------------------------------
summary_rows <- list()
record <- function(sample_id, status, selected_column = NA_character_,
                   composite_score = NA_real_, reason = NA_character_) {
  summary_rows[[length(summary_rows) + 1L]] <<- data.frame(
    sample_id       = sample_id,
    status          = status,
    selected_column = selected_column,
    composite_score = composite_score,
    reason          = reason,
    timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

cat(sprintf("Scanning %d sample(s) under %s\n", length(sample_dirs), output_root))

for (sdir in sample_dirs) {
  sample_id <- sub(".*/Sample_", "", sdir)
  ann_dir   <- file.path(sdir, "07_Annotation")
  json_path <- file.path(ann_dir, "annotation_selection.json")
  ckpt_dir  <- file.path(sdir, "00_Checkpoints", "07_annotation")

  cat(sprintf("\n[%s]\n", sample_id))

  if (file.exists(json_path) && !isTRUE(opts$force)) {
    cat("  Existing JSON present — skipping (pass --force to overwrite).\n")
    existing <- tryCatch(jsonlite::fromJSON(json_path), error = function(e) NULL)
    record(sample_id, "skipped_existing",
           selected_column = existing$selected_annotation_column %||% NA_character_,
           composite_score = existing$composite_score %||% NA_real_,
           reason = "annotation_selection.json already exists")
    next
  }

  if (!dir.exists(ckpt_dir)) {
    cat("  No step-07 checkpoint at ", ckpt_dir, " — skipping.\n", sep = "")
    record(sample_id, "skipped_no_checkpoint",
           reason = paste("missing checkpoint dir:", ckpt_dir))
    next
  }

  dir.create(ann_dir, recursive = TRUE, showWarnings = FALSE)

  result <- tryCatch({
    cat("  Loading checkpoint…\n")
    gobj <- load_giotto_checkpoint(ckpt_dir)
    cat("  Running select_best_annotation()…\n")
    select_best_annotation(
      gobj           = gobj,
      annotation_dir = ann_dir,
      sample_id      = sample_id,
      conf_threshold = opts$conf_threshold,
      score_weights  = score_weights
    )
  }, error = function(e) {
    cat("  FAILED: ", conditionMessage(e), "\n", sep = "")
    record(sample_id, "failed", reason = conditionMessage(e))
    NULL
  })

  if (is.null(result)) next

  if (!file.exists(json_path)) {
    cat("  select_best_annotation() returned but JSON is still missing.\n")
    record(sample_id, "failed",
           reason = "select_best_annotation returned without writing JSON")
    next
  }

  sel <- tryCatch(jsonlite::fromJSON(json_path), error = function(e) NULL)
  cat(sprintf("  ✓ Wrote %s (selected %s, composite=%.3f)\n",
              basename(json_path),
              sel$selected_annotation_column %||% "<unknown>",
              sel$composite_score %||% NA_real_))
  record(sample_id, "repaired",
         selected_column = sel$selected_annotation_column %||% NA_character_,
         composite_score = sel$composite_score %||% NA_real_,
         reason = "regenerated from checkpoint")
}

# ------------------------------------------------------------------------------
# Write summary
# ------------------------------------------------------------------------------
summary_df <- do.call(rbind, summary_rows)
if (!is.null(summary_df) && nrow(summary_df) > 0L) {
  summary_path <- file.path(output_root, "annotation_repair_summary.csv")
  write.csv(summary_df, summary_path, row.names = FALSE)
  cat(sprintf("\nSummary written: %s\n", summary_path))

  tally <- table(summary_df$status)
  cat("Status counts:\n")
  for (s in names(tally)) cat(sprintf("  %-24s %d\n", s, tally[[s]]))
} else {
  cat("\nNo samples processed.\n")
}

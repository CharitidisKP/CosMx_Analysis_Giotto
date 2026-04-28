#!/usr/bin/env Rscript
# ============================================================================
# Make_Test_Subset.R
#
# Build a small, self-contained CosMx sample directory by picking the
# smallest contiguous FOV range that meets a target cell count. Useful for
# end-to-end pipeline smoke tests (~5k cells runs all 13 steps in minutes).
#
# Why FOV-contiguous (not random per-cell):
#   Steps 09 (spatial network), 10 (CCI), and 12 (smiDE) operate on cell
#   neighborhoods. Random subsampling breaks neighbor structure; contiguous
#   FOVs preserve it so the test is a faithful mini-version of a real run.
#
# Usage (inside the Apptainer container so data.table is available):
#
#   apptainer exec "$SIF" Rscript Helper_Scripts/Make_Test_Subset.R \
#     --source-dir   Data/Raw_data/MCFASCIST1250Slide11979 \
#     --target-dir   Data/Raw_data/TEST_5k \
#     --n-cells      5000
#
# Output (4 files, identical structure to a regular sample directory):
#   <target-dir>/<basename>_metadata_file.csv.gz
#   <target-dir>/<basename>_exprMat_file.csv.gz
#   <target-dir>/<basename>-polygons.csv.gz
#   <target-dir>/<basename>_fov_positions_file.csv.gz
#
# After it finishes, the script prints a sample_sheet.csv row to copy in.
# ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required. Install with install.packages('data.table').")
  }
  library(data.table)
})

parse_args <- function(args) {
  opts <- list(source_dir = NULL, target_dir = NULL, n_cells = 5000L)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) stop("Unexpected positional argument: ", arg)
    if (i == length(args)) stop("Missing value for flag: ", arg)
    key <- gsub("-", "_", sub("^--", "", arg))
    if (!(key %in% names(opts))) stop("Unknown option: ", arg)
    opts[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  for (req in c("source_dir", "target_dir")) {
    if (is.null(opts[[req]])) stop("--", gsub("_", "-", req), " is required")
  }
  opts$n_cells <- as.integer(opts$n_cells)
  if (is.na(opts$n_cells) || opts$n_cells <= 0L) stop("--n-cells must be a positive integer")
  opts
}

find_one <- function(dir, pattern, label) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0L) stop(label, " file not found in ", dir, " (pattern: ", pattern, ")")
  if (length(f) > 1L) stop("Multiple ", tolower(label), " files in ", dir, ": ",
                           paste(basename(f), collapse = ", "))
  f
}

pick_colname <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    stop("Could not find ", label, " column (tried ", paste(candidates, collapse = "/"),
         ") in columns: ", paste(names(df), collapse = ", "))
  }
  hit[1]
}

# Pick the smallest contiguous FOV window whose cumulative cell count >= n_cells.
# Sliding window over FOVs sorted ascending. Returns c(fov_min, fov_max, n_cells).
pick_fov_window <- function(fov_counts_dt, n_cells) {
  setorder(fov_counts_dt, fov)
  fovs   <- fov_counts_dt$fov
  counts <- fov_counts_dt$n
  total  <- sum(counts)
  if (total < n_cells) {
    warning(sprintf("Source has only %d cells across all FOVs; returning the full slide.", total))
    return(list(fov_min = min(fovs), fov_max = max(fovs), n_cells = total, n_fovs = length(fovs)))
  }
  best <- list(width = .Machine$integer.max, lo = NA_integer_, hi = NA_integer_, sum = 0L)
  lo <- 1L
  running <- 0L
  for (hi in seq_along(fovs)) {
    running <- running + counts[hi]
    while (running - counts[lo] >= n_cells) {
      running <- running - counts[lo]
      lo <- lo + 1L
    }
    if (running >= n_cells) {
      width <- fovs[hi] - fovs[lo] + 1L
      if (width < best$width) best <- list(width = width, lo = lo, hi = hi, sum = running)
    }
  }
  list(fov_min = fovs[best$lo], fov_max = fovs[best$hi],
       n_cells = best$sum, n_fovs = best$hi - best$lo + 1L)
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))

  if (!dir.exists(opts$source_dir)) stop("Source directory not found: ", opts$source_dir)
  source_base <- basename(normalizePath(opts$source_dir, mustWork = TRUE))
  target_base <- basename(opts$target_dir)
  if (target_base == source_base) stop("Target basename collides with source: ", source_base)

  meta_f <- find_one(opts$source_dir, "_metadata_file\\.csv(\\.gz)?$",     "Metadata")
  expr_f <- find_one(opts$source_dir, "_exprMat_file\\.csv(\\.gz)?$",      "Expression matrix")
  poly_f <- find_one(opts$source_dir, "-polygons\\.csv(\\.gz)?$",          "Polygons")
  fov_f  <- find_one(opts$source_dir, "_fov_positions_file\\.csv(\\.gz)?$", "FOV positions")

  cat("\n=== Reading source files ===\n")
  cat("  ", basename(meta_f), "\n"); meta    <- fread(meta_f)
  cat("  ", basename(expr_f), "\n"); expr    <- fread(expr_f)
  cat("  ", basename(poly_f), "\n"); poly    <- fread(poly_f)
  cat("  ", basename(fov_f),  "\n"); fov_pos <- fread(fov_f)

  meta_fov_c    <- pick_colname(meta,    c("fov", "FOV"), "metadata fov")
  fov_pos_fov_c <- pick_colname(fov_pos, c("fov", "FOV"), "fov_positions fov")

  fov_counts <- meta[, .(n = .N), by = c(meta_fov_c)]
  setnames(fov_counts, meta_fov_c, "fov")

  win <- pick_fov_window(fov_counts, opts$n_cells)
  cat(sprintf("\n=== Selected window ===\n  FOV %d-%d  (%d FOVs, %d cells; target %d)\n",
              win$fov_min, win$fov_max, win$n_fovs, win$n_cells, opts$n_cells))

  meta_keep    <- meta[meta[[meta_fov_c]] >= win$fov_min & meta[[meta_fov_c]] <= win$fov_max]
  expr_keep    <- expr[expr[["fov"]]      >= win$fov_min & expr[["fov"]]      <= win$fov_max]
  poly_keep    <- poly[poly[["fov"]]      >= win$fov_min & poly[["fov"]]      <= win$fov_max]
  fov_pos_keep <- fov_pos[fov_pos[[fov_pos_fov_c]] >= win$fov_min &
                          fov_pos[[fov_pos_fov_c]] <= win$fov_max]

  if (nrow(expr_keep) != nrow(meta_keep)) {
    warning(sprintf("Cell count mismatch: metadata=%d, expression=%d.",
                    nrow(meta_keep), nrow(expr_keep)))
  }

  dir.create(opts$target_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(meta_keep,    file.path(opts$target_dir, paste0(target_base, "_metadata_file.csv.gz")))
  fwrite(expr_keep,    file.path(opts$target_dir, paste0(target_base, "_exprMat_file.csv.gz")))
  fwrite(poly_keep,    file.path(opts$target_dir, paste0(target_base, "-polygons.csv.gz")))
  fwrite(fov_pos_keep, file.path(opts$target_dir, paste0(target_base, "_fov_positions_file.csv.gz")))

  cat(sprintf("\n=== Wrote %d files to %s ===\n", 4L, opts$target_dir))
  cat(sprintf("  cells=%d  fovs=%d  polygon_rows=%d\n",
              nrow(meta_keep), win$n_fovs, nrow(poly_keep)))

  cat("\n=== Add this row to Parameters/sample_sheet.csv ===\n")
  cat(sprintf(
    "%s,,%s,99,batch_TEST,test_pt,,TEST,Test,T0,TRUE,,%d,%d,Smoke-test subset (%d cells, FOV %d-%d) from %s\n",
    target_base, target_base, win$fov_min, win$fov_max, win$n_cells,
    win$fov_min, win$fov_max, source_base
  ))
}

main()

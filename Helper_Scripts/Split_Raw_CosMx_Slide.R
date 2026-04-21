#!/usr/bin/env Rscript
# ============================================================================
# Split_Raw_CosMx_Slide.R
#
# One-time pre-processing utility: split a composite CosMx slide directory
# (multiple biopsies on one physical slide) into self-contained per-biopsy
# sub-folders filtered by FOV range. After running this, each biopsy looks
# identical in structure to a single-sample slide and can be declared as
# its own row in sample_sheet.csv with no composite-slide special case.
#
# The pipeline reads only cell-level files, so the per-transcript
# `*_tx_file.csv.gz` is NOT copied to the sub-folders. The original slide
# folder is left untouched — this script never writes to, renames, or deletes
# anything inside --source-dir. The composite row in sample_sheet.csv keeps
# pointing at the original folder so the whole-slide view remains analysable.
#
# Usage (run inside the Apptainer container so data.table is available):
#
#   Rscript Helper_Scripts/Split_Raw_CosMx_Slide.R \
#     --source-dir  Data/Raw_data/MCFASCIST1250Slide11979876960634322 \
#     --output-root Data/Raw_data \
#     --splits "MCFASCIST1250Slide11979:55-71,\
# MCFASCIST1250Slide18769:42-54,\
# MCFASCIST1250Slide16063:29-41,\
# MCFASCIST1250Slide14320:1-28" \
#     [--dry-run]
#
# --splits is a comma-separated list of <target_folder_name>:<fov_min>-<fov_max>.
# FOV ranges are inclusive on both ends. Ranges must not overlap.
#
# Outputs (per target folder):
#   <target>/<target>_metadata_file.csv.gz
#   <target>/<target>_exprMat_file.csv.gz
#   <target>/<target>-polygons.csv.gz
#   <target>/<target>_fov_positions_file.csv.gz
# Plus one audit file at <output-root>/<source_basename>_split_audit.csv
# summarising cell / FOV / gene counts per target.
# ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required. Install with install.packages('data.table').")
  }
  library(data.table)
})

# ---- CLI -------------------------------------------------------------------
parse_args <- function(args) {
  opts <- list(
    source_dir  = NULL,
    output_root = NULL,
    splits      = NULL,
    dry_run     = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--dry-run") { opts$dry_run <- TRUE; i <- i + 1L; next }
    if (!startsWith(arg, "--")) stop("Unexpected positional argument: ", arg)
    if (i == length(args)) stop("Missing value for flag: ", arg)
    key <- gsub("-", "_", sub("^--", "", arg))
    if (!(key %in% names(opts))) stop("Unknown option: ", arg)
    opts[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  for (req in c("source_dir", "output_root", "splits")) {
    if (is.null(opts[[req]])) stop("--", gsub("_", "-", req), " is required")
  }
  opts
}

parse_splits <- function(splits_str) {
  parts <- trimws(strsplit(splits_str, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) stop("--splits is empty")
  out <- lapply(parts, function(p) {
    m <- regmatches(p, regexec("^([^:]+):(\\d+)-(\\d+)$", p))[[1]]
    if (length(m) != 4L) {
      stop("Invalid split spec '", p, "'. Expected <target>:<fov_min>-<fov_max>.")
    }
    list(target = m[2], fov_min = as.integer(m[3]), fov_max = as.integer(m[4]))
  })
  tbl <- data.table::rbindlist(out)
  if (any(tbl$fov_min > tbl$fov_max)) {
    bad <- tbl[fov_min > fov_max]
    stop("Invalid range in: ", paste(bad$target, collapse = ", "))
  }
  if (anyDuplicated(tbl$target)) {
    stop("Duplicate target directory: ",
         paste(unique(tbl$target[duplicated(tbl$target)]), collapse = ", "))
  }
  # Overlap check
  ord <- order(tbl$fov_min)
  o <- tbl[ord]
  for (j in seq_len(nrow(o) - 1L)) {
    if (o$fov_max[j] >= o$fov_min[j + 1L]) {
      stop("FOV ranges overlap between '", o$target[j], "' and '",
           o$target[j + 1L], "'")
    }
  }
  tbl
}

find_one <- function(dir, pattern, label) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(f) == 0L) stop(label, " file not found in ", dir,
                            " (pattern: ", pattern, ")")
  if (length(f) > 1L) stop("Multiple ", tolower(label), " files in ", dir, ": ",
                           paste(basename(f), collapse = ", "))
  f
}

pick_colname <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    stop("Could not find ", label, " column (tried ",
         paste(candidates, collapse = "/"), ") in columns: ",
         paste(names(df), collapse = ", "))
  }
  hit[1]
}

# ---- Main ------------------------------------------------------------------
main <- function() {
  opts   <- parse_args(commandArgs(trailingOnly = TRUE))
  splits <- parse_splits(opts$splits)

  if (!dir.exists(opts$source_dir)) {
    stop("Source directory not found: ", opts$source_dir)
  }
  if (!dir.exists(opts$output_root)) {
    stop("Output root not found: ", opts$output_root)
  }

  source_base <- basename(normalizePath(opts$source_dir, mustWork = TRUE))
  cat("\n=== Split plan ===\n")
  cat(sprintf("  Source: %s\n", opts$source_dir))
  cat(sprintf("  Output root: %s\n", opts$output_root))
  for (i in seq_len(nrow(splits))) {
    cat(sprintf("    %-40s  FOV %3d-%3d\n",
                splits$target[i], splits$fov_min[i], splits$fov_max[i]))
  }

  # Validate target folders do not collide with the source folder.
  if (any(splits$target == source_base)) {
    stop("Target directory '", source_base,
         "' collides with the source directory basename. ",
         "Choose a different target name.")
  }

  if (isTRUE(opts$dry_run)) {
    cat("\n=== Dry run — no files written ===\n")
    return(invisible(NULL))
  }

  # -- Locate and read source files (one pass) ------------------------------
  meta_f <- find_one(opts$source_dir, "_metadata_file\\.csv(\\.gz)?$", "Metadata")
  expr_f <- find_one(opts$source_dir, "_exprMat_file\\.csv(\\.gz)?$", "Expression matrix")
  poly_f <- find_one(opts$source_dir, "-polygons\\.csv(\\.gz)?$",      "Polygons")
  fov_f  <- find_one(opts$source_dir, "_fov_positions_file\\.csv(\\.gz)?$", "FOV positions")

  cat("\n=== Reading source files ===\n")
  cat("  ", basename(meta_f), "\n"); meta    <- data.table::fread(meta_f)
  cat("  ", basename(expr_f), "\n"); expr    <- data.table::fread(expr_f)
  cat("  ", basename(poly_f), "\n"); poly    <- data.table::fread(poly_f)
  cat("  ", basename(fov_f),  "\n"); fov_pos <- data.table::fread(fov_f)

  meta_fov_c    <- pick_colname(meta,    c("fov", "FOV"),          "metadata fov")
  meta_cell_c   <- pick_colname(meta,    c("cell_ID", "cellID"),   "metadata cell_ID")
  expr_cell_c   <- pick_colname(expr,    c("cell_ID", "cellID"),   "expression cell_ID")
  poly_cell_c   <- pick_colname(poly,    c("cellID", "cell_ID"),   "polygons cellID")
  fov_pos_fov_c <- pick_colname(fov_pos, c("fov", "FOV"),          "fov_positions fov")

  # CosMx exports sometimes store cell_ID as integer in one file and character
  # in another. Coerce both sides of every comparison with as.character() so
  # %chin% never hits a type mismatch. Values written to disk are unchanged
  # because we only coerce inside the comparison — original columns are intact.

  # Cross-file sanity: every cell in the expression matrix should exist in metadata.
  if (!all(as.character(expr[[expr_cell_c]]) %chin% as.character(meta[[meta_cell_c]]))) {
    warning("Some cells in the expression matrix have no metadata row. ",
            "They will be dropped from every split.")
  }

  audit <- data.table::data.table()

  for (i in seq_len(nrow(splits))) {
    tgt  <- splits$target[i]
    fmin <- splits$fov_min[i]
    fmax <- splits$fov_max[i]
    out_dir <- file.path(opts$output_root, tgt)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    cat(sprintf("\n--- %s (FOV %d-%d) ---\n", tgt, fmin, fmax))

    meta_keep    <- meta[meta[[meta_fov_c]] >= fmin & meta[[meta_fov_c]] <= fmax]
    keep_cells   <- as.character(meta_keep[[meta_cell_c]])
    if (length(keep_cells) == 0L) {
      stop("No cells match FOV range ", fmin, "-", fmax, " for target '", tgt,
           "'. Check --splits spec.")
    }
    expr_keep    <- expr[as.character(expr[[expr_cell_c]]) %chin% keep_cells]
    poly_keep    <- poly[as.character(poly[[poly_cell_c]]) %chin% keep_cells]
    fov_pos_keep <- fov_pos[fov_pos[[fov_pos_fov_c]] >= fmin &
                            fov_pos[[fov_pos_fov_c]] <= fmax]

    if (nrow(expr_keep) != nrow(meta_keep)) {
      warning(sprintf(
        "[%s] cell count mismatch: metadata=%d, expression=%d. ",
        tgt, nrow(meta_keep), nrow(expr_keep)))
    }

    fwrite(meta_keep,    file.path(out_dir, paste0(tgt, "_metadata_file.csv.gz")))
    fwrite(expr_keep,    file.path(out_dir, paste0(tgt, "_exprMat_file.csv.gz")))
    fwrite(poly_keep,    file.path(out_dir, paste0(tgt, "-polygons.csv.gz")))
    fwrite(fov_pos_keep, file.path(out_dir, paste0(tgt, "_fov_positions_file.csv.gz")))

    n_fovs_kept      <- length(unique(meta_keep[[meta_fov_c]]))
    n_fovs_requested <- fmax - fmin + 1L
    cat(sprintf("  cells=%d  fovs=%d/%d  polygons_rows=%d\n",
                nrow(meta_keep), n_fovs_kept, n_fovs_requested, nrow(poly_keep)))

    audit <- rbind(audit, data.table::data.table(
      target_directory = tgt,
      fov_min          = fmin,
      fov_max          = fmax,
      n_cells          = nrow(meta_keep),
      n_fovs_kept      = n_fovs_kept,
      n_fovs_requested = n_fovs_requested,
      n_polygon_rows   = nrow(poly_keep)
    ))
  }

  audit_path <- file.path(opts$output_root, paste0(source_base, "_split_audit.csv"))
  fwrite(audit, audit_path)
  cat("\n=== Audit written to: ", audit_path, " ===\n", sep = "")
  print(audit)

  # Cross-split accounting: total cells should match the FOVs we touched.
  covered <- meta[meta[[meta_fov_c]] >= min(splits$fov_min) &
                  meta[[meta_fov_c]] <= max(splits$fov_max)]
  covered <- covered[
    Reduce("|", Map(function(lo, hi) covered[[meta_fov_c]] >= lo & covered[[meta_fov_c]] <= hi,
                    splits$fov_min, splits$fov_max))
  ]
  if (nrow(covered) != sum(audit$n_cells)) {
    warning(sprintf(
      "Cross-split cell total (%d) does not match the union of requested FOVs (%d). ",
      sum(audit$n_cells), nrow(covered)))
  }
}

main()

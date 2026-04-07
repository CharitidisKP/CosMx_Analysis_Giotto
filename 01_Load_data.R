# Load CosMx data and create Giotto object --------------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 01_Load_data.R
# Load CosMx data and create Giotto object with provided polygons
# ==============================================================================

.find_single_cosmx_file <- function(data_dir, pattern, label) {
  matches <- list.files(
    data_dir,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(matches) == 0) {
    stop(label, " file not found in ", data_dir)
  }
  
  if (length(matches) > 1) {
    stop(
      "Expected one ", label, " file in ", data_dir,
      ", found ", length(matches), ": ",
      paste(basename(matches), collapse = ", ")
    )
  }
  
  matches[[1]]
}

.read_cosmx_csv <- function(path, ...) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    readr::read_csv(gzfile(path), show_col_types = FALSE, ...)
  } else {
    readr::read_csv(path, show_col_types = FALSE, ...)
  }
}

.require_columns <- function(df, required, label) {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

.infer_sample_number <- function(cell_ids) {
  cell_ids <- cell_ids[!is.na(cell_ids) & nzchar(cell_ids)]
  if (length(cell_ids) == 0) {
    return(1L)
  }
  
  parsed <- suppressWarnings(as.integer(sub("^c_(\\d+)_.*$", "\\1", cell_ids)))
  parsed <- parsed[!is.na(parsed)]
  
  if (length(parsed) == 0) {
    return(1L)
  }
  
  modal_idx <- which.max(tabulate(match(parsed, unique(parsed))))
  unique(parsed)[modal_idx]
}

.align_expression_and_metadata <- function(expr_data, metadata, sample_number) {
  expr_data <- expr_data %>%
    dplyr::mutate(giotto_cell_ID = paste0("c_", sample_number, "_", fov, "_", cell_ID))
  
  common_ids <- intersect(expr_data$giotto_cell_ID, metadata$cell_id)
  cat("Matching cell IDs:", length(common_ids), "expression/metadata overlaps\n")
  
  if (length(common_ids) == 0) {
    stop("No overlapping cell IDs between expression matrix and metadata.")
  }
  
  if (length(common_ids) < nrow(expr_data) || length(common_ids) < nrow(metadata)) {
    cat(
      "⚠ Retaining the shared cells only:",
      length(common_ids), "of", nrow(expr_data), "expression rows and",
      nrow(metadata), "metadata rows\n"
    )
  } else {
    cat("✓ 100% cell ID match\n")
  }
  
  expr_data <- expr_data %>%
    dplyr::filter(giotto_cell_ID %in% common_ids) %>%
    dplyr::arrange(match(giotto_cell_ID, common_ids))
  
  metadata <- metadata %>%
    dplyr::filter(cell_id %in% common_ids) %>%
    dplyr::arrange(match(cell_id, common_ids)) %>%
    dplyr::mutate(giotto_cell_ID = cell_id)
  
  list(expr_data = expr_data, metadata = metadata)
}

.write_fov_report <- function(metadata, fov_positions, results_folder, sample_id) {
  if (!all(c("fov", "CenterX_global_px", "CenterY_global_px") %in% names(metadata))) {
    return(invisible(NULL))
  }
  
  if (!all(c("FOV", "x_global_px", "y_global_px") %in% names(fov_positions))) {
    return(invisible(NULL))
  }
  
  fov_summary <- metadata %>%
    dplyr::mutate(fov = suppressWarnings(as.integer(fov))) %>%
    dplyr::group_by(fov) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_center_x_px = stats::median(CenterX_global_px, na.rm = TRUE),
      median_center_y_px = stats::median(CenterY_global_px, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      fov_positions %>%
        dplyr::transmute(
          fov = suppressWarnings(as.integer(FOV)),
          fov_origin_x_px = suppressWarnings(as.numeric(x_global_px)),
          fov_origin_y_px = suppressWarnings(as.numeric(y_global_px)),
          fov_origin_x_mm = suppressWarnings(as.numeric(x_global_mm)),
          fov_origin_y_mm = suppressWarnings(as.numeric(y_global_mm))
        ),
      by = "fov"
    ) %>%
    dplyr::mutate(
      delta_center_x_px = median_center_x_px - fov_origin_x_px,
      delta_center_y_px = median_center_y_px - fov_origin_y_px
    ) %>%
    dplyr::arrange(fov)
  
  readr::write_csv(
    fov_summary,
    file.path(results_folder, paste0(sample_id, "_fov_summary.csv"))
  )
  
  invisible(fov_summary)
}

.extract_spatial_locations_from_metadata <- function(metadata) {
  x_col <- c("CenterX_global_px", "x_global_px", "CenterX_local_px")[c("CenterX_global_px", "x_global_px", "CenterX_local_px") %in% names(metadata)][1]
  y_col <- c("CenterY_global_px", "y_global_px", "CenterY_local_px")[c("CenterY_global_px", "y_global_px", "CenterY_local_px") %in% names(metadata)][1]
  
  if (is.na(x_col) || is.na(y_col)) {
    cat("⚠ No coordinate columns found in metadata; Giotto object will be created without explicit spatial_locs\n\n")
    return(NULL)
  }
  
  spatlocs <- metadata %>%
    dplyr::transmute(
      cell_ID = giotto_cell_ID,
      sdimx = suppressWarnings(as.numeric(.data[[x_col]])),
      sdimy = suppressWarnings(as.numeric(.data[[y_col]]))
    ) %>%
    dplyr::filter(!is.na(sdimx) & !is.na(sdimy))
  
  if (nrow(spatlocs) == 0) {
    cat("⚠ Metadata coordinates were empty after numeric conversion; Giotto object will be created without explicit spatial_locs\n\n")
    return(NULL)
  }
  
  as.data.frame(spatlocs, stringsAsFactors = FALSE)
}

#' Internal: Giotto object integrity report
#'
#' Checks polygon <-> cell ID consistency and per-FOV cell counts.
#' Writes two CSVs to results_folder and prints a console summary.
#' Stops with an informative error if overlap falls below min_overlap_pct.
#'
#' @param gobj            Giotto object (post-polygon attachment)
#' @param sample_id       Sample identifier string
#' @param results_folder  Directory to write the CSV reports
#' @param min_overlap_pct Minimum acceptable polygon<->cell ID overlap (0-100).
#'                        Default 80. Set to 0 to warn only.
#' @return Invisibly returns a list: overview, fov_cells, overlap_pct,
#'         n_poly_only, n_cell_only

.giotto_object_report <- function(gobj,
                                  sample_id,
                                  results_folder,
                                  min_overlap_pct = 80) {
  
  cat("\n--- Object integrity report ---\n")
  
  # ── Basic counts ─────────────────────────────────────────────────────────────
  n_cells <- tryCatch(length(gobj@cell_ID$cell), error = function(e) NA_integer_)
  n_feats <- tryCatch(length(gobj@feat_ID$rna),  error = function(e) NA_integer_)
  
  # ── Polygon presence ────────────────────────────────────────────────────────��─
  spat_names <- tryCatch(names(gobj@spatial_info), error = function(e) character(0))
  has_polys  <- "cell" %in% spat_names && !is.null(gobj@spatial_info$cell)
  
  n_polys <- if (has_polys) {
    tryCatch(
      length(gobj@spatial_info$cell@unique_ID_cache),
      error = function(e) NA_integer_
    )
  } else { 0L }
  
  # ── ID overlap ────────────────────────────────────────────────────────────────
  poly_ids <- if (has_polys) {
    tryCatch(gobj@spatial_info$cell@unique_ID_cache, error = function(e) character(0))
  } else { character(0) }
  
  cell_ids <- tryCatch(gobj@cell_ID$cell, error = function(e) character(0))
  
  n_overlap   <- length(intersect(poly_ids, cell_ids))
  n_poly_only <- length(setdiff(poly_ids, cell_ids))   # orphan polygons
  n_cell_only <- length(setdiff(cell_ids, poly_ids))   # cells missing a polygon
  
  overlap_pct <- if (length(cell_ids) > 0) {
    round(100 * n_overlap / length(cell_ids), 1)
  } else { 0 }
  
  # ── Per-FOV cell counts ────────────────────────────────────────────────────────
  fov_tab <- tryCatch({
    fov_str <- sub("^c_(\\d+)_(\\d+)_.*$", "\\2", cell_ids)
    fov_df  <- data.frame(FOV = as.integer(fov_str))
    fov_df  <- aggregate(list(n_cells = fov_df$FOV), list(FOV = fov_df$FOV), length)
    fov_df[order(fov_df$FOV), ]
  }, error = function(e) data.frame(FOV = integer(), n_cells = integer()))
  
  # ── Console output ────────────────────────────────────────────────────────────
  cat(sprintf("  Cells loaded:           %d\n", n_cells))
  cat(sprintf("  RNA features:           %d\n", n_feats))
  cat(sprintf("  Cell polygons:          %d\n", n_polys))
  cat(sprintf("  Polygon \u2229 cell overlap:  %d / %d cells (%.1f%%)\n",
              n_overlap, length(cell_ids), overlap_pct))
  
  if (n_poly_only > 0)
    cat(sprintf("  \u26A0 Orphan polygons (no matching cell): %d\n", n_poly_only))
  if (n_cell_only > 0)
    cat(sprintf("  \u26A0 Cells missing polygon:              %d\n", n_cell_only))
  
  cat(sprintf("  FOVs represented:       %d\n", nrow(fov_tab)))
  if (nrow(fov_tab) > 0)
    cat(sprintf("  Cells/FOV (min/med/max): %d / %g / %d\n",
                min(fov_tab$n_cells),
                stats::median(fov_tab$n_cells),
                max(fov_tab$n_cells)))
  
  # ── Write CSV reports ──────────────────────────────────────────────────────────
  overview_df <- data.frame(
    Metric = c("Cells", "RNA features", "Cell polygons",
               "Polygon-cell overlap (n)", "Polygon-cell overlap (%)",
               "Orphan polygons", "Cells missing polygon", "FOVs"),
    Value  = c(n_cells, n_feats, n_polys,
               n_overlap, overlap_pct,
               n_poly_only, n_cell_only, nrow(fov_tab))
  )
  
  write.csv(overview_df,
            file.path(results_folder, paste0(sample_id, "_object_overview.csv")),
            row.names = FALSE)
  
  if (nrow(fov_tab) > 0)
    write.csv(fov_tab,
              file.path(results_folder, paste0(sample_id, "_cells_per_fov.csv")),
              row.names = FALSE)
  
  # ── Hard stop if overlap is critically low ─────────────────────────────────────
  if (min_overlap_pct > 0 && overlap_pct < min_overlap_pct) {
    stop(sprintf(
      paste0("Polygon<->cell ID overlap is %.1f%% (threshold: %.0f%%).\n",
             "  Check polygon IDs were built with the correct sample_id prefix.\n",
             "  Orphan polygons: %d  |  Cells missing polygon: %d"),
      overlap_pct, min_overlap_pct, n_poly_only, n_cell_only
    ))
  }
  
  cat("\u2713 Object integrity report complete\n")
  
  invisible(list(
    overview    = overview_df,
    fov_cells   = fov_tab,
    overlap_pct = overlap_pct,
    n_poly_only = n_poly_only,
    n_cell_only = n_cell_only
  ))
}


#' Load CosMx Sample with Provided Polygons
#'
#' @param sample_id       Sample identifier
#' @param data_dir        Path to raw data directory
#' @param output_dir      Output directory for this sample
#' @param min_overlap_pct Minimum polygon<->cell ID overlap %. Default 80.
#'                        Set to 0 to warn only (no hard stop).
#' @return Giotto object with polygons

load_cosmx_sample <- function(sample_id, data_dir, output_dir,
                              min_overlap_pct = 80) {
  
  cat("\n========================================\n")
  cat("STEP 01: Loading CosMx Data\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  # Create output directories
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results_folder <- file.path(output_dir, "01_Data_Loading")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Set Giotto instructions
  instructions <- createGiottoInstructions(
    save_dir = results_folder,
    save_plot = TRUE,
    show_plot = FALSE,
    return_plot = FALSE,
    python_path = python_path
  )
  
  # Find required files
  cat("Searching for data files...\n")
  
  expr_file <- .find_single_cosmx_file(data_dir, "_exprMat_file\\.csv(\\.gz)?$", "Expression matrix")
  meta_file <- .find_single_cosmx_file(data_dir, "_metadata_file\\.csv(\\.gz)?$", "Metadata")
  fov_file  <- .find_single_cosmx_file(data_dir, "_fov_positions_file\\.csv(\\.gz)?$", "FOV positions")
  poly_file <- .find_single_cosmx_file(data_dir, "-polygons\\.csv(\\.gz)?$", "Polygons")
  
  required_files <- c(expr_file, meta_file, fov_file, poly_file)
  for (file_path in required_files) {
    cat("✓", basename(file_path), "\n")
  }
  cat("\n")
  
  # Load expression matrix
  cat("Loading expression matrix...\n")
  expr_data <- .read_cosmx_csv(expr_file)
  .require_columns(expr_data, c("fov", "cell_ID"), "Expression matrix")
  
  cat("Initial dimensions:", nrow(expr_data), "rows x", ncol(expr_data), "columns\n")
  
  # Load metadata
  cat("\nLoading metadata...\n")
  metadata <- .read_cosmx_csv(meta_file)
  .require_columns(
    metadata,
    c("cell_id", "fov", "cell_ID"),
    "Metadata"
  )
  
  cat("Metadata dimensions:", nrow(metadata), "rows x", ncol(metadata), "columns\n")
  
  cat("\nLoading FOV positions...\n")
  fov_positions <- .read_cosmx_csv(fov_file)
  .require_columns(fov_positions, c("FOV", "x_global_px", "y_global_px"), "FOV positions")
  cat("FOV positions:", nrow(fov_positions), "rows\n")
  
  # Get sample number from metadata cell_id
  sample_number <- .infer_sample_number(metadata$cell_id)
  cat("Sample number:", sample_number, "\n\n")
  
  cat("Aligning expression and metadata on CosMx cell IDs...\n")
  aligned <- .align_expression_and_metadata(
    expr_data = expr_data,
    metadata = metadata,
    sample_number = sample_number
  )
  expr_data <- aligned$expr_data
  metadata <- aligned$metadata
  cat("\n")
  
  # Extract gene columns
  meta_cols    <- c("fov", "cell_ID", "giotto_cell_ID")
  gene_col_idx <- which(!names(expr_data) %in% meta_cols)
  gene_names   <- names(expr_data)[gene_col_idx]
  
  cat("Genes:", length(gene_names), "\n")
  cat("Cells:", nrow(expr_data), "\n")
  
  non_numeric_genes <- gene_names[!vapply(expr_data[gene_col_idx], is.numeric, logical(1))]
  if (length(non_numeric_genes) > 0) {
    stop(
      "Expression matrix contains non-numeric gene columns: ",
      paste(head(non_numeric_genes, 10), collapse = ", ")
    )
  }
  
  # Create expression matrix (genes x cells)
  expr_matrix_ct <- as.matrix(expr_data[, gene_col_idx])
  expr_matrix    <- Matrix::Matrix(t(expr_matrix_ct), sparse = TRUE)
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- expr_data$giotto_cell_ID
  rm(expr_matrix_ct)
  
  # Handle duplicate gene names
  if (any(duplicated(rownames(expr_matrix)))) {
    n_dups <- sum(duplicated(rownames(expr_matrix)))
    cat("⚠ Making", n_dups, "duplicate gene names unique...\n")
    rownames(expr_matrix) <- make.unique(rownames(expr_matrix), sep = "_")
  }
  
  cat("✓ Expression matrix prepared:", nrow(expr_matrix), "genes x",
      ncol(expr_matrix), "cells\n\n")
  
  # Prepare metadata
  metadata <- metadata %>%
    dplyr::mutate(giotto_cell_ID = cell_id)
  
  spatial_locs <- .extract_spatial_locations_from_metadata(metadata)
  
  # Create Giotto object
  cat("Creating Giotto object...\n")
  
  cosmx <- createGiottoObject(
    expression   = expr_matrix,
    spatial_locs = spatial_locs,
    instructions = instructions
  )
  
  cat("✓ Giotto object created\n\n")
  
  # Add metadata
  cat("Adding metadata...\n")
  
  cosmx <- addCellMetadata(
    gobject        = cosmx,
    new_metadata   = metadata,
    by_column      = TRUE,
    column_cell_ID = "giotto_cell_ID"
  )
  
  cat("✓ Metadata added\n\n")
  
  if (!is.null(spatial_locs)) {
    cat("✓ Spatial locations registered during object creation\n\n")
  }
  
  .write_fov_report(
    metadata = metadata,
    fov_positions = fov_positions,
    results_folder = results_folder,
    sample_id = sample_id
  )
  
  # Add polygons
  cat("Loading polygons...\n")
  cosmx <- add_polygons_from_csv(cosmx, poly_file)
  cat("✓ Polygons added\n\n")
  
  # ── NEW: Object integrity report ──────────────────────────────────────────────
  .giotto_object_report(
    gobj            = cosmx,
    sample_id       = sample_id,
    results_folder  = results_folder,
    min_overlap_pct = min_overlap_pct
  )
  
  cleanup_fn <- get0("cleanup_memory", mode = "function", inherits = TRUE)
  if (is.function(cleanup_fn)) {
    cleanup_fn(
      remove = c("expr_data", "metadata", "aligned", "spatial_locs", "fov_positions"),
      envir = environment(),
      label = "Data loading",
      verbose = TRUE
    )
  } else {
    gc(verbose = FALSE)
  }
  
  cat("\n✓ Data loading complete\n")
  
  return(cosmx)
}


#' Add polygons from CSV (simplified version using existing cell IDs)
add_polygons_from_csv <- function(gobj, polygon_file) {
  
  # Read polygon data
  poly_data <- .read_cosmx_csv(polygon_file)
  .require_columns(
    poly_data,
    c("cell", "x_global_px", "y_global_px"),
    "Polygon file"
  )
  
  cat("✓ Polygon data loaded:", nrow(poly_data), "vertices\n")
  
  # The 'cell' column already has the proper IDs!
  cat("Using existing cell IDs from polygon file\n")
  
  # Group by cell and create polygons
  poly_list <- poly_data %>%
    dplyr::group_by(cell) %>%
    dplyr::summarise(
      coords = list(cbind(x_global_px, y_global_px)),
      .groups = "drop"
    )
  
  cat("Creating", nrow(poly_list), "polygons...\n")
  
  # Convert to terra SpatVector
  vects <- lapply(seq_len(nrow(poly_list)), function(i) {
    coords <- poly_list$coords[[i]]
    
    # Close polygon
    if (nrow(coords) > 2 && !all(coords[1,] == coords[nrow(coords),])) {
      coords <- rbind(coords, coords[1,])
    }
    
    tryCatch({
      terra::vect(coords, type = "polygons")
    }, error = function(e) NULL)
  })
  
  # Remove invalid
  valid_idx <- !sapply(vects, is.null)
  vects     <- vects[valid_idx]
  poly_list <- poly_list[valid_idx, ]
  
  if (length(vects) == 0) {
    stop("No valid polygons created")
  }
  
  cat("✓ Created", length(vects), "valid polygons\n")
  
  # Combine all polygons
  spat_vect <- do.call(rbind, vects)
  terra::values(spat_vect) <- data.frame(
    poly_ID = as.character(poly_list$cell)
  )
  
  # Create GiottoPolygon
  cat("Adding polygons to Giotto object...\n")
  
  gpolygon <- new("giottoPolygon",
                  spatVector          = spat_vect,
                  spatVectorCentroids = terra::centroids(spat_vect),
                  overlaps            = data.table(),
                  name                = "cell",
                  unique_ID_cache     = as.character(poly_list$cell))
  
  if (exists("setPolygonInfo", mode = "function")) {
    gobj <- setPolygonInfo(
      gobject = gobj,
      x = list(cell = gpolygon),
      centroids_to_spatlocs = FALSE
    )
  } else {
    gobj <- addGiottoPolygons(
      gobject   = gobj,
      gpolygons = list(cell = gpolygon)
    )
    
    tryCatch({
      gobj <- addSpatialCentroidLocations(
        gobject   = gobj,
        poly_info = "cell"
      )
    }, error = function(e) {
      cat("⚠ Could not add centroid locations (non-critical)\n")
    })
  }
  
  return(gobj)
}


# Run if sourced directly
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) {
    script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
    script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
    bootstrap_script <- file.path(script_dir, "Helper_Scripts", "Script_Bootstrap.R")
    if (file.exists(bootstrap_script)) {
      source(bootstrap_script, local = .GlobalEnv)
      bootstrap_pipeline_environment(script_dir, load_pipeline_utils = TRUE, verbose = FALSE)
    }
    
    sample_id  <- args[1]
    data_dir   <- args[2]
    output_dir <- args[3]
    cosmx <- load_cosmx_sample(sample_id, data_dir, output_dir)
    if (exists("save_giotto_checkpoint")) {
      save_giotto_checkpoint(
        gobj = cosmx,
        checkpoint_dir = file.path(output_dir, "Giotto_Object_Loaded"),
        metadata = list(stage = "load", sample_id = sample_id)
      )
    }
    saveRDS(cosmx, file.path(output_dir, paste0(sample_id, "_cosmx_loaded.rds")))
  } else {
    stop("Usage: Rscript 01_Load_data.R <sample_id> <data_dir> <output_dir>")
  }
}

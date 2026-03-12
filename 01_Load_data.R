# Load CosMx data and create Giotto object --------------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 01_Load_data.R
# Load CosMx data and create Giotto object with provided polygons
# ==============================================================================

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
  
  expr_file <- list.files(data_dir, pattern = "_exprMat_file\\.csv",
                          full.names = TRUE)[1]
  meta_file <- list.files(data_dir, pattern = "_metadata_file\\.csv",
                          full.names = TRUE)[1]
  fov_file  <- list.files(data_dir, pattern = "_fov_positions_file\\.csv",
                          full.names = TRUE)[1]
  poly_file <- list.files(data_dir, pattern = "-polygons\\.csv",
                          full.names = TRUE)[1]
  
  # Validate files exist
  required_files <- list(
    "Expression matrix" = expr_file,
    "Metadata"          = meta_file,
    "FOV positions"     = fov_file,
    "Polygons"          = poly_file
  )
  
  for (file_type in names(required_files)) {
    file_path <- required_files[[file_type]]
    if (is.na(file_path) || !file.exists(file_path)) {
      stop(file_type, " file not found in ", data_dir)
    }
    cat("✓", file_type, ":", basename(file_path), "\n")
  }
  
  cat("\n")
  
  # Load expression matrix
  cat("Loading expression matrix...\n")
  if (grepl("\\.gz$", expr_file)) {
    expr_data <- read_csv(gzfile(expr_file), show_col_types = FALSE)
  } else {
    expr_data <- read_csv(expr_file, show_col_types = FALSE)
  }
  
  cat("Initial dimensions:", nrow(expr_data), "rows x", ncol(expr_data), "columns\n")
  
  # Load metadata
  cat("\nLoading metadata...\n")
  if (grepl("\\.gz$", meta_file)) {
    metadata <- read_csv(gzfile(meta_file), show_col_types = FALSE)
  } else {
    metadata <- read_csv(meta_file, show_col_types = FALSE)
  }
  
  cat("Metadata dimensions:", nrow(metadata), "rows x", ncol(metadata), "columns\n")
  
  # Get sample number from metadata cell_id
  example_id <- metadata$cell_id[1]
  id_parts <- str_split(example_id, "_")[[1]]
  sample_number <- if (length(id_parts) >= 4 && id_parts[1] == "c") {
    as.integer(id_parts[2])
  } else {
    1
  }
  
  cat("Sample number:", sample_number, "\n\n")
  
  # Create matching cell IDs in expression data
  cat("Creating cell IDs in expression data...\n")
  
  expr_data <- expr_data %>%
    mutate(
      giotto_cell_ID = paste0("c_", sample_number, "_", fov, "_", cell_ID)
    )
  
  # Verify matching
  common_ids <- intersect(expr_data$giotto_cell_ID, metadata$cell_id)
  cat("Matching cell IDs:", length(common_ids), "out of", nrow(expr_data), "\n")
  
  if (length(common_ids) == nrow(expr_data)) {
    cat("✓ 100% cell ID match\n\n")
  }
  
  # Extract gene columns
  meta_cols    <- c("fov", "cell_ID", "giotto_cell_ID")
  gene_col_idx <- which(!names(expr_data) %in% meta_cols)
  gene_names   <- names(expr_data)[gene_col_idx]
  
  cat("Genes:", length(gene_names), "\n")
  cat("Cells:", nrow(expr_data), "\n")
  
  # Create expression matrix (genes x cells)
  expr_matrix_ct <- as.matrix(expr_data[, gene_col_idx])
  expr_matrix    <- t(expr_matrix_ct)
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- expr_data$giotto_cell_ID
  
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
    mutate(giotto_cell_ID = cell_id)
  
  # Create Giotto object
  cat("Creating Giotto object...\n")
  
  cosmx <- createGiottoObject(
    expression   = expr_matrix,
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
  
  cat("\n✓ Data loading complete\n")
  
  return(cosmx)
}


#' Add polygons from CSV (simplified version using existing cell IDs)
add_polygons_from_csv <- function(gobj, polygon_file) {
  
  # Read polygon data
  if (grepl("\\.gz$", polygon_file)) {
    poly_data <- read_csv(gzfile(polygon_file), show_col_types = FALSE)
  } else {
    poly_data <- read_csv(polygon_file, show_col_types = FALSE)
  }
  
  cat("✓ Polygon data loaded:", nrow(poly_data), "vertices\n")
  
  # The 'cell' column already has the proper IDs!
  cat("Using existing cell IDs from polygon file\n")
  
  # Group by cell and create polygons
  poly_list <- poly_data %>%
    group_by(cell) %>%
    summarise(
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
  
  # Add to Giotto object
  gobj <- addGiottoPolygons(
    gobject   = gobj,
    gpolygons = list(cell = gpolygon)
  )
  
  # Calculate and add centroids as spatial locations
  tryCatch({
    gobj <- addSpatialCentroidLocations(
      gobject   = gobj,
      poly_info = "cell"
    )
  }, error = function(e) {
    cat("⚠ Could not add centroid locations (non-critical)\n")
  })
  
  return(gobj)
}


# Run if sourced directly
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) {
    sample_id  <- args[1]
    data_dir   <- args[2]
    output_dir <- args[3]
    cosmx <- load_cosmx_sample(sample_id, data_dir, output_dir)
    saveRDS(cosmx, file.path(output_dir, paste0(sample_id, "_cosmx_loaded.rds")))
  }
}
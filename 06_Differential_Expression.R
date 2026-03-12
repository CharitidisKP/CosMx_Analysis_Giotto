# Find cluster-specific and DEG marker genes -----------------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 06_marker_analysis.R
# Find cluster marker genes
# ==============================================================================

#' Marker Gene Analysis
#'
#' @param gobj Giotto object or path
#' @param sample_id Sample identifier
#' @param output_dir Output directory
#' @param cluster_column Clustering column name
#' @param top_n Top N markers per cluster
#' @return Giotto object with marker results

marker_analysis <- function(gobj,
                            sample_id,
                            output_dir,
                            cluster_column = "leiden_clust",
                            top_n = 25) {
  
  cat("\n========================================\n")
  cat("STEP 06: Marker Analysis\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }
  
  results_folder <- file.path(output_dir, "06_Markers")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  cat("Finding marker genes...\n")
  cat("  Cluster column:", cluster_column, "\n")
  cat("  Top N per cluster:", top_n, "\n\n")
  
  # Find markers using scran
  markers_scran <- findMarkers_one_vs_all(
    gobject = gobj,
    method = "scran",
    expression_values = "normalized",
    cluster_column = cluster_column,
    min_feats = 3
  )
  
  cat("✓ Markers found\n\n")
  
  # Extract top markers per cluster
  top_markers <- markers_scran %>%
    as_tibble() %>%
    group_by(cluster) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # Save all markers
  write_csv(markers_scran, 
            file.path(results_folder, paste0(sample_id, "_all_markers.csv")))
  
  # Save top markers
  write_csv(top_markers,
            file.path(results_folder, paste0(sample_id, "_top_markers.csv")))
  
  cat("✓ Marker tables saved\n\n")
  
  # Create heatmap of top markers
  cat("Creating marker heatmap...\n")
  
  top_genes <- top_markers %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    pull(feats)
  
  plotMetaDataHeatmap(
    gobject = gobj,
    expression_values = "normalized",
    metadata_cols = cluster_column,
    selected_feats = top_genes,
    save_plot = TRUE,
    save_param = list(
      save_name = paste0(sample_id, "_marker_heatmap"),
      save_dir = results_folder,
      base_width = 14,
      base_height = 12
    )
  )
  
  cat("✓ Heatmap saved\n\n")
  
  # Violin plots for top markers
  cat("Creating violin plots for top 5 markers per cluster...\n")
  
  top5_per_cluster <- top_markers %>%
    group_by(cluster) %>%
    slice_head(n = 5) %>%
    pull(feats) %>%
    unique()
  
  violinPlot(
    gobject = gobj,
    feats = head(top5_per_cluster, 20),  # Limit to avoid overcrowding
    cluster_column = cluster_column,
    strip_text = 8,
    strip_position = "right",
    save_plot = TRUE,
    save_param = list(
      save_name = paste0(sample_id, "_marker_violins"),
      save_dir = results_folder,
      base_width = 16,
      base_height = 20
    )
  )
  
  cat("✓ Violin plots saved\n\n")
  
  cat("=== Marker Analysis Summary ===\n")
  cat("Total marker genes:", nrow(markers_scran), "\n")
  cat("Top markers per cluster:", top_n, "\n")
  cat("Clusters analyzed:", length(unique(markers_scran$cluster)), "\n")
  
  cat("\n✓ Marker analysis complete for", sample_id, "\n\n")
  
  return(gobj)
}

# Run if sourced directly
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) {
    sample_id <- args[1]
    input_path <- args[2]
    output_dir <- args[3]
    
    gobj <- marker_analysis(
      gobj = input_path,
      sample_id = sample_id,
      output_dir = output_dir
    )
    
    saveGiotto(gobj, dir = output_dir, foldername = "Giotto_Object_DEG_Markers", overwrite = TRUE)
  }
}

#!/usr/bin/env Rscript
# ==============================================================================
# 05_clustering.R
# Graph-based clustering
# ==============================================================================

#' Perform Clustering
#'
#' @param gobj Giotto object or path
#' @param sample_id Sample identifier
#' @param output_dir Output directory
#' @param k_nn Number of nearest neighbors
#' @param resolution Leiden resolution
#' @param dimensions_to_use PCA dimensions to use
#' @return Clustered Giotto object

perform_clustering <- function(gobj,
                               sample_id,
                               output_dir,
                               k_nn = 15,
                               resolution = 0.5,
                               dimensions_to_use = 1:30,
                               dim_reduction_to_use = "pca",
                               dim_reduction_name = NULL,
                               network_name = NULL) {
  
  cat("\n========================================\n")
  cat("STEP 05: Clustering\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }
  
  results_folder <- file.path(output_dir, "05_Clustering")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Create nearest neighbor network
  cat("Creating nearest neighbor network...\n")
  cat("  k:", k_nn, "\n")
  cat("  Dimensionality reduction:", dim_reduction_to_use, "\n")
  
  gobj <- createNearestNetwork(
    gobject = gobj,
    dim_reduction_to_use = dim_reduction_to_use,
    dim_reduction_name = dim_reduction_name,
    dimensions_to_use = dimensions_to_use,
    k = k_nn,
    name = network_name
  )
  
  cat("✓ NN network created\n\n")
  
  # Leiden clustering
  cat("Running Leiden clustering...\n")
  cat("  Resolution:", resolution, "\n")
  
  gobj <- doLeidenCluster(
    gobject = gobj,
    resolution = resolution,
    n_iterations = 1000,
    name = "leiden_clust"
  )
  
  n_clusters <- length(unique(pDataDT(gobj)$leiden_clust))
  cat("✓ Clustering complete:", n_clusters, "clusters\n\n")
  
  # Create clustering visualizations
  cat("Creating clustering plots...\n")
  
  # UMAP colored by cluster
  dimPlot2D(
    gobject = gobj,
    dim_reduction_to_use = "umap",
    cell_color = "leiden_clust",
    point_size = 0.8,
    show_legend = TRUE,
    save_plot = TRUE,
    save_param = list(
      save_name = paste0(sample_id, "_umap_clusters"),
      save_dir = results_folder,
      base_width = 11,
      base_height = 8
    )
  )
  
  # t-SNE colored by cluster
  dimPlot2D(
    gobject = gobj,
    dim_reduction_to_use = "tsne",
    cell_color = "leiden_clust",
    point_size = 0.8,
    show_legend = TRUE,
    save_plot = TRUE,
    save_param = list(
      save_name = paste0(sample_id, "_tsne_clusters"),
      save_dir = results_folder,
      base_width = 11,
      base_height = 8
    )
  )
  
  # Spatial plot by cluster
  spatPlot2D(
    gobject = gobj,
    cell_color = "leiden_clust",
    point_size = 0.5,
    show_image = FALSE,
    save_plot = TRUE,
    save_param = list(
      save_name = paste0(sample_id, "_spatial_clusters"),
      save_dir = results_folder,
      base_width = 14,
      base_height = 10
    )
  )
  
  cat("✓ Plots saved\n\n")
  
  # Cluster statistics
  cluster_stats <- pDataDT(gobj) %>%
    as_tibble() %>%
    group_by(leiden_clust) %>%
    summarise(
      n_cells = n(),
      mean_genes = mean(nr_feats),
      mean_counts = mean(total_expr),
      .groups = "drop"
    ) %>%
    arrange(leiden_clust)
  
  write_csv(cluster_stats, 
            file.path(results_folder, paste0(sample_id, "_cluster_stats.csv")))
  
  cat("=== Clustering Summary ===\n")
  cat("Number of clusters:", n_clusters, "\n")
  cat("Cluster sizes:\n")
  print(cluster_stats)
  
  cat("\n✓ Clustering complete for", sample_id, "\n\n")
  
  return(gobj)
}

# Run if sourced directly
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) {
    sample_id <- args[1]
    input_path <- args[2]
    output_dir <- args[3]
    
    gobj <- perform_clustering(
      gobj = input_path,
      sample_id = sample_id,
      output_dir = output_dir
    )
    
    saveGiotto(gobj, dir = output_dir, foldername = "Giotto_Object_Clustered", overwrite = TRUE)
  }
}

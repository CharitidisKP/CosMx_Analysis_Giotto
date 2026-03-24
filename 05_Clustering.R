# Graph-based clustering  -------------------------------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 05_Clustering.R
# ==============================================================================

#' Perform Clustering
#'
#' Builds a shared nearest-neighbour (sNN) graph on PCA or Harmony embeddings,
#' optionally inspects network topology via inspect_nn_network(), and
#' runs Leiden clustering. Saves UMAP, t-SNE, and spatial plots coloured
#' by cluster, plus a per-cluster statistics CSV.
#'
#' @param gobj                Giotto object or path to a saved Giotto object
#' @param sample_id           Sample identifier
#' @param output_dir          Root output directory
#' @param k_nn                Number of nearest neighbours for the sNN graph (default: 15)
#' @param resolution          Leiden resolution parameter (default: 0.5)
#' @param dimensions_to_use   PCA/Harmony dimensions to use (default: 1:30)
#' @param dim_reduction_to_use Which reduction to use for NN graph (default: "pca")
#' @param dim_reduction_name  Name of the reduction object (default: NULL)
#' @param network_name        Name to assign to the resulting NN network (default: NULL)
#' @param scripts_dir         Path to the Scripts directory; used to source
#'                            Helper_Scripts/Inspect_sNN_network.R.
#'                            Defaults to the directory of the current script.
#' @param inspect_snn         Run sNN inspection after network construction?
#'                            (default: TRUE)
#' @param n_cells_subgraph    Maximum cells to subsample for the subgraph plot
#'                            inside inspect_nn_network() (default: 100000)
#' @return Clustered Giotto object

perform_clustering <- function(gobj,
                               sample_id,
                               output_dir,
                               k_nn                 = 15,
                               resolution           = 0.5,
                               dimensions_to_use    = 1:30,
                               dim_reduction_to_use = "pca",
                               dim_reduction_name   = NULL,
                               network_name         = NULL,
                               scripts_dir          = NULL,
                               inspect_snn          = TRUE,
                               n_cells_subgraph     = 100000) {
  
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
  
  # ── Create nearest neighbour network ────────────────────────────────────────
  cat("Creating nearest neighbor network...\n")
  cat("  k:", k_nn, "\n")
  cat("  Dimensionality reduction:", dim_reduction_to_use, "\n")
  
  gobj <- createNearestNetwork(
    gobject              = gobj,
    dim_reduction_to_use = dim_reduction_to_use,
    dim_reduction_name   = dim_reduction_name,
    dimensions_to_use    = dimensions_to_use,
    k                    = k_nn,
    name                 = network_name
  )
  
  cat("✓ NN network created\n\n")
  
  # ── sNN network inspection ───────────────────────────────────────────────────
  if (inspect_snn) {
    
    # Resolve scripts_dir if not supplied
    if (is.null(scripts_dir)) {
      scripts_dir <- tryCatch(
        dirname(sys.frame(1)$ofile),
        error = function(e) getwd()
      )
    }
    
    scripts_dir <- normalizePath(scripts_dir, winslash = "/", mustWork = FALSE)
    inspect_candidates <- c(
      file.path(scripts_dir, "Helper_Scripts", "Inspect_sNN_network.R"),
      file.path(scripts_dir, "Scripts", "Helper_Scripts", "Inspect_sNN_network.R"),
      file.path(dirname(scripts_dir), "Helper_Scripts", "Inspect_sNN_network.R")
    )
    inspect_candidates <- unique(inspect_candidates)
    inspect_script <- inspect_candidates[file.exists(inspect_candidates)][1]
    
    if (!is.na(inspect_script) && nzchar(inspect_script) && file.exists(inspect_script)) {
      source(inspect_script)
      
      cat("Running sNN network inspection...\n")
      
      tryCatch({
        nn_results <- inspect_nn_network(
          gobject          = gobj,
          nn_name          = if (!is.null(network_name)) network_name else "sNN.pca",
          output_dir       = results_folder,
          sample_id        = sample_id,
          subfolder_name   = "sNN_inspection",
          n_cells_subgraph = n_cells_subgraph,
          save_plots       = TRUE,
          verbose          = TRUE
        )
        
        cat("\u2713 sNN inspection complete\n")
        cat("  Mean degree:         ", round(nn_results$mean_degree, 2), "\n")
        cat("  Connected components:", nn_results$n_components, "\n")
        cat("  Largest component:   ", nn_results$largest_component_size,
            "cells\n\n")
        
      }, error = function(e) {
        cat("\u26A0 sNN inspection failed:", conditionMessage(e),
            "\n  Continuing with clustering.\n\n")
      })
      
    } else {
      cat("\u26A0 inspect_nn_network.R not found. Tried:\n  ",
          paste(inspect_candidates, collapse = "\n  "),
          "\n  Skipping sNN inspection.\n\n")
    }
  }
  
  # ── Leiden clustering ────────────────────────────────────────────────────────
  cat("Running Leiden clustering...\n")
  cat("  Resolution:", resolution, "\n")
  
  gobj <- doLeidenCluster(
    gobject      = gobj,
    resolution   = resolution,
    n_iterations = 1000,
    name         = "leiden_clust"
  )
  
  n_clusters <- length(unique(pDataDT(gobj)$leiden_clust))
  cat("✓ Clustering complete:", n_clusters, "clusters\n\n")
  
  # ── Clustering visualisations ────────────────────────────────────────────────
  cat("Creating clustering plots...\n")
  
  # UMAP coloured by cluster
  dimPlot2D(
    gobject              = gobj,
    dim_reduction_to_use = "umap",
    cell_color           = "leiden_clust",
    point_size           = 0.8,
    show_legend          = TRUE,
    save_plot            = TRUE,
    save_param           = list(
      save_name   = paste0(sample_id, "_umap_clusters"),
      save_dir    = results_folder,
      base_width  = 11,
      base_height = 8
    )
  )
  
  # t-SNE coloured by cluster
  dimPlot2D(
    gobject              = gobj,
    dim_reduction_to_use = "tsne",
    cell_color           = "leiden_clust",
    point_size           = 0.8,
    show_legend          = TRUE,
    save_plot            = TRUE,
    save_param           = list(
      save_name   = paste0(sample_id, "_tsne_clusters"),
      save_dir    = results_folder,
      base_width  = 11,
      base_height = 8
    )
  )
  
  # Spatial plot by cluster
  spatPlot2D(
    gobject    = gobj,
    cell_color = "leiden_clust",
    point_size = 0.5,
    show_image = FALSE,
    save_plot  = TRUE,
    save_param = list(
      save_name   = paste0(sample_id, "_spatial_clusters"),
      save_dir    = results_folder,
      base_width  = 14,
      base_height = 10
    )
  )
  
  cat("✓ Plots saved\n\n")
  
  # ── Cluster statistics ───────────────────────────────────────────────────────
  cluster_stats <- pDataDT(gobj) %>%
    as_tibble() %>%
    dplyr::group_by(leiden_clust) %>%
    dplyr::summarise(
      n_cells     = dplyr::n(),
      mean_genes  = mean(nr_feats),
      mean_counts = mean(total_expr),
      .groups     = "drop"
    ) %>%
    dplyr::arrange(leiden_clust)
  
  readr::write_csv(
    cluster_stats,
    file.path(results_folder, paste0(sample_id, "_cluster_stats.csv"))
  )
  
  cat("=== Clustering Summary ===\n")
  cat("Number of clusters:", n_clusters, "\n")
  cat("Cluster sizes:\n")
  print(cluster_stats)
  
  cat("\n\u2713 STEP 05 complete for", sample_id, "\n\n")
  
  return(gobj)
}


# Run if sourced directly ---------------------------------------------------
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) {
    script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
    script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
    bootstrap_script <- file.path(script_dir, "Helper_Scripts", "Script_Bootstrap.R")
    if (file.exists(bootstrap_script)) {
      source(bootstrap_script, local = .GlobalEnv)
      bootstrap_pipeline_environment(script_dir, load_pipeline_utils = FALSE, verbose = FALSE)
    }
    
    sample_id  <- args[1]
    input_path <- args[2]
    output_dir <- args[3]
    
    gobj <- perform_clustering(
      gobj       = input_path,
      sample_id  = sample_id,
      output_dir = output_dir
    )
    
    saveGiotto(gobj,
               dir        = output_dir,
               foldername = "Giotto_Object_Clustered",
               overwrite  = TRUE)
  } else {
    stop("Usage: Rscript 05_Clustering.R <sample_id> <input_path> <output_dir>")
  }
}

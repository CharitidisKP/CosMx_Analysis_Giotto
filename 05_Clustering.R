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
#' @param leiden_n_iterations Number of Leiden iterations for the primary
#'                            clustering run (default: 200)
#' @param resolution_sweep    Optional numeric vector of extra Leiden
#'                            resolutions to save as additional metadata columns
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
                               python_path          = NULL,
                               scripts_dir          = NULL,
                               inspect_snn          = TRUE,
                               n_cells_subgraph     = 100000,
                               leiden_n_iterations  = 200,
                               resolution_sweep     = NULL) {
  
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
  
  if (is.null(python_path) || !nzchar(python_path)) {
    python_path <- Sys.getenv("COSMX_PYTHON_PATH", unset = "")
  }
  if ((!nzchar(python_path)) && exists("python_path", envir = .GlobalEnv, inherits = FALSE)) {
    python_path <- get("python_path", envir = .GlobalEnv)
  }
  
  if (nzchar(python_path) && requireNamespace("reticulate", quietly = TRUE)) {
    if (reticulate::py_available(initialize = TRUE)) {
      py_cfg <- reticulate::py_config()
      active_python <- normalizePath(py_cfg$python, winslash = "/", mustWork = FALSE)
      requested_python <- normalizePath(python_path, winslash = "/", mustWork = FALSE)
      
      if (!identical(active_python, requested_python)) {
        stop(
          "reticulate is already initialized with ", active_python,
          " but clustering is configured to use ", requested_python,
          ". Restart the R session before running step 05."
        )
      }
    } else {
      reticulate::use_python(python_path, required = TRUE)
    }
    
    missing_modules <- c("igraph", "leidenalg")[!vapply(
      c("igraph", "leidenalg"),
      reticulate::py_module_available,
      logical(1)
    )]
    
    if (length(missing_modules) > 0) {
      stop(
        "Leiden clustering requires Python module(s) ",
        paste(missing_modules, collapse = ", "),
        " in ", python_path, ". Install them in that environment and restart R."
      )
    }
  }
  
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
        mean_degree <- nn_results$mean_degree
        if (is.null(mean_degree) && !is.null(nn_results$degree_distribution)) {
          mean_degree <- mean(nn_results$degree_distribution)
        }
        n_components <- nn_results$n_components
        if (is.null(n_components) && !is.null(nn_results$network_metrics)) {
          n_components <- nn_results$network_metrics$Value[
            nn_results$network_metrics$Metric == "Number of Components"
          ][1]
        }
        largest_component_size <- nn_results$largest_component_size
        if (is.null(largest_component_size) && !is.null(nn_results$network_metrics)) {
          largest_component_size <- nn_results$network_metrics$Value[
            nn_results$network_metrics$Metric == "Largest Component Size"
          ][1]
        }
        cat("  Mean degree:         ", round(as.numeric(mean_degree), 2), "\n")
        cat("  Connected components:", as.character(n_components), "\n")
        cat("  Largest component:   ", as.character(largest_component_size),
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
  cat("  Iterations:", leiden_n_iterations, "\n")
  
  gobj <- doLeidenCluster(
    gobject      = gobj,
    python_path  = if (nzchar(python_path)) python_path else NULL,
    resolution   = resolution,
    n_iterations = leiden_n_iterations,
    name         = "leiden_clust"
  )
  
  n_clusters <- length(unique(pDataDT(gobj)$leiden_clust))
  cat("✓ Clustering complete:", n_clusters, "clusters\n\n")
  
  resolution_sweep <- unique(as.numeric(resolution_sweep))
  resolution_sweep <- resolution_sweep[is.finite(resolution_sweep) & resolution_sweep > 0]
  resolution_sweep <- setdiff(resolution_sweep, resolution)
  sweep_summary <- list()
  
  if (length(resolution_sweep) > 0) {
    cat("Running auxiliary Leiden resolution sweep...\n")
    for (res in resolution_sweep) {
      sweep_name <- paste0(
        "leiden_clust_res_",
        gsub("(^_+|_+$)", "", gsub("[^A-Za-z0-9]+", "_", format(res, trim = TRUE, scientific = FALSE)))
      )
      gobj <- doLeidenCluster(
        gobject = gobj,
        python_path = if (nzchar(python_path)) python_path else NULL,
        resolution = res,
        n_iterations = leiden_n_iterations,
        name = sweep_name
      )
      sweep_summary[[length(sweep_summary) + 1L]] <- data.frame(
        cluster_column = sweep_name,
        resolution = res,
        n_clusters = length(unique(pDataDT(gobj)[[sweep_name]])),
        stringsAsFactors = FALSE
      )
      cat("  ✓ Saved", sweep_name, "(", tail(sweep_summary, 1)[[1]]$n_clusters, "clusters)\n")
    }
    readr::write_csv(
      dplyr::bind_rows(sweep_summary),
      file.path(results_folder, paste0(sample_id, "_resolution_sweep_summary.csv"))
    )
    cat("\n")
  }
  
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
  
  cluster_script_base <- if (is.null(scripts_dir)) getwd() else scripts_dir
  cluster_vis_candidates <- c(
    file.path(cluster_script_base, "Helper_Scripts", "Cluster_Visualisations.R"),
    file.path(cluster_script_base, "Scripts", "Helper_Scripts", "Cluster_Visualisations.R"),
    file.path(dirname(cluster_script_base), "Helper_Scripts", "Cluster_Visualisations.R")
  )
  cluster_vis_candidates <- unique(cluster_vis_candidates)
  cluster_vis_script <- cluster_vis_candidates[file.exists(cluster_vis_candidates)][1]
  
  if (!is.na(cluster_vis_script) && nzchar(cluster_vis_script) && file.exists(cluster_vis_script)) {
    source(cluster_vis_script)
    tryCatch({
      create_clustering_visualization(
        gobject = gobj,
        cluster_column = "leiden_clust",
        title_suffix = paste0(" - ", sample_id, " Leiden Clusters"),
        save_plots = TRUE,
        save_dir = results_folder,
        prefix = paste0(sample_id, "_custom_clusters")
      )
      cat("✓ Custom presentation clustering plots saved\n\n")
    }, error = function(e) {
      cat("⚠ Custom clustering visualizations failed:", conditionMessage(e), "\n\n")
    })
  }
  
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

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
                               resolution_sweep     = NULL,
                               seed                 = 42) {
  
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
  cat("  Seed:", seed, "\n")

  set.seed(seed)
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
      set.seed(seed)
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
    sweep_df <- dplyr::bind_rows(sweep_summary)
    # Include the primary run for visual context
    sweep_df <- dplyr::bind_rows(
      sweep_df,
      data.frame(
        cluster_column = "leiden_clust",
        resolution     = resolution,
        n_clusters     = n_clusters,
        stringsAsFactors = FALSE
      )
    )
    sweep_df <- sweep_df[order(sweep_df$resolution), ]
    readr::write_csv(
      sweep_df,
      file.path(results_folder, paste0(sample_id, "_resolution_sweep_summary.csv"))
    )
    tryCatch({
      sweep_df$is_primary <- sweep_df$resolution == resolution
      p_sweep <- ggplot2::ggplot(
          sweep_df,
          ggplot2::aes(x = resolution, y = n_clusters)
        ) +
        ggplot2::geom_line(colour = "grey40", linewidth = 0.5) +
        ggplot2::geom_point(
          ggplot2::aes(colour = is_primary, size = is_primary)
        ) +
        ggplot2::scale_colour_manual(
          values = c(`TRUE` = "#E41A1C", `FALSE` = "grey20"),
          guide = "none"
        ) +
        ggplot2::scale_size_manual(
          values = c(`TRUE` = 3.5, `FALSE` = 2.2),
          guide = "none"
        ) +
        ggplot2::scale_x_continuous(breaks = sweep_df$resolution) +
        ggplot2::labs(
          title    = sample_plot_title(sample_id, "Leiden resolution sweep"),
          subtitle = "Red point = primary resolution",
          x = "Resolution", y = "Number of clusters"
        ) +
        presentation_theme(base_size = 12)
      save_presentation_plot(
        plot     = p_sweep,
        filename = file.path(results_folder,
                             paste0(sample_id, "_resolution_sweep.png")),
        width    = 8, height = 6, dpi = 300
      )
    }, error = function(e) {
      cat("\u26A0 Resolution-sweep plot failed:", conditionMessage(e), "\n")
    })
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

  # Cluster-size bar chart
  tryCatch({
    p_size <- ggplot2::ggplot(
        cluster_stats,
        ggplot2::aes(x = reorder(factor(leiden_clust), -n_cells), y = n_cells)
      ) +
      ggplot2::geom_col(fill = "#4C72B0") +
      ggplot2::geom_text(
        ggplot2::aes(label = n_cells),
        vjust = -0.3, size = 3
      ) +
      ggplot2::labs(
        title = sample_plot_title(sample_id, "Cells per Leiden cluster"),
        x = "Cluster (ordered by size)", y = "Number of cells"
      ) +
      presentation_theme(base_size = 12)
    save_presentation_plot(
      plot     = p_size,
      filename = file.path(results_folder,
                           paste0(sample_id, "_cluster_sizes.png")),
      width    = max(8, 0.5 * n_clusters + 4),
      height   = 6, dpi = 300
    )
  }, error = function(e) {
    cat("\u26A0 Cluster-size bar plot failed:", conditionMessage(e), "\n")
  })

  # Silhouette score per cluster (subsampled for tractability on >10k cells)
  tryCatch({
    if (requireNamespace("cluster", quietly = TRUE)) {
      pca_mat <- tryCatch(
        getDimReduction(
          gobject = gobj, spat_unit = "cell", feat_type = "rna",
          reduction = "cells", reduction_method = "pca",
          name = "pca", output = "matrix"
        ),
        error = function(e) NULL
      )
      if (!is.null(pca_mat) && nrow(pca_mat) > 0) {
        meta <- as.data.frame(pDataDT(gobj))
        clust <- as.integer(as.factor(meta$leiden_clust))
        names(clust) <- meta$cell_ID
        rn <- rownames(pca_mat)
        common <- intersect(rn, names(clust))
        pca_mat <- pca_mat[common, seq_len(min(ncol(pca_mat), 20)), drop = FALSE]
        clust <- clust[common]
        # Subsample per cluster up to 400 cells (keeps dist matrix tractable)
        set.seed(seed)
        subs <- unlist(lapply(unique(clust), function(cl) {
          ids <- which(clust == cl)
          if (length(ids) <= 400) ids else sample(ids, 400)
        }))
        pca_sub <- pca_mat[subs, , drop = FALSE]
        clust_sub <- clust[subs]
        if (length(unique(clust_sub)) >= 2 && nrow(pca_sub) >= 50) {
          sil <- cluster::silhouette(clust_sub, stats::dist(pca_sub))
          sil_df <- data.frame(
            cluster = factor(sil[, "cluster"]),
            sil_width = sil[, "sil_width"]
          )
          sil_summary <- dplyr::summarise(
            dplyr::group_by(sil_df, cluster),
            mean_sil = mean(sil_width),
            .groups = "drop"
          )
          p_sil <- ggplot2::ggplot(
              sil_summary,
              ggplot2::aes(x = reorder(cluster, -mean_sil), y = mean_sil,
                           fill = mean_sil > 0)
            ) +
            ggplot2::geom_col() +
            ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
            ggplot2::scale_fill_manual(
              values = c(`TRUE` = "#4C72B0", `FALSE` = "#C44E52"),
              guide = "none"
            ) +
            ggplot2::labs(
              title    = sample_plot_title(sample_id, "Mean silhouette per cluster"),
              subtitle = sprintf("PCA space (first %d PCs); subsampled to %d cells",
                                 ncol(pca_sub), nrow(pca_sub)),
              x = "Cluster", y = "Mean silhouette width"
            ) +
            presentation_theme(base_size = 12)
          save_presentation_plot(
            plot     = p_sil,
            filename = file.path(results_folder,
                                 paste0(sample_id, "_cluster_silhouette.png")),
            width    = max(8, 0.5 * n_clusters + 4),
            height   = 6, dpi = 300
          )
        }
      }
    } else {
      cat("  note: 'cluster' package not available; skipping silhouette plot\n")
    }
  }, error = function(e) {
    cat("\u26A0 Silhouette plot failed:", conditionMessage(e), "\n")
  })

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

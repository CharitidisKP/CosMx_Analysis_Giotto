# Calculate HVGs, PCA, UMAP, and tSNE -------------------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 04_dimensionality_reduction.R
# PCA, UMAP, and t-SNE dimensionality reduction
# Includes support for spatially-informed HVG selection
# ==============================================================================

#' Calculate row variances for sparse matrix
#' @param x Sparse matrix
#' @return Vector of row variances

sparse_row_var <- function(x) {
  row_means <- Matrix::rowMeans(x)
  n <- ncol(x)
  row_sq_diffs <- Matrix::rowSums((x - row_means)^2)
  row_vars <- row_sq_diffs / (n - 1)
  return(row_vars)
}

#' Dimensionality Reduction
#'
#' @param gobj Giotto object or path
#' @param sample_id Sample identifier
#' @param output_dir Output directory
#' @param n_hvgs Number of highly variable genes
#' @param n_pcs Number of principal components
#' @param umap_n_neighbors UMAP neighbors parameter
#' @param umap_min_dist UMAP minimum distance
#' @param spatial_hvg Use spatially-informed HVG selection
#' @return Giotto object with dimensionality reductions

dimensionality_reduction <- function(gobj,
                                     sample_id,
                                     output_dir,
                                     n_hvgs = 500,
                                     n_pcs = 30,
                                     umap_n_neighbors = 30,
                                     umap_min_dist = 0.3,
                                     spatial_hvg = FALSE) {
  
  cat("\n========================================\n")
  cat("STEP 04: Dimensionality Reduction\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }
  
  results_folder <- file.path(output_dir, "04_Dimensionality_Reduction")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Identify highly variable genes
  if (spatial_hvg) {
    cat("Identifying spatially variable genes...\n")
    cat("  Target SVGs:", n_hvgs, "\n")
    
    # Use spatial method (binSpect or spark)
    gobj <- runSpatialEnrich(
      gobject = gobj,
      sign_method = "rank",
      expression_values = "normalized",
      min_cells_per_grid = 4
    )
    
    # Get top spatial genes
    spatial_genes <- getSpatialEnrichment(
      gobj,
      output = "data.table"
    ) %>%
      as_tibble() %>%
      arrange(desc(rank_score)) %>%
      slice_head(n = n_hvgs) %>%
      pull(feats)
    
    hvg_genes <- spatial_genes
    
    cat("✓ Selected", length(hvg_genes), "spatially variable genes\n\n")
    
  } else {
    cat("Identifying highly variable genes...\n")
    cat("  Target HVGs:", n_hvgs, "\n")
    
    gobj <- calculateHVF(
      gobject = gobj,
      method = "cov_loess",
      reverse_log_scale = TRUE,
      difference_in_cov = 0.1,
      show_plot = FALSE,
      save_plot = TRUE,
      save_param = list(
        save_name = paste0(sample_id, "_hvg_plot"),
        save_dir = results_folder,
        base_width = 10,
        base_height = 8
      )
    )
    
    # Get HVG genes
    feat_metadata <- getFeatureMetadata(gobj, output = "data.table") %>%
      as_tibble()
    
    cat("Feature metadata columns:", paste(names(feat_metadata), collapse = ", "), "\n")
    
    # Find HVG column
    hvf_col <- if ("hvf" %in% names(feat_metadata)) {
      "hvf"
    } else if ("hvg" %in% names(feat_metadata)) {
      "hvg"
    } else {
      stop("Cannot find HVF/HVG column in feature metadata")
    }
    
    # Get HVG genes
    hvg_genes <- feat_metadata %>%
      filter(.data[[hvf_col]] == "yes") %>%
      pull(feat_ID)
    
    # If we have more HVGs than requested, select top N by variance
    if (length(hvg_genes) > n_hvgs) {
      cat("Found", length(hvg_genes), "HVGs, selecting top", n_hvgs, "by variance\n")
      
      # Calculate variance for HVG genes
      expr_mat <- getExpression(gobj, values = "normalized", output = "matrix")
      hvg_subset <- expr_mat[hvg_genes, ]
      
      # Use custom sparse row variance function
      hvg_vars <- sparse_row_var(hvg_subset)
      names(hvg_vars) <- hvg_genes
      
      hvg_genes <- names(sort(hvg_vars, decreasing = TRUE)[1:n_hvgs])
    }
    
    cat("✓ Selected", length(hvg_genes), "highly variable genes\n\n")
  }
  
  # Save HVG list
  write_lines(hvg_genes, 
              file.path(results_folder, paste0(sample_id, "_hvg_genes.txt")))
  
  # PCA
  cat("Running PCA...\n")
  cat("  Number of PCs:", n_pcs, "\n")
  cat("  Using", length(hvg_genes), "features\n")
  
  gobj <- runPCA(
    gobject = gobj,
    feats_to_use = hvg_genes,
    scale_unit = TRUE,
    center = TRUE,
    ncp = n_pcs
  )
  
  cat("✓ PCA complete\n\n")
  
  # Scree plot
  tryCatch({
    screePlot(
      gobject = gobj,
      ncp = min(30, n_pcs),
      save_plot = TRUE,
      save_param = list(
        save_name = paste0(sample_id, "_scree_plot"),
        save_dir = results_folder,
        base_width = 10,
        base_height = 6
      )
    )
    cat("✓ Scree plot saved\n\n")
  }, error = function(e) {
    cat("⚠ Scree plot warning:", conditionMessage(e), "\n\n")
  })
  
  # UMAP
  cat("Running UMAP...\n")
  cat("  Neighbors:", umap_n_neighbors, "\n")
  cat("  Min distance:", umap_min_dist, "\n")
  
  gobj <- runUMAP(
    gobject = gobj,
    dimensions_to_use = 1:n_pcs,
    n_neighbors = umap_n_neighbors,
    min_dist = umap_min_dist
  )
  
  cat("✓ UMAP complete\n\n")
  
  # t-SNE
  cat("Running t-SNE...\n")
  
  gobj <- runtSNE(
    gobject = gobj,
    dimensions_to_use = 1:n_pcs,
    perplexity = 30
  )
  
  cat("✓ t-SNE complete\n\n")
  
  # Create dimensionality reduction plots
  cat("Creating dimensionality reduction plots...\n")
  
  tryCatch({
    # UMAP plots
    dimPlot2D(
      gobject = gobj,
      dim_reduction_to_use = "umap",
      cell_color = "nr_feats",
      point_size = 0.5,
      save_plot = TRUE,
      save_param = list(
        save_name = paste0(sample_id, "_umap_by_genes"),
        save_dir = results_folder,
        base_width = 10,
        base_height = 8
      )
    )
    
    dimPlot2D(
      gobject = gobj,
      dim_reduction_to_use = "umap",
      cell_color = "total_expr",
      point_size = 0.5,
      save_plot = TRUE,
      save_param = list(
        save_name = paste0(sample_id, "_umap_by_counts"),
        save_dir = results_folder,
        base_width = 10,
        base_height = 8
      )
    )
    
    # t-SNE plots
    dimPlot2D(
      gobject = gobj,
      dim_reduction_to_use = "tsne",
      cell_color = "nr_feats",
      point_size = 0.5,
      save_plot = TRUE,
      save_param = list(
        save_name = paste0(sample_id, "_tsne_by_genes"),
        save_dir = results_folder,
        base_width = 10,
        base_height = 8
      )
    )
    
    cat("✓ Plots saved\n\n")
  }, error = function(e) {
    cat("⚠ Plotting warning:", conditionMessage(e), "\n\n")
  })
  
  # Summary
  cat("=== Dimensionality Reduction Summary ===\n")
  cat("Method:", ifelse(spatial_hvg, "Spatial HVG", "Standard HVG"), "\n")
  cat("Features selected:", length(hvg_genes), "\n")
  cat("PCs computed:", n_pcs, "\n")
  
  tryCatch({
    umap_dims <- getDimReduction(gobj, "umap")
    tsne_dims <- getDimReduction(gobj, "tsne")
    cat("UMAP dimensions:", ncol(umap_dims), "\n")
    cat("t-SNE dimensions:", ncol(tsne_dims), "\n")
  }, error = function(e) {
    cat("Dimension retrieval: info unavailable\n")
  })
  
  cat("\n✓ Dimensionality reduction complete for", sample_id, "\n\n")
  
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
      bootstrap_pipeline_environment(script_dir, load_pipeline_utils = FALSE, verbose = FALSE)
    }
    
    sample_id <- args[1]
    input_path <- args[2]
    output_dir <- args[3]
    spatial_hvg <- if (length(args) >= 4) as.logical(args[4]) else FALSE
    
    gobj <- dimensionality_reduction(
      gobj = input_path,
      sample_id = sample_id,
      output_dir = output_dir,
      spatial_hvg = spatial_hvg
    )
    
    saveGiotto(gobj, dir = output_dir, foldername = "Giotto_Object_DimReduced", overwrite = TRUE)
  } else {
    stop("Usage: Rscript 04_Dimensionallity_Reduction.R <sample_id> <input_path> <output_dir> [spatial_hvg]")
  }
}

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
  if (!inherits(x, "sparseMatrix")) {
    return(apply(x, 1, stats::var))
  }
  
  x <- methods::as(x, "dgCMatrix")
  n <- ncol(x)
  if (n < 2) {
    return(rep(0, nrow(x)))
  }
  
  row_means <- Matrix::rowMeans(x)
  x_sq <- x
  x_sq@x <- x_sq@x^2
  row_sq_means <- Matrix::rowSums(x_sq) / n
  pmax(0, (n / (n - 1)) * (row_sq_means - (row_means^2)))
}

.run_known_giotto_warning_safe <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Not all expression matrices share the same cell_IDs", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

.select_sparse_hvgs <- function(gobj, n_hvgs, sample_id, results_folder) {
  expr_mat <- getExpression(gobj, values = "normalized", output = "matrix")
  gene_means <- Matrix::rowMeans(expr_mat)
  gene_vars <- sparse_row_var(expr_mat)
  names(gene_means) <- rownames(expr_mat)
  names(gene_vars) <- rownames(expr_mat)
  
  if (length(gene_vars) == 0) {
    stop("No genes were available for HVG selection")
  }
  
  gene_vars[is.na(gene_vars)] <- 0
  gene_means[is.na(gene_means)] <- 0
  ranked_genes <- names(sort(gene_vars, decreasing = TRUE))
  hvg_genes <- head(ranked_genes, min(n_hvgs, length(ranked_genes)))
  
  hvg_table <- tibble::tibble(
    feat_ID = names(gene_means),
    mean_expr = as.numeric(gene_means),
    variance = as.numeric(gene_vars)
  ) %>%
    dplyr::mutate(selected_hvg = feat_ID %in% hvg_genes) %>%
    dplyr::arrange(dplyr::desc(variance))
  
  readr::write_csv(
    hvg_table,
    file.path(results_folder, paste0(sample_id, "_hvg_summary.csv"))
  )
  
  hvg_plot <- ggplot(
    hvg_table,
    aes(
      x = log10(pmax(mean_expr, 1e-8)),
      y = log10(pmax(variance, 1e-8)),
      color = selected_hvg
    )
  ) +
    geom_point(alpha = 0.5, size = 0.7) +
    scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "firebrick")) +
    labs(
      title = paste(sample_id, "- Highly variable genes"),
      subtitle = "Top genes selected by sparse variance on normalized expression",
      x = "Log10 mean expression",
      y = "Log10 variance",
      color = "Selected"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(
    filename = file.path(results_folder, paste0(sample_id, "_hvg_plot.png")),
    plot = hvg_plot,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  hvg_genes
}

.prepare_dim_plot_data <- function(gobj, dim_reduction_name, color_column) {
  dims <- as.data.frame(getDimReduction(gobj, dim_reduction_name))
  if (ncol(dims) < 2) {
    stop("Dimensionality reduction ", dim_reduction_name, " does not contain at least two components")
  }
  
  colnames(dims)[1:2] <- c("dim1", "dim2")
  dims$cell_ID <- if (!is.null(rownames(dims)) && any(nzchar(rownames(dims)))) {
    rownames(dims)
  } else {
    pDataDT(gobj)$cell_ID[seq_len(nrow(dims))]
  }
  
  metadata <- pDataDT(gobj) %>%
    as_tibble() %>%
    dplyr::select(cell_ID, dplyr::all_of(color_column))
  
  dplyr::left_join(dims, metadata, by = "cell_ID")
}

.save_continuous_dim_plot <- function(plot_data,
                                      color_column,
                                      reduction_label,
                                      value_label,
                                      output_path,
                                      sample_id) {
  p <- ggplot(plot_data, aes(x = dim1, y = dim2, color = .data[[color_column]])) +
    geom_point(size = 0.35, alpha = 0.7) +
    scale_color_viridis_c(option = "magma", na.value = "grey80") +
    labs(
      title = paste(sample_id, "-", reduction_label, "colored by", value_label),
      x = paste(reduction_label, "1"),
      y = paste(reduction_label, "2"),
      color = value_label
    ) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(
    filename = output_path,
    plot = p,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
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
    
    hvg_genes <- .select_sparse_hvgs(
      gobj = gobj,
      n_hvgs = n_hvgs,
      sample_id = sample_id,
      results_folder = results_folder
    )
    
    cat("✓ Selected", length(hvg_genes), "highly variable genes\n\n")
  }
  
  # Save HVG list
  readr::write_lines(
    hvg_genes,
    file.path(results_folder, paste0(sample_id, "_hvg_genes.txt"))
  )
  
  # PCA
  cat("Running PCA...\n")
  cat("  Number of PCs:", n_pcs, "\n")
  cat("  Using", length(hvg_genes), "features\n")
  
  gobj <- .run_known_giotto_warning_safe(
    runPCA(
      gobject = gobj,
      feats_to_use = hvg_genes,
      scale_unit = TRUE,
      center = TRUE,
      ncp = n_pcs
    )
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
  
  gobj <- .run_known_giotto_warning_safe(
    runUMAP(
      gobject = gobj,
      dimensions_to_use = 1:n_pcs,
      n_neighbors = umap_n_neighbors,
      min_dist = umap_min_dist
    )
  )
  
  cat("✓ UMAP complete\n\n")
  
  # t-SNE
  cat("Running t-SNE...\n")
  
  gobj <- .run_known_giotto_warning_safe(
    runtSNE(
      gobject = gobj,
      dimensions_to_use = 1:n_pcs,
      perplexity = 30
    )
  )
  
  cat("✓ t-SNE complete\n\n")
  
  # Create dimensionality reduction plots
  cat("Creating dimensionality reduction plots...\n")
  
  tryCatch({
    umap_genes <- .prepare_dim_plot_data(gobj, "umap", "nr_feats")
    .save_continuous_dim_plot(
      plot_data = umap_genes,
      color_column = "nr_feats",
      reduction_label = "UMAP",
      value_label = "genes per cell",
      output_path = file.path(results_folder, paste0(sample_id, "_umap_by_genes.png")),
      sample_id = sample_id
    )
    
    umap_counts <- .prepare_dim_plot_data(gobj, "umap", "total_expr")
    .save_continuous_dim_plot(
      plot_data = umap_counts,
      color_column = "total_expr",
      reduction_label = "UMAP",
      value_label = "total counts",
      output_path = file.path(results_folder, paste0(sample_id, "_umap_by_counts.png")),
      sample_id = sample_id
    )
    
    tsne_genes <- .prepare_dim_plot_data(gobj, "tsne", "nr_feats")
    .save_continuous_dim_plot(
      plot_data = tsne_genes,
      color_column = "nr_feats",
      reduction_label = "t-SNE",
      value_label = "genes per cell",
      output_path = file.path(results_folder, paste0(sample_id, "_tsne_by_genes.png")),
      sample_id = sample_id
    )
    
    cat("✓ Plots saved\n\n")
  }, error = function(e) {
    cat("⚠ Plotting warning:", conditionMessage(e), "\n\n")
  })
  
  # Summary
  cat("=== Dimensionality Reduction Summary ===\n")
  cat("Method:", ifelse(spatial_hvg, "Spatial HVG", "Standard HVG (sparse variance)"), "\n")
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

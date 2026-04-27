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

current_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  ofiles <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else frame$ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles) > 0) {
    return(dirname(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

pipeline_utils <- file.path(current_script_dir(), "Helper_Scripts", "Pipeline_Utils.R")
if ((!exists("presentation_theme") || !exists("sample_plot_title") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
  source(pipeline_utils)
}

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
      title = sample_plot_title(sample_id, "Highly Variable Genes"),
      subtitle = "Top genes selected by sparse variance on normalized expression",
      x = log10_axis_label("Mean Expression"),
      y = log10_axis_label("Variance"),
      color = "Selected\nas HVG"
    ) +
    presentation_theme(base_size = 12)
  
  save_presentation_plot(
    plot = hvg_plot,
    filename = file.path(results_folder, paste0(sample_id, "_hvg_plot.png")),
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  hvg_genes
}

.prepare_dim_plot_data <- function(gobj, dim_reduction_name, color_column) {
  dims <- as.data.frame(
    getDimReduction(
      gobject = gobj,
      spat_unit = "cell",
      feat_type = "rna",
      reduction = "cells",
      reduction_method = dim_reduction_name,
      name = dim_reduction_name,
      output = "matrix"
    )
  )
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
  reduction_label <- display_reduction_name(reduction_label)
  p <- ggplot(plot_data, aes(x = dim1, y = dim2, color = .data[[color_column]])) +
    geom_point(size = 0.35, alpha = 0.7) +
    scale_color_viridis_c(option = "magma", na.value = "grey80") +
    labs(
      title = sample_plot_title(sample_id, paste(reduction_label, "Embedding Colored By", tools::toTitleCase(value_label))),
      x = embedding_axis_label(reduction_label, 1),
      y = embedding_axis_label(reduction_label, 2),
      color = value_label
    ) +
    presentation_theme(base_size = 12)
  
  save_presentation_plot(
    plot = p,
    filename = output_path,
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
                                     spatial_hvg = FALSE,
                                     seed = 42,
                                     sample_row = NULL,
                                     sample_sheet_path = NULL) {

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

  # Dynamic n_hvgs: when NULL or "auto", use min(500, max(100, round(0.1 * n_genes)))
  # to avoid selecting 50% of a 1K panel as "variable". Explicit numeric values
  # are honoured verbatim (still capped at n_genes downstream).
  n_genes_total <- length(gobj@feat_ID$rna)
  n_hvgs_source <- "explicit"
  if (is.null(n_hvgs) ||
      (is.character(n_hvgs) && tolower(n_hvgs) %in% c("auto", "default")) ||
      (is.numeric(n_hvgs) && n_hvgs <= 0)) {
    n_hvgs <- min(500L, max(100L, as.integer(round(0.1 * n_genes_total))))
    n_hvgs_source <- sprintf("auto (0.1 x %d genes, clamped 100..500)", n_genes_total)
    cat(sprintf("Auto n_hvgs = %d [%s]\n", n_hvgs, n_hvgs_source))
  } else if (is.numeric(n_hvgs) && n_hvgs > n_genes_total) {
    cat(sprintf(
      "\u26A0 n_hvgs (%d) exceeds total genes (%d); capping at %d\n",
      n_hvgs, n_genes_total, n_genes_total
    ))
    n_hvgs <- n_genes_total
    n_hvgs_source <- "capped_at_n_genes"
  }
  
  # Identify highly variable genes
  if (spatial_hvg) {
    cat("Identifying spatially informed feature rankings...\n")
    cat("  Target ranked features:", n_hvgs, "\n")
    cat("  Note: this mode uses Giotto spatial enrichment ranking as a proxy, not a dedicated SVG model.\n")
    cat("  For true SVG detection on large cohorts, prefer nnSVG or a dedicated binSpect workflow.\n")
    
    # Use Giotto spatial enrichment ranking as a lightweight spatial proxy
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
    
    cat("✓ Selected", length(hvg_genes), "spatially ranked features\n\n")
    
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
  
  # Scree plot - capture the ggplot so we can attach a title that matches the
  # rest of the pipeline's presentation style, then save via the shared helper.
  # Giotto's screePlot returns a ggplot but applies its own theme at the end,
  # which silently overwrites the title set by + ggtitle(). Re-apply both the
  # title AND presentation_theme() AFTER capturing the object.
  tryCatch({
    scree_p <- screePlot(
      gobject     = gobj,
      ncp         = min(30, n_pcs),
      save_plot   = FALSE,
      return_plot = TRUE,
      show_plot   = FALSE
    )
    scree_p <- scree_p +
      labs(
        title    = sample_plot_title(sample_id, "PCA scree plot"),
        subtitle = NULL,
        x        = "Principal component",
        y        = "Variance explained"
      ) +
      presentation_theme(base_size = 12)
    save_presentation_plot(
      plot     = scree_p,
      filename = file.path(results_folder, paste0(sample_id, "_scree_plot.png")),
      width    = 10,
      height   = 6,
      dpi      = 300,
      bg       = "white"
    )
    cat("✓ Scree plot saved\n\n")
  }, error = function(e) {
    cat("⚠ Scree plot warning:", conditionMessage(e), "\n\n")
  })
  
  # UMAP
  cat("Running UMAP...\n")
  cat("  Neighbors:", umap_n_neighbors, "\n")
  cat("  Min distance:", umap_min_dist, "\n")
  cat("  Seed:", seed, "\n")

  set.seed(seed)
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
  cat("  Seed:", seed, "\n")

  set.seed(seed)
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

    # UMAP coloured by FOV - batch-drift / spatial-batch indicator.
    # Hide legend when there are too many FOVs for a readable key.
    if ("fov" %in% names(pDataDT(gobj))) {
      umap_fov <- .prepare_dim_plot_data(gobj, "umap", "fov")
      umap_fov$fov <- factor(umap_fov$fov)
      n_fov <- length(levels(umap_fov$fov))
      fov_plot <- ggplot(umap_fov, aes(x = dim1, y = dim2, color = fov)) +
        geom_point(size = 0.35, alpha = 0.7) +
        scale_color_viridis_d(option = "turbo") +
        labs(
          title = sample_plot_title(sample_id, "UMAP Embedding Colored By FOV"),
          x = embedding_axis_label("UMAP", 1),
          y = embedding_axis_label("UMAP", 2),
          color = "FOV"
        ) +
        presentation_theme(base_size = 12) +
        (if (n_fov > 20) theme(legend.position = "none")
         else guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))))
      save_presentation_plot(
        plot = fov_plot,
        filename = file.path(results_folder,
                             paste0(sample_id, "_umap_by_fov.png")),
        width = 12, height = 8, dpi = 300, bg = "white"
      )
    }

    # Composite-only: UMAP coloured by sub-biopsy. Each cell's FOV is mapped
    # to the sub-biopsy whose [fov_min, fov_max] range it falls inside. Useful
    # for spotting biopsy-level batch structure inside a merged composite.
    sub_rows <- discover_composite_subsamples(sample_row, sample_sheet_path)
    if (!is.null(sub_rows) && nrow(sub_rows) > 0 &&
        "fov" %in% names(pDataDT(gobj))) {
      umap_sub <- .prepare_dim_plot_data(gobj, "umap", "fov")
      assign_sub <- function(fov_val) {
        if (is.na(fov_val)) return(NA_character_)
        hit <- which(fov_val >= as.integer(sub_rows$fov_min) &
                     fov_val <= as.integer(sub_rows$fov_max))
        if (length(hit) == 0) NA_character_ else as.character(sub_rows$sample_id[hit[1]])
      }
      umap_sub$sub_sample_id <- vapply(umap_sub$fov, assign_sub, character(1))
      umap_sub <- umap_sub[!is.na(umap_sub$sub_sample_id), , drop = FALSE]

      if (nrow(umap_sub) > 0) {
        umap_sub$sub_sample_id <- factor(
          umap_sub$sub_sample_id,
          levels = sort(unique(umap_sub$sub_sample_id))
        )
        sub_colours <- cluster_palette(levels(umap_sub$sub_sample_id))
        sub_plot <- ggplot(
          umap_sub,
          aes(x = dim1, y = dim2, color = sub_sample_id)
        ) +
          geom_point(size = 0.35, alpha = 0.7) +
          scale_color_manual(values = sub_colours) +
          labs(
            title = sample_plot_title(sample_id, "UMAP Embedding Colored By Sample"),
            x     = embedding_axis_label("UMAP", 1),
            y     = embedding_axis_label("UMAP", 2),
            color = "Sample"
          ) +
          presentation_theme(base_size = 12) +
          guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
        save_presentation_plot(
          plot     = sub_plot,
          filename = file.path(results_folder,
                               paste0(sample_id, "_umap_by_sample.png")),
          width    = 12, height = 8, dpi = 300, bg = "white"
        )
        cat("✓ Composite UMAP by sub-biopsy saved (",
            nlevels(umap_sub$sub_sample_id), " sub-samples)\n", sep = "")
      }
    }

    cat("✓ Plots saved\n\n")
  }, error = function(e) {
    cat("⚠ Plotting warning:", conditionMessage(e), "\n\n")
  })
  
  # Summary
  cat("=== Dimensionality Reduction Summary ===\n")
  cat("Method:", ifelse(spatial_hvg, "Spatial proxy ranking", "Standard HVG (sparse variance)"), "\n")
  cat("Features selected:", length(hvg_genes), "\n")
  cat("PCs computed:", n_pcs, "\n")
  
  tryCatch({
    umap_dims <- getDimReduction(
      gobject = gobj,
      spat_unit = "cell",
      feat_type = "rna",
      reduction = "cells",
      reduction_method = "umap",
      name = "umap",
      output = "matrix"
    )
    tsne_dims <- getDimReduction(
      gobject = gobj,
      spat_unit = "cell",
      feat_type = "rna",
      reduction = "cells",
      reduction_method = "tsne",
      name = "tsne",
      output = "matrix"
    )
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

    # Resolve optional composite context so the CLI entrypoint can emit the
    # sub-biopsy UMAP without expanding its argument signature.
    cli_sheet_path <- Sys.getenv("COSMX_SAMPLE_SHEET", unset = "")
    cli_sample_row <- NULL
    if (nzchar(cli_sheet_path) && file.exists(cli_sheet_path) &&
        exists("safe_read_sheet", mode = "function")) {
      sheet_df <- tryCatch(
        as.data.frame(safe_read_sheet(cli_sheet_path), stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (!is.null(sheet_df) && "sample_id" %in% names(sheet_df)) {
        hit <- sheet_df[as.character(sheet_df$sample_id) == sample_id, , drop = FALSE]
        if (nrow(hit) >= 1) cli_sample_row <- hit[1, , drop = FALSE]
      }
    }

    gobj <- dimensionality_reduction(
      gobj              = input_path,
      sample_id         = sample_id,
      output_dir        = output_dir,
      spatial_hvg       = spatial_hvg,
      sample_row        = cli_sample_row,
      sample_sheet_path = if (nzchar(cli_sheet_path)) cli_sheet_path else NULL
    )
    
    saveGiotto(gobj, dir = output_dir, foldername = "Giotto_Object_DimReduced", overwrite = TRUE)
  } else {
    stop("Usage: Rscript 04_Dimensionallity_Reduction.R <sample_id> <input_path> <output_dir> [spatial_hvg]")
  }
}

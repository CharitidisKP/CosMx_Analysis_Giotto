# Create comprehensive visualizations and final summary -------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 08_visualization.R
# ==============================================================================

#' Create Comprehensive Visualizations
#'
#' @param gobj Giotto object or path
#' @param sample_id Sample identifier
#' @param output_dir Output directory
#' @param celltype_columns Vector of cell type annotation columns
#' @param cluster_column Clustering column name
#' @param marker_genes Optional vector of marker genes to highlight
#' @return Giotto object (unchanged)

.muffle_known_giotto_plot_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      if (grepl("`aes_string\\(\\)` was deprecated", msg, fixed = FALSE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

.filter_celltype_columns <- function(metadata, columns) {
  if (length(columns) == 0) {
    return(columns)
  }
  
  present_cols <- columns[columns %in% names(metadata)]
  if (length(present_cols) == 0) {
    return(character(0))
  }
  
  excluded_score_cols <- present_cols[grepl("(^score_|_score$|score_)", present_cols)]
  candidate_cols <- setdiff(present_cols, excluded_score_cols)
  
  excluded_numeric_cols <- candidate_cols[
    vapply(candidate_cols, function(col) is.numeric(metadata[[col]]), logical(1))
  ]
  
  filtered_cols <- setdiff(candidate_cols, excluded_numeric_cols)
  
  excluded_cols <- unique(c(excluded_score_cols, excluded_numeric_cols))
  if (length(excluded_cols) > 0) {
    cat(
      "Excluded non-categorical annotation columns:",
      paste(excluded_cols, collapse = ", "),
      "\n"
    )
  }
  
  filtered_cols
}

create_visualizations <- function(gobj,
                                  sample_id,
                                  output_dir,
                                  celltype_columns = NULL,
                                  cluster_column = "leiden_clust",
                                  marker_genes = NULL) {
  
  cat("\n========================================\n")
  cat("STEP 08: Comprehensive visualization\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }
  
  results_folder <- file.path(output_dir, "08_Visualization")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Auto-detect cell type columns if not provided
  metadata <- pDataDT(gobj)
  
  if (is.null(celltype_columns)) {
    metadata_cols <- names(metadata)
    celltype_columns <- grep("^celltype_", metadata_cols, value = TRUE)
    celltype_columns <- .filter_celltype_columns(metadata, celltype_columns)
    
    if (length(celltype_columns) > 0) {
      cat("Auto-detected cell type columns:", paste(celltype_columns, collapse = ", "), "\n\n")
    }
  } else {
    celltype_columns <- .filter_celltype_columns(metadata, celltype_columns)
  }
  
  # Auto-detect cluster column
  if (!cluster_column %in% names(metadata)) {
    metadata_cols <- names(metadata)
    cluster_cols <- grep("leiden|louvain|cluster", metadata_cols, value = TRUE, ignore.case = TRUE)
    
    if (length(cluster_cols) > 0) {
      cluster_column <- cluster_cols[1]
      cat("Auto-detected cluster column:", cluster_column, "\n\n")
    }
  }
  
  # ============================================================================
  # Section 1: Dimensionality Reduction Plots
  # ============================================================================
  
  cat("Creating dimensionality reduction plots...\n")
  dimred_folder <- file.path(results_folder, "01_DimReduction")
  dir.create(dimred_folder, recursive = TRUE, showWarnings = FALSE)
  
  # UMAP - Clusters
  tryCatch({
    .muffle_known_giotto_plot_warnings(
      dimPlot2D(
        gobject = gobj,
        dim_reduction_to_use = "umap",
        cell_color = cluster_column,
        point_size = 0.8,
        show_legend = TRUE,
        save_plot = TRUE,
        save_param = list(
          save_name = paste0(sample_id, "_umap_clusters"),
          save_dir = dimred_folder,
          base_width = 12,
          base_height = 9
        )
      )
    )
    cat("  ✓ UMAP clusters\n")
  }, error = function(e) {
    cat("  ⚠ UMAP clusters warning\n")
  })
  
  # UMAP - Cell types
  for (ct_col in celltype_columns) {
    tryCatch({
      .muffle_known_giotto_plot_warnings(
        dimPlot2D(
          gobject = gobj,
          dim_reduction_to_use = "umap",
          cell_color = ct_col,
          point_size = 0.8,
          show_legend = TRUE,
          save_plot = TRUE,
          save_param = list(
            save_name = paste0(sample_id, "_umap_", ct_col),
            save_dir = dimred_folder,
            base_width = 12,
            base_height = 9
          )
        )
      )
      cat("  ✓ UMAP", ct_col, "\n")
    }, error = function(e) {
      cat("  ⚠", ct_col, "warning\n")
    })
  }
  
  # t-SNE - Clusters
  tryCatch({
    .muffle_known_giotto_plot_warnings(
      dimPlot2D(
        gobject = gobj,
        dim_reduction_to_use = "tsne",
        cell_color = cluster_column,
        point_size = 0.8,
        show_legend = TRUE,
        save_plot = TRUE,
        save_param = list(
          save_name = paste0(sample_id, "_tsne_clusters"),
          save_dir = dimred_folder,
          base_width = 12,
          base_height = 9
        )
      )
    )
    cat("  ✓ t-SNE clusters\n")
  }, error = function(e) {
    cat("  ⚠ t-SNE warning\n")
  })
  
  cat("✓ Dimensionality reduction plots complete\n\n")
  
  # ============================================================================
  # Section 2: Spatial Plots
  # ============================================================================
  
  cat("Creating spatial plots...\n")
  spatial_folder <- file.path(results_folder, "02_Spatial")
  dir.create(spatial_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Spatial - Clusters
  tryCatch({
    .muffle_known_giotto_plot_warnings(
      spatPlot2D(
        gobject = gobj,
        cell_color = cluster_column,
        point_size = 0.5,
        show_image = FALSE,
        point_alpha = 0.8,
        save_plot = TRUE,
        save_param = list(
          save_name = paste0(sample_id, "_spatial_clusters"),
          save_dir = spatial_folder,
          base_width = 14,
          base_height = 10
        )
      )
    )
    cat("  ✓ Spatial clusters\n")
  }, error = function(e) {
    cat("  ⚠ Spatial clusters warning\n")
  })
  
  # Spatial - Cell types
  for (ct_col in celltype_columns) {
    tryCatch({
      .muffle_known_giotto_plot_warnings(
        spatPlot2D(
          gobject = gobj,
          cell_color = ct_col,
          point_size = 0.5,
          show_image = FALSE,
          point_alpha = 0.8,
          save_plot = TRUE,
          save_param = list(
            save_name = paste0(sample_id, "_spatial_", ct_col),
            save_dir = spatial_folder,
            base_width = 14,
            base_height = 10
          )
        )
      )
      cat("  ✓ Spatial", ct_col, "\n")
    }, error = function(e) {
      cat("  ⚠", ct_col, "warning\n")
    })
  }
  
  # Spatial - QC metrics
  tryCatch({
    .muffle_known_giotto_plot_warnings(
      spatPlot2D(
        gobject = gobj,
        cell_color = "nr_feats",
        color_as_factor = FALSE,
        gradient_style = "sequential",
        point_size = 0.5,
        show_image = FALSE,
        save_plot = TRUE,
        save_param = list(
          save_name = paste0(sample_id, "_spatial_genes_per_cell"),
          save_dir = spatial_folder,
          base_width = 14,
          base_height = 10
        )
      )
    )
    cat("  ✓ Spatial genes per cell\n")
  }, error = function(e) {
    cat("  ⚠ Spatial QC warning\n")
  })
  
  cat("✓ Spatial plots complete\n\n")
  
  # ============================================================================
  # Section 3: Feature Plots
  # ============================================================================
  
  if (!is.null(marker_genes) && length(marker_genes) > 0) {
    
    cat("Creating feature plots for", length(marker_genes), "genes...\n")
    feature_folder <- file.path(results_folder, "03_Features")
    dir.create(feature_folder, recursive = TRUE, showWarnings = FALSE)
    
    # Spatial feature plots
    for (gene in marker_genes) {
      tryCatch({
        .muffle_known_giotto_plot_warnings(
          spatFeatPlot2D(
            gobject = gobj,
            expression_values = "normalized",
            feats = gene,
            point_size = 0.5,
            show_image = FALSE,
            save_plot = TRUE,
            save_param = list(
              save_name = paste0(sample_id, "_spatial_", gene),
              save_dir = feature_folder,
              base_width = 12,
              base_height = 10
            )
          )
        )
        cat("  ✓", gene, "\n")
      }, error = function(e) {
        cat("  ⚠", gene, "warning\n")
      })
    }
    
    cat("✓ Feature plots complete\n\n")
  }
  
  # ============================================================================
  # Section 4: Summary Plots
  # ============================================================================
  
  cat("Creating summary plots...\n")
  summary_folder <- file.path(results_folder, "04_Summary")
  dir.create(summary_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Cell type proportions
  for (ct_col in celltype_columns) {
    tryCatch({
      metadata <- metadata %>% as_tibble()
      
      type_counts <- metadata %>%
        dplyr::group_by(!!rlang::sym(ct_col)) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(n))
      
      p <- ggplot(type_counts, aes(x = reorder(!!rlang::sym(ct_col), n), y = n, fill = !!rlang::sym(ct_col))) +
        geom_bar(stat = "identity") +
        coord_flip() +
        labs(title = paste(sample_id, "-", gsub("celltype_", "", ct_col)),
             x = "Cell type",
             y = "Number of cells") +
        theme_classic() +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text = element_text(size = 10)
        )
      
      ggsave(
        filename = file.path(summary_folder, paste0(sample_id, "_proportions_", ct_col, ".png")),
        plot = p,
        width = 10,
        height = max(6, nrow(type_counts) * 0.3),
        dpi = 300
      )
      
      cat("  ✓ Proportions", ct_col, "\n")
    }, error = function(e) {
      cat("  ⚠ Proportions", ct_col, "warning\n")
    })
  }
  
  # Cluster composition by cell type
  if (length(celltype_columns) > 0) {
    for (ct_col in celltype_columns) {
      tryCatch({
        metadata <- metadata %>% as_tibble()
        
        composition <- metadata %>%
          dplyr::group_by(!!rlang::sym(cluster_column), !!rlang::sym(ct_col)) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(!!rlang::sym(cluster_column)) %>%
          dplyr::mutate(prop = n / sum(n))
        
        p <- ggplot(composition, aes(x = !!rlang::sym(cluster_column), y = prop, fill = !!rlang::sym(ct_col))) +
          geom_bar(stat = "identity", position = "fill") +
          labs(title = paste(sample_id, "- Cluster composition"),
               subtitle = gsub("celltype_", "", ct_col),
               x = "Cluster",
               y = "Proportion",
               fill = "Cell type") +
          theme_classic() +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1)
          )
        
        ggsave(
          filename = file.path(summary_folder, paste0(sample_id, "_composition_", ct_col, ".png")),
          plot = p,
          width = 12,
          height = 8,
          dpi = 300
        )
        
        cat("  ✓ Composition", ct_col, "\n")
      }, error = function(e) {
        cat("  ⚠ Composition", ct_col, "warning\n")
      })
    }
  }
  
  cat("✓ Summary plots complete\n\n")
  
  # ============================================================================
  # Section 5: Generate HTML Report
  # ============================================================================
  
  cat("Generating HTML summary report...\n")
  
  tryCatch({
    # Create simple HTML report
    html_content <- paste0(
      "<!DOCTYPE html>\n<html>\n<head>\n",
      "<title>CosMx Analysis Report - ", sample_id, "</title>\n",
      "<style>\n",
      "body { font-family: Arial, sans-serif; margin: 40px; }\n",
      "h1 { color: #2c3e50; }\n",
      "h2 { color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; }\n",
      "img { max-width: 100%; height: auto; margin: 20px 0; border: 1px solid #ddd; }\n",
      ".summary { background: #ecf0f1; padding: 20px; border-radius: 5px; margin: 20px 0; }\n",
      "</style>\n",
      "</head>\n<body>\n",
      "<h1>CosMx Analysis Report: ", sample_id, "</h1>\n",
      "<div class='summary'>\n",
      "<h3>Analysis Summary</h3>\n",
      "<p><strong>Total Cells:</strong> ", length(gobj@cell_ID$cell), "</p>\n",
      "<p><strong>Total Genes:</strong> ", length(gobj@feat_ID$rna), "</p>\n"
    )
    
    if (cluster_column %in% names(pDataDT(gobj))) {
      n_clusters <- length(unique(pDataDT(gobj)[[cluster_column]]))
      html_content <- paste0(html_content, "<p><strong>Clusters:</strong> ", n_clusters, "</p>\n")
    }
    
    html_content <- paste0(html_content, "</div>\n")
    
    # Add sections with images
    sections <- list(
      list(title = "Dimensionality Reduction", folder = "01_DimReduction"),
      list(title = "Spatial Visualization", folder = "02_Spatial"),
      list(title = "Summary Statistics", folder = "04_Summary")
    )
    
    for (section in sections) {
      section_path <- file.path(results_folder, section$folder)
      if (dir.exists(section_path)) {
        images <- list.files(section_path, pattern = "\\.png$", full.names = FALSE)
        
        if (length(images) > 0) {
          html_content <- paste0(html_content, "<h2>", section$title, "</h2>\n")
          
          for (img in images) {
            img_rel_path <- file.path(section$folder, img)
            html_content <- paste0(
              html_content,
              "<h3>", gsub("_", " ", gsub(paste0(sample_id, "_|.png"), "", img)), "</h3>\n",
              "<img src='", img_rel_path, "' alt='", img, "'>\n"
            )
          }
        }
      }
    }
    
    html_content <- paste0(html_content, "</body>\n</html>")
    
    writeLines(html_content, file.path(results_folder, paste0(sample_id, "_analysis_report.html")))
    
    cat("✓ HTML report created:", file.path(results_folder, paste0(sample_id, "_analysis_report.html")), "\n")
  }, error = function(e) {
    cat("⚠ HTML report generation warning\n")
  })
  
  cat("\n✓ All visualizations complete for", sample_id, "\n\n")
  
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
    
    gobj <- create_visualizations(
      gobj = input_path,
      sample_id = sample_id,
      output_dir = output_dir
    )
  } else {
    stop("Usage: Rscript 08_Visualisation.R <sample_id> <input_path> <output_dir>")
  }
}

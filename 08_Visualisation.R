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
#' @param max_cells_preview Optional cap for preview plotting on very large
#'                          datasets. When set, Giotto plots are generated on a
#'                          random cell subset while composition summaries still
#'                          use the full object.
#' @param preview_seed Random seed for preview-cell downsampling
#' @return Giotto object (unchanged)

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
if ((!exists("presentation_theme") || !exists("sample_plot_title") ||
     !exists("pretty_plot_label") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
  source(pipeline_utils)
}

# Gene-expression polygon renderer for Section 3 feature plots.
# spatInSituPlotPoints() does not reliably accept a gene name as polygon_fill
# across Giotto versions, so we draw directly with ggplot2 using the
# polygon table from .extract_polygon_df() (defined in 07_Annotation.R).
# Returns a ggplot or NULL if polygon data / expression can't be assembled.
.plot_gene_polygons <- function(gobj, gene, title_txt, norm_mat = NULL) {
  if (is.null(norm_mat)) {
    norm_mat <- tryCatch(
      getExpression(gobj, values = "normalized", output = "matrix"),
      error = function(e) NULL
    )
  }
  if (is.null(norm_mat) || !(gene %in% rownames(norm_mat))) return(NULL)

  poly_df <- tryCatch(.extract_polygon_df(gobj), error = function(e) NULL)
  if (is.null(poly_df) || nrow(poly_df) == 0) return(NULL)
  poly_df$poly_group <- paste(poly_df$cell_ID, poly_df$geom,
                               poly_df$part, sep = "_")

  # Map expression values onto the polygon table by cell_ID.
  expr_vec <- norm_mat[gene, , drop = TRUE]
  poly_df$expr <- unname(expr_vec[poly_df$cell_ID])

  ggplot2::ggplot(poly_df,
                  ggplot2::aes(x = x, y = y,
                               group = poly_group, fill = expr)) +
    ggplot2::geom_polygon(colour = "grey30", linewidth = 0.08) +
    ggplot2::scale_fill_gradient(low = "lightgrey", high = "red",
                                 name = gene, na.value = "grey85") +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = title_txt, x = NULL, y = NULL) +
    presentation_theme(base_size = 11, legend_position = "right") +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

# Find the sub-biopsy split rows belonging to a composite sample.
# Returns a data.frame (rows with sample_id, fov_min, fov_max, group_id,
# patient_id, timepoint) or NULL if the current sample is not a composite
# or the sheet is unavailable.
.discover_composite_subsamples <- function(sample_row, sample_sheet_path) {
  if (is.null(sample_row) || is.null(sample_sheet_path)) return(NULL)
  if (!file.exists(sample_sheet_path)) return(NULL)

  role <- tryCatch(as.character(sample_row$split_role)[1],
                   error = function(e) NA_character_)
  if (is.null(role) || is.na(role) || !nzchar(role) || role != "composite") {
    return(NULL)
  }

  slide_num <- tryCatch(as.character(sample_row$slide_num)[1],
                        error = function(e) NA_character_)

  sheet <- tryCatch(safe_read_sheet(sample_sheet_path),
                    error = function(e) NULL)
  if (is.null(sheet)) return(NULL)
  sheet <- as.data.frame(sheet, stringsAsFactors = FALSE)

  matches_slide <- if (!is.na(slide_num) && "slide_num" %in% names(sheet)) {
    as.character(sheet$slide_num) == slide_num
  } else rep(TRUE, nrow(sheet))

  keep <- matches_slide &
    !is.na(sheet$split_role) &
    as.character(sheet$split_role) == "split"

  subs <- sheet[keep, , drop = FALSE]
  if (nrow(subs) == 0) return(NULL)
  # Require fov_min / fov_max to be set so we can actually subset.
  subs <- subs[!is.na(subs$fov_min) & !is.na(subs$fov_max), , drop = FALSE]
  if (nrow(subs) == 0) return(NULL)
  subs
}

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

.prepare_preview_gobject <- function(gobj, max_cells_preview = NULL, preview_seed = 1) {
  ## Changed from: ##
  ## total_cells <- length(gobj@cell_ID$cell) ##
  cell_ids_slot <- tryCatch(gobj@cell_ID$cell, error = function(e) NULL)
  all_cell_ids  <- if (!is.null(cell_ids_slot) && length(cell_ids_slot) > 0) {
    cell_ids_slot
  } else {
    GiottoClass::getCellID(gobj)   # Giotto Suite accessor
  }
  total_cells <- length(all_cell_ids)
  
  if (is.null(max_cells_preview)) {
    return(list(gobj = gobj, is_preview = FALSE, total_cells = total_cells, preview_cells = total_cells))
  }
  
  max_cells_preview <- suppressWarnings(as.integer(max_cells_preview)[1])
  if (!is.finite(max_cells_preview) || max_cells_preview <= 0 || total_cells <= max_cells_preview) {
    return(list(gobj = gobj, is_preview = FALSE, total_cells = total_cells, preview_cells = total_cells))
  }
  
  set.seed(preview_seed)
  
  ## Changed from: ##
  ## preview_ids <- sample(gobj@cell_ID$cell, max_cells_preview) ##
  ## Check again ##
  preview_ids_slot <- tryCatch(sample(gobj@cell_ID$cell, max_cells_preview), error = function(e) NULL)
  preview_ids <- if (!is.null(preview_ids_slot) && length(preview_ids_slot) > 0) {
    preview_ids_slot
  } else {
    sample(GiottoClass::getCellID(gobj), max_cells_preview)
  }
  preview_gobj <- subsetGiotto(gobject = gobj, cell_ids = preview_ids)
  
  list(
    gobj = preview_gobj,
    is_preview = TRUE,
    total_cells = total_cells,
    preview_cells = length(preview_ids)
  )
}

create_visualizations <- function(gobj,
                                  sample_id,
                                  output_dir,
                                  celltype_columns = NULL,
                                  cluster_column = "leiden_clust",
                                  marker_genes = NULL,
                                  max_cells_preview = NULL,
                                  preview_seed = 1,
                                  sample_row = NULL,
                                  sample_sheet_path = NULL) {
  
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
  
  preview_info <- .prepare_preview_gobject(
    gobj = gobj,
    max_cells_preview = max_cells_preview,
    preview_seed = preview_seed
  )
  plot_gobj <- preview_info$gobj
  if (isTRUE(preview_info$is_preview)) {
    cat(
      "Using a ", preview_info$preview_cells,
      "-cell preview subset for Giotto plot rendering (full object retained for summaries).\n\n",
      sep = ""
    )
  }
  
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
        gobject = plot_gobj,
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
    cat("  ⚠ UMAP clusters failed:", conditionMessage(e), "\n") 
    }
  )
  
  # UMAP - Cell types
  for (ct_col in celltype_columns) {
    tryCatch({
      .muffle_known_giotto_plot_warnings(
        dimPlot2D(
          gobject = plot_gobj,
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
      cat("  ⚠", ct_col, "warning", conditionMessage(e), "\n") 
      }
    )
  }
  
  # t-SNE - Clusters
  tryCatch({
    .muffle_known_giotto_plot_warnings(
      dimPlot2D(
        gobject = plot_gobj,
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
    cat("  ⚠ t-SNE failed:", conditionMessage(e), "\n") 
    }
  )
  
  cat("✓ Dimensionality reduction plots complete\n\n")
  
  # ============================================================================
  # Section 2+3: Polygon-based Spatial & Feature Plots
  # ============================================================================
  # Both spatial (celltype / cluster) and feature (gene-expression) plots use
  # the shared polygon renderer .spat_in_situ_outlined() defined in
  # 07_Annotation.R. Falls back to the previous spatPlot2D / spatFeatPlot2D
  # path only if the polygon helper returns NULL (e.g. no polygon data).

  # Genes on the panel but absent from the expression matrix (dropped by QC
  # or never loaded) cause Giotto's internal subscript to fail opaquely. Filter
  # the list up front so the loop reports a clean "skipped" line instead.
  norm_mat <- tryCatch(
    getExpression(plot_gobj, values = "normalized", output = "matrix"),
    error = function(e) NULL
  )
  available_genes <- if (!is.null(norm_mat)) rownames(norm_mat) else character(0)

  .render_spatial_and_features <- function(gobj_local,
                                            spatial_dir,
                                            feature_dir,
                                            title_row,
                                            sample_label) {
    dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(feature_dir, recursive = TRUE, showWarnings = FALSE)

    meta_local <- as.data.frame(pDataDT(gobj_local))

    # --- Spatial: clusters --------------------------------------------------
    if (cluster_column %in% names(meta_local)) {
      tryCatch({
        clust_vals <- as.character(meta_local[[cluster_column]])
        cmap_clust <- .build_colour_map_safe(clust_vals)
        p <- .spat_in_situ_outlined(
          gobj           = gobj_local,
          fill_col       = cluster_column,
          fill_as_factor = TRUE,
          colour_map     = cmap_clust,
          legend_title   = "Cluster",
          title_txt      = plot_title_for_sample(title_row, "Clusters",
                                                 sample_id_fallback = sample_label)
        )
        if (!is.null(p)) {
          save_presentation_plot(
            plot     = p,
            filename = file.path(spatial_dir,
                                 paste0(sample_label, "_spatial_clusters.png")),
            width = 14, height = 10, dpi = 300
          )
          cat("  ✓ Spatial clusters (polygon)\n")
        } else {
          cat("  ⚠ Spatial clusters: polygon renderer unavailable\n")
        }
      }, error = function(e) {
        cat("  ⚠ Spatial clusters failed:", conditionMessage(e), "\n")
      })
    }

    # --- Spatial: cell-type columns ----------------------------------------
    for (ct_col in celltype_columns) {
      if (!ct_col %in% names(meta_local)) next
      tryCatch({
        ct_vals <- as.character(meta_local[[ct_col]])
        cmap_ct <- .build_colour_map_safe(ct_vals)
        title_subtitle <- paste(pretty_plot_label(ct_col))
        p <- .spat_in_situ_outlined(
          gobj           = gobj_local,
          fill_col       = ct_col,
          fill_as_factor = TRUE,
          colour_map     = cmap_ct,
          legend_title   = "Cell Type",
          title_txt      = plot_title_for_sample(title_row, title_subtitle,
                                                 sample_id_fallback = sample_label)
        )
        if (!is.null(p)) {
          save_presentation_plot(
            plot     = p,
            filename = file.path(spatial_dir,
                                 paste0(sample_label, "_spatial_", ct_col, ".png")),
            width = 14, height = 10, dpi = 300
          )
          cat("  ✓ Spatial", ct_col, "(polygon)\n")
        } else {
          cat("  ⚠ Spatial", ct_col, ": polygon renderer unavailable\n")
        }
      }, error = function(e) {
        cat("  ⚠", ct_col, "failed:", conditionMessage(e), "\n")
      })
    }

    # --- Feature plots ------------------------------------------------------
    if (!is.null(marker_genes) && length(marker_genes) > 0) {
      norm_local <- tryCatch(
        getExpression(gobj_local, values = "normalized", output = "matrix"),
        error = function(e) NULL
      )
      local_genes   <- if (!is.null(norm_local)) rownames(norm_local) else available_genes
      genes_kept    <- intersect(marker_genes, local_genes)
      genes_skipped <- setdiff(marker_genes, local_genes)
      if (length(genes_skipped) > 0) {
        cat("  ℹ Skipping", length(genes_skipped),
            "gene(s) absent from expression matrix:",
            paste(genes_skipped, collapse = ", "), "\n")
      }
      if (length(genes_kept) > 0) {
        cat("Creating feature plots for", length(genes_kept), "genes...\n")
        for (gene in genes_kept) {
          tryCatch({
            p <- .plot_gene_polygons(
              gobj      = gobj_local,
              gene      = gene,
              norm_mat  = norm_local,
              title_txt = plot_title_for_sample(
                title_row, paste(gene, "Expression"),
                sample_id_fallback = sample_label
              )
            )
            if (!is.null(p)) {
              save_presentation_plot(
                plot     = p,
                filename = file.path(feature_dir,
                                     paste0(sample_label, "_spatial_", gene, ".png")),
                width = 12, height = 10, dpi = 300
              )
              cat("  ✓", gene, "\n")
            } else {
              cat("  ⚠", gene, ": polygon data / expression unavailable\n")
            }
          }, error = function(e) {
            cat("  ⚠", gene, "failed:", conditionMessage(e), "\n")
          })
        }
      }
    }
    invisible(NULL)
  }

  # Defensive colour-map builder: prefer 07's .build_colour_map if present in
  # the sourced env, else synthesise something usable locally.
  if (!exists(".build_colour_map_safe", inherits = FALSE)) {
    .build_colour_map_safe <- function(vals) {
      if (exists(".build_colour_map", mode = "function", inherits = TRUE)) {
        return(.build_colour_map(vals))
      }
      lev <- sort(unique(stats::na.omit(as.character(vals))))
      stats::setNames(
        grDevices::hcl.colors(max(length(lev), 1), palette = "Set3"),
        lev
      )
    }
  }

  # --- Main call: full object -------------------------------------------------
  cat("Creating spatial + feature plots (main)...\n")
  .render_spatial_and_features(
    gobj_local   = plot_gobj,
    spatial_dir  = file.path(results_folder, "02_Spatial"),
    feature_dir  = file.path(results_folder, "03_Features"),
    title_row    = sample_row,
    sample_label = sample_id
  )
  cat("✓ Main spatial + feature plots complete\n\n")

  # --- Per-sub-biopsy plots (CART composite) ---------------------------------
  # If this sample is a composite with sub-biopsy split rows in the sample
  # sheet, emit the same plots filtered by FOV range into per-biopsy subfolders.
  sub_rows <- .discover_composite_subsamples(sample_row, sample_sheet_path)
  if (!is.null(sub_rows) && nrow(sub_rows) > 0) {
    cat("Composite sample detected — rendering", nrow(sub_rows),
        "per-sub-biopsy subfolder(s)...\n")
    meta_all <- as.data.frame(pDataDT(plot_gobj))
    if (!"fov" %in% names(meta_all)) {
      cat("  ⚠ No 'fov' column on Giotto object; skipping sub-biopsy split\n")
    } else {
      for (k in seq_len(nrow(sub_rows))) {
        sub_r  <- sub_rows[k, , drop = FALSE]
        sub_id <- as.character(sub_r$sample_id)
        fmin   <- as.integer(sub_r$fov_min)
        fmax   <- as.integer(sub_r$fov_max)
        if (anyNA(c(fmin, fmax))) {
          cat("  ⚠ ", sub_id, ": fov_min/fov_max missing, skipped\n", sep = "")
          next
        }
        cell_ids <- meta_all$cell_ID[
          !is.na(meta_all$fov) & meta_all$fov >= fmin & meta_all$fov <= fmax
        ]
        if (length(cell_ids) == 0) {
          cat("  ⚠ ", sub_id, ": no cells in FOV ", fmin, "-", fmax,
              ", skipped\n", sep = "")
          next
        }
        sub_gobj <- tryCatch(
          subsetGiotto(plot_gobj, cell_ids = cell_ids),
          error = function(e) {
            cat("  ⚠ ", sub_id, ": subsetGiotto failed: ",
                conditionMessage(e), "\n", sep = "")
            NULL
          }
        )
        if (is.null(sub_gobj)) next
        cat("  → ", sub_id, " (FOV ", fmin, "-", fmax,
            ", ", length(cell_ids), " cells)\n", sep = "")
        .render_spatial_and_features(
          gobj_local   = sub_gobj,
          spatial_dir  = file.path(results_folder, "02_Spatial", sub_id),
          feature_dir  = file.path(results_folder, "03_Features", sub_id),
          title_row    = sub_r,
          sample_label = sub_id
        )
      }
      cat("✓ Sub-biopsy spatial + feature plots complete\n\n")
    }
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
      
      ct_display <- pretty_plot_label(gsub("celltype_", "", ct_col), width = 24)
      type_counts$label <- factor(
        pretty_plot_label(type_counts[[ct_col]], width = 24),
        levels = pretty_plot_label(type_counts[[ct_col]], width = 24)
      )
      
      p <- ggplot(type_counts, aes(x = reorder(label, n), y = n, fill = label)) +
        geom_col() +
        coord_flip() +
        labs(
          title = sample_plot_title(sample_id, paste(ct_display, "Cell Type Composition")),
          subtitle = "Cell counts ranked from most abundant to least abundant annotation.",
          x = "Cell Type",
          y = "Number of Cells"
        ) +
        presentation_theme(base_size = 12) +
        theme(legend.position = "none")
      
      save_presentation_plot(
        plot = p,
        filename = file.path(summary_folder, paste0(sample_id, "_proportions_", ct_col, ".png")),
        width = 10,
        height = max(6, nrow(type_counts) * 0.3),
        dpi = 300
      )
      
      cat("  ✓ Proportions", ct_col, "\n")
    }, error = function(e) {
      cat("  ⚠ Proportions", ct_col, "failed:", conditionMessage(e), "\n") 
      }
    )
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
        
        ct_display <- pretty_plot_label(gsub("celltype_", "", ct_col), width = 24)
        composition$celltype_label <- pretty_plot_label(composition[[ct_col]], width = 24)
        
        p <- ggplot(composition, aes(x = !!rlang::sym(cluster_column), y = prop, fill = celltype_label)) +
          geom_bar(stat = "identity", position = "fill") +
          labs(
            title = sample_plot_title(sample_id, "Cluster Composition by Annotation"),
            subtitle = ct_display,
            x = "Cluster",
            y = "Cell Type Proportion",
            fill = "Cell Type"
          ) +
          presentation_theme(base_size = 12, x_angle = 45)
        
        save_presentation_plot(
          plot = p,
          filename = file.path(summary_folder, paste0(sample_id, "_composition_", ct_col, ".png")),
          width = 12,
          height = 8,
          dpi = 300
        )
        
        cat("  ✓ Composition", ct_col, "\n")
      }, error = function(e) {
        cat("  ⚠ Composition", ct_col, "failed:", conditionMessage(e), "\n") 
        }
      )
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
  
  # Multi-panel per-sample overview (patchwork: UMAP + spatial + composition)
  tryCatch({
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      cat("\u26A0 patchwork not installed; skipping multi-panel overview\n")
    } else {
      meta_ov <- as.data.frame(pDataDT(gobj))
      meta_cols <- names(meta_ov)
      # Prefer selected celltype, then best supervised, then leiden_clust
      colour_col <- NULL
      if ("celltype" %in% meta_cols &&
          sum(!is.na(meta_ov$celltype) & nzchar(meta_ov$celltype)) > 0) {
        colour_col <- "celltype"
      } else {
        sup_cols <- grep("^celltype_.*_supervised$", meta_cols, value = TRUE)
        if (length(sup_cols)) {
          colour_col <- sup_cols[1]
        } else if ("leiden_clust" %in% meta_cols) {
          colour_col <- "leiden_clust"
        }
      }
      have_umap <- tryCatch({
        !is.null(getDimReduction(gobject = gobj, spat_unit = "cell",
                                 feat_type = "rna", reduction = "cells",
                                 reduction_method = "umap", name = "umap",
                                 output = "matrix"))
      }, error = function(e) FALSE)
      have_spat <- all(c("sdimx", "sdimy") %in%
                       names(getSpatialLocations(gobj, output = "data.table")))

      if (!is.null(colour_col) && have_umap && have_spat) {
        umap_mat <- getDimReduction(gobject = gobj, spat_unit = "cell",
                                    feat_type = "rna", reduction = "cells",
                                    reduction_method = "umap", name = "umap",
                                    output = "matrix")
        spat_dt <- as.data.frame(getSpatialLocations(gobj, output = "data.table"))
        ids <- rownames(umap_mat)
        if (is.null(ids) || !length(ids)) ids <- meta_ov$cell_ID[seq_len(nrow(umap_mat))]
        ov <- data.frame(
          cell_ID = ids,
          umap1 = umap_mat[, 1], umap2 = umap_mat[, 2],
          stringsAsFactors = FALSE
        )
        ov <- merge(ov, spat_dt[, c("cell_ID", "sdimx", "sdimy")], by = "cell_ID")
        ov <- merge(ov, meta_ov[, c("cell_ID", colour_col)], by = "cell_ID")
        ov[[colour_col]] <- factor(ov[[colour_col]])

        p_umap <- ggplot2::ggplot(ov,
            ggplot2::aes(x = umap1, y = umap2,
                         colour = .data[[colour_col]])) +
          ggplot2::geom_point(size = 0.35, alpha = 0.8) +
          ggplot2::labs(title = "UMAP", x = "UMAP 1", y = "UMAP 2") +
          presentation_theme(base_size = 10) +
          ggplot2::guides(colour = ggplot2::guide_legend(
            ncol = 1, override.aes = list(size = 2.5), title = NULL
          ))

        p_spat <- ggplot2::ggplot(ov,
            ggplot2::aes(x = sdimx, y = sdimy,
                         colour = .data[[colour_col]])) +
          ggplot2::geom_point(size = 0.25, alpha = 0.85) +
          ggplot2::coord_fixed() +
          ggplot2::labs(title = "Spatial", x = NULL, y = NULL) +
          presentation_theme(base_size = 10) +
          ggplot2::theme(legend.position = "none",
                         axis.text = ggplot2::element_blank(),
                         axis.ticks = ggplot2::element_blank(),
                         panel.grid = ggplot2::element_blank())

        comp_tbl <- as.data.frame(table(ov[[colour_col]]))
        names(comp_tbl) <- c("label", "n")
        comp_tbl <- comp_tbl[order(-comp_tbl$n), ]
        comp_tbl$label <- factor(comp_tbl$label, levels = comp_tbl$label)
        p_comp <- ggplot2::ggplot(comp_tbl,
            ggplot2::aes(x = label, y = n, fill = label)) +
          ggplot2::geom_col() +
          ggplot2::labs(title = "Composition", x = NULL, y = "Cells") +
          presentation_theme(base_size = 10) +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "none"
          )

        overview <- patchwork::wrap_plots(
          p_umap, p_spat, p_comp,
          design = "AAB\nAAB\nCCC", heights = c(1, 1, 0.8)
        ) +
          patchwork::plot_annotation(
            title = sample_plot_title(sample_id, "Sample overview"),
            subtitle = paste0("Coloured by: ", colour_col)
          )
        save_presentation_plot(
          plot     = overview,
          filename = file.path(results_folder,
                               paste0(sample_id, "_overview_multipanel.png")),
          width    = 16, height = 12, dpi = 300
        )
        cat("  \u2713 Multi-panel overview saved\n")
      }
    }
  }, error = function(e) {
    cat("\u26A0 Multi-panel overview failed:", conditionMessage(e), "\n")
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

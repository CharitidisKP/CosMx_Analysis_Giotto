# Functions for creating clustering visualization dataframes and plots --------
# Designed to work with any clustering method stored in Giotto metadata

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

pipeline_utils <- file.path(current_script_dir(), "Pipeline_Utils.R")
if ((!exists("presentation_theme") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
  source(pipeline_utils)
}

## Make sure cell type annotations dont ##
## convert characters to NA instead of numeric ##
.make_cluster_factor <- function(vals) {
  char_vals <- as.character(vals)
  int_vals  <- suppressWarnings(as.integer(char_vals))
  if (all(!is.na(int_vals))) {
    factor(char_vals, levels = as.character(sort(unique(int_vals))))
  } else {
    factor(char_vals, levels = sort(unique(char_vals)))
  }
}

# Function 1: Extract and prepare visualization dataframes ----------------
prepare_clustering_dataframes <- function(gobject, 
                                          cluster_column,
                                          spat_unit = "cell",
                                          umap_name = "umap",
                                          tsne_name = "tsne") {
  
  cat("Extracting coordinates and metadata...\n")
  
  # Get UMAP coordinates
  umap_coords <- getDimReduction(
    gobject = gobject,
    reduction = "cells",
    reduction_method = "umap",
    name = umap_name,
    output = "matrix"
  )
  
  # Get tSNE coordinates
  tsne_coords <- getDimReduction(
    gobject = gobject,
    reduction = "cells",
    reduction_method = "tsne",
    name = tsne_name,
    output = "matrix"
  )
  
  # Get spatial locations
  spat_locs <- getSpatialLocations(
    gobject = gobject,
    spat_unit = spat_unit,
    output = "data.table"
  )
  
  # Get metadata
  metadata <- pDataDT(gobject) %>%
    as_tibble()
  
  # Check if cluster column exists
  if (!cluster_column %in% colnames(metadata)) {
    stop(paste("Cluster column", cluster_column, "not found in metadata"))
  }
  
  # Create UMAP dataframe with ordered factors
  umap_df <- umap_coords %>%
    as.data.frame() %>%
    setNames(c("UMAP_1", "UMAP_2")) %>%
    rownames_to_column(var = "cell_ID") %>%
    inner_join(
      metadata %>% select(cell_ID, all_of(cluster_column)),
      by = "cell_ID"
    ) %>%
    mutate(cluster = .make_cluster_factor(.data[[cluster_column]]))
  ## Fixed from: ##
    # mutate(
    #   cluster = factor(
    #     as.character(.data[[cluster_column]]),
    #     levels = as.character(sort(unique(as.integer(.data[[cluster_column]]))))
    #   )
    # )
  
  # Create tSNE dataframe with ordered factors
  tsne_df <- tsne_coords %>%
    as.data.frame() %>%
    setNames(c("tSNE_1", "tSNE_2")) %>%
    rownames_to_column(var = "cell_ID") %>%
    inner_join(
      metadata %>% select(cell_ID, all_of(cluster_column)),
      by = "cell_ID"
    ) %>%
    mutate(cluster = .make_cluster_factor(.data[[cluster_column]]))
  ## Fixed from: ##
  # mutate(
  #   cluster = factor(
  #     as.character(.data[[cluster_column]]),
  #     levels = as.character(sort(unique(as.integer(.data[[cluster_column]]))))
  #   )
  # )
  
  # Create spatial dataframe with ordered factors
  spatial_df <- spat_locs %>%
    as_tibble() %>%
    select(cell_ID, sdimx, sdimy) %>%
    inner_join(
      metadata %>% select(cell_ID, all_of(cluster_column)),
      by = "cell_ID"
    ) %>%
    mutate(cluster = .make_cluster_factor(.data[[cluster_column]]))
  ## Fixed from: ##
  # mutate(
  #   cluster = factor(
  #     as.character(.data[[cluster_column]]),
  #     levels = as.character(sort(unique(as.integer(.data[[cluster_column]]))))
  #   )
  # )
  
  # Get cluster information
  cluster_levels <- levels(umap_df$cluster)
  n_clusters <- length(cluster_levels)
  
  cat("Number of clusters:", n_clusters, "\n")
  cat("Cluster IDs:", paste(cluster_levels, collapse = ", "), "\n")
  
  return(list(
    umap_df = umap_df,
    tsne_df = tsne_df,
    spatial_df = spatial_df,
    cluster_levels = cluster_levels,
    n_clusters = n_clusters
  ))
}

# Function 2: Generate cluster colors -------------------------------------

generate_cluster_colors <- function(cluster_levels, 
                                    palette = "Paired",
                                    custom_colors = NULL) {
  
  n_clusters <- length(cluster_levels)
  
  if (!is.null(custom_colors) && length(custom_colors) >= n_clusters) {
    clus_colors <- custom_colors[1:n_clusters]
  } else {
    ## Check if this is properly running or if I need to change the syntax ##
    if (requireNamespace("RColorBrewer", quietly = TRUE) &&
        palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      clus_colors <- colorRampPalette(RColorBrewer::brewer.pal(max_n, palette))(n_clusters)
    } else {
      clus_colors <- scales::hue_pal()(n_clusters)   # scales is already used in the repo
    }
  }
  
  names(clus_colors) <- as.character(cluster_levels)
  
  cat("\nCluster color mapping:\n")
  print(data.frame(
    Cluster = names(clus_colors),
    Color = clus_colors,
    row.names = NULL
  ))
  
  return(clus_colors)
}

# Function 3: Create individual clustering plots --------------------------

create_clustering_plots <- function(umap_df, 
                                    tsne_df, 
                                    spatial_df,
                                    clus_colors,
                                    point_size_umap = 1,
                                    point_size_tsne = 1,
                                    point_size_spatial = 0.5,
                                    alpha = 0.8,
                                    title_suffix = "") {
  
  # Subtitles on the individual panels just say the projection type - the
  # main panel title (passed via title_suffix on the combined plot) already
  # carries the "Leiden" framing, so we don't repeat it here.
  cat("Creating UMAP plot...\n")
  custom_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size = point_size_umap, alpha = alpha, stroke = 0) +
    scale_color_manual(
      values = clus_colors,
      name = "Leiden cluster"
    ) +
    labs(
      title    = paste0("UMAP projection", title_suffix),
      subtitle = "UMAP",
      x        = "UMAP dimension 1",
      y        = "UMAP dimension 2"
    ) +
    presentation_theme(base_size = 12, legend_position = "right") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    add_centroid_labels(umap_df, "cluster", x_col = "UMAP_1", y_col = "UMAP_2")

  cat("Creating tSNE plot...\n")
  custom_tsne <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
    geom_point(size = point_size_tsne, alpha = alpha, stroke = 0) +
    scale_color_manual(
      values = clus_colors,
      name = "Leiden cluster"
    ) +
    labs(
      title    = paste0("t-SNE projection", title_suffix),
      subtitle = "t-SNE",
      x        = "t-SNE dimension 1",
      y        = "t-SNE dimension 2"
    ) +
    presentation_theme(base_size = 12, legend_position = "right") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    add_centroid_labels(tsne_df, "cluster", x_col = "tSNE_1", y_col = "tSNE_2")

  # The custom spatial ggplot has been retired from this helper: per-pipeline
  # rule, spatial cell plots must render actual cell polygons (not points).
  # 05_Clustering.R invokes plot_cells_polygon() directly after the dim-reduction
  # panels are produced. The combined figure therefore shows UMAP | t-SNE only.
  # On the combined plot we strip the per-panel title (the patchwork-level
  # title already says "Leiden clusters") and keep just the projection-type
  # subtitle so each subplot is clearly labelled.
  cat("Creating combined plot (UMAP | t-SNE)...\n")
  combined_umap <- custom_umap + labs(title = NULL)
  combined_tsne <- custom_tsne + labs(title = NULL)
  combined_plot <- (combined_umap + combined_tsne) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("Clustering results", title_suffix),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
    )

  return(list(
    custom_umap   = custom_umap,
    custom_tsne   = custom_tsne,
    combined_plot = combined_plot
  ))
}

# Function 4: Complete workflow wrapper -----------------------------------

create_clustering_visualization <- function(gobject,
                                            cluster_column,
                                            spat_unit = "cell",
                                            umap_name = "umap",
                                            tsne_name = "tsne",
                                            palette = "Paired",
                                            custom_colors = NULL,
                                            point_size_umap = 1,
                                            point_size_tsne = 1,
                                            point_size_spatial = 0.5,
                                            alpha = 0.8,
                                            title_suffix = "",
                                            save_plots = FALSE,
                                            save_dir = NULL,
                                            prefix = "clustering") {
  
  cat("\n=== Creating clustering visualization for:", cluster_column, "===\n")
  
  # Step 1: Prepare dataframes
  data_list <- prepare_clustering_dataframes(
    gobject = gobject,
    cluster_column = cluster_column,
    spat_unit = spat_unit,
    umap_name = umap_name,
    tsne_name = tsne_name
  )
  
  # Step 2: Generate colors
  clus_colors <- generate_cluster_colors(
    cluster_levels = data_list$cluster_levels,
    palette = palette,
    custom_colors = custom_colors
  )
  
  # Step 3: Create plots
  plots_list <- create_clustering_plots(
    umap_df = data_list$umap_df,
    tsne_df = data_list$tsne_df,
    spatial_df = data_list$spatial_df,
    clus_colors = clus_colors,
    point_size_umap = point_size_umap,
    point_size_tsne = point_size_tsne,
    point_size_spatial = point_size_spatial,
    alpha = alpha,
    title_suffix = title_suffix
  )
  
  # Step 4: Save plots if requested
  if (save_plots && !is.null(save_dir)) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    
    save_presentation_plot(
      plot = plots_list$custom_umap,
      filename = file.path(save_dir, paste0(prefix, "_umap.png")),
      width = 11,
      height = 8,
      dpi = 300,
      bg = "white"
    )

    save_presentation_plot(
      plot = plots_list$custom_tsne,
      filename = file.path(save_dir, paste0(prefix, "_tsne.png")),
      width = 11,
      height = 8,
      dpi = 300,
      bg = "white"
    )

    # `_spatial.png` removed. The spatial panel is now a polygon plot emitted
    # by 05_Clustering.R via plot_cells_polygon().

    save_presentation_plot(
      plot = plots_list$combined_plot,
      filename = file.path(save_dir, paste0(prefix, "_combined.png")),
      width = 18,
      height = 8,
      dpi = 300,
      bg = "white"
    )

    cat("\nPlots saved to:", save_dir, "\n")
  }
  
  # Return everything
  return(list(
    dataframes = data_list,
    colors = clus_colors,
    plots = plots_list
  ))
}

# Example usage:
# leiden_viz <- create_clustering_visualization(
#   gobject = cosmx_dim,
#   cluster_column = "leiden_clus",
#   title_suffix = " - Leiden",
#   save_plots = TRUE,
#   save_dir = results_folder,
#   prefix = "leiden_clustering"
# )
# 
# other_viz <- create_clustering_visualization(
#   gobject = cosmx_dim,
#   cluster_column = "other_clus",
#   title_suffix = " - Other Method",
#   save_plots = TRUE,
#   save_dir = results_folder,
#   prefix = "other_clustering"
# )
# 
# final_comparison <- (leiden_viz$plots$combined_plot) / 
#                     (other_viz$plots$combined_plot) +
#   plot_annotation(
#     title = "Leiden vs Other Clustering Method",
#     theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
#   )

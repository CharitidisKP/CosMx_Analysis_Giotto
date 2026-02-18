# Functions for creating clustering visualization dataframes and plots --------
# Designed to work with any clustering method stored in Giotto metadata


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
    mutate(
      cluster = factor(
        as.character(.data[[cluster_column]]),
        levels = as.character(sort(unique(as.integer(.data[[cluster_column]]))))
      )
    )
  
  # Create tSNE dataframe with ordered factors
  tsne_df <- tsne_coords %>%
    as.data.frame() %>%
    setNames(c("tSNE_1", "tSNE_2")) %>%
    rownames_to_column(var = "cell_ID") %>%
    inner_join(
      metadata %>% select(cell_ID, all_of(cluster_column)),
      by = "cell_ID"
    ) %>%
    mutate(
      cluster = factor(
        as.character(.data[[cluster_column]]),
        levels = as.character(sort(unique(as.integer(.data[[cluster_column]]))))
      )
    )
  
  # Create spatial dataframe with ordered factors
  spatial_df <- spat_locs %>%
    as_tibble() %>%
    select(cell_ID, sdimx, sdimy) %>%
    inner_join(
      metadata %>% select(cell_ID, all_of(cluster_column)),
      by = "cell_ID"
    ) %>%
    mutate(
      cluster = factor(
        as.character(.data[[cluster_column]]),
        levels = as.character(sort(unique(as.integer(.data[[cluster_column]]))))
      )
    )
  
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
    clus_colors <- getColors(palette, n_clusters)
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
  
  cat("Creating UMAP plot...\n")
  custom_umap <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
    geom_point(size = point_size_umap, alpha = alpha, stroke = 0) +
    scale_color_manual(
      values = clus_colors,
      name = "Cluster"
    ) +
    labs(
      title = paste0("UMAP Projection", title_suffix),
      x = "UMAP Dimension 1",
      y = "UMAP Dimension 2"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  cat("Creating tSNE plot...\n")
  custom_tsne <- ggplot(tsne_df, aes(x = tSNE_1, y = tSNE_2, color = cluster)) +
    geom_point(size = point_size_tsne, alpha = alpha, stroke = 0) +
    scale_color_manual(
      values = clus_colors,
      name = "Cluster"
    ) +
    labs(
      title = paste0("t-SNE Projection", title_suffix),
      x = "t-SNE Dimension 1",
      y = "t-SNE Dimension 2"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  cat("Creating spatial plot...\n")
  custom_spatial <- ggplot(spatial_df, aes(x = sdimx, y = sdimy, color = cluster)) +
    geom_point(size = point_size_spatial, alpha = alpha, stroke = 0) +
    scale_color_manual(
      values = clus_colors,
      name = "Cluster"
    ) +
    labs(
      title = paste0("Spatial Distribution", title_suffix),
      x = "Spatial X (microns)",
      y = "Spatial Y (microns)"
    ) +
    coord_fixed() +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  cat("Creating combined plot...\n")
  combined_plot <- (custom_umap + custom_tsne + custom_spatial) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("Clustering Results", title_suffix),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
    )
  
  return(list(
    custom_umap = custom_umap,
    custom_tsne = custom_tsne,
    custom_spatial = custom_spatial,
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
    
    ggsave(
      filename = file.path(save_dir, paste0(prefix, "_umap.png")),
      plot = plots_list$custom_umap,
      width = 11,
      height = 8,
      dpi = 300,
      bg = "white"
    )
    
    ggsave(
      filename = file.path(save_dir, paste0(prefix, "_tsne.png")),
      plot = plots_list$custom_tsne,
      width = 11,
      height = 8,
      dpi = 300,
      bg = "white"
    )
    
    ggsave(
      filename = file.path(save_dir, paste0(prefix, "_spatial.png")),
      plot = plots_list$custom_spatial,
      width = 14,
      height = 12,
      dpi = 300,
      bg = "white"
    )
    
    ggsave(
      filename = file.path(save_dir, paste0(prefix, "_combined.png")),
      plot = plots_list$combined_plot,
      width = 24,
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
# SOLO PROCESSING FUNCTION FOR SINGLE SAMPLES ------------------------

#' Inspect and Visualize Nearest Neighbor Network
#'
#' @param gobject Giotto object containing the nearest neighbor network
#' @param nn_name Name of the nearest neighbor network (default: "sNN.pca")
#' @param output_dir Directory to save output files (default: "nn_analysis")
#' @param sample_id Sample identifier (e.g., "19H827D")
#' @param subfolder_name Subfolder name within sample directory (default: "sNN_results")
#' @param n_cells_subgraph Number of cells to sample for network visualization (default: 1000)
#' @param save_plots Logical, whether to save plots (default: TRUE)
#' @param verbose Logical, whether to print messages (default: TRUE)
#'
#' @return A list containing summary statistics and plots
#' @export

inspect_nn_network <- function(gobject,
                               nn_name = "sNN.pca",
                               output_dir = "nn_analysis",
                               sample_id = "sample",
                               subfolder_name = "sNN_results",
                               n_cells_subgraph = 1000,
                               save_plots = TRUE,
                               verbose = TRUE) {
  
  # Load required libraries
  required_packages <- c("ggplot2", "patchwork", "pheatmap", "igraph", 
                         "ggraph", "tidygraph", "gridExtra", "grid", "knitr", "kableExtra")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (verbose) message(paste("Installing package:", pkg))
      install.packages(pkg, quiet = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
  
  # Create full output path: output_dir/sample_id/subfolder_name
  full_output_path <- file.path(output_dir, sample_id, subfolder_name)
  
  # Expand tilde if present
  full_output_path <- path.expand(full_output_path)
  
  # Create output directory
  if (save_plots && !dir.exists(full_output_path)) {
    dir.create(full_output_path, recursive = TRUE)
    if (verbose) message(paste("Created output directory:", full_output_path))
  }
  
  if (verbose) message(paste("\n========== Analyzing Network for Sample:", sample_id, "==========\n"))
  
  # 1. EXTRACT NETWORK DATA -------------------------------------------------

  if (verbose) message("Step 1: Extracting nearest neighbor network...")
  
  nn <- getNearestNetwork(
    gobject = gobject,
    nn_type = "sNN",
    name = nn_name
  )
  
  nn_igraph <- nn@igraph
  
  # Store original edge attributes before conversion
  original_edge_attrs <- list(
    weight = E(nn_igraph)$weight,
    distance = E(nn_igraph)$distance,
    shared = E(nn_igraph)$shared,
    rank = E(nn_igraph)$rank
  )
  
  # Convert to undirected graph if needed (FIX FOR LOUVAIN ERROR)
  if (is_directed(nn_igraph)) {
    if (verbose) message("  - Converting directed graph to undirected...")
    
    # Use mode = "each" to keep all edges and preserve attributes better
    nn_igraph <- as_undirected(nn_igraph, mode = "each", edge.attr.comb = "first")
    
    # Check if attributes were preserved, if not, restore them
    if (is.null(E(nn_igraph)$weight)) {
      if (verbose) message("  - Restoring edge attributes after conversion...")
      E(nn_igraph)$weight <- original_edge_attrs$weight
      E(nn_igraph)$distance <- original_edge_attrs$distance
      E(nn_igraph)$shared <- original_edge_attrs$shared
      E(nn_igraph)$rank <- original_edge_attrs$rank
    }
  }
  
  # Basic network properties
  n_vertices <- vcount(nn_igraph)
  n_edges <- ecount(nn_igraph)
  
  if (verbose) {
    message(paste("  - Number of cells/vertices:", n_vertices))
    message(paste("  - Number of edges:", n_edges))
    message(paste("  - Graph is directed:", is_directed(nn_igraph)))
  }
  

  # 2. EXTRACT EDGE ATTRIBUTES ----------------------------------------------

  if (verbose) message("\nStep 2: Extracting edge attributes...")
  
  edge_weights <- E(nn_igraph)$weight
  edge_distances <- E(nn_igraph)$distance
  edge_shared <- E(nn_igraph)$shared
  edge_ranks <- E(nn_igraph)$rank
  
  # Check if attributes exist and have correct length
  if (is.null(edge_weights) || length(edge_weights) == 0) {
    stop("Edge weights are missing or empty. Check the network object.")
  }
  
  if (verbose) {
    message(paste("  - Edge weights: ", length(edge_weights), "values"))
    message(paste("  - Edge distances: ", length(edge_distances), "values"))
    message(paste("  - Edge shared: ", length(edge_shared), "values"))
    message(paste("  - Edge ranks: ", length(edge_ranks), "values"))
  }
  
  # Handle missing attributes by creating defaults
  if (is.null(edge_distances) || length(edge_distances) == 0) {
    if (verbose) message("  - Warning: Edge distances missing, using NA")
    edge_distances <- rep(NA, length(edge_weights))
  }
  
  if (is.null(edge_shared) || length(edge_shared) == 0) {
    if (verbose) message("  - Warning: Edge shared neighbors missing, using NA")
    edge_shared <- rep(NA, length(edge_weights))
  }
  
  if (is.null(edge_ranks) || length(edge_ranks) == 0) {
    if (verbose) message("  - Warning: Edge ranks missing, using NA")
    edge_ranks <- rep(NA, length(edge_weights))
  }
  
  # Create edge data frame
  edge_data <- data.frame(
    weight = edge_weights,
    distance = edge_distances,
    shared = edge_shared,
    rank = edge_ranks
  )
  
  # Remove rows with all NAs (if any attribute was completely missing)
  edge_data_complete <- edge_data[complete.cases(edge_data), ]
  
  if (nrow(edge_data_complete) < nrow(edge_data) && verbose) {
    message(paste("  - Warning: Removed", nrow(edge_data) - nrow(edge_data_complete), 
                  "edges with missing attributes"))
  }
  
  
  # 3. CALCULATE SUMMARY STATISTICS -----------------------------------------

  
  if (verbose) message("\nStep 3: Calculating summary statistics...")
  
  # Edge attribute summaries (use complete data for statistics)
  edge_summaries <- data.frame(
    Attribute = c("Weight", "Distance", "Shared Neighbors", "Rank"),
    Min = c(min(edge_data$weight, na.rm = TRUE),
            min(edge_data$distance, na.rm = TRUE),
            min(edge_data$shared, na.rm = TRUE),
            min(edge_data$rank, na.rm = TRUE)),
    Q1 = c(quantile(edge_data$weight, 0.25, na.rm = TRUE),
           quantile(edge_data$distance, 0.25, na.rm = TRUE),
           quantile(edge_data$shared, 0.25, na.rm = TRUE),
           quantile(edge_data$rank, 0.25, na.rm = TRUE)),
    Median = c(median(edge_data$weight, na.rm = TRUE),
               median(edge_data$distance, na.rm = TRUE),
               median(edge_data$shared, na.rm = TRUE),
               median(edge_data$rank, na.rm = TRUE)),
    Mean = c(mean(edge_data$weight, na.rm = TRUE),
             mean(edge_data$distance, na.rm = TRUE),
             mean(edge_data$shared, na.rm = TRUE),
             mean(edge_data$rank, na.rm = TRUE)),
    Q3 = c(quantile(edge_data$weight, 0.75, na.rm = TRUE),
           quantile(edge_data$distance, 0.75, na.rm = TRUE),
           quantile(edge_data$shared, 0.75, na.rm = TRUE),
           quantile(edge_data$rank, 0.75, na.rm = TRUE)),
    Max = c(max(edge_data$weight, na.rm = TRUE),
            max(edge_data$distance, na.rm = TRUE),
            max(edge_data$shared, na.rm = TRUE),
            max(edge_data$rank, na.rm = TRUE))
  )
  
  # Round numeric columns
  edge_summaries[, 2:7] <- round(edge_summaries[, 2:7], 3)
  
  # Degree distribution
  degree_dist <- degree(nn_igraph)
  
  # Network connectivity
  is_conn <- is_connected(nn_igraph)
  comp <- components(nn_igraph)
  
  # Network metrics
  network_metrics <- data.frame(
    Metric = c("Sample ID",
               "Number of Cells", 
               "Number of Edges", 
               "Average Degree",
               "Min Degree",
               "Max Degree",
               "Network Density",
               "Is Connected?",
               "Number of Components",
               "Largest Component Size",
               "Clustering Coefficient"),
    Value = as.character(c(
      sample_id,
      n_vertices,
      n_edges,
      round(mean(degree_dist), 2),
      min(degree_dist),
      max(degree_dist),
      round(edge_density(nn_igraph), 6),
      is_conn,
      comp$no,
      max(comp$csize),
      round(transitivity(nn_igraph, type = "global"), 3)
    ))
  )
  

  # 4. CREATE VISUALIZATIONS ------------------------------------------------
  
  if (verbose) message("\nStep 4: Creating visualizations...")
  
  # Use complete edge data for plotting
  plot_data <- edge_data_complete
  
  ## 4.1 Edge Attribute Distributions
  p1 <- ggplot(plot_data, aes(x = weight)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    labs(title = "Edge Weights", 
         x = "Weight", 
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p2 <- ggplot(plot_data, aes(x = distance)) +
    geom_histogram(bins = 30, fill = "coral", color = "black", alpha = 0.7) +
    labs(title = "Edge Distances", 
         x = "Distance", 
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p3 <- ggplot(plot_data, aes(x = shared)) +
    geom_histogram(bins = 15, fill = "seagreen", color = "black", alpha = 0.7) +
    labs(title = "Shared Neighbors", 
         x = "Number Shared", 
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p4 <- ggplot(plot_data, aes(x = rank)) +
    geom_histogram(bins = 15, fill = "mediumpurple", color = "black", alpha = 0.7) +
    scale_x_continuous(breaks = seq(min(plot_data$rank, na.rm = TRUE), 
                                    max(plot_data$rank, na.rm = TRUE), by = 1)) +
    labs(title = "Neighbor Ranks", 
         x = "Rank", 
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Combine edge distribution plots
  combined_plot <- (p1 | p2) / (p3 | p4)
  combined_plot <- combined_plot + 
    plot_annotation(title = paste("sNN Network Edge Attributes -", sample_id),
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  
  ## 4.2 Degree Distribution
  degree_df <- data.frame(degree = degree_dist)
  
  p_degree <- ggplot(degree_df, aes(x = degree)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
    geom_vline(aes(xintercept = mean(degree)), 
               color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = mean(degree_dist) + 2, y = Inf, 
             label = paste("Mean =", round(mean(degree_dist), 2)), 
             vjust = 2, color = "red", fontface = "bold") +
    labs(title = paste("Degree Distribution -", sample_id),
         x = "Number of Neighbors (Degree)",
         y = "Frequency") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  
  ## 4.3 Network Graph Visualization (Subgraph) - ROBUST VERSION
  if (verbose) message(paste("  - Creating network subgraph with", n_cells_subgraph, "cells..."))
  
  set.seed(123)
  n_sample <- min(n_cells_subgraph, n_vertices)
  sample_vertices <- sample(vcount(nn_igraph), n_sample)
  nn_subgraph <- induced_subgraph(nn_igraph, sample_vertices)
  
  # Detect communities
  communities <- cluster_louvain(nn_subgraph)
  if (verbose) message(paste("  - Detected", length(unique(communities$membership)), "communities"))
  
  # Get node degree
  node_degrees <- degree(nn_subgraph)
  
  # METHOD 1: ggraph with manual layout
  if (verbose) message("  - Computing graph layout (Fruchterman-Reingold)...")
  layout_coords <- layout_with_fr(nn_subgraph, niter = 500)
  
  # Create data frame for plotting
  node_data <- data.frame(
    x = layout_coords[, 1],
    y = layout_coords[, 2],
    community = as.factor(communities$membership),
    degree = node_degrees,
    vertex_id = V(nn_subgraph)$name
  )
  
  # Get edge list with coordinates
  edge_list <- as_edgelist(nn_subgraph, names = FALSE)
  edge_data <- data.frame(
    x = layout_coords[edge_list[, 1], 1],
    y = layout_coords[edge_list[, 1], 2],
    xend = layout_coords[edge_list[, 2], 1],
    yend = layout_coords[edge_list[, 2], 2]
  )
  
  # Get edge weights if available
  if (!is.null(E(nn_subgraph)$weight)) {
    edge_data$weight <- E(nn_subgraph)$weight
  } else {
    edge_data$weight <- 1
  }
  
  # Debug: Check layout coordinates
  if (verbose) {
    message(paste("  - Layout coordinate range X:", 
                  round(min(layout_coords[,1]), 2), "to", 
                  round(max(layout_coords[,1]), 2)))
    message(paste("  - Layout coordinate range Y:", 
                  round(min(layout_coords[,2]), 2), "to", 
                  round(max(layout_coords[,2]), 2)))
    message(paste("  - Any infinite values:", any(!is.finite(layout_coords))))
  }
  
  # Create ggplot version (more reliable)
  p_network <- ggplot() +
    geom_segment(data = edge_data,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 alpha = 0.1,
                 color = "gray70",
                 linewidth = 0.3) +
    geom_point(data = node_data,
               aes(x = x, y = y, color = community, size = degree),
               alpha = 0.7) +
    scale_size_continuous(range = c(0.5, 3), name = "Degree") +
    scale_color_discrete(name = "Community") +
    labs(title = paste("sNN Network with Community Detection -", sample_id),
         subtitle = paste("Showing", n_sample, 
                          "random cells | Node size = degree | Colour = community | Number of communities: ", 
                          length(unique(communities$membership)))) +
    theme_void() +
    theme(legend.position = "none",
          # legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          plot.margin = margin(10, 10, 10, 10),
          aspect.ratio = 1)
  
  # METHOD 2: Alternative using base igraph plot (save separately)
  if (save_plots) {
    if (verbose) message("  - Creating alternative igraph visualization...")
    
    # Set up colors for communities
    n_communities <- length(unique(communities$membership))
    community_colors <- rainbow(n_communities, alpha = 0.7)[communities$membership]
    
    # Create igraph plot
    png(file.path(full_output_path, paste0(sample_id, "_network_igraph.png")),
        width = 3000, height = 3000, res = 300)
    
    plot(nn_subgraph,
         layout = layout_coords,
         vertex.color = community_colors,
         vertex.size = sqrt(node_degrees) * 2,
         vertex.label = NA,
         vertex.frame.color = NA,
         edge.arrow.size = 0,
         edge.width = 0.3,
         edge.color = adjustcolor("gray50", alpha = 0.2),
         main = paste("sNN Network -", sample_id))
    
    # # Add legend
    # legend("topright",
    #        legend = paste("Community", 1:n_communities),
    #        fill = rainbow(n_communities, alpha = 0.7),
    #        border = NA,
    #        bty = "n",
    #        cex = 0.8)
    
    dev.off()
    
    if (verbose) message("  - Alternative igraph plot saved")
  }
  
  if (verbose) message("  - Network graph created successfully")
  
  ## 4.4 Adjacency Matrix Heatmap (first 100 cells)
  if (verbose) message("  - Creating adjacency matrix heatmap...")
  
  adj_matrix <- as_adjacency_matrix(nn_igraph, 
                                    attr = "weight", 
                                    sparse = TRUE)
  
  n_heatmap <- min(100, n_vertices)
  adj_subset <- as.matrix(adj_matrix[1:n_heatmap, 1:n_heatmap])
  

  # 5. PRINT SUMMARY TABLES -------------------------------------------------
  
  
  if (verbose) {
    message("\n========== NETWORK SUMMARY STATISTICS ==========\n")
    
    cat("\n--- Network Metrics ---\n")
    print(kable(network_metrics, format = "simple", align = 'lr'))
    
    cat("\n--- Edge Attribute Summaries ---\n")
    print(kable(edge_summaries, format = "simple", align = 'lrrrrrr'))
    
    cat("\n--- Rank Frequencies ---\n")
    
    # Safe rank table creation with error handling
    tryCatch({
      rank_freq <- table(edge_data$rank, useNA = "ifany")
      
      if (length(rank_freq) > 0) {
        rank_table <- data.frame(
          Rank = names(rank_freq),
          Frequency = as.vector(rank_freq),
          stringsAsFactors = FALSE
        )
        print(kable(rank_table, format = "simple", align = 'lr', row.names = FALSE))
      } else {
        message("  - No rank data available")
        rank_table <- data.frame(Rank = character(0), Frequency = numeric(0))
      }
    }, error = function(e) {
      message(paste("  - Error creating rank table:", e$message))
      rank_table <- data.frame(Rank = character(0), Frequency = numeric(0))
    })
  }
  

  # 6. SAVE OUTPUTS -----------------------------------------------------

  
  if (save_plots) {
    if (verbose) message("\nStep 5: Saving outputs...")
    
    # Save combined edge distributions
    ggsave(file.path(full_output_path, paste0(sample_id, "_edge_distributions.png")), 
           plot = combined_plot, 
           width = 12, 
           height = 10, 
           dpi = 300)
    
    # Save degree distribution
    ggsave(file.path(full_output_path, paste0(sample_id, "_degree_distribution.png")), 
           plot = p_degree, 
           width = 8, 
           height = 6, 
           dpi = 300)
    
    # Save network visualization
    ggsave(file.path(full_output_path, paste0(sample_id, "_network_graph.png")), 
           plot = p_network, 
           width = 12, 
           height = 10, 
           dpi = 300)
    
    # Save adjacency heatmap
    png(file.path(full_output_path, paste0(sample_id, "_adjacency_heatmap.png")),
        width = 2400, height = 2400, res = 300)
    pheatmap(adj_subset, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE,
             show_rownames = FALSE,
             show_colnames = FALSE,
             main = paste("sNN Connectivity (first", n_heatmap, "cells) -", sample_id))
    dev.off()
    
    # Save summary tables
    write.csv(network_metrics, 
              file.path(full_output_path, paste0(sample_id, "_network_metrics.csv")),
              row.names = FALSE)
    
    write.csv(edge_summaries, 
              file.path(full_output_path, paste0(sample_id, "_edge_summaries.csv")),
              row.names = FALSE)
    
    # Safe rank table save
    tryCatch({
      rank_freq <- table(edge_data$rank, useNA = "ifany")
      
      if (length(rank_freq) > 0) {
        rank_table_save <- data.frame(
          Rank = names(rank_freq),
          Frequency = as.vector(rank_freq),
          stringsAsFactors = FALSE
        )
        write.csv(rank_table_save, 
                  file.path(full_output_path, paste0(sample_id, "_rank_frequencies.csv")),
                  row.names = FALSE)
      } else {
        if (verbose) message("  - Skipping rank frequencies (no data available)")
      }
    }, error = function(e) {
      if (verbose) message(paste("  - Warning: Could not save rank frequencies:", e$message))
    })
    
    if (verbose) message(paste("  - All outputs saved to:", full_output_path))
  }
  


  # 7. RETURN RESULTS -----------------------------------------------------

  
  # Create rank table for return object (with error handling)
  rank_table <- tryCatch({
    rank_freq <- table(edge_data$rank, useNA = "ifany")
    
    if (length(rank_freq) > 0) {
      data.frame(
        Rank = names(rank_freq),
        Frequency = as.vector(rank_freq),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(Rank = character(0), Frequency = numeric(0))
    }
  }, error = function(e) {
    data.frame(Rank = character(0), Frequency = numeric(0))
  })
  
  results <- list(
    sample_id = sample_id,
    output_path = full_output_path,
    network_metrics = network_metrics,
    edge_summaries = edge_summaries,
    rank_frequencies = rank_table,
    degree_distribution = degree_dist,
    plots = list(
      edge_distributions = combined_plot,
      degree_distribution = p_degree,
      network_graph = p_network
    ),
    igraph_object = nn_igraph,
    nn_object = nn
  )
  
  if (verbose) message("\n========== Analysis Complete ==========\n")
  
  return(invisible(results))
  
}



# BATCH PROCESSING FUNCTION FOR MULTIPLE SAMPLES ------------------------

#' Batch Inspect Multiple Nearest Neighbor Networks
#'
#' @param gobject_list Named list of Giotto objects (names should be sample IDs)
#' @param nn_name Name of the nearest neighbor network (default: "sNN.pca")
#' @param output_dir Base output directory (default: current working directory)
#' @param subfolder_name Subfolder name within each sample directory (default: "sNN_results")
#' @param n_cells_subgraph Number of cells to sample for network visualization (default: 1000)
#' @param save_plots Logical, whether to save plots (default: TRUE)
#' @param verbose Logical, whether to print messages (default: TRUE)
#'
#' @return A named list of results for each sample
#' @export

batch_inspect_nn_network <- function(gobject_list,
                                     nn_name = "sNN.pca",
                                     output_dir = getwd(),
                                     subfolder_name = "sNN_results",
                                     n_cells_subgraph = 1000,
                                     save_plots = TRUE,
                                     verbose = TRUE) {
  
  if (!is.list(gobject_list)) {
    stop("gobject_list must be a named list of Giotto objects")
  }
  
  if (is.null(names(gobject_list))) {
    names(gobject_list) <- paste0("sample_", seq_along(gobject_list))
    warning("gobject_list was not named. Assigning default names: ", 
            paste(names(gobject_list), collapse = ", "))
  }
  
  # Expand tilde in output_dir
  output_dir <- path.expand(output_dir)
  
  if (verbose) {
    message("========================================================")
    message("BATCH ANALYSIS: Processing ", length(gobject_list), " samples")
    message("Base output directory: ", output_dir)
    message("========================================================\n")
  }
  
  # Initialize results list
  all_results <- list()
  
  # Process each sample
  for (sample_id in names(gobject_list)) {
    
    tryCatch({
      results <- inspect_nn_network(
        gobject = gobject_list[[sample_id]],
        nn_name = nn_name,
        output_dir = output_dir,
        sample_id = sample_id,
        subfolder_name = subfolder_name,
        n_cells_subgraph = n_cells_subgraph,
        save_plots = save_plots,
        verbose = verbose
      )
      
      all_results[[sample_id]] <- results
      
    }, error = function(e) {
      message(paste("\nERROR processing sample:", sample_id))
      message(paste("Error message:", e$message, "\n"))
      all_results[[sample_id]] <- NULL
    })
  }
  
  # Create comparison summary across all samples
  if (verbose) message("\nCreating cross-sample comparison summary...")
  
  comparison_metrics <- do.call(rbind, lapply(names(all_results), function(s) {
    if (!is.null(all_results[[s]])) {
      metrics <- all_results[[s]]$network_metrics
      data.frame(Sample = s, 
                 Metric = metrics$Metric, 
                 Value = metrics$Value,
                 stringsAsFactors = FALSE)
    }
  }))
  
  # Save comparison table
  if (save_plots && nrow(comparison_metrics) > 0) {
    comparison_dir <- file.path(output_dir, "cross_sample_comparison")
    if (!dir.exists(comparison_dir)) {
      dir.create(comparison_dir, recursive = TRUE)
    }
    
    write.csv(comparison_metrics, 
              file.path(comparison_dir, "all_samples_comparison.csv"),
              row.names = FALSE)
    
    # Create wide format for easier comparison
    comparison_wide <- reshape(comparison_metrics, 
                               idvar = "Metric", 
                               timevar = "Sample", 
                               direction = "wide")
    
    write.csv(comparison_wide, 
              file.path(comparison_dir, "all_samples_comparison_wide.csv"),
              row.names = FALSE)
    
    if (verbose) message(paste("  - Comparison tables saved to:", comparison_dir))
  }
  
  if (verbose) {
    message("\n========================================================")
    message("BATCH ANALYSIS COMPLETE")
    message("Successfully processed: ", length(all_results), "/", length(gobject_list), " samples")
    message("Results structure: ", output_dir, "/<sample_id>/", subfolder_name)
    message("========================================================\n")
  }
  
  return(invisible(all_results))
}
#!/usr/bin/env Rscript
# ==============================================================================
# 09_merged_batch_correction.R
# Merge Giotto objects and correct slide-level batch effects with Harmony
# ==============================================================================

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
if (!exists("save_giotto_checkpoint") && file.exists(pipeline_utils)) {
  source(pipeline_utils)
}

parse_dimension_vector <- function(x) {
  if (is.null(x)) {
    return(1:30)
  }
  if (is.numeric(x)) {
    return(as.integer(x))
  }
  eval(parse(text = x))
}

embedding_to_tibble <- function(gobj, reduction_method, name = NULL, color_columns = character()) {
  coords <- get_dimReduction(
    gobject = gobj,
    reduction = "cells",
    reduction_method = reduction_method,
    name = name,
    output = "data.table"
  )
  
  coord_df <- as_tibble(coords)
  if (!"cell_ID" %in% names(coord_df)) {
    names(coord_df)[1] <- "cell_ID"
  }
  dim_cols <- setdiff(names(coord_df), "cell_ID")
  dim_cols <- dim_cols[seq_len(min(2, length(dim_cols)))]
  coord_df <- coord_df %>%
    select(cell_ID, all_of(dim_cols))
  names(coord_df)[2:ncol(coord_df)] <- c("dim_1", "dim_2")[seq_len(ncol(coord_df) - 1)]
  
  metadata <- pDataDT(gobj) %>%
    as_tibble() %>%
    select(any_of(c("cell_ID", color_columns)))
  
  coord_df %>%
    left_join(metadata, by = "cell_ID")
}

save_embedding_plot <- function(df, color_column, title, output_file) {
  p <- ggplot(df, aes(x = dim_1, y = dim_2, color = .data[[color_column]])) +
    geom_point(size = 0.4, alpha = 0.8) +
    labs(
      title = title,
      x = "Dimension 1",
      y = "Dimension 2",
      color = color_column
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
}

merge_giotto_samples <- function(gobject_list,
                                 sample_table,
                                 output_dir,
                                 join_method = "shift",
                                 x_padding = 1000,
                                 save_object = FALSE) {
  if (length(gobject_list) == 0) {
    stop("No Giotto objects were provided for merging.")
  }
  
  if (length(gobject_list) == 1) {
    merged_gobj <- gobject_list[[1]]
  } else {
    merged_gobj <- joinGiottoObjects(
      gobject_list = gobject_list,
      gobject_names = sample_table$sample_id,
      join_method = join_method,
      x_padding = x_padding
    )
  }
  
  results_dir <- ensure_dir(file.path(output_dir, "09_Merged"))
  readr::write_csv(sample_table, file.path(results_dir, "merged_sample_manifest.csv"))
  readr::write_csv(as_tibble(pDataDT(merged_gobj)), file.path(results_dir, "merged_cell_metadata.csv"))
  
  if (save_object) {
    save_giotto_checkpoint(
      gobj = merged_gobj,
      checkpoint_dir = file.path(output_dir, "Giotto_Object_Merged"),
      metadata = list(stage = "merge")
    )
  }
  
  merged_gobj
}

batch_correct_merged_object <- function(gobj,
                                        sample_id = "merged",
                                        output_dir,
                                        batch_column = "slide_num",
                                        dimensions_to_use = 1:30,
                                        umap_n_neighbors = 30,
                                        umap_min_dist = 0.3,
                                        k_nn = 15,
                                        resolution = 0.3,
                                        create_plots = TRUE,
                                        save_object = FALSE) {
  cat("\n========================================\n")
  cat("STEP 09: Merged Batch Correction\n")
  cat("Run:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    manifest_path <- file.path(gobj, "manifest.json")
    if (file.exists(manifest_path)) {
      gobj <- load_giotto_checkpoint(gobj)
    } else {
      gobj <- loadGiotto(gobj)
    }
  }
  
  results_dir <- ensure_dir(file.path(output_dir, "09_Merged", "Batch_Correction"))
  dimensions_to_use <- parse_dimension_vector(dimensions_to_use)
  
  metadata <- pDataDT(gobj) %>% as_tibble()
  if (!batch_column %in% names(metadata)) {
    fallback_column <- c("slide_num", "list_ID", "sample_id")[c("slide_num", "list_ID", "sample_id") %in% names(metadata)][1]
    if (is.na(fallback_column) || !nzchar(fallback_column)) {
      stop("Batch column '", batch_column, "' is not present in merged metadata.")
    }
    batch_column <- fallback_column
  }
  
  cat("Running Harmony using batch column:", batch_column, "\n")
  cat("Dimensions:", paste(range(dimensions_to_use), collapse = ":"), "\n\n")
  
  gobj <- runGiottoHarmony(
    gobject = gobj,
    vars_use = batch_column,
    dim_reduction_to_use = "pca",
    dimensions_to_use = dimensions_to_use,
    name = "harmony"
  )
  
  gobj <- runUMAP(
    gobject = gobj,
    dim_reduction_to_use = "harmony",
    dim_reduction_name = "harmony",
    dimensions_to_use = dimensions_to_use,
    n_neighbors = umap_n_neighbors,
    min_dist = umap_min_dist,
    name = "umap"
  )
  
  gobj <- runtSNE(
    gobject = gobj,
    dim_reduction_to_use = "harmony",
    dim_reduction_name = "harmony",
    dimensions_to_use = dimensions_to_use,
    perplexity = 30,
    name = "tsne"
  )
  
  gobj <- createNearestNetwork(
    gobject = gobj,
    dim_reduction_to_use = "harmony",
    dim_reduction_name = "harmony",
    dimensions_to_use = dimensions_to_use,
    k = k_nn,
    name = "NN.harmony"
  )
  
  gobj <- doLeidenCluster(
    gobject = gobj,
    network_name = "NN.harmony",
    resolution = resolution,
    n_iterations = 1000,
    name = "leiden_clust"
  )
  
  readr::write_csv(as_tibble(pDataDT(gobj)), file.path(results_dir, paste0(sample_id, "_batch_corrected_metadata.csv")))
  
  if (create_plots) {
    batch_plot_cols <- c(batch_column, "leiden_clust")
    
    pca_df <- embedding_to_tibble(gobj, reduction_method = "pca", name = NULL, color_columns = batch_plot_cols)
    harmony_df <- embedding_to_tibble(gobj, reduction_method = "harmony", name = "harmony", color_columns = batch_plot_cols)
    umap_df <- embedding_to_tibble(gobj, reduction_method = "umap", name = "umap", color_columns = batch_plot_cols)
    
    save_embedding_plot(
      pca_df,
      color_column = batch_column,
      title = paste(sample_id, "- PCA by", batch_column),
      output_file = file.path(results_dir, paste0(sample_id, "_pca_by_", batch_column, ".png"))
    )
    save_embedding_plot(
      harmony_df,
      color_column = batch_column,
      title = paste(sample_id, "- Harmony by", batch_column),
      output_file = file.path(results_dir, paste0(sample_id, "_harmony_by_", batch_column, ".png"))
    )
    save_embedding_plot(
      umap_df,
      color_column = batch_column,
      title = paste(sample_id, "- Harmony UMAP by", batch_column),
      output_file = file.path(results_dir, paste0(sample_id, "_umap_by_", batch_column, ".png"))
    )
    save_embedding_plot(
      umap_df,
      color_column = "leiden_clust",
      title = paste(sample_id, "- Harmony UMAP by leiden_clust"),
      output_file = file.path(results_dir, paste0(sample_id, "_umap_by_leiden_clust.png"))
    )
  }
  
  if (save_object) {
    save_giotto_checkpoint(
      gobj = gobj,
      checkpoint_dir = file.path(output_dir, "Giotto_Object_BatchCorrected"),
      metadata = list(stage = "batch_correction", batch_column = batch_column)
    )
  }
  
  gobj
}

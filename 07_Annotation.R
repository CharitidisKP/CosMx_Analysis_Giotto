# Cell type annotation using InSituType with multiple reference profiles -------

#!/usr/bin/env Rscript
# ==============================================================================
# 07_Annotation.R
# Cell type annotation using InSituType with multiple reference profiles
# ==============================================================================

#' Load reference profile based on configuration
#' @param profile_config Profile configuration list
#' @return Reference profile matrix
load_reference_profile <- function(profile_config) {
  
  if (profile_config$type == "url") {
    cat("  Downloading from URL...\n")
    ref_data <- read_profile_csv(url = profile_config$source)
    
  } else if (profile_config$type == "hca") {
    cat("  Downloading HCA profile via SpatialDecon...\n")
    ref_data <- SpatialDecon::download_profile_matrix(
      species = profile_config$species,
      age_group = profile_config$age_group,
      matrixname = profile_config$matrixname
    )
    ref_data <- as.matrix(ref_data)
    
  } else if (profile_config$type == "file") {
    cat("  Loading from local file...\n")
    ref_data <- read_profile_csv(url = profile_config$source)
    
  } else {
    stop("Unknown profile type: ", profile_config$type)
  }
  
  return(ref_data)
}

#' Annotate Cells Using InSituType
#'
#' @param gobj Giotto object or path
#' @param sample_id Sample identifier
#' @param output_dir Output directory
#' @param profiles List of profile configurations
#' @param default_profile Name of default profile to use (if NULL, uses first)
#' @param align_genes Align genes between reference and data
#' @param n_clusts_semi Number of novel clusters for semi-supervised
#' @param n_starts Number of random starts (semi-supervised only)
#' @param cohort_column Column to use as cohort for semi-supervised
#' @param min_gene_overlap Minimum gene overlap required (default: 100)
#' @param create_plots Create visualization plots
#' @return Giotto object with annotations
annotate_cells <- function(gobj,
                           sample_id,
                           output_dir,
                           profiles = NULL,
                           default_profile = NULL,
                           align_genes = TRUE,
                           n_clusts_semi = 5,
                           n_starts = 10,
                           cohort_column = "leiden_clus",
                           min_gene_overlap = 100,
                           create_plots = TRUE) {
  
  cat("\n========================================\n")
  cat("STEP 07: Cell Type Annotation\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }
  
  results_folder <- file.path(output_dir, "07_Annotation")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Prepare data for InSituType
  cat("Preparing data for InSituType...\n")
  
  # Get raw counts (genes x cells)
  counts_raw <- getExpression(
    gobject = gobj,
    feat_type = "rna",
    spat_unit = "cell",
    values = "raw",
    output = "matrix"
  )
  
  # Transpose to cells x genes (InSituType format)
  counts_mat <- t(as.matrix(counts_raw))
  
  cat("  Count matrix:", nrow(counts_mat), "cells x", ncol(counts_mat), "genes\n")
  
  # === CORRECT BACKGROUND CALCULATION ===
  # Extract negative control genes from ORIGINAL matrix (before transpose)
  all_feats <- rownames(counts_raw)
  neg_genes <- grep("^Negative", all_feats, value = TRUE, ignore.case = TRUE)
  
  if (length(neg_genes) > 0) {
    cat("  Using", length(neg_genes), "negative control genes for background\n")
    
    # Get negative probe counts (genes x cells)
    neg_counts_mat <- counts_raw[neg_genes, , drop = FALSE]
    
    # Calculate MEAN across negative probes per cell
    bg_per_cell <- colMeans(as.matrix(neg_counts_mat))
    
    cat("  Background per cell - Mean:", round(mean(bg_per_cell), 3), 
        "Median:", round(median(bg_per_cell), 3), "\n")
  } else {
    cat("  No negative probes found, using 1% of mean expression\n")
    mean_counts_per_cell <- Matrix::rowMeans(counts_mat)
    bg_per_cell <- mean_counts_per_cell * 0.01
  }
  
  # Ensure bg_per_cell matches counts_mat order
  names(bg_per_cell) <- rownames(counts_mat)
  
  cat("✓ Data prepared\n\n")
  
  # Store original cell order
  cell_order <- rownames(counts_mat)
  
  # Determine which profile to use as default
  if (!is.null(default_profile)) {
    default_idx <- which(sapply(profiles, function(p) p$name == default_profile))
    if (length(default_idx) == 0) {
      cat("⚠ Default profile '", default_profile, "' not found, using first profile\n")
      default_idx <- 1
    }
  } else {
    default_idx <- 1
  }
  
  cat("Default profile:", profiles[[default_idx]]$name, "\n\n")
  
  # Process each reference profile
  for (i in seq_along(profiles)) {
    
    profile_config <- profiles[[i]]
    profile_name <- profile_config$name
    is_default <- (i == default_idx)
    
    cat("================================================================================\n")
    cat("Running annotation with:", profile_name)
    if (is_default) cat(" [DEFAULT]")
    cat("\n================================================================================\n\n")
    
    profile_folder <- file.path(results_folder, profile_name)
    dir.create(profile_folder, recursive = TRUE, showWarnings = FALSE)
    
    tryCatch({
      
      # Load reference profile
      cat("Loading reference profile...\n")
      ref_profiles <- load_reference_profile(profile_config)
      
      cat("✓ Profile loaded\n")
      cat("  Cell types:", ncol(ref_profiles), "\n")
      cat("  Genes in reference:", nrow(ref_profiles), "\n")
      
      # Check gene overlap
      ref_genes <- rownames(ref_profiles)
      data_genes <- colnames(counts_mat)
      common_genes <- intersect(ref_genes, data_genes)
      
      cat("  Overlapping genes:", length(common_genes), "\n")
      
      if (length(common_genes) < min_gene_overlap) {
        cat("✗ Insufficient overlap (", length(common_genes), " < ", min_gene_overlap, ")\n")
        cat("  Skipping", profile_name, "\n\n")
        next
      }
      
      cat("  Sufficient overlap - proceeding\n\n")
      
      # === SUPERVISED ANNOTATION ===
      cat("Running InSituType (supervised)...\n")
      cat("  Align genes:", align_genes, "\n\n")
      
      insitu_supervised <- InSituType::insitutypeML(
        x = counts_mat,
        neg = bg_per_cell,
        reference_profiles = ref_profiles,
        align_genes = align_genes
      )
      
      cat("✓ Supervised complete\n")
      cat("  Cell types found:", length(unique(insitu_supervised$clust)), "\n\n")
      
      # Create annotation columns
      celltype_col_sup <- paste0("celltype_", profile_name, "_supervised")
      score_col_sup <- paste0("score_", profile_name, "_supervised")
      
      # Also create default annotation columns if this is the default profile
      if (is_default) {
        default_celltype_col <- "celltype"
        default_score_col <- "celltype_score"
      }
      
      # Build annotation data frame
      annot_sup <- data.frame(
        cell_ID = cell_order,
        celltype = as.character(insitu_supervised$clust),
        max_score = apply(insitu_supervised$prob, 1, max),
        stringsAsFactors = FALSE
      )
      
      names(annot_sup)[2:3] <- c(celltype_col_sup, score_col_sup)
      
      # Add default columns if this is default profile
      if (is_default) {
        annot_sup[[default_celltype_col]] <- annot_sup[[celltype_col_sup]]
        annot_sup[[default_score_col]] <- annot_sup[[score_col_sup]]
      }
      
      cat("Adding annotations to Giotto object...\n")
      gobj <- addCellMetadata(gobj, annot_sup, by_column = TRUE, column_cell_ID = "cell_ID")
      
      # Save results
      write_csv(annot_sup, 
                file.path(profile_folder, paste0(sample_id, "_supervised_celltypes.csv")))
      
      cat("✓ Annotations added and saved\n\n")
      
      # === SEMI-SUPERVISED ANNOTATION ===
      if (cohort_column %in% names(pDataDT(gobj)) && n_clusts_semi > 0) {
        
        cat("Running InSituType (semi-supervised)...\n")
        cat("  Novel clusters:", n_clusts_semi, "\n")
        cat("  Cohort column:", cohort_column, "\n")
        cat("  Random starts:", n_starts, "\n\n")
        
        metadata_dt <- pDataDT(gobj)
        cohort_vec <- as.character(metadata_dt[[cohort_column]])
        
        insitu_semi <- InSituType::insitutype(
          x = counts_mat,
          neg = bg_per_cell,
          reference_profiles = ref_profiles,
          n_clusts = n_clusts_semi,
          cohort = cohort_vec,
          align_genes = align_genes,
          n_starts = n_starts
        )
        
        cat("✓ Semi-supervised complete\n\n")
        
        celltype_col_semi <- paste0("celltype_", profile_name, "_semi")
        
        annot_semi <- data.frame(
          cell_ID = cell_order,
          celltype = as.character(insitu_semi$clust),
          stringsAsFactors = FALSE
        )
        names(annot_semi)[2] <- celltype_col_semi
        
        gobj <- addCellMetadata(gobj, annot_semi, by_column = TRUE, column_cell_ID = "cell_ID")
        
        write_csv(annot_semi,
                  file.path(profile_folder, paste0(sample_id, "_semi_celltypes.csv")))
        
        cat("✓ Semi-supervised annotations added\n\n")
      }
      
      # === SUMMARIES ===
      sup_summary <- annot_sup %>%
        group_by(!!sym(celltype_col_sup)) %>%
        summarise(
          n_cells = n(),
          mean_score = mean(!!sym(score_col_sup), na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(desc(n_cells))
      
      write_csv(sup_summary,
                file.path(profile_folder, paste0(sample_id, "_supervised_summary.csv")))
      
      cat("=== Supervised Summary ===\n")
      print(sup_summary)
      cat("\n")
      
      # === VISUALIZATIONS ===
      if (create_plots) {
        cat("Creating visualizations...\n")
        
        # Flight plot
        tryCatch({
          pdf(file.path(profile_folder, paste0(sample_id, "_flightplot.pdf")), 
              width = 12, height = 10)
          InSituType::flightplot_info(
            mat = counts_mat,
            neg = bg_per_cell,
            bg = insitu_supervised$bg,
            clust = insitu_supervised$clust,
            profiles = ref_profiles
          )
          dev.off()
          cat("  ✓ Flight plot\n")
        }, error = function(e) cat("  ⚠ Flight plot warning\n"))
        
        # UMAP
        tryCatch({
          suppressWarnings({
            dimPlot2D(gobj, dim_reduction_to_use = "umap", cell_color = celltype_col_sup,
                      point_size = 0.8, save_plot = TRUE,
                      save_param = list(save_name = paste0(sample_id, "_umap_", profile_name),
                                        save_dir = profile_folder, base_width = 12, base_height = 9))
          })
          cat("  ✓ UMAP\n")
        }, error = function(e) cat("  ⚠ UMAP warning\n"))
        
        # Spatial
        tryCatch({
          suppressWarnings({
            spatPlot2D(gobj, cell_color = celltype_col_sup, point_size = 0.5, 
                       show_image = FALSE, save_plot = TRUE,
                       save_param = list(save_name = paste0(sample_id, "_spatial_", profile_name),
                                         save_dir = profile_folder, base_width = 14, base_height = 10))
          })
          cat("  ✓ Spatial\n")
        }, error = function(e) cat("  ⚠ Spatial warning\n"))
        
        # Proportions
        tryCatch({
          p <- ggplot(sup_summary, aes(x = reorder(!!sym(celltype_col_sup), n_cells), 
                                       y = n_cells, fill = !!sym(celltype_col_sup))) +
            geom_bar(stat = "identity") + coord_flip() +
            labs(title = paste(sample_id, "-", profile_name), x = "Cell Type", y = "Cells") +
            theme_classic() + 
            theme(legend.position = "none", 
                  plot.title = element_text(hjust = 0.5, face = "bold"))
          
          ggsave(file.path(profile_folder, paste0(sample_id, "_proportions.png")),
                 p, width = 10, height = max(6, nrow(sup_summary) * 0.4), dpi = 300)
          cat("  ✓ Proportions\n")
        }, error = function(e) cat("  ⚠ Proportions warning\n"))
        
        cat("✓ Visualizations complete\n\n")
      }
      
      cat("✓ Annotation complete:", profile_name, "\n\n")
      
    }, error = function(e) {
      cat("✗ Error with", profile_name, ":\n")
      cat("  ", conditionMessage(e), "\n\n")
    })
  }
  
  cat("✓ All annotations complete for", sample_id, "\n\n")
  
  return(gobj)
}

# Command-line interface
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) {
    sample_id <- args[1]
    input_path <- args[2]
    output_dir <- args[3]
    config_file <- if (length(args) >= 4) args[4] else NULL
    
    if (!is.null(config_file) && file.exists(config_file)) {
      config <- yaml::read_yaml(config_file)
      
      # Source helper functions
      helper_path <- file.path(config$paths$scripts_dir, "Helper_Scripts/Helper_Functions.R")
      if (file.exists(helper_path)) source(helper_path)
      
      profiles <- config$parameters$annotation$profiles
      default_profile <- config$parameters$annotation$default_profile
      align_genes <- config$parameters$annotation$align_genes
      n_starts <- config$parameters$annotation$n_starts
      n_clusts_semi <- config$parameters$annotation$n_clusts_semi
      cohort_column <- config$parameters$annotation$cohort_column
      min_gene_overlap <- config$parameters$annotation$min_gene_overlap %||% 100
    } else {
      stop("Config file required")
    }
    
    gobj <- annotate_cells(
      gobj = input_path, sample_id = sample_id, output_dir = output_dir,
      profiles = profiles, default_profile = default_profile,
      align_genes = align_genes, n_starts = n_starts,
      n_clusts_semi = n_clusts_semi, cohort_column = cohort_column,
      min_gene_overlap = min_gene_overlap
    )
    
    saveGiotto(gobj, dir = output_dir, foldername = "Giotto_Object_Annotated", overwrite = TRUE)
  }
}
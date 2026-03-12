# Filter cells and genes based on QC metrics ------------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 02_quality_control.R
# Quality control and filtering of CosMx data
# ==============================================================================

#' Quality Control for CosMx Sample
#'
#' @param gobj Giotto object (or path to saved object)
#' @param sample_id Sample identifier
#' @param output_dir Output directory for this sample
#' @param gene_min_cells Minimum cells expressing a gene
#' @param cell_min_genes Minimum genes per cell
#' @param cell_max_genes Maximum genes per cell (NULL = no max)
#' @param min_count Minimum total counts per cell
#' @param max_mito_pct Maximum mitochondrial percentage (NULL = no filter)
#' @return Filtered Giotto object

quality_control <- function(gobj, 
                            sample_id, 
                            output_dir,
                            gene_min_cells = 10,
                            cell_min_genes = 50,
                            cell_max_genes = NULL,
                            min_count = 100,
                            max_mito_pct = NULL) {
  
  cat("\n========================================\n")
  cat("STEP 02: Quality Control\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  # Load Giotto object if path provided
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }
  
  # Create output directory
  results_folder <- file.path(output_dir, "02_QC")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Get initial stats
  n_cells_initial <- length(gobj@cell_ID$cell)
  n_genes_initial <- length(gobj@feat_ID$rna)
  
  cat("=== Initial Statistics ===\n")
  cat("Cells:", n_cells_initial, "\n")
  cat("Genes:", n_genes_initial, "\n\n")
  
  # Add statistics if not already present
  cat("Calculating cell/gene statistics...\n")
  
  gobj <- addStatistics(
    gobject = gobj,
    expression_values = "raw"
  )
  
  cat("✓ Statistics calculated\n\n")
  
  # Calculate mitochondrial percentage if MT genes present
  cat("Checking for mitochondrial genes...\n")
  
  mt_genes <- grep("^MT-|^Mt-|^mt-", fDataDT(gobj)$feat_ID, value = TRUE)
  
  if (length(mt_genes) > 0) {
    cat("Found", length(mt_genes), "mitochondrial genes\n")
    cat("  Examples:", paste(head(mt_genes, 5), collapse = ", "), "\n")
    
    gobj <- addFeatsPerc(
      gobject = gobj,
      feats = mt_genes,
      vector_name = "mito_pct"
    )
    
    cat("✓ Mitochondrial percentage calculated\n\n")
  } else {
    cat("⚠ No mitochondrial genes found (pattern: ^MT-)\n\n")
  }
  
  # Get cell metadata for QC plots
  metadata <- pDataDT(gobj) %>% as_tibble()
  
  cat("=== QC Metrics Summary ===\n")
  cat("Total counts per cell:\n")
  cat("  Min:", min(metadata$total_expr, na.rm = TRUE), "\n")
  cat("  Median:", median(metadata$total_expr, na.rm = TRUE), "\n")
  cat("  Max:", max(metadata$total_expr, na.rm = TRUE), "\n")
  cat("\nGenes per cell:\n")
  cat("  Min:", min(metadata$nr_feats, na.rm = TRUE), "\n")
  cat("  Median:", median(metadata$nr_feats, na.rm = TRUE), "\n")
  cat("  Max:", max(metadata$nr_feats, na.rm = TRUE), "\n")
  
  if ("mito_pct" %in% names(metadata)) {
    cat("\nMitochondrial %:\n")
    cat("  Median:", round(median(metadata$mito_pct, na.rm = TRUE), 2), "%\n")
    cat("  95th percentile:", round(quantile(metadata$mito_pct, 0.95, na.rm = TRUE), 2), "%\n")
  }
  cat("\n")
  
  # Create QC plots before filtering
  cat("Creating pre-filtering QC plots...\n")
  
  # Distribution plots
  p1 <- ggplot(metadata, aes(x = nr_feats)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    geom_vline(xintercept = cell_min_genes, color = "red", linetype = "dashed", linewidth = 1) +
    labs(title = "Genes per Cell (Pre-filtering)",
         x = "Number of Genes",
         y = "Cell Count") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(cell_max_genes)) {
    p1 <- p1 + geom_vline(xintercept = cell_max_genes, color = "red", linetype = "dashed", linewidth = 1)
  }
  
  p2 <- ggplot(metadata, aes(x = total_expr)) +
    geom_histogram(bins = 50, fill = "darkgreen", color = "white") +
    geom_vline(xintercept = min_count, color = "red", linetype = "dashed", linewidth = 1) +
    scale_x_log10() +
    labs(title = "Total Counts per Cell (Pre-filtering)",
         x = "Total Counts (log10)",
         y = "Cell Count") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if ("mito_pct" %in% names(metadata) && !is.null(max_mito_pct)) {
    p3 <- ggplot(metadata, aes(x = mito_pct)) +
      geom_histogram(bins = 50, fill = "coral", color = "white") +
      geom_vline(xintercept = max_mito_pct, color = "red", linetype = "dashed", linewidth = 1) +
      labs(title = "Mitochondrial % (Pre-filtering)",
           x = "Mitochondrial %",
           y = "Cell Count") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    combined <- (p1 | p2 | p3) +
      plot_annotation(
        title = paste(sample_id, "- QC Metrics Before Filtering"),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      )
  } else {
    combined <- (p1 | p2) +
      plot_annotation(
        title = paste(sample_id, "- QC Metrics Before Filtering"),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
      )
  }
  
  ggsave(
    filename = file.path(results_folder, paste0(sample_id, "_qc_pre_filtering.png")),
    plot = combined,
    width = 18,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✓ Pre-filtering plots saved\n\n")
  
  # Spatial QC plots
  cat("Creating spatial QC plots...\n")
  
  tryCatch({
    spatPlot2D(
      gobject = gobj,
      show_image = FALSE,
      point_alpha = 0.7,
      point_size = 0.5,
      cell_color = "nr_feats",
      color_as_factor = FALSE,
      gradient_style = "sequential",
      save_plot = TRUE,
      save_param = list(
        save_name = paste0(sample_id, "_spatial_genes_per_cell"),
        save_dir = results_folder,
        base_width = 12,
        base_height = 10
      )
    )
    
    spatPlot2D(
      gobject = gobj,
      show_image = FALSE,
      point_alpha = 0.7,
      point_size = 0.5,
      cell_color = "total_expr",
      color_as_factor = FALSE,
      gradient_style = "sequential",
      save_plot = TRUE,
      save_param = list(
        save_name = paste0(sample_id, "_spatial_total_counts"),
        save_dir = results_folder,
        base_width = 12,
        base_height = 10
      )
    )
    
    if ("mito_pct" %in% names(metadata)) {
      spatPlot2D(
        gobject = gobj,
        show_image = FALSE,
        point_alpha = 0.7,
        point_size = 0.5,
        cell_color = "mito_pct",
        color_as_factor = FALSE,
        gradient_style = "sequential",
        save_plot = TRUE,
        save_param = list(
          save_name = paste0(sample_id, "_spatial_mito_pct"),
          save_dir = results_folder,
          base_width = 12,
          base_height = 10
        )
      )
    }
    
    cat("✓ Spatial plots saved\n\n")
  }, error = function(e) {
    cat("⚠ Spatial plotting warning:", conditionMessage(e), "\n\n")
  })
  
  # Filter genes
  cat("Filtering genes...\n")
  cat("  Minimum cells per gene:", gene_min_cells, "\n")
  
  gobj <- filterGiotto(
    gobject = gobj,
    feat_type = "rna",
    expression_threshold = 1,
    feat_det_in_min_cells = gene_min_cells,
    min_det_feats_per_cell = 0  # Don't filter cells yet
  )
  
  n_genes_after <- length(gobj@feat_ID$rna)
  cat("✓ Genes after filtering:", n_genes_after, 
      "(removed", n_genes_initial - n_genes_after, ")\n\n")
  
  # Filter cells
  cat("Filtering cells...\n")
  cat("  Minimum genes per cell:", cell_min_genes, "\n")
  if (!is.null(cell_max_genes)) {
    cat("  Maximum genes per cell:", cell_max_genes, "\n")
  }
  cat("  Minimum total counts:", min_count, "\n")
  if (!is.null(max_mito_pct)) {
    cat("  Maximum mitochondrial %:", max_mito_pct, "\n")
  }
  
  # Build cell filter criteria
  cell_metadata <- pDataDT(gobj)
  
  cells_to_keep <- cell_metadata$cell_ID[
    cell_metadata$nr_feats >= cell_min_genes &
      cell_metadata$total_expr >= min_count
  ]
  
  if (!is.null(cell_max_genes)) {
    cells_to_keep <- intersect(
      cells_to_keep,
      cell_metadata$cell_ID[cell_metadata$nr_feats <= cell_max_genes]
    )
  }
  
  if (!is.null(max_mito_pct) && "mito_pct" %in% names(cell_metadata)) {
    cells_to_keep <- intersect(
      cells_to_keep,
      cell_metadata$cell_ID[cell_metadata$mito_pct <= max_mito_pct]
    )
  }
  
  # Apply cell filter
  gobj <- subsetGiotto(
    gobject = gobj,
    cell_ids = cells_to_keep
  )
  
  n_cells_after <- length(gobj@cell_ID$cell)
  cat("✓ Cells after filtering:", n_cells_after,
      "(removed", n_cells_initial - n_cells_after, ")\n\n")
  
  # Post-filtering plots
  cat("Creating post-filtering QC plots...\n")
  
  metadata_post <- pDataDT(gobj) %>% as_tibble()
  
  p1_post <- ggplot(metadata_post, aes(x = nr_feats)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    labs(title = "Genes per Cell (Post-filtering)",
         x = "Number of Genes",
         y = "Cell Count") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  p2_post <- ggplot(metadata_post, aes(x = total_expr)) +
    geom_histogram(bins = 50, fill = "darkgreen", color = "white") +
    scale_x_log10() +
    labs(title = "Total Counts per Cell (Post-filtering)",
         x = "Total Counts (log10)",
         y = "Cell Count") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  combined_post <- (p1_post | p2_post) +
    plot_annotation(
      title = paste(sample_id, "- QC Metrics After Filtering"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    )
  
  ggsave(
    filename = file.path(results_folder, paste0(sample_id, "_qc_post_filtering.png")),
    plot = combined_post,
    width = 14,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✓ Post-filtering plots saved\n\n")
  
  # Save QC summary
  qc_summary <- tibble(
    metric = c("Cells (initial)", "Cells (final)", "Cells (removed)",
               "Genes (initial)", "Genes (final)", "Genes (removed)",
               "Median genes/cell", "Median counts/cell"),
    value = c(n_cells_initial, n_cells_after, n_cells_initial - n_cells_after,
              n_genes_initial, n_genes_after, n_genes_initial - n_genes_after,
              median(metadata_post$nr_feats),
              median(metadata_post$total_expr))
  )
  
  write_csv(qc_summary, file.path(results_folder, paste0(sample_id, "_qc_summary.csv")))
  
  cat("=== Final Statistics ===\n")
  print(qc_summary)
  cat("\n")
  
  cat("✓ Quality control complete for", sample_id, "\n\n")
  
  return(gobj)
}

# Run if sourced directly
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) {
    sample_id <- args[1]
    input_path <- args[2]
    output_dir <- args[3]
    
    gobj <- quality_control(
      gobj = input_path,
      sample_id = sample_id,
      output_dir = output_dir
    )
    
    # Save filtered object
    saveGiotto(
      gobject = gobj,
      dir = output_dir,
      foldername = "Giotto_Object_Filtered",
      overwrite = TRUE
    )
  }
}

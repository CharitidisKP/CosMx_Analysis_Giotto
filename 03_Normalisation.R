# Normalize expression data and add statistics ----------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 03_normalization.R
# Normalize expression data
# ==============================================================================

#' Normalize Expression Data
#'
#' @param gobj Giotto object or path to saved object
#' @param sample_id Sample identifier
#' @param output_dir Output directory
#' @param scalefactor Scale factor for normalization (default: 6000)
#' @param log_transform Apply log transformation
#' @return Normalized Giotto object

normalize_expression <- function(gobj,
                                 sample_id,
                                 output_dir,
                                 scalefactor = 6000,
                                 log_transform = TRUE) {
  
  cat("\n========================================\n")
  cat("STEP 03: Normalization\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  # Load if path
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }
  
  results_folder <- file.path(output_dir, "03_Normalization")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  cat("Normalizing expression data...\n")
  cat("  Scale factor:", scalefactor, "\n")
  cat("  Log transform:", log_transform, "\n\n")
  
  # Normalize
  gobj <- normalizeGiotto(
    gobject = gobj,
    spat_unit = "cell",
    feat_type = "rna",
    scalefactor = scalefactor,
    verbose = TRUE,
    norm_methods = "standard",
    library_size_norm = TRUE,
    log_norm = log_transform,
    log_offset = 1,
    scale_feats = FALSE,
    scale_cells = FALSE
  )
  
  cat("\n✓ Normalization complete\n\n")
  
  # Check normalized data
  cat("=== Normalization Summary ===\n")
  cat("Expression slots available:\n")
  print(names(gobj@expression$cell$rna))
  
  # Compare raw vs normalized
  if ("normalized" %in% names(gobj@expression$cell$rna)) {
    
    raw_expr <- gobj@expression$cell$rna$raw[]
    norm_expr <- gobj@expression$cell$rna$normalized[]
    
    cat("\nRaw expression:\n")
    cat("  Mean:", round(mean(raw_expr), 2), "\n")
    cat("  Median:", round(median(raw_expr), 2), "\n")
    cat("  Max:", round(max(raw_expr), 2), "\n")
    
    cat("\nNormalized expression:\n")
    cat("  Mean:", round(mean(norm_expr), 2), "\n")
    cat("  Median:", round(median(norm_expr), 2), "\n")
    cat("  Max:", round(max(norm_expr), 2), "\n")
    
    # Create comparison plots
    comparison_data <- tibble(
      raw_mean = Matrix::rowMeans(raw_expr),
      norm_mean = Matrix::rowMeans(norm_expr)
    )
    
    p <- ggplot(comparison_data, aes(x = raw_mean, y = norm_mean)) +
      geom_point(alpha = 0.3, size = 0.5) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      scale_x_log10() +
      scale_y_log10() +
      labs(title = paste(sample_id, "- Raw vs Normalized Expression"),
           x = "Raw Mean Expression (log10)",
           y = "Normalized Mean Expression (log10)") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(
      filename = file.path(results_folder, paste0(sample_id, "_raw_vs_normalized.png")),
      plot = p,
      width = 8,
      height = 8,
      dpi = 300,
      bg = "white"
    )
    
    cat("\n✓ Comparison plots saved\n")
  }
  
  cat("\n✓ Normalization complete for", sample_id, "\n\n")
  
  return(gobj)
}

# Run if sourced directly
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2) {
    sample_id <- args[1]
    input_path <- args[2]
    output_dir <- args[3]
    
    gobj <- normalize_expression(
      gobj = input_path,
      sample_id = sample_id,
      output_dir = output_dir
    )
    
    saveGiotto(gobj, dir = output_dir, foldername = "Giotto_Object_Normalized", overwrite = TRUE)
  }
}

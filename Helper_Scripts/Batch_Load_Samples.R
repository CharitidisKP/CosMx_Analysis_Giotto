#!/usr/bin/env Rscript
# ==============================================================================
# batch_load_samples.R
# Batch load CosMx samples using smart sample detection
# ==============================================================================

# Setup
setwd("~/P_lab/CosMx_analysis")
source("Scripts/00_Setup.R")
source("Scripts/01_Load_data.R")
source("Scripts/Helper_Scripts/Get_samples.R")

# Get samples
# Uses CSV if exists, auto-detects otherwise
samples <- get_samples(
  method = "auto",  # or "csv" or "detect"
  # filter_ids = c("18H12037", "1281")  # Optional: process only these samples
)

cat("=== Samples to Process ===\n")
print(samples %>% select(sample_id, slide_num, directory_name))
cat("\n")

# Track results
results <- tibble(
  sample_id = character(),
  slide_num = character(),
  status = character(),
  cells = integer(),
  genes = integer(),
  polygons = integer(),
  time_seconds = numeric(),
  error_msg = character()
)

# Process each sample
cat("\n")
cat("================================================================================\n")
cat("                    BATCH LOADING CosMx SAMPLES                                \n")
cat("================================================================================\n\n")

for (i in seq_len(nrow(samples))) {
  
  sample <- samples[i, ]
  
  cat("################################################################################\n")
  cat("## SAMPLE", i, "of", nrow(samples), ":", sample$sample_id, "(Slide", sample$slide_num, ")\n")
  cat("################################################################################\n\n")
  
  start_time <- Sys.time()
  
  result <- tryCatch({
    
    # Load sample
    cosmx <- load_cosmx_sample(
      sample_id = sample$sample_id,
      data_dir = sample$data_dir,
      output_dir = sample$output_dir
    )
    
    # Get counts
    n_cells <- length(cosmx@cell_ID$cell)
    n_genes <- length(cosmx@feat_ID$rna)
    n_poly <- if ("cell" %in% names(cosmx@spatial_info)) {
      nrow(cosmx@spatial_info$cell)
    } else {
      0
    }
    
    # Save RDS
    output_file <- file.path(sample$output_dir, paste0(sample$sample_id, "_cosmx_loaded.rds"))
    saveRDS(cosmx, output_file)
    cat("\n✓ Saved to:", output_file, "\n")
    
    list(
      status = "SUCCESS",
      cells = n_cells,
      genes = n_genes,
      polygons = n_poly,
      error = NA_character_
    )
    
  }, error = function(e) {
    cat("\n✗ ERROR:", conditionMessage(e), "\n")
    list(
      status = "FAILED",
      cells = NA_integer_,
      genes = NA_integer_,
      polygons = NA_integer_,
      error = conditionMessage(e)
    )
  })
  
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Record result
  results <- bind_rows(
    results,
    tibble(
      sample_id = sample$sample_id,
      slide_num = sample$slide_num,
      status = result$status,
      cells = result$cells,
      genes = result$genes,
      polygons = result$polygons,
      time_seconds = elapsed,
      error_msg = result$error
    )
  )
  
  if (result$status == "SUCCESS") {
    cat("\n✓✓✓ Sample", sample$sample_id, "COMPLETE ✓✓✓\n")
    cat("    Cells:", result$cells, "| Genes:", result$genes, 
        "| Polygons:", result$polygons, "| Time:", round(elapsed, 1), "sec\n\n")
  } else {
    cat("\n✗✗✗ Sample", sample$sample_id, "FAILED ✗✗✗\n\n")
  }
  
  gc()
}

# Save results
results_file <- file.path("~/P_lab/CosMx_analysis/Output", 
                          paste0("batch_load_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"))
write_csv(results, results_file)

# Summary
cat("\n")
cat("================================================================================\n")
cat("                         BATCH LOADING SUMMARY                                  \n")
cat("================================================================================\n\n")

print(results)

cat("\nSuccess rate:", sum(results$status == "SUCCESS"), "/", nrow(results), "\n")
cat("Total time:", round(sum(results$time_seconds, na.rm = TRUE) / 60, 1), "minutes\n")
cat("Average time per sample:", round(mean(results$time_seconds, na.rm = TRUE), 1), "seconds\n")
cat("\nResults saved to:", results_file, "\n\n")

cat("================================================================================\n")
cat("                           BATCH COMPLETE                                       \n")
cat("================================================================================\n\n")
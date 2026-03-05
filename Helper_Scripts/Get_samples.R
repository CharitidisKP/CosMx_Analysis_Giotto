#!/usr/bin/env Rscript
# ==============================================================================
# Get_samples.R
# Smart sample loader: uses CSV if exists, auto-detects otherwise
# ==============================================================================

library(tidyverse)

source("~/P_lab/CosMx_analysis/Scripts/auto_detect_samples.R")

#' Get samples using best available method
#'
#' @param method "auto", "csv", or "detect"
#' @param csv_path Path to sample registry CSV (if using CSV)
#' @param filter_ids Optional vector of sample IDs to filter
#' @return Validated tibble with sample information
get_samples <- function(method = "auto",
                        csv_path = "~/P_lab/CosMx_analysis/Data/sample_registry.csv",
                        filter_ids = NULL) {
  
  cat("\n=== Getting Sample Information ===\n\n")
  
  # Determine method
  if (method == "auto") {
    if (file.exists(csv_path)) {
      cat("✓ Found sample registry CSV, using that\n")
      method <- "csv"
    } else {
      cat("⚠ No sample registry found, auto-detecting\n")
      method <- "detect"
    }
  }
  
  # Load samples
  if (method == "csv") {
    samples <- load_samples_from_csv(csv_path)
  } else {
    samples <- auto_detect_samples()
  }
  
  # Validate
  samples <- validate_samples(samples)
  
  # Filter to valid samples only
  valid_samples <- samples %>% filter(is_valid)
  
  if (nrow(valid_samples) == 0) {
    stop("No valid samples found!")
  }
  
  # Apply ID filter if provided
  if (!is.null(filter_ids)) {
    cat("Filtering to specified sample IDs:", paste(filter_ids, collapse = ", "), "\n")
    valid_samples <- valid_samples %>%
      filter(sample_id %in% filter_ids)
    
    if (nrow(valid_samples) == 0) {
      stop("No samples match the filter criteria")
    }
  }
  
  cat("\n✓ Ready to process", nrow(valid_samples), "samples\n\n")
  
  return(valid_samples)
}
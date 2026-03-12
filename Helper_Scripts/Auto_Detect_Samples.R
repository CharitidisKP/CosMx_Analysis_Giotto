#!/usr/bin/env Rscript
# ==============================================================================
# auto_detect_samples.R
# Automatically detect and parse CosMx sample information
# ==============================================================================

library(tidyverse)

default_project_dir <- function() {
  project_dir <- Sys.getenv("COSMX_PROJECT_DIR", unset = "")
  if (nzchar(project_dir)) {
    return(path.expand(project_dir))
  }
  getwd()
}

#' Auto-detect samples from raw data directory
#'
#' @param raw_data_dir Path to raw data directory
#' @return Tibble with sample information

auto_detect_samples <- function(raw_data_dir = file.path(default_project_dir(), "Data", "Raw_data")) {
  
  cat("Scanning directory:", raw_data_dir, "\n")
  
  # List all directories
  dirs <- list.dirs(raw_data_dir, recursive = FALSE, full.names = FALSE)
  
  if (length(dirs) == 0) {
    stop("No sample directories found in ", raw_data_dir)
  }
  
  cat("Found", length(dirs), "directories\n\n")
  
  # Parse each directory name
  samples <- tibble(directory_name = dirs) %>%
    mutate(
      # Extract slide number
      # Pattern: Slide followed by 1-2 digits, then sample ID
      slide_num = str_extract(directory_name, "(?<=Slide)\\d{1,2}"),
      
      # Extract sample ID (everything after slide number)
      # Pattern: after "Slide[1-2 digits]", get remaining alphanumeric
      sample_id_raw = str_extract(directory_name, "(?<=Slide\\d{1,2}).*$"),
      
      # Clean sample ID (remove leading underscores)
      sample_id = str_remove(sample_id_raw, "^_+"),
      
      # Create full paths
      data_dir = file.path(raw_data_dir, directory_name),
      output_dir = file.path(default_project_dir(), "Output", paste0("Sample_", sample_id)),
      
      # Detect file type prefix (everything before "Slide")
      prefix = str_extract(directory_name, "^.*(?=Slide)")
    ) %>%
    select(sample_id, slide_num, directory_name, data_dir, output_dir, prefix)
  
  return(samples)
}

#' Validate detected samples
#'
#' @param samples Tibble from auto_detect_samples
#' @return Validated samples with file checks
validate_samples <- function(samples) {
  
  cat("Validating samples...\n\n")
  
  samples <- samples %>%
    rowwise() %>%
    mutate(
      # Check if directory exists
      dir_exists = dir.exists(data_dir),
      
      # Check for required files
      has_expr = if (dir_exists) {
        length(list.files(data_dir, pattern = "_exprMat_file\\.csv", full.names = TRUE)) > 0
      } else {
        FALSE
      },
      has_meta = if (dir_exists) {
        length(list.files(data_dir, pattern = "_metadata_file\\.csv", full.names = TRUE)) > 0
      } else {
        FALSE
      },
      has_poly = if (dir_exists) {
        length(list.files(data_dir, pattern = "-polygons\\.csv", full.names = TRUE)) > 0
      } else {
        FALSE
      },
      has_fov = if (dir_exists) {
        length(list.files(data_dir, pattern = "_fov_positions_file\\.csv", full.names = TRUE)) > 0
      } else {
        FALSE
      },
      
      # Overall validity
      is_valid = dir_exists & has_expr & has_meta & has_poly & has_fov
    ) %>%
    ungroup()
  
  # Print validation summary
  cat("=== Validation Summary ===\n")
  cat("Total samples:", nrow(samples), "\n")
  cat("Valid samples:", sum(samples$is_valid), "\n")
  cat("Invalid samples:", sum(!samples$is_valid), "\n\n")
  
  if (any(!samples$is_valid)) {
    cat("Invalid samples:\n")
    invalid <- samples %>% 
      filter(!is_valid) %>%
      select(sample_id, directory_name, dir_exists, has_expr, has_meta, has_poly, has_fov)
    print(invalid)
    cat("\n")
  }
  
  return(samples)
}

#' Load samples from CSV registry
#'
#' @param csv_path Path to sample registry CSV
#' @param raw_data_dir Base directory for raw data
#' @return Tibble with sample information
load_samples_from_csv <- function(csv_path = file.path(default_project_dir(), "Data", "sample_registry.csv"),
                                  raw_data_dir = file.path(default_project_dir(), "Data", "Raw_data")) {
  
  if (!file.exists(csv_path)) {
    stop("Sample registry CSV not found: ", csv_path)
  }
  
  cat("Loading samples from:", csv_path, "\n")
  
  samples <- read_csv(csv_path, show_col_types = FALSE) %>%
    mutate(
      directory_name = Filename,
      sample_id = Sample_ID,
      slide_num = as.character(Slide_num),
      data_dir = file.path(raw_data_dir, directory_name),
      output_dir = file.path(default_project_dir(), "Output", paste0("Sample_", sample_id))
    ) %>%
    select(sample_id, slide_num, directory_name, data_dir, output_dir, Notes)
  
  cat("✓ Loaded", nrow(samples), "samples from CSV\n\n")
  
  return(samples)
}

#' Create sample registry CSV from auto-detection
#'
#' @param output_path Where to save the CSV
create_sample_registry <- function(output_path = file.path(default_project_dir(), "Data", "sample_registry.csv")) {
  
  cat("Creating sample registry...\n\n")
  
  # Auto-detect samples
  samples <- auto_detect_samples()
  samples <- validate_samples(samples)
  
  # Create registry format
  registry <- samples %>%
    select(
      Filename = directory_name,
      Sample_ID = sample_id,
      Slide_num = slide_num
    ) %>%
    mutate(
      Notes = "",
      Valid = samples$is_valid
    )
  
  # Save
  write_csv(registry, output_path)
  
  cat("✓ Sample registry saved to:", output_path, "\n")
  cat("  Review and edit this file, then use it for batch loading\n\n")
  
  # Print preview
  cat("=== Registry Preview ===\n")
  print(registry)
  
  return(registry)
}

# If run directly, create registry
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  registry <- create_sample_registry()
}

# Main wrapper script to run complete CosMx analysis pipeline -------------

#!/usr/bin/env Rscript
# ==============================================================================
# Combines automatic sample detection with flexible YAML configuration
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(optparse)
})

#' Load and validate configuration
#' @param config_path Path to YAML configuration file
#' @return Configuration list
load_config <- function(config_path) {
  
  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }
  
  cat("Loading configuration from:", config_path, "\n")
  config <- yaml::read_yaml(config_path)
  
  # Validate required paths
  if (!dir.exists(config$paths$main_dir)) {
    stop("Main directory not found: ", config$paths$main_dir)
  }
  
  if (!dir.exists(config$paths$raw_data_dir)) {
    stop("Raw data directory not found: ", config$paths$raw_data_dir)
  }
  
  cat("✓ Configuration validated\n\n")
  
  return(config)
}

#' Setup environment and source all scripts
#' @param config Configuration list
setup_environment <- function(config) {
  
  cat("Setting up environment...\n")
  
  # Set working directory
  setwd(config$paths$main_dir)
  
  # Source environment setup
  env_script <- file.path(config$paths$scripts_dir, "00_setup_environment.R")
  if (!file.exists(env_script)) {
    stop("Environment setup script not found: ", env_script)
  }
  source(env_script)
  
  # Source all pipeline scripts
  scripts <- c(
    "01_load_data.R",
    "02_quality_control.R",
    "03_normalization.R",
    "04_dimensionality_reduction.R",
    "05_clustering.R",
    "06_marker_analysis.R",
    "07_annotation.R"
  )
  
  for (script in scripts) {
    script_path <- file.path(config$paths$scripts_dir, script)
    if (!file.exists(script_path)) {
      warning("Script not found: ", script_path)
    } else {
      source(script_path)
    }
  }
  
  # Source helper scripts
  helper_script <- file.path(config$paths$scripts_dir, "Helper_Scripts/auto_detect_samples.R")
  if (file.exists(helper_script)) {
    source(helper_script)
  }
  
  cat("✓ Environment setup complete\n\n")
}

#' Get samples to process
#' @param config Configuration list
#' @param samples_override Optional sample IDs to override config
#' @return Tibble of samples to process
get_samples_to_process <- function(config, samples_override = NULL) {
  
  cat("Determining samples to process...\n")
  
  # Determine sample IDs
  if (!is.null(samples_override)) {
    sample_ids <- samples_override
    cat("Using command-line specified samples:", paste(sample_ids, collapse = ", "), "\n")
  } else if (config$samples$sample_ids != "all") {
    sample_ids <- config$samples$sample_ids
    cat("Using config specified samples:", paste(sample_ids, collapse = ", "), "\n")
  } else {
    sample_ids <- NULL
    cat("Processing all available samples\n")
  }
  
  # Get and validate samples
  samples <- get_samples_smart(
    csv_path = config$paths$sample_registry,
    filter_ids = sample_ids
  )
  
  cat("✓ Found", nrow(samples), "samples to process\n\n")
  
  return(samples)
}

#' Process a single sample through the pipeline
#' @param sample Sample information (single row tibble)
#' @param config Configuration list
#' @return Result row with status for each step
process_sample <- function(sample, config) {
  
  sample_id <- sample$sample_id
  
  cat("\n")
  cat("################################################################################\n")
  cat("## PROCESSING SAMPLE:", sample_id, "\n")
  cat("################################################################################\n\n")
  
  sample_start <- Sys.time()
  
  # Initialize result row
  result_row <- tibble(
    sample_id = sample_id,
    step_01_load = "SKIP",
    step_02_qc = "SKIP",
    step_03_norm = "SKIP",
    step_04_dimred = "SKIP",
    step_05_cluster = "SKIP",
    step_06_markers = "SKIP",
    step_07_annotate = "SKIP",
    total_time_min = NA_real_,
    final_cells = NA_integer_,
    n_clusters = NA_integer_
  )
  
  # Track Giotto object through pipeline
  gobj <- NULL
  
  tryCatch({
    
    # STEP 01: Load Data
    if (config$pipeline$steps$step_01_load) {
      cat(">>> STEP 01: LOADING DATA <<<\n\n")
      
      gobj <- load_cosmx_sample(
        sample_id = sample_id,
        data_dir = sample$data_dir,
        output_dir = sample$output_dir
      )
      
      result_row$step_01_load <- "SUCCESS"
      
      if (config$pipeline$save_intermediates) {
        saveGiotto(gobj, sample$output_dir, "Giotto_Step01_Loaded", overwrite = TRUE)
      }
    } else {
      # Load existing
      gobj <- loadGiotto(file.path(sample$output_dir, "Giotto_Object"))
    }
    
    # STEP 02: Quality Control
    if (config$pipeline$steps$step_02_qc) {
      cat("\n>>> STEP 02: QUALITY CONTROL <<<\n\n")
      
      gobj <- quality_control(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = sample$output_dir,
        gene_min_cells = config$parameters$qc$gene_min_cells,
        cell_min_genes = config$parameters$qc$cell_min_genes,
        cell_max_genes = config$parameters$qc$cell_max_genes,
        min_count = config$parameters$qc$min_count,
        max_mito_pct = config$parameters$qc$max_mito_pct
      )
      
      result_row$step_02_qc <- "SUCCESS"
      
      if (config$pipeline$save_intermediates) {
        saveGiotto(gobj, sample$output_dir, "Giotto_Step02_QC", overwrite = TRUE)
      }
    }
    
    # STEP 03: Normalization
    if (config$pipeline$steps$step_03_norm) {
      cat("\n>>> STEP 03: NORMALIZATION <<<\n\n")
      
      gobj <- normalize_expression(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = sample$output_dir,
        scalefactor = config$parameters$normalization$scalefactor,
        log_transform = config$parameters$normalization$log_transform
      )
      
      result_row$step_03_norm <- "SUCCESS"
      
      if (config$pipeline$save_intermediates) {
        saveGiotto(gobj, sample$output_dir, "Giotto_Step03_Normalized", overwrite = TRUE)
      }
    }
    
    # STEP 04: Dimensionality Reduction
    if (config$pipeline$steps$step_04_dimred) {
      cat("\n>>> STEP 04: DIMENSIONALITY REDUCTION <<<\n\n")
      
      gobj <- dimensionality_reduction(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = sample$output_dir,
        n_hvgs = config$parameters$dim_reduction$n_hvgs,
        n_pcs = config$parameters$dim_reduction$n_pcs,
        umap_n_neighbors = config$parameters$dim_reduction$umap_n_neighbors,
        umap_min_dist = config$parameters$dim_reduction$umap_min_dist,
        spatial_hvg = config$parameters$dim_reduction$spatial_hvg
      )
      
      result_row$step_04_dimred <- "SUCCESS"
      
      if (config$pipeline$save_intermediates) {
        saveGiotto(gobj, sample$output_dir, "Giotto_Step04_DimReduced", overwrite = TRUE)
      }
    }
    
    # STEP 05: Clustering
    if (config$pipeline$steps$step_05_cluster) {
      cat("\n>>> STEP 05: CLUSTERING <<<\n\n")
      
      # Parse dimensions_to_use
      dims_to_use <- if (is.character(config$parameters$clustering$dimensions_to_use)) {
        eval(parse(text = config$parameters$clustering$dimensions_to_use))
      } else {
        config$parameters$clustering$dimensions_to_use
      }
      
      gobj <- perform_clustering(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = sample$output_dir,
        k_nn = config$parameters$clustering$k_nn,
        resolution = config$parameters$clustering$resolution,
        dimensions_to_use = dims_to_use
      )
      
      result_row$step_05_cluster <- "SUCCESS"
      
      if (config$pipeline$save_intermediates) {
        saveGiotto(gobj, sample$output_dir, "Giotto_Step05_Clustered", overwrite = TRUE)
      }
    }
    
    # STEP 06: Marker Analysis
    if (config$pipeline$steps$step_06_markers) {
      cat("\n>>> STEP 06: MARKER ANALYSIS <<<\n\n")
      
      gobj <- marker_analysis(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = sample$output_dir,
        cluster_column = config$parameters$markers$cluster_column,
        top_n = config$parameters$markers$top_n
      )
      
      result_row$step_06_markers <- "SUCCESS"
      
      if (config$pipeline$save_intermediates) {
        saveGiotto(gobj, sample$output_dir, "Giotto_Step06_WithMarkers", overwrite = TRUE)
      }
    }
    
    # STEP 07: Annotation
    if (config$pipeline$steps$step_07_annotate) {
      cat("\n>>> STEP 07: ANNOTATION <<<\n\n")
      
      gobj <- annotate_cells(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = sample$output_dir,
        profiles = config$parameters$annotation$profiles,
        default_profile = config$parameters$annotation$default_profile,
        align_genes = config$parameters$annotation$align_genes,
        n_starts = config$parameters$annotation$n_starts,
        n_clusts_semi = config$parameters$annotation$n_clusts_semi,
        cohort_column = config$parameters$annotation$cohort_column,
        min_gene_overlap = config$parameters$annotation$min_gene_overlap,
        create_plots = config$parameters$annotation$create_flightplots
      )
      
      result_row$step_07_annotate <- "SUCCESS"
      
      if (config$pipeline$save_intermediates) {
        saveGiotto(gobj, sample$output_dir, "Giotto_Step07_Annotated", overwrite = TRUE)
      }
    }
    
    # Save final object
    saveGiotto(gobj, sample$output_dir, "Giotto_Object_Final", overwrite = TRUE)
    
    # Collect final stats
    result_row$final_cells <- length(gobj@cell_ID$cell)
    
    if ("leiden_clus" %in% names(pDataDT(gobj))) {
      result_row$n_clusters <- length(unique(pDataDT(gobj)$leiden_clus))
    }
    
    sample_end <- Sys.time()
    result_row$total_time_min <- as.numeric(difftime(sample_end, sample_start, units = "mins"))
    
    cat("\n✓✓✓ SAMPLE", sample_id, "COMPLETE ✓✓✓\n")
    cat("    Time:", round(result_row$total_time_min, 1), "minutes\n\n")
    
  }, error = function(e) {
    cat("\n✗✗✗ SAMPLE", sample_id, "FAILED ✗✗✗\n")
    cat("    Error:", conditionMessage(e), "\n\n")
    
    # Mark all non-SUCCESS steps as FAILED
    for (col in grep("^step_", names(result_row), value = TRUE)) {
      if (result_row[[col]] == "SKIP") {
        result_row[[col]] <- "FAILED"
      }
    }
  })
  
  return(result_row)
}

#' Main pipeline function
#' @param config_path Path to configuration YAML
#' @param samples_override Optional sample IDs
#' @param steps_override Optional steps to run
#' @return Pipeline results tibble
run_pipeline <- function(config_path = "~/P_lab/CosMx_analysis/config/pipeline_config.yaml",
                         samples_override = NULL,
                         steps_override = NULL) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("                      CosMx ANALYSIS PIPELINE                                   \n")
  cat("================================================================================\n\n")
  
  # Load configuration
  config <- load_config(config_path)
  
  # Setup environment
  setup_environment(config)
  
  # Get samples
  samples <- get_samples_to_process(config, samples_override)
  
  # Initialize results tracking
  pipeline_results <- tibble(
    sample_id = character(),
    step_01_load = character(),
    step_02_qc = character(),
    step_03_norm = character(),
    step_04_dimred = character(),
    step_05_cluster = character(),
    step_06_markers = character(),
    step_07_annotate = character(),
    total_time_min = numeric(),
    final_cells = integer(),
    n_clusters = integer()
  )
  
  # Process each sample
  for (i in seq_len(nrow(samples))) {
    sample <- samples[i, ]
    result_row <- process_sample(sample, config)
    pipeline_results <- bind_rows(pipeline_results, result_row)
    gc()
  }
  
  # Save results
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_file <- file.path(config$paths$output_dir, 
                            paste0("pipeline_results_", timestamp, ".csv"))
  write_csv(pipeline_results, results_file)
  
  # Print summary
  cat("\n")
  cat("================================================================================\n")
  cat("                         PIPELINE SUMMARY                                       \n")
  cat("================================================================================\n\n")
  
  print(pipeline_results)
  
  cat("\nResults saved to:", results_file, "\n\n")
  
  cat("================================================================================\n")
  cat("                         PIPELINE COMPLETE                                      \n")
  cat("================================================================================\n\n")
  
  return(invisible(pipeline_results))
}

# Command-line interface
if (!interactive()) {
  
  option_list <- list(
    make_option(c("-c", "--config"), type = "character",
                default = "~/P_lab/CosMx_analysis/config/pipeline_config.yaml",
                help = "Path to configuration YAML [default: %default]"),
    make_option(c("-s", "--samples"), type = "character", default = NULL,
                help = "Comma-separated sample IDs [default: from config]"),
    make_option(c("--steps"), type = "character", default = NULL,
                help = "Steps to run (not implemented yet) [default: from config]")
  )
  
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  args <- parse_args(parser)
  
  # Parse samples
  samples_override <- if (!is.null(args$samples)) {
    strsplit(args$samples, ",")[[1]]
  } else {
    NULL
  }
  
  # Run pipeline
  results <- run_pipeline(
    config_path = args$config,
    samples_override = samples_override
  )
  
  # Exit with status
  if (any(grepl("FAILED", unlist(results[, grep("^step_", names(results))])))) {
    quit(status = 1)
  } else {
    quit(status = 0)
  }
}
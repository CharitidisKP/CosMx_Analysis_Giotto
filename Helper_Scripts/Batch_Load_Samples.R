#!/usr/bin/env Rscript
# ==============================================================================
# batch_load_samples.R
# Smart batch loading with automatic registry creation and validation
# ==============================================================================

# Setup
setwd("~/P_lab/CosMx_analysis")
source("Scripts/00_Setup.R")
source("Scripts/01_Load_data.R")

# Source helper functions
source("Scripts/Helper_Scripts/Auto_Detect_Samples.R")

#' Smart sample getter - creates registry if needed
#'
#' @param csv_path Path to registry CSV
#' @param force_detect Force re-detection even if CSV exists
#' @param filter_ids Optional sample IDs to process
#' @return Validated sample information

get_samples_smart <- function(csv_path = "~/P_lab/CosMx_analysis/Data/sample_registry.csv",
                              force_detect = FALSE,
                              filter_ids = NULL) {
  
  cat("\n=== Loading Sample Information ===\n\n")
  
  if (file.exists(csv_path) && !force_detect) {
    cat("✓ Using existing registry:", csv_path, "\n")
    samples <- load_samples_from_csv(csv_path)
  } else {
    if (force_detect) {
      cat("⚠ Force detect enabled\n")
    } else {
      cat("⚠ No registry found, auto-detecting\n")
    }
    
    samples <- auto_detect_samples()
    
    registry <- samples %>%
      select(Filename = directory_name, Sample_ID = sample_id, Slide_num = slide_num) %>%
      mutate(Notes = "")
    
    dir.create(dirname(csv_path), recursive = TRUE, showWarnings = FALSE)
    write_csv(registry, csv_path)
    cat("✓ Created registry:", csv_path, "\n")
  }
  
  samples <- validate_samples(samples)
  valid_samples <- samples %>% filter(is_valid)
  
  if (nrow(valid_samples) == 0) {
    stop("No valid samples found!")
  }
  
  if (!is.null(filter_ids)) {
    cat("Filtering to:", paste(filter_ids, collapse = ", "), "\n")
    valid_samples <- valid_samples %>%
      filter(sample_id %in% filter_ids)
    
    if (nrow(valid_samples) == 0) {
      stop("No samples match filter")
    }
  }
  
  cat("✓ Ready to process", nrow(valid_samples), "samples\n\n")
  
  return(valid_samples)
}

#' Save Giotto object safely (avoids C stack issues)
#'
#' @param gobj Giotto object
#' @param output_dir Output directory
#' @param sample_id Sample identifier
#' @return TRUE if successful
save_giotto_safe <- function(gobj, output_dir, sample_id) {
  
  cat("Saving Giotto object...\n")
  
  # Try saveGiotto first (Giotto's native format)
  tryCatch({
    saveGiotto(
      gobject = gobj,
      dir = output_dir,
      foldername = "Giotto_Object",
      overwrite = TRUE
    )
    cat("✓ Saved using saveGiotto\n")
    return(TRUE)
  }, error = function(e) {
    cat("⚠ saveGiotto failed, trying qs package...\n")
    
    # Try qs package (faster and handles large objects better)
    if (requireNamespace("qs", quietly = TRUE)) {
      tryCatch({
        output_file <- file.path(output_dir, paste0(sample_id, "_cosmx_loaded.qs"))
        qs::qsave(gobj, output_file, preset = "fast")
        cat("✓ Saved using qs package:", output_file, "\n")
        return(TRUE)
      }, error = function(e2) {
        cat("⚠ qs save also failed\n")
      })
    }
    
    # Last resort: save as RDS with reduced recursion
    tryCatch({
      output_file <- file.path(output_dir, paste0(sample_id, "_cosmx_loaded.rds"))
      # Increase stack limit temporarily
      old_limit <- Cstack_info()["size"]
      options(expressions = 500000)
      saveRDS(gobj, output_file, compress = "xz")
      cat("✓ Saved as RDS:", output_file, "\n")
      return(TRUE)
    }, error = function(e3) {
      cat("✗ All save methods failed\n")
      cat("  Error:", conditionMessage(e3), "\n")
      return(FALSE)
    })
  })
}

#' Main batch loading function
batch_load_cosmx <- function(sample_ids = NULL,
                             force_detect = FALSE,
                             save_intermediates = TRUE) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("                    BATCH LOADING CosMx SAMPLES                                \n")
  cat("================================================================================\n\n")
  
  # Get samples
  samples <- get_samples_smart(
    force_detect = force_detect,
    filter_ids = sample_ids
  )
  
  cat("=== Samples to Process ===\n")
  print(samples %>% select(sample_id, slide_num, directory_name))
  cat("\n")
  
  # Initialize results
  results <- tibble(
    sample_id = character(),
    slide_num = character(),
    status = character(),
    cells = integer(),
    genes = integer(),
    polygons = integer(),
    saved = logical(),
    time_seconds = numeric(),
    error_msg = character()
  )
  
  # Process each sample
  for (i in seq_len(nrow(samples))) {
    
    sample <- samples[i, ]
    
    cat("################################################################################\n")
    cat("## SAMPLE", i, "of", nrow(samples), ":", sample$sample_id, "\n")
    cat("##", "Slide", sample$slide_num, "|", sample$directory_name, "\n")
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
      
      # Save if requested
      saved <- FALSE
      if (save_intermediates) {
        saved <- save_giotto_safe(cosmx, sample$output_dir, sample$sample_id)
      }
      
      # Clean up
      rm(cosmx)
      gc()
      
      list(
        status = "SUCCESS",
        cells = n_cells,
        genes = n_genes,
        polygons = n_poly,
        saved = saved,
        error = NA_character_
      )
      
    }, error = function(e) {
      cat("\n✗ ERROR:", conditionMessage(e), "\n")
      list(
        status = "FAILED",
        cells = NA_integer_,
        genes = NA_integer_,
        polygons = NA_integer_,
        saved = FALSE,
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
        saved = result$saved,
        time_seconds = elapsed,
        error_msg = result$error
      )
    )
    
    if (result$status == "SUCCESS") {
      cat("\n✓✓✓ COMPLETE ✓✓✓\n")
      cat("    Cells:", format(result$cells, big.mark = ","), 
          "| Genes:", result$genes,
          "| Polygons:", format(result$polygons, big.mark = ","),
          "| Saved:", result$saved,
          "| Time:", round(elapsed, 1), "sec\n\n")
    } else {
      cat("\n✗✗✗ FAILED ✗✗✗\n\n")
    }
    
    # Force garbage collection
    gc(verbose = FALSE)
  }
  
  # Save results summary
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_file <- file.path("~/P_lab/CosMx_analysis/Output", 
                            paste0("batch_load_results_", timestamp, ".csv"))
  
  # Create Output directory if needed
  dir.create("~/P_lab/CosMx_analysis/Output", recursive = TRUE, showWarnings = FALSE)
  write_csv(results, results_file)
  
  # Print summary
  cat("\n")
  cat("================================================================================\n")
  cat("                         BATCH LOADING SUMMARY                                  \n")
  cat("================================================================================\n\n")
  
  print(results %>% select(sample_id, slide_num, status, cells, genes, saved, time_seconds))
  
  n_success <- sum(results$status == "SUCCESS")
  n_failed <- sum(results$status == "FAILED")
  n_saved <- sum(results$saved, na.rm = TRUE)
  
  cat("\nSuccess:", n_success, "/", nrow(results), "\n")
  cat("Saved:", n_saved, "/", n_success, "\n")
  
  if (n_failed > 0) {
    cat("Failed:", n_failed, "\n")
    cat("\nFailed samples:\n")
    failed <- results %>% filter(status == "FAILED")
    for (j in seq_len(nrow(failed))) {
      cat("  -", failed$sample_id[j], ":", failed$error_msg[j], "\n")
    }
  }
  
  cat("\nTotal time:", round(sum(results$time_seconds, na.rm = TRUE) / 60, 1), "minutes\n")
  if (n_success > 0) {
    cat("Average per sample:", round(mean(results$time_seconds[results$status == "SUCCESS"], na.rm = TRUE), 1), "seconds\n")
  }
  cat("\nResults saved to:", results_file, "\n\n")
  
  cat("================================================================================\n")
  cat("                           BATCH COMPLETE                                       \n")
  cat("================================================================================\n\n")
  
  return(results)
}

# Command-line interface
if (!interactive()) {
  
  library(optparse)
  
  option_list <- list(
    make_option(c("-s", "--samples"), type = "character", default = NULL,
                help = "Comma-separated sample IDs [default: all]"),
    make_option(c("-f", "--force-detect"), action = "store_true", default = FALSE,
                help = "Force re-detection [default: FALSE]"),
    make_option(c("--no-save"), action = "store_true", default = FALSE,
                help = "Don't save intermediate files [default: FALSE]")
  )
  
  parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
  args <- parse_args(parser)
  
  sample_ids <- if (!is.null(args$samples)) {
    strsplit(args$samples, ",")[[1]]
  } else {
    NULL
  }
  
  results <- batch_load_cosmx(
    sample_ids = sample_ids,
    force_detect = args$`force-detect`,
    save_intermediates = !args$`no-save`
  )
  
  if (any(results$status == "FAILED")) {
    quit(status = 1)
  } else {
    quit(status = 0)
  }
}
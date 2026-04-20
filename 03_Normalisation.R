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
if ((!exists("presentation_theme") || !exists("sample_plot_title") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
  source(pipeline_utils)
}

.run_known_giotto_warning_safe <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Not all expression matrices share the same cell_IDs", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

.sparse_matrix_median <- function(x) {
  if (!inherits(x, "sparseMatrix")) {
    return(stats::median(as.numeric(x), na.rm = TRUE))
  }
  
  x <- methods::as(x, "dgCMatrix")
  n_total <- nrow(x) * ncol(x)
  if (n_total == 0) {
    return(NA_real_)
  }
  
  zero_count <- n_total - Matrix::nnzero(x)
  nonzero_values <- sort(x@x)
  
  get_order_stat <- function(k) {
    if (k <= zero_count) {
      return(0)
    }
    nonzero_values[[k - zero_count]]
  }
  
  if (n_total %% 2 == 1) {
    return(get_order_stat((n_total + 1) / 2))
  }
  
  lower <- get_order_stat(n_total / 2)
  upper <- get_order_stat((n_total / 2) + 1)
  (lower + upper) / 2
}

.summarize_expression_matrix <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    n_total <- nrow(x) * ncol(x)
    x_c <- methods::as(x, "dgCMatrix")
    return(list(
      mean = as.numeric(sum(x_c) / n_total),
      median = .sparse_matrix_median(x_c),
      max = if (Matrix::nnzero(x_c) == 0) 0 else max(x_c@x)
    ))
  }
  
  list(
    mean = mean(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

.sanitize_expression_slots <- function(gobj) {
  expr_slots <- names(gobj@expression$cell$rna)
  if (!"raw" %in% expr_slots) {
    return(gobj)
  }
  
  reference_ids <- tryCatch(
    colnames(getExpression(gobj, values = "raw", output = "matrix")),
    error = function(e) NULL
  )
  
  if (is.null(reference_ids)) {
    return(gobj)
  }
  
  for (slot_name in setdiff(expr_slots, "raw")) {
    slot_matrix <- tryCatch(
      getExpression(gobj, values = slot_name, output = "matrix"),
      error = function(e) NULL
    )
    if (is.null(slot_matrix)) {
      next
    }
    
    slot_ids <- colnames(slot_matrix)
    if (is.null(slot_ids) || !identical(slot_ids, reference_ids)) {
      cat("⚠ Removing misaligned expression slot:", slot_name, "\n")
      gobj@expression$cell$rna[[slot_name]] <- NULL
    }
  }
  
  if ("scaled" %in% names(gobj@expression$cell$rna)) {
    cat("Removing unused scaled expression slot after normalization\n")
    gobj@expression$cell$rna$scaled <- NULL
  }
  
  gobj
}

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

  # Dynamic scalefactor: if NULL / 0 / negative, use the per-sample median
  # library size. Fixed 6000 does not fit mixed panel sizes (1K vs 6K).
  scalefactor_source <- "config"
  if (is.null(scalefactor) || !is.numeric(scalefactor) || scalefactor <= 0) {
    raw_mat <- getExpression(gobj, values = "raw", output = "matrix")
    lib_sizes <- Matrix::colSums(raw_mat)
    dyn_sf <- as.numeric(stats::median(lib_sizes, na.rm = TRUE))
    if (!is.finite(dyn_sf) || dyn_sf <= 0) dyn_sf <- 6000
    scalefactor <- round(dyn_sf)
    scalefactor_source <- "dynamic_median_lib_size"
    cat(sprintf(
      "Using dynamic scalefactor = %d (median library size across %d cells)\n",
      scalefactor, length(lib_sizes)
    ))
  }

  cat("Normalizing expression data...\n")
  cat("  Scale factor:", scalefactor, " (", scalefactor_source, ")\n", sep = "")
  cat("  Log transform:", log_transform, "\n\n")
  
  # Normalize
  gobj <- .run_known_giotto_warning_safe(
    normalizeGiotto(
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
  )
  
  gobj <- .sanitize_expression_slots(gobj)
  
  cat("\n✓ Normalization complete\n\n")
  
  # Check normalized data
  cat("=== Normalization Summary ===\n")
  cat("Expression slots available:\n")
  print(names(gobj@expression$cell$rna))
  
  # Compare raw vs normalized
  if ("normalized" %in% names(gobj@expression$cell$rna)) {
    
    raw_expr <- getExpression(gobj, values = "raw", output = "matrix")
    norm_expr <- getExpression(gobj, values = "normalized", output = "matrix")
    raw_summary <- .summarize_expression_matrix(raw_expr)
    norm_summary <- .summarize_expression_matrix(norm_expr)
    
    cat("\nRaw expression:\n")
    cat("  Mean:", round(raw_summary$mean, 2), "\n")
    cat("  Median:", round(raw_summary$median, 2), "\n")
    cat("  Max:", round(raw_summary$max, 2), "\n")
    
    cat("\nNormalized expression:\n")
    cat("  Mean:", round(norm_summary$mean, 2), "\n")
    cat("  Median:", round(norm_summary$median, 2), "\n")
    cat("  Max:", round(norm_summary$max, 2), "\n")
    
    # Create comparison plots
    comparison_data <- tibble(
      raw_mean = pmax(Matrix::rowMeans(raw_expr), 1e-8),
      norm_mean = pmax(Matrix::rowMeans(norm_expr), 1e-8)
    )
    
    p <- ggplot(comparison_data, aes(x = raw_mean, y = norm_mean)) +
      geom_point(alpha = 0.3, size = 0.5) +
      geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
      scale_x_log10() +
      scale_y_log10() +
      labs(
        title = sample_plot_title(sample_id, "Raw Vs Normalized Gene Expression"),
        subtitle = "Each point represents one gene. The dashed line marks equal mean expression before and after normalization.",
        x = log10_axis_label("Raw Mean Expression"),
        y = log10_axis_label("Normalized Mean Expression")
      ) +
      presentation_theme(base_size = 12)
    
    save_presentation_plot(
      plot = p,
      filename = file.path(results_folder, paste0(sample_id, "_raw_vs_normalized.png")),
      width = 8,
      height = 8,
      dpi = 300,
      bg = "white"
    )
    
    cat("\n✓ Comparison plots saved\n")
  }

  # Audit trail: record normalization parameters + summary statistics.
  tryCatch({
    raw_summary  <- if (exists("raw_summary",  inherits = FALSE)) raw_summary
                    else list(mean = NA_real_, median = NA_real_, max = NA_real_)
    norm_summary <- if (exists("norm_summary", inherits = FALSE)) norm_summary
                    else list(mean = NA_real_, median = NA_real_, max = NA_real_)
    norm_params <- tibble(
      parameter = c("scalefactor", "scalefactor_source",
                    "log_transform", "library_size_norm",
                    "norm_method", "log_offset",
                    "raw_mean", "raw_median", "raw_max",
                    "normalized_mean", "normalized_median", "normalized_max",
                    "timestamp"),
      value = c(as.character(scalefactor),
                scalefactor_source,
                as.character(log_transform),
                "TRUE",
                "standard",
                "1",
                as.character(round(raw_summary$mean,  3)),
                as.character(round(raw_summary$median, 3)),
                as.character(round(raw_summary$max,   3)),
                as.character(round(norm_summary$mean,   3)),
                as.character(round(norm_summary$median, 3)),
                as.character(round(norm_summary$max,    3)),
                format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
    )
    write_csv(
      norm_params,
      file.path(results_folder, paste0(sample_id, "_normalization_parameters.csv"))
    )
  }, error = function(e) {
    cat("  \u26A0 Could not write normalization_parameters.csv:",
        conditionMessage(e), "\n")
  })

  cat("\n✓ Normalization complete for", sample_id, "\n\n")

  return(gobj)
}

# Run if sourced directly
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) {
    script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
    script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
    bootstrap_script <- file.path(script_dir, "Helper_Scripts", "Script_Bootstrap.R")
    if (file.exists(bootstrap_script)) {
      source(bootstrap_script, local = .GlobalEnv)
      bootstrap_pipeline_environment(script_dir, load_pipeline_utils = FALSE, verbose = FALSE)
    }

    sample_id <- args[1]
    input_path <- args[2]
    output_dir <- args[3]
    
    gobj <- normalize_expression(
      gobj = input_path,
      sample_id = sample_id,
      output_dir = output_dir
    )
    
    saveGiotto(gobj, dir = output_dir, foldername = "Giotto_Object_Normalized", overwrite = TRUE)
  } else {
    stop("Usage: Rscript 03_Normalisation.R <sample_id> <input_path> <output_dir>")
  }
}

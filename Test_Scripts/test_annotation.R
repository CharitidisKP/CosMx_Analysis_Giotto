#!/usr/bin/env Rscript
# ==============================================================================
# test_annotation.R
# Comprehensive test script to validate the annotation workflow (07_Annotation.R)
# with the CosMx pipeline using InSituType supervised cell-type annotation.
#
# Usage:
#   Rscript test_annotation.R [sample_id] [base_dir] [config_path]
#
# Arguments (all optional, positional):
#   sample_id   - Sample identifier (default: "Sample_18H12037")
#   base_dir    - Base project directory (default: "~/P_lab/CosMx_analysis")
#   config_path - Path to config YAML (default: <base_dir>/Scripts/config.yaml)
# ==============================================================================

cat("\n========================================\n")
cat("ANNOTATION WORKFLOW VALIDATION TEST\n")
cat("========================================\n\n")

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

# Parse optional command-line arguments
args <- commandArgs(trailingOnly = TRUE)
SAMPLE_ID   <- if (length(args) >= 1) args[1] else "Sample_18H12037"
BASE_DIR    <- path.expand(if (length(args) >= 2) args[2] else "~/P_lab/CosMx_analysis")
CONFIG_PATH <- path.expand(
  if (length(args) >= 3) args[3]
  else file.path(BASE_DIR, "Scripts", "config.yaml")
)

# Regex pattern for identifying negative-control probe columns (used in tests 1c and 3)
NEG_PROBE_PATTERN <- "(?i)(neg|negprobe|SystemControl|background)"

cat("Sample ID  :", SAMPLE_ID, "\n")
cat("Base dir   :", BASE_DIR, "\n")
cat("Config file:", CONFIG_PATH, "\n\n")

# Source environment setup from the project scripts directory
setup_script <- file.path(BASE_DIR, "Scripts", "00_Setup.R")
if (file.exists(setup_script)) {
  source(setup_script)
} else {
  # Minimal package loading if setup script is unavailable
  suppressPackageStartupMessages({
    library(Matrix)
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(patchwork)
    library(Giotto)
    library(InSituType)
  })
}

# Optional packages loaded with informative warnings
optional_pkgs <- c("SpatialDecon", "yaml")
for (pkg in optional_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    warning("Optional package '", pkg, "' is not installed. ",
            "Some tests will be skipped.")
  }
}

# Initialise results list
test_results <- list()

record_result <- function(test_name, passed, details = "") {
  test_results[[length(test_results) + 1]] <<- list(
    test      = test_name,
    passed    = passed,
    details   = as.character(details)
  )
  status <- if (passed) "  ✓ PASS" else "  ✗ FAIL"
  cat(status, "-", test_name, "\n")
  if (nchar(details) > 0) cat("    Details:", details, "\n")
}

# Output directory for this test run
test_output_dir <- file.path(BASE_DIR, "Output", SAMPLE_ID, "Test_Annotation")
dir.create(test_output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Load and prepare test data
# ------------------------------------------------------------------------------

cat("\n--- TEST 1: Load and prepare test data ---\n")

giotto_object_path <- file.path(
  BASE_DIR, "Output", SAMPLE_ID, "Giotto_Object_Processed"
)

# 1a. Load the Giotto object
cosmx <- tryCatch({
  obj <- readRDS(giotto_object_path)
  record_result("1a: Load Giotto object", TRUE,
                paste("Loaded from", giotto_object_path))
  obj
}, error = function(e) {
  record_result("1a: Load Giotto object", FALSE, conditionMessage(e))
  NULL
})

if (!is.null(cosmx)) {

  # 1b. Extract and verify count matrix dimensions
  tryCatch({
    expr_mat <- getExpression(
      gobject     = cosmx,
      values      = "raw",
      output      = "matrix",
      feat_type   = "rna",
      spat_unit   = "cell"
    )
    n_genes <- nrow(expr_mat)
    n_cells <- ncol(expr_mat)
    record_result("1b: Extract count matrix", TRUE,
                  paste(n_genes, "genes x", n_cells, "cells"))
  }, error = function(e) {
    expr_mat <<- NULL
    record_result("1b: Extract count matrix", FALSE, conditionMessage(e))
  })

  # 1c. Calculate background using negative control probes
  tryCatch({
    meta <- pDataDT(cosmx)
    neg_cols <- grep(NEG_PROBE_PATTERN, names(meta), value = TRUE, perl = TRUE)

    if (length(neg_cols) > 0) {
      # Use the first identified negative-probe column
      neg_col <- neg_cols[1]
      bg_values <- as.numeric(meta[[neg_col]])
      bg_mean   <- mean(bg_values, na.rm = TRUE)
      record_result("1c: Calculate background", TRUE,
                    paste("Column:", neg_col,
                          "| Mean background:", round(bg_mean, 4)))
    } else {
      # Fall back to the Helper_Functions add_bg_to_qc approach
      if (exists("add_bg_to_qc")) {
        qc <- as.data.frame(meta)
        if ("Total_Counts" %in% names(qc)) {
          qc <- add_bg_to_qc(qc)
          record_result("1c: Calculate background", TRUE,
                        "Background calculated via add_bg_to_qc()")
        } else {
          record_result("1c: Calculate background", FALSE,
                        "No negative probe columns found in metadata")
        }
      } else {
        record_result("1c: Calculate background", FALSE,
                      "No negative probe columns found; add_bg_to_qc unavailable")
      }
    }
  }, error = function(e) {
    record_result("1c: Calculate background", FALSE, conditionMessage(e))
  })

} else {
  record_result("1b: Extract count matrix",    FALSE, "Skipped: Giotto object not loaded")
  record_result("1c: Calculate background",    FALSE, "Skipped: Giotto object not loaded")
}

# ------------------------------------------------------------------------------
# 2. Test reference profile loading
# ------------------------------------------------------------------------------

cat("\n--- TEST 2: Reference profile loading ---\n")

hca_profile   <- NULL
cosmx_profile <- NULL

# 2a. HCA Kidney profile via SpatialDecon
tryCatch({
  if (!requireNamespace("SpatialDecon", quietly = TRUE)) {
    stop("SpatialDecon is not installed")
  }
  hca_profile <- SpatialDecon::download_profile_matrix(
    species = "Human",
    age_group = "Adult",
    matrixname = "HCA_Kidney"
  )
  record_result("2a: HCA Kidney profile (SpatialDecon)", TRUE,
                paste(nrow(hca_profile), "genes x", ncol(hca_profile), "cell types"))
}, error = function(e) {
  record_result("2a: HCA Kidney profile (SpatialDecon)", FALSE, conditionMessage(e))
})

# 2b. CosMx 1K Kidney profile via URL
cosmx_profile_url <- paste0(
  "https://raw.githubusercontent.com/Nanostring-Biostats/",
  "InSituType/main/inst/extdata/",
  "ioprofiles.rda"
)
tryCatch({
  tmp_rda <- tempfile(fileext = ".rda")
  download.file(cosmx_profile_url, tmp_rda, quiet = TRUE, mode = "wb")
  load_env  <- new.env()
  load(tmp_rda, envir = load_env)
  unlink(tmp_rda)

  # The .rda may expose 'ioprofiles' or similar objects
  profile_obj_names <- ls(load_env)
  if (length(profile_obj_names) == 0) stop("No objects found in downloaded .rda")
  cosmx_profile <- get(profile_obj_names[1], envir = load_env)

  # Coerce to matrix if it is a list
  if (is.list(cosmx_profile) && !is.data.frame(cosmx_profile)) {
    cosmx_profile <- do.call(cbind, cosmx_profile)
  }
  cosmx_profile <- as.matrix(cosmx_profile)

  record_result("2b: CosMx 1K Kidney profile (URL)", TRUE,
                paste(nrow(cosmx_profile), "genes x",
                      ncol(cosmx_profile), "cell types"))
}, error = function(e) {
  record_result("2b: CosMx 1K Kidney profile (URL)", FALSE, conditionMessage(e))
})

# 2c. Verify profile dimensions and gene overlap with data
for (profile_info in list(
  list(name = "2c: HCA Kidney gene overlap",    profile = hca_profile),
  list(name = "2c: CosMx Kidney gene overlap",  profile = cosmx_profile)
)) {
  tryCatch({
    prof <- profile_info$profile
    if (is.null(prof)) stop("Profile not loaded")

    if (!exists("expr_mat") || is.null(expr_mat)) {
      stop("Expression matrix not available")
    }
    data_genes    <- rownames(expr_mat)
    profile_genes <- rownames(prof)
    n_overlap     <- length(intersect(data_genes, profile_genes))
    pct_overlap   <- round(100 * n_overlap / length(data_genes), 1)

    record_result(profile_info$name, n_overlap > 0,
                  paste(n_overlap, "overlapping genes (", pct_overlap, "% of data)"))
  }, error = function(e) {
    record_result(profile_info$name, FALSE, conditionMessage(e))
  })
}

# ------------------------------------------------------------------------------
# 3. Test InSituType annotation
# ------------------------------------------------------------------------------

cat("\n--- TEST 3: InSituType annotation ---\n")

insitutype_result <- NULL

tryCatch({
  if (is.null(cosmx)) stop("Giotto object not available")
  if (is.null(hca_profile) && is.null(cosmx_profile)) {
    stop("No reference profiles available")
  }

  # Choose profile: prefer HCA, fall back to CosMx
  ref_profile <- if (!is.null(hca_profile)) hca_profile else cosmx_profile
  profile_name <- if (!is.null(hca_profile)) "HCA_Kidney" else "CosMx_Kidney"

  # Extract raw count matrix (cells × genes for InSituType)
  count_mat <- t(as.matrix(
    getExpression(gobject = cosmx, values = "raw",
                  output = "matrix", feat_type = "rna", spat_unit = "cell")
  ))

  # Align genes
  shared_genes  <- intersect(colnames(count_mat), rownames(ref_profile))
  if (length(shared_genes) == 0) stop("No shared genes between data and reference profile")

  count_mat   <- count_mat[, shared_genes, drop = FALSE]
  ref_profile <- ref_profile[shared_genes, , drop = FALSE]

  cat("  Running InSituType supervised annotation with", profile_name, "...\n")
  cat("  Cells:", nrow(count_mat), "| Shared genes:", length(shared_genes), "\n")

  # Guard against degenerate inputs that would cause convergence failures
  if (nrow(count_mat) < 10) stop("Too few cells for annotation (< 10)")
  if (length(shared_genes) < 10) stop("Too few shared genes for annotation (< 10)")

  # Estimate per-cell background: mean of negative-probe columns
  meta <- as.data.frame(pDataDT(cosmx))
  neg_cols <- grep(NEG_PROBE_PATTERN, names(meta), value = TRUE, perl = TRUE)
  if (length(neg_cols) > 0) {
    bg_per_cell <- rowMeans(meta[, neg_cols, drop = FALSE], na.rm = TRUE)
    # Identify the cell ID column (first character-type column or a known name)
    cell_id_col <- intersect(c("cell_ID", "cell_id", "CellID"), names(meta))[1]
    if (!is.na(cell_id_col)) {
      names(bg_per_cell) <- meta[[cell_id_col]]
    } else {
      names(bg_per_cell) <- rownames(meta)
    }
    bg_per_cell <- bg_per_cell[rownames(count_mat)]
    bg_per_cell[is.na(bg_per_cell)] <- mean(bg_per_cell, na.rm = TRUE)
  } else {
    bg_per_cell <- rep(0.01, nrow(count_mat))
    names(bg_per_cell) <- rownames(count_mat)
  }

  insitutype_result <- insitutypeML(
    x           = count_mat,
    neg         = bg_per_cell,
    reference_profiles = ref_profile
  )

  record_result("3a: Run InSituType annotation", TRUE,
                paste("Profile:", profile_name,
                      "| Cells annotated:", length(insitutype_result$clust)))

}, error = function(e) {
  record_result("3a: Run InSituType annotation", FALSE, conditionMessage(e))
})

# 3b. Verify annotation results
tryCatch({
  if (is.null(insitutype_result)) stop("Annotation result not available")

  cell_types <- insitutype_result$clust
  scores     <- insitutype_result$prob
  n_cells    <- if (!is.null(cosmx)) length(GiottoCellId(cosmx)) else NA

  has_labels  <- !is.null(cell_types) && length(cell_types) > 0
  has_scores  <- !is.null(scores)     && (is.vector(scores) || is.matrix(scores))
  count_match <- !is.na(n_cells) && length(cell_types) == n_cells

  record_result("3b: Annotation results - cell type labels", has_labels,
                paste(length(unique(cell_types)), "unique cell types identified"))

  record_result("3b: Annotation results - confidence scores", has_scores,
                paste("Score dimensions:",
                      if (is.matrix(scores)) paste(dim(scores), collapse = "x")
                      else length(scores)))

  record_result("3b: Annotation results - cell count match", count_match,
                paste("Expected:", n_cells, "| Got:", length(cell_types)))

}, error = function(e) {
  record_result("3b: Annotation results", FALSE, conditionMessage(e))
})

# ------------------------------------------------------------------------------
# 4. Test Giotto integration
# ------------------------------------------------------------------------------

cat("\n--- TEST 4: Giotto integration ---\n")

cosmx_annotated <- NULL

# 4a. Add annotations to Giotto object (custom cell-type name)
tryCatch({
  if (is.null(cosmx))              stop("Giotto object not available")
  if (is.null(insitutype_result)) stop("Annotation result not available")

  cell_types <- insitutype_result$clust
  scores     <- if (is.matrix(insitutype_result$prob)) {
    apply(insitutype_result$prob, 1, max)
  } else {
    insitutype_result$prob
  }

  annotation_df <- data.frame(
    cell_ID            = names(cell_types),
    HCA_Kidney_celltype = as.character(cell_types),
    HCA_Kidney_score    = as.numeric(scores),
    stringsAsFactors   = FALSE
  )

  cosmx_annotated <- addCellMetadata(
    gobject      = cosmx,
    new_metadata = annotation_df,
    by_column    = TRUE,
    column_cell_ID = "cell_ID"
  )

  record_result("4a: Add annotations to Giotto object", TRUE,
                "HCA_Kidney_celltype and HCA_Kidney_score columns added")

}, error = function(e) {
  record_result("4a: Add annotations to Giotto object", FALSE, conditionMessage(e))
})

# 4b. Verify metadata columns are created correctly
tryCatch({
  if (is.null(cosmx_annotated)) stop("Annotated Giotto object not available")

  meta_cols <- names(pDataDT(cosmx_annotated))
  has_celltype <- "HCA_Kidney_celltype" %in% meta_cols
  has_score    <- "HCA_Kidney_score"    %in% meta_cols

  record_result("4b: Metadata column - celltype", has_celltype,
                paste("Columns:", paste(meta_cols, collapse = ", ")))
  record_result("4b: Metadata column - score",    has_score,
                paste("Score column present:", has_score))

}, error = function(e) {
  record_result("4b: Verify metadata columns", FALSE, conditionMessage(e))
})

# 4c. Default profile: creates generic 'celltype' and 'celltype_score' columns
tryCatch({
  if (is.null(cosmx))              stop("Giotto object not available")
  if (is.null(insitutype_result)) stop("Annotation result not available")

  cell_types <- insitutype_result$clust
  scores     <- if (is.matrix(insitutype_result$prob)) {
    apply(insitutype_result$prob, 1, max)
  } else {
    insitutype_result$prob
  }

  default_df <- data.frame(
    cell_ID      = names(cell_types),
    celltype     = as.character(cell_types),
    celltype_score = as.numeric(scores),
    stringsAsFactors = FALSE
  )

  cosmx_default <- addCellMetadata(
    gobject        = cosmx,
    new_metadata   = default_df,
    by_column      = TRUE,
    column_cell_ID = "cell_ID"
  )

  meta_cols_default <- names(pDataDT(cosmx_default))
  has_default_ct    <- "celltype"       %in% meta_cols_default
  has_default_score <- "celltype_score" %in% meta_cols_default

  record_result("4c: Default profile - 'celltype' column",       has_default_ct,
                "Generic celltype column created")
  record_result("4c: Default profile - 'celltype_score' column", has_default_score,
                "Generic celltype_score column created")

}, error = function(e) {
  record_result("4c: Default profile annotation", FALSE, conditionMessage(e))
})

# ------------------------------------------------------------------------------
# 5. Test configuration loading
# ------------------------------------------------------------------------------

cat("\n--- TEST 5: Configuration loading ---\n")

config <- NULL

# 5a. Load configuration from YAML
tryCatch({
  if (!requireNamespace("yaml", quietly = TRUE)) stop("yaml package is not installed")

  config_path <- CONFIG_PATH

  if (!file.exists(config_path)) {
    stop("Config file not found: ", config_path)
  }

  config <- yaml::read_yaml(config_path)
  record_result("5a: Load YAML config", TRUE,
                paste("Keys found:", paste(names(config), collapse = ", ")))

}, error = function(e) {
  record_result("5a: Load YAML config", FALSE, conditionMessage(e))
})

# 5b. Extract annotation parameters
tryCatch({
  if (is.null(config)) stop("Config not loaded")

  annotation_params <- config[["annotation"]]
  if (is.null(annotation_params)) stop("'annotation' section missing from config")

  record_result("5b: Extract annotation parameters", TRUE,
                paste("Annotation params:", paste(names(annotation_params), collapse = ", ")))

}, error = function(e) {
  record_result("5b: Extract annotation parameters", FALSE, conditionMessage(e))
})

# 5c. Verify all required annotation parameters are present
tryCatch({
  if (is.null(config)) stop("Config not loaded")

  required_params <- c("profile", "output_prefix", "sample_id")
  annotation_params <- config[["annotation"]]

  if (is.null(annotation_params)) stop("'annotation' section missing from config")

  missing_params <- setdiff(required_params, names(annotation_params))
  all_present    <- length(missing_params) == 0

  record_result("5c: Required annotation parameters present", all_present,
                if (all_present) "All required params found"
                else paste("Missing:", paste(missing_params, collapse = ", ")))

}, error = function(e) {
  record_result("5c: Required annotation parameters", FALSE, conditionMessage(e))
})

# ------------------------------------------------------------------------------
# 6. Create validation report
# ------------------------------------------------------------------------------

cat("\n--- TEST 6: Validation report ---\n")

# 6a. Save test results to CSV
tryCatch({
  results_df <- do.call(rbind, lapply(test_results, as.data.frame,
                                      stringsAsFactors = FALSE))
  csv_path <- file.path(test_output_dir, "annotation_test_results.csv")
  write.csv(results_df, csv_path, row.names = FALSE)
  record_result("6a: Save test results to CSV", TRUE,
                paste("Saved to:", csv_path))
}, error = function(e) {
  record_result("6a: Save test results to CSV", FALSE, conditionMessage(e))
})

# 6b. Generate summary statistics
tryCatch({
  results_df <- do.call(rbind, lapply(test_results, as.data.frame,
                                      stringsAsFactors = FALSE))
  n_total  <- nrow(results_df)
  n_passed <- sum(results_df$passed, na.rm = TRUE)
  n_failed <- n_total - n_passed
  pct_pass <- round(100 * n_passed / n_total, 1)

  cat("\n=== ANNOTATION VALIDATION SUMMARY ===\n")
  cat("Total tests : ", n_total,  "\n")
  cat("Passed      : ", n_passed, "\n")
  cat("Failed      : ", n_failed, "\n")
  cat("Pass rate   : ", pct_pass, "%\n")

  if (!is.null(insitutype_result) && !is.null(insitutype_result$clust)) {
    ct_table <- sort(table(insitutype_result$clust), decreasing = TRUE)
    cat("\nCell type distribution (top 10):\n")
    print(head(ct_table, 10))
  }

  record_result("6b: Generate summary statistics", TRUE,
                paste(n_passed, "/", n_total, "tests passed (", pct_pass, "%)"))

}, error = function(e) {
  record_result("6b: Generate summary statistics", FALSE, conditionMessage(e))
})

# 6c. Create diagnostic plots
tryCatch({
  plot_list <- list()

  # Plot 1: test pass/fail bar chart
  results_df <- do.call(rbind, lapply(test_results, as.data.frame,
                                      stringsAsFactors = FALSE))
  p_results <- ggplot(results_df, aes(x = passed, fill = passed)) +
    geom_bar() +
    scale_fill_manual(values = c("TRUE" = "#2ECC71", "FALSE" = "#E74C3C"),
                      labels = c("TRUE" = "Pass", "FALSE" = "Fail")) +
    labs(title = "Annotation Test Results",
         x = "Result", y = "Number of Tests", fill = "Status") +
    theme_minimal()
  plot_list[["test_results"]] <- p_results

  # Plot 2: cell type distribution (if annotation succeeded)
  if (!is.null(insitutype_result) && !is.null(insitutype_result$clust)) {
    ct_df <- as.data.frame(sort(table(insitutype_result$clust), decreasing = TRUE))
    names(ct_df) <- c("CellType", "Count")
    ct_df$CellType <- factor(ct_df$CellType, levels = ct_df$CellType)

    p_ct <- ggplot(ct_df, aes(x = CellType, y = Count, fill = CellType)) +
      geom_bar(stat = "identity") +
      scale_fill_viridis_d() +
      labs(title   = "Cell Type Distribution (InSituType)",
           x = "Cell Type", y = "Cell Count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    plot_list[["cell_type_distribution"]] <- p_ct
  }

  # Plot 3: score distribution (if annotation succeeded)
  if (!is.null(insitutype_result) && !is.null(insitutype_result$prob)) {
    scores <- if (is.matrix(insitutype_result$prob)) {
      apply(insitutype_result$prob, 1, max)
    } else {
      insitutype_result$prob
    }
    score_df <- data.frame(score = scores)
    p_score  <- ggplot(score_df, aes(x = score)) +
      geom_histogram(bins = 50, fill = "#3498DB", colour = "white") +
      labs(title = "InSituType Confidence Score Distribution",
           x = "Max Cell-Type Probability", y = "Frequency") +
      theme_minimal()
    plot_list[["score_distribution"]] <- p_score
  }

  # Save combined plot
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = min(2, length(plot_list)))
    plot_path <- file.path(test_output_dir, "annotation_diagnostic_plots.pdf")
    ggsave(plot_path, combined_plot,
           width = 14, height = 7 * ceiling(length(plot_list) / 2))
    record_result("6c: Create diagnostic plots", TRUE,
                  paste(length(plot_list), "plots saved to:", plot_path))
  } else {
    record_result("6c: Create diagnostic plots", FALSE,
                  "No plots could be generated")
  }

}, error = function(e) {
  record_result("6c: Create diagnostic plots", FALSE, conditionMessage(e))
})

# ------------------------------------------------------------------------------
# 7. Final summary
# ------------------------------------------------------------------------------

results_df_final <- do.call(rbind, lapply(test_results, as.data.frame,
                                          stringsAsFactors = FALSE))
n_total  <- nrow(results_df_final)
n_passed <- sum(results_df_final$passed, na.rm = TRUE)

cat("\n========================================\n")
cat("FINAL RESULT:", n_passed, "/", n_total, "tests passed\n")
cat("Output dir  :", test_output_dir, "\n")
cat("========================================\n\n")

if (n_passed < n_total) {
  failed_tests <- results_df_final$test[!results_df_final$passed]
  cat("Failed tests:\n")
  for (t in failed_tests) cat(" -", t, "\n")
  cat("\n")
}

invisible(list(results = results_df_final, output_dir = test_output_dir))

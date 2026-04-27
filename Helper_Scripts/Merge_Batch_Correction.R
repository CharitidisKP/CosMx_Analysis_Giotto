#!/usr/bin/env Rscript
# ==============================================================================
# Merge_Batch_Correction.R
# Merge Giotto objects and correct slide-level batch effects with Harmony
# ==============================================================================

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

pipeline_utils <- file.path(current_script_dir(), "Pipeline_Utils.R")
if ((!exists("save_giotto_checkpoint") || !exists("presentation_theme") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
  source(pipeline_utils)
}

.giotto_pdata_dt <- function(gobj) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("pDataDT", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("pDataDT", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("pDataDT", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("pDataDT", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("pDataDT", mode = "function")
  }
  
  accessor(gobj)
}

.giotto_load <- function(path) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("loadGiotto", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("loadGiotto", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("loadGiotto", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("loadGiotto", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("loadGiotto", mode = "function")
  }
  
  accessor(path)
}

.giotto_get_dim_reduction <- function(gobj, reduction_method, name = NULL, output = "data.table") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("get_dimReduction", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("get_dimReduction", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("get_dimReduction", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("get_dimReduction", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("get_dimReduction", mode = "function")
  }
  
  accessor(
    gobject = gobj,
    reduction = "cells",
    reduction_method = reduction_method,
    name = name,
    output = output
  )
}

.giotto_join_objects <- function(gobject_list, gobject_names, join_method = "shift", x_padding = 1000) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("joinGiottoObjects", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("joinGiottoObjects", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("joinGiottoObjects", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("joinGiottoObjects", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("joinGiottoObjects", mode = "function")
  }
  
  accessor(
    gobject_list = gobject_list,
    gobject_names = gobject_names,
    join_method = join_method,
    x_padding = x_padding
  )
}

.giotto_run_pca <- function(gobj, feats_to_use = NULL, n_pcs = 30) {
  .run_known_giotto_warning_safe(
    runPCA(
      gobject   = gobj,
      feats_to_use = feats_to_use,
      scale_unit = TRUE,
      center     = TRUE,
      ncp        = n_pcs
    )
  )
}

.giotto_run_harmony <- function(gobj, vars_use, dimensions_to_use, name = "harmony") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("runGiottoHarmony", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("runGiottoHarmony", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("runGiottoHarmony", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("runGiottoHarmony", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("runGiottoHarmony", mode = "function")
  }
  
  accessor(
    gobject = gobj,
    vars_use = vars_use,
    dim_reduction_to_use = "pca",
    dimensions_to_use = dimensions_to_use,
    name = name
  )
}

.giotto_run_umap <- function(gobj,
                             dim_reduction_to_use,
                             dim_reduction_name,
                             dimensions_to_use,
                             n_neighbors,
                             min_dist,
                             name = "umap") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("runUMAP", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("runUMAP", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("runUMAP", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("runUMAP", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("runUMAP", mode = "function")
  }
  
  accessor(
    gobject = gobj,
    dim_reduction_to_use = dim_reduction_to_use,
    dim_reduction_name = dim_reduction_name,
    dimensions_to_use = dimensions_to_use,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    name = name
  )
}

.giotto_run_tsne <- function(gobj,
                             dim_reduction_to_use,
                             dim_reduction_name,
                             dimensions_to_use,
                             perplexity = 30,
                             name = "tsne") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("runtSNE", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("runtSNE", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("runtSNE", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("runtSNE", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("runtSNE", mode = "function")
  }
  
  accessor(
    gobject = gobj,
    dim_reduction_to_use = dim_reduction_to_use,
    dim_reduction_name = dim_reduction_name,
    dimensions_to_use = dimensions_to_use,
    perplexity = perplexity,
    name = name
  )
}

.giotto_create_nearest_network <- function(gobj,
                                           dim_reduction_to_use,
                                           dim_reduction_name,
                                           dimensions_to_use,
                                           k,
                                           name = "NN.harmony") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("createNearestNetwork", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("createNearestNetwork", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("createNearestNetwork", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("createNearestNetwork", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("createNearestNetwork", mode = "function")
  }
  
  accessor(
    gobject = gobj,
    dim_reduction_to_use = dim_reduction_to_use,
    dim_reduction_name = dim_reduction_name,
    dimensions_to_use = dimensions_to_use,
    k = k,
    name = name
  )
}

.giotto_do_leiden_cluster <- function(gobj,
                                      network_name,
                                      resolution = 0.3,
                                      n_iterations = 1000,
                                      name = "leiden_clust") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("doLeidenCluster", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("doLeidenCluster", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("doLeidenCluster", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("doLeidenCluster", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("doLeidenCluster", mode = "function")
  }
  
  accessor(
    gobject = gobj,
    network_name = network_name,
    resolution = resolution,
    n_iterations = n_iterations,
    name = name
  )
}

.sub_ckpt_exists <- function(dir) {
  dir.exists(file.path(dir, "giotto")) ||
    file.exists(file.path(dir, "object.qs")) ||
    file.exists(file.path(dir, "object.rds"))
}

parse_dimension_vector <- function(x) {
  if (is.null(x)) {
    return(1:30)
  }
  if (is.numeric(x)) {
    return(as.integer(x))
  }
  eval(parse(text = x))
}

embedding_to_tibble <- function(gobj, reduction_method, name = NULL, color_columns = character()) {
  coords <- .giotto_get_dim_reduction(
    gobj = gobj,
    reduction_method = reduction_method,
    name = name,
    output = "data.table"
  )
  
  coord_df <- tibble::as_tibble(coords)
  if (!"cell_ID" %in% names(coord_df)) {
    names(coord_df)[1] <- "cell_ID"
  }
  dim_cols <- setdiff(names(coord_df), "cell_ID")
  dim_cols <- dim_cols[seq_len(min(2, length(dim_cols)))]
  coord_df <- dplyr::select(coord_df, cell_ID, dplyr::all_of(dim_cols))
  names(coord_df)[2:ncol(coord_df)] <- c("dim_1", "dim_2")[seq_len(ncol(coord_df) - 1)]
  
  metadata <- tibble::as_tibble(.giotto_pdata_dt(gobj))
  metadata <- dplyr::select(metadata, dplyr::any_of(c("cell_ID", color_columns)))

  coord_df$cell_ID <- as.character(coord_df$cell_ID)
  metadata$cell_ID <- as.character(metadata$cell_ID)

  dplyr::left_join(coord_df, metadata, by = "cell_ID")
}

save_embedding_plot <- function(df, color_column, title, output_file) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = dim_1, y = dim_2, color = rlang::.data[[color_column]])) +
    ggplot2::geom_point(size = 0.4, alpha = 0.8) +
    ggplot2::labs(
      title = title,
      x = "Embedding Dimension 1",
      y = "Embedding Dimension 2",
      color = pretty_plot_label(color_column)
    ) +
    presentation_theme(base_size = 12, legend_position = "right") +
    ggplot2::theme(
      legend.key.height = grid::unit(0.45, "cm")
    )
  
  save_presentation_plot(p, output_file, width = 10, height = 8)
}

merge_giotto_samples <- function(gobject_list,
                                 sample_table,
                                 output_dir,
                                 join_method = "shift",
                                 x_padding = 1000,
                                 renormalize = TRUE,
                                 scalefactor = NULL,
                                 log_transform = TRUE,
                                 save_object = FALSE) {
  if (length(gobject_list) == 0) {
    stop("No Giotto objects were provided for merging.")
  }

  if (length(gobject_list) == 1) {
    merged_gobj <- gobject_list[[1]]
  } else {
    merged_gobj <- .giotto_join_objects(
      gobject_list = gobject_list,
      gobject_names = sample_table$sample_id,
      join_method = join_method,
      x_padding = x_padding
    )
  }

  # Re-normalize on the merged raw counts so per-sample scale-factor mismatches
  # don't leak into HVG selection / Harmony / any downstream consumer of the
  # normalized matrix. Pseudobulk DE engines use raw counts and are unaffected.
  if (isTRUE(renormalize)) {
    merged_gobj <- .renormalize_merged_object(
      merged_gobj,
      scalefactor   = scalefactor,
      log_transform = log_transform
    )
  } else {
    cat("Skipping merged-object renormalization (renormalize = FALSE).\n",
        "  HVG/PCA/Harmony will run on the per-sample normalized matrix as-joined.\n",
        sep = "")
  }

  results_dir <- ensure_dir(file.path(output_dir, "10_Merged"))
  readr::write_csv(sample_table, file.path(results_dir, "merged_sample_manifest.csv"))
  readr::write_csv(
    tibble::as_tibble(.giotto_pdata_dt(merged_gobj)),
    file.path(results_dir, "merged_cell_metadata.csv")
  )

  if (save_object) {
    save_giotto_checkpoint(
      gobj = merged_gobj,
      checkpoint_dir = file.path(output_dir, "Giotto_Object_Merged"),
      metadata = list(stage = "merge",
                      renormalize = renormalize,
                      scalefactor = scalefactor %||% "dynamic_median_lib_size")
    )
  }

  merged_gobj
}

# Re-runs library-size + log normalization on the merged object using a
# global scale factor. Mirrors 03_Normalisation.R's call so the merged
# normalized values follow the same definition as the per-sample ones.
.renormalize_merged_object <- function(merged_gobj,
                                       scalefactor = NULL,
                                       log_transform = TRUE) {
  raw_mat <- tryCatch({
    if (requireNamespace("Giotto", quietly = TRUE) &&
        exists("getExpression", envir = asNamespace("Giotto"),
               inherits = FALSE)) {
      Giotto::getExpression(merged_gobj, values = "raw", output = "matrix")
    } else if (requireNamespace("GiottoClass", quietly = TRUE)) {
      GiottoClass::getExpression(merged_gobj, values = "raw",
                                 output = "matrix")
    } else NULL
  }, error = function(e) NULL)

  if (is.null(raw_mat)) {
    cat("⚠ Merged object has no raw expression slot — skipping ",
        "re-normalization.\n", sep = "")
    return(merged_gobj)
  }

  lib_sizes <- Matrix::colSums(raw_mat)
  median_lib <- stats::median(lib_sizes[lib_sizes > 0])
  if (!is.finite(median_lib) || median_lib <= 0) median_lib <- 6000

  source_tag <- "config"
  if (is.null(scalefactor) || !is.numeric(scalefactor) || scalefactor <= 0) {
    scalefactor <- round(median_lib)
    source_tag <- "dynamic_median_lib_size"
  }

  cat("Re-normalizing merged object on raw counts...\n")
  cat("  Scale factor   : ", scalefactor, " (", source_tag, ")\n", sep = "")
  cat("  Median lib size: ", round(median_lib), " across ",
      length(lib_sizes), " cells\n", sep = "")
  cat("  Log transform  : ", log_transform, "\n", sep = "")

  merged_gobj <- .run_known_giotto_warning_safe(
    normalizeGiotto(
      gobject           = merged_gobj,
      spat_unit         = "cell",
      feat_type         = "rna",
      scalefactor       = scalefactor,
      verbose           = FALSE,
      norm_methods      = "standard",
      library_size_norm = TRUE,
      log_norm          = log_transform,
      log_offset        = 1,
      scale_feats       = FALSE,
      scale_cells       = FALSE
    )
  )

  cat("✓ Merged-object re-normalization complete\n\n")
  merged_gobj
}

batch_correct_merged_object <- function(gobj,
                                        sample_id = "merged",
                                        output_dir,
                                        batch_column = "slide_id",
                                        annotation_columns = NULL,
                                        n_hvgs = 500,
                                        n_pcs = 30,
                                        dimensions_to_use = 1:30,
                                        umap_n_neighbors = 30,
                                        umap_min_dist = 0.3,
                                        k_nn = 15,
                                        resolution = 0.3,
                                        create_plots = TRUE,
                                        save_object = FALSE,
                                        force_recompute = FALSE) {
  cat("\n========================================\n")
  cat("MERGED: Batch Correction\n")
  cat("Run:", sample_id, "\n")
  cat("========================================\n\n")

  if (is.character(gobj)) {
    manifest_path <- file.path(gobj, "manifest.json")
    if (file.exists(manifest_path)) {
      gobj <- load_giotto_checkpoint(gobj)
    } else {
      gobj <- .giotto_load(gobj)
    }
  }

  results_dir       <- ensure_dir(file.path(output_dir, "10_Merged", "Batch_Correction"))
  dimensions_to_use <- parse_dimension_vector(dimensions_to_use)

  # Read metadata now for batch_column validation and end-of-run plot detection.
  # This is refreshed after each sub-checkpoint load below.
  metadata <- tibble::as_tibble(.giotto_pdata_dt(gobj))
  if (!batch_column %in% names(metadata)) {
    fallback_column <- c("slide_id", "slide_num", "list_ID", "sample_id")[
      c("slide_id", "slide_num", "list_ID", "sample_id") %in% names(metadata)
    ][1]
    if (is.na(fallback_column) || !nzchar(fallback_column)) {
      stop("Batch column '", batch_column, "' is not present in merged metadata.")
    }
    batch_column <- fallback_column
  }

  # Hard guard: Harmony must never regress out biological variables.
  forbidden_batch <- c("group_id", "treatment", "group", "arm")
  if (tolower(batch_column) %in% forbidden_batch) {
    stop(sprintf(
      "batch_column '%s' is a biological variable; Harmony would erase CART/Control/Conventional signal. Use slide_id / slide_num / sample_id instead.",
      batch_column))
  }

  # Advisory: warn if batch_column is strongly correlated with group_id.
  if ("group_id" %in% names(metadata) && batch_column != "group_id") {
    tab <- table(metadata[[batch_column]], metadata[["group_id"]])
    if (all(dim(tab) >= 2L)) {
      chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
      n   <- sum(tab)
      k   <- min(dim(tab)) - 1L
      cramers_v <- if (n > 0 && k > 0) {
        sqrt(as.numeric(chi$statistic) / (n * k))
      } else NA_real_
      if (is.finite(cramers_v) && cramers_v > 0.5) {
        warning(sprintf(
          "batch_column '%s' has Cramer's V = %.2f with group_id (> 0.5); Harmony may partially erase biological signal.",
          batch_column, cramers_v))
      }
    }
  }

  # Sub-checkpoints — saved inside output_dir alongside the main checkpoints.
  # Each represents one completed stage so a crash mid-run can resume from
  # the last successful save rather than restarting from Giotto_Object_Merged.
  #
  #   Stage 0  Giotto_Object_Merged           (input — created by merge step)
  #   Stage 1  Giotto_Object_Merged_PCA       (after HVG selection + PCA)
  #   Stage 2  Giotto_Object_Merged_Harmony   (after Harmony batch correction)
  #   Stage 3  Giotto_Object_Merged_Embeddings (after UMAP + tSNE)
  #   Stage 4  Giotto_Object_BatchCorrected   (after NN + Leiden — final output)

  subckpt_pca        <- file.path(output_dir, "Giotto_Object_Merged_PCA")
  subckpt_harmony    <- file.path(output_dir, "Giotto_Object_Merged_Harmony")
  subckpt_embeddings <- file.path(output_dir, "Giotto_Object_Merged_Embeddings")

  resume_stage <- 0L
  if (!isTRUE(force_recompute)) {
    if (.sub_ckpt_exists(subckpt_embeddings)) {
      resume_stage <- 3L
      cat("\u21a9 Sub-checkpoint found — resuming after UMAP/tSNE\n")
      cat("  Loading:", subckpt_embeddings, "\n\n")
      gobj     <- load_giotto_checkpoint(subckpt_embeddings)
      metadata <- tibble::as_tibble(.giotto_pdata_dt(gobj))
    } else if (.sub_ckpt_exists(subckpt_harmony)) {
      resume_stage <- 2L
      cat("\u21a9 Sub-checkpoint found — resuming after Harmony\n")
      cat("  Loading:", subckpt_harmony, "\n\n")
      gobj     <- load_giotto_checkpoint(subckpt_harmony)
      metadata <- tibble::as_tibble(.giotto_pdata_dt(gobj))
    } else if (.sub_ckpt_exists(subckpt_pca)) {
      resume_stage <- 1L
      cat("\u21a9 Sub-checkpoint found — resuming after PCA\n")
      cat("  Loading:", subckpt_pca, "\n\n")
      gobj     <- load_giotto_checkpoint(subckpt_pca)
      metadata <- tibble::as_tibble(.giotto_pdata_dt(gobj))
    }
  }

  if (resume_stage > 0L) {
    cat("Skipping completed stages (1-", resume_stage, "). ",
        "Pass force_recompute = TRUE to re-run from scratch.\n\n", sep = "")
  }

  # ── Stage 1: HVG selection + PCA ──────────────────────────────────────────
  if (resume_stage < 1L) {
    cat("Selecting HVGs for merged PCA...\n")
    hvg_genes <- tryCatch(
      .select_sparse_hvgs(
        gobj           = gobj,
        n_hvgs         = n_hvgs,
        sample_id      = sample_id,
        results_folder = results_dir
      ),
      error = function(e) {
        cat("\u26A0 HVG selection failed (", conditionMessage(e), ") \u2014 using all features\n")
        NULL
      }
    )
    cat("Running PCA on merged object (",
        length(hvg_genes) %||% "all", " features, ", n_pcs, " PCs)...\n", sep = "")
    gobj <- .giotto_run_pca(gobj, feats_to_use = hvg_genes, n_pcs = n_pcs)
    cat("\u2713 PCA complete\n\n")
    save_giotto_checkpoint(gobj, subckpt_pca,
                           metadata = list(stage = "pca", sample_id = sample_id))
    cat("\u2713 Sub-checkpoint saved:", subckpt_pca, "\n\n")
  }

  # ── Stage 2: Harmony ──────────────────────────────────────────────────────
  if (resume_stage < 2L) {
    cat("Running Harmony using batch column:", batch_column, "\n")
    cat("Dimensions:", paste(range(dimensions_to_use), collapse = ":"), "\n\n")
    gobj <- runGiottoHarmony(
      gobject              = gobj,
      vars_use             = batch_column,
      dim_reduction_to_use = "pca",
      dimensions_to_use    = dimensions_to_use,
      name                 = "harmony"
    )
    cat("\u2713 Harmony complete\n\n")
    save_giotto_checkpoint(gobj, subckpt_harmony,
                           metadata = list(stage = "harmony", sample_id = sample_id,
                                           batch_column = batch_column))
    cat("\u2713 Sub-checkpoint saved:", subckpt_harmony, "\n\n")
  }

  # ── Stage 3: UMAP + tSNE ─────────────────────────────────────────────────
  # Call Giotto dim-reduction functions DIRECTLY (no wrappers) — they use
  # sys.frame(-2) internally; an extra wrapper frame shifts the lookup to the
  # wrong environment and causes "accessor not found". (Do-Not-Repeat rule)
  if (resume_stage < 3L) {
    cat("Running UMAP on Harmony embedding...\n")
    gobj <- .run_known_giotto_warning_safe(
      runUMAP(
        gobject              = gobj,
        dim_reduction_to_use = "harmony",
        dim_reduction_name   = "harmony",
        dimensions_to_use    = dimensions_to_use,
        n_neighbors          = umap_n_neighbors,
        min_dist             = umap_min_dist,
        name                 = "umap"
      )
    )
    cat("\u2713 UMAP complete\n\n")

    cat("Running t-SNE on Harmony embedding...\n")
    gobj <- .run_known_giotto_warning_safe(
      runtSNE(
        gobject              = gobj,
        dim_reduction_to_use = "harmony",
        dim_reduction_name   = "harmony",
        dimensions_to_use    = dimensions_to_use,
        perplexity           = 30,
        name                 = "tsne"
      )
    )
    cat("\u2713 t-SNE complete\n\n")
    save_giotto_checkpoint(gobj, subckpt_embeddings,
                           metadata = list(stage = "embeddings", sample_id = sample_id))
    cat("\u2713 Sub-checkpoint saved:", subckpt_embeddings, "\n\n")
  }

  # ── Stage 4: NN network + Leiden ─────────────────────────────────────────
  cat("Building nearest-neighbour network on Harmony embedding...\n")
  gobj <- createNearestNetwork(
    gobject              = gobj,
    dim_reduction_to_use = "harmony",
    dim_reduction_name   = "harmony",
    dimensions_to_use    = dimensions_to_use,
    k                    = k_nn,
    name                 = "NN.harmony"
  )
  cat("\u2713 NN network complete\n\n")

  cat("Running Leiden clustering...\n")
  gobj <- doLeidenCluster(
    gobject           = gobj,
    nn_network_to_use = "sNN",
    network_name      = "NN.harmony",
    resolution        = resolution,
    n_iterations      = 1000,
    name              = "leiden_clust"
  )
  cat("\u2713 Leiden clustering complete\n\n")

  readr::write_csv(
    tibble::as_tibble(.giotto_pdata_dt(gobj)),
    file.path(results_dir, paste0(sample_id, "_batch_corrected_metadata.csv"))
  )

  if (create_plots) {
    batch_plot_cols <- c(batch_column, "leiden_clust")

    pca_df     <- embedding_to_tibble(gobj, reduction_method = "pca",
                                      name = NULL,      color_columns = batch_plot_cols)
    harmony_df <- embedding_to_tibble(gobj, reduction_method = "harmony",
                                      name = "harmony", color_columns = batch_plot_cols)
    umap_df    <- embedding_to_tibble(gobj, reduction_method = "umap",
                                      name = "umap",    color_columns = batch_plot_cols)

    save_embedding_plot(pca_df,
      color_column = batch_column,
      title = paste(sample_id, "- PCA by", batch_column, "(before Harmony)"),
      output_file  = file.path(results_dir, paste0(sample_id, "_pca_by_", batch_column, ".png")))
    save_embedding_plot(harmony_df,
      color_column = batch_column,
      title = paste(sample_id, "- Harmony dims by", batch_column),
      output_file  = file.path(results_dir, paste0(sample_id, "_harmony_by_", batch_column, ".png")))
    save_embedding_plot(umap_df,
      color_column = batch_column,
      title = paste(sample_id, "- UMAP by", batch_column),
      output_file  = file.path(results_dir, paste0(sample_id, "_umap_by_", batch_column, ".png")))
    save_embedding_plot(umap_df,
      color_column = "leiden_clust",
      title = paste(sample_id, "- UMAP by Leiden cluster"),
      output_file  = file.path(results_dir, paste0(sample_id, "_umap_by_leiden_clust.png")))

    # Validation: UMAP coloured by cell type — confirms Harmony preserved biology.
    if (is.null(annotation_columns)) {
      annotation_columns <- grep("^celltype_", names(metadata), value = TRUE)
    }
    ann_cols_present <- annotation_columns[annotation_columns %in% names(metadata)]

    if (length(ann_cols_present) > 0) {
      cat("Generating annotation validation plots for:",
          paste(ann_cols_present, collapse = ", "), "\n")
      umap_ann <- embedding_to_tibble(gobj, reduction_method = "umap", name = "umap",
                                      color_columns = ann_cols_present)
      for (ann_col in ann_cols_present) {
        save_embedding_plot(umap_ann,
          color_column = ann_col,
          title = paste(sample_id, "- UMAP by", ann_col, "(biology check)"),
          output_file  = file.path(results_dir,
                                   paste0(sample_id, "_umap_by_", ann_col, ".png")))
      }
    }

    if ("sample_id" %in% names(metadata) && "sample_id" != batch_column) {
      umap_sample <- embedding_to_tibble(gobj, reduction_method = "umap", name = "umap",
                                         color_columns = "sample_id")
      save_embedding_plot(umap_sample,
        color_column = "sample_id",
        title = paste(sample_id, "- UMAP by sample_id (mixing check)"),
        output_file  = file.path(results_dir,
                                 paste0(sample_id, "_umap_by_sample_id.png")))
    }
  }

  if (save_object) {
    save_giotto_checkpoint(
      gobj           = gobj,
      checkpoint_dir = file.path(output_dir, "Giotto_Object_BatchCorrected"),
      metadata       = list(stage = "batch_correction", batch_column = batch_column)
    )
  }

  gobj
}

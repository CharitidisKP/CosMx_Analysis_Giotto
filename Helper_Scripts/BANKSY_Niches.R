#!/usr/bin/env Rscript
# ==============================================================================
# BANKSY_Niches.R
# Spatial-domain niche labels via BANKSY (Singhal et al., Nat Genet 2024).
#
# Pipeline (Stage 1, Phase 7):
#   1. compute_banksy_per_sample() - per-sample Giotto -> SpatialExperiment ->
#      Banksy::computeBanksy() with k_geom + lambda. Returns the augmented
#      matrix list keyed on sample_id (each cell has neighborhood-mean +
#      AGF-augmented features alongside the per-cell expression).
#   2. cluster_banksy_niches() - cbind augmented matrices in the order of
#      the merged Giotto's cell_IDs -> PCA -> Harmony on sample_id -> Leiden.
#      Writes the resulting `niche_id` column on the merged Giotto's cell
#      metadata and a niche_summary.csv (cells per niche per sample).
#
# Required (optional) packages:
#   - Banksy        (Bioconductor)
#   - SpatialExperiment (Bioconductor)
#   - SummarizedExperiment / SingleCellExperiment (transitive deps)
#   - Giotto (already in renv)
#   - harmony, igraph (already used elsewhere in pipeline)
#
# All Banksy / Harmony calls are wrapped in tryCatch so a missing dep or
# computational failure produces a clear log line and skips the niche step
# without stopping the pipeline.
# ==============================================================================

if (!exists(".giotto_pdata_dt")) {
  .giotto_pdata_dt <- function(gobj) {
    if (requireNamespace("Giotto", quietly = TRUE) &&
        exists("pDataDT", envir = asNamespace("Giotto"), inherits = FALSE)) {
      get("pDataDT", envir = asNamespace("Giotto"))(gobj)
    } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
               exists("pDataDT", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
      get("pDataDT", envir = asNamespace("GiottoClass"))(gobj)
    } else {
      get("pDataDT", mode = "function")(gobj)
    }
  }
}

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}

# Pulls the normalized expression matrix (genes by cells), BANKSY's neighbour-mean
# smoothing distorts on raw counts because library-size differences dominate.
.banksy_get_norm_counts <- function(gobj) {
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("getExpression", envir = asNamespace("Giotto"), inherits = FALSE)) {
    get("getExpression", envir = asNamespace("Giotto"))(
      gobj, values = "normalized", output = "matrix")
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getExpression", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    get("getExpression", envir = asNamespace("GiottoClass"))(
      gobj, values = "normalized", output = "matrix")
  } else NULL
}

# Pulls per-cell spatial coordinates (sdimx, sdimy) from a Giotto object.
.banksy_get_spatial_locs <- function(gobj) {
  acc <- if (requireNamespace("Giotto", quietly = TRUE) &&
             exists("getSpatialLocations", envir = asNamespace("Giotto"), inherits = FALSE)) {
    get("getSpatialLocations", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getSpatialLocations", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    get("getSpatialLocations", envir = asNamespace("GiottoClass"))
  } else return(NULL)
  tryCatch(acc(gobj, output = "data.table"), error = function(e) NULL)
}

# Builds a SpatialExperiment for one sample from a Giotto object.
.banksy_giotto_to_spe <- function(gobj, sample_id) {
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("SpatialExperiment is required for BANKSY. ",
         "BiocManager::install('SpatialExperiment').")
  }
  counts <- .banksy_get_norm_counts(gobj)
  locs   <- .banksy_get_spatial_locs(gobj)
  if (is.null(counts) || is.null(locs)) return(NULL)
  # Match cell_ID column order across counts and locs.
  cell_ids <- colnames(counts)
  locs_dt <- as.data.frame(locs)
  locs_dt <- locs_dt[match(cell_ids, locs_dt$cell_ID), , drop = FALSE]
  coords <- as.matrix(locs_dt[, c("sdimx", "sdimy")])
  rownames(coords) <- cell_ids
  # Keep colData minimal - only cell_ID + sample_id. Per-sample Giotto
  # checkpoints carry independent metadata columns (annotation/cluster
  # outputs differ between samples), and SummarizedExperiment::cbind
  # requires identical colData column names across samples. Rich metadata
  # is re-attached from the merged Giotto object after clustering.
  meta <- data.frame(
    cell_ID   = cell_ids,
    sample_id = sample_id,
    stringsAsFactors = FALSE
  )
  rownames(meta) <- cell_ids
  SpatialExperiment::SpatialExperiment(
    assays      = list(normalized = counts),
    colData     = S4Vectors::DataFrame(meta),
    spatialCoords = coords,
    sample_id   = sample_id
  )
}

#' Per-sample BANKSY augmentation. Returns a named list of augmented
#' SpatialExperiment objects (one per sample) ready for cbind in
#' `cluster_banksy_niches()`.
#'
#' @param per_sample_gobjs Named list of per-sample Giotto objects (the
#'   same pre-merge checkpoints used by `load_merge_inputs()`).
#' @param k_geom Spatial-neighborhood size (BANKSY default 18).
#' @param lambda Spatial weight (0 = expression only, 1 = spatial only;
#'   BANKSY default 0.2 for cell-typing-style niches).
compute_banksy_per_sample <- function(per_sample_gobjs,
                                      k_geom = 18,
                                      lambda = 0.2) {
  if (!requireNamespace("Banksy", quietly = TRUE)) {
    cat("Warning: Banksy not installed - skipping BANKSY niche computation.\n")
    return(NULL)
  }

  spes <- list()
  for (nm in names(per_sample_gobjs)) {
    cat("  BANKSY: per-sample augmentation for ", nm, "\n", sep = "")
    spe <- tryCatch(
      .banksy_giotto_to_spe(per_sample_gobjs[[nm]], nm),
      error = function(e) {
        cat(sprintf("    Warning: SPE conversion failed: %s\n",
                    conditionMessage(e)))
        NULL
      }
    )
    if (is.null(spe)) next
    aug <- tryCatch(
      Banksy::computeBanksy(spe,
                            assay_name = "normalized",
                            k_geom = k_geom,
                            compute_agf = TRUE),
      error = function(e) {
        cat(sprintf("    Warning: computeBanksy failed: %s\n",
                    conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(aug)) spes[[nm]] <- aug
  }
  if (length(spes) == 0L) return(NULL)
  spes
}

#' Concatenate augmented SPEs, run PCA -> Harmony -> Leiden, write
#' `niche_id` back to the merged Giotto cell metadata.
#'
#' @param merged_gobj Merged Giotto object (post batch_correct_merged_object).
#' @param banksy_aug_list Output of compute_banksy_per_sample().
#' @param batch_col Column for Harmony correction (default "sample_id").
#' @param resolution Leiden resolution.
#' @param output_dir Merged-pipeline output dir; niche_summary.csv and the
#'   niche composition heatmap are written to <output_dir>/14_BANKSY/BANKSY/.
cluster_banksy_niches <- function(merged_gobj,
                                  banksy_aug_list,
                                  batch_col = "sample_id",
                                  resolution = 0.4,
                                  lambda = 0.2,
                                  output_dir,
                                  n_pcs = 30,
                                  k_nn = 30) {
  if (is.null(banksy_aug_list)) {
    cat("  BANKSY: no augmented data - skipping niche clustering.\n")
    return(invisible(NULL))
  }
  if (!requireNamespace("Banksy", quietly = TRUE)) return(invisible(NULL))

  # cbind the per-sample augmented SPEs, keep only cells present in the
  # merged Giotto (defensive: joinGiottoObjects may have dropped some).
  aug <- do.call(SummarizedExperiment::cbind, unname(banksy_aug_list))
  merged_meta <- as.data.frame(.giotto_pdata_dt(merged_gobj))
  merged_ids <- as.character(merged_meta$cell_ID)
  keep <- colnames(aug) %in% merged_ids
  aug <- aug[, keep]
  # Drop NAs from the match index BEFORE subsetting. Sparse matrices reject NA indices.
  idx <- match(merged_ids, colnames(aug))
  idx <- idx[!is.na(idx)]
  aug <- aug[, idx, drop = FALSE]
  if (ncol(aug) == 0L) {
    cat("  BANKSY: no overlapping cells with merged object - skipping.\n")
    return(invisible(NULL))
  }

  # PCA on the BANKSY-augmented assay. Banksy::runBanksyPCA wraps this.
  pca <- tryCatch(
    Banksy::runBanksyPCA(aug, lambda = lambda, npcs = n_pcs),
    error = function(e) {
      cat(sprintf("    Warning: runBanksyPCA failed: %s\n",
                  conditionMessage(e)))
      NULL
    }
  )
  if (is.null(pca)) return(invisible(NULL))
  expected_slot <- paste0("PCA_M0_lam", lambda)
  pca_mat <- SingleCellExperiment::reducedDim(pca, expected_slot)
  if (is.null(pca_mat)) {
    rd_names <- SingleCellExperiment::reducedDimNames(pca)
    cat(sprintf("    Warning: BANKSY PCA slot '%s' not found; using fallback '%s' (available: %s).\n",
                expected_slot, rd_names[1], paste(rd_names, collapse = ", ")))
    pca_mat <- SingleCellExperiment::reducedDim(pca, rd_names[1])
  }

  # Harmony on sample_id. Use the harmony R package directly.
  if (!requireNamespace("harmony", quietly = TRUE)) {
    cat("    Warning: harmony not installed - skipping batch correction step.\n")
    harmony_mat <- pca_mat
  } else {
    harmony_mat <- tryCatch(
      harmony::HarmonyMatrix(
        data_mat = pca_mat,
        meta_data = data.frame(batch = merged_meta[[batch_col]]),
        vars_use = "batch",
        do_pca = FALSE,
        verbose = FALSE
      ),
      error = function(e) {
        cat(sprintf("    Warning: Harmony failed: %s - using uncorrected PCs.\n",
                    conditionMessage(e)))
        pca_mat
      }
    )
  }

  # Leiden on a kNN graph in Harmony space.
  if (!requireNamespace("igraph", quietly = TRUE)) {
    cat("    Warning: igraph not installed - cannot cluster niches.\n")
    return(invisible(NULL))
  }
  knn <- tryCatch(
    BiocNeighbors::findKNN(harmony_mat, k = k_nn),
    error = function(e) {
      # Fallback: build kNN with stats::dist (slower but always available).
      d <- as.matrix(stats::dist(harmony_mat))
      idx <- t(apply(d, 1, function(r) order(r)[2:(k_nn + 1L)]))
      list(index = idx)
    }
  )
  el <- do.call(rbind, lapply(seq_len(nrow(knn$index)), function(i) {
    cbind(i, knn$index[i, ])
  }))
  g <- igraph::graph_from_edgelist(el, directed = FALSE)
  g <- igraph::simplify(g)
  comm <- tryCatch(
    igraph::cluster_leiden(g,
                           objective_function = "modularity",
                           resolution_parameter = resolution),
    error = function(e) {
      cat("    Warning: Leiden failed (", conditionMessage(e),
          "); falling back to Louvain.\n", sep = "")
      igraph::cluster_louvain(g)
    }
  )
  niche_id <- paste0("niche_", igraph::membership(comm))

  # Push niche_id to the merged Giotto cell metadata.
  acc <- if (requireNamespace("Giotto", quietly = TRUE) &&
             exists("addCellMetadata", envir = asNamespace("Giotto"), inherits = FALSE)) {
    get("addCellMetadata", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("addCellMetadata", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    get("addCellMetadata", envir = asNamespace("GiottoClass"))
  } else {
    get("addCellMetadata", mode = "function")
  }
  niche_df <- data.frame(cell_ID = colnames(aug), niche_id = niche_id,
                         stringsAsFactors = FALSE)
  merged_gobj <- tryCatch(
    acc(merged_gobj, new_metadata = niche_df,
        by_column = TRUE, column_cell_ID = "cell_ID"),
    error = function(e) {
      cat(sprintf("    Warning: addCellMetadata failed for niche_id: %s\n",
                  conditionMessage(e)))
      merged_gobj
    }
  )

  # niche_summary.csv + composition heatmap: cells per niche per sample.
  out_dir <- file.path(output_dir, "14_BANKSY")
  banksy_dir <- file.path(out_dir, "BANKSY")
  dir.create(banksy_dir, recursive = TRUE, showWarnings = FALSE)
  summary_df <- as.data.frame(table(
    sample_id = merged_meta[match(niche_df$cell_ID, merged_ids), batch_col],
    niche_id  = niche_df$niche_id
  ))
  utils::write.csv(summary_df,
                   file.path(banksy_dir, "niche_summary.csv"),
                   row.names = FALSE)
  cat(sprintf("  BANKSY niche labels written; %d niches across %d samples.\n",
              length(unique(niche_id)),
              length(unique(merged_meta[[batch_col]]))))

  merged_meta_with_id <- merged_meta
  if (!"cell_ID" %in% names(merged_meta_with_id)) {
    merged_meta_with_id$cell_ID <- merged_ids
  }
  .banksy_plot_niche_composition(merged_meta_with_id, niche_df,
                                  sample_id = "merged", out_dir = out_dir)

  invisible(merged_gobj)
}

#' Single-sample BANKSY clustering. Mirrors `cluster_banksy_niches()` but
#' skips the cbind + Harmony steps (no batch in single-sample mode) and
#' writes `niche_id` to the per-sample Giotto object.
cluster_banksy_single_sample <- function(aug_spe,
                                         sample_id,
                                         per_sample_gobj,
                                         lambda = 0.2,
                                         resolution = 0.4,
                                         n_pcs = 30,
                                         k_nn = 30,
                                         output_dir) {
  if (is.null(aug_spe)) {
    cat("  BANKSY: no augmented SPE for ", sample_id, ", skipping.\n",
        sep = "")
    return(invisible(NULL))
  }
  if (!requireNamespace("Banksy", quietly = TRUE)) return(invisible(NULL))

  pca <- tryCatch(
    Banksy::runBanksyPCA(aug_spe, lambda = lambda, npcs = n_pcs),
    error = function(e) {
      cat("    Warning: runBanksyPCA failed: ", conditionMessage(e), "\n",
          sep = "")
      NULL
    }
  )
  if (is.null(pca)) return(invisible(NULL))
  expected_slot <- paste0("PCA_M0_lam", lambda)
  pca_mat <- SingleCellExperiment::reducedDim(pca, expected_slot)
  if (is.null(pca_mat)) {
    rd_names <- SingleCellExperiment::reducedDimNames(pca)
    cat(sprintf("    Warning: BANKSY PCA slot '%s' not found; using fallback '%s' (available: %s).\n",
                expected_slot, rd_names[1], paste(rd_names, collapse = ", ")))
    pca_mat <- SingleCellExperiment::reducedDim(pca, rd_names[1])
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    cat("    Warning: igraph not installed, cannot cluster niches.\n")
    return(invisible(NULL))
  }
  knn <- tryCatch(
    BiocNeighbors::findKNN(pca_mat, k = k_nn),
    error = function(e) {
      d <- as.matrix(stats::dist(pca_mat))
      idx <- t(apply(d, 1, function(r) order(r)[2:(k_nn + 1L)]))
      list(index = idx)
    }
  )
  el <- do.call(rbind, lapply(seq_len(nrow(knn$index)), function(i) {
    cbind(i, knn$index[i, ])
  }))
  g <- igraph::simplify(igraph::graph_from_edgelist(el, directed = FALSE))
  comm <- tryCatch(
    igraph::cluster_leiden(g,
                           objective_function = "modularity",
                           resolution_parameter = resolution),
    error = function(e) {
      cat("    Warning: Leiden failed (", conditionMessage(e),
          "), falling back to Louvain.\n", sep = "")
      igraph::cluster_louvain(g)
    }
  )
  niche_id <- paste0("niche_", igraph::membership(comm))

  acc <- if (requireNamespace("Giotto", quietly = TRUE) &&
             exists("addCellMetadata", envir = asNamespace("Giotto"),
                    inherits = FALSE)) {
    get("addCellMetadata", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("addCellMetadata", envir = asNamespace("GiottoClass"),
                    inherits = FALSE)) {
    get("addCellMetadata", envir = asNamespace("GiottoClass"))
  } else {
    get("addCellMetadata", mode = "function")
  }
  niche_df <- data.frame(cell_ID = colnames(aug_spe), niche_id = niche_id,
                         stringsAsFactors = FALSE)
  per_sample_gobj <- tryCatch(
    acc(per_sample_gobj, new_metadata = niche_df,
        by_column = TRUE, column_cell_ID = "cell_ID"),
    error = function(e) {
      cat("    Warning: addCellMetadata failed for niche_id: ",
          conditionMessage(e), "\n", sep = "")
      per_sample_gobj
    }
  )

  out_dir <- file.path(output_dir, "14_BANKSY")
  banksy_dir <- file.path(out_dir, "BANKSY")
  dir.create(banksy_dir, recursive = TRUE, showWarnings = FALSE)
  summary_df <- as.data.frame(table(niche_id = niche_df$niche_id))
  utils::write.csv(summary_df,
                   file.path(banksy_dir,
                             paste0(sample_id, "_niche_summary.csv")),
                   row.names = FALSE)
  cat(sprintf("  BANKSY niche labels (%d niches) written for %s.\n",
              length(unique(niche_id)), sample_id))

  .banksy_plot_niche_polygons(per_sample_gobj, niche_df, sample_id, out_dir)
  per_sample_meta <- tryCatch(as.data.frame(.giotto_pdata_dt(per_sample_gobj)),
                              error = function(e) NULL)
  .banksy_plot_niche_composition(per_sample_meta, niche_df, sample_id, out_dir)

  invisible(per_sample_gobj)
}

# Polygon-coloured niche plot. Reuses .extract_polygon_df + presentation_theme + save_presentation_plot from the canonical Pipeline_Utils / Plot_Helpers stack. Falls back to a centroid scatter if polygons are unavailable.
.banksy_plot_niche_polygons <- function(gobj, niche_df, sample_id, out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(NULL))
  poly_df <- tryCatch(.extract_polygon_df(gobj), error = function(e) NULL)
  banksy_dir <- file.path(out_dir, "BANKSY")
  ensure_dir(banksy_dir)
  out_path <- file.path(banksy_dir, paste0(sample_id, "_niche_spatial.png"))

  niche_lookup <- setNames(niche_df$niche_id, niche_df$cell_ID)
  niche_levels <- paste0("niche_", sort(as.integer(sub("^niche_", "",
                                                        unique(niche_df$niche_id)))))
  n_niches <- length(niche_levels)

  p <- if (!is.null(poly_df) && nrow(poly_df) > 0) {
    poly_df$poly_group <- paste(poly_df$cell_ID, poly_df$geom,
                                 poly_df$part, sep = "_")
    poly_df$niche_id <- factor(niche_lookup[as.character(poly_df$cell_ID)],
                                levels = niche_levels)
    poly_df <- poly_df[!is.na(poly_df$niche_id), , drop = FALSE]
    ggplot2::ggplot(poly_df,
                    ggplot2::aes(x = x, y = y,
                                 group = poly_group, fill = niche_id)) +
      ggplot2::geom_polygon(colour = "grey30", linewidth = 0.05) +
      ggplot2::scale_fill_discrete(labels = function(x) pretty_plot_label(x))
  } else {
    locs <- tryCatch(.banksy_get_spatial_locs(gobj), error = function(e) NULL)
    if (is.null(locs)) return(invisible(NULL))
    locs_df <- as.data.frame(locs)
    locs_df$niche_id <- factor(niche_lookup[as.character(locs_df$cell_ID)],
                                levels = niche_levels)
    locs_df <- locs_df[!is.na(locs_df$niche_id), , drop = FALSE]
    ggplot2::ggplot(locs_df,
                    ggplot2::aes(x = sdimx, y = sdimy, colour = niche_id)) +
      ggplot2::geom_point(size = 0.4) +
      ggplot2::scale_colour_discrete(labels = function(x) pretty_plot_label(x))
  }

  fallback_tag <- if (is.null(poly_df) || !nrow(poly_df)) " (centroid fallback)" else ""
  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = sample_plot_title(sample_id,
        paste0("BANKSY niches (", n_niches, ")", fallback_tag)),
      x = "Global X Coordinate", y = "Global Y Coordinate",
      fill = "Niche", colour = "Niche"
    ) +
    presentation_theme(base_size = 12, legend_position = "right") +
    ggplot2::theme(
      legend.key.height = grid::unit(0.45, "cm"),
      axis.title.x      = element_markdown_safe(margin = ggplot2::margin(t = 10)),
      axis.title.y      = element_markdown_safe(margin = ggplot2::margin(r = 10)),
      plot.margin       = ggplot2::margin(t = 10, r = 20, b = 20, l = 20)
    )

  tryCatch({
    save_presentation_plot(p, out_path, width = 10, height = 8, dpi = 150)
    cat(sprintf("  BANKSY niche spatial plot saved: %s\n", out_path))
  }, error = function(e) {
    cat(sprintf("    Warning: niche spatial plot failed: %s\n",
                conditionMessage(e)))
  })
  invisible(out_path)
}

# Auto-resolves a celltype column from a Giotto pData frame. Prefers refined supervised, then supervised, then semi-supervised, then leiden_clust.
.banksy_resolve_celltype_col <- function(meta) {
  if (is.null(meta) || !length(meta)) return(NA_character_)
  candidates <- c(
    "celltype",
    grep("^celltype_.+_supervised_refined$", names(meta), value = TRUE),
    grep("^celltype_.+_supervised$", names(meta), value = TRUE),
    grep("^celltype_.+_semi$", names(meta), value = TRUE),
    "leiden_clust"
  )
  for (col in candidates) if (col %in% names(meta)) return(col)
  NA_character_
}

# Niche x celltype composition heatmap. Builds a log2 fold-enrichment matrix (observed niche-level fraction / sample-wide fraction) and delegates rendering to the shared .plot_proximity_heatmap_complex(). Cells with zero counts become NA so the helper greys them via na_col.
.banksy_plot_niche_composition <- function(meta, niche_df, sample_id, out_dir,
                                            celltype_col = NULL) {
  if (is.null(meta) || !nrow(meta) || is.null(niche_df) || !nrow(niche_df)) {
    return(invisible(NULL))
  }
  if (is.null(celltype_col) || !nzchar(celltype_col)) {
    celltype_col <- .banksy_resolve_celltype_col(meta)
  }
  if (is.na(celltype_col) || !celltype_col %in% names(meta)) {
    cat("  BANKSY: no celltype column found; skipping composition heatmap.\n")
    return(invisible(NULL))
  }
  if (!exists(".plot_proximity_heatmap_complex", mode = "function")) {
    cat("  BANKSY: .plot_proximity_heatmap_complex unavailable; skipping composition heatmap.\n")
    return(invisible(NULL))
  }

  niche_lookup <- setNames(as.character(niche_df$niche_id),
                            as.character(niche_df$cell_ID))
  meta <- as.data.frame(meta)
  meta$niche_id <- niche_lookup[as.character(meta$cell_ID)]
  meta <- meta[!is.na(meta$niche_id) & !is.na(meta[[celltype_col]]), , drop = FALSE]
  if (!nrow(meta)) return(invisible(NULL))

  count_mat <- as.matrix(unclass(table(
    niche_id = as.character(meta$niche_id),
    celltype = as.character(meta[[celltype_col]])
  )))
  if (!nrow(count_mat) || !ncol(count_mat)) return(invisible(NULL))

  niche_totals <- rowSums(count_mat)
  obs_frac <- sweep(count_mat, 1, niche_totals, FUN = "/")
  global_totals <- colSums(count_mat)
  exp_frac <- global_totals / sum(global_totals)
  log2fe <- sweep(log2(obs_frac), 2, log2(exp_frac), FUN = "-")
  log2fe[!is.finite(log2fe)] <- NA_real_

  niche_order <- order(suppressWarnings(as.integer(sub("^niche_", "",
                                                       rownames(log2fe)))),
                        rownames(log2fe))
  log2fe <- log2fe[niche_order, , drop = FALSE]

  banksy_dir <- file.path(out_dir, "BANKSY")
  ensure_dir(banksy_dir)
  out_path <- file.path(banksy_dir,
                        paste0(sample_id, "_niche_composition_heatmap.png"))

  tryCatch({
    .plot_proximity_heatmap_complex(
      mat       = log2fe,
      title     = sample_plot_title(sample_id, "BANKSY niche composition"),
      subtitle  = "log2 fold-enrichment vs sample-wide cell-type fractions",
      filename  = out_path,
      width     = max(10, 0.6 * ncol(log2fe) + 4),
      height    = max(5, 0.5 * nrow(log2fe) + 3),
      dpi       = 300
    )
    cat(sprintf("  BANKSY niche composition heatmap saved: %s\n", out_path))
  }, error = function(e) {
    cat(sprintf("    Warning: niche composition heatmap failed: %s\n",
                conditionMessage(e)))
  })
  invisible(out_path)
}

#' Mode-dispatch wrapper. `mode = "per_sample"` (default) loops
#' `cluster_banksy_single_sample()` per sample. `mode = "merged"` calls the
#' existing `cluster_banksy_niches()` after concatenating the augmented SPEs.
run_banksy <- function(per_sample_gobjs,
                       mode = c("per_sample", "merged"),
                       merged_gobj = NULL,
                       output_dir,
                       k_geom = 18,
                       lambda = 0.2,
                       resolution = 0.4,
                       n_pcs = 30,
                       k_nn = 30,
                       batch_col = "sample_id") {
  mode <- match.arg(mode)
  aug <- compute_banksy_per_sample(per_sample_gobjs, k_geom = k_geom,
                                    lambda = lambda)
  if (is.null(aug)) return(invisible(NULL))

  if (mode == "per_sample") {
    out <- list()
    for (nm in names(aug)) {
      out[[nm]] <- cluster_banksy_single_sample(
        aug_spe          = aug[[nm]],
        sample_id        = nm,
        per_sample_gobj  = per_sample_gobjs[[nm]],
        lambda           = lambda,
        resolution       = resolution,
        n_pcs            = n_pcs,
        k_nn             = k_nn,
        output_dir       = output_dir
      )
    }
    invisible(out)
  } else {
    if (is.null(merged_gobj))
      stop("run_banksy(mode = 'merged') requires merged_gobj.")
    cluster_banksy_niches(
      merged_gobj     = merged_gobj,
      banksy_aug_list = aug,
      batch_col       = batch_col,
      resolution      = resolution,
      lambda          = lambda,
      output_dir      = output_dir,
      n_pcs           = n_pcs,
      k_nn            = k_nn
    )
  }
}

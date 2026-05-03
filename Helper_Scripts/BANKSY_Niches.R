#!/usr/bin/env Rscript
# ==============================================================================
# BANKSY_Niches.R
# Spatial-domain niche labels via BANKSY (Singhal et al., Nat Genet 2024).
#
# Pipeline (Stage 1, Phase 7):
#   1. compute_banksy_per_sample() — per-sample Giotto → SpatialExperiment →
#      Banksy::computeBanksy() with k_geom + lambda. Returns the augmented
#      matrix list keyed on sample_id (each cell has neighborhood-mean +
#      AGF-augmented features alongside the per-cell expression).
#   2. cluster_banksy_niches() — cbind augmented matrices in the order of
#      the merged Giotto's cell_IDs → PCA → Harmony on sample_id → Leiden.
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

# Pulls the raw expression matrix (genes × cells) from a Giotto object via
# whichever accessor is available (Giotto vs GiottoClass).
.banksy_get_raw_counts <- function(gobj) {
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("getExpression", envir = asNamespace("Giotto"), inherits = FALSE)) {
    get("getExpression", envir = asNamespace("Giotto"))(
      gobj, values = "raw", output = "matrix")
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getExpression", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    get("getExpression", envir = asNamespace("GiottoClass"))(
      gobj, values = "raw", output = "matrix")
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
  counts <- .banksy_get_raw_counts(gobj)
  locs   <- .banksy_get_spatial_locs(gobj)
  if (is.null(counts) || is.null(locs)) return(NULL)
  # Match cell_ID column order across counts and locs.
  cell_ids <- colnames(counts)
  locs_dt <- as.data.frame(locs)
  locs_dt <- locs_dt[match(cell_ids, locs_dt$cell_ID), , drop = FALSE]
  coords <- as.matrix(locs_dt[, c("sdimx", "sdimy")])
  rownames(coords) <- cell_ids
  # Keep colData minimal — only cell_ID + sample_id. Per-sample Giotto
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
    assays      = list(counts = counts),
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
    cat("⚠ Banksy not installed — skipping BANKSY niche computation.\n")
    return(NULL)
  }

  spes <- list()
  for (nm in names(per_sample_gobjs)) {
    cat("  BANKSY: per-sample augmentation for ", nm, "\n", sep = "")
    spe <- tryCatch(
      .banksy_giotto_to_spe(per_sample_gobjs[[nm]], nm),
      error = function(e) {
        cat("    ⚠ SPE conversion failed: ", conditionMessage(e), "\n",
            sep = "")
        NULL
      }
    )
    if (is.null(spe)) next
    aug <- tryCatch(
      Banksy::computeBanksy(spe,
                            assay_name = "counts",
                            k_geom = k_geom,
                            compute_agf = TRUE),
      error = function(e) {
        cat("    ⚠ computeBanksy failed: ", conditionMessage(e), "\n",
            sep = "")
        NULL
      }
    )
    if (!is.null(aug)) spes[[nm]] <- aug
  }
  if (length(spes) == 0L) return(NULL)
  spes
}

#' Concatenate augmented SPEs, run PCA → Harmony → Leiden, write
#' `niche_id` back to the merged Giotto cell metadata.
#'
#' @param merged_gobj Merged Giotto object (post batch_correct_merged_object).
#' @param banksy_aug_list Output of compute_banksy_per_sample().
#' @param batch_col Column for Harmony correction (default "sample_id").
#' @param resolution Leiden resolution.
#' @param output_dir Merged-pipeline output dir; niche_summary.csv is written
#'   to <output_dir>/14_BANKSY/.
cluster_banksy_niches <- function(merged_gobj,
                                  banksy_aug_list,
                                  batch_col = "sample_id",
                                  resolution = 0.4,
                                  lambda = 0.2,
                                  output_dir,
                                  n_pcs = 30,
                                  k_nn = 30) {
  if (is.null(banksy_aug_list)) {
    cat("  BANKSY: no augmented data — skipping niche clustering.\n")
    return(invisible(NULL))
  }
  if (!requireNamespace("Banksy", quietly = TRUE)) return(invisible(NULL))

  # cbind the per-sample augmented SPEs, keep only cells present in the
  # merged Giotto (defensive — joinGiottoObjects may have dropped some).
  aug <- do.call(SummarizedExperiment::cbind, unname(banksy_aug_list))
  merged_meta <- as.data.frame(.giotto_pdata_dt(merged_gobj))
  merged_ids <- as.character(merged_meta$cell_ID)
  keep <- colnames(aug) %in% merged_ids
  aug <- aug[, keep]
  # Re-order to match merged Giotto cell order.
  aug <- aug[, match(merged_ids, colnames(aug))]
  aug <- aug[, !is.na(colnames(aug))]
  if (ncol(aug) == 0L) {
    cat("  BANKSY: no overlapping cells with merged object — skipping.\n")
    return(invisible(NULL))
  }

  # PCA on the BANKSY-augmented assay. Banksy::runBanksyPCA wraps this.
  pca <- tryCatch(
    Banksy::runBanksyPCA(aug, lambda = lambda, npcs = n_pcs),
    error = function(e) {
      cat("    ⚠ runBanksyPCA failed: ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  if (is.null(pca)) return(invisible(NULL))
  pca_mat <- SingleCellExperiment::reducedDim(
    pca, paste0("PCA_M0_lam", lambda))
  if (is.null(pca_mat)) {
    # Older Banksy versions name the slot differently.
    rd_names <- SingleCellExperiment::reducedDimNames(pca)
    pca_mat <- SingleCellExperiment::reducedDim(pca, rd_names[1])
  }

  # Harmony on sample_id. Use the harmony R package directly.
  if (!requireNamespace("harmony", quietly = TRUE)) {
    cat("    ⚠ harmony not installed — skipping batch correction step.\n")
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
        cat("    ⚠ Harmony failed: ", conditionMessage(e),
            " — using uncorrected PCs.\n", sep = "")
        pca_mat
      }
    )
  }

  # Leiden on a kNN graph in Harmony space.
  if (!requireNamespace("igraph", quietly = TRUE)) {
    cat("    ⚠ igraph not installed — cannot cluster niches.\n")
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
    igraph::cluster_leiden(g, resolution_parameter = resolution),
    error = function(e) {
      cat("    ⚠ Leiden failed (", conditionMessage(e),
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
      cat("    ⚠ addCellMetadata failed for niche_id: ",
          conditionMessage(e), "\n", sep = "")
      merged_gobj
    }
  )

  # niche_summary.csv: cells per niche per sample.
  out_dir <- file.path(output_dir, "14_BANKSY")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  summary_df <- as.data.frame(table(
    sample_id = merged_meta[match(niche_df$cell_ID, merged_ids), batch_col],
    niche_id  = niche_df$niche_id
  ))
  utils::write.csv(summary_df,
                   file.path(out_dir, "niche_summary.csv"),
                   row.names = FALSE)
  cat("  ✓ BANKSY niche labels written; ", length(unique(niche_id)),
      " niches across ", length(unique(merged_meta[[batch_col]])),
      " samples.\n", sep = "")

  invisible(merged_gobj)
}

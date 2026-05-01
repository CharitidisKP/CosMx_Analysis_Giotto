#!/usr/bin/env Rscript
# ==============================================================================
# Cluster_Marker_Review.R
# Marker-review CSV per Leiden cluster on the merged Giotto object.
#
# Stage 1, Phase 8: after the global merged-mode annotation runs (annotate_cells
# on the merged object), this helper writes a quick top-N marker table per
# Leiden cluster so mis-labeled rare clusters can be recovered manually.
#
# Output: <output_dir>/10_Merged/global_annotation/marker_review.csv with
# columns: cluster | rank | gene | summary.logFC | summary.AUC | n_cells |
#          top_celltype_match | top_celltype_fraction
#
# top_celltype_match is a plurality vote of `celltype_<profile>_supervised`
# values within the cluster (taken from the per-sample annotations carried
# in by joinGiottoObjects). top_celltype_fraction is the share of cells in
# that cluster whose label matches the plurality.
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

.cmr_get_norm_expression <- function(gobj) {
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

.cmr_resolve_celltype_col <- function(meta, explicit = NULL) {
  if (!is.null(explicit) && nzchar(explicit) && explicit %in% names(meta))
    return(explicit)
  cand <- grep("^celltype_.*_supervised$", names(meta), value = TRUE)
  if (length(cand) > 0L) return(cand[1])
  if ("celltype" %in% names(meta)) return("celltype")
  NA_character_
}

#' Write a per-Leiden-cluster marker review CSV.
#'
#' @param gobj Merged Giotto object (after merged_annotate).
#' @param cluster_col Cluster column to summarize (default "leiden_clust").
#' @param celltype_col Per-sample cell-type column whose plurality vote is
#'   reported as `top_celltype_match` (default: auto-detect first
#'   "celltype_*_supervised" column).
#' @param output_csv Path to write the marker_review.csv.
#' @param top_n Top-N markers per cluster (default 25).
#' @return Invisibly, the marker_review data.frame.
write_cluster_marker_review <- function(gobj,
                                        cluster_col = "leiden_clust",
                                        celltype_col = NULL,
                                        output_csv,
                                        top_n = 25) {
  if (!requireNamespace("scran", quietly = TRUE)) {
    cat("⚠ scran not installed — skipping marker review CSV.\n")
    return(invisible(NULL))
  }
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (!cluster_col %in% names(meta)) {
    cat("⚠ Cluster column '", cluster_col,
        "' not found — skipping marker review.\n", sep = "")
    return(invisible(NULL))
  }
  ctc <- .cmr_resolve_celltype_col(meta, celltype_col)

  norm_mat <- .cmr_get_norm_expression(gobj)
  if (is.null(norm_mat)) {
    cat("⚠ No normalized expression available — skipping marker review.\n")
    return(invisible(NULL))
  }
  # Align columns; scran::findMarkers wants cells × clusters.
  cell_ids <- as.character(meta$cell_ID)
  norm_mat <- norm_mat[, cell_ids, drop = FALSE]
  cluster_vec <- as.factor(meta[[cluster_col]])

  cat("  Running scran::findMarkers on ", length(unique(cluster_vec)),
      " clusters × ", nrow(norm_mat), " genes...\n", sep = "")
  markers <- tryCatch(
    scran::findMarkers(
      x = norm_mat,
      groups = cluster_vec,
      direction = "up",
      pval.type = "any"
    ),
    error = function(e) {
      cat("⚠ scran::findMarkers failed: ", conditionMessage(e), "\n", sep = "")
      NULL
    }
  )
  if (is.null(markers)) return(invisible(NULL))

  # Per-cluster celltype plurality.
  if (!is.na(ctc)) {
    plurality <- vapply(levels(cluster_vec), function(cl) {
      sub_ct <- as.character(meta[[ctc]][cluster_vec == cl])
      sub_ct <- sub_ct[!is.na(sub_ct) & nzchar(sub_ct)]
      if (length(sub_ct) == 0L) return(NA_character_)
      tab <- sort(table(sub_ct), decreasing = TRUE)
      names(tab)[1]
    }, character(1))
    plurality_frac <- vapply(levels(cluster_vec), function(cl) {
      sub_ct <- as.character(meta[[ctc]][cluster_vec == cl])
      sub_ct <- sub_ct[!is.na(sub_ct) & nzchar(sub_ct)]
      if (length(sub_ct) == 0L) return(NA_real_)
      tab <- sort(table(sub_ct), decreasing = TRUE)
      as.numeric(tab[1] / sum(tab))
    }, numeric(1))
  } else {
    plurality <- setNames(rep(NA_character_, nlevels(cluster_vec)),
                          levels(cluster_vec))
    plurality_frac <- setNames(rep(NA_real_, nlevels(cluster_vec)),
                               levels(cluster_vec))
  }

  rows <- list()
  for (cl in names(markers)) {
    df <- as.data.frame(markers[[cl]])
    df <- df[, intersect(c("p.value", "FDR", "summary.logFC", "summary.AUC"),
                         names(df)), drop = FALSE]
    df$gene <- rownames(df)
    # findMarkers already orders by p.value; take top_n.
    df <- head(df, top_n)
    df$rank <- seq_len(nrow(df))
    df$cluster <- cl
    df$n_cells <- sum(cluster_vec == cl)
    df$top_celltype_match <- plurality[cl]
    df$top_celltype_fraction <- plurality_frac[cl]
    rows[[cl]] <- df
  }
  out_df <- do.call(rbind, rows)
  cols_order <- c("cluster", "rank", "gene",
                  intersect(c("summary.logFC", "summary.AUC", "p.value", "FDR"),
                            names(out_df)),
                  "n_cells", "top_celltype_match", "top_celltype_fraction")
  out_df <- out_df[, cols_order, drop = FALSE]

  dir.create(dirname(output_csv), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out_df, output_csv, row.names = FALSE)
  cat("  ✓ Marker review written: ", output_csv, " (",
      nrow(out_df), " rows, ", length(unique(out_df$cluster)),
      " clusters).\n", sep = "")
  invisible(out_df)
}

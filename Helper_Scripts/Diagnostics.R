#!/usr/bin/env Rscript
# ==============================================================================
# Diagnostics.R
# Post-Harmony batch-correction diagnostics for the merged Giotto object.
#
# Non-gating: writes per-metric CSVs and a one-page summary.md alongside the
# merged outputs. Threshold breaches are flagged in summary.md and warnings
# are echoed to stdout, but nothing here ever stops the pipeline.
#
# Required (optional) packages, wrapped in requireNamespace() so the
# pipeline still runs before the Stage-1 deps are installed:
#   - kBET   (theislab/kBET)       , batch-mixing acceptance rate
#   - lisi   (immunogenomics/LISI) , iLISI / cLISI scores
#
# Always available:
#   - Delaunay-edge audit (no external deps; uses Giotto accessors)
# ==============================================================================

current_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]),
                                 winslash = "/", mustWork = FALSE)))
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

# Reuse Giotto accessor shims from Merge_Batch_Correction.R (loaded earlier
# in the pipeline). If sourced standalone, define a minimal pDataDT.
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

# Pulls the Harmony-corrected embedding matrix (cells × dims) from the
# merged Giotto object. Returns NULL with a clear message if the slot is
# missing (so diagnostics fall back to the raw PCA embedding).
.diag_get_embedding <- function(gobj,
                                reduction_name = "harmony",
                                reduction_method = "harmony") {
  acc <- if (requireNamespace("Giotto", quietly = TRUE) &&
             exists("get_dimReduction", envir = asNamespace("Giotto"), inherits = FALSE)) {
    get("get_dimReduction", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("get_dimReduction", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    get("get_dimReduction", envir = asNamespace("GiottoClass"))
  } else {
    get("get_dimReduction", mode = "function")
  }
  emb <- tryCatch(
    acc(gobject = gobj, reduction = "cells",
        reduction_method = reduction_method,
        name = reduction_name,
        output = "matrix"),
    error = function(e) NULL
  )
  if (is.null(emb)) {
    # Fallback to PCA if Harmony slot is absent.
    emb <- tryCatch(
      acc(gobject = gobj, reduction = "cells",
          reduction_method = "pca",
          name = "pca",
          output = "matrix"),
      error = function(e) NULL
    )
    if (!is.null(emb)) {
      cat("  Diagnostics: Harmony embedding not found; falling back to PCA.\n")
    }
  }
  emb
}

# Pulls the Delaunay spatial-network edge table.
.diag_get_delaunay_edges <- function(gobj,
                                     network_name = "Delaunay_network") {
  acc <- if (requireNamespace("Giotto", quietly = TRUE) &&
             exists("getSpatialNetwork", envir = asNamespace("Giotto"), inherits = FALSE)) {
    get("getSpatialNetwork", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getSpatialNetwork", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    get("getSpatialNetwork", envir = asNamespace("GiottoClass"))
  } else {
    return(NULL)
  }
  tryCatch(
    acc(gobject = gobj, name = network_name, output = "data.table"),
    error = function(e) NULL
  )
}

#' Run post-Harmony diagnostics on a merged Giotto object.
#'
#' Writes kbet.csv, lisi.csv, edge_audit.txt, and summary.md to
#' `<output_dir>/10_Merged/diagnostics/`. Non-gating; warns on threshold
#' breaches but never stops the pipeline.
#'
#' @param gobj Merged Giotto object after batch_correct_merged_object().
#' @param batch_col Metadata column treated as the batch identifier
#'   (default "sample_id"; matches cfg$merged$batch_column).
#' @param biology_cols Character vector of biological metadata columns
#'   that should NOT be mixed by Harmony, cLISI on these reports whether
#'   biological signal survived (default c("treatment", "timepoint")).
#' @param output_dir Merged-pipeline output directory (the parent of the
#'   "10_Merged" folder).
#' @param thresholds Named list of pass/fail thresholds for summary.md.
#' @return Invisible list with raw scores + summary text.
run_merge_diagnostics <- function(gobj,
                                  batch_col = "sample_id",
                                  biology_cols = c("treatment", "timepoint"),
                                  output_dir,
                                  thresholds = list(
                                    kbet_acc_min = 0.5,
                                    ilisi_batch_min = 1.5,
                                    clisi_biology_min = 1.2,
                                    cross_sample_edges_max = 0L
                                  )) {
  diag_dir <- file.path(output_dir, "10_Merged", "diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (!batch_col %in% names(meta)) {
    msg <- sprintf("Diagnostics: batch column '%s' not in cell metadata; skipping.",
                   batch_col)
    warning(msg, immediate. = TRUE, call. = FALSE)
    writeLines(msg, file.path(diag_dir, "summary.md"))
    return(invisible(NULL))
  }
  emb <- .diag_get_embedding(gobj)

  kbet_result <- .run_kbet(emb, meta[[batch_col]], diag_dir)
  lisi_result <- .run_lisi(emb, meta, batch_col, biology_cols, diag_dir)
  edge_result <- .run_edge_audit(gobj, meta, batch_col, diag_dir)

  .write_summary_md(diag_dir, batch_col, biology_cols, thresholds,
                    kbet_result, lisi_result, edge_result)

  invisible(list(kbet = kbet_result, lisi = lisi_result, edges = edge_result))
}

# ------------------------------------------------------------------------------
# kBET, batch-mixing acceptance rate.
# ------------------------------------------------------------------------------
.run_kbet <- function(emb, batch, diag_dir) {
  if (is.null(emb)) {
    return(list(status = "skipped", reason = "no embedding available",
                acc_rate = NA_real_))
  }
  if (!requireNamespace("kBET", quietly = TRUE)) {
    return(list(status = "skipped",
                reason = "package 'kBET' not installed",
                acc_rate = NA_real_))
  }
  res <- tryCatch(
    kBET::kBET(df = emb, batch = batch, plot = FALSE, do.pca = FALSE,
               verbose = FALSE),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(list(status = "error", reason = conditionMessage(res),
                acc_rate = NA_real_))
  }
  acc_rate <- 1 - mean(res$summary$kBET.observed, na.rm = TRUE)
  out <- data.frame(
    metric = c("kBET_acceptance_rate", "kBET_observed_mean",
               "kBET_expected_mean", "n_cells"),
    value = c(acc_rate,
              mean(res$summary$kBET.observed, na.rm = TRUE),
              mean(res$summary$kBET.expected, na.rm = TRUE),
              length(batch)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(out, file.path(diag_dir, "kbet.csv"), row.names = FALSE)
  list(status = "ok", acc_rate = acc_rate, raw = res$summary)
}

# ------------------------------------------------------------------------------
# LISI, iLISI on batch (mixing), cLISI on biology cols (preservation).
# ------------------------------------------------------------------------------
.run_lisi <- function(emb, meta, batch_col, biology_cols, diag_dir) {
  if (is.null(emb)) {
    return(list(status = "skipped", reason = "no embedding available",
                scores = data.frame()))
  }
  if (!requireNamespace("lisi", quietly = TRUE)) {
    return(list(status = "skipped",
                reason = "package 'lisi' not installed",
                scores = data.frame()))
  }
  cols <- c(batch_col, biology_cols)
  cols <- cols[cols %in% names(meta)]
  if (length(cols) == 0L) {
    return(list(status = "skipped",
                reason = "no requested columns present in metadata",
                scores = data.frame()))
  }
  scores <- tryCatch(
    lisi::compute_lisi(X = emb, meta_data = meta, label_colnames = cols),
    error = function(e) e
  )
  if (inherits(scores, "error")) {
    return(list(status = "error", reason = conditionMessage(scores),
                scores = data.frame()))
  }
  summary_df <- data.frame(
    column = cols,
    role = ifelse(cols == batch_col, "batch_iLISI", "biology_cLISI"),
    mean_lisi = vapply(cols, function(co) mean(scores[[co]], na.rm = TRUE), numeric(1)),
    median_lisi = vapply(cols, function(co) stats::median(scores[[co]], na.rm = TRUE), numeric(1)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary_df, file.path(diag_dir, "lisi.csv"), row.names = FALSE)
  list(status = "ok", scores = summary_df)
}

# ------------------------------------------------------------------------------
# Cross-sample-edge audit on the Delaunay spatial network.
# ------------------------------------------------------------------------------
.run_edge_audit <- function(gobj, meta, batch_col, diag_dir,
                            network_name = "Delaunay_network") {
  edges <- .diag_get_delaunay_edges(gobj, network_name)
  out_path <- file.path(diag_dir, "edge_audit.txt")
  if (is.null(edges) || nrow(edges) == 0) {
    writeLines(c(
      sprintf("Network '%s' not found on merged object, audit skipped.",
              network_name),
      "If the network was rebuilt post-merge, the cross-sample edge count",
      "should be re-checked manually."
    ), out_path)
    return(list(status = "skipped",
                reason = sprintf("network '%s' missing", network_name),
                cross_sample_edges = NA_integer_))
  }
  # Edge endpoints stored as `from` / `to` with cell_ID values; some
  # Giotto versions name them `cell_ID_1` / `cell_ID_2`.
  from_col <- intersect(c("from", "cell_ID_1"), names(edges))[1]
  to_col   <- intersect(c("to",   "cell_ID_2"), names(edges))[1]
  if (is.na(from_col) || is.na(to_col)) {
    writeLines(c(
      sprintf("Network '%s' present but has no recognisable from/to columns.",
              network_name),
      sprintf("Available columns: %s", paste(names(edges), collapse = ", "))
    ), out_path)
    return(list(status = "skipped",
                reason = "unrecognised edge column layout",
                cross_sample_edges = NA_integer_))
  }
  meta_lookup <- setNames(as.character(meta[[batch_col]]),
                          as.character(meta$cell_ID))
  from_batch <- meta_lookup[as.character(edges[[from_col]])]
  to_batch   <- meta_lookup[as.character(edges[[to_col]])]
  cross <- !is.na(from_batch) & !is.na(to_batch) & from_batch != to_batch
  n_cross <- sum(cross)
  n_total <- length(cross)
  writeLines(c(
    sprintf("Network: %s", network_name),
    sprintf("Total edges: %d", n_total),
    sprintf("Cross-%s edges: %d (%.4f%%)",
            batch_col, n_cross,
            100 * n_cross / max(n_total, 1L)),
    "(Expected: 0, Delaunay network is built per-sample before joinGiottoObjects.)"
  ), out_path)
  list(status = "ok", cross_sample_edges = n_cross, total_edges = n_total)
}

# ------------------------------------------------------------------------------
# Synthesise summary.md.
# ------------------------------------------------------------------------------
.write_summary_md <- function(diag_dir, batch_col, biology_cols, thresh,
                              kbet, lisi, edges) {
  flag_pass <- function(ok) if (isTRUE(ok)) "PASS" else if (is.na(ok)) "SKIP" else "FAIL"
  lines <- c("# Merge Diagnostics Summary",
             "",
             sprintf("Batch column: `%s`", batch_col),
             sprintf("Biology columns checked for cLISI: %s",
                     paste(sprintf("`%s`", biology_cols), collapse = ", ")),
             "",
             "| Metric | Value | Threshold | Status |",
             "|---|---|---|---|")

  # kBET
  kbet_ok <- if (!identical(kbet$status, "ok")) NA
              else kbet$acc_rate >= thresh$kbet_acc_min
  lines <- c(lines, sprintf("| kBET acceptance rate | %s | ≥ %.2f | %s |",
                            if (is.na(kbet$acc_rate)) "," else sprintf("%.3f", kbet$acc_rate),
                            thresh$kbet_acc_min, flag_pass(kbet_ok)))
  if (!identical(kbet$status, "ok")) {
    lines <- c(lines, sprintf("> kBET %s: %s", kbet$status, kbet$reason %||% ""))
  }

  # iLISI / cLISI
  if (identical(lisi$status, "ok") && nrow(lisi$scores) > 0) {
    for (i in seq_len(nrow(lisi$scores))) {
      r <- lisi$scores[i, ]
      thr <- if (r$role == "batch_iLISI") thresh$ilisi_batch_min
             else thresh$clisi_biology_min
      ok <- r$mean_lisi >= thr
      lines <- c(lines, sprintf("| %s(`%s`) mean | %.3f | ≥ %.2f | %s |",
                                r$role, r$column, r$mean_lisi, thr, flag_pass(ok)))
    }
  } else {
    lines <- c(lines, "| LISI |, |, | SKIP |",
                       sprintf("> LISI %s: %s", lisi$status, lisi$reason %||% ""))
  }

  # Edge audit
  edge_ok <- if (!identical(edges$status, "ok")) NA
              else edges$cross_sample_edges <= thresh$cross_sample_edges_max
  lines <- c(lines,
             sprintf("| Cross-`%s` Delaunay edges | %s | ≤ %d | %s |",
                     batch_col,
                     if (is.na(edges$cross_sample_edges)) ","
                       else as.character(edges$cross_sample_edges),
                     thresh$cross_sample_edges_max,
                     flag_pass(edge_ok)))
  if (!identical(edges$status, "ok")) {
    lines <- c(lines, sprintf("> Edge audit %s: %s",
                              edges$status, edges$reason %||% ""))
  }

  lines <- c(lines, "",
             "_Non-gating: any failure here is logged and surfaced in the merged-run summary report; the pipeline continues._")

  writeLines(lines, file.path(diag_dir, "summary.md"))

  # Echo failures to stdout so they're visible in the run log.
  if (isTRUE(!is.na(kbet_ok) && !kbet_ok))
    cat(sprintf("  [diagnostics] kBET acceptance %.3f below threshold %.2f\n",
                kbet$acc_rate, thresh$kbet_acc_min))
  if (identical(lisi$status, "ok")) {
    bad <- lisi$scores[lisi$scores$role == "batch_iLISI" &
                       lisi$scores$mean_lisi < thresh$ilisi_batch_min, , drop = FALSE]
    for (i in seq_len(nrow(bad)))
      cat(sprintf("  [diagnostics] iLISI(%s) mean %.3f below threshold %.2f\n",
                  bad$column[i], bad$mean_lisi[i], thresh$ilisi_batch_min))
    bad <- lisi$scores[lisi$scores$role == "biology_cLISI" &
                       lisi$scores$mean_lisi < thresh$clisi_biology_min, , drop = FALSE]
    for (i in seq_len(nrow(bad)))
      cat(sprintf("  [diagnostics] cLISI(%s) mean %.3f below threshold %.2f (biology may be over-corrected)\n",
                  bad$column[i], bad$mean_lisi[i], thresh$clisi_biology_min))
  }
  if (isTRUE(!is.na(edge_ok) && !edge_ok))
    cat(sprintf("  [diagnostics] %d cross-%s edges in Delaunay_network (expected 0)\n",
                edges$cross_sample_edges, batch_col))
}

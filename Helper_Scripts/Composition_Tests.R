#!/usr/bin/env Rscript
# ==============================================================================
# Composition_Tests.R
# Cell-type composition differential abundance via speckle::propeller.
#
# For each entry in `cfg$pathway$comparisons` with design ∈ {paired,
# unpaired}, builds a per-sample cell-type fraction matrix and runs
# propeller paired by `block` (default patient_id). Skips descriptive
# comparisons (no inferential test makes sense at n=1).
#
# Output: <output_dir>/15_Composition/<comparison_label>.csv with columns
#   cell_type | logFC | tstat | p | padj | n_per_side | note
# `note = "low_power"` is set whenever either side has < 3 samples.
#
# `niche_id` differential abundance (post Phase 7) is computed by the
# parallel helper run_niche_composition_tests() which shares 90% of the
# implementation.
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

# Resolves a comparison's group_a / group_b cell-id sets and the per-sample
# metadata frame. Returns NULL when the contrast can't be filled in (e.g.
# group_b empty, n=1 per side and not enough block levels).
.composition_resolve_groups <- function(meta, comparison_spec, sample_col = "sample_id") {
  ga <- comparison_spec$group_a %||% list()
  gb <- comparison_spec$group_b %||% list()
  match_rows <- function(spec) {
    if (length(spec) == 0L) return(rep(FALSE, nrow(meta)))
    keep <- rep(TRUE, nrow(meta))
    for (col in names(spec)) {
      if (!col %in% names(meta)) return(rep(FALSE, nrow(meta)))
      keep <- keep & as.character(meta[[col]]) == as.character(spec[[col]])
    }
    keep
  }
  rows_a <- match_rows(ga)
  rows_b <- match_rows(gb)
  if (sum(rows_a) == 0L || sum(rows_b) == 0L) return(NULL)
  list(
    samples_a = unique(as.character(meta[[sample_col]][rows_a])),
    samples_b = unique(as.character(meta[[sample_col]][rows_b]))
  )
}

# Builds a per-sample cell-type fraction matrix (rows = cell types, cols
# = samples). Counts are computed from `meta`'s cell-type column; samples
# with zero cells in a given cell type get a zero count (no NaN).
.composition_build_fractions <- function(meta, sample_col, celltype_col, samples) {
  sub <- meta[as.character(meta[[sample_col]]) %in% samples, , drop = FALSE]
  if (nrow(sub) == 0L) return(NULL)
  counts <- table(
    factor(sub[[celltype_col]]),
    factor(sub[[sample_col]], levels = samples)
  )
  fractions <- sweep(counts, 2, pmax(colSums(counts), 1L), "/")
  list(
    counts = unclass(counts),
    fractions = unclass(fractions),
    samples = samples
  )
}

# Tags low-power rows in the propeller result.
.composition_tag_power <- function(df, n_a, n_b, threshold = 3L) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  df$n_per_side <- sprintf("%d / %d", n_a, n_b)
  df$note <- ifelse(min(n_a, n_b) < threshold, "low_power", "")
  df
}

#' Run cell-type composition tests via propeller for every paired/unpaired
#' comparison in `cfg$pathway$comparisons`.
#'
#' @param gobj Merged Giotto object after batch_correct_merged_object().
#' @param comparisons List read from `cfg$pathway$comparisons` (each element
#'   has `label`, `group_a`, `group_b`, `design`, optional `block`).
#' @param output_dir Merged-pipeline output directory (the parent of
#'   "15_Composition").
#' @param cfg Full config list (used for column names and defaults).
#' @return Invisible list, one element per comparison, with the propeller
#'   table or NULL if the contrast was skipped.
run_composition_tests <- function(gobj,
                                  comparisons,
                                  output_dir,
                                  cfg) {
  if (!requireNamespace("speckle", quietly = TRUE)) {
    cat("⚠ speckle not installed — skipping composition tests.\n")
    return(invisible(NULL))
  }

  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  sample_col <- cfg$spatial_de$sample_column %||% "sample_id"
  celltype_col <- cfg$pathway$celltype_column %||%
                  cfg$interaction$annotation_column %||%
                  "celltype_HCA_Kidney_supervised"
  if (!celltype_col %in% names(meta)) {
    # Fall back to any celltype_*_supervised column present.
    cand <- grep("^celltype_.*_supervised$", names(meta), value = TRUE)
    if (length(cand) == 0L) {
      cat("⚠ Composition tests: no cell-type column found.\n")
      return(invisible(NULL))
    }
    celltype_col <- cand[1]
    cat("  Composition: falling back to cell-type column '",
        celltype_col, "'\n", sep = "")
  }

  out_dir <- file.path(output_dir, "15_Composition")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n========================================\n")
  cat("Composition tests (propeller, paired by patient)\n")
  cat("Cell-type column: ", celltype_col, "\n", sep = "")
  cat("========================================\n\n")

  results <- list()
  for (cp in comparisons) {
    label <- cp$label %||% "<unnamed>"
    design_kind <- cp$design %||% "unpaired"
    if (design_kind == "descriptive") {
      cat("  [", label, "] design=descriptive — skipped (no inferential test).\n",
          sep = "")
      next
    }
    block_col <- cp$block %||% "patient_id"

    grp <- .composition_resolve_groups(meta, cp, sample_col = sample_col)
    if (is.null(grp)) {
      cat("  [", label, "] group_a or group_b empty — skipped.\n", sep = "")
      next
    }
    samples <- c(grp$samples_a, grp$samples_b)
    if (length(samples) < 2L) {
      cat("  [", label, "] only ", length(samples), " sample(s) — skipped.\n",
          sep = "")
      next
    }

    fractions <- .composition_build_fractions(
      meta, sample_col, celltype_col, samples
    )
    if (is.null(fractions)) next

    counts_mat <- fractions$counts          # rows = celltypes, cols = samples
    samples_ordered <- colnames(counts_mat)

    # Per-sample group + (optionally) block lookups.
    sample_meta <- unique(meta[meta[[sample_col]] %in% samples_ordered,
                               c(sample_col, block_col)])
    sample_meta <- sample_meta[match(samples_ordered, sample_meta[[sample_col]]),
                                , drop = FALSE]
    group_vec <- ifelse(samples_ordered %in% grp$samples_a, "A", "B")

    n_a <- sum(group_vec == "A")
    n_b <- sum(group_vec == "B")

    propeller_design <- if (design_kind == "paired" &&
                            block_col %in% names(sample_meta) &&
                            length(unique(sample_meta[[block_col]])) >= 2L) "paired"
                        else "unpaired"

    res <- tryCatch({
      if (propeller_design == "paired") {
        # propeller.ttest with blocking on patient_id.
        block <- factor(sample_meta[[block_col]])
        speckle::propeller.ttest(
          prop.list = speckle::getTransformedProps(
            counts = counts_mat, transform = "asin"),
          design = stats::model.matrix(~ block + group_vec),
          coef = "group_vecB",
          robust = TRUE, trend = FALSE, sort = FALSE
        )
      } else {
        speckle::propeller.ttest(
          prop.list = speckle::getTransformedProps(
            counts = counts_mat, transform = "asin"),
          design = stats::model.matrix(~ group_vec),
          coef = "group_vecB",
          robust = TRUE, trend = FALSE, sort = FALSE
        )
      }
    }, error = function(e) {
      cat("  [", label, "] propeller failed: ", conditionMessage(e), "\n",
          sep = "")
      NULL
    })

    if (is.null(res)) next

    out_df <- data.frame(
      cell_type   = rownames(res),
      logFC       = res$PropMean.B - res$PropMean.A,
      tstat       = res$Tstatistic,
      p           = res$P.Value,
      padj        = res$FDR,
      n_per_side  = sprintf("%d / %d", n_a, n_b),
      note        = ifelse(min(n_a, n_b) < 3L, "low_power", ""),
      design_used = propeller_design,
      stringsAsFactors = FALSE
    )
    rownames(out_df) <- NULL
    out_path <- file.path(out_dir, paste0(label, ".csv"))
    utils::write.csv(out_df, out_path, row.names = FALSE)
    cat("  ✓ [", label, "] design=", propeller_design,
        " n=", n_a, "/", n_b,
        " — wrote ", out_path, "\n", sep = "")
    results[[label]] <- out_df
  }

  invisible(results)
}

#' Niche-id composition tests. Identical machinery to the cell-type version
#' but stratified by `niche_id` (post BANKSY). No-op when `niche_id` is
#' absent from cell metadata (Phase 7 hasn't run yet).
run_niche_composition_tests <- function(gobj,
                                        comparisons,
                                        output_dir,
                                        cfg,
                                        niche_col = "niche_id") {
  if (!requireNamespace("speckle", quietly = TRUE)) {
    cat("⚠ speckle not installed — skipping niche composition tests.\n")
    return(invisible(NULL))
  }
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (!niche_col %in% names(meta)) {
    cat("  Niche composition: '", niche_col,
        "' column missing — run BANKSY (Phase 7) first.\n", sep = "")
    return(invisible(NULL))
  }
  # Re-use run_composition_tests by renaming niche_col → celltype_<niche>
  # in a shadow object so the existing path picks it up. To avoid mutating
  # the merged Giotto object, we re-implement the small inner loop here
  # against a tweaked celltype_col.
  cfg2 <- cfg
  cfg2$pathway$celltype_column <- niche_col
  out_dir <- file.path(output_dir, "15_Composition", "by_niche")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # The existing runner expects 15_Composition/<label>.csv; redirect by
  # passing a child output_dir that resolves to 15_Composition/by_niche.
  results <- run_composition_tests(
    gobj         = gobj,
    comparisons  = comparisons,
    output_dir   = out_dir,        # → out_dir/15_Composition/<label>.csv
    cfg          = cfg2
  )
  invisible(results)
}

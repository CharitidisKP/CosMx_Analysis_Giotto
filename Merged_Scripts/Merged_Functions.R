# =============================================================================
# Merged_Functions.R
#
# Reusable helpers for Merged_Scripts/Paired_Merged_Analysis.Rmd.
# Populated incrementally — only multi-chunk helpers live here.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Pull cfg$merged$wtx$<key>, falling back to `default`.
.wtx <- function(cfg, key, default) {
  v <- cfg$merged$wtx[[key]]
  if (is.null(v) || (is.character(v) && !nzchar(v))) default else v
}

# Raw counts (genes x cells) regardless of accessor namespace.
.gobj_raw_counts <- function(gobj) {
  if (requireNamespace("Giotto", quietly = TRUE)) {
    Giotto::getExpression(gobj, values = "raw", output = "matrix")
  } else {
    GiottoClass::getExpression(gobj, values = "raw", output = "matrix")
  }
}

# Push a per-cell metadata data.frame back onto a Giotto object.
.gobj_add_meta <- function(gobj, df, by = "cell_ID") {
  acc <- if (requireNamespace("Giotto", quietly = TRUE) &&
             exists("addCellMetadata", asNamespace("Giotto"), inherits = FALSE)) {
    Giotto::addCellMetadata
  } else {
    GiottoClass::addCellMetadata
  }
  acc(gobj, new_metadata = df, by_column = TRUE, column_cell_ID = by)
}

# Mirrors 03_Normalisation.R:40 — Merge_Batch_Correction.R needs this in scope.
if (!exists(".run_known_giotto_warning_safe", inherits = FALSE)) {
  .run_known_giotto_warning_safe <- function(expr) {
    withCallingHandlers(
      expr,
      warning = function(w) {
        if (grepl("Not all expression matrices share the same cell_IDs",
                  conditionMessage(w), fixed = TRUE)) {
          invokeRestart("muffleWarning")
        }
      }
    )
  }
}

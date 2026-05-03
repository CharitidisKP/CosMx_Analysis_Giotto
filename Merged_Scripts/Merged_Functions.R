# Merged_Functions.R — helpers for Paired_Merged_Analysis.Rmd.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

# Resolve cfg$merged$wtx$<key>, defaulting to `default`.
.wtx <- function(cfg, key, default) {
  v <- cfg$merged$wtx[[key]]
  if (is.null(v) || (is.character(v) && !nzchar(v))) default else v
}

# Null out cached polygon centroids so Giotto::joinGiottoObjects() can rebuild them.
.strip_polygon_centroids <- function(gobj) {
  tryCatch({
    poly_slot <- gobj@spatial_info
    if (!is.null(poly_slot) && length(poly_slot) > 0) {
      for (nm in names(poly_slot)) {
        poly <- poly_slot[[nm]]
        if (inherits(poly, "giottoPolygon")) {
          poly@spatVectorCentroids <- NULL
          gobj@spatial_info[[nm]] <- poly
        }
      }
    }
    gobj
  }, error = function(e) {
    message("  ⚠ Could not strip polygon centroids (non-fatal): ", conditionMessage(e))
    gobj
  })
}

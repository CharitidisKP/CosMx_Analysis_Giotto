# =============================================================================
# Merged_Functions.R
#
# Reusable helpers for the paired-merged analysis notebook
# (Merged_Scripts/Paired_Merged_Analysis.Rmd).
#
# Populated incrementally — only functions that are actually called from more
# than one chunk live here. One-off chunk logic stays inline in the Rmd.
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

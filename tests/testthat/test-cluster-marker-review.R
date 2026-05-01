# Phase 8 — Cluster_Marker_Review.R. Skipped if scran isn't installed
# (Stage 1 helpers wrap their package calls in requireNamespace fallbacks
# so the rest of the pipeline keeps working in that case).

test_that("write_cluster_marker_review returns NULL when scran missing", {
  skip_if(requireNamespace("scran", quietly = TRUE),
          "scran is installed; this test only verifies the skip-and-log fallback path")
  env <- .source_helper(file.path("Helper_Scripts", "Cluster_Marker_Review.R"))
  # Pass dummy gobj — function should bail early on requireNamespace check.
  res <- expect_invisible(env$write_cluster_marker_review(
    gobj = NULL, output_csv = tempfile()))
  expect_null(res)
})

test_that("Cluster_Marker_Review.R parses without sourcing external deps", {
  env <- .source_helper(file.path("Helper_Scripts", "Cluster_Marker_Review.R"))
  expect_true(exists("write_cluster_marker_review", envir = env, inherits = FALSE))
  expect_true(exists(".cmr_resolve_celltype_col", envir = env, inherits = FALSE))
})

# Phase 6 — Composition_Tests.R. Uses speckle when available; otherwise
# tests the group-resolution helper only.

test_that("Composition_Tests.R parses and exposes the public runner", {
  env <- .source_helper(file.path("Helper_Scripts", "Composition_Tests.R"))
  expect_true(exists("run_composition_tests", envir = env, inherits = FALSE))
  expect_true(exists("run_niche_composition_tests", envir = env, inherits = FALSE))
})

test_that(".composition_resolve_groups picks correct samples", {
  env <- .source_helper(file.path("Helper_Scripts", "Composition_Tests.R"))
  meta <- .synthetic_full_selection()
  meta <- meta[rep(seq_len(nrow(meta)), each = 5), , drop = FALSE]  # 5 cells/sample
  cp <- list(
    group_a = list(treatment = "CART", timepoint = "T0"),
    group_b = list(treatment = "CART", timepoint = "T12")
  )
  res <- env$.composition_resolve_groups(meta, cp, sample_col = "sample_id")
  expect_equal(sort(res$samples_a), sort(c("CART_T0_S1", "CART_T0_S2")))
  expect_equal(sort(res$samples_b), sort(c("CART_T12_S1", "CART_T12_S2")))
})

test_that("speckle-dependent path skips cleanly when speckle missing", {
  skip_if(requireNamespace("speckle", quietly = TRUE),
          "speckle installed; this test only verifies the skip path")
  env <- .source_helper(file.path("Helper_Scripts", "Composition_Tests.R"))
  res <- expect_invisible(env$run_composition_tests(
    gobj = NULL, comparisons = list(), output_dir = tempdir(), cfg = list()))
  expect_null(res)
})

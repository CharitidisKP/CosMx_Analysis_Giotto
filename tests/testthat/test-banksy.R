# Phase 7 — BANKSY_Niches.R. Most paths require Banksy + SpatialExperiment;
# here we test parsing + the missing-package skip path.

test_that("BANKSY_Niches.R parses and exposes the two public functions", {
  env <- .source_helper(file.path("Helper_Scripts", "BANKSY_Niches.R"))
  expect_true(exists("compute_banksy_per_sample", envir = env, inherits = FALSE))
  expect_true(exists("cluster_banksy_niches", envir = env, inherits = FALSE))
})

test_that("compute_banksy_per_sample returns NULL when Banksy missing", {
  skip_if(requireNamespace("Banksy", quietly = TRUE),
          "Banksy installed; this test only verifies the skip path")
  env <- .source_helper(file.path("Helper_Scripts", "BANKSY_Niches.R"))
  res <- env$compute_banksy_per_sample(per_sample_gobjs = list())
  expect_null(res)
})

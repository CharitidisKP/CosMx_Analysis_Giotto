# Phase 2 — build_smide_formula() helper. Pure R, no smiDE / Giotto deps.

test_that("multi-patient stratum injects + patient_id", {
  env <- .source_helper("12_Spatial_Differential_Expression.R")
  m <- data.frame(
    patient_id = rep(c("CART_pt1", "CART_pt2", "Conv_pt1", "Conv_pt2"),
                     each = 25)
  )
  fb <- env$build_smide_formula(
    base_terms = c("RankNorm(otherct_expr)", "smide_niche_treatment", "sample_id"),
    metadata = m
  )
  expect_equal(fb$paired_blocking, "patient_id")
  expect_match(fb$formula_str, "\\+ patient_id")
})

test_that("single-patient stratum skips block; no error", {
  env <- .source_helper("12_Spatial_Differential_Expression.R")
  m <- data.frame(patient_id = rep("CART_pt1", 100))
  fb <- env$build_smide_formula(
    base_terms = c("RankNorm(otherct_expr)", "treatment"),
    metadata = m,
    quiet = TRUE
  )
  expect_equal(fb$paired_blocking, "skipped_single_patient")
  expect_false(grepl("\\+ patient_id", fb$formula_str))
})

test_that("missing block column is reported, not fatal", {
  env <- .source_helper("12_Spatial_Differential_Expression.R")
  m <- data.frame(other = 1:10)
  fb <- env$build_smide_formula(
    base_terms = c("RankNorm(otherct_expr)", "treatment"),
    metadata = m,
    quiet = TRUE
  )
  expect_equal(fb$paired_blocking, "skipped_block_col_missing")
})

test_that("enable_blocking = FALSE returns 'disabled'", {
  env <- .source_helper("12_Spatial_Differential_Expression.R")
  m <- data.frame(patient_id = rep(c("a", "b"), 50))
  fb <- env$build_smide_formula(
    base_terms = c("RankNorm(otherct_expr)", "treatment"),
    metadata = m,
    enable_blocking = FALSE
  )
  expect_equal(fb$paired_blocking, "disabled")
})

test_that("formula always ends with offset(log(nCount_RNA))", {
  env <- .source_helper("12_Spatial_Differential_Expression.R")
  m <- data.frame(patient_id = rep(c("a", "b"), 50))
  fb <- env$build_smide_formula(
    base_terms = c("RankNorm(otherct_expr)", "treatment"),
    metadata = m
  )
  expect_match(fb$formula_str, "offset\\(log\\(nCount_RNA\\)\\)$")
})

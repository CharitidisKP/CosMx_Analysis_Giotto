# Phase 11 — verify step orders + alias resolution after Stage 1 wiring.

test_that("MERGED_STEP_ORDER includes all five new steps in correct order", {
  env <- new.env(parent = globalenv())
  options(cosmx.disable_cli = TRUE)
  source(file.path(.repo_dir(), "CosMx_pipeline.R"), local = env)
  expected <- c("merge", "merge_batch", "merged_annotate", "banksy",
                "composition", "12_spatial_de", "cci_paired", "13_pathway",
                "summary")
  expect_identical(env$MERGED_STEP_ORDER, expected)
})

test_that("SAMPLE_STEP_ORDER ends with 13_pathway", {
  env <- new.env(parent = globalenv())
  options(cosmx.disable_cli = TRUE)
  source(file.path(.repo_dir(), "CosMx_pipeline.R"), local = env)
  expect_equal(tail(env$SAMPLE_STEP_ORDER, 1), "13_pathway")
})

test_that("merged-side aliases resolve all new step IDs", {
  env <- new.env(parent = globalenv())
  options(cosmx.disable_cli = TRUE)
  source(file.path(.repo_dir(), "CosMx_pipeline.R"), local = env)

  expect_equal(env$canonical_step_ids("banksy", "merged"), "banksy")
  expect_equal(env$canonical_step_ids("14_banksy", "merged"), "banksy")
  expect_equal(env$canonical_step_ids("composition", "merged"), "composition")
  expect_equal(env$canonical_step_ids("propeller", "merged"), "composition")
  expect_equal(env$canonical_step_ids("merged_annotate", "merged"),
               "merged_annotate")
  expect_equal(env$canonical_step_ids("annotate_merged", "merged"),
               "merged_annotate")
  expect_equal(env$canonical_step_ids("cci_paired", "merged"), "cci_paired")
  expect_equal(env$canonical_step_ids("paired_cci", "merged"), "cci_paired")
  expect_equal(env$canonical_step_ids("summary", "merged"), "summary")
  expect_equal(env$canonical_step_ids("report", "merged"), "summary")
})

test_that("sample-side bare numeric '13' resolves to 13_pathway", {
  env <- new.env(parent = globalenv())
  options(cosmx.disable_cli = TRUE)
  source(file.path(.repo_dir(), "CosMx_pipeline.R"), local = env)
  expect_equal(env$canonical_step_ids("13", "sample"), "13_pathway")
})

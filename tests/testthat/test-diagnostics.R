# Phase 3 — Diagnostics.R. Most paths require kBET / lisi which aren't
# installed locally; here we test the parts that don't.

test_that("Diagnostics.R parses and exposes run_merge_diagnostics", {
  env <- .source_helper(file.path("Helper_Scripts", "Diagnostics.R"))
  expect_true(exists("run_merge_diagnostics", envir = env, inherits = FALSE))
})

test_that("kBET / LISI sub-runners gracefully skip when packages missing", {
  skip_if(requireNamespace("kBET", quietly = TRUE),
          "kBET installed; skip test only checks the missing-package path")
  env <- .source_helper(file.path("Helper_Scripts", "Diagnostics.R"))
  res <- env$.run_kbet(emb = matrix(rnorm(100), 10, 10),
                       batch = rep(c("a", "b"), 5),
                       diag_dir = tempdir())
  expect_equal(res$status, "skipped")
  expect_match(res$reason, "kBET")
})

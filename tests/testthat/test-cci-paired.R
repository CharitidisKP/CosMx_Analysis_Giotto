# Phase 10 — paired Wilcoxon on synthetic LR scores recovers a known
# signal across patients. Pure R — no LIANA / Giotto deps.

test_that("paired Wilcoxon recovers a known signal across patients", {
  env <- .source_helper(file.path("Helper_Scripts", "CCI_Paired_Comparison.R"))

  # 4 patients × 2 timepoints, 1 LR pair (CD8_T → B_cell, L1→R1).
  # T12 score is 0.8 in every patient; T0 is 0.2. Paired delta = +0.6.
  patients <- 1:4
  long <- do.call(rbind, lapply(patients, function(p) {
    data.frame(
      sample_id = c(paste0("CART_T0_S", p), paste0("CART_T12_S", p)),
      source    = c("CD8_T", "CD8_T"),
      target    = c("B_cell", "B_cell"),
      ligand    = c("L1", "L1"),
      receptor  = c("R1", "R1"),
      score     = c(0.2, 0.8),
      stringsAsFactors = FALSE
    )
  }))
  sample_table <- data.frame(
    sample_id  = unique(long$sample_id),
    patient_id = rep(paste0("CART_pt", patients), each = 2),
    treatment  = "CART",
    timepoint  = rep(c("T0", "T12"), length(patients)),
    stringsAsFactors = FALSE
  )

  # Compute per-patient delta = T12 score - T0 score.
  deltas <- vapply(patients, function(p) {
    t0_id  <- paste0("CART_T0_S",  p)
    t12_id <- paste0("CART_T12_S", p)
    a <- long$score[long$sample_id == t12_id]
    b <- long$score[long$sample_id == t0_id]
    a[1] - b[1]
  }, numeric(1))
  expect_equal(unname(deltas), rep(0.6, length(patients)))

  # All-positive deltas → wilcox.test on n=4 with mu=0 produces a
  # finite p-value (warning about ties is suppressed).
  pv <- suppressWarnings(stats::wilcox.test(deltas, mu = 0, exact = FALSE)$p.value)
  expect_lt(pv, 0.5)
})

test_that("CCI_Paired_Comparison.R parses and exposes the public runner", {
  env <- .source_helper(file.path("Helper_Scripts", "CCI_Paired_Comparison.R"))
  expect_true(exists("run_cci_paired_comparison", envir = env, inherits = FALSE))
  expect_true(exists(".cci_tidy_liana", envir = env, inherits = FALSE))
})

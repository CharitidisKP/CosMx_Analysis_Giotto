# Phase 9 — paired GSEA-score sub-step. Synthetic per-sample NES tables;
# verify paired Wilcoxon on per-patient deltas detects a known signal.

test_that("paired Wilcoxon on synthetic NES deltas recovers a known shift", {
  # 4 patients × 2 timepoints. T12 NES − T0 NES = +1.0 for pathway P1
  # (always positive shift), 0 for pathway P2 (null).
  patients <- paste0("CART_pt", 1:4)
  pathways <- c("P1", "P2")
  per_sample <- list()
  for (i in seq_along(patients)) {
    p <- patients[i]
    sid_t0  <- paste0("CART_T0_S",  i)
    sid_t12 <- paste0("CART_T12_S", i)
    per_sample[[sid_t0]]  <- data.frame(pathway = pathways,
                                         NES = c(0, 0))
    per_sample[[sid_t12]] <- data.frame(pathway = pathways,
                                         NES = c(1, 0))
  }
  sids <- names(per_sample)

  # Build pathway × sample matrix.
  mat <- do.call(cbind, lapply(sids, function(s) {
    per_sample[[s]]$NES[match(pathways, per_sample[[s]]$pathway)]
  }))
  colnames(mat) <- sids; rownames(mat) <- pathways

  # Per-block deltas (T12 − T0).
  delta_mat <- vapply(patients, function(p) {
    a <- paste0("CART_T12_S", which(patients == p))
    b <- paste0("CART_T0_S",  which(patients == p))
    mat[, a] - mat[, b]
  }, numeric(length(pathways)))

  # P1: all +1 → strong rejection of H0
  v1 <- delta_mat["P1", ]
  expect_true(all(v1 == 1))
  pv1 <- suppressWarnings(stats::wilcox.test(v1, mu = 0, exact = FALSE)$p.value)
  expect_lt(pv1, 0.5)

  # P2: all 0 → null deltas; wilcox handles the all-tied case without crashing
  v2 <- delta_mat["P2", ]
  expect_true(all(v2 == 0))
})

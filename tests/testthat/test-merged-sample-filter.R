# Phase 1 — merged-side sample filter. Verifies:
#   - Default (no --samples): treatment %in% c("CART", "Conventional")
#   - Explicit --samples (any list): pass-through, Control allowed.

test_that("default merged filter excludes Control", {
  selected <- .synthetic_full_selection()
  default_treatments <- c("CART", "Conventional")
  explicit_samples <- FALSE

  merged <- if (isTRUE(explicit_samples)) {
    selected
  } else {
    keep <- selected$treatment %in% default_treatments
    selected[keep, , drop = FALSE]
  }

  expect_equal(nrow(merged), 8L)
  expect_false("1281" %in% merged$sample_id)
  expect_true(all(merged$treatment %in% c("CART", "Conventional")))
})

test_that("explicit --samples with Control is honored", {
  selected <- .synthetic_full_selection()
  # User asked for Control + 4 CART biopsies = 5 samples.
  ex <- c("1281", "CART_T0_S1", "CART_T12_S1", "CART_T0_S2", "CART_T12_S2")
  selected <- selected[selected$sample_id %in% ex, , drop = FALSE]
  explicit_samples <- TRUE

  merged <- if (isTRUE(explicit_samples)) {
    selected
  } else {
    keep <- selected$treatment %in% c("CART", "Conventional")
    selected[keep, , drop = FALSE]
  }

  expect_equal(nrow(merged), 5L)
  expect_true("1281" %in% merged$sample_id)
})

test_that("no-match fallback returns full selection (auto-skip handles <2)", {
  # User passes some weird selection where no treatment matches the default.
  selected <- data.frame(
    sample_id = c("X1", "X2"),
    treatment = c("Other", "Other"),
    stringsAsFactors = FALSE
  )
  default_treatments <- c("CART", "Conventional")
  keep <- selected$treatment %in% default_treatments
  merged <- if (!any(keep)) selected else selected[keep, , drop = FALSE]
  expect_equal(nrow(merged), 2L)
})

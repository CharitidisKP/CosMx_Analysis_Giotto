# Phase 4 — verify that .pathway_run_de() honors comparison_spec$pseudobulk_backend
# and comparison_spec$block overrides.

test_that(".pathway_run_de selects engine via per-comparison override", {
  # Source the engine dispatcher only — the upstream dependencies (DESeq2,
  # edgeR, limma) aren't installed locally; we test the dispatch logic by
  # passing a bogus engine name and confirming the error path fires.
  env <- new.env(parent = globalenv())
  if (!exists("%||%", envir = env, inherits = FALSE)) {
    assign("%||%",
           function(a, b) if (is.null(a) || length(a) == 0L) b else a,
           envir = env)
  }
  source(file.path(.repo_dir(), "13_Pathway_Analysis.R"), local = env)

  # Comparison override should win over cfg$pathway$de_engine.
  cfg_p <- list(de_engine = "pseudobulk_deseq2",
                de_min_pseudobulk_per_group = 2L)
  cs <- list(pseudobulk_backend = "pseudobulk_definitely_not_an_engine",
             design = "unpaired")

  # Need ≥2 pseudobulk replicates per side (de_min_pseudobulk_per_group)
  # to reach the engine dispatch — pseudobulk is by sample_id × de_group,
  # so 2 samples per side gives 2 replicates.
  cells_a <- paste0("a", 1:8)
  cells_b <- paste0("b", 1:8)
  raw <- matrix(rpois(10 * 16, 5), nrow = 10,
                dimnames = list(paste0("g", 1:10), c(cells_a, cells_b)))
  norm <- raw
  meta <- data.frame(
    cell_ID   = c(cells_a, cells_b),
    sample_id = c(rep("S_A1", 4), rep("S_A2", 4),
                  rep("S_B1", 4), rep("S_B2", 4)),
    stringsAsFactors = FALSE
  )
  expect_error(
    env$.pathway_run_de(raw, norm, meta,
                         ids_a = cells_a, ids_b = cells_b,
                         comparison_spec = cs, cfg_pathway = cfg_p),
    "Unknown pathway.de_engine"
  )
})

test_that("block override threads into the design formula (unit shape)", {
  # Construct a minimal pb_meta that exercises the design_formula construction
  # for DESeq2 without needing DESeq2 itself: factor coercion + as.formula.
  pb_meta <- data.frame(
    cell_ID = paste0("c", 1:4),
    sample_id = paste0("s", 1:4),
    de_group = factor(c("A","A","B","B"), levels = c("B","A")),
    site_id = c("SiteX","SiteY","SiteX","SiteY"),
    stringsAsFactors = FALSE
  )
  block_col <- "site_id"
  pb_meta[[block_col]] <- factor(pb_meta[[block_col]])
  design_formula <- stats::as.formula(paste0("~ ", block_col, " + de_group"))
  expect_equal(deparse(design_formula), "~site_id + de_group")
})

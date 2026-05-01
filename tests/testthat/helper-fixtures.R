# Synthetic fixtures used by Stage 1 helper tests. Pure R / Matrix only;
# does NOT depend on Giotto / Banksy / kBET / scran being installed.

# Resolve the project root from the testthat invocation directory.
.repo_dir <- function() {
  here <- normalizePath(getwd(), mustWork = FALSE)
  candidate <- here
  for (i in 1:5) {
    if (file.exists(file.path(candidate, "CosMx_pipeline.R"))) {
      return(candidate)
    }
    candidate <- dirname(candidate)
  }
  here  # best effort; tests will fail loudly if path is wrong
}

# Source a Helper_Scripts file relative to repo root, isolated in `env`.
# Sets cosmx.disable_cli = TRUE for the duration of the source so that
# top-level CLI dispatch blocks (which check
# `!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))`) skip
# their stop("Usage: ...") branch when sourced under Rscript.
.source_helper <- function(rel_path, env = new.env(parent = globalenv())) {
  path <- file.path(.repo_dir(), rel_path)
  if (!file.exists(path)) skip(paste0("helper not found: ", rel_path))
  # Define %||% in env so helpers that reference it before defining it work.
  if (!exists("%||%", envir = env, inherits = FALSE)) {
    assign("%||%", function(a, b) if (is.null(a) || length(a) == 0L) b else a,
           envir = env)
  }
  prev <- options(cosmx.disable_cli = TRUE)
  on.exit(options(prev), add = TRUE)
  source(path, local = env)
  env
}

# Build a synthetic per-sample cell metadata frame for Composition / paired
# tests. Eight paired samples, 4 patients × 2 timepoints, 2 arms.
.synthetic_paired_table <- function() {
  data.frame(
    sample_id = c("CART_T0_S1", "CART_T12_S1", "CART_T0_S2", "CART_T12_S2",
                  "CONV_T0_P1", "CONV_T12_P1", "CONV_T0_P2", "CONV_T12_P2"),
    patient_id = c("CART_pt1", "CART_pt1", "CART_pt2", "CART_pt2",
                   "Conv_pt1", "Conv_pt1", "Conv_pt2", "Conv_pt2"),
    treatment = c(rep("CART", 4), rep("Conventional", 4)),
    timepoint = c("T0", "T12", "T0", "T12", "T0", "T12", "T0", "T12"),
    stringsAsFactors = FALSE
  )
}

# Build a synthetic 9-sample selection (the 8 paired + 1 Control) modelling
# the user's actual sample sheet for the Phase 1 filter test.
.synthetic_full_selection <- function() {
  rbind(
    .synthetic_paired_table(),
    data.frame(
      sample_id = "1281", patient_id = "Control_pt1",
      treatment = "Control", timepoint = "T0",
      stringsAsFactors = FALSE
    )
  )
}

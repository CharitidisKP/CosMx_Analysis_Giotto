#!/usr/bin/env Rscript
# ==============================================================================
# Install_Stage1_Dependencies.R
# Installs the new R packages introduced by Stage 1 of the merged-run
# enhancement: 8 runtime packages (kBET, lisi, Banksy, speckle, limma,
# UCell, scran, rmarkdown) + testthat for the unit-test suite.
#
# Run from inside the project's Apptainer container so the project's
# .Rprofile activates the host-mounted renv library (packages persist
# across container runs):
#
#   apptainer exec --cleanenv --pwd "$PROJECT_DIR" \
#     --bind "$HOME":"$HOME" --bind "$TMPDIR_HOST":/tmp \
#     "$SIF" Rscript "$PROJECT_DIR/Scripts/Parameters/Install_Stage1_Dependencies.R"
#
# Notes:
#   - Banksy / kBET / lisi may pull in C++ deps; if a compile fails inside
#     the existing .sif, the most likely fix is to rebuild the image with
#     the missing system library. All other Stage-1 helpers wrap their
#     calls in requireNamespace() so the rest of the pipeline remains
#     functional even if one package fails to install.
#   - Pandoc availability is asserted at the end so the Phase 13 summary
#     report won't fail to render later.
# ==============================================================================

cat("=================================================\n")
cat("Installing Stage 1 dependencies (8 packages)\n")
cat("=================================================\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Bioconductor, pure R for the most part; should compile cleanly inside
# the existing image. update = FALSE keeps already-pinned versions intact.
cat("-> Bioconductor: Banksy, speckle, limma, UCell, scran\n")
BiocManager::install(c("Banksy", "speckle", "limma", "UCell", "scran"),
                     update = FALSE, ask = FALSE)

# CRAN, rmarkdown is needed for the Phase 13 summary report; testthat is
# needed to run the Stage 1 unit-test suite inside the container. Both are
# idempotent.
cat("\n-> CRAN: rmarkdown, testthat\n")
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  install.packages("rmarkdown")
}
if (!requireNamespace("testthat", quietly = TRUE)) {
  install.packages("testthat")
}

# GitHub-only, kBET (theislab) and lisi (immunogenomics).
cat("\n-> GitHub: kBET (theislab), lisi (immunogenomics)\n")
remotes::install_github("theislab/kBET",       upgrade = "never")
remotes::install_github("immunogenomics/LISI", upgrade = "never")

# Snapshot to renv.lock so subsequent containers reproduce the install.
# Non-fatal: pre-existing renv-state issues (e.g. unrelated packages missing
# or BioC-release skew) should not block the Stage 1 install. The packages
# are installed regardless; the user can run `renv::status()` and
# `renv::snapshot()` interactively later to update the lockfile.
if (requireNamespace("renv", quietly = TRUE)) {
  cat("\n-> renv::snapshot()\n")
  snap_ok <- tryCatch({ renv::snapshot(prompt = FALSE); TRUE },
                      error = function(e) {
                        cat("  Warning: renv::snapshot() aborted:\n    ",
                            conditionMessage(e), "\n",
                            "    The install itself succeeded. Run\n",
                            "    `renv::status()` to inspect, then\n",
                            "    `renv::snapshot()` interactively when\n",
                            "    the unrelated issues are resolved.\n",
                            sep = "")
                        FALSE
                      })
}

# Verify all nine packages load cleanly (eight runtime + testthat for tests).
cat("\n=== Version report ===\n")
for (pkg in c("Banksy", "speckle", "limma", "UCell", "scran",
              "rmarkdown", "kBET", "lisi", "testthat")) {
  v <- tryCatch(as.character(packageVersion(pkg)),
                error = function(e) "MISSING")
  cat(sprintf("  %-12s %s\n", pkg, v))
}

# Pandoc is needed by rmarkdown::render() for the Phase 13 summary report.
# Stage 1 phases 1-11 are unaffected if pandoc is absent; only Phase 13
# will fail.
cat("\n=== Pandoc availability ===\n")
pandoc_ok <- tryCatch(rmarkdown::pandoc_available(),
                      error = function(e) FALSE)
cat("  rmarkdown::pandoc_available() = ", pandoc_ok, "\n", sep = "")
if (!isTRUE(pandoc_ok)) {
  cat("  Warning: Phase 13 (merged-run summary) will fail to render until\n",
      "    pandoc is available inside the .sif (system-level install).\n",
      sep = "")
}

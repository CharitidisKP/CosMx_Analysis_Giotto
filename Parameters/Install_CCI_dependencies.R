#!/usr/bin/env Rscript
# ==============================================================================
# Install_CCI_Dependencies.R
# One-time installer for optional dependencies used by 10_CCI_Analysis.R
# and optional spatial DE backends such as smiDE.
# ==============================================================================

cci_dependency_registry <- function() {
  list(
    insitucor = list(
      pkg = "InSituCor",
      source = "github",
      install = function() {
        ensure_installer_package("remotes")
        remotes::install_github(
          "Nanostring-Biostats/InSituCor",
          upgrade = "never"
        )
      }
    ),
    liana = list(
      pkg = "liana",
      source = "github",
      install = function() {
        ensure_installer_package("remotes")
        remotes::install_github(
          "saezlab/liana",
          upgrade = "never"
        )
      }
    ),
    seurat = list(
      pkg = "Seurat",
      source = "cran",
      install = function() install.packages("Seurat", repos = "https://cloud.r-project.org")
    ),
    nichenet = list(
      pkg = "nichenetr",
      source = "github",
      install = function() {
        ensure_installer_package("remotes")
        remotes::install_github(
          "saeyslab/nichenetr",
          upgrade = "never"
        )
      }
    ),
    misty = list(
      pkg = "mistyR",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("mistyR", ask = FALSE, update = FALSE)
      }
    ),
    nnsvg = list(
      pkg = "nnSVG",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("nnSVG", ask = FALSE, update = FALSE)
      }
    ),
    seraster = list(
      pkg = "SEraster",
      source = "github",
      install = function() {
        ensure_installer_package("remotes")
        remotes::install_github(
          "JEFworks-Lab/SEraster",
          upgrade = "never"
        )
      }
    ),
    spatialexperiment = list(
      pkg = "SpatialExperiment",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("SpatialExperiment", ask = FALSE, update = FALSE)
      }
    ),
    smide = list(
      pkg = "smiDE",
      source = "github",
      install = function() {
        ensure_installer_package("remotes")
        remotes::install_github(
          "Nanostring-Biostats/CosMx-Analysis-Scratch-Space",
          subdir = "_code/smiDE",
          ref = "Main"
        )
      }
    ),

    # ------------------------------------------------------------------
    # Pathway enrichment + GSEA stack used by 13_Pathway_Analysis.R.
    # Install with: Rscript Parameters/Install_CCI_dependencies.R \
    #   msigdbr,clusterprofiler,fgsea,enrichplot,reactomepa,orghs,annotationdbi,\
    #   decoupler,gseabase,progeny,presto,pheatmap,upsetr,complexheatmap
    # ------------------------------------------------------------------
    msigdbr = list(
      pkg = "msigdbr",
      source = "cran",
      install = function() install.packages("msigdbr", repos = "https://cloud.r-project.org")
    ),
    clusterprofiler = list(
      pkg = "clusterProfiler",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
      }
    ),
    fgsea = list(
      pkg = "fgsea",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("fgsea", ask = FALSE, update = FALSE)
      }
    ),
    enrichplot = list(
      pkg = "enrichplot",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("enrichplot", ask = FALSE, update = FALSE)
      }
    ),
    reactomepa = list(
      pkg = "ReactomePA",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("ReactomePA", ask = FALSE, update = FALSE)
      }
    ),
    orghs = list(
      pkg = "org.Hs.eg.db",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)
      }
    ),
    annotationdbi = list(
      pkg = "AnnotationDbi",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
      }
    ),
    decoupler = list(
      pkg = "decoupleR",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("decoupleR", ask = FALSE, update = FALSE)
      }
    ),
    gseabase = list(
      pkg = "GSEABase",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("GSEABase", ask = FALSE, update = FALSE)
      }
    ),
    progeny = list(
      pkg = "progeny",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("progeny", ask = FALSE, update = FALSE)
      }
    ),
    presto = list(
      pkg = "presto",
      source = "github",
      install = function() {
        ensure_installer_package("remotes")
        remotes::install_github("immunogenomics/presto", upgrade = "never")
      }
    ),
    pheatmap = list(
      pkg = "pheatmap",
      source = "cran",
      install = function() install.packages("pheatmap", repos = "https://cloud.r-project.org")
    ),
    upsetr = list(
      pkg = "UpSetR",
      source = "cran",
      install = function() install.packages("UpSetR", repos = "https://cloud.r-project.org")
    ),
    complexheatmap = list(
      pkg = "ComplexHeatmap",
      source = "bioc",
      install = function() {
        ensure_installer_package("BiocManager")
        BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
      }
    ),

    # ------------------------------------------------------------------
    # Parallelism stack used by sample-level parallelism (CosMx_pipeline.R)
    # and smiDE per-cell-type parallelism (12_Spatial_Differential_Expression.R).
    # All three are pure CRAN, no Bioconductor dependency.
    # Install with: Rscript Parameters/Install_CCI_dependencies.R \
    #   future,futureapply,progressr
    # ------------------------------------------------------------------
    future = list(
      pkg = "future",
      source = "cran",
      install = function() install.packages("future", repos = "https://cloud.r-project.org")
    ),
    futureapply = list(
      pkg = "future.apply",
      source = "cran",
      install = function() install.packages("future.apply", repos = "https://cloud.r-project.org")
    ),
    progressr = list(
      pkg = "progressr",
      source = "cran",
      install = function() install.packages("progressr", repos = "https://cloud.r-project.org")
    )
  )
}

ensure_installer_package <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  install.packages(pkg, repos = "https://cloud.r-project.org")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Failed to install required installer package: ", pkg)
  }
  invisible(TRUE)
}

report_cci_dependency_status <- function(registry = cci_dependency_registry()) {
  status <- vapply(registry, function(entry) {
    requireNamespace(entry$pkg, quietly = TRUE)
  }, logical(1))
  
  cat("\n=== Optional CCI dependency status ===\n")
  for (nm in names(registry)) {
    icon <- if (status[[nm]]) "\u2713" else "\u2717"
    cat(sprintf("  %s %-18s (%s)\n", icon, registry[[nm]]$pkg, nm))
  }
  cat("\n")
  
  invisible(status)
}

install_cci_dependencies <- function(targets = names(cci_dependency_registry()),
                                     force = FALSE,
                                     registry = cci_dependency_registry()) {
  unknown <- setdiff(targets, names(registry))
  if (length(unknown) > 0) {
    stop("Unknown targets: ", paste(unknown, collapse = ", "))
  }
  
  for (target in targets) {
    entry <- registry[[target]]
    installed <- requireNamespace(entry$pkg, quietly = TRUE)
    
    if (installed && !force) {
      cat("\u2713", entry$pkg, "already installed; skipping\n")
      next
    }
    
    cat("Installing", entry$pkg, "from", entry$source, "...\n")
    entry$install()
    
    if (!requireNamespace(entry$pkg, quietly = TRUE)) {
      stop("Installation did not make package available: ", entry$pkg)
    }
    
    cat("\u2713", entry$pkg, "installed\n")
  }
  
  invisible(TRUE)
}

parse_targets <- function(args) {
  if (length(args) == 0) {
    return(names(cci_dependency_registry()))
  }
  unique(trimws(unlist(strsplit(args[[1]], ",", fixed = TRUE))))
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  targets <- parse_targets(args)
  report_cci_dependency_status()
  install_cci_dependencies(targets = targets)
  report_cci_dependency_status()
}

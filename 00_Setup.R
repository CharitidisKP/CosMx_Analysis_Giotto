# Load all required packages for CosMx analysis pipeline ------------------
#!/usr/bin/env Rscript
# ==============================================================================
# 00_setup_environment.R
# Session-aware: won't reload packages that are already loaded
# ==============================================================================

current_setup_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE)))
  }
  
  ofiles <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else frame$ofile
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles) > 0) {
    return(dirname(normalizePath(tail(ofiles, 1), winslash = "/", mustWork = FALSE)))
  }
  
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

setup_environment <- function(verbose = TRUE) {
  
  if (verbose) {
    cat("\n========================================\n")
    cat("STEP 00: Setting Up Environment\n")
    cat("========================================\n\n")
  }
  
  start_time <- Sys.time()
  
  project_dir <- Sys.getenv("COSMX_PROJECT_DIR", unset = current_setup_dir())
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = FALSE)
  
  python_candidates <- unique(c(
    getOption("cosmx.python_path", NULL),
    Sys.which("python3"),
    Sys.which("python")
  ))
  python_candidates <- python_candidates[
    nzchar(python_candidates) & vapply(
      python_candidates,
      function(path) file.exists(path) || nzchar(Sys.which(path)),
      logical(1)
    )
  ]
  
  python_path <- Sys.getenv("COSMX_PYTHON_PATH", unset = "")
  if (!nzchar(python_path)) {
    python_path <- if (length(python_candidates) > 0) python_candidates[[1]] else "python3"
  }
  Sys.setenv(
    COSMX_PYTHON_PATH = python_path,
    RETICULATE_PYTHON = python_path
  )

  # Pin reticulate to the chosen python BEFORE other packages load, so any
  # package calling reticulate::use_python / use_condaenv in its .onLoad
  # (Giotto, SpatialDecon, etc.) finds python already initialised and cannot
  # override it back to /usr/bin/python3.
  if (nzchar(python_path) && file.exists(python_path) &&
      requireNamespace("reticulate", quietly = TRUE)) {
    tryCatch({
      reticulate::use_python(python_path, required = TRUE)
      reticulate::py_run_string("pass")
    }, error = function(e) {
      if (verbose) cat("  ⚠ early reticulate pin failed:", conditionMessage(e), "\n")
    })
  }

  # Define package list in STRICT load order
  packages <- c(
    "Matrix", "data.table", "tibble", "dplyr", "tidyr", "readr",
    "purrr", "stringr", "ggplot2", "terra", "sf", "viridis",
    "RColorBrewer", "scales", "patchwork", "cowplot", "ggpubr",
    "pheatmap", "gridExtra", "grid", "knitr", "kableExtra",
    "igraph", "tidygraph", "ggraph", "scran", "SingleR",
    "celldex", "InSituType", "reticulate", "SpatialDecon", "Giotto",
    "yaml", "optparse", "readxl", "jsonlite", "processx", "ggrepel",
    "edgeR", "smiDE",
    "harmony", "shiny", "bslib", "DT"
  )
  
  # Check which packages are already loaded
  already_loaded <- sapply(packages, function(pkg) {
    paste0("package:", pkg) %in% search()
  })
  
  packages_to_load <- packages[!already_loaded]
  
  if (verbose) {
    if (sum(already_loaded) > 0) {
      cat("Already loaded:", sum(already_loaded), "packages\n")
      cat("  Loaded:", paste(names(already_loaded)[already_loaded], collapse = ", "), "\n\n")
    }
    if (length(packages_to_load) > 0) {
      cat("Need to load:", length(packages_to_load), "packages\n")
      cat("  To load:", paste(packages_to_load, collapse = ", "), "\n\n")
    }
  }
  
  allow_install <- tolower(Sys.getenv("COSMX_ALLOW_INSTALL", unset = "false")) %in%
    c("1", "true", "yes", "y")
  
  manual_install_packages <- c("InSituType", "SpatialDecon", "smiDE", "liana", "nichenetr", "mistyR")
  
  # Function to safely load package
  load_package <- function(pkg) {
    # First check if already loaded
    if (paste0("package:", pkg) %in% search()) {
      return(TRUE)
    }
    
    # Check if installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (allow_install && !pkg %in% manual_install_packages) {
        if (verbose) cat("  Installing", pkg, "...\n")
        suppressWarnings(
          install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
        )
      } else {
        if (verbose) {
          if (pkg %in% manual_install_packages) {
            cat("  ✗", pkg, "- not installed (install manually; non-CRAN dependency)\n")
          } else {
            cat("  ✗", pkg, "- not installed\n")
          }
        }
        return(FALSE)
      }
    }
    
    # Try to load
    success <- suppressPackageStartupMessages({
      require(pkg, character.only = TRUE, quietly = TRUE)
    })
    
    return(success)
  }
  
  # Load packages
  if (length(packages_to_load) > 0) {
    if (verbose) cat("Loading new packages...\n")
    
    failed_packages <- c()
    for (pkg in packages_to_load) {
      success <- load_package(pkg)
      if (success) {
        if (verbose) cat("  ✓", pkg, "\n")
      } else {
        if (verbose) cat("  ✗", pkg, "- FAILED\n")
        failed_packages <- c(failed_packages, pkg)
      }
    }
    
    if (verbose) cat("\n")
    
    # Only stop if critical packages failed
    critical_failed <- intersect(failed_packages,
                                 c("Giotto", "Matrix", "terra", "ggplot2"))
    if (length(critical_failed) > 0) {
      stop("Critical packages failed to load: ",
           paste(critical_failed, collapse = ", "))
    }
    
    if (length(failed_packages) > 0) {
      warning("Some packages failed to load: ",
              paste(failed_packages, collapse = ", "))
    }
  } else {
    if (verbose) cat("✓ All required packages already loaded\n\n")
  }
  
  # Verify critical packages
  critical_packages <- c("Giotto", "Matrix", "terra", "ggplot2")
  loaded_check <- sapply(critical_packages, function(pkg) {
    paste0("package:", pkg) %in% search()
  })
  
  if (!all(loaded_check)) {
    missing <- names(loaded_check)[!loaded_check]
    stop("Critical packages not available: ", paste(missing, collapse = ", "))
  }
  
  if (verbose) cat("✓ All critical packages verified\n\n")
  
  # Configure Python
  if (verbose) cat("Configuring Python...\n")
  tryCatch({
    if (!reticulate::py_available(initialize = FALSE)) {
      reticulate::use_python(python_path, required = TRUE)
    }
    py_config <- reticulate::py_config()
    if (verbose) {
      cat("✓ Python:", as.character(py_config$version), "\n")
      cat("  Path:", py_config$python, "\n")
    }
    
    required_python_modules <- c("umap", "igraph", "leidenalg")
    missing_python_modules <- required_python_modules[!vapply(
      required_python_modules,
      reticulate::py_module_available,
      logical(1)
    )]
    if (length(missing_python_modules) > 0) {
      warning(
        "Python module(s) missing in ", py_config$python, ": ",
        paste(missing_python_modules, collapse = ", "),
        ". Giotto functions that rely on them may fail."
      )
      if (verbose) {
        cat("  ⚠ missing python modules:", paste(missing_python_modules, collapse = ", "), "\n")
      }
    }
    if (verbose) cat("\n")
  }, error = function(e) {
    if (verbose) cat("⚠ Python warning, continuing...\n\n")
  })
  
  # Source helpers
  if (verbose) cat("Loading helper functions...\n")
  helper_dir_env <- Sys.getenv("COSMX_HELPER_DIR", unset = "")
  helper_candidates <- c(
    helper_dir_env,
    file.path(project_dir, "Helper_Scripts"),
    file.path(project_dir, "Scripts", "Helper_Scripts"),
    "~/P_lab/CosMx_analysis/Scripts/Helper_Scripts"
  )
  helper_candidates <- unique(helper_candidates[nzchar(helper_candidates)])
  helper_dir <- helper_candidates[file.exists(file.path(helper_candidates, "Helper_Functions.R"))][1]
  
  if (is.na(helper_dir) || !nzchar(helper_dir)) {
    helper_dir <- helper_candidates[1]
  }
  
  helper_files <- c(
    "Helper_Functions.R",
    "Merge_Batch_Correction.R",
    "Inspect_sNN_network.R",
    "Remove_HVF_duplicates.R",
    "Cluster_Visualisations.R",
    "Arrange_Feature_plots.R",
    "Feature_plots_panel.R",
    "CCI_Summary.R",
    "Plot_Helpers.R"
  )
  
  sourced_count <- 0
  for (hf in helper_files) {
    path <- file.path(helper_dir, hf)
    if (file.exists(path)) {
      tryCatch({
        source(path)
        sourced_count <- sourced_count + 1
        if (verbose) cat("  ✓", hf, "\n")
      }, error = function(e) {
        if (verbose) cat("  ⚠", hf, "\n")
      })
    } else {
      if (verbose) cat("  ⚠", hf, "- not found\n")
    }
  }
  
  if (verbose) {
    cat("\n✓ Environment setup complete\n")
    cat("  Project dir:", project_dir, "\n")
    cat("  Packages checked:", length(packages), "\n")
    cat("  Helpers loaded:", sourced_count, "/", length(helper_files), "\n\n")
  }
  
  return(invisible(list(
    project_dir = project_dir,
    python_path = python_path,
    helper_dir = helper_dir,
    packages_loaded = sum(sapply(packages, function(p) paste0("package:", p) %in% search()))
  )))
}

# Auto-execute
if (!exists(".cosmx_env_setup_complete")) {
  env_config <- setup_environment()
  .cosmx_project_dir <<- env_config$project_dir
  python_path <<- env_config$python_path
  .cosmx_env_setup_complete <<- TRUE
} else {
  cat("\n✓ Environment already initialized\n\n")
  env_config <- list(
    project_dir = get(".cosmx_project_dir", envir = .GlobalEnv),
    python_path = get("python_path", envir = .GlobalEnv),
    helper_dir = Sys.getenv("COSMX_HELPER_DIR", unset = file.path(getwd(), "Helper_Scripts"))
  )
}

invisible(env_config)

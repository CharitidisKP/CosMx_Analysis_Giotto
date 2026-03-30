#!/usr/bin/env Rscript
# ==============================================================================
# 10_CCI_Analysis.R
# Cell-cell interaction (Layer 2) and spatially-resolved DE (Layer 3)
#
# Optional dependencies for this script should be installed separately via
# Parameters/Install_CCI_Dependencies.R. This analysis script does not mutate the R
# environment while running.
#
# INPUT:  Giotto object saved by 09_Spatial_Network.R
#         (contains annotation + spatial networks + PAGE/rank enrichment)
# OUTPUT: output_dir/10_CCI_Analysis/
#           ├── insitucor/    Spatially co-expressed gene modules (NanoString)
#           ├── liana/        Consensus L-R scores (all methods)
#           ├── nichenet/     Ligand activity scores
#           ├── misty/        Intercellular spatial modelling
#           ├── svg/          Spatially variable genes (nnSVG)
#           └── Summary/      Step-level summary tables and overview plots
# ==============================================================================


# ==============================================================================
# SECTION 0 — InSituCor (NanoString-native spatially co-expressed modules)
# ==============================================================================

.giotto_get_expression <- function(gobj, values, output = "matrix") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("getExpression", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("getExpression", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getExpression", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("getExpression", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("getExpression", mode = "function")
  }
  
  accessor(gobj, values = values, output = output)
}

.giotto_pdata_dt <- function(gobj) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("pDataDT", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("pDataDT", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("pDataDT", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("pDataDT", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("pDataDT", mode = "function")
  }
  
  accessor(gobj)
}

.giotto_get_spatial_locations <- function(gobj, output = "data.table") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("getSpatialLocations", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("getSpatialLocations", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getSpatialLocations", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("getSpatialLocations", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("getSpatialLocations", mode = "function")
  }
  
  accessor(gobj, output = output)
}

.sanitize_feature_names <- function(x) {
  safe <- make.names(x, unique = TRUE)
  names(safe) <- x
  safe
}

.first_existing_column <- function(df, candidates) {
  cols <- intersect(candidates, colnames(df))
  if (length(cols) == 0) {
    return(NULL)
  }
  cols[1]
}

.resolve_nichenet_prior_files <- function(network_dir) {
  required_sets <- list(
    v2 = c(
      ligand_target_matrix = "ligand_target_matrix_nsga2r_final.rds",
      lr_network = "lr_network_human_21122021.rds",
      weighted_networks = "weighted_networks_nsga2r_final.rds"
    ),
    legacy = c(
      ligand_target_matrix = "ligand_target_matrix.rds",
      lr_network = "lr_network.rds",
      weighted_networks = "weighted_networks.rds"
    )
  )
  
  for (set_name in names(required_sets)) {
    candidate_paths <- file.path(network_dir, unname(required_sets[[set_name]]))
    if (all(file.exists(candidate_paths))) {
      names(candidate_paths) <- names(required_sets[[set_name]])
      return(list(
        version = set_name,
        paths = candidate_paths
      ))
    }
  }
  
  expected_paths <- lapply(required_sets, function(files) {
    stats::setNames(file.path(network_dir, unname(files)), names(files))
  })
  
  missing_paths <- lapply(expected_paths, function(files) {
    files[!file.exists(files)]
  })
  
  list(
    version = NULL,
    paths = NULL,
    missing = missing_paths
  )
}

.subset_expression_columns <- function(expr_mat, cell_ids) {
  matched_ids <- intersect(as.character(cell_ids), colnames(expr_mat))
  if (length(matched_ids) == 0) {
    return(NULL)
  }
  
  subset_mat <- expr_mat[, matched_ids, drop = FALSE]
  
  if (is.null(dim(subset_mat))) {
    subset_mat <- matrix(
      subset_mat,
      nrow = nrow(expr_mat),
      dimnames = list(rownames(expr_mat), matched_ids)
    )
  }
  
  subset_mat
}

.row_means_matrix_like <- function(x) {
  if (inherits(x, "Matrix")) {
    return(Matrix::rowMeans(x))
  }
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  rowMeans(as.matrix(x))
}

.expressed_genes_from_matrix <- function(x, pct = 0.10, min_expr = 0) {
  if (is.null(x) || ncol(x) == 0) {
    return(character(0))
  }
  
  detected_fraction <- if (inherits(x, "Matrix")) {
    Matrix::rowMeans(x > min_expr)
  } else {
    rowMeans(as.matrix(x) > min_expr)
  }
  
  names(detected_fraction)[detected_fraction >= pct]
}

.load_cci_summary_helper <- function() {
  if (exists("create_cci_summary", mode = "function", inherits = TRUE)) {
    return(TRUE)
  }
  
  helper_candidates <- unique(c(
    Sys.getenv("COSMX_HELPER_DIR", unset = ""),
    file.path(Sys.getenv("COSMX_PROJECT_DIR", unset = getwd()), "Helper_Scripts"),
    file.path(Sys.getenv("COSMX_PROJECT_DIR", unset = getwd()), "Scripts", "Helper_Scripts"),
    file.path(getwd(), "Helper_Scripts")
  ))
  helper_candidates <- helper_candidates[nzchar(helper_candidates)]
  
  helper_file <- file.path(helper_candidates, "CCI_Summary.R")
  helper_file <- helper_file[file.exists(helper_file)][1]
  if (is.na(helper_file) || !nzchar(helper_file)) {
    return(FALSE)
  }
  
  source(helper_file, local = .GlobalEnv)
  exists("create_cci_summary", mode = "function", inherits = TRUE)
}

#' Run InSituCor: spatially co-expressed gene modules
#'
#' Conditions on local cell-type composition, so co-expression found here is
#' truly spatial (above and beyond what cell type alone explains).  This is the
#' natural first step before L-R analysis — it tells you which gene programmes
#' are active in each spatial niche.
#'
#' @param gobj          Annotated Giotto object
#' @param sample_id     Sample identifier
#' @param output_dir    Root output directory
#' @param celltype_col  Cell type annotation column
#' @param n_cores       Parallel cores (default: 4)
#' @return List with insitucor results; gobj unchanged

run_insitucor <- function(gobj,
                          sample_id,
                          output_dir,
                          celltype_col = NULL,
                          n_cores      = 4) {
  
  if (!requireNamespace("InSituCor", quietly = TRUE))
    stop("InSituCor not installed.\n",
         'Install: devtools::install_github("Nanostring-Biostats/InSituCor")')
  
  cat("\n--- InSituCor: Spatially co-expressed gene modules ---\n")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "insitucor")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Resolve celltype column
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (is.null(celltype_col)) {
    sup_cols <- grep("celltype_.*_supervised$", names(meta), value = TRUE)
    celltype_col <- if (length(sup_cols) > 0) sup_cols[1] else "leiden_clust"
  }
  
  # Expression matrix (cells x genes)
  # InSituCor models neighborhood expression from count-like data, so use raw counts.
  counts <- t(as.matrix(
    .giotto_get_expression(gobj, values = "raw", output = "matrix")
  ))
  
  # Spatial coordinates
  spat <- as.data.frame(
    .giotto_get_spatial_locations(gobj, output = "data.table")
  )
  rownames(spat) <- spat$cell_ID
  
  # Cell type vector — must align with rows of counts
  ct_vec <- setNames(meta[[celltype_col]], meta$cell_ID)
  ct_vec <- ct_vec[rownames(counts)]
  valid_cells <- !is.na(ct_vec) &
    rownames(counts) %in% rownames(spat) &
    stats::complete.cases(spat[rownames(counts), c("sdimx", "sdimy"), drop = FALSE])
  counts <- counts[valid_cells, , drop = FALSE]
  ct_vec <- as.character(ct_vec[valid_cells])
  xy <- spat[rownames(counts), c("sdimx", "sdimy"), drop = FALSE]
  if (nrow(counts) == 0) {
    stop("No valid cells with both cell type labels and spatial coordinates for InSituCor.")
  }
  
  meta_aligned <- meta[match(rownames(counts), meta$cell_ID), , drop = FALSE]
  total_counts <- rowSums(counts)
  conditionon_df <- data.frame(
    celltype = ct_vec,
    total_counts = as.numeric(total_counts),
    stringsAsFactors = FALSE,
    row.names = rownames(counts)
  )
  
  tissue_col <- intersect(
    c("fov", "fov_ID", "fov_id", "fov_num", "FOV", "slide", "slide_id"),
    colnames(meta_aligned)
  )
  tissue_vec <- NULL
  if (length(tissue_col) > 0) {
    candidate_tissue <- meta_aligned[[tissue_col[1]]]
    if (!all(is.na(candidate_tissue)) && length(unique(stats::na.omit(candidate_tissue))) > 1) {
      tissue_vec <- candidate_tissue
      conditionon_df$tissue <- as.character(candidate_tissue)
    }
  }
  
  cat("  Cells:", nrow(counts), "  Genes:", ncol(counts), "\n")
  cat("  Cell type column:", celltype_col, "\n")
  
  insitucor_fun <- get("insitucor", envir = asNamespace("InSituCor"))
  insitucor_formals <- names(formals(insitucor_fun))
  insitucor_args <- list(
    counts = counts,
    conditionon = conditionon_df,
    celltype = ct_vec,
    xy = xy
  )
  
  if ("n_cores" %in% insitucor_formals) {
    insitucor_args$n_cores <- n_cores
  }
  if ("n.cores" %in% insitucor_formals) {
    insitucor_args$n.cores <- n_cores
  }
  if ("ncores" %in% insitucor_formals) {
    insitucor_args$ncores <- n_cores
  }
  if ("k" %in% insitucor_formals) {
    insitucor_args$k <- 6
  } else if ("radius" %in% insitucor_formals) {
    insitucor_args$radius <- 50
  }
  if ("tissue" %in% insitucor_formals && !is.null(tissue_vec)) {
    insitucor_args$tissue <- tissue_vec
  }
  
  if ("verbose" %in% insitucor_formals) {
    insitucor_args$verbose <- TRUE
  }
  
  cor_results <- tryCatch(
    do.call(insitucor_fun, insitucor_args),
    error = function(e) {
      cat("\u26A0 InSituCor failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(cor_results)) {
    # Save module assignments
    saveRDS(cor_results,
            file.path(out_dir, paste0(sample_id, "_insitucor.rds")))
    
    if (!is.null(cor_results$modules)) {
      write.csv(cor_results$modules,
                file.path(out_dir, paste0(sample_id, "_insitucor_modules.csv")),
                row.names = FALSE)
    }
    cat("\u2713 InSituCor complete. Results saved to:", out_dir, "\n")
  }
  
  invisible(cor_results)
}


# ==============================================================================
# SECTION 1 — LIANA: consensus ligand-receptor scoring
# ==============================================================================

#' Run LIANA consensus ligand-receptor analysis
#'
#' Runs multiple LR methods (CellPhoneDB, CellChat, NATMI, Connectome, etc.)
#' and aggregates into a consensus rank.
#'
#' @param gobj          Annotated Giotto object
#' @param sample_id     Sample identifier
#' @param output_dir    Root output directory
#' @param celltype_col  Cell type annotation column
#' @param methods       LIANA methods to run (NULL = all available)
#' @return liana_result data frame

run_liana <- function(gobj,
                      sample_id,
                      output_dir,
                      celltype_col = NULL,
                      methods      = NULL) {
  
  if (!requireNamespace("liana", quietly = TRUE))
    stop("liana not installed.\n  Install: remotes::install_github('saezlab/liana')")
  
  cat("\n--- LIANA: Consensus ligand-receptor analysis ---\n")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "liana")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Resolve celltype column
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (is.null(celltype_col)) {
    sup_cols <- grep("celltype_.*_supervised$", names(meta), value = TRUE)
    celltype_col <- if (length(sup_cols) > 0) sup_cols[1] else "leiden_clust"
  }
  
  cat("  Converting Giotto to SingleCellExperiment for LIANA...\n")
  sce_obj <- tryCatch({
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
      stop("SingleCellExperiment required for LIANA conversion.")
    if (!requireNamespace("S4Vectors", quietly = TRUE))
      stop("S4Vectors required for LIANA conversion.")
    
    counts_mat <- .giotto_get_expression(gobj, values = "raw", output = "matrix")
    logcounts_mat <- .giotto_get_expression(gobj, values = "normalized", output = "matrix")
    cell_meta <- meta[match(colnames(counts_mat), meta$cell_ID), , drop = FALSE]
    
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(
        counts = counts_mat,
        logcounts = logcounts_mat
      ),
      colData = S4Vectors::DataFrame(cell_meta)
    )
    SingleCellExperiment::colLabels(sce) <- factor(cell_meta[[celltype_col]])
    sce
  }, error = function(e) {
    cat("\u26A0 SingleCellExperiment conversion failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(sce_obj)) return(invisible(NULL))
  
  if (is.null(methods)) {
    preferred_methods <- c("natmi", "connectome", "logfc", "sca", "cellphonedb")
    available_methods <- tryCatch(liana::show_methods(), error = function(e) preferred_methods)
    methods <- intersect(preferred_methods, available_methods)
    if (length(methods) == 0) methods <- available_methods
  }
  
  liana_res <- tryCatch(
    liana::liana_wrap(
      sce_obj,
      method    = methods,
      resource  = "Consensus",
      idents_col = celltype_col,
      assay.type = "logcounts",
      verbose   = TRUE
    ),
    error = function(e) {
      cat("\u26A0 LIANA failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(liana_res)) {
    liana_agg <- liana::liana_aggregate(liana_res)
    
    write.csv(liana_agg,
              file.path(out_dir, paste0(sample_id, "_liana_aggregate.csv")),
              row.names = FALSE)
    
    # Dot plot of top interactions
    tryCatch({
      specificity_col <- .first_existing_column(
        liana_agg,
        c(
          "natmi.edge_specificity",
          "connectome.edge_specificity",
          "cellphonedb.pvalue",
          "lr_means"
        )
      )
      magnitude_col <- .first_existing_column(
        liana_agg,
        c(
          "sca.LRscore",
          "connectome.scaled_weight",
          "logfc.logfc_comb",
          "natmi.edge_average",
          "lr_means"
        )
      )
      
      if (is.null(specificity_col) || is.null(magnitude_col)) {
        stop("Could not identify LIANA plot columns for the selected methods.")
      }
      
      p <- liana::liana_dotplot(
        liana_agg,
        source_groups = NULL,   # all senders
        target_groups = NULL,   # all receivers
        ntop          = 20,
        specificity   = specificity_col,
        magnitude     = magnitude_col
      )
      ggplot2::ggsave(
        file.path(out_dir, paste0(sample_id, "_liana_dotplot.png")),
        p, width = 14, height = 10, dpi = 150
      )
      cat("\u2713 LIANA dotplot saved\n")
    }, error = function(e) {
      cat("\u26A0 LIANA dotplot failed:", conditionMessage(e), "\n")
    })
    
    cat("\u2713 LIANA complete. Results saved to:", out_dir, "\n")
    invisible(liana_agg)
  } else {
    invisible(NULL)
  }
}


# ==============================================================================
# SECTION 2 — NicheNet: targeted ligand activity scoring
# ==============================================================================

#' Run NicheNet: ligand activity prediction for a sender-receiver pair
#'
#' Best used for hypothesis-driven questions, e.g.:
#'   sender   = "Myofibroblast"
#'   receiver = "Proximal.tubule"
#'
#' @param gobj          Annotated Giotto object
#' @param sample_id     Sample identifier
#' @param output_dir    Root output directory
#' @param celltype_col  Cell type annotation column
#' @param sender_celltypes  Character vector of sender cell types
#' @param receiver_celltype Single receiver cell type
#' @param top_n_ligands     Number of top ligands to report (default: 20)
#' @return NicheNet results list

run_nichenet <- function(gobj,
                         sample_id,
                         output_dir,
                         celltype_col       = NULL,
                         sender_celltypes   = NULL,
                         receiver_celltype  = NULL,
                         target_genes       = NULL,
                         top_n_ligands      = 20,
                         network_dir        = NULL) {
  
  if (!requireNamespace("nichenetr", quietly = TRUE))
    stop("nichenetr not installed.\n",
         'Install: devtools::install_github("saeyslab/nichenetr")')
  
  cat("\n--- NicheNet: Ligand activity prediction ---\n")
  cat("  Sender(s):  ", paste(sender_celltypes,  collapse = ", "), "\n")
  cat("  Receiver:   ", receiver_celltype, "\n\n")
  
  if (is.null(sender_celltypes) || is.null(receiver_celltype))
    stop("sender_celltypes and receiver_celltype must both be specified.")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "nichenet",
                       paste0(receiver_celltype, "_from_",
                              paste(sender_celltypes, collapse = "_")))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Resolve celltype column
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (is.null(celltype_col)) {
    sup_cols <- grep("celltype_.*_supervised$", names(meta), value = TRUE)
    celltype_col <- if (length(sup_cols) > 0) sup_cols[1] else "leiden_clust"
  }
  
  # NOTE: NicheNet requires pre-built networks (ligand-target, lr-network,
  # weighted networks). Download from the official NicheNet 2.0 Zenodo record:
  # https://doi.org/10.5281/zenodo.5884439
  cat("  \u26A0 NicheNet requires pre-built prior networks.\n")
  cat("    Download from: https://doi.org/10.5281/zenodo.5884439\n")
  cat("    Then set network_dir in run_nichenet() or load them manually.\n\n")
  
  network_dir <- network_dir %||%
    Sys.getenv("COSMX_NICHENET_DIR", unset = "")
  if (!nzchar(network_dir)) {
    network_dir <- file.path(
      Sys.getenv("HOME"),
      "P_lab", "CosMx_analysis",
      "Reference", "NicheNet_networks"
    )
  }
  
  if (!dir.exists(network_dir)) {
    cat("\u26A0 NicheNet network directory not found:", network_dir, "\n")
    cat("  Skipping NicheNet.\n\n")
    return(invisible(NULL))
  }
  
  prior_files <- .resolve_nichenet_prior_files(network_dir)
  if (is.null(prior_files$paths)) {
    cat("\u26A0 Required NicheNet prior files were not found in:", network_dir, "\n")
    cat("  Accepted filename sets:\n")
    for (set_name in names(prior_files$missing)) {
      cat("   -", set_name, "\n")
      for (path in unname(prior_files$missing[[set_name]])) {
        cat("     ", path, "\n", sep = "")
      }
    }
    cat("  Skipping NicheNet.\n\n")
    return(invisible(NULL))
  }
  
  cat("  Using", prior_files$version, "NicheNet priors from:", network_dir, "\n")
  
  ligand_target_matrix <- readRDS(
    prior_files$paths[["ligand_target_matrix"]])
  lr_network            <- readRDS(
    prior_files$paths[["lr_network"]])
  weighted_networks     <- readRDS(
    prior_files$paths[["weighted_networks"]])
  
  # Expression matrix
  expr_mat <- .giotto_get_expression(gobj, values = "normalized", output = "matrix")
  
  # Sender / receiver cell IDs
  sender_ids   <- meta$cell_ID[meta[[celltype_col]] %in% sender_celltypes]
  receiver_ids <- meta$cell_ID[meta[[celltype_col]] == receiver_celltype]
  sender_mat <- .subset_expression_columns(expr_mat, sender_ids)
  receiver_mat <- .subset_expression_columns(expr_mat, receiver_ids)
  
  if (is.null(sender_mat) || is.null(receiver_mat)) {
    stop("No cells matched the requested sender/receiver cell types for NicheNet.")
  }
  
  sender_ids <- colnames(sender_mat)
  receiver_ids <- colnames(receiver_mat)
  
  if (length(sender_ids) == 0 || length(receiver_ids) == 0) {
    stop("No cells matched the requested sender/receiver cell types for NicheNet.")
  }
  
  sender_expr   <- .row_means_matrix_like(sender_mat)
  receiver_expr <- .row_means_matrix_like(receiver_mat)
  
  # Background and expressed genes
  background_genes  <- rownames(expr_mat)
  expressed_ligands <- intersect(
    .expressed_genes_from_matrix(sender_mat, pct = 0.10),
    lr_network$from
  )
  expressed_receptors <- intersect(
    .expressed_genes_from_matrix(receiver_mat, pct = 0.10),
    lr_network$to
  )
  
  potential_ligands <- expressed_ligands[
    expressed_ligands %in% lr_network$from[lr_network$to %in% expressed_receptors]
  ]
  
  # DE genes in receiver as target genes.
  # Prefer an explicit target gene set; otherwise fall back to a placeholder proxy.
  receiver_de <- if (!is.null(target_genes) && length(target_genes) > 0) {
    intersect(unique(as.character(target_genes)), rownames(expr_mat))
  } else {
    tryCatch({
      if (!"HAVCR1" %in% rownames(expr_mat)) {
        character(0)
      } else {
        high_cells <- receiver_ids[expr_mat["HAVCR1", receiver_ids] > 0]
        low_cells  <- setdiff(receiver_ids, high_cells)
        if (length(high_cells) < 5 || length(low_cells) < 5) character(0) else {
          fc <- rowMeans(expr_mat[, high_cells, drop = FALSE]) -
            rowMeans(expr_mat[, low_cells,  drop = FALSE])
          names(fc)[fc > 0.5]
        }
      }
    }, error = function(e) character(0))
  }
  
  if (length(receiver_de) < 10) {
    cat("\u26A0 Insufficient DE genes for NicheNet target gene set (<10).\n")
    cat("  Supply target_genes from 06_Differential_Expression.R for best results.\n\n")
    return(invisible(NULL))
  }
  
  # Ligand activity
  ligand_activities <- nichenetr::predict_ligand_activities(
    geneset               = receiver_de,
    background_expressed_genes = background_genes,
    ligand_target_matrix  = ligand_target_matrix,
    potential_ligands     = potential_ligands
  )
  
  top_ligands <- head(
    ligand_activities[order(-ligand_activities$pearson), "test_ligand"],
    top_n_ligands
  )
  
  write.csv(ligand_activities,
            file.path(out_dir, paste0(sample_id, "_ligand_activities.csv")),
            row.names = FALSE)
  
  cat("\u2713 NicheNet complete. Top ligands:\n")
  print(top_ligands)
  cat("  Results saved to:", out_dir, "\n\n")
  
  invisible(list(
    ligand_activities     = ligand_activities,
    top_ligands           = top_ligands,
    ligand_target_matrix  = ligand_target_matrix
  ))
}


# ==============================================================================
# SECTION 3 — MISTy: intercellular spatial modelling
# ==============================================================================

#' Run MISTy: multi-view intercellular spatial modelling
#'
#' Decomposes the variance of each target gene into:
#'   - intraview  (explained by the cell itself)
#'   - juxtaview  (explained by direct neighbours)
#'   - paraview   (explained by a broader spatial context)
#'
#' @param gobj          Annotated Giotto object
#' @param sample_id     Sample identifier
#' @param output_dir    Root output directory
#' @param target_genes  Genes to model as targets (default: HVGs from object)
#' @param juxta_radius  Radius for juxtaview (in coordinate units, default: 50)
#' @param para_radius   Radius for paraview (in coordinate units, default: 200)
#' @param n_cores       Parallel cores (default: 4)
#' @return MISTy results list

run_misty <- function(gobj,
                      sample_id,
                      output_dir,
                      target_genes  = NULL,
                      juxta_radius  = 50,
                      para_radius   = 200,
                      n_cores       = 4) {
  
  if (!requireNamespace("mistyR", quietly = TRUE))
    stop("mistyR not installed.\n",
         'Install: BiocManager::install("mistyR")')
  
  cat("\n--- MISTy: Intercellular spatial modelling ---\n")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "misty")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Expression matrix (cells x genes) — log-normalised
  expr_mat <- t(as.matrix(
    .giotto_get_expression(gobj, values = "normalized", output = "matrix")
  ))
  
  # Target genes: use HVGs if available, else top variable genes
  if (is.null(target_genes)) {
    target_genes <- tryCatch({
      hvgs <- getHVFInfo(gobj, var_type = "hvf")
      head(hvgs[order(-hvgs$hvf_value), "feat_ID"], 200)
    }, error = function(e) {
      rv <- apply(expr_mat, 2, var)
      names(sort(rv, decreasing = TRUE))[1:200]
    })
    cat("  Target genes: top", length(target_genes),
        "highly variable genes\n")
  }
  
  target_genes <- intersect(target_genes, colnames(expr_mat))
  feature_name_map <- data.frame(
    original = colnames(expr_mat),
    safe = make.names(colnames(expr_mat), unique = TRUE),
    stringsAsFactors = FALSE
  )
  colnames(expr_mat) <- feature_name_map$safe
  target_genes <- feature_name_map$safe[match(target_genes, feature_name_map$original)]
  target_genes <- unique(target_genes[!is.na(target_genes)])
  write.csv(
    feature_name_map,
    file.path(out_dir, paste0(sample_id, "_misty_feature_name_map.csv")),
    row.names = FALSE
  )
  
  # Spatial coordinates
  spat <- as.data.frame(.giotto_get_spatial_locations(gobj, output = "data.table"))
  rownames(spat) <- spat$cell_ID
  xy <- spat[rownames(expr_mat), c("sdimx", "sdimy")]
  
  cat("  Cells:", nrow(expr_mat), "  Targets:", length(target_genes), "\n")
  cat("  Juxtaview radius:", juxta_radius, "  Paraview radius:", para_radius, "\n\n")
  
  misty_res <- tryCatch({
    misty_views <- mistyR::create_initial_view(
      as.data.frame(expr_mat[, target_genes, drop = FALSE]),
      unique.id = paste0(sample_id, "_misty")
    )
    misty_views <- mistyR::add_juxtaview(
      misty_views,
      xy,
      neighbor.thr = juxta_radius,
      cached = FALSE,
      verbose = TRUE
    )
    misty_views <- mistyR::add_paraview(
      misty_views,
      xy,
      l = para_radius,
      zoi = juxta_radius,
      cached = FALSE,
      verbose = TRUE
    )
    
    mistyR::run_misty(
      views = misty_views,
      results.folder = out_dir,
      target.subset = target_genes,
      num.threads = n_cores
    )
  }, error = function(e) {
    cat("\u26A0 MISTy failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (!is.null(misty_res)) {
    misty_results <- mistyR::collect_results(out_dir)
    
    # Save improvement summary
    if (!is.null(misty_results$improvements)) {
      write.csv(misty_results$improvements,
                file.path(out_dir, paste0(sample_id, "_misty_improvements.csv")),
                row.names = FALSE)
    }
    
    cat("\u2713 MISTy complete. Results saved to:", out_dir, "\n\n")
    invisible(misty_results)
  } else {
    invisible(NULL)
  }
}


# ==============================================================================
# SECTION 4 — nnSVG: spatially variable gene detection
# ==============================================================================

#' Run nnSVG: spatially variable gene detection
#'
#' @param gobj          Annotated Giotto object
#' @param sample_id     Sample identifier
#' @param output_dir    Root output directory
#' @param n_top_svgs    Number of top SVGs to report (default: 100)
#' @param n_cores       Parallel cores (default: 4)
#' @return nnSVG results data frame

run_nnsvg <- function(gobj,
                      sample_id,
                      output_dir,
                      n_top_svgs = 100,
                      n_cores    = 4) {
  
  if (!requireNamespace("nnSVG", quietly = TRUE))
    stop("nnSVG not installed.\n",
         'Install: BiocManager::install("nnSVG")')
  
  cat("\n--- nnSVG: Spatially variable gene detection ---\n")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "svg")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build SpatialExperiment from Giotto
  counts_mat <- .giotto_get_expression(gobj, values = "raw", output = "matrix")
  if (inherits(counts_mat, "data.frame")) {
    counts_mat <- as.matrix(counts_mat)
  }
  if (is.null(dim(counts_mat)) || length(dim(counts_mat)) != 2) {
    stop("Giotto raw expression accessor did not return a 2D matrix-like object for nnSVG.")
  }
  if (is.null(rownames(counts_mat)) || is.null(colnames(counts_mat))) {
    stop("Giotto raw expression matrix is missing row or column names required for nnSVG.")
  }
  
  spat <- as.data.frame(.giotto_get_spatial_locations(gobj, output = "data.table"))
  rownames(spat) <- spat$cell_ID
  matched_cells <- intersect(colnames(counts_mat), rownames(spat))
  if (length(matched_cells) == 0) {
    stop("No overlapping cell IDs between raw expression and spatial coordinates for nnSVG.")
  }
  counts_mat <- counts_mat[, matched_cells, drop = FALSE]
  logcounts_mat <- .giotto_get_expression(gobj, values = "normalized", output = "matrix")
  if (inherits(logcounts_mat, "data.frame")) {
    logcounts_mat <- as.matrix(logcounts_mat)
  }
  if (is.null(dim(logcounts_mat)) || length(dim(logcounts_mat)) != 2) {
    stop("Giotto normalized expression accessor did not return a 2D matrix-like object for nnSVG.")
  }
  logcounts_mat <- logcounts_mat[, matched_cells, drop = FALSE]
  
  coords <- as.matrix(spat[matched_cells, c("sdimx", "sdimy"), drop = FALSE])
  
  # Follow the nnSVG preprocessing guidance more closely:
  # retain genes with >= 3 counts in >= 0.5% of spots, then drop all-zero spots.
  min_spots <- max(1L, ceiling(0.005 * ncol(counts_mat)))
  keep_genes <- if (inherits(counts_mat, "Matrix")) {
    Matrix::rowSums(counts_mat >= 3) >= min_spots
  } else {
    base::rowSums(counts_mat >= 3) >= min_spots
  }
  keep_spots <- if (inherits(counts_mat, "Matrix")) {
    Matrix::colSums(counts_mat) > 0
  } else {
    base::colSums(counts_mat) > 0
  }
  counts_mat <- counts_mat[keep_genes, keep_spots, drop = FALSE]
  logcounts_mat <- logcounts_mat[keep_genes, keep_spots, drop = FALSE]
  coords <- coords[keep_spots, , drop = FALSE]
  
  if (nrow(counts_mat) == 0 || ncol(counts_mat) == 0) {
    stop("nnSVG filtering removed all genes or cells.")
  }
  
  spe <- tryCatch({
    if (!requireNamespace("SpatialExperiment", quietly = TRUE))
      stop("SpatialExperiment required.\n",
           'Install: BiocManager::install("SpatialExperiment")')
    
    SpatialExperiment::SpatialExperiment(
      assays        = list(counts = counts_mat, logcounts = logcounts_mat),
      spatialCoords = coords
    )
  }, error = function(e) {
    cat("\u26A0 SpatialExperiment creation failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(spe)) return(invisible(NULL))
  
  cat("  Genes after filtering:", nrow(spe), "\n")
  cat("  Cells:", ncol(spe), "\n\n")
  
  nnsvg_res <- tryCatch(
    nnSVG::nnSVG(spe, assay_name = "logcounts", n_threads = n_cores),
    error = function(e) {
      cat("\u26A0 nnSVG failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(nnsvg_res)) {
    result_df <- as.data.frame(
      SummarizedExperiment::rowData(nnsvg_res)
    )
    result_df <- result_df[order(result_df$rank), ]
    
    write.csv(result_df,
              file.path(out_dir,
                        paste0(sample_id, "_nnSVG_results.csv")))
    
    top_svgs <- head(rownames(result_df), n_top_svgs)
    cat("\u2713 nnSVG complete.\n")
    cat("  Top 10 SVGs:", paste(head(top_svgs, 10), collapse = ", "), "\n")
    cat("  Results saved to:", out_dir, "\n\n")
    
    invisible(result_df)
  } else {
    invisible(NULL)
  }
}


# ==============================================================================
# MAIN WRAPPER
# ==============================================================================

#' Run full Layer 2/3 CCI analysis
#'
#' Calls each section in order. Each section is wrapped in tryCatch so a
#' failure in one does not abort the others.
#'
#' @param gobj              Annotated Giotto object or path
#' @param sample_id         Sample identifier
#' @param output_dir        Root output directory
#' @param celltype_col      Cell type column (auto-detected if NULL)
#' @param sender_celltypes  For NicheNet: sender cell type vector
#' @param receiver_celltype For NicheNet: single receiver cell type
#' @param run_sections      Which sections to run. Named logical vector with
#'                          names: insitucor, liana, nichenet, misty, nnsvg.
#'                          Default: all TRUE.
#' @return Named list of results from each section

.normalize_run_sections <- function(run_sections) {
  default_sections <- c(
    insitucor = TRUE,
    liana = TRUE,
    nichenet = TRUE,
    misty = TRUE,
    nnsvg = TRUE
  )
  
  if (is.null(run_sections) || length(run_sections) == 0) {
    return(default_sections)
  }
  
  if (is.list(run_sections)) {
    run_sections <- unlist(run_sections, use.names = TRUE)
  }
  
  run_section_names <- names(run_sections)
  run_sections <- as.logical(run_sections)
  names(run_sections) <- run_section_names %||% character(length(run_sections))
  
  normalized <- default_sections
  matched_names <- intersect(names(default_sections), names(run_sections))
  if (length(matched_names) > 0) {
    normalized[matched_names] <- run_sections[matched_names]
  }
  
  normalized[is.na(normalized)] <- FALSE
  normalized
}

.run_cci_section <- function(label, runner) {
  tryCatch(
    {
      result <- runner()
      list(
        result = result,
        status = if (is.null(result)) "error" else "ok",
        message = if (is.null(result)) "Section returned no result" else NULL
      )
    },
    error = function(e) {
      msg <- conditionMessage(e)
      cat("\u26A0", label, "error:", msg, "\n")
      list(result = NULL, status = "error", message = msg)
    }
  )
}

run_cci_analysis <- function(gobj,
                             sample_id,
                             output_dir,
                             celltype_col       = NULL,
                             sender_celltypes   = NULL,
                             receiver_celltype  = NULL,
                             target_genes       = NULL,
                             nichenet_network_dir = NULL,
                             run_sections       = c(insitucor = TRUE,
                                                    liana     = TRUE,
                                                    nichenet  = TRUE,
                                                    misty     = TRUE,
                                                    nnsvg     = TRUE)) {
  
  cat("\n========================================\n")
  cat("STEP 10: CCI Analysis (Layer 2 & 3)\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  run_sections <- .normalize_run_sections(run_sections)
  enabled_sections <- names(run_sections)[run_sections]
  cat(
    "Sections enabled:",
    if (length(enabled_sections) > 0) paste(enabled_sections, collapse = ", ") else "<none>",
    "\n\n"
  )
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("\u2713 Loaded\n\n")
  }
  
  results <- list()
  section_status <- setNames(rep("disabled", length(run_sections)), names(run_sections))
  section_messages <- setNames(as.list(rep(NA_character_, length(run_sections))), names(run_sections))
  
  if (!any(run_sections)) {
    cat("\u26A0 No step 10 sections were enabled. Nothing to run.\n")
    cat("\n\u2713 STEP 10 complete for", sample_id, "\n\n")
    return(invisible(results))
  }
  
  if (isTRUE(run_sections["insitucor"])) {
    section_out <- .run_cci_section(
      "InSituCor",
      function() run_insitucor(gobj, sample_id, output_dir, celltype_col)
    )
    results$insitucor <- section_out$result
    section_status["insitucor"] <- section_out$status
    section_messages[["insitucor"]] <- section_out$message
  }
  
  if (isTRUE(run_sections["liana"])) {
    section_out <- .run_cci_section(
      "LIANA",
      function() run_liana(gobj, sample_id, output_dir, celltype_col)
    )
    results$liana <- section_out$result
    section_status["liana"] <- section_out$status
    section_messages[["liana"]] <- section_out$message
  }
  
  if (isTRUE(run_sections["nichenet"])) {
    section_out <- .run_cci_section(
      "NicheNet",
      function() run_nichenet(
        gobj,
        sample_id,
        output_dir,
        celltype_col,
        sender_celltypes,
        receiver_celltype,
        target_genes = target_genes,
        network_dir = nichenet_network_dir
      )
    )
    results$nichenet <- section_out$result
    section_status["nichenet"] <- section_out$status
    section_messages[["nichenet"]] <- section_out$message
  }
  
  if (isTRUE(run_sections["misty"])) {
    section_out <- .run_cci_section(
      "MISTy",
      function() run_misty(gobj, sample_id, output_dir)
    )
    results$misty <- section_out$result
    section_status["misty"] <- section_out$status
    section_messages[["misty"]] <- section_out$message
  }
  
  if (isTRUE(run_sections["nnsvg"])) {
    section_out <- .run_cci_section(
      "nnSVG",
      function() run_nnsvg(gobj, sample_id, output_dir)
    )
    results$nnsvg <- section_out$result
    section_status["nnsvg"] <- section_out$status
    section_messages[["nnsvg"]] <- section_out$message
  }
  attr(results, "section_status") <- section_status
  attr(results, "section_messages") <- section_messages
  
  summary_outputs <- NULL
  if (.load_cci_summary_helper()) {
    summary_outputs <- tryCatch(
      create_cci_summary(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = output_dir,
        cci_results = results
      ),
      error = function(e) {
        cat("\u26A0 CCI summary error:", conditionMessage(e), "\n")
        NULL
      }
    )
  } else {
    cat("\u26A0 CCI summary helper not found; skipping summary outputs.\n")
  }
  results$summary <- summary_outputs
  
  cat("\n=== Step 10 summary ===\n")
  for (section_name in names(section_status)) {
    icon <- switch(
      section_status[[section_name]],
      ok = "\u2713",
      error = "\u26A0",
      disabled = "-",
      "-"
    )
    cat("  ", icon, " ", section_name, ": ", section_status[[section_name]], "\n", sep = "")
  }
  
  cat("\n\u2713 STEP 10 complete for", sample_id, "\n\n")
  invisible(results)
}


# Run if sourced directly ---------------------------------------------------
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) {
    script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
    script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
    bootstrap_script <- file.path(script_dir, "Helper_Scripts", "Script_Bootstrap.R")
    if (file.exists(bootstrap_script)) {
      source(bootstrap_script, local = .GlobalEnv)
      bootstrap_pipeline_environment(script_dir, load_pipeline_utils = FALSE, verbose = FALSE)
    }
    
    run_cci_analysis(
      gobj       = args[2],
      sample_id  = args[1],
      output_dir = args[3]
    )
  } else {
    stop("Usage: Rscript 10_CCI_Analysis.R <sample_id> <input_path> <output_dir>")
  }
}

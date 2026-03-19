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
#           └── svg/          Spatially variable genes (nnSVG)
# ==============================================================================


# ==============================================================================
# SECTION 0 — InSituCor (NanoString-native spatially co-expressed modules)
# ==============================================================================

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
  meta <- as.data.frame(pDataDT(gobj))
  if (is.null(celltype_col)) {
    sup_cols <- grep("celltype_.*_supervised$", names(meta), value = TRUE)
    celltype_col <- if (length(sup_cols) > 0) sup_cols[1] else "leiden_clust"
  }
  
  # Expression matrix (cells x genes)
  counts <- t(as.matrix(
    getExpression(gobj, values = "normalized", output = "matrix")
  ))
  
  # Spatial coordinates
  spat <- as.data.frame(
    getSpatialLocations(gobj, output = "data.table")
  )
  rownames(spat) <- spat$cell_ID
  
  # Cell type vector — must align with rows of counts
  ct_vec <- setNames(meta[[celltype_col]], meta$cell_ID)
  ct_vec <- ct_vec[rownames(counts)]
  
  cat("  Cells:", nrow(counts), "  Genes:", ncol(counts), "\n")
  cat("  Cell type column:", celltype_col, "\n")
  
  cor_results <- tryCatch(
    InSituCor::insitucor(
      counts         = counts,
      celltype       = ct_vec,
      xy             = spat[rownames(counts), c("sdimx", "sdimy")],
      n_cores        = n_cores,
      verbose        = TRUE
    ),
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
                file.path(out_dir, paste0(sample_id, "_insitucor_modules.csv")))
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
    stop("liana not installed.\n  Install: install.packages('liana')")
  
  cat("\n--- LIANA: Consensus ligand-receptor analysis ---\n")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "liana")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Resolve celltype column
  meta <- as.data.frame(pDataDT(gobj))
  if (is.null(celltype_col)) {
    sup_cols <- grep("celltype_.*_supervised$", names(meta), value = TRUE)
    celltype_col <- if (length(sup_cols) > 0) sup_cols[1] else "leiden_clust"
  }
  
  # Convert Giotto to Seurat for LIANA (liana expects Seurat or SCE)
  # NOTE: if SeuratDisk is available, use SaveH5Seurat / LoadH5Seurat for
  # large objects to avoid memory duplication.
  cat("  Converting Giotto to Seurat for LIANA...\n")
  seurat_obj <- tryCatch({
    # Build a minimal Seurat object from the count matrix + cell types
    if (!requireNamespace("Seurat", quietly = TRUE))
      stop("Seurat required for LIANA conversion.")
    
    counts_mat <- getExpression(gobj, values = "raw", output = "matrix")
    seu <- Seurat::CreateSeuratObject(counts = counts_mat)
    seu <- Seurat::NormalizeData(seu, verbose = FALSE)
    seu$celltype <- meta[[celltype_col]][match(colnames(seu), meta$cell_ID)]
    Seurat::Idents(seu) <- "celltype"
    seu
  }, error = function(e) {
    cat("\u26A0 Seurat conversion failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(seurat_obj)) return(invisible(NULL))
  
  if (is.null(methods)) methods <- liana::show_methods()
  
  liana_res <- tryCatch(
    liana::liana_wrap(
      seurat_obj,
      method    = methods,
      resource  = "Consensus",
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
      p <- liana::liana_dotplot(
        liana_agg,
        source_groups = NULL,   # all senders
        target_groups = NULL,   # all receivers
        ntop          = 20
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
  meta <- as.data.frame(pDataDT(gobj))
  if (is.null(celltype_col)) {
    sup_cols <- grep("celltype_.*_supervised$", names(meta), value = TRUE)
    celltype_col <- if (length(sup_cols) > 0) sup_cols[1] else "leiden_clust"
  }
  
  # NOTE: NicheNet requires pre-built networks (ligand-target, lr-network,
  # weighted networks). Download from:
  # https://zenodo.org/record/7074291
  # and set paths below.
  cat("  \u26A0 NicheNet requires pre-built prior networks.\n")
  cat("    Download from: https://zenodo.org/record/7074291\n")
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
  
  ligand_target_matrix <- readRDS(
    file.path(network_dir, "ligand_target_matrix_nsga2r_final.rds"))
  lr_network            <- readRDS(
    file.path(network_dir, "lr_network_human_21122021.rds"))
  weighted_networks     <- readRDS(
    file.path(network_dir, "weighted_networks_nsga2r_final.rds"))
  
  # Expression matrix
  expr_mat <- getExpression(gobj, values = "normalized", output = "matrix")
  
  # Sender / receiver cell IDs
  sender_ids   <- meta$cell_ID[meta[[celltype_col]] %in% sender_celltypes]
  receiver_ids <- meta$cell_ID[meta[[celltype_col]] == receiver_celltype]
  if (length(sender_ids) == 0 || length(receiver_ids) == 0) {
    stop("No cells matched the requested sender/receiver cell types for NicheNet.")
  }
  
  sender_expr   <- rowMeans(expr_mat[, sender_ids,   drop = FALSE])
  receiver_expr <- rowMeans(expr_mat[, receiver_ids, drop = FALSE])
  
  # Background and expressed genes
  background_genes  <- rownames(expr_mat)
  expressed_ligands <- nichenetr::intersect_vectors(
    nichenetr::get_expressed_genes(sender_ids,   expr_mat, pct = 0.10),
    lr_network$from
  )
  expressed_receptors <- nichenetr::intersect_vectors(
    nichenetr::get_expressed_genes(receiver_ids, expr_mat, pct = 0.10),
    lr_network$to
  )
  
  potential_ligands <- expressed_ligands[
    expressed_ligands %in% lr_network$from[lr_network$to %in% expressed_receptors]
  ]
  
  # DE genes in receiver as target genes
  # (requires a condition column; using high-vs-low expression split as proxy)
  # For real use: supply DE results from 06_Differential_Expression.R
  receiver_de <- tryCatch({
    high_cells <- receiver_ids[expr_mat["HAVCR1", receiver_ids] > 0]  # example
    low_cells  <- setdiff(receiver_ids, high_cells)
    if (length(high_cells) < 5 || length(low_cells) < 5) character(0) else {
      # Minimal DE: log2FC > 0.5 between high/low
      fc <- rowMeans(expr_mat[, high_cells, drop = FALSE]) -
        rowMeans(expr_mat[, low_cells,  drop = FALSE])
      names(fc)[fc > 0.5]
    }
  }, error = function(e) character(0))
  
  if (length(receiver_de) < 10) {
    cat("\u26A0 Insufficient DE genes for NicheNet target gene set (<10).\n")
    cat("  Supply DE results from 06_Differential_Expression.R for best results.\n\n")
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
         'Install: devtools::install_github("saezlab/mistyR")')
  
  cat("\n--- MISTy: Intercellular spatial modelling ---\n")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "misty")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Expression matrix (cells x genes) — log-normalised
  expr_mat <- t(as.matrix(
    getExpression(gobj, values = "normalized", output = "matrix")
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
  
  # Spatial coordinates
  spat <- as.data.frame(getSpatialLocations(gobj, output = "data.table"))
  rownames(spat) <- spat$cell_ID
  xy <- spat[rownames(expr_mat), c("sdimx", "sdimy")]
  
  cat("  Cells:", nrow(expr_mat), "  Targets:", length(target_genes), "\n")
  cat("  Juxtaview radius:", juxta_radius, "  Paraview radius:", para_radius, "\n\n")
  
  misty_res <- tryCatch({
    mistyR::run_misty(
      views = mistyR::create_views(
        intra_view  = mistyR::create_intra_view(
          data = as.data.frame(expr_mat[, target_genes])
        ),
        juxta_view  = mistyR::create_juxta_view(
          data   = as.data.frame(expr_mat[, target_genes]),
          positions = xy,
          l      = juxta_radius
        ),
        para_view   = mistyR::create_para_view(
          data      = as.data.frame(expr_mat),
          positions = xy,
          l         = para_radius,
          zoi       = juxta_radius
        )
      ),
      results_folder = out_dir,
      n_jobs         = n_cores
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
  counts_mat <- getExpression(gobj, values = "raw", output = "matrix")
  spat <- as.data.frame(getSpatialLocations(gobj, output = "data.table"))
  rownames(spat) <- spat$cell_ID
  
  coords <- as.matrix(spat[colnames(counts_mat), c("sdimx", "sdimy")])
  
  spe <- tryCatch({
    if (!requireNamespace("SpatialExperiment", quietly = TRUE))
      stop("SpatialExperiment required.\n",
           'Install: BiocManager::install("SpatialExperiment")')
    
    SpatialExperiment::SpatialExperiment(
      assays        = list(counts = counts_mat),
      spatialCoords = coords
    )
  }, error = function(e) {
    cat("\u26A0 SpatialExperiment creation failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(spe)) return(invisible(NULL))
  
  # Normalise within SPE
  spe <- scuttle::computeLibraryFactors(spe)
  spe <- scuttle::logNormCounts(spe)
  
  # Filter low-expressed genes
  spe <- spe[rowSums(SummarizedExperiment::assay(spe, "counts") > 0) >=
               ceiling(0.05 * ncol(spe)), ]
  
  cat("  Genes after filtering:", nrow(spe), "\n")
  cat("  Cells:", ncol(spe), "\n\n")
  
  nnsvg_res <- tryCatch(
    nnSVG::nnSVG(spe, n_threads = n_cores),
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

run_cci_analysis <- function(gobj,
                             sample_id,
                             output_dir,
                             celltype_col       = NULL,
                             sender_celltypes   = NULL,
                             receiver_celltype  = NULL,
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
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("\u2713 Loaded\n\n")
  }
  
  results <- list()
  
  if (isTRUE(run_sections["insitucor"])) {
    results$insitucor <- tryCatch(
      run_insitucor(gobj, sample_id, output_dir, celltype_col),
      error = function(e) { cat("\u26A0 InSituCor error:", conditionMessage(e), "\n"); NULL }
    )
  }
  
  if (isTRUE(run_sections["liana"])) {
    results$liana <- tryCatch(
      run_liana(gobj, sample_id, output_dir, celltype_col),
      error = function(e) { cat("\u26A0 LIANA error:", conditionMessage(e), "\n"); NULL }
    )
  }
  
  if (isTRUE(run_sections["nichenet"])) {
    results$nichenet <- tryCatch(
      run_nichenet(gobj, sample_id, output_dir, celltype_col,
                   sender_celltypes, receiver_celltype,
                   network_dir = nichenet_network_dir),
      error = function(e) { cat("\u26A0 NicheNet error:", conditionMessage(e), "\n"); NULL }
    )
  }
  
  if (isTRUE(run_sections["misty"])) {
    results$misty <- tryCatch(
      run_misty(gobj, sample_id, output_dir),
      error = function(e) { cat("\u26A0 MISTy error:", conditionMessage(e), "\n"); NULL }
    )
  }
  
  if (isTRUE(run_sections["nnsvg"])) {
    results$nnsvg <- tryCatch(
      run_nnsvg(gobj, sample_id, output_dir),
      error = function(e) { cat("\u26A0 nnSVG error:", conditionMessage(e), "\n"); NULL }
    )
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

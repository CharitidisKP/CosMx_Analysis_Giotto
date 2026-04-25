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

current_script_dir <- function() {
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

pipeline_utils <- file.path(current_script_dir(), "Helper_Scripts", "Pipeline_Utils.R")
if ((!exists("presentation_theme") || !exists("sample_plot_title") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
  source(pipeline_utils)
}


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

.normalize_named_character_list <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.list(x)) {
    return(NULL)
  }
  
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    return(NULL)
  }
  
  out <- lapply(x, function(v) unique(as.character(v[!is.na(v)])))
  names(out) <- nm
  out
}

.normalize_nichenet_spatial_option <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return("none")
  }
  if (is.logical(x)) {
    return(if (isTRUE(x[[1]])) "filtered" else "none")
  }
  
  value <- tolower(as.character(x[[1]]))
  if (value %in% c("true", "filtered", "spatial", "yes")) {
    return("filtered")
  }
  if (value %in% c("false", "none", "no", "unfiltered")) {
    return("none")
  }
  if (value %in% c("both", "all")) {
    return("both")
  }
  
  stop("Unsupported nichenet spatial option: ", x[[1]],
       ". Use FALSE, TRUE, 'filtered', 'none', or 'both'.")
}

.resolve_target_genes_for_receiver <- function(receiver_celltype,
                                               target_genes = NULL,
                                               target_genes_by_receiver = NULL) {
  target_genes_by_receiver <- .normalize_named_character_list(target_genes_by_receiver)
  if (!is.null(target_genes_by_receiver) &&
      receiver_celltype %in% names(target_genes_by_receiver)) {
    return(target_genes_by_receiver[[receiver_celltype]])
  }
  
  if (!is.null(target_genes) && length(target_genes) > 0) {
    return(unique(as.character(target_genes)))
  }
  
  NULL
}

.find_proximity_enrichment_path <- function(output_dir,
                                            sample_id,
                                            proximity_enrichment_path = NULL) {
  candidates <- unique(c(
    proximity_enrichment_path %||% "",
    file.path(output_dir, "09_Spatial_Network", "proximity",
              paste0(sample_id, "_proximity_enrichment.csv"))
  ))
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit) || !nzchar(hit)) {
    return(NULL)
  }
  hit
}

.load_proximity_enrichment <- function(output_dir,
                                       sample_id,
                                       proximity_enrichment_path = NULL) {
  prox_path <- .find_proximity_enrichment_path(
    output_dir = output_dir,
    sample_id = sample_id,
    proximity_enrichment_path = proximity_enrichment_path
  )
  if (is.null(prox_path)) {
    return(NULL)
  }
  
  prox_tbl <- tryCatch(
    data.table::fread(prox_path, data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(prox_tbl) || !"unified_int" %in% names(prox_tbl)) {
    return(NULL)
  }
  
  prox_tbl$unified_int <- as.character(prox_tbl$unified_int)
  prox_tbl
}

.pair_is_spatially_supported <- function(sender,
                                         receiver,
                                         proximity_tbl,
                                         spatial_padj_threshold = 0.05) {
  if (is.null(proximity_tbl) || nrow(proximity_tbl) == 0) {
    return(FALSE)
  }
  
  pair_labels <- c(
    paste(sender, receiver, sep = "--"),
    paste(receiver, sender, sep = "--")
  )
  pair_rows <- proximity_tbl[proximity_tbl$unified_int %in% pair_labels, , drop = FALSE]
  if (nrow(pair_rows) == 0) {
    return(FALSE)
  }
  
  if ("type_int" %in% names(pair_rows)) {
    pair_rows <- pair_rows[is.na(pair_rows$type_int) | pair_rows$type_int == "hetero", , drop = FALSE]
  }
  if (nrow(pair_rows) == 0) {
    return(FALSE)
  }
  
  enrich_ok <- if ("enrichm" %in% names(pair_rows)) pair_rows$enrichm > 0 else TRUE
  padj_ok <- if ("p.adj_higher" %in% names(pair_rows)) {
    !is.na(pair_rows$p.adj_higher) & pair_rows$p.adj_higher < spatial_padj_threshold
  } else {
    TRUE
  }
  
  any(enrich_ok & padj_ok, na.rm = TRUE)
}

.build_nichenet_comparison_table <- function(metadata,
                                             celltype_col,
                                             mode = c("single", "all_senders_to_receiver", "all_pairs"),
                                             sender_celltypes = NULL,
                                             receiver_celltype = NULL,
                                             min_cells_per_celltype = 5,
                                             include_self_pairs = FALSE) {
  mode <- match.arg(mode)
  
  celltype_counts <- sort(table(stats::na.omit(metadata[[celltype_col]])), decreasing = TRUE)
  valid_celltypes <- names(celltype_counts[celltype_counts >= min_cells_per_celltype])
  
  if (mode == "single") {
    if (is.null(sender_celltypes) || is.null(receiver_celltype)) {
      stop("NicheNet mode 'single' requires sender_celltypes and receiver_celltype.")
    }
    return(data.frame(
      sender_label = paste(unique(as.character(sender_celltypes)), collapse = " + "),
      receiver = as.character(receiver_celltype),
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(sender_list = I(list(unique(as.character(sender_celltypes))))))
  }
  
  if (mode == "all_senders_to_receiver") {
    if (is.null(receiver_celltype) || !nzchar(receiver_celltype)) {
      stop("NicheNet mode 'all_senders_to_receiver' requires receiver_celltype.")
    }
    senders <- valid_celltypes
    if (!is.null(sender_celltypes) && length(sender_celltypes) > 0) {
      senders <- intersect(senders, unique(as.character(sender_celltypes)))
    }
    if (!include_self_pairs) {
      senders <- setdiff(senders, receiver_celltype)
    }
    return(data.frame(
      sender_label = senders,
      receiver = rep(receiver_celltype, length(senders)),
      stringsAsFactors = FALSE
    ) |>
      dplyr::mutate(sender_list = as.list(sender_label)))
  }
  
  pair_grid <- expand.grid(
    sender_label = valid_celltypes,
    receiver = valid_celltypes,
    stringsAsFactors = FALSE
  )
  if (!include_self_pairs) {
    pair_grid <- pair_grid[pair_grid$sender_label != pair_grid$receiver, , drop = FALSE]
  }
  pair_grid |>
    dplyr::mutate(sender_list = as.list(sender_label))
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
#' @param expr_cache    Optional pre-extracted expression cache list with
#'                      $raw and $normalized matrices (genes x cells).
#'                      Avoids redundant dense materializations.
#' @return List with insitucor results; gobj unchanged

run_insitucor <- function(gobj,
                          sample_id,
                          output_dir,
                          celltype_col      = NULL,
                          n_cores           = 4,
                          expr_cache        = NULL,
                          sample_row        = NULL,
                          sample_sheet_path = NULL) {
  
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
  # FIX #6: Use pre-extracted cache when available to avoid repeated ~10 GB
  # dense materializations across CCI sections.
  raw_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$raw)) {
    expr_cache$raw
  } else {
    .giotto_get_expression(gobj, values = "raw", output = "matrix")
  }
  counts <- t(as.matrix(raw_mat))
  
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

    # --- InSituCor visualisations --------------------------------------
    plot_insitucor_results(
      cor_results       = cor_results,
      gobj              = gobj,
      meta              = meta,
      counts            = counts,
      celltype_col      = celltype_col,
      out_dir           = out_dir,
      sample_id         = sample_id,
      sample_row        = sample_row,
      sample_sheet_path = sample_sheet_path
    )

    cat("\u2713 InSituCor complete. Results saved to:", out_dir, "\n")
  }

  invisible(cor_results)
}

# ------------------------------------------------------------------
# InSituCor visualisation helper
# ------------------------------------------------------------------
#
# Produces the following plots + CSVs for an InSituCor result:
#   Plots
#   1. {id}_insitucor_module_sizes.png                 — module size + mean weight bars
#   2. {id}_insitucor_celltype_involvement.png         — cell-type x module heatmap
#   3. spatial_modules/{id}_insitucor_spatial_<m>.png  — per-module polygon spatial map
#   4. {id}_insitucor_module_network.png               — module co-expression network
#   5. {id}_insitucor_module_correlation_heatmap.png   — module x module correlation heatmap
#   6. {id}_insitucor_gene_weights_panel.png           — per-module gene-weight facets
#   7. {id}_insitucor_spatial_modules_combined.png     — combined spatial overlay (patchwork)
#   CSVs (in addition to {id}_insitucor_modules.csv and _module_stats.csv written by caller)
#   - {id}_insitucor_celltype_involvement.csv          — long celltype/module/involvement
#   - {id}_insitucor_cell_module_scores.csv            — wide cell_ID x module weighted scores
#   - {id}_insitucor_module_cor.csv                    — module correlation matrix
# Safe to call whenever cor_results contains `modules` (others optional).

plot_insitucor_results <- function(cor_results,
                                   gobj,
                                   meta,
                                   counts,           # cells x genes (from caller)
                                   celltype_col,
                                   out_dir,
                                   sample_id,
                                   top_n_modules     = 12,
                                   sample_row        = NULL,
                                   sample_sheet_path = NULL) {
  modules_df <- tryCatch(tibble::as_tibble(cor_results$modules),
                         error = function(e) NULL)
  if (is.null(modules_df) || nrow(modules_df) == 0 ||
      !all(c("module", "gene", "weight") %in% colnames(modules_df))) {
    cat("  \u26A0 InSituCor: no module table \u2014 skipping plots.\n")
    return(invisible(NULL))
  }

  module_stats <- modules_df %>%
    dplyr::group_by(module) %>%
    dplyr::summarise(
      n_genes      = dplyr::n(),
      mean_weight  = mean(weight, na.rm = TRUE),
      max_weight   = max(weight,  na.rm = TRUE),
      top_genes    = paste(
        gene[order(-weight)][seq_len(min(3, dplyr::n()))],
        collapse = ", "
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_genes))

  write.csv(module_stats,
            file.path(out_dir, paste0(sample_id, "_insitucor_module_stats.csv")),
            row.names = FALSE)

  top_modules <- head(module_stats, top_n_modules)

  # -- Plot 6: per-module gene-weight panel -------------------------------
  # Complements plot 1 (which only annotates the top-3 genes) by showing
  # the full weight distribution inside each top module.
  tryCatch({
    gw_df <- modules_df[modules_df$module %in% top_modules$module, , drop = FALSE]
    gw_df <- as.data.frame(gw_df, stringsAsFactors = FALSE)
    gw_df$module <- factor(gw_df$module, levels = top_modules$module)
    # Per-facet gene ordering (highest weight first)
    gw_df <- gw_df[order(gw_df$module, -gw_df$weight), ]
    gw_df$gene_ord <- factor(
      paste(gw_df$module, gw_df$gene, sep = "__"),
      levels = unique(paste(gw_df$module, gw_df$gene, sep = "__"))
    )
    n_mod <- length(levels(gw_df$module))
    ncol_facet <- min(3L, n_mod)
    p_gw <- ggplot2::ggplot(
        gw_df,
        ggplot2::aes(x = gene_ord, y = weight, fill = weight)
      ) +
      ggplot2::geom_col(color = "grey30", linewidth = 0.15) +
      ggplot2::scale_x_discrete(labels = function(x) sub(".*__", "", x)) +
      ggplot2::scale_fill_viridis_c(option = "C", name = "Weight") +
      ggplot2::facet_wrap(~ module, scales = "free", ncol = ncol_facet) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title    = sample_plot_title(sample_id,
                     "InSituCor Gene Weights per Module"),
        subtitle = paste0("All genes in each of the top ", n_mod,
                          " modules (by gene count)"),
        x = NULL, y = "Weight"
      ) +
      presentation_theme(base_size = 10) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = 7),
        strip.text  = ggplot2::element_text(face = "bold", size = 9),
        panel.spacing = grid::unit(0.6, "lines")
      )
    n_row_facet <- ceiling(n_mod / ncol_facet)
    save_presentation_plot(
      plot     = p_gw,
      filename = file.path(out_dir,
                           paste0(sample_id, "_insitucor_gene_weights_panel.png")),
      width    = max(12, 4.2 * ncol_facet),
      height   = max(8,  3.0 * n_row_facet),
      dpi      = 150
    )
    cat("  \u2713 InSituCor gene-weight panel saved\n")
  }, error = function(e) {
    cat("  \u26A0 InSituCor gene-weight panel failed:",
        conditionMessage(e), "\n")
  })

  # -- Plot 1: module sizes with top-3 genes annotated --------------------
  tryCatch({
    p1 <- ggplot2::ggplot(
      top_modules,
      ggplot2::aes(x = stats::reorder(module, n_genes),
                   y = n_genes, fill = mean_weight)
    ) +
      ggplot2::geom_col(color = "grey30", linewidth = 0.2) +
      ggplot2::geom_text(
        ggplot2::aes(label = top_genes),
        hjust = -0.05, size = 3.2, color = "grey25"
      ) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::scale_fill_viridis_c(option = "C", name = "Mean weight") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.4))) +
      ggplot2::labs(
        title    = sample_plot_title(sample_id, "InSituCor Modules"),
        subtitle = paste0("Top ", nrow(top_modules),
                          " modules by gene count; top-3 weighted genes annotated"),
        x = NULL, y = "Number of genes"
      ) +
      presentation_theme(base_size = 12)
    save_presentation_plot(
      plot     = p1,
      filename = file.path(out_dir,
                           paste0(sample_id, "_insitucor_module_sizes.png")),
      width    = 12, height = max(6, 0.5 * nrow(top_modules) + 2), dpi = 150
    )
    cat("  \u2713 InSituCor module sizes plot saved\n")
  }, error = function(e) {
    cat("  \u26A0 InSituCor module sizes failed:", conditionMessage(e), "\n")
  })

  # -- Plot 2: celltype involvement as ggplot heatmap ---------------------
  tryCatch({
    inv <- cor_results$celltypeinvolvement
    if (!is.null(inv)) {
      inv_mat <- as.matrix(inv)
      # Orient so modules are columns, celltypes are rows
      if (!is.null(colnames(inv_mat)) && all(top_modules$module %in% colnames(inv_mat))) {
        inv_mat <- inv_mat[, top_modules$module, drop = FALSE]
      } else if (!is.null(rownames(inv_mat)) &&
                 all(top_modules$module %in% rownames(inv_mat))) {
        inv_mat <- t(inv_mat[top_modules$module, , drop = FALSE])
      }
      # Order by max involvement, then take top 25 celltypes — but always
      # preserve the B-cell row if it exists in the full matrix, even when
      # its max involvement puts it outside the top 25 (otherwise B cells
      # can silently drop out of the heatmap in samples where they are rare).
      row_max  <- apply(inv_mat, 1, max, na.rm = TRUE)
      inv_mat  <- inv_mat[order(-row_max), , drop = FALSE]
      top_rows <- rownames(inv_mat)[seq_len(min(nrow(inv_mat), 25))]
      # Match any row name that normalises to "B cell" (handles "B.cell",
      # "B_cell", "B cell" variants across reference profiles).
      .norm_nm <- function(x) trimws(tolower(gsub("[._]+", " ", x)))
      bcell_row <- rownames(inv_mat)[.norm_nm(rownames(inv_mat)) == "b cell"]
      keep_rows <- unique(c(top_rows, bcell_row))
      inv_mat   <- inv_mat[keep_rows, , drop = FALSE]

      inv_long <- as.data.frame(inv_mat)
      inv_long$celltype <- rownames(inv_mat)
      inv_long <- tidyr::pivot_longer(inv_long,
                                       cols      = -celltype,
                                       names_to  = "module",
                                       values_to = "involvement")

      inv_long$celltype <- factor(inv_long$celltype, levels = rownames(inv_mat))
      inv_long$module   <- factor(inv_long$module,   levels = colnames(inv_mat))

      # Persist the numeric values used to build the heatmap so downstream
      # scripts can do module-celltype querying without re-loading the RDS.
      tryCatch({
        inv_out <- data.frame(
          celltype    = as.character(inv_long$celltype),
          module      = as.character(inv_long$module),
          involvement = as.numeric(inv_long$involvement),
          stringsAsFactors = FALSE
        )
        write.csv(
          inv_out,
          file.path(out_dir,
                    paste0(sample_id, "_insitucor_celltype_involvement.csv")),
          row.names = FALSE
        )
      }, error = function(e) {
        cat("  \u26A0 celltype_involvement CSV failed:",
            conditionMessage(e), "\n")
      })

      p2 <- ggplot2::ggplot(inv_long,
               ggplot2::aes(x = module, y = celltype, fill = involvement)) +
        ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
        ggplot2::scale_fill_viridis_c(option = "magma",
                                      name = "Involvement") +
        ggplot2::labs(
          title    = sample_plot_title(sample_id,
                       "InSituCor Cell-type Involvement"),
          subtitle = "Row = cell type, column = module",
          x = "Module", y = "Cell type"
        ) +
        presentation_theme(base_size = 11, legend_position = "right") +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y  = ggplot2::element_text(size = 9),
          axis.title.x = element_markdown_safe(margin = ggplot2::margin(t = 12)),
          axis.title.y = element_markdown_safe(margin = ggplot2::margin(r = 12)),
          plot.margin  = ggplot2::margin(t = 10, r = 24, b = 24, l = 20)
        )
      save_presentation_plot(
        plot     = p2,
        filename = file.path(out_dir,
                             paste0(sample_id,
                                    "_insitucor_celltype_involvement.png")),
        width    = max(12, ncol(inv_mat) * 0.60 + 5),
        height   = max(9,  nrow(inv_mat) * 0.34 + 3),
        dpi      = 150
      )
      cat("  \u2713 InSituCor celltype involvement plot saved\n")
    } else {
      cat("  \u26A0 InSituCor: no celltypeinvolvement slot; heatmap skipped.\n")
    }
  }, error = function(e) {
    cat("  \u26A0 InSituCor celltype involvement failed:",
        conditionMessage(e), "\n")
  })

  # -- Per-cell module scores CSV (numeric, no polygons required) --------
  # Writes wide cell_ID x module score matrix using the same log1p-normalised
  # weighted aggregate as the spatial plots below. Saved regardless of
  # whether polygons are available.
  tryCatch({
    all_mod_ids <- unique(module_stats$module)
    cell_total  <- rowSums(counts)
    cell_total[cell_total == 0] <- 1
    score_mat <- matrix(0, nrow = nrow(counts), ncol = length(all_mod_ids),
                        dimnames = list(rownames(counts), all_mod_ids))
    for (mod_id in all_mod_ids) {
      genes_in_mod   <- modules_df$gene[modules_df$module == mod_id]
      weights_in_mod <- modules_df$weight[modules_df$module == mod_id]
      keep <- genes_in_mod %in% colnames(counts)
      genes_in_mod   <- genes_in_mod[keep]
      weights_in_mod <- weights_in_mod[keep]
      if (length(genes_in_mod) == 0) next
      sub_mat  <- counts[, genes_in_mod, drop = FALSE]
      sub_norm <- log1p(sub_mat / cell_total * 1e4)
      score_mat[, mod_id] <- as.numeric(sub_norm %*% weights_in_mod)
    }
    cell_scores_df <- data.frame(
      cell_ID = rownames(score_mat),
      score_mat,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    write.csv(
      cell_scores_df,
      file.path(out_dir,
                paste0(sample_id, "_insitucor_cell_module_scores.csv")),
      row.names = FALSE
    )
    cat("  \u2713 InSituCor per-cell module scores CSV saved (",
        nrow(score_mat), " cells x ", ncol(score_mat), " modules)\n", sep = "")
  }, error = function(e) {
    cat("  \u26A0 InSituCor per-cell module scores CSV failed:",
        conditionMessage(e), "\n")
  })

  # -- Plot 3: spatial polygon maps of aggregate module expression --------
  spatial_plots <- list()
  tryCatch({
    poly_df <- .cci_extract_polygon_df(gobj)
    if (is.null(poly_df) || nrow(poly_df) == 0) {
      cat("  \u26A0 InSituCor spatial module maps skipped: polygon data unavailable.\n")
    } else {
      mod_spatial_dir <- file.path(out_dir, "spatial_modules")
      dir.create(mod_spatial_dir, showWarnings = FALSE, recursive = TRUE)
      module_list <- head(module_stats$module, min(top_n_modules, 8))
      n_saved <- 0L
      for (mod_id in module_list) {
        tryCatch({
          genes_in_mod <- modules_df$gene[modules_df$module == mod_id]
          weights_in_mod <- modules_df$weight[modules_df$module == mod_id]
          genes_in_mod <- genes_in_mod[genes_in_mod %in% colnames(counts)]
          if (length(genes_in_mod) == 0) return(invisible(NULL))
          weights_in_mod <- weights_in_mod[match(genes_in_mod,
                                                 modules_df$gene[modules_df$module == mod_id])]
          sub_mat <- counts[, genes_in_mod, drop = FALSE]
          # Weighted aggregate expression per cell (sum weight_i * norm_expr_i)
          # Using counts -> convert to log1p-normalised per cell
          cell_total <- rowSums(counts)
          cell_total[cell_total == 0] <- 1
          sub_norm   <- log1p(sub_mat / cell_total * 1e4)
          mod_score  <- as.numeric(sub_norm %*% weights_in_mod)

          score_df <- data.frame(cell_ID = rownames(counts),
                                 module_score = mod_score,
                                 stringsAsFactors = FALSE)
          plot_df <- merge(poly_df, score_df, by = "cell_ID", all.x = TRUE)
          plot_df$module_score[is.na(plot_df$module_score)] <- 0

          p3 <- ggplot2::ggplot(plot_df,
                   ggplot2::aes(x = x, y = y,
                                group = interaction(geom, part),
                                fill  = module_score)) +
            ggplot2::geom_polygon(color = NA) +
            ggplot2::scale_fill_viridis_c(option = "magma",
                                          name = "Module\nscore") +
            ggplot2::coord_equal() +
            ggplot2::labs(
              title    = sample_plot_title(sample_id,
                           paste0("InSituCor Module: ", mod_id)),
              subtitle = paste0(length(genes_in_mod),
                                " genes; weighted log1p aggregate expression")
            ) +
            presentation_theme(base_size = 11, legend_position = "right") +
            ggplot2::theme(
              panel.background = ggplot2::element_rect(fill = "grey15")
            )

          fname <- paste0(sample_id, "_insitucor_spatial_",
                          gsub("[^A-Za-z0-9]", "_", mod_id), ".png")
          save_presentation_plot(
            plot     = p3,
            filename = file.path(mod_spatial_dir, fname),
            width    = 12, height = 10, dpi = 180
          )
          spatial_plots[[as.character(mod_id)]] <- p3
          n_saved <- n_saved + 1L
        }, error = function(e) {
          cat("  \u26A0 InSituCor spatial module ", mod_id,
              " failed: ", conditionMessage(e), "\n", sep = "")
        })
      }
      cat("  \u2713 InSituCor spatial module maps saved (", n_saved,
          " plots) \u2192 ", mod_spatial_dir, "\n", sep = "")

      # Composite CART slides: emit per-sub-biopsy variants of the top
      # modules into spatial_modules/subsamples/. Iterates the shared
      # discover_composite_subsamples() helper; silently skips for
      # non-composite samples.
      sub_rows <- tryCatch(
        discover_composite_subsamples(sample_row, sample_sheet_path),
        error = function(e) NULL
      )
      if (!is.null(sub_rows) && nrow(sub_rows) > 0 && "fov" %in% names(meta)) {
        sub_dir <- file.path(mod_spatial_dir, "subsamples")
        dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
        cat("  Composite sample - rendering per-sub-biopsy spatial module variants\n")
        # Cap to top 4 modules per sub-biopsy to keep the subsample output volume bounded.
        sub_module_list <- head(module_list, min(length(module_list), 4))
        for (k in seq_len(nrow(sub_rows))) {
          sub_r  <- sub_rows[k, , drop = FALSE]
          sub_id <- as.character(sub_r$sample_id)
          fmin   <- as.integer(sub_r$fov_min)
          fmax   <- as.integer(sub_r$fov_max)
          if (anyNA(c(fmin, fmax))) next
          sub_cells <- meta$cell_ID[
            !is.na(meta$fov) & meta$fov >= fmin & meta$fov <= fmax
          ]
          if (length(sub_cells) == 0) next
          sub_poly <- poly_df[poly_df$cell_ID %in% sub_cells, , drop = FALSE]
          if (nrow(sub_poly) == 0) next
          for (mod_id in sub_module_list) {
            tryCatch({
              genes_in_mod   <- modules_df$gene[modules_df$module == mod_id]
              weights_in_mod <- modules_df$weight[modules_df$module == mod_id]
              genes_in_mod   <- genes_in_mod[genes_in_mod %in% colnames(counts)]
              if (length(genes_in_mod) == 0) return(invisible(NULL))
              weights_in_mod <- weights_in_mod[match(genes_in_mod,
                modules_df$gene[modules_df$module == mod_id])]
              sub_counts  <- counts[sub_cells, , drop = FALSE]
              sub_mat     <- sub_counts[, genes_in_mod, drop = FALSE]
              cell_total  <- rowSums(sub_counts)
              cell_total[cell_total == 0] <- 1
              sub_norm    <- log1p(sub_mat / cell_total * 1e4)
              mod_score_s <- as.numeric(sub_norm %*% weights_in_mod)
              score_df_s  <- data.frame(cell_ID = rownames(sub_counts),
                                        module_score = mod_score_s,
                                        stringsAsFactors = FALSE)
              plot_df_s   <- merge(sub_poly, score_df_s, by = "cell_ID", all.x = TRUE)
              plot_df_s$module_score[is.na(plot_df_s$module_score)] <- 0
              p_sub <- ggplot2::ggplot(plot_df_s,
                       ggplot2::aes(x = x, y = y,
                                    group = interaction(geom, part),
                                    fill  = module_score)) +
                ggplot2::geom_polygon(color = NA) +
                ggplot2::scale_fill_viridis_c(option = "magma",
                                              name = "Module\nscore") +
                ggplot2::coord_equal() +
                ggplot2::labs(
                  title    = sample_plot_title(sub_id,
                               paste0("InSituCor Module: ", mod_id)),
                  subtitle = paste0(length(genes_in_mod),
                                    " genes; weighted log1p aggregate expression")
                ) +
                presentation_theme(base_size = 11, legend_position = "right") +
                ggplot2::theme(
                  panel.background = ggplot2::element_rect(fill = "grey15")
                )
              fname_s <- paste0(sub_id, "_insitucor_spatial_",
                                gsub("[^A-Za-z0-9]", "_", mod_id), ".png")
              save_presentation_plot(
                plot     = p_sub,
                filename = file.path(sub_dir, fname_s),
                width    = 12, height = 10, dpi = 180
              )
            }, error = function(e) {
              cat("    \u26A0 sub-biopsy ", sub_id,
                  ", module ", mod_id,
                  " failed: ", conditionMessage(e), "\n", sep = "")
            })
          }
        }
      }
    }
  }, error = function(e) {
    cat("  \u26A0 InSituCor spatial module maps failed:",
        conditionMessage(e), "\n")
  })

  # -- Plot 7: combined spatial modules overlay ---------------------------
  # Stitches the per-module spatial maps (plot 3) into a single multi-panel
  # PNG. Requires `patchwork`; silently skipped otherwise.
  tryCatch({
    if (length(spatial_plots) >= 2) {
      if (!requireNamespace("patchwork", quietly = TRUE)) {
        cat("  \u26A0 InSituCor combined spatial overlay skipped: 'patchwork' not installed.\n")
      } else {
        # Patchwork layout — tighten title / subtitle spacing, constrain
        # each panel to square aspect so composite slides (long, thin
        # tissue) don't stretch the grid.
        ncol_panel <- if (length(spatial_plots) <= 4) 2L else min(3L, length(spatial_plots))
        n_row_panel <- ceiling(length(spatial_plots) / ncol_panel)
        combined <- patchwork::wrap_plots(spatial_plots, ncol = ncol_panel) +
          patchwork::plot_annotation(
            title    = sample_plot_title(sample_id,
                         "InSituCor Modules - Combined Spatial Overlay"),
            subtitle = paste0(length(spatial_plots),
                              " top modules; weighted log1p aggregate expression"),
            theme    = ggplot2::theme(
              plot.title    = ggplot2::element_text(face = "bold", size = 16,
                                                    margin = ggplot2::margin(b = 8)),
              plot.subtitle = ggplot2::element_text(size = 11,
                                                    margin = ggplot2::margin(b = 14)),
              plot.margin   = ggplot2::margin(t = 16, r = 20, b = 16, l = 20)
            )
          ) &
          ggplot2::theme(
            plot.margin = ggplot2::margin(t = 6, r = 8, b = 6, l = 8),
            plot.title  = element_markdown_safe(size = 12, margin = ggplot2::margin(b = 4))
          )
        save_presentation_plot(
          plot     = combined,
          filename = file.path(out_dir,
                               paste0(sample_id,
                                      "_insitucor_spatial_modules_combined.png")),
          width    = max(16, 6.5 * ncol_panel + 2),
          height   = max(11, 5.5 * n_row_panel + 2),
          dpi      = 150
        )
        cat("  \u2713 InSituCor combined spatial modules overlay saved\n")
      }
    }
  }, error = function(e) {
    cat("  \u26A0 InSituCor combined spatial overlay failed:",
        conditionMessage(e), "\n")
  })

  # Resolve module-correlation matrix once; used by both the CSV export,
  # the new correlation heatmap (plot 5), and the existing network (plot 4).
  mod_cor <- cor_results$modulecor %||% cor_results$moduleCor %||%
             cor_results$module_cor %||% NULL

  # -- Module correlation CSV + heatmap (plot 5) --------------------------
  tryCatch({
    if (is.null(mod_cor)) {
      cat("  \u26A0 InSituCor module correlation heatmap/CSV skipped: no modulecor slot.\n")
    } else {
      mc_full <- as.matrix(mod_cor)
      write.csv(
        mc_full,
        file.path(out_dir,
                  paste0(sample_id, "_insitucor_module_cor.csv")),
        row.names = TRUE
      )

      mc_long <- as.data.frame(mc_full)
      mc_long$module_row <- rownames(mc_full)
      mc_long <- tidyr::pivot_longer(
        mc_long,
        cols      = -module_row,
        names_to  = "module_col",
        values_to = "cor"
      )
      mc_long$module_row <- factor(mc_long$module_row, levels = rownames(mc_full))
      mc_long$module_col <- factor(mc_long$module_col, levels = colnames(mc_full))

      p_cor <- ggplot2::ggplot(
          mc_long,
          ggplot2::aes(x = module_col, y = module_row, fill = cor)
        ) +
        ggplot2::geom_tile(color = "grey90", linewidth = 0.2) +
        ggplot2::scale_fill_gradient2(
          low = "#2166AC", mid = "white", high = "#B2182B",
          midpoint = 0, limits = c(-1, 1), name = "Pearson r",
          na.value = "grey85"
        ) +
        ggplot2::coord_fixed() +
        ggplot2::labs(
          title    = sample_plot_title(sample_id,
                       "InSituCor Module Correlation Matrix"),
          subtitle = "Pairwise module-score correlations (diverging palette)",
          x = NULL, y = NULL
        ) +
        presentation_theme(base_size = 11, legend_position = "right") +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
          axis.text.y = ggplot2::element_text(size = 9)
        )
      n_mods_cor <- nrow(mc_full)
      save_presentation_plot(
        plot     = p_cor,
        filename = file.path(out_dir,
                             paste0(sample_id,
                                    "_insitucor_module_correlation_heatmap.png")),
        width    = max(9, n_mods_cor * 0.45 + 3),
        height   = max(8, n_mods_cor * 0.42 + 2),
        dpi      = 150
      )
      cat("  \u2713 InSituCor module correlation heatmap saved\n")
    }
  }, error = function(e) {
    cat("  \u26A0 InSituCor module correlation heatmap/CSV failed:",
        conditionMessage(e), "\n")
  })

  # -- Plot 1b/6b: B-cell-only module size + gene weight panels ------------
  # Ranks modules by B-cell involvement (from plot 2's involvement matrix)
  # and re-uses the same rendering style as plots 1 and 6. Silent no-op when
  # no B-cell row is present or the top B-cell involvement is zero.
  tryCatch({
    inv_all <- cor_results$celltypeinvolvement
    if (!is.null(inv_all) && nrow(modules_df) > 0) {
      inv_all_mat <- as.matrix(inv_all)
      # Orient so modules are columns, celltypes are rows.
      if (!is.null(colnames(inv_all_mat)) &&
          any(modules_df$module %in% colnames(inv_all_mat))) {
        # already module-cols; leave as is
      } else if (!is.null(rownames(inv_all_mat)) &&
                 any(modules_df$module %in% rownames(inv_all_mat))) {
        inv_all_mat <- t(inv_all_mat)
      }
      .norm_nm2 <- function(x) trimws(tolower(gsub("[._]+", " ", x)))
      bcell_row <- rownames(inv_all_mat)[.norm_nm2(rownames(inv_all_mat)) == "b cell"]
      if (length(bcell_row) > 0) {
        bcell_inv <- inv_all_mat[bcell_row[1], , drop = TRUE]
        bcell_inv <- bcell_inv[!is.na(bcell_inv)]
        bcell_inv <- sort(bcell_inv, decreasing = TRUE)
        top_bc_modules <- head(names(bcell_inv[bcell_inv > 0]),
                               min(top_n_modules, 8))
        if (length(top_bc_modules) > 0) {
          bc_stats <- module_stats[match(top_bc_modules, module_stats$module), ]
          # (i) B-cell module sizes — same style as Plot 1 ---------------
          tryCatch({
            p1b <- ggplot2::ggplot(
              bc_stats,
              ggplot2::aes(x = stats::reorder(module, n_genes),
                           y = n_genes, fill = mean_weight)
            ) +
              ggplot2::geom_col(color = "grey30", linewidth = 0.2) +
              ggplot2::geom_text(
                ggplot2::aes(label = top_genes),
                hjust = -0.05, size = 3.2, color = "grey25"
              ) +
              ggplot2::coord_flip(clip = "off") +
              ggplot2::scale_fill_viridis_c(option = "C", name = "Mean weight") +
              ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0, 0.4))
              ) +
              ggplot2::labs(
                title    = sample_plot_title(sample_id,
                             "InSituCor Modules (B-cell-ranked)"),
                subtitle = paste0("Top ", nrow(bc_stats),
                                  " modules by B-cell involvement"),
                x = NULL, y = "Number of genes"
              ) +
              presentation_theme(base_size = 12)
            save_presentation_plot(
              plot     = p1b,
              filename = file.path(out_dir,
                                   paste0(sample_id, "_insitucor_bcell_module_sizes.png")),
              width = 12,
              height = max(6, 0.5 * nrow(bc_stats) + 2),
              dpi = 150
            )
            cat("  \u2713 InSituCor B-cell module sizes plot saved\n")
          }, error = function(e) {
            cat("  \u26A0 B-cell module sizes failed:",
                conditionMessage(e), "\n")
          })

          # (ii) B-cell gene-weight panel — same style as Plot 6 ---------
          tryCatch({
            gw_b <- modules_df[modules_df$module %in% top_bc_modules, , drop = FALSE]
            gw_b <- as.data.frame(gw_b, stringsAsFactors = FALSE)
            gw_b$module <- factor(gw_b$module, levels = top_bc_modules)
            gw_b <- gw_b[order(gw_b$module, -gw_b$weight), ]
            gw_b$gene_ord <- factor(
              paste(gw_b$module, gw_b$gene, sep = "__"),
              levels = unique(paste(gw_b$module, gw_b$gene, sep = "__"))
            )
            n_mod_b   <- length(levels(gw_b$module))
            ncol_fb   <- min(3L, n_mod_b)
            p_gwb <- ggplot2::ggplot(
              gw_b,
              ggplot2::aes(x = gene_ord, y = weight, fill = weight)
            ) +
              ggplot2::geom_col(color = "grey30", linewidth = 0.15) +
              ggplot2::scale_x_discrete(labels = function(x) sub(".*__", "", x)) +
              ggplot2::scale_fill_viridis_c(option = "C", name = "Weight") +
              ggplot2::facet_wrap(~ module, scales = "free", ncol = ncol_fb) +
              ggplot2::coord_flip() +
              ggplot2::labs(
                title    = sample_plot_title(sample_id,
                             "InSituCor Gene Weights (B-cell-ranked)"),
                subtitle = paste0("All genes in the top ", n_mod_b,
                                  " B-cell-involved modules"),
                x = NULL, y = "Weight"
              ) +
              presentation_theme(base_size = 10) +
              ggplot2::theme(
                axis.text.y = ggplot2::element_text(size = 7),
                strip.text  = ggplot2::element_text(face = "bold", size = 9),
                panel.spacing = grid::unit(0.6, "lines")
              )
            n_row_fb <- ceiling(n_mod_b / ncol_fb)
            save_presentation_plot(
              plot     = p_gwb,
              filename = file.path(out_dir,
                                   paste0(sample_id, "_insitucor_bcell_gene_weights_panel.png")),
              width = max(12, 4.2 * ncol_fb),
              height = max(8, 3.0 * n_row_fb),
              dpi = 150
            )
            cat("  \u2713 InSituCor B-cell gene-weight panel saved\n")
          }, error = function(e) {
            cat("  \u26A0 B-cell gene-weight panel failed:",
                conditionMessage(e), "\n")
          })
        } else {
          cat("  \u2139 No modules with positive B-cell involvement; skipping B-cell panels\n")
        }
      } else {
        cat("  \u2139 No B-cell row in celltype involvement; skipping B-cell panels\n")
      }
    }
  }, error = function(e) {
    cat("  \u26A0 InSituCor B-cell module/gene panels failed:",
        conditionMessage(e), "\n")
  })

  # -- Plot 4: module co-expression network -------------------------------
  tryCatch({
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      cat("  \u26A0 InSituCor module network skipped: 'ggraph' not installed.\n")
    } else {
      # Edges: between modules that share genes (weighted by shared count),
      # or — if an insitucor $condcor / $modulecor matrix exists — by its
      # off-diagonal correlations.
      mod_net_edges <- NULL
      if (!is.null(mod_cor)) {
        mc <- as.matrix(mod_cor)
        diag(mc) <- NA
        mc[lower.tri(mc)] <- NA
        idx <- which(!is.na(mc) & abs(mc) > 0.15, arr.ind = TRUE)
        if (nrow(idx) > 0) {
          mod_net_edges <- data.frame(
            from   = rownames(mc)[idx[, 1]],
            to     = colnames(mc)[idx[, 2]],
            weight = mc[idx],
            stringsAsFactors = FALSE
          )
        }
      }
      if (is.null(mod_net_edges)) {
        # Fallback: shared-gene edges
        mods <- split(modules_df$gene, modules_df$module)
        mods <- mods[top_modules$module]
        pairs <- utils::combn(names(mods), 2, simplify = FALSE)
        mod_net_edges <- do.call(rbind, lapply(pairs, function(pr) {
          shared <- length(intersect(mods[[pr[1]]], mods[[pr[2]]]))
          if (shared > 0) {
            data.frame(from = pr[1], to = pr[2], weight = shared,
                       stringsAsFactors = FALSE)
          } else NULL
        }))
      }

      if (is.null(mod_net_edges) || nrow(mod_net_edges) == 0) {
        cat("  \u26A0 InSituCor module network skipped: no inter-module edges.\n")
      } else {
        g <- igraph::graph_from_data_frame(mod_net_edges, directed = FALSE,
               vertices = data.frame(
                 name = top_modules$module,
                 n_genes = top_modules$n_genes,
                 stringsAsFactors = FALSE
               ))
        p4 <- ggraph::ggraph(g, layout = "fr") +
          ggraph::geom_edge_link(
            ggplot2::aes(edge_width = abs(weight), edge_alpha = abs(weight)),
            color = "steelblue"
          ) +
          ggraph::geom_node_point(
            ggplot2::aes(size = n_genes),
            color = "firebrick", alpha = 0.85
          ) +
          ggraph::geom_node_label(
            ggplot2::aes(label = name),
            size = 3.5, fill = "white",
            label.padding = grid::unit(0.15, "lines")
          ) +
          ggraph::scale_edge_width_continuous(range = c(0.3, 2.5),
                                              name = "Edge weight") +
          ggraph::scale_edge_alpha_continuous(range = c(0.25, 0.85),
                                              guide = "none") +
          ggplot2::scale_size_continuous(range = c(3, 12),
                                         name = "Module\nsize (genes)") +
          ggplot2::labs(
            title    = sample_plot_title(sample_id,
                         "InSituCor Module Co-expression Network"),
            subtitle = "Nodes = modules (sized by gene count); edges weighted by module correlation or shared genes"
          ) +
          ggraph::theme_graph(background = "white", base_size = 11) +
          ggplot2::theme(
            plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 10, color = "grey20")
          )
        save_presentation_plot(
          plot     = p4,
          filename = file.path(out_dir,
                               paste0(sample_id,
                                      "_insitucor_module_network.png")),
          width    = 12, height = 10, dpi = 150
        )
        cat("  \u2713 InSituCor module co-expression network saved\n")
      }
    }
  }, error = function(e) {
    cat("  \u26A0 InSituCor module network failed:",
        conditionMessage(e), "\n")
  })

  invisible(NULL)
}


# ==============================================================================
# SECTION 1 helper — Extended LIANA visualizations
# ==============================================================================

# Extract cell polygon vertices from a Giotto object as a flat data.frame
# with columns: geom, part, x, y, hole, cell_ID. Returns NULL when polygon
# data are unavailable. Mirrors the helper in 07_Annotation.R so this file
# can render polygon-backed spatial plots without cross-script sourcing.
.cci_extract_polygon_df <- function(gobj) {
  poly_sv <- tryCatch({
    p <- GiottoClass::getPolygonInfo(gobject = gobj, polygon_name = "cell",
                                     return_giottoPolygon = FALSE)
    if (inherits(p, "giottoPolygon")) p@spatVector else p
  }, error = function(e) {
    tryCatch({
      gp <- GiottoClass::getPolygonInfo(gobject = gobj, polygon_name = "cell")
      if (inherits(gp, "giottoPolygon")) gp@spatVector else NULL
    }, error = function(e2) NULL)
  })
  if (is.null(poly_sv)) {
    poly_sv <- tryCatch(gobj@polygon$cell@spatVector, error = function(e) NULL)
  }
  if (is.null(poly_sv)) return(NULL)

  poly_attr <- tryCatch(as.data.frame(poly_sv), error = function(e) NULL)
  if (is.null(poly_attr) || nrow(poly_attr) == 0) return(NULL)
  poly_coords <- tryCatch(terra::geom(poly_sv, df = TRUE), error = function(e) NULL)
  if (is.null(poly_coords) || nrow(poly_coords) == 0) return(NULL)

  id_col <- intersect(c("poly_ID", "cell_ID", "id"), names(poly_attr))
  if (length(id_col) == 0) return(NULL)
  poly_coords$cell_ID <- poly_attr[[id_col[1]]][poly_coords$geom]
  poly_coords
}

# ==============================================================================
# Shared B-cell polygon-overlay helpers used by run_nichenet / run_misty /
# run_nnsvg. Reuse .cci_extract_polygon_df() above. Each plot outlines B cells
# (thick blue border) on top of a base polygon map coloured by `value`.
# ==============================================================================

.resolve_bcell_ids <- function(gobj, celltype_col, focus_celltype = "^B cell$") {
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (is.null(celltype_col) || !celltype_col %in% names(meta)) return(character(0))
  hit <- grepl(focus_celltype, as.character(meta[[celltype_col]]),
               ignore.case = TRUE)
  unique(as.character(meta$cell_ID[hit]))
}

.bcell_polygon_base <- function(poly_df, bcell_ids) {
  poly_df$poly_group <- paste(poly_df$cell_ID, poly_df$geom, poly_df$part,
                              sep = "_")
  poly_df$.is_bcell  <- poly_df$cell_ID %in% bcell_ids
  poly_df
}

#' Single-panel B-cell polygon overlay.
#' @param values named numeric vector keyed by cell_ID (or gene-expression row).
.plot_bcell_polygon_panel <- function(poly_df, values, title = "", subtitle = "",
                                      fill_label = "value",
                                      highlight_colour = "mediumblue") {
  poly_df$.val <- unname(values[poly_df$cell_ID])
  bcell_df <- poly_df[poly_df$.is_bcell, , drop = FALSE]
  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly_df,
      mapping = ggplot2::aes(x = x, y = y, group = poly_group, fill = .val),
      colour = "grey40", linewidth = 0.05
    ) +
    ggplot2::geom_polygon(
      data = bcell_df,
      mapping = ggplot2::aes(x = x, y = y, group = poly_group),
      fill = NA, colour = highlight_colour, linewidth = 0.35
    ) +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90",
                                  name = fill_label) +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 9,  colour = "grey30"),
      legend.key.height = grid::unit(10, "pt")
    )
}

#' Multi-panel grid: one panel per gene, B cells outlined everywhere.
#' `genes` is a character vector of gene names present in `expr_mat` rows.
.plot_bcell_polygon_grid <- function(poly_df, expr_mat, genes, bcell_ids,
                                     outfile, ncol_grid = 3, width = 12,
                                     height = 10, title = NULL,
                                     highlight_colour = "mediumblue") {
  genes <- intersect(genes, rownames(expr_mat))
  if (length(genes) == 0) return(FALSE)
  panels <- lapply(genes, function(g) {
    vals <- expr_mat[g, ]
    .plot_bcell_polygon_panel(
      poly_df = poly_df, values = vals,
      title = g, fill_label = "expr",
      highlight_colour = highlight_colour
    )
  })
  if (!requireNamespace("patchwork", quietly = TRUE)) return(FALSE)
  combined <- patchwork::wrap_plots(panels, ncol = ncol_grid)
  if (!is.null(title)) {
    combined <- combined + patchwork::plot_annotation(title = title)
  }
  ggplot2::ggsave(outfile, combined, width = width, height = height, dpi = 150)
  TRUE
}

#' Colour each cell polygon by a categorical celltype; outline B cells.
.plot_bcell_polygon_category <- function(poly_df, ct_vec, focus_labels,
                                         title, outfile, width = 9, height = 8,
                                         highlight_colour = "mediumblue") {
  poly_df$.ct <- unname(ct_vec[poly_df$cell_ID])
  poly_df$.ct <- ifelse(poly_df$.ct %in% focus_labels, poly_df$.ct, "Other")
  bcell_df <- poly_df[poly_df$.is_bcell, , drop = FALSE]
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = poly_df,
      mapping = ggplot2::aes(x = x, y = y, group = poly_group, fill = .ct),
      colour = "grey40", linewidth = 0.05
    ) +
    ggplot2::geom_polygon(
      data = bcell_df,
      mapping = ggplot2::aes(x = x, y = y, group = poly_group),
      fill = NA, colour = highlight_colour, linewidth = 0.35
    ) +
    ggplot2::scale_fill_brewer(palette = "Set1", na.value = "grey90",
                               name = "sender") +
    ggplot2::coord_equal() +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 10)
  ggplot2::ggsave(outfile, p, width = width, height = height, dpi = 150)
  TRUE
}

plot_liana_extended <- function(liana_agg,
                                meta,
                                gobj,
                                expr_cache,
                                out_dir,
                                sample_id,
                                celltype_col,
                                top_n             = 20,
                                focus_celltype    = NULL,
                                sample_row        = NULL,
                                sample_sheet_path = NULL) {

  if (is.null(liana_agg) || nrow(liana_agg) == 0) return(invisible(NULL))

  source_col   <- .first_existing_column(liana_agg, c("source", "sender", "Sender"))
  target_col   <- .first_existing_column(liana_agg, c("target", "receiver", "Receiver"))
  ligand_col   <- .first_existing_column(liana_agg, c("ligand_complex",   "ligand.complex",   "ligand",   "Ligand"))
  receptor_col <- .first_existing_column(liana_agg, c("receptor_complex", "receptor.complex", "receptor", "Receptor"))
  rank_col     <- .first_existing_column(liana_agg, c("aggregate_rank", "rank", "mean_rank"))

  if (is.null(source_col) || is.null(target_col) ||
      is.null(ligand_col) || is.null(receptor_col) || is.null(rank_col)) {
    cat("\u26A0 plot_liana_extended: could not resolve required LIANA columns; skipping extended plots.\n")
    return(invisible(NULL))
  }

  agg <- liana_agg
  agg$source           <- agg[[source_col]]
  agg$target           <- agg[[target_col]]
  agg$ligand_complex   <- agg[[ligand_col]]
  agg$receptor_complex <- agg[[receptor_col]]
  agg$aggregate_rank   <- agg[[rank_col]]

  count_df <- tryCatch(
    dplyr::count(agg, source, target, name = "n_interactions"),
    error = function(e) NULL
  )

  # Dedicated subfolder for all focus-cell-type (B cell) specific plots
  bcell_dir <- file.path(out_dir, "B_cell_specific")
  if (!is.null(focus_celltype) && length(focus_celltype) > 0) {
    dir.create(bcell_dir, showWarnings = FALSE, recursive = TRUE)
  }

  # -- Plot 1: LR ranking bar chart ------------------------------------------
  tryCatch({
    agg_sorted <- agg[order(agg$aggregate_rank, na.last = TRUE), ]
    top_lr     <- head(agg_sorted, top_n)

    # Ensure every source cell type gets at least one entry even if not in top_n
    represented <- unique(top_lr$source)
    extra <- do.call(rbind, lapply(setdiff(unique(agg$source), represented), function(ct) {
      sub <- agg_sorted[agg_sorted$source == ct, ]
      if (nrow(sub) > 0) sub[1, ] else NULL
    }))
    if (!is.null(extra) && nrow(extra) > 0) top_lr <- rbind(top_lr, extra)
    top_lr$interaction_label  <- paste(top_lr$ligand_complex, "\u2192", top_lr$receptor_complex)
    top_lr$neg_log10_rank     <- -log10(pmax(top_lr$aggregate_rank, 1e-10))
    top_lr$is_focus           <- !is.null(focus_celltype) & top_lr$source %in% focus_celltype

    n_rows    <- nrow(top_lr)
    n_senders <- length(unique(top_lr$source))
    p1_w      <- max(14, ceiling(n_rows * 0.32) + 6)
    p1_h      <- max(8,  ceiling(n_rows * 0.28))

    # Uniform divider thickness: map the outline colour via aes so focus
    # rows stand out by colour (black vs grey) but every bar has the same
    # linewidth — the old two-layer overlay gave focus bars a visibly
    # thicker outline, which TASKS.md flagged as inconsistent.
    p1 <- ggplot2::ggplot(top_lr,
             ggplot2::aes(
               x     = reorder(interaction_label, neg_log10_rank),
               y     = neg_log10_rank,
               fill  = source,
               color = is_focus
             )) +
      ggplot2::geom_col(linewidth = 0.25) +
      ggplot2::scale_color_manual(
        values = c(`FALSE` = "grey40", `TRUE` = "black"),
        guide  = "none"
      ) +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::labs(
        title    = sample_plot_title(sample_id, "Top Ligand-Receptor Interactions"),
        subtitle = paste0("Top ", top_n,
                          " interactions + best per sender (\u2212log\u2081\u2080; taller = more significant)"),
        x        = "Interaction (Ligand \u2192 Receptor)",
        y        = "\u2212log\u2081\u2080(Aggregate Rank)",
        fill     = "Sender"
      ) +
      presentation_theme(base_size = 11, legend_position = "right") +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
      ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, max(80, max(nchar(top_lr$interaction_label)) * 5)))

    save_presentation_plot(
      plot     = p1,
      filename = file.path(out_dir, paste0(sample_id, "_liana_lr_ranking.png")),
      width    = p1_w,
      height   = p1_h,
      dpi      = 150
    )
    cat("\u2713 LIANA LR ranking plot saved\n")
  }, error = function(e) {
    cat("\u26A0 LIANA LR ranking plot failed:", conditionMessage(e), "\n")
  })

  # -- Plot 1b: LIANA dotplot (main + B-cell-specific) -----------------------
  # Capped at 8 targets for legibility; x-axis rotated 45° via
  # presentation_theme(x_angle = 45) so rotation survives liana_dotplot's
  # internal theme.
  render_liana_dotplot <- function(source_groups, target_groups,
                                   title_text, subtitle_text, out_path,
                                   base_size   = 11,
                                   axis_size   = 10,
                                   strip_size  = 11,
                                   width_per_target = 2.8) {
    if (!requireNamespace("liana", quietly = TRUE)) {
      cat("\u26A0 LIANA dotplot skipped: 'liana' not installed.\n")
      return(invisible(FALSE))
    }
    specificity_col <- .first_existing_column(
      liana_agg,
      c("natmi.edge_specificity", "connectome.edge_specificity",
        "cellphonedb.pvalue", "lr_means")
    )
    magnitude_col <- .first_existing_column(
      liana_agg,
      c("sca.LRscore", "connectome.scaled_weight",
        "logfc.logfc_comb", "natmi.edge_average", "lr_means")
    )
    if (is.null(specificity_col) || is.null(magnitude_col)) {
      stop("Could not identify LIANA plot columns for the selected methods.")
    }

    # Pre-filter to requested source/target groups so liana_dotplot()'s
    # internal ntop cut doesn't leave the facet variable empty (which
    # aborts ggplot2's facet_wrap with "Faceting variables must have at
    # least one value").
    agg_sub <- liana_agg
    if (!is.null(source_groups) && length(source_groups) > 0) {
      agg_sub <- agg_sub[agg_sub$source %in% source_groups, , drop = FALSE]
    }
    if (!is.null(target_groups) && length(target_groups) > 0) {
      agg_sub <- agg_sub[agg_sub$target %in% target_groups, , drop = FALSE]
    }
    if (nrow(agg_sub) == 0) {
      cat("\u26A0 LIANA dotplot skipped: no rows after source/target filter.\n")
      return(invisible(FALSE))
    }

    n_targets <- length(target_groups %||% unique(agg_sub$target))
    dp_width  <- max(22, n_targets * width_per_target + 10)
    dp_height <- max(16, 20 * 0.6 + 8)

    p <- liana::liana_dotplot(
      agg_sub,
      source_groups = source_groups,
      target_groups = target_groups,
      ntop          = 20,
      specificity   = specificity_col,
      magnitude     = magnitude_col
    ) +
      ggplot2::labs(title = title_text, subtitle = subtitle_text) +
      presentation_theme(base_size = base_size, legend_position = "right",
                         x_angle = 45) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(size = axis_size),
        axis.text.y  = ggplot2::element_text(size = axis_size),
        # Facet strip (cell-type) titles: no rotation, centered
        strip.text.x = ggplot2::element_text(size = strip_size,
                                             angle = 0,
                                             hjust = 0.5, vjust = 0.5,
                                             face = "bold",
                                             margin = ggplot2::margin(t = 4, b = 4)),
        strip.text.y = ggplot2::element_text(size = strip_size),
        strip.clip   = "off",
        plot.margin  = ggplot2::margin(t = 40, r = 20, b = 40, l = 20)
      )
    save_presentation_plot(plot = p, filename = out_path,
                           width = dp_width, height = dp_height, dpi = 150)
    TRUE
  }

  # Main dotplot: top 8 senders x top 8 receivers (by interaction count).
  tryCatch({
    src_col_dp <- .first_existing_column(liana_agg, c("source", "sender", "Sender"))
    tgt_col_dp <- .first_existing_column(liana_agg, c("target", "receiver", "Receiver"))
    top_n_groups <- 8L

    dp_sources <- names(sort(table(liana_agg[[src_col_dp]]), decreasing = TRUE))[
      seq_len(min(top_n_groups, length(unique(liana_agg[[src_col_dp]]))))]
    dp_targets <- names(sort(table(liana_agg[[tgt_col_dp]]), decreasing = TRUE))[
      seq_len(min(top_n_groups, length(unique(liana_agg[[tgt_col_dp]]))))]

    render_liana_dotplot(
      source_groups = dp_sources,
      target_groups = dp_targets,
      title_text    = sample_plot_title(sample_id,
                        "LIANA Ligand-Receptor Interactions"),
      subtitle_text = paste0("Top ", top_n_groups,
                             " senders \u00D7 top ", top_n_groups,
                             " receivers by interaction count"),
      out_path      = file.path(out_dir,
                       paste0(sample_id, "_liana_dotplot.png"))
    )
    cat("\u2713 LIANA dotplot saved\n")
  }, error = function(e) {
    cat("\u26A0 LIANA dotplot failed:", conditionMessage(e), "\n")
  })

  # B-cell dotplot: focus_celltype as sender, top 10 receivers.
  tryCatch({
    if (is.null(focus_celltype) || length(focus_celltype) == 0) {
      cat("\u26A0 LIANA B-cell dotplot skipped: no focus_celltype resolved.\n")
    } else {
      focus_pairs <- count_df[count_df$source %in% focus_celltype, ]
      if (nrow(focus_pairs) == 0) {
        cat("\u26A0 LIANA B-cell dotplot skipped: no ",
            paste(focus_celltype, collapse = "/"),
            " sender edges.\n", sep = "")
      } else {
        bcell_targets <- head(
          focus_pairs[order(-focus_pairs$n_interactions), "target"], 10
        )
        render_liana_dotplot(
          source_groups = focus_celltype,
          target_groups = bcell_targets,
          title_text    = sample_plot_title(sample_id,
                            "B-cell LIANA Ligand-Receptor Interactions"),
          subtitle_text = paste0(paste(focus_celltype, collapse = "/"),
                                 " as sender \u2014 top 20 L-R pairs, top ",
                                 length(bcell_targets),
                                 " receivers by interaction count"),
          out_path      = file.path(bcell_dir,
                            paste0(sample_id, "_liana_dotplot_bcell.png")),
          # Larger type for the B-cell version \u2014 fewer panels = more room
          base_size        = 14,
          axis_size        = 13,
          strip_size       = 15,
          width_per_target = 3.6
        )
        cat("\u2713 LIANA B-cell dotplot saved \u2192", bcell_dir, "\n")
      }
    }
  }, error = function(e) {
    cat("\u26A0 LIANA B-cell dotplot failed:", conditionMessage(e), "\n")
  })

  # -- Plot 2: CCI heatmap ---------------------------------------------------
  tryCatch({
    if (is.null(count_df) || nrow(count_df) == 0) stop("No interaction count data.")

    n_ct      <- length(unique(c(count_df$source, count_df$target)))
    tile_size <- max(10, ceiling(n_ct * 0.45))

    p2 <- ggplot2::ggplot(count_df,
             ggplot2::aes(x = source, y = target, fill = n_interactions)) +
      ggplot2::geom_tile(color = "grey90", linewidth = 0.3) +
      ggplot2::geom_text(ggplot2::aes(label = n_interactions), size = 2.5) +
      ggplot2::scale_fill_viridis_c(option = "plasma", name = "Interactions (n)") +
      ggplot2::labs(
        title    = sample_plot_title(sample_id, "Cell-Cell Interaction Heatmap"),
        subtitle = "Number of significant L-R pairs per sender-receiver combination",
        x        = "Sender",
        y        = "Receiver"
      ) +
      presentation_theme(base_size = 10, legend_position = "right") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = ggplot2::element_text(size = 7)
      )

    save_presentation_plot(
      plot     = p2,
      filename = file.path(out_dir, paste0(sample_id, "_liana_cci_heatmap.png")),
      width    = tile_size,
      height   = tile_size,
      dpi      = 150
    )
    cat("\u2713 LIANA CCI heatmap saved\n")
  }, error = function(e) {
    cat("\u26A0 LIANA CCI heatmap failed:", conditionMessage(e), "\n")
    message("  Detail: ", conditionMessage(e))
  })

  # -- Plot 3: CCI network graph (requires ggraph) ---------------------------
  # Shared renderer for circular CCI network plots. Takes a source-target
  # edge table (net_df) and renders it with the standard ggraph style;
  # highlights any vertices in `highlight_celltype` in firebrick.
  render_cci_network <- function(net_df, highlight_celltype,
                                 title_text, subtitle_text, out_path) {
    net_df <- net_df[!duplicated(paste(net_df$source, net_df$target)), ]
    net_df$edge_color <- net_df$source

    out_deg <- stats::setNames(
      tapply(net_df$n_interactions, net_df$source, sum, default = 0),
      unique(net_df$source)
    )
    in_deg <- stats::setNames(
      tapply(net_df$n_interactions, net_df$target, sum, default = 0),
      unique(net_df$target)
    )
    all_nodes <- unique(c(net_df$source, net_df$target))
    node_size <- vapply(all_nodes, function(n)
      (out_deg[n] %||% 0) + (in_deg[n] %||% 0), numeric(1))

    g <- igraph::graph_from_data_frame(net_df, directed = TRUE)
    if (!is.null(highlight_celltype)) {
      missing_focus <- setdiff(highlight_celltype, igraph::V(g)$name)
      if (length(missing_focus) > 0) {
        g <- igraph::add_vertices(g, length(missing_focus), name = missing_focus)
        node_size <- c(node_size,
                       stats::setNames(rep(2, length(missing_focus)), missing_focus))
      }
    }
    igraph::V(g)$total_degree <- node_size[igraph::V(g)$name]
    igraph::V(g)$is_focus     <- igraph::V(g)$name %in% (highlight_celltype %||% character(0))

    p <- ggraph::ggraph(g, layout = "circle") +
      ggraph::geom_edge_arc(
        ggplot2::aes(edge_width = n_interactions, edge_colour = edge_color,
                     edge_alpha = n_interactions),
        arrow   = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
        end_cap = ggraph::circle(5, "mm")
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(size = total_degree),
        color = "grey30", alpha = 0.6
      ) +
      ggraph::geom_node_point(
        data  = function(x) x[x$is_focus, ],
        ggplot2::aes(size = total_degree),
        color = "firebrick", alpha = 0.9
      ) +
      ggraph::geom_node_label(
        ggplot2::aes(label = name,
                     fontface = ifelse(is_focus, "bold", "plain")),
        size = 2.5, fill = "white", label.padding = grid::unit(0.15, "lines")
      ) +
      ggraph::scale_edge_width_continuous(range = c(0.4, 2.5),
                                          name = "Interactions (n)") +
      ggraph::scale_edge_alpha_continuous(range = c(0.25, 0.85), guide = "none") +
      ggraph::scale_edge_colour_discrete(guide = "none") +
      ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
      ggplot2::labs(title = title_text, subtitle = subtitle_text) +
      ggraph::theme_graph(background = "white", base_size = 11) +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13,
                                                margin = ggplot2::margin(b = 6)),
        plot.subtitle   = ggplot2::element_text(hjust = 0.5, size = 10, color = "grey20"),
        legend.position = "right",
        legend.title    = ggplot2::element_text(face = "bold", size = 10),
        legend.text     = ggplot2::element_text(size = 9),
        plot.margin     = ggplot2::margin(60, 80, 60, 80)
      ) +
      ggplot2::coord_cartesian(clip = "off")

    save_presentation_plot(plot = p, filename = out_path,
                           width = 14, height = 14, dpi = 150)
  }

  # Plot 3a: pure top-25 edges (no focus union) -----------------------------
  tryCatch({
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      cat("\u26A0 LIANA CCI network skipped: 'ggraph' is not installed.\n")
    } else if (is.null(count_df) || nrow(count_df) == 0) {
      stop("No interaction count data.")
    } else {
      net_top25 <- head(count_df[order(-count_df$n_interactions), ], 25)
      render_cci_network(
        net_df             = net_top25,
        highlight_celltype = focus_celltype,
        title_text         = sample_plot_title(sample_id, "CCI Network"),
        subtitle_text      = paste0("Top 25 sender-receiver pairs; ",
                                    "edge width proportional to N significant L-R pairs"),
        out_path           = file.path(out_dir,
                               paste0(sample_id, "_liana_cci_network.png"))
      )
      cat("\u2713 LIANA CCI network saved (top-25)\n")
    }
  }, error = function(e) {
    cat("\u26A0 LIANA CCI network failed:", conditionMessage(e), "\n")
    message("  Detail: ", conditionMessage(e))
  })

  # Plot 3a-withBcell: top 25 + every focus-cell-type edge ------------------
  tryCatch({
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      # already warned in 3a
    } else if (is.null(focus_celltype) || length(focus_celltype) == 0) {
      cat("\u26A0 LIANA CCI network (top25 + B-cell) skipped: no focus_celltype.\n")
    } else if (is.null(count_df) || nrow(count_df) == 0) {
      stop("No interaction count data.")
    } else {
      net_top25  <- head(count_df[order(-count_df$n_interactions), ], 25)
      focus_edges <- count_df[count_df$source %in% focus_celltype |
                                count_df$target %in% focus_celltype, ]
      net_merged <- rbind(net_top25, focus_edges)
      render_cci_network(
        net_df             = net_merged,
        highlight_celltype = focus_celltype,
        title_text         = sample_plot_title(sample_id,
                               "CCI Network (top 25 + B-cell edges)"),
        subtitle_text      = paste0("Top 25 sender-receiver pairs plus every edge ",
                                    "involving ", paste(focus_celltype, collapse = "/"),
                                    "; edge width \u221D N significant L-R pairs"),
        out_path           = file.path(out_dir,
                               paste0(sample_id,
                                      "_liana_cci_network_top25_with_bcell.png"))
      )
      cat("\u2713 LIANA CCI network saved (top-25 + B-cell)\n")
    }
  }, error = function(e) {
    cat("\u26A0 LIANA CCI network (top25 + B-cell) failed:",
        conditionMessage(e), "\n")
  })

  # -- Plot 3b: B-cell-focused CCI network (requires ggraph) ----------------
  tryCatch({
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      cat("\u26A0 LIANA B-cell CCI network skipped: 'ggraph' is not installed.\n")
    } else if (is.null(focus_celltype) || length(focus_celltype) == 0) {
      cat("\u26A0 LIANA B-cell CCI network skipped: no focus_celltype resolved.\n")
    } else if (is.null(count_df) || nrow(count_df) == 0) {
      stop("No interaction count data.")
    } else {
      bcell_df <- count_df[count_df$source %in% focus_celltype |
                             count_df$target %in% focus_celltype, ]
      if (nrow(bcell_df) == 0) {
        cat("\u26A0 LIANA B-cell CCI network skipped: no ",
            paste(focus_celltype, collapse = "/"),
            " edges in count_df.\n", sep = "")
      } else {
        bcell_df$edge_color <- bcell_df$source

        out_deg <- stats::setNames(
          tapply(bcell_df$n_interactions, bcell_df$source, sum, default = 0),
          unique(bcell_df$source)
        )
        in_deg  <- stats::setNames(
          tapply(bcell_df$n_interactions, bcell_df$target, sum, default = 0),
          unique(bcell_df$target)
        )
        all_nodes <- unique(c(bcell_df$source, bcell_df$target))
        node_size <- vapply(all_nodes, function(n)
          (out_deg[n] %||% 0) + (in_deg[n] %||% 0), numeric(1))

        gB <- igraph::graph_from_data_frame(bcell_df, directed = TRUE)
        igraph::V(gB)$total_degree <- node_size[igraph::V(gB)$name]
        igraph::V(gB)$is_focus     <- igraph::V(gB)$name %in% focus_celltype

        p3b <- ggraph::ggraph(gB, layout = "circle") +
          ggraph::geom_edge_arc(
            ggplot2::aes(edge_width = n_interactions,
                         edge_colour = edge_color,
                         edge_alpha  = n_interactions),
            arrow   = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
            end_cap = ggraph::circle(5, "mm")
          ) +
          ggraph::geom_node_point(
            ggplot2::aes(size = total_degree),
            color = "grey30", alpha = 0.6
          ) +
          ggraph::geom_node_point(
            data  = function(x) x[x$is_focus, ],
            ggplot2::aes(size = total_degree),
            color = "firebrick", alpha = 0.9
          ) +
          ggraph::geom_node_label(
            ggplot2::aes(label = name,
                         fontface = ifelse(is_focus, "bold", "plain")),
            size = 2.5, fill = "white",
            label.padding = grid::unit(0.15, "lines")
          ) +
          ggraph::scale_edge_width_continuous(range = c(0.4, 2.5),
                                              name  = "Interactions (n)") +
          ggraph::scale_edge_alpha_continuous(range = c(0.25, 0.85),
                                              guide = "none") +
          ggraph::scale_edge_colour_discrete(guide = "none") +
          ggplot2::scale_size_continuous(range = c(2, 8), guide = "none") +
          ggplot2::labs(
            title    = sample_plot_title(sample_id, "B-cell CCI Network"),
            subtitle = paste0("All edges where ",
                              paste(focus_celltype, collapse = "/"),
                              " is sender or receiver; edge width \u221D ",
                              "N significant L-R pairs")
          ) +
          ggraph::theme_graph(background = "white", base_size = 11) +
          ggplot2::theme(
            plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                    size = 13,
                                                    margin = ggplot2::margin(b = 6)),
            plot.subtitle   = ggplot2::element_text(hjust = 0.5, size = 10,
                                                    color = "grey20"),
            legend.position = "right",
            legend.title    = ggplot2::element_text(face = "bold", size = 10),
            legend.text     = ggplot2::element_text(size = 9),
            plot.margin     = ggplot2::margin(60, 80, 60, 80)
          ) +
          ggplot2::coord_cartesian(clip = "off")

        save_presentation_plot(
          plot     = p3b,
          filename = file.path(bcell_dir,
                               paste0(sample_id, "_liana_cci_network_bcell.png")),
          width    = 14,
          height   = 14,
          dpi      = 150
        )
        cat("\u2713 LIANA B-cell CCI network saved \u2192", bcell_dir, "\n")
      }
    }
  }, error = function(e) {
    cat("\u26A0 LIANA B-cell CCI network failed:", conditionMessage(e), "\n")
    message("  Detail: ", conditionMessage(e))
  })

  # -- Plot 3c: Top-N B-cell CCI network (filtered to most significant L-R) --
  # Uses aggregate_rank (LIANA consensus rank across methods) to keep the
  # top N most significant B-cell-involved L-R pairs, then aggregates to
  # source-target counts for the graph. More readable than Plot 3b when the
  # focus cell type interacts with many partners.
  tryCatch({
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      cat("\u26A0 LIANA B-cell top-N CCI network skipped: 'ggraph' not installed.\n")
    } else if (is.null(focus_celltype) || length(focus_celltype) == 0) {
      cat("\u26A0 LIANA B-cell top-N CCI network skipped: no focus_celltype resolved.\n")
    } else {
      top_n_bcell_pairs <- 50L
      bcell_agg <- agg[agg$source %in% focus_celltype |
                         agg$target %in% focus_celltype, ]
      if (nrow(bcell_agg) == 0) {
        cat("\u26A0 LIANA B-cell top-N CCI network skipped: no B-cell L-R pairs.\n")
      } else {
        bcell_agg_top <- head(
          bcell_agg[order(bcell_agg$aggregate_rank, na.last = TRUE), ],
          top_n_bcell_pairs
        )
        bcell_top_df <- dplyr::count(
          bcell_agg_top, source, target, name = "n_interactions"
        )
        bcell_top_df$edge_color <- bcell_top_df$source

        out_deg <- stats::setNames(
          tapply(bcell_top_df$n_interactions, bcell_top_df$source,
                 sum, default = 0),
          unique(bcell_top_df$source)
        )
        in_deg <- stats::setNames(
          tapply(bcell_top_df$n_interactions, bcell_top_df$target,
                 sum, default = 0),
          unique(bcell_top_df$target)
        )
        all_nodes <- unique(c(bcell_top_df$source, bcell_top_df$target))
        node_size <- vapply(all_nodes, function(n)
          (out_deg[n] %||% 0) + (in_deg[n] %||% 0), numeric(1))

        gC <- igraph::graph_from_data_frame(bcell_top_df, directed = TRUE)
        igraph::V(gC)$total_degree <- node_size[igraph::V(gC)$name]
        igraph::V(gC)$is_focus     <- igraph::V(gC)$name %in% focus_celltype

        p3c <- ggraph::ggraph(gC, layout = "circle") +
          ggraph::geom_edge_arc(
            ggplot2::aes(edge_width = n_interactions,
                         edge_colour = edge_color,
                         edge_alpha  = n_interactions),
            arrow   = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
            end_cap = ggraph::circle(5, "mm")
          ) +
          ggraph::geom_node_point(
            ggplot2::aes(size = total_degree),
            color = "grey30", alpha = 0.6
          ) +
          ggraph::geom_node_point(
            data  = function(x) x[x$is_focus, ],
            ggplot2::aes(size = total_degree),
            color = "firebrick", alpha = 0.9
          ) +
          ggraph::geom_node_label(
            ggplot2::aes(label = name,
                         fontface = ifelse(is_focus, "bold", "plain")),
            size = 3, fill = "white",
            label.padding = grid::unit(0.15, "lines")
          ) +
          ggraph::scale_edge_width_continuous(range = c(0.6, 3),
                                              name  = "Top-N L-R pairs") +
          ggraph::scale_edge_alpha_continuous(range = c(0.4, 0.9),
                                              guide = "none") +
          ggraph::scale_edge_colour_discrete(guide = "none") +
          ggplot2::scale_size_continuous(range = c(3, 10), guide = "none") +
          ggplot2::labs(
            title    = sample_plot_title(sample_id,
                         "B-cell CCI Network (top interactions)"),
            subtitle = paste0("Top ", top_n_bcell_pairs,
                              " ", paste(focus_celltype, collapse = "/"),
                              " L-R pairs by LIANA aggregate_rank; ",
                              "edge width = N pairs per partner")
          ) +
          ggraph::theme_graph(background = "white", base_size = 11) +
          ggplot2::theme(
            plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold",
                                                    size = 13,
                                                    margin = ggplot2::margin(b = 6)),
            plot.subtitle   = ggplot2::element_text(hjust = 0.5, size = 10,
                                                    color = "grey20"),
            legend.position = "right",
            legend.title    = ggplot2::element_text(face = "bold", size = 10),
            legend.text     = ggplot2::element_text(size = 9),
            plot.margin     = ggplot2::margin(60, 80, 60, 80)
          ) +
          ggplot2::coord_cartesian(clip = "off")

        save_presentation_plot(
          plot     = p3c,
          filename = file.path(bcell_dir,
                               paste0(sample_id, "_liana_cci_network_bcell_topN.png")),
          width    = 14,
          height   = 14,
          dpi      = 150
        )
        cat("\u2713 LIANA B-cell top-N CCI network saved \u2192", bcell_dir, "\n")
      }
    }
  }, error = function(e) {
    cat("\u26A0 LIANA B-cell top-N CCI network failed:", conditionMessage(e), "\n")
    message("  Detail: ", conditionMessage(e))
  })

  # -- Plot 4: Chord diagram (requires circlize) -----------------------------
  tryCatch({
    if (!requireNamespace("circlize", quietly = TRUE)) {
      cat("\u26A0 LIANA chord diagram skipped: 'circlize' is not installed.\n")
    } else if (is.null(count_df) || nrow(count_df) == 0) {
      stop("No interaction count data.")
    } else {
      # If focus_celltype set: weight by top-N B-cell L-R importance; else top 20 pairs
      if (!is.null(focus_celltype)) {
        # Weight chord thickness by top-N aggregate_rank importance —
        # count the number of "significant" L-R pairs (top 100 by rank)
        # per source-target, so highly-ranked interactions drive width.
        top_n_chord <- 100L
        bcell_agg_full <- agg[agg$source %in% focus_celltype |
                                agg$target %in% focus_celltype, ]
        bcell_agg_full <- bcell_agg_full[order(bcell_agg_full$aggregate_rank,
                                               na.last = TRUE), ]
        top_bcell <- head(bcell_agg_full, top_n_chord)
        chord_df <- if (nrow(top_bcell) > 0) {
          dplyr::count(top_bcell, source, target, name = "n_interactions")
        } else {
          count_df[count_df$source %in% focus_celltype |
                     count_df$target %in% focus_celltype, ]
        }
        if (nrow(chord_df) == 0) chord_df <- count_df  # fallback
        chord_title <- paste0(sample_id, " \u2014 B-cell CCI Chord (top ",
                              min(top_n_chord, nrow(bcell_agg_full)),
                              " L-R pairs)")
        chord_path <- file.path(bcell_dir,
                                paste0(sample_id, "_liana_chord_bcell.png"))
      } else {
        top_ct_pairs <- head(count_df[order(-count_df$n_interactions), ], 20)
        keep_cts     <- unique(c(top_ct_pairs$source, top_ct_pairs$target))
        chord_df     <- count_df[count_df$source %in% keep_cts &
                                    count_df$target %in% keep_cts, ]
        chord_title  <- paste0(sample_id, " \u2014 CCI Chord Diagram")
        chord_path   <- file.path(out_dir, paste0(sample_id, "_liana_chord.png"))
      }
      mat <- stats::xtabs(n_interactions ~ source + target, data = chord_df)
      # Wider canvas + reserve outer ring for sector labels so long
      # cell-type names ("MNP.b.non.classical.monocyte.derived") aren't
      # cut off at the plot margin.
      grDevices::png(chord_path, width = 3200, height = 3200, res = 220)
      circlize::circos.clear()
      circlize::circos.par(
        gap.degree   = 3,
        start.degree = 90,
        track.margin = c(0.01, 0.01),
        canvas.xlim  = c(-1.4, 1.4),
        canvas.ylim  = c(-1.4, 1.4)
      )
      chord_args <- list(
        x                 = mat,
        transparency      = 0.35,
        annotationTrack   = "grid",
        preAllocateTracks = list(track.height = circlize::mm_h(14))
      )
      if (!is.null(focus_celltype)) {
        focus_in_mat <- intersect(focus_celltype, unique(c(rownames(mat), colnames(mat))))
        if (length(focus_in_mat) > 0) {
          grid_col <- stats::setNames(
            rep("grey60", length(unique(c(rownames(mat), colnames(mat))))),
            unique(c(rownames(mat), colnames(mat)))
          )
          grid_col[focus_in_mat] <- "firebrick"
          chord_args$grid.col <- grid_col
        }
      }
      do.call(circlize::chordDiagram, chord_args)
      circlize::circos.trackPlotRegion(
        track.index = 1,
        panel.fun   = function(x, y) {
          xlim <- circlize::get.cell.meta.data("xlim")
          ylim <- circlize::get.cell.meta.data("ylim")
          sector.name <- circlize::get.cell.meta.data("sector.index")
          circlize::circos.text(mean(xlim), ylim[1] + 0.4,
                                sector.name, facing = "clockwise",
                                niceFacing = TRUE,
                                adj = c(0, 0.5), cex = 0.8)
        },
        bg.border = NA
      )
      graphics::title(main = chord_title, cex.main = 1.4, line = -1)
      grDevices::dev.off()
      circlize::circos.clear()
      cat("\u2713 LIANA chord diagram saved \u2192", chord_path, "\n")
    }
  }, error = function(e) {
    tryCatch(grDevices::dev.off(), error = function(e2) NULL)
    cat("\u26A0 LIANA chord diagram failed:", conditionMessage(e), "\n")
  })

  # -- Plot 5: Spatial LR expression maps (polygon-based) ---------------------
  # Renders real cell shapes via getPolygonInfo(). Emits two directories:
  #   out_dir/spatial_lr/                         — default top-10 LR pairs
  #   bcell_dir/spatial_lr/ (if focus set)        — EVERY B-cell L-R pair
  tryCatch({
    poly_df <- .cci_extract_polygon_df(gobj)
    spat    <- tryCatch(
      as.data.frame(.giotto_get_spatial_locations(gobj, output = "data.table")),
      error = function(e) NULL
    )
    use_polygons <- !is.null(poly_df) && nrow(poly_df) > 0
    if (!use_polygons) {
      if (is.null(spat) || !all(c("cell_ID", "sdimx", "sdimy") %in% names(spat))) {
        stop("Neither polygon nor centroid spatial data available.")
      }
      cat("  (polygon info unavailable; falling back to centroid points)\n")
    }

    expr_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$normalized)) {
      expr_cache$normalized
    } else {
      tryCatch(
        .giotto_get_expression(gobj, values = "normalized", output = "matrix"),
        error = function(e) NULL
      )
    }
    if (is.null(expr_mat)) stop("Expression matrix unavailable.")

    extract_expr_vec <- function(gene, cell_ids) {
      ids <- intersect(as.character(cell_ids), colnames(expr_mat))
      if (length(ids) == 0 || !gene %in% rownames(expr_mat)) return(NULL)
      vals <- expr_mat[gene, ids, drop = TRUE]
      data.frame(cell_ID = ids, expression = as.numeric(vals),
                 stringsAsFactors = FALSE)
    }

    # Renders one ligand/receptor pair to the given directory. Reusable
    # closure so the default (top-10) and B-cell-all loops share logic.
    render_lr_pair <- function(pair, out_path,
                               cell_ids_subset = NULL,
                               title_sample_id = NULL) {
      ligand_gene   <- strsplit(as.character(pair$ligand_complex),   "_")[[1]][1]
      receptor_gene <- strsplit(as.character(pair$receptor_complex), "_")[[1]][1]
      sender_ct     <- as.character(pair$source)
      receiver_ct   <- as.character(pair$target)

      # Composite sub-biopsy support: when cell_ids_subset is supplied,
      # restrict sender / receiver / background layers to those cells.
      if (!is.null(cell_ids_subset) && length(cell_ids_subset) > 0) {
        subset_meta <- meta[meta$cell_ID %in% cell_ids_subset, , drop = FALSE]
        sender_ids   <- subset_meta$cell_ID[!is.na(subset_meta[[celltype_col]]) &
                                              subset_meta[[celltype_col]] == sender_ct]
        receiver_ids <- subset_meta$cell_ID[!is.na(subset_meta[[celltype_col]]) &
                                              subset_meta[[celltype_col]] == receiver_ct]
      } else {
        sender_ids   <- meta$cell_ID[!is.na(meta[[celltype_col]]) &
                                       meta[[celltype_col]] == sender_ct]
        receiver_ids <- meta$cell_ID[!is.na(meta[[celltype_col]]) &
                                       meta[[celltype_col]] == receiver_ct]
      }
      lig_df <- extract_expr_vec(ligand_gene,   sender_ids)
      rec_df <- extract_expr_vec(receptor_gene, receiver_ids)
      if (is.null(lig_df) && is.null(rec_df)) return(FALSE)

      panels <- Filter(Negate(is.null), list(
        if (!is.null(lig_df)) {
          lig_df$panel <- paste0("Ligand: ", ligand_gene, "\n(", sender_ct, ")"); lig_df
        },
        if (!is.null(rec_df)) {
          rec_df$panel <- paste0("Receptor: ", receptor_gene, "\n(", receiver_ct, ")"); rec_df
        }
      ))
      panel_expr <- do.call(rbind, panels)

      # When subsetting, restrict background polygons + centroid frame too.
      .poly_df_use <- if (!is.null(cell_ids_subset)) {
        poly_df[poly_df$cell_ID %in% cell_ids_subset, , drop = FALSE]
      } else poly_df
      .spat_use <- if (!is.null(cell_ids_subset) && !is.null(spat)) {
        spat[spat$cell_ID %in% cell_ids_subset, , drop = FALSE]
      } else spat

      if (use_polygons) {
        # Every polygon vertex gets the matching cell's expression so
        # geom_polygon fills by expression. Cells outside the current
        # panel's sender/receiver get fill = NA (drawn light grey).
        bg_df <- .poly_df_use  # all cell outlines (possibly subset)
        panel_df <- do.call(rbind, lapply(unique(panel_expr$panel),
          function(pnl) {
            sub <- panel_expr[panel_expr$panel == pnl, ]
            m   <- merge(.poly_df_use, sub[, c("cell_ID", "expression")],
                         by = "cell_ID", all.x = FALSE)
            m$panel <- pnl
            m
          }))
        if (nrow(panel_df) == 0) return(FALSE)

        # Duplicate bg outlines per panel so facet_wrap shows them in each
        bg_df_panels <- do.call(rbind, lapply(unique(panel_expr$panel),
          function(pnl) { b <- bg_df; b$panel <- pnl; b }))

        p5 <- ggplot2::ggplot() +
          ggplot2::geom_polygon(
            data = bg_df_panels,
            ggplot2::aes(x = x, y = y,
                         group = interaction(panel, geom, part)),
            fill = "grey88", color = NA
          ) +
          ggplot2::geom_polygon(
            data = panel_df,
            ggplot2::aes(x = x, y = y,
                         group = interaction(panel, geom, part),
                         fill  = expression),
            color = NA
          ) +
          ggplot2::facet_wrap(~panel) +
          ggplot2::scale_fill_viridis_c(option = "magma", name = "Expression") +
          ggplot2::coord_equal()
      } else {
        panel_df <- merge(panel_expr,
                          .spat_use[, c("cell_ID", "sdimx", "sdimy")],
                          by = "cell_ID")
        p5 <- ggplot2::ggplot() +
          ggplot2::geom_point(data = .spat_use,
            ggplot2::aes(x = sdimx, y = sdimy),
            color = "grey80", size = 0.15, alpha = 0.35) +
          ggplot2::geom_point(data = panel_df,
            ggplot2::aes(x = sdimx, y = sdimy, color = expression),
            size = 0.8, alpha = 0.9) +
          ggplot2::facet_wrap(~panel) +
          ggplot2::scale_color_viridis_c(option = "magma", name = "Expression") +
          ggplot2::coord_equal()
      }

      .title_id <- if (!is.null(title_sample_id) && nzchar(title_sample_id)) title_sample_id else sample_id
      p5 <- p5 +
        ggplot2::labs(
          title    = sample_plot_title(.title_id,
                       paste0("Spatial Expression: ",
                              ligand_gene, " \u2192 ", receptor_gene)),
          subtitle = paste0(sender_ct, " (ligand) \u2192 ",
                            receiver_ct, " (receptor)"),
          x = "x", y = "y"
        ) +
        presentation_theme(base_size = 12, legend_position = "right") +
        ggplot2::theme(panel.background =
                         ggplot2::element_rect(fill = "grey15"))

      save_presentation_plot(
        plot = p5, filename = out_path,
        width = 16, height = 8, dpi = 200
      )
      TRUE
    }

    # Default selection: top 8 focus + top 2 others (or top 10 overall)
    if (!is.null(focus_celltype)) {
      focus_rows <- agg[agg$source %in% focus_celltype |
                          agg$target %in% focus_celltype, ]
      focus_rows <- focus_rows[order(focus_rows$aggregate_rank,
                                     na.last = TRUE), ]
      other_rows <- agg[!(agg$source %in% focus_celltype |
                            agg$target %in% focus_celltype), ]
      other_rows <- other_rows[order(other_rows$aggregate_rank,
                                     na.last = TRUE), ]
      pairs_to_plot <- rbind(head(focus_rows, 8), head(other_rows, 2))
    } else {
      pairs_to_plot <- head(agg[order(agg$aggregate_rank, na.last = TRUE), ], 10)
    }
    if (nrow(pairs_to_plot) == 0) stop("No interactions to map.")

    sp_lr_dir <- file.path(out_dir, "spatial_lr")
    dir.create(sp_lr_dir, showWarnings = FALSE, recursive = TRUE)

    n_saved <- 0L
    for (i in seq_len(nrow(pairs_to_plot))) {
      pair <- pairs_to_plot[i, ]
      tryCatch({
        ligand_gene   <- strsplit(as.character(pair$ligand_complex),   "_")[[1]][1]
        receptor_gene <- strsplit(as.character(pair$receptor_complex), "_")[[1]][1]
        fname <- paste0(sample_id, "_spatial_lr_",
                        sprintf("%02d", i), "_",
                        gsub("[^A-Za-z0-9]", "_", ligand_gene), "_",
                        gsub("[^A-Za-z0-9]", "_", receptor_gene), ".png")
        ok <- render_lr_pair(pair, file.path(sp_lr_dir, fname))
        if (isTRUE(ok)) n_saved <- n_saved + 1L
      }, error = function(e) {
        cat("  \u26A0 Spatial LR pair", i, "failed:", conditionMessage(e), "\n")
      })
    }
    cat("\u2713 LIANA spatial LR maps saved (", n_saved, "plots) \u2192",
        sp_lr_dir, "\n")

    # Composite CART slides: split the main spatial_lr panels per sub-
    # biopsy into spatial_lr/subsamples/{sub_id}/. Per-FOV variants are
    # intentionally NOT emitted here (would produce hundreds of files
    # per sample); the per-FOV asks in TASKS.md item (spatial_lr) remain
    # a deferred enhancement pending user confirmation.
    sub_rows_liana <- tryCatch(
      discover_composite_subsamples(sample_row, sample_sheet_path),
      error = function(e) NULL
    )
    if (!is.null(sub_rows_liana) && nrow(sub_rows_liana) > 0 &&
        "fov" %in% names(meta)) {
      sub_root <- file.path(sp_lr_dir, "subsamples")
      dir.create(sub_root, recursive = TRUE, showWarnings = FALSE)
      # Cap to top 5 pairs per sub-biopsy to keep file counts bounded.
      sub_pairs <- head(pairs_to_plot, min(nrow(pairs_to_plot), 5))
      for (k in seq_len(nrow(sub_rows_liana))) {
        sub_r  <- sub_rows_liana[k, , drop = FALSE]
        sub_id <- as.character(sub_r$sample_id)
        fmin   <- as.integer(sub_r$fov_min)
        fmax   <- as.integer(sub_r$fov_max)
        if (anyNA(c(fmin, fmax))) next
        sub_cells <- meta$cell_ID[
          !is.na(meta$fov) & meta$fov >= fmin & meta$fov <= fmax
        ]
        if (length(sub_cells) == 0) next
        sub_dir_k <- file.path(sub_root, sub_id)
        dir.create(sub_dir_k, recursive = TRUE, showWarnings = FALSE)
        for (i in seq_len(nrow(sub_pairs))) {
          pair <- sub_pairs[i, ]
          tryCatch({
            lig_g <- strsplit(as.character(pair$ligand_complex),   "_")[[1]][1]
            rec_g <- strsplit(as.character(pair$receptor_complex), "_")[[1]][1]
            fn <- paste0(sub_id, "_spatial_lr_",
                         sprintf("%02d", i), "_",
                         gsub("[^A-Za-z0-9]", "_", lig_g), "_",
                         gsub("[^A-Za-z0-9]", "_", rec_g), ".png")
            render_lr_pair(pair, file.path(sub_dir_k, fn),
                           cell_ids_subset = sub_cells,
                           title_sample_id = sub_id)
          }, error = function(e) {
            cat("    \u26A0 ", sub_id, " pair ", i,
                " failed: ", conditionMessage(e), "\n", sep = "")
          })
        }
      }
      cat("\u2713 LIANA spatial LR per-sub-biopsy variants saved \u2192 ",
          sub_root, "\n", sep = "")
    }

    # Additionally: EVERY B-cell L-R pair in its own B_cell_specific subfolder
    if (!is.null(focus_celltype) && length(focus_celltype) > 0) {
      bcell_all <- agg[agg$source %in% focus_celltype |
                         agg$target %in% focus_celltype, ]
      bcell_all <- bcell_all[order(bcell_all$aggregate_rank, na.last = TRUE), ]
      # Cap at 60 to keep runtime reasonable; ranked by aggregate_rank.
      bcell_cap <- min(60L, nrow(bcell_all))
      bcell_all <- head(bcell_all, bcell_cap)

      if (nrow(bcell_all) > 0) {
        bcell_sp_dir <- file.path(bcell_dir, "spatial_lr")
        dir.create(bcell_sp_dir, showWarnings = FALSE, recursive = TRUE)
        n_bcell_saved <- 0L
        for (i in seq_len(nrow(bcell_all))) {
          pair <- bcell_all[i, ]
          tryCatch({
            ligand_gene   <- strsplit(as.character(pair$ligand_complex),   "_")[[1]][1]
            receptor_gene <- strsplit(as.character(pair$receptor_complex), "_")[[1]][1]
            fname <- paste0(sample_id, "_spatial_lr_bcell_",
                            sprintf("%02d", i), "_",
                            gsub("[^A-Za-z0-9]", "_", as.character(pair$source)), "_to_",
                            gsub("[^A-Za-z0-9]", "_", as.character(pair$target)), "_",
                            gsub("[^A-Za-z0-9]", "_", ligand_gene), "_",
                            gsub("[^A-Za-z0-9]", "_", receptor_gene), ".png")
            ok <- render_lr_pair(pair, file.path(bcell_sp_dir, fname))
            if (isTRUE(ok)) n_bcell_saved <- n_bcell_saved + 1L
          }, error = function(e) {
            cat("  \u26A0 B-cell spatial LR pair", i,
                "failed:", conditionMessage(e), "\n")
          })
        }
        cat("\u2713 LIANA B-cell spatial LR maps saved (",
            n_bcell_saved, "plots) \u2192", bcell_sp_dir, "\n")

        # Composite sub-biopsy split for B-cell LR pairs.
        sub_rows_bcell <- tryCatch(
          discover_composite_subsamples(sample_row, sample_sheet_path),
          error = function(e) NULL
        )
        if (!is.null(sub_rows_bcell) && nrow(sub_rows_bcell) > 0 &&
            "fov" %in% names(meta)) {
          sub_root_b <- file.path(bcell_sp_dir, "subsamples")
          dir.create(sub_root_b, recursive = TRUE, showWarnings = FALSE)
          sub_pairs_b <- head(bcell_all, min(nrow(bcell_all), 5))
          for (k in seq_len(nrow(sub_rows_bcell))) {
            sub_r  <- sub_rows_bcell[k, , drop = FALSE]
            sub_id <- as.character(sub_r$sample_id)
            fmin   <- as.integer(sub_r$fov_min)
            fmax   <- as.integer(sub_r$fov_max)
            if (anyNA(c(fmin, fmax))) next
            sub_cells <- meta$cell_ID[
              !is.na(meta$fov) & meta$fov >= fmin & meta$fov <= fmax
            ]
            if (length(sub_cells) == 0) next
            sub_dir_k <- file.path(sub_root_b, sub_id)
            dir.create(sub_dir_k, recursive = TRUE, showWarnings = FALSE)
            for (i in seq_len(nrow(sub_pairs_b))) {
              pair <- sub_pairs_b[i, ]
              tryCatch({
                lig_g <- strsplit(as.character(pair$ligand_complex),   "_")[[1]][1]
                rec_g <- strsplit(as.character(pair$receptor_complex), "_")[[1]][1]
                fn <- paste0(sub_id, "_spatial_lr_bcell_",
                             sprintf("%02d", i), "_",
                             gsub("[^A-Za-z0-9]", "_", as.character(pair$source)), "_to_",
                             gsub("[^A-Za-z0-9]", "_", as.character(pair$target)), "_",
                             gsub("[^A-Za-z0-9]", "_", lig_g), "_",
                             gsub("[^A-Za-z0-9]", "_", rec_g), ".png")
                render_lr_pair(pair, file.path(sub_dir_k, fn),
                               cell_ids_subset = sub_cells,
                               title_sample_id = sub_id)
              }, error = function(e) {
                cat("    \u26A0 ", sub_id, " b-cell pair ", i,
                    " failed: ", conditionMessage(e), "\n", sep = "")
              })
            }
          }
          cat("\u2713 LIANA B-cell spatial LR per-sub-biopsy variants saved \u2192 ",
              sub_root_b, "\n", sep = "")
        }
      }
    }
  }, error = function(e) {
    cat("\u26A0 LIANA spatial LR map failed:", conditionMessage(e), "\n")
  })

  # -- Plot 6a: UMAP CCI overview -----------------------------------------------
  tryCatch({
    if (is.null(count_df) || nrow(count_df) == 0) stop("No interaction count data.")

    # Robust UMAP extraction: try several stored names, handle both
    # data.table-with-cell_ID-column and matrix-with-rownames outputs,
    # and keep two dimension columns regardless of how Giotto returns them.
    .extract_umap_df <- function(gobj) {
      candidate_names <- c("umap", "umap.harmony", "umap_harmony",
                           "umap.pca", "umap_pca")
      for (nm in candidate_names) {
        raw <- tryCatch(
          .giotto_get_dim_reduction(gobj, reduction_method = "umap",
                                    name = nm, output = "data.table"),
          error = function(e) NULL
        )
        if (is.null(raw)) next
        df <- as.data.frame(raw, stringsAsFactors = FALSE)
        # Case 1: matrix-like with cell IDs in rownames
        if (!"cell_ID" %in% names(df) && !is.null(rownames(df)) &&
            length(rownames(df)) == nrow(df)) {
          rn <- rownames(df)
          # Only promote rownames if they look like cell IDs (non-numeric)
          if (any(suppressWarnings(is.na(as.numeric(rn))))) {
            df$cell_ID <- rn
            df <- df[, c("cell_ID",
                         setdiff(names(df), "cell_ID")), drop = FALSE]
          }
        }
        dim_cols <- setdiff(names(df), "cell_ID")
        if (length(dim_cols) >= 2 && "cell_ID" %in% names(df)) {
          df <- df[, c("cell_ID", dim_cols[1], dim_cols[2]), drop = FALSE]
          names(df)[2:3] <- c("umap_1", "umap_2")
          return(df)
        }
      }
      NULL
    }

    umap_coords <- .extract_umap_df(gobj)
    if (is.null(umap_coords)) {
      stop("UMAP dim reduction unavailable under names: umap, umap.harmony, umap.pca.")
    }

    umap_df <- merge(
      umap_coords,
      meta[, c("cell_ID", celltype_col)],
      by = "cell_ID"
    )
    names(umap_df)[names(umap_df) == celltype_col] <- "celltype"
    umap_df$cell_ID <- as.character(umap_df$cell_ID)

    # Per-cell L-R activity score:
    # For each of the top N L-R pairs we add the cell's ligand expression
    # when its celltype is the pair's source, and its receptor expression
    # when its celltype is the pair's target. The final score is the sum,
    # log1p-scaled. Cells from non-participating cell types stay at 0.
    top_n_lr <- 20L
    top_lr_agg <- head(agg[order(agg$aggregate_rank, na.last = TRUE), ],
                      top_n_lr)

    expr_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$normalized)) {
      expr_cache$normalized
    } else {
      tryCatch(
        .giotto_get_expression(gobj, values = "normalized", output = "matrix"),
        error = function(e) NULL
      )
    }
    if (is.null(expr_mat)) stop("Expression matrix unavailable for L-R activity.")

    cell_ids_all <- as.character(umap_df$cell_ID)
    valid_cells  <- intersect(cell_ids_all, colnames(expr_mat))
    ct_by_cell   <- stats::setNames(umap_df$celltype, umap_df$cell_ID)
    lr_score     <- stats::setNames(rep(0, length(cell_ids_all)), cell_ids_all)

    for (i in seq_len(nrow(top_lr_agg))) {
      pair          <- top_lr_agg[i, ]
      ligand_gene   <- strsplit(as.character(pair$ligand_complex),   "_")[[1]][1]
      receptor_gene <- strsplit(as.character(pair$receptor_complex), "_")[[1]][1]
      src_ct        <- as.character(pair$source)
      tgt_ct        <- as.character(pair$target)

      if (ligand_gene %in% rownames(expr_mat)) {
        src_cells <- intersect(valid_cells,
                               cell_ids_all[ct_by_cell[cell_ids_all] == src_ct])
        if (length(src_cells) > 0) {
          lr_score[src_cells] <- lr_score[src_cells] +
            as.numeric(expr_mat[ligand_gene, src_cells, drop = TRUE])
        }
      }
      if (receptor_gene %in% rownames(expr_mat)) {
        tgt_cells <- intersect(valid_cells,
                               cell_ids_all[ct_by_cell[cell_ids_all] == tgt_ct])
        if (length(tgt_cells) > 0) {
          lr_score[tgt_cells] <- lr_score[tgt_cells] +
            as.numeric(expr_mat[receptor_gene, tgt_cells, drop = TRUE])
        }
      }
    }

    umap_df$lr_activity <- log1p(lr_score[umap_df$cell_ID])
    umap_df$is_focus    <- umap_df$celltype %in% (focus_celltype %||% character(0))
    umap_df             <- umap_df[order(umap_df$lr_activity), ]

    p6 <- ggplot2::ggplot(umap_df,
             ggplot2::aes(x = umap_1, y = umap_2, color = lr_activity)) +
      ggplot2::geom_point(size = 0.35, alpha = 0.85) +
      ggplot2::scale_color_viridis_c(
        option = "magma", name = "L-R activity\n(log1p)"
      ) +
      ggplot2::labs(
        title    = sample_plot_title(sample_id,
                     "CCI L-R Activity - UMAP"),
        subtitle = paste0("Per-cell score = sum of ligand expression (when ",
                          "sender) + receptor expression (when receiver) ",
                          "across top ", top_n_lr, " L-R pairs"),
        x = "UMAP 1", y = "UMAP 2"
      ) +
      presentation_theme(base_size = 11, legend_position = "right")

    save_presentation_plot(
      plot = p6,
      filename = file.path(out_dir,
                           paste0(sample_id, "_liana_umap_lr_activity.png")),
      width = 12, height = 10, dpi = 150
    )
    cat("\u2713 LIANA UMAP L-R activity saved\n")

    # B-cell-only variant: same score but only keep cells whose type
    # participates in any B-cell L-R pair; others dimmed. Goes to
    # B_cell_specific/.
    if (!is.null(focus_celltype) && length(focus_celltype) > 0) {
      bcell_agg <- agg[agg$source %in% focus_celltype |
                         agg$target %in% focus_celltype, , drop = FALSE]
      bcell_agg <- bcell_agg[order(bcell_agg$aggregate_rank,
                                   na.last = TRUE), , drop = FALSE]
      bcell_lr  <- head(bcell_agg, top_n_lr)
      lr_score_b <- stats::setNames(rep(0, length(cell_ids_all)), cell_ids_all)
      for (i in seq_len(nrow(bcell_lr))) {
        pair          <- bcell_lr[i, ]
        ligand_gene   <- strsplit(as.character(pair$ligand_complex),   "_")[[1]][1]
        receptor_gene <- strsplit(as.character(pair$receptor_complex), "_")[[1]][1]
        src_ct        <- as.character(pair$source)
        tgt_ct        <- as.character(pair$target)
        if (ligand_gene %in% rownames(expr_mat)) {
          src_cells <- intersect(valid_cells,
                                 cell_ids_all[ct_by_cell[cell_ids_all] == src_ct])
          if (length(src_cells) > 0) {
            lr_score_b[src_cells] <- lr_score_b[src_cells] +
              as.numeric(expr_mat[ligand_gene, src_cells, drop = TRUE])
          }
        }
        if (receptor_gene %in% rownames(expr_mat)) {
          tgt_cells <- intersect(valid_cells,
                                 cell_ids_all[ct_by_cell[cell_ids_all] == tgt_ct])
          if (length(tgt_cells) > 0) {
            lr_score_b[tgt_cells] <- lr_score_b[tgt_cells] +
              as.numeric(expr_mat[receptor_gene, tgt_cells, drop = TRUE])
          }
        }
      }
      umap_df$lr_activity_bcell <- log1p(lr_score_b[umap_df$cell_ID])
      umap_df <- umap_df[order(umap_df$lr_activity_bcell), ]
      p6b <- ggplot2::ggplot(umap_df,
               ggplot2::aes(x = umap_1, y = umap_2, color = lr_activity_bcell)) +
        ggplot2::geom_point(size = 0.35, alpha = 0.85) +
        ggplot2::scale_color_viridis_c(
          option = "magma", name = "B-cell L-R\nactivity (log1p)"
        ) +
        ggplot2::labs(
          title    = sample_plot_title(sample_id,
                       "B-cell CCI L-R Activity - UMAP"),
          subtitle = paste0("Per-cell score across top ", nrow(bcell_lr),
                            " ", paste(focus_celltype, collapse = "/"),
                            " L-R pairs"),
          x = "UMAP 1", y = "UMAP 2"
        ) +
        presentation_theme(base_size = 11, legend_position = "right")
      save_presentation_plot(
        plot = p6b,
        filename = file.path(bcell_dir,
                             paste0(sample_id, "_liana_umap_lr_activity_bcell.png")),
        width = 12, height = 10, dpi = 150
      )
      cat("\u2713 LIANA B-cell UMAP L-R activity saved \u2192", bcell_dir, "\n")
    }
  }, error = function(e) {
    cat("\u26A0 LIANA CCI UMAP failed:", conditionMessage(e), "\n")
    message("  Detail: ", conditionMessage(e))
  })

  # Plot 6b (centroid-level spatial CCI overlay) removed — the centroid-
  # arrow rendering was not interpretable on tissue-scale data.

  # -- Plot 7: Dominance scatter (outgoing vs incoming per cell type) ---------
  tryCatch({
    if (!requireNamespace("ggrepel", quietly = TRUE)) {
      cat("\u26A0 LIANA dominance scatter skipped: 'ggrepel' is not installed.\n")
    } else if (is.null(count_df) || nrow(count_df) == 0) {
      stop("No interaction count data.")
    } else {

    out_per <- stats::aggregate(n_interactions ~ source, data = count_df, FUN = sum)
    in_per  <- stats::aggregate(n_interactions ~ target, data = count_df, FUN = sum)
    dom_df  <- merge(out_per, in_per, by.x = "source", by.y = "target", all = TRUE)
    names(dom_df) <- c("celltype", "outgoing", "incoming")
    dom_df[is.na(dom_df)] <- 0
    dom_df$total    <- dom_df$outgoing + dom_df$incoming
    dom_df$is_focus <- dom_df$celltype %in% (focus_celltype %||% character(0))

    p7 <- ggplot2::ggplot(dom_df,
             ggplot2::aes(x = outgoing, y = incoming,
                          size = total, color = celltype,
                          label = celltype)) +
      ggplot2::geom_point(alpha = 0.75) +
      ggplot2::geom_point(
        data  = dom_df[dom_df$is_focus, ],
        color = "firebrick", shape = 21, stroke = 2, fill = NA
      ) +
      ggrepel::geom_label_repel(
        data  = {
          top8   <- head(dom_df[order(-dom_df$total), ], 8)
          focus  <- dom_df[dom_df$is_focus, ]
          unique(rbind(top8, focus))
        },
        size            = 3,
        max.overlaps    = 30,
        show.legend     = FALSE,
        force           = 2,
        min.segment.length = 0.2,
        fill            = ggplot2::alpha("white", 0.85),
        label.size      = 0.15,
        label.padding   = grid::unit(0.18, "lines"),
        label.r         = grid::unit(0.12, "lines"),
        segment.color   = "grey40",
        segment.size    = 0.25
      ) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "grey50") +
      ggplot2::scale_size_continuous(range = c(2, 10), guide = "none") +
      ggplot2::scale_color_discrete(guide = "none") +
      ggplot2::labs(
        title    = sample_plot_title(sample_id, "CCI Sender/Receiver Dominance"),
        subtitle = "x = total outgoing L-R pairs as sender; y = total incoming as receiver; diagonal = balanced",
        x        = "Outgoing interactions (as sender)",
        y        = "Incoming interactions (as receiver)"
      ) +
      presentation_theme(base_size = 12, legend_position = "none")

    save_presentation_plot(
      plot = p7,
      filename = file.path(out_dir, paste0(sample_id, "_liana_dominance.png")),
      width = 10, height = 9, dpi = 150
    )
    cat("\u2713 LIANA dominance scatter saved\n")
    }
  }, error = function(e) {
    cat("\u26A0 LIANA dominance scatter failed:", conditionMessage(e), "\n")
  })

  # -- Plot 8: Top sender-receiver pairs (information flow bar chart) ---------
  tryCatch({
    if (is.null(count_df) || nrow(count_df) == 0) stop("No interaction count data.")

    n_pairs  <- min(25L, nrow(count_df))
    flow_df  <- head(count_df[order(-count_df$n_interactions), ], n_pairs)
    flow_df$pair_label <- paste0(flow_df$source, " \u2192 ", flow_df$target)
    flow_df$is_focus   <- flow_df$source %in% (focus_celltype %||% character(0)) |
                          flow_df$target %in% (focus_celltype %||% character(0))

    p8 <- ggplot2::ggplot(flow_df,
             ggplot2::aes(
               x    = reorder(pair_label, n_interactions),
               y    = n_interactions,
               fill = source,
               color = is_focus
             )) +
      ggplot2::geom_col(linewidth = 0.4) +
      ggplot2::scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA),
                                  guide = "none") +
      ggplot2::coord_flip(clip = "off") +
      ggplot2::labs(
        title    = sample_plot_title(sample_id, "CCI Information Flow"),
        subtitle = paste0("Top ", n_pairs, " sender-receiver pairs by total L-R interactions"),
        x        = "Sender \u2192 Receiver",
        y        = "N significant L-R pairs",
        fill     = "Sender"
      ) +
      presentation_theme(base_size = 11, legend_position = "right") +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
      ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 160))

    save_presentation_plot(
      plot = p8,
      filename = file.path(out_dir, paste0(sample_id, "_liana_information_flow.png")),
      width = 14, height = 10, dpi = 150
    )
    cat("\u2713 LIANA information flow saved\n")
  }, error = function(e) {
    cat("\u26A0 LIANA information flow failed:", conditionMessage(e), "\n")
  })

  # -- Plot 8b: B-cell-only information flow ---------------------------------
  tryCatch({
    if (is.null(focus_celltype) || length(focus_celltype) == 0) {
      cat("\u26A0 LIANA B-cell information flow skipped: no focus_celltype resolved.\n")
    } else if (is.null(count_df) || nrow(count_df) == 0) {
      stop("No interaction count data.")
    } else {
      bcell_flow <- count_df[count_df$source %in% focus_celltype |
                               count_df$target %in% focus_celltype, ]
      if (nrow(bcell_flow) == 0) {
        cat("\u26A0 LIANA B-cell information flow skipped: no B-cell edges.\n")
      } else {
        n_bcell <- min(30L, nrow(bcell_flow))
        bcell_flow <- head(bcell_flow[order(-bcell_flow$n_interactions), ], n_bcell)
        bcell_flow$pair_label <- paste0(bcell_flow$source, " \u2192 ", bcell_flow$target)
        bcell_flow$direction  <- ifelse(
          bcell_flow$source %in% focus_celltype,
          paste0(paste(focus_celltype, collapse = "/"), " as sender"),
          paste0(paste(focus_celltype, collapse = "/"), " as receiver")
        )

        p8b <- ggplot2::ggplot(bcell_flow,
                 ggplot2::aes(
                   x    = reorder(pair_label, n_interactions),
                   y    = n_interactions,
                   fill = direction
                 )) +
          ggplot2::geom_col(color = "grey30", linewidth = 0.2) +
          ggplot2::scale_fill_manual(
            values = stats::setNames(
              c("firebrick", "steelblue"),
              c(paste0(paste(focus_celltype, collapse = "/"), " as sender"),
                paste0(paste(focus_celltype, collapse = "/"), " as receiver"))
            ),
            name = NULL
          ) +
          ggplot2::coord_flip(clip = "off") +
          ggplot2::labs(
            title    = sample_plot_title(sample_id,
                         "B-cell CCI Information Flow"),
            subtitle = paste0("Top ", n_bcell, " ",
                              paste(focus_celltype, collapse = "/"),
                              " sender/receiver pairs by total L-R interactions"),
            x        = "Sender \u2192 Receiver",
            y        = "N significant L-R pairs"
          ) +
          presentation_theme(base_size = 11, legend_position = "top") +
          ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 160))

        save_presentation_plot(
          plot     = p8b,
          filename = file.path(bcell_dir,
                               paste0(sample_id, "_liana_information_flow_bcell.png")),
          width    = 14, height = 10, dpi = 150
        )
        cat("\u2713 LIANA B-cell information flow saved \u2192", bcell_dir, "\n")
      }
    }
  }, error = function(e) {
    cat("\u26A0 LIANA B-cell information flow failed:", conditionMessage(e), "\n")
  })

  invisible(NULL)
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
#' @param expr_cache    Optional pre-extracted expression cache list.
#' @return liana_result data frame

run_liana <- function(gobj,
                      sample_id,
                      output_dir,
                      celltype_col      = NULL,
                      methods           = NULL,
                      expr_cache        = NULL,
                      focus_celltype    = "B cell",
                      sample_row        = NULL,
                      sample_sheet_path = NULL) {
  
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
    
    # FIX #6: Use pre-extracted cache when available.
    counts_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$raw)) {
      expr_cache$raw
    } else {
      .giotto_get_expression(gobj, values = "raw", output = "matrix")
    }
    logcounts_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$normalized)) {
      expr_cache$normalized
    } else {
      .giotto_get_expression(gobj, values = "normalized", output = "matrix")
    }
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

    # Resolve focus cell type BEFORE the dotplot tryCatch so it survives any dotplot failure
    {
      src_col_early <- .first_existing_column(liana_agg, c("source", "sender", "Sender"))
      tgt_col_early <- .first_existing_column(liana_agg, c("target", "receiver", "Receiver"))
      all_cts_early <- unique(c(liana_agg[[src_col_early]], liana_agg[[tgt_col_early]]))
      resolved_focus <- if (!is.null(focus_celltype) && focus_celltype %in% all_cts_early) {
        focus_celltype
      } else if (!is.null(focus_celltype)) {
        hits <- all_cts_early[grepl(focus_celltype, all_cts_early, ignore.case = TRUE)]
        if (length(hits) > 0) hits else NULL
      } else {
        NULL
      }
      cat("  Focus cell type resolved to:", paste(resolved_focus %||% "none", collapse = ", "), "\n")
    }

    # Dotplots (main + B-cell) are rendered inside plot_liana_extended()
    # so that manual re-runs on a cached liana_aggregate.csv also regenerate
    # them without re-running the full LIANA pipeline.

    plot_liana_extended(
      liana_agg         = liana_agg,
      meta              = meta,
      gobj              = gobj,
      expr_cache        = expr_cache,
      out_dir           = out_dir,
      sample_id         = sample_id,
      celltype_col      = celltype_col,
      focus_celltype    = resolved_focus,
      sample_row        = sample_row,
      sample_sheet_path = sample_sheet_path
    )

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
#' @param expr_cache    Optional pre-extracted expression cache list.
#' @return NicheNet results list

run_nichenet <- function(gobj,
                         sample_id,
                         output_dir,
                         celltype_col       = NULL,
                         sender_celltypes   = NULL,
                         receiver_celltype  = NULL,
                         target_genes       = NULL,
                         top_n_ligands      = 20,
                         network_dir        = NULL,
                         result_root_dir    = NULL,
                         expr_cache         = NULL) {
  
  if (!requireNamespace("nichenetr", quietly = TRUE))
    stop("nichenetr not installed.\n",
         'Install: devtools::install_github("saeyslab/nichenetr")')
  
  cat("\n--- NicheNet: Ligand activity prediction ---\n")
  cat("  Sender(s):  ", paste(sender_celltypes,  collapse = ", "), "\n")
  cat("  Receiver:   ", receiver_celltype, "\n\n")
  
  if (is.null(sender_celltypes) || is.null(receiver_celltype))
    stop("sender_celltypes and receiver_celltype must both be specified.")
  
  if (is.null(result_root_dir) || !nzchar(result_root_dir)) {
    result_root_dir <- file.path(output_dir, "10_CCI_Analysis", "nichenet")
  }
  
  out_dir <- file.path(result_root_dir,
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
  
  # FIX #6: Use pre-extracted cache when available.
  expr_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$normalized)) {
    expr_cache$normalized
  } else {
    .giotto_get_expression(gobj, values = "normalized", output = "matrix")
  }
  
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
  
  # FIX #4: Removed kidney-specific HAVCR1 proxy. NicheNet requires explicit
  # target genes — a tissue-specific proxy like HAVCR1 fails silently for
  # immune cell receivers (B cells, T cells, etc.) and produces misleading
  # "insufficient DE genes" warnings. Users must supply target genes via
  # config.yaml cci.target_genes or cci.target_genes_by_receiver.
  receiver_de <- if (!is.null(target_genes) && length(target_genes) > 0) {
    intersect(unique(as.character(target_genes)), rownames(expr_mat))
  } else {
    # No target genes supplied. Rather than using a tissue-specific proxy
    # that would fail silently for immune cell types, we require explicit
    # target genes. These can be derived from step 06 (cluster markers) or
    # step 12 (spatial DE) and configured via cci.target_genes or
    # cci.target_genes_by_receiver in config.yaml.
    cat("\u26A0 No target genes provided for NicheNet receiver '",
        receiver_celltype, "'.\n", sep = "")
    cat("  NicheNet requires a set of DE/response genes in the receiver.\n")
    cat("  Options:\n")
    cat("    1. Set cci.target_genes in config.yaml (applies to all receivers)\n")
    cat("    2. Set cci.target_genes_by_receiver.<celltype> in config.yaml\n")
    cat("    3. Use marker genes from step 06_Differential_Expression.R\n\n")
    character(0)
  }
  
  if (length(receiver_de) < 10) {
    cat("\u26A0 Insufficient target genes for NicheNet (<10 after intersection with expression matrix).\n")
    cat("  Skipping NicheNet for receiver '", receiver_celltype, "'.\n\n", sep = "")
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
  
  # --- P1: top-20 ligand activity bar plot ---
  tryCatch({
    la_sorted <- ligand_activities[order(-ligand_activities$pearson), , drop = FALSE]
    la_top    <- head(la_sorted, 20)
    p_la <- ggplot2::ggplot(
      la_top,
      ggplot2::aes(x = stats::reorder(test_ligand, pearson), y = pearson)
    ) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL, y = "Pearson score",
        title = sprintf("%s \u2014 NicheNet top-20 ligands \u2192 %s",
                        sample_id, receiver_celltype)
      ) +
      ggplot2::theme_minimal(base_size = 11)
    ggplot2::ggsave(
      file.path(out_dir, paste0(sample_id, "_nichenet_top_ligand_activity.png")),
      p_la, width = 6, height = 7, dpi = 150
    )
  }, error = function(e) {
    cat("\u26A0 NicheNet ligand activity plot failed:", conditionMessage(e), "\n")
  })
  
  # --- P2: ligand -> target heatmap (top 20 ligands x top 30 targets) ---
  tryCatch({
    top_ligs_p2 <- head(ligand_activities[order(-ligand_activities$pearson),
                                          "test_ligand"], 20)
    top_ligs_p2 <- intersect(top_ligs_p2, rownames(ligand_target_matrix))
    if (length(top_ligs_p2) > 0) {
      lt_sub <- as.matrix(ligand_target_matrix[top_ligs_p2, , drop = FALSE])
      col_sums <- colSums(lt_sub, na.rm = TRUE)
      top_targs_p2 <- names(sort(col_sums, decreasing = TRUE))[
        seq_len(min(30L, length(col_sums)))
      ]
      lt_sub <- lt_sub[, top_targs_p2, drop = FALSE]
      heat_df <- as.data.frame(as.table(lt_sub), stringsAsFactors = FALSE)
      names(heat_df) <- c("ligand", "target", "weight")
      heat_df$ligand <- factor(heat_df$ligand, levels = rev(top_ligs_p2))
      heat_df$target <- factor(heat_df$target, levels = top_targs_p2)
      p_lt <- ggplot2::ggplot(heat_df,
                              ggplot2::aes(target, ligand, fill = weight)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient(low = "white", high = "#440154",
                                     name = "regulatory\npotential") +
        ggplot2::labs(
          x = "Target gene", y = "Ligand",
          title = sprintf("%s \u2014 NicheNet ligand-target potential",
                          sample_id)
        ) +
        ggplot2::theme_minimal(base_size = 9) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(
          angle = 90, hjust = 1, vjust = 0.5))
      ggplot2::ggsave(
        file.path(out_dir,
                  paste0(sample_id, "_nichenet_ligand_target_heatmap.png")),
        p_lt, width = 10, height = 7, dpi = 150
      )
    }
  }, error = function(e) {
    cat("\u26A0 NicheNet ligand-target heatmap failed:", conditionMessage(e), "\n")
  })
  
  # --- P3, P4, P5: B-cell polygon overlays ---
  tryCatch({
    bcell_ids <- meta$cell_ID[meta[[celltype_col]] == receiver_celltype]
    poly_df   <- .cci_extract_polygon_df(gobj)
    if (!is.null(poly_df) && length(bcell_ids) > 0) {
      poly_df <- .bcell_polygon_base(poly_df, bcell_ids)
      
      # P3: sender map (top-5 senders coloured, B cells outlined)
      ct_vec <- stats::setNames(as.character(meta[[celltype_col]]), meta$cell_ID)
      senders_p3 <- utils::head(sender_celltypes, 5)
      .plot_bcell_polygon_category(
        poly_df = poly_df,
        ct_vec  = ct_vec,
        focus_labels = senders_p3,
        title   = sprintf("%s \u2014 NicheNet senders around %s",
                          sample_id, receiver_celltype),
        outfile = file.path(out_dir,
                            paste0(sample_id, "_nichenet_bcell_sender_spatial.png"))
      )
      
      # P4: top-6 ligand expression on polygons (senders show, receivers outlined)
      top_ligs_spatial <- utils::head(
        ligand_activities[order(-ligand_activities$pearson), "test_ligand"], 6
      )
      .plot_bcell_polygon_grid(
        poly_df  = poly_df,
        expr_mat = expr_mat,
        genes    = as.character(top_ligs_spatial),
        bcell_ids = bcell_ids,
        outfile  = file.path(out_dir,
                             paste0(sample_id, "_nichenet_bcell_top_ligand_expr.png")),
        ncol_grid = 3, width = 12, height = 8,
        title = sprintf("%s \u2014 top-6 NicheNet ligands", sample_id)
      )
      
      # P5: top-6 target gene expression on polygons
      top_targs_spatial <- {
        top_ligs_for_targets <- utils::head(
          ligand_activities[order(-ligand_activities$pearson), "test_ligand"], 20
        )
        top_ligs_for_targets <- intersect(as.character(top_ligs_for_targets),
                                          rownames(ligand_target_matrix))
        if (length(top_ligs_for_targets) > 0) {
          lt_sub   <- as.matrix(ligand_target_matrix[top_ligs_for_targets, , drop = FALSE])
          col_sums <- colSums(lt_sub, na.rm = TRUE)
          utils::head(names(sort(col_sums, decreasing = TRUE)), 6)
        } else character(0)
      }
      if (length(top_targs_spatial) > 0) {
        .plot_bcell_polygon_grid(
          poly_df   = poly_df,
          expr_mat  = expr_mat,
          genes     = as.character(top_targs_spatial),
          bcell_ids = bcell_ids,
          outfile   = file.path(out_dir,
                                paste0(sample_id, "_nichenet_bcell_top_target_expr.png")),
          ncol_grid = 3, width = 12, height = 8,
          title = sprintf("%s \u2014 top-6 NicheNet targets in %s",
                          sample_id, receiver_celltype)
        )
      }
    }
  }, error = function(e) {
    cat("\u26A0 NicheNet B-cell polygon overlays failed:", conditionMessage(e), "\n")
  })
  
  cat("\u2713 NicheNet complete. Top ligands:\n")
  print(top_ligands)
  cat("  Results saved to:", out_dir, "\n\n")
  
  invisible(list(
    ligand_activities     = ligand_activities,
    top_ligands           = top_ligands,
    ligand_target_matrix  = ligand_target_matrix
  ))
}

run_nichenet_batch <- function(gobj,
                               sample_id,
                               output_dir,
                               celltype_col             = NULL,
                               sender_celltypes         = NULL,
                               receiver_celltype        = NULL,
                               target_genes             = NULL,
                               target_genes_by_receiver = NULL,
                               network_dir              = NULL,
                               mode                     = c("single", "all_senders_to_receiver", "all_pairs"),
                               use_spatial_filter       = FALSE,
                               proximity_enrichment_path = NULL,
                               spatial_padj_threshold   = 0.05,
                               min_cells_per_celltype   = 5,
                               include_self_pairs       = FALSE,
                               result_root_dir          = NULL,
                               expr_cache               = NULL) {
  mode <- match.arg(mode)
  cat("NicheNet batch mode:", mode, "\n")
  
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (is.null(celltype_col)) {
    sup_cols <- grep("celltype_.*_supervised$", names(meta), value = TRUE)
    celltype_col <- if (length(sup_cols) > 0) sup_cols[1] else "leiden_clust"
  }
  
  comparisons <- .build_nichenet_comparison_table(
    metadata = meta,
    celltype_col = celltype_col,
    mode = mode,
    sender_celltypes = sender_celltypes,
    receiver_celltype = receiver_celltype,
    min_cells_per_celltype = min_cells_per_celltype,
    include_self_pairs = include_self_pairs
  )
  
  if (nrow(comparisons) == 0) {
    cat("\u26A0 No NicheNet comparisons passed the initial filters.\n")
    return(invisible(NULL))
  }
  
  if (is.null(result_root_dir) || !nzchar(result_root_dir)) {
    result_root_dir <- file.path(output_dir, "10_CCI_Analysis", "nichenet")
  }
  nichenet_root <- result_root_dir
  dir.create(nichenet_root, recursive = TRUE, showWarnings = FALSE)
  
  utils::write.csv(
    comparisons[, c("sender_label", "receiver"), drop = FALSE],
    file.path(nichenet_root, paste0(sample_id, "_nichenet_requested_comparisons.csv")),
    row.names = FALSE
  )
  cat("Requested NicheNet comparisons:", nrow(comparisons), "\n")
  
  proximity_tbl <- if (isTRUE(use_spatial_filter)) {
    .load_proximity_enrichment(
      output_dir = output_dir,
      sample_id = sample_id,
      proximity_enrichment_path = proximity_enrichment_path
    )
  } else {
    NULL
  }
  
  if (isTRUE(use_spatial_filter) && is.null(proximity_tbl)) {
    cat("\u26A0 Spatial filter requested for NicheNet, but no proximity enrichment table was found.\n")
    cat("  Proceeding without the spatial filter.\n\n")
  }
  
  if (isTRUE(use_spatial_filter) && !is.null(proximity_tbl)) {
    keep_idx <- vapply(seq_len(nrow(comparisons)), function(i) {
      .pair_is_spatially_supported(
        sender = comparisons$sender_label[[i]],
        receiver = comparisons$receiver[[i]],
        proximity_tbl = proximity_tbl,
        spatial_padj_threshold = spatial_padj_threshold
      )
    }, logical(1))
    comparisons <- comparisons[keep_idx, , drop = FALSE]
    if (nrow(comparisons) == 0) {
      cat("\u26A0 No NicheNet comparisons passed the spatial proximity filter.\n")
      return(invisible(NULL))
    }
    cat("NicheNet spatial filter retained", nrow(comparisons), "comparison(s).\n\n")
    utils::write.csv(
      comparisons[, c("sender_label", "receiver"), drop = FALSE],
      file.path(nichenet_root, paste0(sample_id, "_nichenet_spatially_retained_comparisons.csv")),
      row.names = FALSE
    )
  }
  
  if (mode == "all_pairs" && is.null(.normalize_named_character_list(target_genes_by_receiver)) &&
      (is.null(target_genes) || length(target_genes) == 0)) {
    stop("NicheNet mode 'all_pairs' requires target_genes_by_receiver or a fallback target_genes vector.")
  }
  
  summary_rows <- list()
  results <- list()
  
  for (i in seq_len(nrow(comparisons))) {
    sender_label <- comparisons$sender_label[[i]]
    sender_list <- unlist(comparisons$sender_list[[i]], use.names = FALSE)
    receiver_label <- comparisons$receiver[[i]]
    receiver_targets <- .resolve_target_genes_for_receiver(
      receiver_celltype = receiver_label,
      target_genes = target_genes,
      target_genes_by_receiver = target_genes_by_receiver
    )
    
    if (is.null(receiver_targets) || length(receiver_targets) < 10) {
      cat("\u26A0 Skipping NicheNet", sender_label, "->", receiver_label,
          ": fewer than 10 receiver target genes available.\n")
      next
    }
    
    nichenet_out <- tryCatch(
      run_nichenet(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = output_dir,
        celltype_col = celltype_col,
        sender_celltypes = sender_list,
        receiver_celltype = receiver_label,
        target_genes = receiver_targets,
        network_dir = network_dir,
        result_root_dir = nichenet_root,
        expr_cache = expr_cache
      ),
      error = function(e) {
        cat("\u26A0 NicheNet comparison failed for", sender_label, "->", receiver_label, ":", conditionMessage(e), "\n")
        NULL
      }
    )
    
    comparison_key <- paste(receiver_label, "from", sender_label, sep = "__")
    results[[comparison_key]] <- nichenet_out
    
    if (is.null(nichenet_out) || is.null(nichenet_out$ligand_activities)) {
      next
    }
    
    top_tbl <- as.data.frame(nichenet_out$ligand_activities, stringsAsFactors = FALSE)
    if ("pearson" %in% names(top_tbl)) {
      top_tbl <- top_tbl[order(-top_tbl$pearson), , drop = FALSE]
    }
    top_tbl <- utils::head(top_tbl, 10)
    if (nrow(top_tbl) == 0) {
      next
    }
    
    top_tbl$sender <- sender_label
    top_tbl$receiver <- receiver_label
    summary_rows[[length(summary_rows) + 1L]] <- top_tbl
  }
  
  summary_df <- if (length(summary_rows) > 0) dplyr::bind_rows(summary_rows) else data.frame()
  if (nrow(summary_df) > 0) {
    utils::write.csv(
      summary_df,
      file.path(nichenet_root, paste0(sample_id, "_nichenet_comparison_summary.csv")),
      row.names = FALSE
    )
  }
  
  invisible(list(
    comparisons = comparisons,
    summary = summary_df,
    results = results
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
#' @param expr_cache    Optional pre-extracted expression cache list.
#' @return MISTy results list

run_misty <- function(gobj,
                      sample_id,
                      output_dir,
                      target_genes  = NULL,
                      juxta_radius  = 50,
                      para_radius   = 200,
                      n_cores       = 4,
                      expr_cache    = NULL,
                      celltype_col   = NULL,
                      focus_celltype = "^B cell$") {
  
  if (!requireNamespace("mistyR", quietly = TRUE))
    stop("mistyR not installed.\n",
         'Install: BiocManager::install("mistyR")')
  # Guard: .misty_runtime_ready() should have already checked this, but we
  # verify here in case run_misty() is called directly.  The full diagnostic
  # message (with remediation steps) is printed by .misty_runtime_ready().
  ridge_ready <- tryCatch(requireNamespace("ridge", quietly = TRUE), error = function(e) FALSE)
  if (!isTRUE(ridge_ready)) {
    misty_check <- .misty_runtime_ready()
    cat("\u26A0", misty_check$reason, "\n")
    return(invisible(NULL))
  }
  
  cat("\n--- MISTy: Intercellular spatial modelling ---\n")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "misty")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Expression matrix (cells x genes) — log-normalised
  # FIX #6: Use pre-extracted cache when available.
  norm_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$normalized)) {
    expr_cache$normalized
  } else {
    .giotto_get_expression(gobj, values = "normalized", output = "matrix")
  }
  expr_mat <- t(as.matrix(norm_mat))
  
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
  
  # FIX #8: Guard against running MISTy on datasets too large for pairwise
  # distance computation. add_juxtaview() and add_paraview() build O(n^2)
  # distance matrices which at 200k cells would require ~300 GB of RAM.
  MAX_CELLS_MISTY <- 50000
  if (nrow(expr_mat) > MAX_CELLS_MISTY) {
    cat("\u26A0 MISTy: Dataset has", nrow(expr_mat), "cells, which exceeds the",
        MAX_CELLS_MISTY, "cell limit for pairwise distance computation.\n")
    cat("  Downsampling to", MAX_CELLS_MISTY, "cells.\n")
    set.seed(42)
    keep_idx <- sample(nrow(expr_mat), MAX_CELLS_MISTY)
    expr_mat <- expr_mat[keep_idx, , drop = FALSE]
    xy <- xy[keep_idx, , drop = FALSE]
    cat("  Cells after downsampling:", nrow(expr_mat), "\n\n")
  }
  
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
    
    # --- P6, P7, P8: mistyR native plots ---
    tryCatch({
      grDevices::png(
        file.path(out_dir, paste0(sample_id, "_misty_improvement_stats.png")),
        width = 900, height = 600, res = 150
      )
      print(mistyR::plot_improvement_stats(misty_results, "gain.R2"))
      grDevices::dev.off()
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      cat("\u26A0 MISTy improvement_stats plot failed:",
          conditionMessage(e), "\n")
    })
    tryCatch({
      grDevices::png(
        file.path(out_dir, paste0(sample_id, "_misty_view_contributions.png")),
        width = 900, height = 600, res = 150
      )
      print(mistyR::plot_view_contributions(misty_results))
      grDevices::dev.off()
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      cat("\u26A0 MISTy view_contributions plot failed:",
          conditionMessage(e), "\n")
    })
    for (misty_view in c("juxtaview_50", "paraview_200")) {
      tryCatch({
        grDevices::png(
          file.path(out_dir, sprintf("%s_misty_interaction_heatmap_%s.png",
                                     sample_id, misty_view)),
          width = 900, height = 900, res = 150
        )
        print(mistyR::plot_interaction_heatmap(misty_results, misty_view))
        grDevices::dev.off()
      }, error = function(e) {
        try(grDevices::dev.off(), silent = TRUE)
        cat("\u26A0 MISTy interaction_heatmap (", misty_view,
            ") failed: ", conditionMessage(e), "\n", sep = "")
      })
    }
    
    # --- P9, P10: B-cell polygon overlays ---
    tryCatch({
      celltype_col_local <- .resolve_celltype_column_auto(gobj, celltype_col)
      bcell_ids <- .resolve_bcell_ids(gobj, celltype_col_local, focus_celltype)
      poly_df   <- .cci_extract_polygon_df(gobj)
      if (!is.null(poly_df) && length(bcell_ids) > 0 &&
          !is.null(misty_results$improvements)) {
        poly_df <- .bcell_polygon_base(poly_df, bcell_ids)
        imp_df  <- as.data.frame(misty_results$improvements)
        
        # P9: top-9 paraview-boosted targets, polygon grid with expression
        gain_col <- intersect(c("gain.R2", "gain"), names(imp_df))[1]
        tgt_col  <- intersect(c("target", "Target"), names(imp_df))[1]
        if (!is.na(gain_col) && !is.na(tgt_col)) {
          imp_ord <- imp_df[order(-imp_df[[gain_col]]), , drop = FALSE]
          top_tgts_safe <- utils::head(unique(as.character(imp_ord[[tgt_col]])), 9L)
          # Map safe names back to original panel names
          name_map_path <- file.path(out_dir,
                                     paste0(sample_id, "_misty_feature_name_map.csv"))
          top_tgts <- top_tgts_safe
          if (file.exists(name_map_path)) {
            nm <- tryCatch(read.csv(name_map_path, stringsAsFactors = FALSE),
                           error = function(e) NULL)
            if (!is.null(nm) && all(c("original", "safe") %in% names(nm))) {
              idx <- match(top_tgts_safe, nm$safe)
              top_tgts <- ifelse(is.na(idx), top_tgts_safe, nm$original[idx])
            }
          }
          # Original expression matrix (before the make.names rename).
          expr_mat_orig <- if (!is.null(expr_cache) &&
                               !is.null(expr_cache$normalized)) {
            expr_cache$normalized
          } else {
            .giotto_get_expression(gobj, values = "normalized", output = "matrix")
          }
          .plot_bcell_polygon_grid(
            poly_df   = poly_df,
            expr_mat  = expr_mat_orig,
            genes     = as.character(top_tgts),
            bcell_ids = bcell_ids,
            outfile   = file.path(
              out_dir, paste0(sample_id, "_misty_bcell_top_targets_spatial.png")
            ),
            ncol_grid = 3, width = 12, height = 12,
            title = sprintf("%s \u2014 MISTy top-9 paraview-boosted targets",
                            sample_id)
          )
        }
        
        # P10: per-cell paraview importance score
        imp_detail <- if (!is.null(misty_results$importances))
          as.data.frame(misty_results$importances) else NULL
        if (!is.null(imp_detail)) {
          view_col <- intersect(c("view", "View"), names(imp_detail))[1]
          tgt2     <- intersect(c("Target", "target"), names(imp_detail))[1]
          pred_col <- intersect(c("Predictor", "predictor"), names(imp_detail))[1]
          imp_col  <- intersect(c("Importance", "importance"), names(imp_detail))[1]
          para_view_name <- grep("^para", unique(imp_detail[[view_col]]),
                                 value = TRUE, ignore.case = TRUE)[1]
          if (!is.na(para_view_name) && !is.na(tgt2) && !is.na(imp_col)) {
            imp_para <- imp_detail[imp_detail[[view_col]] == para_view_name, ,
                                    drop = FALSE]
            per_target_mean <- tapply(abs(imp_para[[imp_col]]),
                                      imp_para[[tgt2]], mean, na.rm = TRUE)
            genes_in_mat <- intersect(names(per_target_mean),
                                      rownames(expr_mat_orig))
            if (length(genes_in_mat) > 0) {
              expr_sub <- expr_mat_orig[genes_in_mat, , drop = FALSE]
              weights  <- per_target_mean[genes_in_mat]
              per_cell_score <- as.numeric(weights %*% as.matrix(expr_sub))
              names(per_cell_score) <- colnames(expr_sub)
              bcell_mean <- mean(per_cell_score[intersect(names(per_cell_score),
                                                          bcell_ids)],
                                 na.rm = TRUE)
              p_score <- .plot_bcell_polygon_panel(
                poly_df = poly_df, values = per_cell_score,
                title = sprintf("%s \u2014 MISTy paraview score", sample_id),
                subtitle = sprintf("Mean B-cell score: %.3f", bcell_mean),
                fill_label = "score"
              )
              ggplot2::ggsave(
                file.path(out_dir,
                          paste0(sample_id, "_misty_bcell_neighbourhood_score.png")),
                p_score, width = 9, height = 8, dpi = 150
              )
            }
          }
        }
      }
    }, error = function(e) {
      cat("\u26A0 MISTy B-cell polygon overlays failed:", conditionMessage(e), "\n")
    })
    
    cat("\u2713 MISTy complete. Results saved to:", out_dir, "\n\n")
    invisible(misty_results)
  } else {
    invisible(NULL)
  }
}


# ==============================================================================
# SECTION 4 — nnSVG: spatially variable gene detection
# ==============================================================================

# Pick a SEraster pixel side length that aggregates ~target cells per pixel.
# Unit-agnostic: extent and resolution share whatever unit spatialCoords() uses,
# so this works regardless of whether the export is in micrometers, mm, or
# native imaging units (CosMx ships with ~0.12 um/unit on this cluster, which
# silently no-ops a hard-coded resolution = 50).
.autopick_raster_resolution <- function(coords,
                                        target_cells_per_pixel = 30,
                                        floor_min = 10) {
  if (is.null(coords) || nrow(coords) < 100) return(NULL)
  rng_x <- diff(range(coords[, 1], na.rm = TRUE))
  rng_y <- diff(range(coords[, 2], na.rm = TRUE))
  if (!is.finite(rng_x) || !is.finite(rng_y) || rng_x <= 0 || rng_y <= 0) {
    return(NULL)
  }
  R <- sqrt(target_cells_per_pixel * rng_x * rng_y / nrow(coords))
  max(R, floor_min)
}

#' Run nnSVG: spatially variable gene detection
#'
#' @param gobj          Annotated Giotto object
#' @param sample_id     Sample identifier
#' @param output_dir    Root output directory
#' @param n_top_svgs    Number of top SVGs to report (default: 100)
#' @param n_cores       Parallel cores (default: 4)
#' @param expr_cache    Optional pre-extracted expression cache list.
#' @return nnSVG results data frame

run_nnsvg <- function(gobj,
                      sample_id,
                      output_dir,
                      n_top_svgs = 100,
                      n_cores    = 4,
                      expr_cache = NULL,
                      celltype_col   = NULL,
                      focus_celltype = "^B cell$",
                      raster_resolution = "auto") {
  
  if (!requireNamespace("nnSVG", quietly = TRUE))
    stop("nnSVG not installed.\n",
         'Install: BiocManager::install("nnSVG")')
  
  cat("\n--- nnSVG: Spatially variable gene detection ---\n")
  cat("  \u26A0 nnSVG is computationally heavy (~15-60 min per sample at 20k cells).\n")
  
  out_dir <- file.path(output_dir, "10_CCI_Analysis", "svg")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Build SpatialExperiment from Giotto
  # FIX #6: Use pre-extracted cache when available.
  counts_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$raw)) {
    expr_cache$raw
  } else {
    .giotto_get_expression(gobj, values = "raw", output = "matrix")
  }
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
  logcounts_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$normalized)) {
    expr_cache$normalized
  } else {
    .giotto_get_expression(gobj, values = "normalized", output = "matrix")
  }
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

  # Resolve celltype column once and attach to colData so SEraster::rasterizeCellType
  # aggregates onto the same pixel grid as gene expression. Used by P14 / P15.
  celltype_col_local <- .resolve_celltype_column_auto(gobj, celltype_col)
  if (!is.null(celltype_col_local)) {
    pdata_full <- as.data.frame(.giotto_pdata_dt(gobj))
    rownames(pdata_full) <- pdata_full$cell_ID
    ct_vec <- pdata_full[colnames(spe), celltype_col_local, drop = TRUE]
    SummarizedExperiment::colData(spe)$celltype <- as.character(ct_vec)
  }

  # ---- B-cell-focused pass (uses full spe, before tissue-wide downsampling) --
  # Runs nnSVG on B cells only to answer: "what varies spatially WITHIN the
  # B-cell compartment?". Complements the tissue-wide pass below, which asks
  # "what is spatially structured across the whole tissue, and which of those
  # top hits are B-cell enriched?". Output prefix: <sample>_bcell_nnSVG_*
  bcell_nnsvg_df <- tryCatch({
    celltype_col_local <- .resolve_celltype_column_auto(gobj, celltype_col)
    bcell_ids_all <- .resolve_bcell_ids(gobj, celltype_col_local, focus_celltype)
    bcell_cols    <- intersect(bcell_ids_all, colnames(spe))
    if (length(bcell_cols) < 50) {
      cat("\n--- nnSVG B-cell-focused pass skipped: only ",
          length(bcell_cols), " B cells (need >= 50).\n", sep = "")
      NULL
    } else {
      cat("\n--- nnSVG: B-cell-focused pass ---\n")
      cat("  B cells:", length(bcell_cols), "\n")
      spe_b <- spe[, bcell_cols]
      # Re-filter: drop genes with too few counts among B cells specifically.
      counts_b    <- SummarizedExperiment::assay(spe_b, "counts")
      min_spots_b <- max(1L, ceiling(0.005 * ncol(spe_b)))
      keep_b <- if (inherits(counts_b, "Matrix")) {
        Matrix::rowSums(counts_b >= 3) >= min_spots_b
      } else {
        base::rowSums(counts_b >= 3) >= min_spots_b
      }
      spe_b <- spe_b[keep_b, ]
      cat("  Genes after B-cell-level filtering:", nrow(spe_b), "\n\n")
      if (nrow(spe_b) < 5) {
        cat("\u26A0 B-cell pass skipped: < 5 genes survived B-cell-level filtering.\n")
        NULL
      } else {
        res_b <- nnSVG::nnSVG(spe_b, assay_name = "logcounts",
                              n_threads = n_cores)
        df_b  <- as.data.frame(SummarizedExperiment::rowData(res_b))
        df_b  <- df_b[order(df_b$rank), ]
        write.csv(df_b,
                  file.path(out_dir,
                            paste0(sample_id, "_bcell_nnSVG_results.csv")))
        # Rank plot
        lr_col_b <- intersect(c("LR_stat", "stat"), names(df_b))[1]
        if (!is.na(lr_col_b)) {
          df_plot_b <- df_b
          df_plot_b$.gene <- rownames(df_plot_b)
          p_rank_b <- ggplot2::ggplot(
            df_plot_b, ggplot2::aes(x = rank, y = .data[[lr_col_b]])
          ) +
            ggplot2::geom_point(size = 0.6, alpha = 0.5, colour = "grey40") +
            ggplot2::geom_point(
              data = utils::head(df_plot_b, 20L),
              colour = "firebrick", size = 1.2
            ) +
            ggplot2::labs(
              x = "nnSVG rank", y = lr_col_b,
              title = sprintf("%s \u2014 B-cell-only nnSVG: rank vs %s",
                              sample_id, lr_col_b)
            ) +
            ggplot2::theme_minimal(base_size = 11)
          ggplot2::ggsave(
            file.path(out_dir,
                      paste0(sample_id, "_bcell_nnSVG_rank_plot.png")),
            p_rank_b, width = 7, height = 5, dpi = 150
          )
        }
        cat("\u2713 B-cell-focused nnSVG complete.\n\n")
        df_b
      }
    }
  }, error = function(e) {
    cat("\u26A0 B-cell-focused nnSVG pass failed:", conditionMessage(e), "\n")
    NULL
  })
  
  # SEraster rasterization preprocessing. nnSVG's nearest-neighbor GP does
  # not scale to im-SRT cell counts (CosMx slides routinely exceed 50k cells,
  # and the paper's own guidance is Visium-scale ~5k spots). Pattanayak et
  # al. 2023 (PMC10491191) show that aggregating cells into square pixels
  # with SEraster before nnSVG is the recommended route for MERFISH/CosMx:
  # it preserves every cell's signal by averaging into bins rather than
  # discarding cells, and the pixel-level SpatialExperiment fits nnSVG's
  # design scale. Falls back to random subsampling if SEraster is not
  # installed.
  MAX_CELLS_NNSVG        <- 20000
  TARGET_CELLS_PER_PIXEL <- 30
  MIN_REDUCTION_OK       <- 3   # sub-3x reduction = bins are sub-cellular

  # Resolve "auto" / NULL / numeric into a single numeric R (or NULL = disabled).
  # "auto" picks resolution from coord extent so it is robust to whatever unit
  # the SPE coords use (CosMx exports often use ~0.12 um/unit, which silently
  # no-ops a hard-coded resolution = 50).
  resolved_resolution <- raster_resolution
  if (is.character(resolved_resolution) &&
      identical(tolower(resolved_resolution), "auto")) {
    resolved_resolution <- .autopick_raster_resolution(
      coords = SpatialExperiment::spatialCoords(spe),
      target_cells_per_pixel = TARGET_CELLS_PER_PIXEL
    )
    if (!is.null(resolved_resolution)) {
      cat("  Auto-picked SEraster resolution = ",
          sprintf("%.2f", resolved_resolution),
          " coord units (target ", TARGET_CELLS_PER_PIXEL,
          " cells/pixel)\n", sep = "")
    } else {
      cat("  \u26A0 Auto-pick failed (fewer than 100 cells or zero coord ",
          "extent); rasterization disabled.\n", sep = "")
    }
  }

  use_raster <- !is.null(resolved_resolution) &&
                is.numeric(resolved_resolution) &&
                resolved_resolution > 0 &&
                ncol(spe) > MAX_CELLS_NNSVG
  nnsvg_assay <- "logcounts"
  # Cell-level snapshot (used by rasterizeCellType on the same grid as the
  # gene-expression raster, AND as the source for the random-subsample
  # fallback below) and slot for the celltype pixel counts that P14/P15
  # consume.
  spe_cells <- spe
  ct_raster <- NULL
  if (use_raster) {
    if (!requireNamespace("SEraster", quietly = TRUE)) {
      cat("\u26A0 nnSVG: dataset has", ncol(spe),
          "cells but SEraster is not installed.\n",
          "  Falling back to random subsample (",
          MAX_CELLS_NNSVG, " cells).\n",
          "  Install: Rscript Parameters/Install_CCI_dependencies.R seraster\n",
          sep = "")
      set.seed(42)
      keep_idx <- sample(ncol(spe), MAX_CELLS_NNSVG)
      spe <- spe[, keep_idx]
      cat("  Cells after downsampling:", ncol(spe), "\n\n")
    } else {
      n_cells_in <- ncol(spe)
      cat("  Rasterizing with SEraster (resolution = ",
          sprintf("%.2f", resolved_resolution),
          " coord units, square pixels, fun = mean)...\n", sep = "")
      spe_raster <- tryCatch(
        SEraster::rasterizeGeneExpression(
          spe,
          assay_name = "logcounts",
          resolution = resolved_resolution,
          fun        = "mean",
          square     = TRUE
        ),
        error = function(e) {
          cat("\u26A0 SEraster::rasterizeGeneExpression failed: ",
              conditionMessage(e), "\n", sep = "")
          NULL
        }
      )

      reduction <- if (!is.null(spe_raster) && ncol(spe_raster) > 0) {
        n_cells_in / ncol(spe_raster)
      } else 0
      acceptable <- !is.null(spe_raster) &&
                    ncol(spe_raster) >= 50 &&
                    reduction >= MIN_REDUCTION_OK

      if (!acceptable) {
        cat("\u26A0 Rasterization rejected: pixels = ",
            if (is.null(spe_raster)) 0 else ncol(spe_raster),
            ", reduction = ", sprintf("%.2fx", reduction),
            " (need pixels >= 50 AND reduction >= ", MIN_REDUCTION_OK,
            "x).\n",
            "  Bins are at sub-cellular scale for these coord units.\n",
            "  Falling back to random subsample (", MAX_CELLS_NNSVG,
            " cells).\n",
            "  Override by setting cci.nnsvg_raster_resolution: ",
            "<larger number> in config.yaml.\n", sep = "")
        if (exists("spe_raster", inherits = FALSE)) rm(spe_raster)
        set.seed(42)
        keep_idx <- sample(ncol(spe_cells), MAX_CELLS_NNSVG)
        spe <- spe_cells[, keep_idx]
      } else {
        cat("  \u2713 Rasterized: ", n_cells_in, " cells \u2192 ",
            ncol(spe_raster), " pixels (",
            sprintf("%.1fx", reduction),
            " reduction)\n", sep = "")
        spe <- spe_raster
        # SEraster emits the aggregated assay under the name "pixelval".
        nnsvg_assay <- "pixelval"

        # Rasterize cell-type counts onto the same pixel grid. Pixel IDs
        # align with spe_raster because both rasterizations use the same
        # input coordinates + resolution. Required for P14 / P15 overlays.
        if (!is.null(celltype_col_local)) {
          ct_raster <- tryCatch(
            SEraster::rasterizeCellType(
              spe_cells,
              col_name   = "celltype",
              resolution = resolved_resolution,
              square     = TRUE
            ),
            error = function(e) {
              cat("\u26A0 SEraster::rasterizeCellType failed: ",
                  conditionMessage(e), "\n", sep = "")
              NULL
            }
          )
          if (!is.null(ct_raster)) {
            ct_counts_mat <- as.matrix(
              SummarizedExperiment::assay(ct_raster, 1)
            )
            tryCatch(
              write.csv(
                ct_counts_mat,
                file.path(out_dir,
                  paste0(sample_id, "_celltype_pixel_counts.csv"))
              ),
              error = function(e) {
                cat("\u26A0 celltype pixel counts CSV write failed: ",
                    conditionMessage(e), "\n", sep = "")
              }
            )
            cat("  \u2713 Celltype rasterized: ",
                nrow(ct_counts_mat), " types \u00D7 ",
                ncol(ct_counts_mat), " pixels\n", sep = "")
          }
        }
      }
    }
  }

  cat("  Genes after filtering:", nrow(spe), "\n")
  cat("  Spots (",
      if (nnsvg_assay == "pixelval") "pixels" else "cells",
      "): ", ncol(spe), "\n\n", sep = "")

  set.seed(42)
  nnsvg_res <- tryCatch(
    nnSVG::nnSVG(spe, assay_name = nnsvg_assay, n_threads = n_cores),
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
    
    # --- P11: rank plot (LR statistic vs rank, label top-20) ---
    tryCatch({
      lr_col <- intersect(c("LR_stat", "stat"), names(result_df))[1]
      rank_df <- result_df
      rank_df$.gene <- rownames(rank_df)
      if (!is.na(lr_col)) {
        p_rank <- ggplot2::ggplot(
          rank_df, ggplot2::aes(x = rank, y = .data[[lr_col]])
        ) +
          ggplot2::geom_point(size = 0.6, alpha = 0.5, colour = "grey40") +
          ggplot2::geom_point(
            data = utils::head(rank_df, 20L),
            colour = "firebrick", size = 1.2
          ) +
          ggplot2::labs(
            x = "nnSVG rank", y = lr_col,
            title = sprintf("%s \u2014 top-20 spatially variable genes", sample_id)
          ) +
          ggplot2::theme_minimal(base_size = 11)
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p_rank <- p_rank + ggrepel::geom_text_repel(
            data = utils::head(rank_df, 20L),
            ggplot2::aes(label = .gene), size = 3, max.overlaps = 30L
          )
        }
        ggplot2::ggsave(
          file.path(out_dir, paste0(sample_id, "_nnSVG_rank_plot.png")),
          p_rank, width = 8, height = 6, dpi = 150
        )
      }
    }, error = function(e) {
      cat("\u26A0 nnSVG rank plot failed:", conditionMessage(e), "\n")
    })
    
    # --- P12, P13: B-cell polygon overlays ---
    tryCatch({
      celltype_col_local <- .resolve_celltype_column_auto(gobj, celltype_col)
      bcell_ids <- .resolve_bcell_ids(gobj, celltype_col_local, focus_celltype)
      poly_df   <- .cci_extract_polygon_df(gobj)
      top20 <- utils::head(rownames(result_df), 20L)
      expr_mat_nnsvg <- if (!is.null(expr_cache) &&
                            !is.null(expr_cache$normalized)) {
        expr_cache$normalized
      } else {
        .giotto_get_expression(gobj, values = "normalized", output = "matrix")
      }
      
      # P12: 4x5 grid of top-20 SVGs on polygons, B cells outlined
      if (!is.null(poly_df) && length(top20) > 0) {
        poly_df <- .bcell_polygon_base(poly_df, bcell_ids)
        .plot_bcell_polygon_grid(
          poly_df   = poly_df,
          expr_mat  = expr_mat_nnsvg,
          genes     = as.character(top20),
          bcell_ids = bcell_ids,
          outfile   = file.path(out_dir,
                                paste0(sample_id, "_nnSVG_top20_spatial_grid.png")),
          ncol_grid = 5, width = 18, height = 14,
          title = sprintf("%s \u2014 nnSVG top-20 SVGs", sample_id)
        )
      }
      
      # P13: per-SVG B-cell enrichment (mean in B cells / mean overall)
      if (length(bcell_ids) > 5 && length(top20) > 0) {
        genes_present <- intersect(top20, rownames(expr_mat_nnsvg))
        if (length(genes_present) > 0) {
          expr_sub  <- expr_mat_nnsvg[genes_present, , drop = FALSE]
          mean_all  <- Matrix::rowMeans(expr_sub)
          bcell_in  <- intersect(bcell_ids, colnames(expr_sub))
          mean_b    <- Matrix::rowMeans(expr_sub[, bcell_in, drop = FALSE])
          enr       <- mean_b / pmax(1e-6, mean_all)
          enr_df    <- data.frame(gene = names(enr), fold = as.numeric(enr))
          enr_df    <- enr_df[order(-enr_df$fold), ]
          enr_df$gene <- factor(enr_df$gene, levels = rev(enr_df$gene))
          p_enr <- ggplot2::ggplot(enr_df,
                                    ggplot2::aes(x = gene, y = fold)) +
            ggplot2::geom_col(fill = "mediumblue") +
            ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                                colour = "grey40") +
            ggplot2::coord_flip() +
            ggplot2::labs(
              x = NULL, y = "B-cell / overall expression",
              title = sprintf("%s \u2014 top-20 SVGs: B-cell enrichment",
                              sample_id)
            ) +
            ggplot2::theme_minimal(base_size = 11)
          ggplot2::ggsave(
            file.path(out_dir,
                      paste0(sample_id, "_nnSVG_bcell_enrichment_bar.png")),
            p_enr, width = 7, height = 7, dpi = 150
          )
        }
      }
    }, error = function(e) {
      cat("\u26A0 nnSVG B-cell polygon overlays failed:", conditionMessage(e), "\n")
    })

    # --- P14: pixel-level gene x celltype correlation heatmap ---
    # Pearson correlation between each top-20 SVG gene's pixel-mean expression
    # and each celltype's pixel count. High r means that gene's spatial
    # pattern tracks where that celltype lives. Requires rasterization.
    tryCatch({
      if (!is.null(ct_raster) && exists("spe_raster", inherits = FALSE) &&
          nnsvg_assay == "pixelval") {
        top20_p14 <- utils::head(rownames(result_df), 20L)
        expr_mat  <- SummarizedExperiment::assay(spe_raster, "pixelval")
        ct_mat    <- as.matrix(SummarizedExperiment::assay(ct_raster, 1))
        common_px <- intersect(colnames(expr_mat), colnames(ct_mat))
        expr_sub  <- expr_mat[intersect(top20_p14, rownames(expr_mat)),
                              common_px, drop = FALSE]
        ct_sub    <- ct_mat[, common_px, drop = FALSE]
        # Drop celltypes with zero variance (correlation undefined).
        ct_var    <- apply(ct_sub, 1, stats::var, na.rm = TRUE)
        ct_sub    <- ct_sub[is.finite(ct_var) & ct_var > 0, , drop = FALSE]

        if (nrow(expr_sub) > 0 && nrow(ct_sub) > 0) {
          cor_mat <- cor(t(as.matrix(expr_sub)), t(ct_sub),
                         method = "pearson",
                         use    = "pairwise.complete.obs")
          write.csv(
            cor_mat,
            file.path(out_dir,
              paste0(sample_id,
                     "_nnSVG_top20_x_celltype_correlation.csv"))
          )

          cor_long <- data.frame(
            gene     = rep(rownames(cor_mat), ncol(cor_mat)),
            celltype = rep(colnames(cor_mat), each = nrow(cor_mat)),
            r        = as.vector(cor_mat)
          )
          cor_long$gene     <- factor(cor_long$gene,
                                      levels = rev(rownames(cor_mat)))
          cor_long$celltype <- factor(cor_long$celltype,
                                      levels = colnames(cor_mat))

          p_cor <- ggplot2::ggplot(cor_long,
                     ggplot2::aes(x = celltype, y = gene, fill = r)) +
            ggplot2::geom_tile(colour = "grey92", linewidth = 0.2) +
            ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white",
                                          high = "#b2182b", midpoint = 0,
                                          limits = c(-1, 1),
                                          name = "Pearson r") +
            ggplot2::labs(
              title    = sprintf(
                "%s \u2014 top-20 SVGs \u00D7 celltype pixel correlation",
                sample_id),
              subtitle = sprintf(
                "Resolution: %s coord units; pixels: %d",
                raster_resolution, ncol(expr_sub)),
              x = NULL, y = NULL
            ) +
            presentation_theme(base_size = 11, legend_position = "right",
                               x_angle = 45)
          save_presentation_plot(
            plot     = p_cor,
            filename = file.path(out_dir,
              paste0(sample_id,
                     "_nnSVG_top20_x_celltype_correlation.png")),
            width    = max(10, ncol(cor_mat) * 0.5 + 5),
            height   = max(8,  nrow(cor_mat) * 0.34 + 3),
            dpi      = 150
          )
          cat("  \u2713 P14 gene \u00D7 celltype correlation heatmap saved\n")
        }
      }
    }, error = function(e) {
      cat("\u26A0 nnSVG gene x celltype correlation failed: ",
          conditionMessage(e), "\n", sep = "")
    })

    # --- P15: focus-celltype side-by-side pixel maps ---
    # For each top-20 SVG gene, render two pixel maps: (a) gene pixel-mean
    # expression; (b) focus-celltype pixel count, both on the same grid.
    tryCatch({
      if (!is.null(ct_raster) && exists("spe_raster", inherits = FALSE) &&
          nnsvg_assay == "pixelval") {
        ct_mat    <- as.matrix(SummarizedExperiment::assay(ct_raster, 1))
        expr_mat  <- SummarizedExperiment::assay(spe_raster, "pixelval")
        focus_rows <- grep(focus_celltype, rownames(ct_mat),
                           value = TRUE, ignore.case = TRUE, perl = TRUE)
        if (length(focus_rows) == 0) {
          cat("  \u2139 P15 skipped: focus_celltype (", focus_celltype,
              ") matched no rows in celltype raster.\n", sep = "")
        } else {
          focus_label <- paste(focus_rows, collapse = " / ")
          focus_count <- colSums(ct_mat[focus_rows, , drop = FALSE])
          coords_mat  <- as.data.frame(
            SpatialExperiment::spatialCoords(spe_raster)
          )
          names(coords_mat)[1:2] <- c("px_x", "px_y")
          coords_mat$pixel_id <- rownames(coords_mat)
          if (is.null(coords_mat$pixel_id) || all(is.na(coords_mat$pixel_id))) {
            coords_mat$pixel_id <- colnames(spe_raster)
          }

          top20_p15 <- intersect(utils::head(rownames(result_df), 20L),
                                 rownames(expr_mat))
          if (length(top20_p15) > 0) {
            # SEraster's pixelval assay is a dgCMatrix; reshape2::melt
            # cannot coerce sparse matrices to data.frame, so densify the
            # top-20 x n_pixel subset (small; ~100k cells at 437 pixels).
            expr_mat_p15 <- as.matrix(expr_mat[top20_p15, , drop = FALSE])
            expr_long <- reshape2::melt(
              expr_mat_p15,
              varnames = c("gene", "pixel_id"),
              value.name = "value"
            )
            expr_long$layer <- "Gene expression"
            # Align focus count to pixel ids referenced in expr_long.
            fc_aligned <- focus_count[as.character(expr_long$pixel_id)]
            ct_long <- data.frame(
              gene     = expr_long$gene,
              pixel_id = expr_long$pixel_id,
              value    = as.numeric(fc_aligned),
              layer    = sprintf("%s count", focus_label)
            )
            pl_df <- rbind(
              data.frame(gene = expr_long$gene,
                         pixel_id = expr_long$pixel_id,
                         value = expr_long$value,
                         layer = expr_long$layer,
                         stringsAsFactors = FALSE),
              ct_long
            )
            pl_df <- merge(pl_df, coords_mat, by = "pixel_id", all.x = TRUE)
            pl_df$gene  <- factor(pl_df$gene, levels = top20_p15)
            pl_df$layer <- factor(pl_df$layer,
                                  levels = c("Gene expression",
                                             sprintf("%s count", focus_label)))

            p_sbs <- ggplot2::ggplot(pl_df,
                      ggplot2::aes(x = px_x, y = px_y, fill = value)) +
              ggplot2::geom_tile() +
              ggplot2::coord_fixed() +
              ggplot2::facet_grid(gene ~ layer, scales = "free") +
              ggplot2::scale_fill_viridis_c(option = "magma",
                                            na.value = "grey92") +
              ggplot2::labs(
                title = sprintf(
                  "%s \u2014 top-20 SVGs vs %s pixel counts",
                  sample_id, focus_label),
                subtitle = sprintf(
                  "Pixel grid at resolution = %s coord units",
                  raster_resolution),
                x = NULL, y = NULL
              ) +
              presentation_theme(base_size = 10, legend_position = "right") +
              ggplot2::theme(
                axis.text = ggplot2::element_blank(),
                axis.ticks = ggplot2::element_blank(),
                panel.grid = ggplot2::element_blank()
              )
            save_presentation_plot(
              plot     = p_sbs,
              filename = file.path(out_dir,
                paste0(sample_id,
                       "_nnSVG_top20_focus_pixel_sidebyside.png")),
              width    = 9,
              height   = min(24, 1.0 * length(top20_p15) + 2),
              dpi      = 120
            )
            cat("  \u2713 P15 focus-celltype side-by-side panel saved\n")
          }
        }
      }
    }, error = function(e) {
      cat("\u26A0 nnSVG focus-celltype side-by-side failed: ",
          conditionMessage(e), "\n", sep = "")
    })

    top_svgs <- head(rownames(result_df), n_top_svgs)
    cat("\u2713 nnSVG complete.\n")
    cat("  Top 10 SVGs:", paste(head(top_svgs, 10), collapse = ", "), "\n")
    cat("  Results saved to:", out_dir, "\n\n")
    
    attr(result_df, "bcell_nnsvg") <- bcell_nnsvg_df
    invisible(result_df)
  } else {
    invisible(bcell_nnsvg_df)
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

.maybe_cleanup_between_cci_sections <- function(enabled = TRUE, label = NULL) {
  if (!isTRUE(enabled)) {
    return(invisible(NULL))
  }
  
  cleanup_fn <- get0("cleanup_memory", mode = "function", inherits = TRUE)
  if (is.function(cleanup_fn)) {
    cleanup_fn(label = paste("CCI", label), verbose = TRUE)
  } else {
    gc(verbose = FALSE)
  }
}

.normalize_run_sections <- function(run_sections) {
  default_sections <- c(
    insitucor = TRUE,
    liana = TRUE,
    nichenet = FALSE,
    misty = TRUE,
    nnsvg = FALSE
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

.resolve_celltype_column_auto <- function(gobj, celltype_col = NULL) {
  if (!is.null(celltype_col) && nzchar(celltype_col)) return(celltype_col)
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  sup_cols <- grep("celltype_.*_supervised$", names(meta), value = TRUE)
  if (length(sup_cols) > 0) return(sup_cols[1])
  if ("leiden_clust" %in% names(meta)) return("leiden_clust")
  NA_character_
}

.resolve_nichenet_senders_auto <- function(gobj,
                                           output_dir,
                                           sample_id,
                                           celltype_col,
                                           receiver_celltype,
                                           min_cells = 5L,
                                           top_n     = 5L) {
  # Step-09 proximity enrichment, if available.
  prox_csv <- file.path(output_dir, "09_Spatial_Network", "proximity",
                        paste0(sample_id, "_proximity_enrichment.csv"))
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  celltype_col <- .resolve_celltype_column_auto(gobj, celltype_col)
  if (is.na(celltype_col) || !celltype_col %in% names(meta)) return(NULL)
  ct_counts <- table(meta[[celltype_col]])
  eligible  <- setdiff(names(ct_counts)[ct_counts >= min_cells], receiver_celltype)

  if (file.exists(prox_csv)) {
    prox <- tryCatch(read.csv(prox_csv, stringsAsFactors = FALSE),
                     error = function(e) NULL)
    if (!is.null(prox) && all(c("unified_int", "enrichm") %in% names(prox))) {
      # Rows involving the receiver; pick significant enrichm > 0.
      pair_re <- paste0("(^|--)", gsub("([.+*?^$()\\[\\]])", "\\\\\\1", receiver_celltype), "(--|$)")
      rows <- prox[grepl(pair_re, prox$unified_int), , drop = FALSE]
      padj <- if ("p.adj_higher" %in% names(rows)) rows$p.adj_higher else rep(NA_real_, nrow(rows))
      ok   <- !is.na(rows$enrichm) & rows$enrichm > 0 &
              (is.na(padj) | padj < 0.05)
      rows <- rows[ok, , drop = FALSE]
      rows <- rows[order(-rows$enrichm), , drop = FALSE]
      partners <- vapply(strsplit(as.character(rows$unified_int), "--", fixed = TRUE),
                         function(x) setdiff(x, receiver_celltype)[1],
                         character(1))
      partners <- partners[!is.na(partners) & partners %in% eligible]
      partners <- unique(partners)
      if (length(partners) > 0) return(head(partners, top_n))
    }
  }
  # Fallback: use the top_n most abundant non-receiver cell types.
  abundant <- names(sort(ct_counts[eligible], decreasing = TRUE))
  if (length(abundant) > 0) return(head(abundant, top_n))
  NULL
}

.resolve_nichenet_targets_auto <- function(gobj,
                                           output_dir,
                                           sample_id,
                                           celltype_col,
                                           receiver_celltype,
                                           top_n = 50L) {
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  celltype_col <- .resolve_celltype_column_auto(gobj, celltype_col)

  # First try: step 12 per-sample DE if it has run for this receiver.
  de_glob <- list.files(file.path(output_dir, "12_Spatial_DE"),
                        pattern = paste0("^", sample_id, "_.*sample_de\\.csv$"),
                        full.names = TRUE)
  rec_token <- gsub("[^A-Za-z0-9]+", ".", receiver_celltype)
  de_glob   <- de_glob[grepl(rec_token, basename(de_glob), ignore.case = TRUE)]
  if (length(de_glob) > 0) {
    de <- tryCatch(read.csv(de_glob[1], stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (!is.null(de) && "gene" %in% names(de)) {
      lfc_col <- intersect(c("logFC", "log2FC", "estimate"), names(de))[1]
      if (!is.na(lfc_col)) {
        de <- de[order(-abs(de[[lfc_col]])), , drop = FALSE]
        out <- head(unique(de$gene), top_n)
        if (length(out) > 0) return(out)
      }
    }
  }

  # Second try: step 06 cluster markers, filtered to the cluster that is
  # dominated by the receiver cell type.
  mk_csv <- file.path(output_dir, "06_Markers",
                      paste0(sample_id, "_all_markers.csv"))
  if (file.exists(mk_csv) && !is.na(celltype_col) && celltype_col %in% names(meta) &&
      "leiden_clust" %in% names(meta)) {
    tab <- table(meta[[celltype_col]], meta$leiden_clust)
    rec_row <- tab[rownames(tab) == receiver_celltype, , drop = FALSE]
    if (nrow(rec_row) > 0 && sum(rec_row) > 0) {
      # Clusters where receiver is the plurality celltype.
      cluster_frac <- rec_row[1, , drop = TRUE] / pmax(1, colSums(tab))
      receiver_clusters <- names(cluster_frac)[cluster_frac >= 0.5]
      if (length(receiver_clusters) > 0) {
        mk <- tryCatch(read.csv(mk_csv, stringsAsFactors = FALSE),
                       error = function(e) NULL)
        if (!is.null(mk)) {
          gene_col <- intersect(c("feats", "feature", "gene"), names(mk))[1]
          clust_col <- intersect(c("cluster", "cluster.1"), names(mk))[1]
          lfc_col   <- intersect(c("summary.logFC", "logFC", "avg_log2FC"), names(mk))[1]
          if (!is.na(gene_col) && !is.na(clust_col) && !is.na(lfc_col)) {
            mk <- mk[as.character(mk[[clust_col]]) %in% receiver_clusters, , drop = FALSE]
            mk <- mk[order(-abs(mk[[lfc_col]])), , drop = FALSE]
            out <- head(unique(mk[[gene_col]]), top_n)
            if (length(out) > 0) return(out)
          }
        }
      }
    }
  }

  character(0)
}

.nichenet_section_ready <- function(mode, sender_celltypes = NULL, receiver_celltype = NULL) {
  mode <- match.arg(mode, c("single", "all_senders_to_receiver", "all_pairs"))
  
  if (mode == "single") {
    if (is.null(sender_celltypes) || length(sender_celltypes) == 0 ||
        is.null(receiver_celltype) || !nzchar(as.character(receiver_celltype)[1])) {
      return(list(
        ok = FALSE,
        reason = "NicheNet mode 'single' requires sender_celltypes and receiver_celltype."
      ))
    }
  }
  
  if (mode == "all_senders_to_receiver" &&
      (is.null(receiver_celltype) || !nzchar(as.character(receiver_celltype)[1]))) {
    return(list(
      ok = FALSE,
      reason = "NicheNet mode 'all_senders_to_receiver' requires receiver_celltype."
    ))
  }
  
  list(ok = TRUE, reason = NULL)
}

.misty_runtime_ready <- function() {
  if (!requireNamespace("mistyR", quietly = TRUE)) {
    return(list(ok = FALSE, reason = "mistyR is not installed."))
  }

  ridge_ok <- tryCatch(requireNamespace("ridge", quietly = TRUE), error = function(e) FALSE)

  if (!isTRUE(ridge_ok)) {
    skip_msg <- paste0(
      "MISTy skipped: the 'ridge' package could not load (needs libgsl.so.27).\n",
      "  The container image must ship Ubuntu's libgsl27 package; conda-forge\n",
      "  GSL uses a different SONAME (libgsl.so.25 / .so.28) and will not work.\n",
      "  Fix: rebuild the .sif from an Apptainer .def that includes 'libgsl27'\n",
      "  in the apt install list of the %post section."
    )
    return(list(ok = FALSE, reason = skip_msg))
  }

  list(ok = TRUE, reason = NULL)
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
                             focus_celltype     = "^B\\.cell$",
                             sender_celltypes   = NULL,
                             receiver_celltype  = NULL,
                             target_genes       = NULL,
                             target_genes_by_receiver = NULL,
                             nichenet_network_dir = NULL,
                             nichenet_mode      = c("single", "all_senders_to_receiver", "all_pairs"),
                             nichenet_spatial_filter = FALSE,
                             nichenet_proximity_enrichment_path = NULL,
                             nichenet_spatial_padj_threshold = 0.05,
                             nichenet_min_cells_per_celltype = 5,
                             nichenet_include_self_pairs = FALSE,
                             run_sections       = c(insitucor = TRUE,
                                                    liana     = TRUE,
                                                    nichenet  = FALSE,
                                                    misty     = TRUE,
                                                    nnsvg     = FALSE),
                             nnsvg_raster_resolution = "auto",
                             cleanup_between_sections = TRUE,
                             sample_row         = NULL,
                             sample_sheet_path  = NULL) {
  
  cat("\n========================================\n")
  cat("STEP 10: CCI Analysis (Layer 2 & 3)\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  nichenet_mode <- match.arg(nichenet_mode)
  nichenet_spatial_option <- .normalize_nichenet_spatial_option(nichenet_spatial_filter)
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
  
  # FIX #6: Pre-extract expression matrices ONCE to avoid repeated dense
  # materializations across CCI sections. At ~200k cells x ~6k genes, each
  # call to getExpression() + as.matrix() allocates ~10 GB. Extracting once
  # and passing as expr_cache saves ~50% peak memory for step 10.
  cat("Pre-extracting expression matrices for CCI sections...\n")
  .cci_expr_cache <- list(
    raw = .giotto_get_expression(gobj, values = "raw", output = "matrix"),
    normalized = .giotto_get_expression(gobj, values = "normalized", output = "matrix")
  )
  cat("\u2713 Expression cache ready (raw:", paste(dim(.cci_expr_cache$raw), collapse = "x"),
      ", normalized:", paste(dim(.cci_expr_cache$normalized), collapse = "x"), ")\n\n")
  
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
      function() run_insitucor(gobj, sample_id, output_dir, celltype_col,
                               expr_cache        = .cci_expr_cache,
                               sample_row        = sample_row,
                               sample_sheet_path = sample_sheet_path)
    )
    results$insitucor <- section_out$result
    section_status["insitucor"] <- section_out$status
    section_messages[["insitucor"]] <- section_out$message
    .maybe_cleanup_between_cci_sections(cleanup_between_sections, "InSituCor")
  }
  
  if (isTRUE(run_sections["liana"])) {
    section_out <- .run_cci_section(
      "LIANA",
      function() run_liana(gobj, sample_id, output_dir, celltype_col,
                           expr_cache        = .cci_expr_cache,
                           focus_celltype    = focus_celltype,
                           sample_row        = sample_row,
                           sample_sheet_path = sample_sheet_path)
    )
    results$liana <- section_out$result
    section_status["liana"] <- section_out$status
    section_messages[["liana"]] <- section_out$message
    .maybe_cleanup_between_cci_sections(cleanup_between_sections, "LIANA")
  }
  
  if (isTRUE(run_sections["nichenet"])) {
    # Auto-wire senders / targets from step 09 proximity and step 06/12 DE
    # when the user leaves the config keys null.
    resolved_celltype_col <- .resolve_celltype_column_auto(gobj, celltype_col)
    if (is.null(sender_celltypes) || length(sender_celltypes) == 0) {
      auto_senders <- .resolve_nichenet_senders_auto(
        gobj, output_dir, sample_id,
        celltype_col      = resolved_celltype_col,
        receiver_celltype = receiver_celltype,
        min_cells         = nichenet_min_cells_per_celltype
      )
      if (!is.null(auto_senders) && length(auto_senders) > 0) {
        sender_celltypes <- auto_senders
        cat("  \u2139 NicheNet senders (auto): ",
            paste(sender_celltypes, collapse = ", "), "\n", sep = "")
      }
    }
    if ((is.null(target_genes) || length(target_genes) == 0) &&
        !is.null(receiver_celltype) && nzchar(receiver_celltype)) {
      auto_targets <- .resolve_nichenet_targets_auto(
        gobj, output_dir, sample_id,
        celltype_col      = resolved_celltype_col,
        receiver_celltype = receiver_celltype
      )
      if (length(auto_targets) > 0) {
        target_genes <- auto_targets
        cat("  \u2139 NicheNet targets (auto): ", length(target_genes),
            " genes derived from step 06/12 for receiver '",
            receiver_celltype, "'\n", sep = "")
      }
    }
    nichenet_ready <- .nichenet_section_ready(nichenet_mode, sender_celltypes, receiver_celltype)
    if (!isTRUE(nichenet_ready$ok)) {
      cat("\u26A0 Skipping NicheNet:", nichenet_ready$reason, "\n")
      section_status["nichenet"] <- "skipped"
      section_messages[["nichenet"]] <- nichenet_ready$reason
    } else {
      section_out <- .run_cci_section(
        "NicheNet",
        function() {
          if (identical(nichenet_spatial_option, "both")) {
            list(
              unfiltered = run_nichenet_batch(
                gobj,
                sample_id,
                output_dir,
                celltype_col,
                sender_celltypes,
                receiver_celltype,
                target_genes = target_genes,
                target_genes_by_receiver = target_genes_by_receiver,
                network_dir = nichenet_network_dir,
                mode = nichenet_mode,
                use_spatial_filter = FALSE,
                proximity_enrichment_path = nichenet_proximity_enrichment_path,
                spatial_padj_threshold = nichenet_spatial_padj_threshold,
                min_cells_per_celltype = nichenet_min_cells_per_celltype,
                include_self_pairs = nichenet_include_self_pairs,
                result_root_dir = file.path(output_dir, "10_CCI_Analysis", "nichenet", "unfiltered"),
                expr_cache = .cci_expr_cache
              ),
              spatial_filtered = run_nichenet_batch(
                gobj,
                sample_id,
                output_dir,
                celltype_col,
                sender_celltypes,
                receiver_celltype,
                target_genes = target_genes,
                target_genes_by_receiver = target_genes_by_receiver,
                network_dir = nichenet_network_dir,
                mode = nichenet_mode,
                use_spatial_filter = TRUE,
                proximity_enrichment_path = nichenet_proximity_enrichment_path,
                spatial_padj_threshold = nichenet_spatial_padj_threshold,
                min_cells_per_celltype = nichenet_min_cells_per_celltype,
                include_self_pairs = nichenet_include_self_pairs,
                result_root_dir = file.path(output_dir, "10_CCI_Analysis", "nichenet", "spatial_filtered"),
                expr_cache = .cci_expr_cache
              )
            )
          } else {
            run_nichenet_batch(
              gobj,
              sample_id,
              output_dir,
              celltype_col,
              sender_celltypes,
              receiver_celltype,
              target_genes = target_genes,
              target_genes_by_receiver = target_genes_by_receiver,
              network_dir = nichenet_network_dir,
              mode = nichenet_mode,
              use_spatial_filter = identical(nichenet_spatial_option, "filtered"),
              proximity_enrichment_path = nichenet_proximity_enrichment_path,
              spatial_padj_threshold = nichenet_spatial_padj_threshold,
              min_cells_per_celltype = nichenet_min_cells_per_celltype,
              include_self_pairs = nichenet_include_self_pairs,
              result_root_dir = file.path(output_dir, "10_CCI_Analysis", "nichenet"),
              expr_cache = .cci_expr_cache
            )
          }
        }
      )
      results$nichenet <- section_out$result
      section_status["nichenet"] <- section_out$status
      section_messages[["nichenet"]] <- section_out$message
      .maybe_cleanup_between_cci_sections(cleanup_between_sections, "NicheNet")
    }
  }
  
  if (isTRUE(run_sections["misty"])) {
    misty_ready <- .misty_runtime_ready()
    if (!isTRUE(misty_ready$ok)) {
      cat("\n", strrep("=", 72), "\n", sep = "")
      cat("  \u26A0  MISTy SKIPPED \u2014 sample: ", sample_id, "\n", sep = "")
      cat(strrep("=", 72), "\n")
      cat(misty_ready$reason, "\n")
      cat("  Recommended substitutes for cell-cell interaction inference:\n")
      cat("    - InSituCor (spatial co-expression modules)\n")
      cat("    - LIANA (L-R consensus)\n")
      cat(strrep("=", 72), "\n\n")
      section_status["misty"] <- "skipped"
      section_messages[["misty"]] <- misty_ready$reason

      # Persist skip record so analysts can audit post-hoc without grepping logs.
      tryCatch({
        misty_dir <- file.path(output_dir, "10_CCI_Analysis", "misty")
        dir.create(misty_dir, recursive = TRUE, showWarnings = FALSE)
        skip_path <- file.path(misty_dir, paste0(sample_id, "_misty_skipped.txt"))
        writeLines(
          c(
            paste0("sample_id: ", sample_id),
            paste0("timestamp: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
            "status: skipped",
            "reason:",
            misty_ready$reason,
            "",
            "recommended_substitutes: InSituCor, LIANA"
          ),
          skip_path
        )
      }, error = function(e) {
        cat("\u26A0 Could not write misty_skipped.txt:", conditionMessage(e), "\n")
      })
    } else {
      section_out <- .run_cci_section(
        "MISTy",
        function() run_misty(gobj, sample_id, output_dir,
                             expr_cache     = .cci_expr_cache,
                             celltype_col   = celltype_col,
                             focus_celltype = focus_celltype)
      )
      results$misty <- section_out$result
      section_status["misty"] <- section_out$status
      section_messages[["misty"]] <- section_out$message
      .maybe_cleanup_between_cci_sections(cleanup_between_sections, "MISTy")
    }
  }
  
  if (isTRUE(run_sections["nnsvg"])) {
    section_out <- .run_cci_section(
      "nnSVG",
      function() run_nnsvg(gobj, sample_id, output_dir,
                           expr_cache        = .cci_expr_cache,
                           celltype_col      = celltype_col,
                           focus_celltype    = focus_celltype,
                           raster_resolution = nnsvg_raster_resolution)
    )
    results$nnsvg <- section_out$result
    section_status["nnsvg"] <- section_out$status
    section_messages[["nnsvg"]] <- section_out$message
    .maybe_cleanup_between_cci_sections(cleanup_between_sections, "nnSVG")
  }
  
  # Release the expression cache now that all sections are done.
  rm(.cci_expr_cache)
  gc(verbose = FALSE)
  
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
      skipped = "-",
      disabled = "-",
      "-"
    )
    cat("  ", icon, " ", section_name, ": ", section_status[[section_name]], "\n", sep = "")
    if (!is.null(section_messages[[section_name]]) && nzchar(section_messages[[section_name]])) {
      cat("      ", section_messages[[section_name]], "\n", sep = "")
    }
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
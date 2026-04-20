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
                          celltype_col = NULL,
                          n_cores      = 4,
                          expr_cache   = NULL) {
  
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
    cat("\u2713 InSituCor complete. Results saved to:", out_dir, "\n")
  }
  
  invisible(cor_results)
}


# ==============================================================================
# SECTION 1 helper — Extended LIANA visualizations
# ==============================================================================

plot_liana_extended <- function(liana_agg,
                                meta,
                                gobj,
                                expr_cache,
                                out_dir,
                                sample_id,
                                celltype_col,
                                top_n = 20) {

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

  # -- Plot 1: LR ranking bar chart ------------------------------------------
  tryCatch({
    top_lr <- head(agg[order(agg$aggregate_rank, na.last = TRUE), ], top_n)
    top_lr$interaction_label  <- paste(top_lr$ligand_complex, "\u2192", top_lr$receptor_complex)
    top_lr$neg_log10_rank     <- -log10(pmax(top_lr$aggregate_rank, 1e-10))

    p1 <- ggplot2::ggplot(top_lr,
             ggplot2::aes(
               x    = reorder(interaction_label, neg_log10_rank),
               y    = neg_log10_rank,
               fill = source
             )) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title    = sample_plot_title(sample_id, "Top Ligand-Receptor Interactions"),
        subtitle = paste0("Top ", top_n,
                          " interactions by LIANA aggregate rank (\u2212log\u2081\u2080; taller = more significant)"),
        x        = "Interaction (Ligand \u2192 Receptor)",
        y        = "\u2212log\u2081\u2080(Aggregate Rank)",
        fill     = "Sender"
      ) +
      presentation_theme(base_size = 12, legend_position = "right")

    save_presentation_plot(
      plot     = p1,
      filename = file.path(out_dir, paste0(sample_id, "_liana_lr_ranking.png")),
      width    = 12,
      height   = 8,
      dpi      = 150
    )
    cat("\u2713 LIANA LR ranking plot saved\n")
  }, error = function(e) {
    cat("\u26A0 LIANA LR ranking plot failed:", conditionMessage(e), "\n")
  })

  # -- Plot 2: CCI heatmap ---------------------------------------------------
  tryCatch({
    if (is.null(count_df) || nrow(count_df) == 0) stop("No interaction count data.")

    p2 <- ggplot2::ggplot(count_df,
             ggplot2::aes(x = source, y = target, fill = n_interactions)) +
      ggplot2::geom_tile(color = "white") +
      ggplot2::geom_text(ggplot2::aes(label = n_interactions), color = "white", size = 3) +
      ggplot2::scale_fill_viridis_c(option = "plasma", name = "# Interactions") +
      ggplot2::labs(
        title    = sample_plot_title(sample_id, "Cell-Cell Interaction Heatmap"),
        subtitle = "Number of significant L-R pairs per sender-receiver combination",
        x        = "Sender",
        y        = "Receiver"
      ) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      presentation_theme(base_size = 12, legend_position = "right")

    save_presentation_plot(
      plot     = p2,
      filename = file.path(out_dir, paste0(sample_id, "_liana_cci_heatmap.png")),
      width    = 10,
      height   = 9,
      dpi      = 150
    )
    cat("\u2713 LIANA CCI heatmap saved\n")
  }, error = function(e) {
    cat("\u26A0 LIANA CCI heatmap failed:", conditionMessage(e), "\n")
  })

  # -- Plot 3: CCI network graph (requires ggraph) ---------------------------
  tryCatch({
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      cat("\u26A0 LIANA CCI network skipped: 'ggraph' is not installed.\n")
    } else if (is.null(count_df) || nrow(count_df) == 0) {
      stop("No interaction count data.")
    } else {
      g <- igraph::graph_from_data_frame(count_df, directed = TRUE)

      p3 <- ggraph::ggraph(g, layout = "circle") +
        ggraph::geom_edge_arc(
          ggplot2::aes(width = n_interactions),
          arrow  = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
          end_cap = ggraph::circle(4, "mm"),
          alpha  = 0.7,
          color  = "steelblue"
        ) +
        ggraph::geom_node_label(ggplot2::aes(label = name), size = 3) +
        ggraph::scale_edge_width_continuous(range = c(0.5, 3), name = "# Interactions") +
        ggplot2::labs(
          title    = sample_plot_title(sample_id, "CCI Network"),
          subtitle = "Edge width proportional to number of significant L-R pairs"
        ) +
        presentation_theme(base_size = 12, legend_position = "right")

      save_presentation_plot(
        plot     = p3,
        filename = file.path(out_dir, paste0(sample_id, "_liana_cci_network.png")),
        width    = 10,
        height   = 10,
        dpi      = 150
      )
      cat("\u2713 LIANA CCI network saved\n")
    }
  }, error = function(e) {
    cat("\u26A0 LIANA CCI network failed:", conditionMessage(e), "\n")
  })

  # -- Plot 4: Chord diagram (requires circlize) -----------------------------
  tryCatch({
    if (!requireNamespace("circlize", quietly = TRUE)) {
      cat("\u26A0 LIANA chord diagram skipped: 'circlize' is not installed.\n")
    } else if (is.null(count_df) || nrow(count_df) == 0) {
      stop("No interaction count data.")
    } else {
      mat <- stats::xtabs(n_interactions ~ source + target, data = count_df)
      chord_path <- file.path(out_dir, paste0(sample_id, "_liana_chord.png"))
      grDevices::png(chord_path, width = 2400, height = 2400, res = 200)
      circlize::chordDiagram(
        mat,
        transparency     = 0.3,
        annotationTrack  = c("grid", "name"),
        preAllocateTracks = 1
      )
      graphics::title(main = paste0(sample_id, " \u2014 CCI Chord Diagram"))
      grDevices::dev.off()
      cat("\u2713 LIANA chord diagram saved\n")
    }
  }, error = function(e) {
    tryCatch(grDevices::dev.off(), error = function(e2) NULL)
    cat("\u26A0 LIANA chord diagram failed:", conditionMessage(e), "\n")
  })

  # -- Plot 5: Spatial LR expression map ------------------------------------
  tryCatch({
    spat <- tryCatch(
      as.data.frame(.giotto_get_spatial_locations(gobj, output = "data.table")),
      error = function(e) NULL
    )
    if (is.null(spat) || !all(c("cell_ID", "sdimx", "sdimy") %in% names(spat))) {
      stop("Spatial locations unavailable.")
    }

    top_pair <- head(agg[order(agg$aggregate_rank, na.last = TRUE), ], 1)
    if (nrow(top_pair) == 0) stop("No interactions to map.")

    ligand_gene   <- strsplit(as.character(top_pair$ligand_complex[1]),   "_")[[1]][1]
    receptor_gene <- strsplit(as.character(top_pair$receptor_complex[1]), "_")[[1]][1]
    sender_ct     <- as.character(top_pair$source[1])
    receiver_ct   <- as.character(top_pair$target[1])

    expr_mat <- if (!is.null(expr_cache) && !is.null(expr_cache$normalized)) {
      expr_cache$normalized
    } else {
      tryCatch(
        .giotto_get_expression(gobj, values = "normalized", output = "matrix"),
        error = function(e) NULL
      )
    }
    if (is.null(expr_mat)) stop("Expression matrix unavailable.")

    sender_ids   <- meta$cell_ID[!is.na(meta[[celltype_col]]) &
                                    meta[[celltype_col]] == sender_ct]
    receiver_ids <- meta$cell_ID[!is.na(meta[[celltype_col]]) &
                                    meta[[celltype_col]] == receiver_ct]

    extract_expr_vec <- function(gene, cell_ids) {
      ids <- intersect(as.character(cell_ids), colnames(expr_mat))
      if (length(ids) == 0 || !gene %in% rownames(expr_mat)) return(NULL)
      vals <- expr_mat[gene, ids, drop = TRUE]
      data.frame(cell_ID = ids, expression = as.numeric(vals), stringsAsFactors = FALSE)
    }

    lig_df <- extract_expr_vec(ligand_gene,   sender_ids)
    rec_df <- extract_expr_vec(receptor_gene, receiver_ids)
    if (is.null(lig_df) && is.null(rec_df)) stop("No expression data for top L-R pair.")

    panels <- Filter(Negate(is.null), list(
      if (!is.null(lig_df)) {
        lig_df$panel <- paste0("Ligand: ", ligand_gene, " (", sender_ct, ")"); lig_df
      },
      if (!is.null(rec_df)) {
        rec_df$panel <- paste0("Receptor: ", receptor_gene, " (", receiver_ct, ")"); rec_df
      }
    ))
    panel_df <- merge(do.call(rbind, panels),
                      spat[, c("cell_ID", "sdimx", "sdimy")],
                      by = "cell_ID")

    p5 <- ggplot2::ggplot(panel_df,
             ggplot2::aes(x = sdimx, y = sdimy, color = expression)) +
      ggplot2::geom_point(size = 0.8, alpha = 0.8) +
      ggplot2::facet_wrap(~panel) +
      ggplot2::scale_color_viridis_c(option = "magma", name = "Expression") +
      ggplot2::labs(
        title    = sample_plot_title(sample_id,
                     paste0("Spatial Expression: ", ligand_gene, " \u2192 ", receptor_gene)),
        subtitle = "Normalized expression of top-ranked L-R pair in sender/receiver cells",
        x        = "x",
        y        = "y"
      ) +
      presentation_theme(base_size = 12, legend_position = "right")

    save_presentation_plot(
      plot     = p5,
      filename = file.path(out_dir, paste0(sample_id, "_liana_spatial_lr.png")),
      width    = 14,
      height   = 7,
      dpi      = 150
    )
    cat("\u2713 LIANA spatial LR map saved\n")
  }, error = function(e) {
    cat("\u26A0 LIANA spatial LR map failed:", conditionMessage(e), "\n")
  })

  # -- Plot 6: Spatial CCI map -----------------------------------------------
  tryCatch({
    spat <- tryCatch(
      as.data.frame(.giotto_get_spatial_locations(gobj, output = "data.table")),
      error = function(e) NULL
    )
    if (is.null(spat) || !all(c("cell_ID", "sdimx", "sdimy") %in% names(spat))) {
      stop("Spatial locations unavailable.")
    }
    if (is.null(count_df) || nrow(count_df) == 0) stop("No interaction count data.")

    cell_spat <- merge(
      spat[, c("cell_ID", "sdimx", "sdimy")],
      meta[, c("cell_ID", celltype_col)],
      by = "cell_ID"
    )
    names(cell_spat)[names(cell_spat) == celltype_col] <- "celltype"

    centroids <- stats::aggregate(
      cbind(sdimx, sdimy) ~ celltype,
      data  = cell_spat,
      FUN   = mean,
      na.rm = TRUE
    )

    top_pairs <- head(count_df[order(-count_df$n_interactions), ], 10)
    seg_df <- merge(top_pairs,  centroids, by.x = "source", by.y = "celltype")
    names(seg_df)[names(seg_df) %in% c("sdimx", "sdimy")] <- c("x_start", "y_start")
    seg_df <- merge(seg_df, centroids, by.x = "target", by.y = "celltype")
    names(seg_df)[names(seg_df) %in% c("sdimx", "sdimy")] <- c("x_end", "y_end")

    p6 <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = cell_spat,
        ggplot2::aes(x = sdimx, y = sdimy, color = celltype),
        size = 0.5, alpha = 0.5
      ) +
      ggplot2::geom_segment(
        data = seg_df,
        ggplot2::aes(
          x = x_start, y = y_start, xend = x_end, yend = y_end,
          linewidth = n_interactions
        ),
        color = "black", alpha = 0.7,
        arrow = grid::arrow(length = grid::unit(3, "mm"), type = "open")
      ) +
      ggplot2::scale_linewidth_continuous(range = c(0.5, 3), name = "# Interactions") +
      ggplot2::labs(
        title    = sample_plot_title(sample_id, "Spatial CCI Map"),
        subtitle = "Arrows connect cell-type centroids; width proportional to # significant L-R pairs",
        x        = "x",
        y        = "y",
        color    = "Cell Type"
      ) +
      presentation_theme(base_size = 12, legend_position = "right")

    save_presentation_plot(
      plot     = p6,
      filename = file.path(out_dir, paste0(sample_id, "_liana_spatial_cci_map.png")),
      width    = 12,
      height   = 10,
      dpi      = 150
    )
    cat("\u2713 LIANA spatial CCI map saved\n")
  }, error = function(e) {
    cat("\u26A0 LIANA spatial CCI map failed:", conditionMessage(e), "\n")
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
                      celltype_col = NULL,
                      methods      = NULL,
                      expr_cache   = NULL) {
  
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
      
      top_n_groups  <- 15L
      src_col_dp    <- .first_existing_column(liana_agg, c("source", "sender", "Sender"))
      tgt_col_dp    <- .first_existing_column(liana_agg, c("target", "receiver", "Receiver"))
      top_sources   <- names(sort(table(liana_agg[[src_col_dp]]), decreasing = TRUE))[
        seq_len(min(top_n_groups, length(unique(liana_agg[[src_col_dp]]))))]
      top_targets   <- names(sort(table(liana_agg[[tgt_col_dp]]), decreasing = TRUE))[
        seq_len(min(top_n_groups, length(unique(liana_agg[[tgt_col_dp]]))))]

      p <- liana::liana_dotplot(
        liana_agg,
        source_groups = top_sources,
        target_groups = top_targets,
        ntop          = 20,
        specificity   = specificity_col,
        magnitude     = magnitude_col
      ) +
        ggplot2::labs(
          title    = sample_plot_title(sample_id, "Top LIANA Ligand-Receptor Interactions"),
          subtitle = paste0("Top ", top_n_groups, " senders/receivers by interaction count.")
        ) +
        presentation_theme(base_size = 11, legend_position = "right") +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = ggplot2::element_text(size = 7)
        )
      save_presentation_plot(
        plot     = p,
        filename = file.path(out_dir, paste0(sample_id, "_liana_dotplot.png")),
        width    = 22,
        height   = 14,
        dpi      = 150
      )
      cat("\u2713 LIANA dotplot saved\n")
    }, error = function(e) {
      cat("\u26A0 LIANA dotplot failed:", conditionMessage(e), "\n")
    })

    plot_liana_extended(
      liana_agg    = liana_agg,
      meta         = meta,
      gobj         = gobj,
      expr_cache   = expr_cache,
      out_dir      = out_dir,
      sample_id    = sample_id,
      celltype_col = celltype_col
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
                      expr_cache    = NULL) {
  
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
#' @param expr_cache    Optional pre-extracted expression cache list.
#' @return nnSVG results data frame

run_nnsvg <- function(gobj,
                      sample_id,
                      output_dir,
                      n_top_svgs = 100,
                      n_cores    = 4,
                      expr_cache = NULL) {
  
  if (!requireNamespace("nnSVG", quietly = TRUE))
    stop("nnSVG not installed.\n",
         'Install: BiocManager::install("nnSVG")')
  
  cat("\n--- nnSVG: Spatially variable gene detection ---\n")
  
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
  
  # FIX #7: Guard against running nnSVG on datasets too large for GP-based
  # SVG detection. nnSVG was designed for Visium-scale data (~5k spots); at
  # >50k cells the nearest-neighbor GP approximation becomes prohibitively
  # slow (hours to days).
  MAX_CELLS_NNSVG <- 50000
  if (ncol(spe) > MAX_CELLS_NNSVG) {
    cat("\u26A0 nnSVG: Dataset has", ncol(spe), "cells, which exceeds the recommended",
        MAX_CELLS_NNSVG, "cell limit for GP-based SVG detection.\n")
    cat("  Downsampling to", MAX_CELLS_NNSVG, "cells (random subsample).\n")
    set.seed(42)
    keep_idx <- sample(ncol(spe), MAX_CELLS_NNSVG)
    spe <- spe[, keep_idx]
    cat("  Cells after downsampling:", ncol(spe), "\n\n")
  }
  
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
    # ── Auto-detection: try to locate libgsl on common cluster / conda paths ──
    # The pipeline runs inside an Apptainer image so system paths inside the
    # container are checked first, then the conda env that supplies Python.
    conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
    conda_python  <- Sys.getenv("COSMX_PYTHON_PATH", unset = "")   # set by Run_Giotto_Pipeline.sh
    # Derive conda env lib dir from COSMX_PYTHON_PATH if CONDA_PREFIX is unset
    conda_lib_from_python <- if (nzchar(conda_python)) {
      file.path(dirname(dirname(conda_python)), "lib")
    } else ""

    gsl_candidates <- unique(c(
      "/usr/lib", "/usr/lib/x86_64-linux-gnu", "/usr/local/lib",
      conda_lib_from_python,
      if (nzchar(conda_prefix)) conda_prefix else "",
      if (nzchar(conda_prefix)) file.path(conda_prefix, "lib") else ""
    ))
    gsl_found <- Filter(function(p) {
      nzchar(p) && file.exists(p) &&
        any(grepl("libgsl\\.so", list.files(p, full.names = FALSE)))
    }, gsl_candidates)

    if (length(gsl_found) > 0) {
      existing_ldpath <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
      Sys.setenv(LD_LIBRARY_PATH = paste(
        c(gsl_found, existing_ldpath)[nzchar(c(gsl_found, existing_ldpath))],
        collapse = ":"
      ))
      message("  \u2139 Found libgsl in: ", paste(gsl_found, collapse = ", "),
              " \u2014 updated LD_LIBRARY_PATH, retrying ridge load...")
      ridge_ok <- tryCatch(requireNamespace("ridge", quietly = TRUE), error = function(e) FALSE)
    }
  }

  if (!isTRUE(ridge_ok)) {
    skip_msg <- paste0(
      "MISTy skipped: the 'ridge' package requires libgsl.so.27 (GNU Scientific Library).\n",
      "  To fix (no sudo required):\n",
      "    Option 1 \u2014 conda env used by this pipeline:\n",
      "        conda activate giotto_Py_3_11\n",
      "        conda install -c conda-forge gsl\n",
      "      then add to Run_Giotto_Pipeline.sh:\n",
      "        --env LD_LIBRARY_PATH=\"$CONDA_PREFIX/lib\"\n",
      "    Option 2 \u2014 bind host GSL libs into Apptainer:\n",
      "        --bind /path/to/gsl/lib:/usr/local/lib:ro\n",
      "      (add to the 'apptainer exec' call in Run_Giotto_Pipeline.sh)\n",
      "    Option 3 \u2014 rebuild the .sif image with GSL included\n",
      "    Option 4 \u2014 ask sysadmin to install libgsl27 / libgsl-dev on the host\n",
      "  MISTy outputs will be absent from this run."
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
                             cleanup_between_sections = TRUE) {
  
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
                               expr_cache = .cci_expr_cache)
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
                           expr_cache = .cci_expr_cache)
    )
    results$liana <- section_out$result
    section_status["liana"] <- section_out$status
    section_messages[["liana"]] <- section_out$message
    .maybe_cleanup_between_cci_sections(cleanup_between_sections, "LIANA")
  }
  
  if (isTRUE(run_sections["nichenet"])) {
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
      cat("\u26A0 Skipping MISTy:", misty_ready$reason, "\n")
      section_status["misty"] <- "skipped"
      section_messages[["misty"]] <- misty_ready$reason
    } else {
      section_out <- .run_cci_section(
        "MISTy",
        function() run_misty(gobj, sample_id, output_dir,
                             expr_cache = .cci_expr_cache)
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
                           expr_cache = .cci_expr_cache)
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
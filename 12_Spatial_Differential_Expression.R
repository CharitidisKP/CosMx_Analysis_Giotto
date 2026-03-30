#!/usr/bin/env Rscript
# ==============================================================================
# 12_Spatial_Differential_Expression.R
#
# Spatially resolved differential expression for CosMx / Giotto objects.
#
# Supported scopes:
# - sample: single-sample niche DE within one sample, comparing the same cell
#   type across spatial niches using FOV-level pseudobulk replicates
# - merged: cohort-level niche DE on merged objects, comparing treatments
#   within the same cell type and niche across samples
#
# The preferred sample-mode backend is smiDE, which is designed for spatially
# resolved DE with explicit attention to segmentation bias and spatial
# correlation. edgeR pseudobulk remains available as a fallback backend and is
# retained for merged multi-sample contrasts.
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
if (!exists("save_giotto_checkpoint") && file.exists(pipeline_utils)) {
  source(pipeline_utils)
}

if (!exists("%||%")) {
  `%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0) {
      return(y)
    }
    if (length(x) == 1 && is.atomic(x) && is.na(x)) {
      return(y)
    }
    x
  }
}

inspect_smide_namespace <- function(verbose = TRUE) {
  if (!requireNamespace("smiDE", quietly = TRUE)) {
    stop("The smiDE package is not installed in the current R library.")
  }
  
  exports <- tryCatch(sort(getNamespaceExports("smiDE")), error = function(e) character())
  extdata_dir <- system.file("extdata", package = "smiDE")
  extdata_files <- if (nzchar(extdata_dir) && dir.exists(extdata_dir)) {
    sort(list.files(extdata_dir))
  } else {
    character()
  }
  
  info <- list(
    version = as.character(utils::packageVersion("smiDE")),
    exports = exports,
    extdata_dir = extdata_dir,
    extdata_files = extdata_files,
    has_spatial_model = "spatial_model" %in% exports,
    has_overlap_ratio_metric = "overlap_ratio_metric" %in% exports
  )
  
  if (verbose) {
    cat("\n=== smiDE namespace inspection ===\n")
    cat("Version:", info$version, "\n")
    cat("Exports:", if (length(exports) > 0) paste(exports, collapse = ", ") else "<none>", "\n")
    cat("Extdata:", if (length(extdata_files) > 0) paste(extdata_files, collapse = ", ") else "<none>", "\n\n")
  }
  
  invisible(info)
}

align_metadata_with_spatial_locations <- function(metadata, spatial_locations) {
  if (all(c("sdimx", "sdimy") %in% names(metadata))) {
    return(metadata)
  }
  
  loc_cols <- c("cell_ID", "sdimx", "sdimy")
  if (!all(loc_cols %in% names(spatial_locations))) {
    stop("Spatial locations must contain cell_ID, sdimx, and sdimy.")
  }
  
  merge(
    metadata,
    spatial_locations[, loc_cols],
    by = "cell_ID",
    all.x = TRUE,
    sort = FALSE
  )
}

get_counts_cells_by_genes <- function(expr_mat, metadata) {
  keep_ids <- metadata$cell_ID[metadata$cell_ID %in% colnames(expr_mat)]
  if (length(keep_ids) == 0) {
    stop("No cell IDs from metadata were found in the expression matrix.")
  }
  
  counts <- t(as.matrix(expr_mat[, keep_ids, drop = FALSE]))
  storage.mode(counts) <- "numeric"
  counts
}

coerce_overlap_metric_table <- function(overlap_obj) {
  if (is.null(overlap_obj)) {
    return(NULL)
  }
  
  if (is.numeric(overlap_obj) && !is.null(names(overlap_obj))) {
    return(data.frame(
      gene = names(overlap_obj),
      overlap_metric = as.numeric(overlap_obj),
      stringsAsFactors = FALSE
    ))
  }
  
  if (is.matrix(overlap_obj)) {
    overlap_obj <- as.data.frame(overlap_obj, stringsAsFactors = FALSE)
  }
  
  if (is.data.frame(overlap_obj)) {
    overlap_df <- as.data.frame(overlap_obj, stringsAsFactors = FALSE)
    if (!"gene" %in% names(overlap_df)) {
      overlap_df$gene <- rownames(overlap_df)
    }
    rownames(overlap_df) <- NULL
    return(overlap_df)
  }
  
  NULL
}

extract_overlap_targets <- function(overlap_df, threshold = NULL) {
  if (is.null(overlap_df) || is.null(threshold)) {
    return(NULL)
  }
  
  numeric_cols <- names(overlap_df)[vapply(overlap_df, is.numeric, logical(1))]
  if (length(numeric_cols) == 0) {
    return(NULL)
  }
  
  metric_col <- numeric_cols[1]
  keep <- overlap_df[[metric_col]] <= threshold
  targets <- overlap_df$gene[keep & !is.na(overlap_df$gene)]
  unique(targets[nzchar(targets)])
}

coerce_smide_results_table <- function(res_obj) {
  if (is.null(res_obj)) {
    return(NULL)
  }
  
  if (is.data.frame(res_obj)) {
    return(as.data.frame(res_obj, stringsAsFactors = FALSE))
  }
  
  if (is.matrix(res_obj)) {
    return(as.data.frame(res_obj, stringsAsFactors = FALSE))
  }
  
  if (is.list(res_obj)) {
    if ("results" %in% names(res_obj)) {
      return(coerce_smide_results_table(res_obj$results))
    }
    
    df_items <- vapply(res_obj, function(x) is.data.frame(x) || is.matrix(x), logical(1))
    if (any(df_items)) {
      out <- lapply(names(res_obj)[df_items], function(nm) {
        df <- as.data.frame(res_obj[[nm]], stringsAsFactors = FALSE)
        df$result_component <- nm
        df
      })
      return(do.call(rbind, out))
    }
  }
  
  NULL
}

guess_significance_column <- function(result_df) {
  if (is.null(result_df) || nrow(result_df) == 0) {
    return(NULL)
  }
  
  candidates <- c(
    grep("^fdr$|adj.*p|padj|qval", names(result_df), value = TRUE, ignore.case = TRUE),
    grep("^p.value$|^p_val$|^pval$|p.value|p_val|pval", names(result_df), value = TRUE, ignore.case = TRUE)
  )
  candidates <- unique(candidates[candidates %in% names(result_df)])
  if (length(candidates) == 0) {
    return(NULL)
  }
  candidates[1]
}

run_smide_sample_backend <- function(expr_mat,
                                     metadata,
                                     run_label,
                                     annotation_column,
                                     tables_dir,
                                     sample_replicate_column = "fov",
                                     sample_contrast = c("one_vs_rest", "pairwise"),
                                     min_cells_per_niche = 30,
                                     spatial_locations = NULL,
                                     smide_family = "nbinom2",
                                     smide_radius = 0.05,
                                     smide_ncores = 1,
                                     smide_overlap_threshold = NULL,
                                     smide_save_raw = TRUE) {
  sample_contrast <- match.arg(sample_contrast)
  
  if (!requireNamespace("smiDE", quietly = TRUE)) {
    stop("smiDE is not installed, but backend = 'smiDE' was requested.")
  }
  
  metadata <- align_metadata_with_spatial_locations(metadata, spatial_locations)
  counts_cells_genes <- get_counts_cells_by_genes(expr_mat, metadata)
  metadata <- metadata[match(rownames(counts_cells_genes), metadata$cell_ID), , drop = FALSE]
  
  sample_ids <- unique(stats::na.omit(metadata$sample_id))
  if (length(sample_ids) > 1) {
    stop("smiDE sample-mode analysis expects exactly one sample in the input object.")
  }
  
  replicate_column <- detect_replicate_column(metadata, sample_replicate_column)
  split_col <- if (!is.null(replicate_column) && length(unique(stats::na.omit(metadata[[replicate_column]]))) > 1) {
    replicate_column
  } else {
    NULL
  }
  
  utils::write.csv(
    metadata,
    file.path(tables_dir, paste0(run_label, "_cell_metadata_with_niches.csv")),
    row.names = FALSE
  )
  
  summary_rows <- list()
  all_results <- list()
  idx <- 1L
  
  for (cell_type in sort(unique(stats::na.omit(metadata[[annotation_column]])))) {
    cell_meta <- metadata[metadata[[annotation_column]] == cell_type, , drop = FALSE]
    niche_sizes <- table(cell_meta$spatial_niche)
    valid_niches <- names(niche_sizes[niche_sizes >= min_cells_per_niche])
    cell_meta <- cell_meta[cell_meta$spatial_niche %in% valid_niches, , drop = FALSE]
    if (length(valid_niches) < 2) {
      next
    }
    
    cell_ids <- cell_meta$cell_ID
    cell_counts <- counts_cells_genes[cell_ids, , drop = FALSE]
    cell_meta$spatial_niche <- factor(cell_meta$spatial_niche, levels = valid_niches)
    
    cat("Running smiDE:", run_label, "|", cell_type, "|", paste(valid_niches, collapse = ", "), "\n")
    
    pre_obj <- tryCatch(
      smiDE::pre_de(
        adjacencies_only = FALSE,
        metadata = metadata,
        ref_celltype = cell_type,
        cell_type_metadata_colname = annotation_column,
        cellid_colname = "cell_ID",
        sdimx_colname = "sdimx",
        sdimy_colname = "sdimy",
        split_neighbors_by_colname = split_col,
        mm_radius = smide_radius,
        counts = counts_cells_genes,
        verbose = FALSE
      ),
      error = function(e) {
        message("Skipping smiDE pre_de for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(pre_obj)) {
      next
    }
    
    file_stub <- paste(
      gsub("[^A-Za-z0-9]+", "_", run_label),
      "sample",
      gsub("[^A-Za-z0-9]+", "_", cell_type),
      "smide",
      sep = "_"
    )
    
    if (isTRUE(smide_save_raw)) {
      saveRDS(pre_obj, file.path(tables_dir, paste0(file_stub, "_pre_de.rds")))
    }
    
    overlap_df <- tryCatch(
      coerce_overlap_metric_table(
        smiDE::overlap_ratio_metric(
          assay_matrix = counts_cells_genes,
          metadata = metadata,
          cluster_col = annotation_column,
          cellid_col = "cell_ID",
          sdimx_col = "sdimx",
          sdimy_col = "sdimy",
          split_neighbors_by_colname = split_col,
          verbose = FALSE
        )
      ),
      error = function(e) {
        message("Could not compute overlap ratio metric for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (!is.null(overlap_df)) {
      utils::write.csv(
        overlap_df,
        file.path(tables_dir, paste0(file_stub, "_overlap_ratio_metric.csv")),
        row.names = FALSE
      )
    }
    
    targets <- extract_overlap_targets(overlap_df, threshold = smide_overlap_threshold)
    if (!is.null(targets) && length(targets) == 0) {
      message("All targets were removed by the overlap threshold for ", run_label, " / ", cell_type)
      next
    }
    
    fit <- tryCatch(
      smiDE::smi_de(
        assay_matrix = cell_counts,
        metadata = cell_meta,
        formula = stats::as.formula("~ spatial_niche"),
        pre_de_obj = pre_obj,
        groupVar = "spatial_niche",
        groupVar_levels = valid_niches,
        nCores = smide_ncores,
        family = smide_family,
        targets = targets,
        cellid_colname = "cell_ID"
      ),
      error = function(e) {
        message("Skipping smiDE model for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) {
      next
    }
    
    if (isTRUE(smide_save_raw)) {
      saveRDS(fit, file.path(tables_dir, paste0(file_stub, "_smide_fit.rds")))
    }
    
    comparison_mode <- if (sample_contrast == "pairwise") "pairwise" else "one.vs.rest"
    res_obj <- tryCatch(
      smiDE::results(
        smide_results = fit,
        comparisons = comparison_mode,
        variable = "spatial_niche",
        targets = "all"
      ),
      error = function(e) {
        message("Could not extract smiDE results for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    
    if (isTRUE(smide_save_raw)) {
      saveRDS(res_obj, file.path(tables_dir, paste0(file_stub, "_smide_results_raw.rds")))
    }
    
    result_df <- coerce_smide_results_table(res_obj)
    if (is.null(result_df) || nrow(result_df) == 0) {
      summary_rows[[idx]] <- data.frame(
        analysis_scope = "sample",
        backend = "smiDE",
        run_label = run_label,
        sample_id = sample_ids[1] %||% run_label,
        annotation = cell_type,
        spatial_niche = paste(valid_niches, collapse = " / "),
        comparison = comparison_mode,
        replicate_column = split_col %||% "<none>",
        n_cells = nrow(cell_meta),
        n_genes_tested = if (is.null(targets)) ncol(cell_counts) else length(targets),
        n_fdr_005 = NA_integer_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
      next
    }
    
    result_df$analysis_scope <- "sample"
    result_df$backend <- "smiDE"
    result_df$run_label <- run_label
    result_df$sample_id <- sample_ids[1] %||% run_label
    result_df$annotation <- cell_type
    result_df$replicate_column <- split_col %||% "<none>"
    all_results[[idx]] <- result_df
    
    utils::write.csv(
      result_df,
      file.path(tables_dir, paste0(file_stub, "_spatial_de.csv")),
      row.names = FALSE
    )
    
    sig_col <- guess_significance_column(result_df)
    n_sig <- if (!is.null(sig_col)) sum(result_df[[sig_col]] < 0.05, na.rm = TRUE) else NA_integer_
    
    summary_rows[[idx]] <- data.frame(
      analysis_scope = "sample",
      backend = "smiDE",
      run_label = run_label,
      sample_id = sample_ids[1] %||% run_label,
      annotation = cell_type,
      spatial_niche = paste(valid_niches, collapse = " / "),
      comparison = comparison_mode,
      replicate_column = split_col %||% "<none>",
      n_cells = nrow(cell_meta),
      n_genes_tested = nrow(result_df),
      n_fdr_005 = n_sig,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
  
  list(
    metadata = metadata,
    summary = if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame(),
    results = if (length(all_results) > 0) do.call(rbind, all_results) else data.frame()
  )
}

detect_annotation_column <- function(metadata, preferred = NULL) {
  if (!is.null(preferred) && preferred %in% names(metadata)) {
    return(preferred)
  }
  
  candidates <- c(
    "celltype",
    grep("^celltype", names(metadata), value = TRUE),
    grep("annotation|annot", names(metadata), value = TRUE, ignore.case = TRUE),
    grep("leiden|cluster", names(metadata), value = TRUE, ignore.case = TRUE)
  )
  candidates <- unique(candidates[candidates %in% names(metadata)])
  if (length(candidates) == 0) {
    stop("No annotation column was found for spatial differential expression.")
  }
  candidates[1]
}

detect_sample_column <- function(metadata, preferred = "sample_id", required = TRUE) {
  candidates <- c(preferred, "sample_id", "list_ID", "slide_num")
  candidates <- unique(candidates[candidates %in% names(metadata)])
  if (length(candidates) == 0) {
    if (required) {
      stop("No sample identifier column was found in the input object.")
    }
    return(NULL)
  }
  candidates[1]
}

detect_treatment_column <- function(metadata, preferred = "treatment", required = TRUE) {
  candidates <- c(preferred, grep("treatment|condition|group", names(metadata), value = TRUE, ignore.case = TRUE))
  candidates <- unique(candidates[candidates %in% names(metadata)])
  if (length(candidates) == 0) {
    if (required) {
      stop("No treatment column was found in the input object.")
    }
    return(NULL)
  }
  candidates[1]
}

detect_patient_column <- function(metadata, preferred = "patient_id") {
  candidates <- c(preferred, grep("patient", names(metadata), value = TRUE, ignore.case = TRUE))
  candidates <- unique(candidates[candidates %in% names(metadata)])
  if (length(candidates) == 0) {
    return(NULL)
  }
  candidates[1]
}

detect_replicate_column <- function(metadata, preferred = "fov") {
  candidates <- c(
    preferred,
    "fov",
    "fov_ID",
    "fov_id",
    "fov_num",
    "FOV",
    grep("^fov", names(metadata), value = TRUE, ignore.case = TRUE)
  )
  candidates <- unique(candidates[candidates %in% names(metadata)])
  if (length(candidates) == 0) {
    return(NULL)
  }
  candidates[1]
}

infer_analysis_scope <- function(metadata,
                                 requested_scope = c("auto", "sample", "merged"),
                                 sample_column = NULL,
                                 treatment_column = NULL) {
  requested_scope <- match.arg(requested_scope)
  if (requested_scope != "auto") {
    return(requested_scope)
  }
  
  n_samples <- if (!is.null(sample_column) && sample_column %in% names(metadata)) {
    length(unique(stats::na.omit(metadata[[sample_column]])))
  } else {
    1L
  }
  n_treatments <- if (!is.null(treatment_column) && treatment_column %in% names(metadata)) {
    length(unique(stats::na.omit(metadata[[treatment_column]])))
  } else {
    0L
  }
  
  if (n_samples > 1 && n_treatments > 1) {
    return("merged")
  }
  "sample"
}

ensure_spatial_network <- function(gobj, spatial_network_name = "Delaunay_network") {
  network_exists <- tryCatch({
    getSpatialNetwork(gobject = gobj, name = spatial_network_name, output = "networkDT")
    TRUE
  }, error = function(e) {
    tryCatch({
      get_spatialNetwork(gobject = gobj, name = spatial_network_name, output = "networkDT")
      TRUE
    }, error = function(e2) FALSE)
  })
  
  if (network_exists) {
    return(gobj)
  }
  
  createSpatialNetwork(
    gobject = gobj,
    name = spatial_network_name,
    method = "Delaunay",
    minimum_k = 2
  )
}

get_network_dt <- function(gobj, spatial_network_name) {
  tryCatch({
    getSpatialNetwork(gobject = gobj, name = spatial_network_name, output = "networkDT")
  }, error = function(e) {
    get_spatialNetwork(gobject = gobj, name = spatial_network_name, output = "networkDT")
  })
}

get_spatial_locations_df <- function(gobj) {
  as.data.frame(getSpatialLocations(gobj, output = "data.table"), stringsAsFactors = FALSE)
}

normalise_rows <- function(mat) {
  rs <- rowSums(mat)
  rs[rs == 0] <- 1
  mat / rs
}

paired_design_possible <- function(meta_df, treatment_column, patient_column) {
  if (is.null(patient_column) || !patient_column %in% names(meta_df)) {
    return(FALSE)
  }
  patient_table <- table(meta_df[[patient_column]], meta_df[[treatment_column]])
  nrow(patient_table) >= 2 && all(rowSums(patient_table > 0) >= 2)
}

build_neighbourhood_matrix <- function(gobj,
                                       metadata,
                                       annotation_column,
                                       spatial_network_name) {
  net_dt <- as.data.frame(get_network_dt(gobj, spatial_network_name))
  if (!all(c("from", "to") %in% names(net_dt))) {
    stop("Spatial network table must contain 'from' and 'to' columns.")
  }
  
  cell_types <- setNames(metadata[[annotation_column]], metadata$cell_ID)
  edges <- rbind(
    data.frame(cell_ID = net_dt$from, neighbor_ID = net_dt$to, stringsAsFactors = FALSE),
    data.frame(cell_ID = net_dt$to, neighbor_ID = net_dt$from, stringsAsFactors = FALSE)
  )
  edges$neighbor_type <- unname(cell_types[edges$neighbor_ID])
  edges <- edges[!is.na(edges$neighbor_type), , drop = FALSE]
  
  niche_mat <- xtabs(~ cell_ID + neighbor_type, data = edges)
  niche_mat <- as.matrix(niche_mat)
  
  missing_cells <- setdiff(metadata$cell_ID, rownames(niche_mat))
  if (length(missing_cells) > 0) {
    filler <- matrix(
      0,
      nrow = length(missing_cells),
      ncol = ncol(niche_mat),
      dimnames = list(missing_cells, colnames(niche_mat))
    )
    niche_mat <- rbind(niche_mat, filler)
  }
  
  niche_mat[metadata$cell_ID, , drop = FALSE]
}

assign_spatial_niches <- function(niche_props, n_niches = 6) {
  if (nrow(niche_props) == 0) {
    return(character())
  }
  
  max_centers <- min(
    n_niches,
    nrow(niche_props),
    nrow(unique(round(niche_props, 6)))
  )
  
  if (max_centers <= 1) {
    return(rep("niche_1", nrow(niche_props)))
  }
  
  set.seed(42)
  niche_fit <- stats::kmeans(niche_props, centers = max_centers, nstart = 25)
  paste0("niche_", niche_fit$cluster)
}

aggregate_pseudobulk <- function(expr_mat, meta_df, group_columns) {
  keep_cols <- c("cell_ID", group_columns)
  meta_df <- meta_df[, keep_cols, drop = FALSE]
  meta_df <- meta_df[meta_df$cell_ID %in% colnames(expr_mat), , drop = FALSE]
  if (nrow(meta_df) == 0) {
    return(NULL)
  }
  
  meta_df$group_id <- do.call(
    paste,
    c(meta_df[group_columns], sep = "__")
  )
  group_factor <- factor(meta_df$group_id, levels = unique(meta_df$group_id))
  design_mat <- Matrix::sparse.model.matrix(~ 0 + group_factor)
  colnames(design_mat) <- levels(group_factor)
  
  counts <- expr_mat[, meta_df$cell_ID, drop = FALSE] %*% design_mat
  counts <- as.matrix(counts)
  
  pb_meta <- unique(meta_df[, c("group_id", group_columns), drop = FALSE])
  pb_meta <- pb_meta[match(colnames(counts), pb_meta$group_id), , drop = FALSE]
  pb_meta$n_cells <- as.integer(table(meta_df$group_id)[pb_meta$group_id])
  rownames(pb_meta) <- pb_meta$group_id
  
  list(counts = counts, meta = pb_meta)
}

run_pairwise_edgeR <- function(count_mat,
                               meta_df,
                               treatment_column,
                               patient_column = NULL,
                               comparison_label) {
  treatment_factor <- factor(meta_df[[treatment_column]])
  if (nlevels(treatment_factor) != 2) {
    stop("Exactly two treatment groups are required for a pairwise DE test.")
  }
  
  paired <- paired_design_possible(meta_df, treatment_column, patient_column)
  if (paired) {
    meta_df[[patient_column]] <- factor(meta_df[[patient_column]])
    design <- model.matrix(
      stats::as.formula(paste("~", patient_column, "+", treatment_column)),
      data = meta_df
    )
  } else {
    design <- model.matrix(
      stats::as.formula(paste("~", treatment_column)),
      data = meta_df
    )
  }
  
  dge <- edgeR::DGEList(counts = count_mat)
  keep <- edgeR::filterByExpr(dge, design = design)
  if (!any(keep)) {
    return(NULL)
  }
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  
  coef_name <- grep(paste0("^", treatment_column), colnames(design), value = TRUE)[1]
  if (is.na(coef_name) || !nzchar(coef_name)) {
    stop("Could not identify a treatment coefficient in the design matrix.")
  }
  
  test <- edgeR::glmQLFTest(fit, coef = coef_name)
  de_tbl <- edgeR::topTags(test, n = Inf)$table
  de_tbl$gene <- rownames(de_tbl)
  de_tbl$comparison <- comparison_label
  de_tbl$paired_design <- paired
  rownames(de_tbl) <- NULL
  de_tbl
}

run_group_edgeR <- function(count_mat,
                            group_vector,
                            target_group,
                            reference_group,
                            comparison_label) {
  group_factor <- factor(group_vector, levels = c(reference_group, target_group))
  if (nlevels(group_factor) != 2 || any(table(group_factor) == 0)) {
    stop("Exactly two non-empty groups are required for this DE test.")
  }
  
  design <- model.matrix(~ group_factor)
  dge <- edgeR::DGEList(counts = count_mat)
  keep <- edgeR::filterByExpr(dge, group = group_factor)
  if (!any(keep)) {
    return(NULL)
  }
  
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)
  fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
  
  coef_name <- grep("^group_factor", colnames(design), value = TRUE)[1]
  test <- edgeR::glmQLFTest(fit, coef = coef_name)
  de_tbl <- edgeR::topTags(test, n = Inf)$table
  de_tbl$gene <- rownames(de_tbl)
  de_tbl$comparison <- comparison_label
  de_tbl$target_group <- target_group
  de_tbl$reference_group <- reference_group
  rownames(de_tbl) <- NULL
  de_tbl
}

create_spatial_patches <- function(metadata,
                                   spatial_locations,
                                   n_spatial_patches = 8) {
  xy <- spatial_locations[match(metadata$cell_ID, spatial_locations$cell_ID), c("sdimx", "sdimy"), drop = FALSE]
  xy <- as.matrix(xy)
  if (nrow(xy) < 2 || anyNA(xy)) {
    stop("Could not create spatial patches because spatial coordinates are missing.")
  }
  
  n_centers <- min(n_spatial_patches, nrow(unique(round(xy, 3))))
  if (n_centers < 2) {
    metadata$spatial_patch <- "patch_1"
    return(list(metadata = metadata, replicate_column = "spatial_patch"))
  }
  
  set.seed(42)
  km <- stats::kmeans(xy, centers = n_centers, nstart = 25)
  metadata$spatial_patch <- paste0("patch_", km$cluster)
  list(metadata = metadata, replicate_column = "spatial_patch")
}

ensure_replicate_units <- function(metadata,
                                   spatial_locations,
                                   sample_replicate_column = "fov",
                                   n_spatial_patches = 8) {
  replicate_column <- detect_replicate_column(metadata, sample_replicate_column)
  if (!is.null(replicate_column)) {
    return(list(metadata = metadata, replicate_column = replicate_column, generated = FALSE))
  }
  
  patched <- create_spatial_patches(
    metadata = metadata,
    spatial_locations = spatial_locations,
    n_spatial_patches = n_spatial_patches
  )
  list(
    metadata = patched$metadata,
    replicate_column = patched$replicate_column,
    generated = TRUE
  )
}

save_niche_plots <- function(metadata,
                             spatial_locations,
                             output_dir,
                             run_label,
                             annotation_column) {
  merged_df <- merge(
    spatial_locations[, c("cell_ID", "sdimx", "sdimy")],
    metadata[, c("cell_ID", "spatial_niche", annotation_column)],
    by = "cell_ID",
    all.x = TRUE
  )
  
  p1 <- ggplot2::ggplot(
    merged_df,
    ggplot2::aes(x = sdimx, y = sdimy, colour = spatial_niche)
  ) +
    ggplot2::geom_point(size = 0.35, alpha = 0.8) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste(run_label, "- Spatial niche assignments"),
      x = "x",
      y = "y",
      colour = "Spatial niche"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(run_label, "_spatial_niches.png")),
    plot = p1,
    width = 13,
    height = 10,
    dpi = 300,
    bg = "white"
  )
  
  niche_comp <- as.data.frame(table(metadata$spatial_niche, metadata[[annotation_column]]), stringsAsFactors = FALSE)
  names(niche_comp) <- c("spatial_niche", "annotation", "n_cells")
  
  p2 <- ggplot2::ggplot(
    niche_comp,
    ggplot2::aes(x = spatial_niche, y = n_cells, fill = annotation)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::labs(
      title = paste(run_label, "- Cell-type composition per spatial niche"),
      x = "Spatial niche",
      y = "Proportion",
      fill = "Annotation"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(run_label, "_niche_composition.png")),
    plot = p2,
    width = 12,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

run_sample_scope_spatial_de <- function(expr_mat,
                                        metadata,
                                        run_label,
                                        annotation_column,
                                        tables_dir,
                                        sample_replicate_column = "fov",
                                        sample_contrast = c("one_vs_rest", "pairwise"),
                                        min_cells_per_niche = 30,
                                        min_cells_per_replicate = 5,
                                        min_replicates_per_group = 2,
                                        n_spatial_patches = 8,
                                        spatial_locations = NULL) {
  sample_contrast <- match.arg(sample_contrast)
  
  sample_ids <- unique(stats::na.omit(metadata$sample_id))
  if (length(sample_ids) > 1) {
    stop("Sample-scope spatial DE expects a single-sample object.")
  }
  
  replicate_info <- ensure_replicate_units(
    metadata = metadata,
    spatial_locations = spatial_locations,
    sample_replicate_column = sample_replicate_column,
    n_spatial_patches = n_spatial_patches
  )
  metadata <- replicate_info$metadata
  replicate_column <- replicate_info$replicate_column
  
  utils::write.csv(
    metadata,
    file.path(tables_dir, paste0(run_label, "_cell_metadata_with_niches.csv")),
    row.names = FALSE
  )
  
  summary_rows <- list()
  all_results <- list()
  idx <- 1L
  
  for (cell_type in sort(unique(stats::na.omit(metadata[[annotation_column]])))) {
    cell_meta <- metadata[metadata[[annotation_column]] == cell_type, , drop = FALSE]
    if (nrow(cell_meta) < min_cells_per_niche) {
      next
    }
    
    replicate_counts <- as.data.frame(
      table(cell_meta$spatial_niche, cell_meta[[replicate_column]]),
      stringsAsFactors = FALSE
    )
    names(replicate_counts) <- c("spatial_niche", replicate_column, "n_cells")
    keep_groups <- replicate_counts[replicate_counts$n_cells >= min_cells_per_replicate, , drop = FALSE]
    if (nrow(keep_groups) == 0) {
      next
    }
    
    keep_keys <- paste(keep_groups$spatial_niche, keep_groups[[replicate_column]], sep = "__")
    cell_meta$keep_key <- paste(cell_meta$spatial_niche, cell_meta[[replicate_column]], sep = "__")
    cell_meta <- cell_meta[cell_meta$keep_key %in% keep_keys, , drop = FALSE]
    if (nrow(cell_meta) < min_cells_per_niche) {
      next
    }
    
    niche_sizes <- table(cell_meta$spatial_niche)
    valid_niches <- names(niche_sizes[niche_sizes >= min_cells_per_niche])
    cell_meta <- cell_meta[cell_meta$spatial_niche %in% valid_niches, , drop = FALSE]
    if (length(valid_niches) < 2) {
      next
    }
    
    pb <- aggregate_pseudobulk(
      expr_mat = expr_mat,
      meta_df = cell_meta,
      group_columns = c("spatial_niche", replicate_column)
    )
    if (is.null(pb)) {
      next
    }
    
    pb_meta <- pb$meta
    pb_meta[[replicate_column]] <- as.character(pb_meta[[replicate_column]])
    pb_meta$spatial_niche <- as.character(pb_meta$spatial_niche)
    
    if (sample_contrast == "one_vs_rest") {
      for (niche in sort(unique(pb_meta$spatial_niche))) {
        group_vector <- ifelse(pb_meta$spatial_niche == niche, niche, "other_niches")
        group_sizes <- table(group_vector)
        if (length(group_sizes) != 2 || any(group_sizes < min_replicates_per_group)) {
          next
        }
        
        comparison_label <- paste0(niche, "_vs_rest")
        cat("Running sample spatial DE:", run_label, "|", cell_type, "|", comparison_label, "\n")
        
        de_tbl <- tryCatch(
          run_group_edgeR(
            count_mat = pb$counts,
            group_vector = group_vector,
            target_group = niche,
            reference_group = "other_niches",
            comparison_label = comparison_label
          ),
          error = function(e) {
            message("Skipping ", run_label, " / ", cell_type, " / ", comparison_label, ": ", conditionMessage(e))
            NULL
          }
        )
        
        if (is.null(de_tbl) || nrow(de_tbl) == 0) {
          next
        }
        
        de_tbl$analysis_scope <- "sample"
        de_tbl$run_label <- run_label
        de_tbl$sample_id <- sample_ids[1] %||% run_label
        de_tbl$annotation <- cell_type
        de_tbl$spatial_niche <- niche
        de_tbl$replicate_column <- replicate_column
        all_results[[idx]] <- de_tbl
        
        file_stub <- paste(
          gsub("[^A-Za-z0-9]+", "_", run_label),
          "sample",
          gsub("[^A-Za-z0-9]+", "_", cell_type),
          gsub("[^A-Za-z0-9]+", "_", comparison_label),
          sep = "_"
        )
        utils::write.csv(
          de_tbl,
          file.path(tables_dir, paste0(file_stub, "_spatial_de.csv")),
          row.names = FALSE
        )
        
        summary_rows[[idx]] <- data.frame(
          analysis_scope = "sample",
          run_label = run_label,
          sample_id = sample_ids[1] %||% run_label,
          annotation = cell_type,
          spatial_niche = niche,
          comparison = comparison_label,
          replicate_column = replicate_column,
          n_pseudobulk_profiles = length(group_vector),
          n_genes_tested = nrow(de_tbl),
          n_fdr_005 = sum(de_tbl$FDR < 0.05, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    } else {
      niche_levels <- names(which(table(pb_meta$spatial_niche) >= min_replicates_per_group))
      if (length(niche_levels) < 2) {
        next
      }
      
      for (pair in combn(sort(niche_levels), 2, simplify = FALSE)) {
        keep_idx <- pb_meta$spatial_niche %in% pair
        pair_meta <- pb_meta[keep_idx, , drop = FALSE]
        pair_counts <- pb$counts[, keep_idx, drop = FALSE]
        group_sizes <- table(pair_meta$spatial_niche)
        if (any(group_sizes < min_replicates_per_group)) {
          next
        }
        
        comparison_label <- paste(pair, collapse = "_vs_")
        cat("Running sample spatial DE:", run_label, "|", cell_type, "|", comparison_label, "\n")
        
        de_tbl <- tryCatch(
          run_group_edgeR(
            count_mat = pair_counts,
            group_vector = pair_meta$spatial_niche,
            target_group = pair[2],
            reference_group = pair[1],
            comparison_label = comparison_label
          ),
          error = function(e) {
            message("Skipping ", run_label, " / ", cell_type, " / ", comparison_label, ": ", conditionMessage(e))
            NULL
          }
        )
        
        if (is.null(de_tbl) || nrow(de_tbl) == 0) {
          next
        }
        
        de_tbl$analysis_scope <- "sample"
        de_tbl$run_label <- run_label
        de_tbl$sample_id <- sample_ids[1] %||% run_label
        de_tbl$annotation <- cell_type
        de_tbl$spatial_niche <- paste(pair, collapse = " / ")
        de_tbl$replicate_column <- replicate_column
        all_results[[idx]] <- de_tbl
        
        file_stub <- paste(
          gsub("[^A-Za-z0-9]+", "_", run_label),
          "sample",
          gsub("[^A-Za-z0-9]+", "_", cell_type),
          gsub("[^A-Za-z0-9]+", "_", comparison_label),
          sep = "_"
        )
        utils::write.csv(
          de_tbl,
          file.path(tables_dir, paste0(file_stub, "_spatial_de.csv")),
          row.names = FALSE
        )
        
        summary_rows[[idx]] <- data.frame(
          analysis_scope = "sample",
          run_label = run_label,
          sample_id = sample_ids[1] %||% run_label,
          annotation = cell_type,
          spatial_niche = paste(pair, collapse = " / "),
          comparison = comparison_label,
          replicate_column = replicate_column,
          n_pseudobulk_profiles = nrow(pair_meta),
          n_genes_tested = nrow(de_tbl),
          n_fdr_005 = sum(de_tbl$FDR < 0.05, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  
  list(
    metadata = metadata,
    summary = if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame(),
    results = if (length(all_results) > 0) do.call(rbind, all_results) else data.frame()
  )
}

run_merged_scope_spatial_de <- function(expr_mat,
                                        metadata,
                                        run_label,
                                        annotation_column,
                                        sample_column,
                                        treatment_column,
                                        patient_column,
                                        tables_dir,
                                        min_cells_per_niche = 30,
                                        min_cells_per_sample = 10,
                                        min_samples_per_group = 2) {
  utils::write.csv(
    metadata,
    file.path(tables_dir, paste0(run_label, "_cell_metadata_with_niches.csv")),
    row.names = FALSE
  )
  
  summary_rows <- list()
  all_results <- list()
  idx <- 1L
  
  for (cell_type in sort(unique(stats::na.omit(metadata[[annotation_column]])))) {
    cell_type_meta <- metadata[metadata[[annotation_column]] == cell_type, , drop = FALSE]
    for (niche in sort(unique(stats::na.omit(cell_type_meta$spatial_niche)))) {
      niche_meta <- cell_type_meta[cell_type_meta$spatial_niche == niche, , drop = FALSE]
      if (nrow(niche_meta) < min_cells_per_niche) {
        next
      }
      
      sample_counts <- as.data.frame(table(niche_meta[[sample_column]]), stringsAsFactors = FALSE)
      names(sample_counts) <- c(sample_column, "n_cells")
      keep_samples <- sample_counts[sample_counts$n_cells >= min_cells_per_sample, sample_column]
      niche_meta <- niche_meta[niche_meta[[sample_column]] %in% keep_samples, , drop = FALSE]
      if (nrow(niche_meta) < min_cells_per_niche) {
        next
      }
      
      pb <- aggregate_pseudobulk(expr_mat, niche_meta, group_columns = sample_column)
      if (is.null(pb)) {
        next
      }
      
      pb_meta <- unique(niche_meta[, c(sample_column, treatment_column, patient_column), drop = FALSE])
      pb_meta <- pb_meta[match(colnames(pb$counts), pb_meta[[sample_column]]), , drop = FALSE]
      
      treatment_levels <- unique(stats::na.omit(pb_meta[[treatment_column]]))
      if (length(treatment_levels) < 2) {
        next
      }
      
      for (pair in combn(sort(treatment_levels), 2, simplify = FALSE)) {
        keep_idx <- pb_meta[[treatment_column]] %in% pair
        pair_meta <- pb_meta[keep_idx, , drop = FALSE]
        pair_counts <- pb$counts[, keep_idx, drop = FALSE]
        
        group_sizes <- table(pair_meta[[treatment_column]])
        if (length(group_sizes) != 2 || any(group_sizes < min_samples_per_group)) {
          next
        }
        
        comparison_label <- paste(pair, collapse = "_vs_")
        cat("Running merged spatial DE:", run_label, "|", cell_type, "|", niche, "|", comparison_label, "\n")
        
        de_tbl <- tryCatch(
          run_pairwise_edgeR(
            count_mat = pair_counts,
            meta_df = pair_meta,
            treatment_column = treatment_column,
            patient_column = patient_column,
            comparison_label = comparison_label
          ),
          error = function(e) {
            message("Skipping ", run_label, " / ", cell_type, " / ", niche, " / ", comparison_label, ": ", conditionMessage(e))
            NULL
          }
        )
        
        if (is.null(de_tbl) || nrow(de_tbl) == 0) {
          next
        }
        
        de_tbl$analysis_scope <- "merged"
        de_tbl$run_label <- run_label
        de_tbl$annotation <- cell_type
        de_tbl$spatial_niche <- niche
        all_results[[idx]] <- de_tbl
        
        file_stub <- paste(
          gsub("[^A-Za-z0-9]+", "_", run_label),
          "merged",
          gsub("[^A-Za-z0-9]+", "_", cell_type),
          niche,
          gsub("[^A-Za-z0-9]+", "_", comparison_label),
          sep = "_"
        )
        utils::write.csv(
          de_tbl,
          file.path(tables_dir, paste0(file_stub, "_spatial_de.csv")),
          row.names = FALSE
        )
        
        summary_rows[[idx]] <- data.frame(
          analysis_scope = "merged",
          run_label = run_label,
          annotation = cell_type,
          spatial_niche = niche,
          comparison = comparison_label,
          n_samples = nrow(pair_meta),
          n_genes_tested = nrow(de_tbl),
          n_fdr_005 = sum(de_tbl$FDR < 0.05, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  
  list(
    metadata = metadata,
    summary = if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame(),
    results = if (length(all_results) > 0) do.call(rbind, all_results) else data.frame()
  )
}

run_spatial_differential_expression <- function(gobj,
                                                run_label = "spatial_de",
                                                output_dir,
                                                analysis_scope = c("auto", "sample", "merged"),
                                                backend = c("smiDE", "edgeR"),
                                                annotation_column = NULL,
                                                sample_column = "sample_id",
                                                treatment_column = "treatment",
                                                patient_column = "patient_id",
                                                spatial_network_name = "Delaunay_network",
                                                n_niches = 6,
                                                min_cells_per_niche = 30,
                                                min_cells_per_sample = 10,
                                                min_samples_per_group = 2,
                                                sample_replicate_column = "fov",
                                                sample_contrast = c("one_vs_rest", "pairwise"),
                                                min_cells_per_replicate = 5,
                                                min_replicates_per_group = 2,
                                                n_spatial_patches = 8,
                                                smide_family = "nbinom2",
                                                smide_radius = 0.05,
                                                smide_ncores = 1,
                                                smide_overlap_threshold = NULL,
                                                smide_save_raw = TRUE,
                                                save_object = FALSE) {
  analysis_scope <- match.arg(analysis_scope)
  backend <- match.arg(backend)
  sample_contrast <- match.arg(sample_contrast)
  
  cat("\n========================================\n")
  cat("STEP 12: Spatial Differential Expression\n")
  cat("Run:", run_label, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    manifest_path <- file.path(gobj, "manifest.json")
    if (file.exists(manifest_path)) {
      gobj <- load_giotto_checkpoint(gobj)
    } else {
      gobj <- loadGiotto(gobj)
    }
  }
  
  results_dir <- ensure_dir(file.path(output_dir, "12_Spatial_Differential_Expression"))
  tables_dir <- ensure_dir(file.path(results_dir, "tables"))
  plots_dir <- ensure_dir(file.path(results_dir, "plots"))
  
  metadata <- as.data.frame(pDataDT(gobj), stringsAsFactors = FALSE)
  if (!"cell_ID" %in% names(metadata)) {
    stop("Cell metadata must contain a cell_ID column.")
  }
  
  annotation_column <- detect_annotation_column(metadata, annotation_column)
  sample_col_guess <- detect_sample_column(metadata, sample_column, required = FALSE)
  treatment_col_guess <- detect_treatment_column(metadata, treatment_column, required = FALSE)
  analysis_scope <- infer_analysis_scope(
    metadata = metadata,
    requested_scope = analysis_scope,
    sample_column = sample_col_guess,
    treatment_column = treatment_col_guess
  )
  
  if (is.null(sample_col_guess)) {
    metadata$sample_id <- run_label
    sample_column <- "sample_id"
  } else {
    sample_column <- sample_col_guess
  }
  metadata$sample_id <- metadata[[sample_column]]
  
  if (analysis_scope == "merged") {
    treatment_column <- detect_treatment_column(metadata, treatment_column, required = TRUE)
    patient_column <- detect_patient_column(metadata, patient_column)
    if (length(unique(stats::na.omit(metadata[[sample_column]]))) < 2) {
      stop("Merged-scope spatial DE requires multiple samples in the input object.")
    }
  } else {
    patient_column <- NULL
    unique_samples <- unique(stats::na.omit(metadata[[sample_column]]))
    if (length(unique_samples) > 1) {
      stop("Sample-scope spatial DE requires exactly one sample in the input object.")
    }
  }
  
  cat("Analysis scope:", analysis_scope, "\n")
  cat("Backend:", backend, "\n")
  cat("Annotation column:", annotation_column, "\n")
  cat("Sample column:", sample_column, "\n")
  if (analysis_scope == "merged") {
    cat("Treatment column:", treatment_column, "\n")
    cat("Patient column:", patient_column %||% "<none>", "\n")
  } else {
    cat("Replicate column:", sample_replicate_column, "(auto-detected if available)\n")
    cat("Sample contrast:", sample_contrast, "\n")
  }
  cat("Spatial network:", spatial_network_name, "\n\n")
  
  gobj <- ensure_spatial_network(gobj, spatial_network_name = spatial_network_name)
  
  niche_counts <- build_neighbourhood_matrix(
    gobj = gobj,
    metadata = metadata,
    annotation_column = annotation_column,
    spatial_network_name = spatial_network_name
  )
  niche_props <- normalise_rows(niche_counts)
  metadata$spatial_niche <- assign_spatial_niches(niche_props, n_niches = n_niches)
  
  spatial_locations <- get_spatial_locations_df(gobj)
  save_niche_plots(
    metadata = metadata,
    spatial_locations = spatial_locations,
    output_dir = plots_dir,
    run_label = run_label,
    annotation_column = annotation_column
  )
  
  niche_summary <- as.data.frame(table(metadata[[annotation_column]], metadata$spatial_niche), stringsAsFactors = FALSE)
  names(niche_summary) <- c("annotation", "spatial_niche", "n_cells")
  utils::write.csv(
    niche_summary,
    file.path(tables_dir, paste0(run_label, "_niche_celltype_summary.csv")),
    row.names = FALSE
  )
  
  expr_mat <- getExpression(gobj, values = "raw", output = "matrix")
  
  scope_out <- if (analysis_scope == "sample") {
    if (backend == "smiDE") {
      run_smide_sample_backend(
        expr_mat = expr_mat,
        metadata = metadata,
        run_label = run_label,
        annotation_column = annotation_column,
        tables_dir = tables_dir,
        sample_replicate_column = sample_replicate_column,
        sample_contrast = sample_contrast,
        min_cells_per_niche = min_cells_per_niche,
        spatial_locations = spatial_locations,
        smide_family = smide_family,
        smide_radius = smide_radius,
        smide_ncores = smide_ncores,
        smide_overlap_threshold = smide_overlap_threshold,
        smide_save_raw = smide_save_raw
      )
    } else {
      run_sample_scope_spatial_de(
        expr_mat = expr_mat,
        metadata = metadata,
        run_label = run_label,
        annotation_column = annotation_column,
        tables_dir = tables_dir,
        sample_replicate_column = sample_replicate_column,
        sample_contrast = sample_contrast,
        min_cells_per_niche = min_cells_per_niche,
        min_cells_per_replicate = min_cells_per_replicate,
        min_replicates_per_group = min_replicates_per_group,
        n_spatial_patches = n_spatial_patches,
        spatial_locations = spatial_locations
      )
    }
  } else {
    if (backend != "edgeR") {
      stop("Merged-scope spatial DE currently supports only the edgeR backend.")
    }
    run_merged_scope_spatial_de(
      expr_mat = expr_mat,
      metadata = metadata,
      run_label = run_label,
      annotation_column = annotation_column,
      sample_column = sample_column,
      treatment_column = treatment_column,
      patient_column = patient_column,
      tables_dir = tables_dir,
      min_cells_per_niche = min_cells_per_niche,
      min_cells_per_sample = min_cells_per_sample,
      min_samples_per_group = min_samples_per_group
    )
  }
  
  metadata <- scope_out$metadata
  combined_results <- scope_out$results
  summary_df <- scope_out$summary
  
  if (nrow(summary_df) == 0 || nrow(combined_results) == 0) {
    warning("No spatial DE comparisons met the current thresholds.")
  } else {
    utils::write.csv(
      combined_results,
      file.path(results_dir, paste0(run_label, "_all_spatial_de_results.csv")),
      row.names = FALSE
    )
    utils::write.csv(
      summary_df,
      file.path(results_dir, paste0(run_label, "_spatial_de_summary.csv")),
      row.names = FALSE
    )
  }
  
  if (save_object) {
    save_giotto_checkpoint(
      gobj = gobj,
      checkpoint_dir = file.path(output_dir, "Giotto_Object_Spatial_DE"),
      metadata = list(stage = "spatial_differential_expression", analysis_scope = analysis_scope)
    )
  }
  
  attr(gobj, "spatial_de_scope") <- analysis_scope
  attr(gobj, "spatial_de_summary") <- summary_df
  attr(gobj, "spatial_de_results") <- combined_results
  invisible(gobj)
}

if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) {
    bootstrap_script <- file.path(current_script_dir(), "Helper_Scripts", "Script_Bootstrap.R")
    if (file.exists(bootstrap_script)) {
      source(bootstrap_script, local = .GlobalEnv)
      bootstrap_pipeline_environment(current_script_dir(), load_pipeline_utils = TRUE, verbose = FALSE)
    }
    
    run_spatial_differential_expression(
      gobj = args[2],
      run_label = args[1],
      output_dir = args[3],
      analysis_scope = if (length(args) >= 4) args[4] else "auto",
      backend = if (length(args) >= 5) args[5] else "smiDE"
    )
  } else {
    stop("Usage: Rscript 12_Spatial_Differential_Expression.R <run_label> <input_path> <output_dir> [analysis_scope] [backend]")
  }
}

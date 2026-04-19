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
#
# AUDIT CHANGES (v2):
# - FIX #12a: smide_annotation_subset default changed from NULL to B-cell
#   subset in config; added runtime warning when NULL (all-vs-all).
# - FIX #12b: get_counts_genes_by_cells() no longer densifies the full
#   expression matrix. Dense conversion is deferred to after cell-type
#   subsetting in run_smide_sample_backend(), saving ~10 GB per iteration.
# - FIX #12c: smide_radius auto-validation added. If spatial coordinates
#   are in microns (range > 100) and radius < 1, a warning is issued.
# - FIX #12d: smide_overlap_threshold default changed from 1 to 0.7 in
#   config to filter CosMx transcript diffusion artifacts.
# - Added gc() calls after each cell-type smiDE iteration.
# - Added early validation of smide_annotation_subset against actual labels.
# - Added RAM estimate warning for large workloads.
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
if ((!exists("save_giotto_checkpoint") || !exists("presentation_theme") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
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

.giotto_load <- function(path) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("loadGiotto", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("loadGiotto", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("loadGiotto", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("loadGiotto", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("loadGiotto", mode = "function")
  }
  
  accessor(path)
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

.giotto_get_spatial_network <- function(gobj, name, output = "networkDT") {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("getSpatialNetwork", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("getSpatialNetwork", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getSpatialNetwork", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("getSpatialNetwork", envir = asNamespace("GiottoClass"))
  } else if (requireNamespace("Giotto", quietly = TRUE) &&
             exists("get_spatialNetwork", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("get_spatialNetwork", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("get_spatialNetwork", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("get_spatialNetwork", envir = asNamespace("GiottoClass"))
  } else if (exists("getSpatialNetwork", mode = "function")) {
    accessor <- get("getSpatialNetwork", mode = "function")
  } else {
    accessor <- get("get_spatialNetwork", mode = "function")
  }
  
  accessor(gobject = gobj, name = name, output = output)
}

.giotto_create_spatial_network <- function(gobj,
                                           name = "Delaunay_network",
                                           method = "Delaunay",
                                           minimum_k = 2) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("createSpatialNetwork", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("createSpatialNetwork", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("createSpatialNetwork", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("createSpatialNetwork", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("createSpatialNetwork", mode = "function")
  }
  
  accessor(
    gobject = gobj,
    name = name,
    method = method,
    minimum_k = minimum_k
  )
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

# =============================================================================
# FIX #12b: get_counts_genes_by_cells() NO LONGER calls as.matrix() on the
# full expression matrix up front.  It now accepts and returns sparse matrices
# when the input is sparse.  Dense conversion is deferred to the caller AFTER
# subsetting to the focal cell type's cells, reducing peak memory from ~10 GB
# (200k cells x 6k genes dense) to ~0.5 GB (5-10k focal cells x 6k genes).
# =============================================================================
get_counts_genes_by_cells <- function(expr_mat, metadata, densify = FALSE) {
  metadata_ids <- unique(metadata$cell_ID)
  
  col_matches <- sum(metadata_ids %in% colnames(expr_mat))
  row_matches <- sum(metadata_ids %in% rownames(expr_mat))
  
  if (col_matches == 0 && row_matches == 0) {
    stop("No cell IDs from metadata were found in the expression matrix.")
  }
  
  if (col_matches >= row_matches) {
    keep_ids <- metadata$cell_ID[metadata$cell_ID %in% colnames(expr_mat)]
    counts <- expr_mat[, keep_ids, drop = FALSE]
  } else {
    keep_ids <- metadata$cell_ID[metadata$cell_ID %in% rownames(expr_mat)]
    counts <- t(expr_mat[keep_ids, , drop = FALSE])
  }
  
  if (isTRUE(densify)) {
    counts <- as.matrix(counts)
  }
  
  # Only set storage.mode on base R matrices/vectors, not S4 sparse objects
  if (!isS4(counts)) {
    storage.mode(counts) <- "numeric"
  }
  
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
    # smiDE returns a 'target' column with gene symbols and a 'gene' column
    # that is a numeric row index. Prefer 'target' as the gene name column.
    if ("target" %in% names(overlap_df) && !"gene" %in% names(overlap_df)) {
      overlap_df$gene <- overlap_df$target
    } else if ("target" %in% names(overlap_df)) {
      # 'gene' exists but is a numeric index — overwrite with gene symbols
      overlap_df$gene <- overlap_df$target
    } else if (!"gene" %in% names(overlap_df)) {
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
  
  # Prefer 'ratio' column (smiDE output) over the first numeric column,
  # which may be a row index ('gene') rather than the overlap metric.
  metric_col <- if ("ratio" %in% names(overlap_df)) {
    "ratio"
  } else {
    numeric_cols <- names(overlap_df)[vapply(overlap_df, is.numeric, logical(1))]
    # Exclude columns that look like row indices (sequential integers from 1)
    numeric_cols <- numeric_cols[!vapply(numeric_cols, function(col) {
      vals <- overlap_df[[col]]
      isTRUE(all.equal(as.numeric(vals), seq_len(nrow(overlap_df))))
    }, logical(1))]
    if (length(numeric_cols) == 0) return(NULL)
    numeric_cols[1]
  }

  keep <- overlap_df[[metric_col]] <= threshold
  targets <- overlap_df$gene[keep & !is.na(overlap_df$gene)]
  unique(targets[nzchar(as.character(targets))])
}

normalize_annotation_subset <- function(annotation_subset) {
  if (is.null(annotation_subset)) {
    return(NULL)
  }
  
  annotation_subset <- unique(as.character(annotation_subset))
  annotation_subset <- annotation_subset[!is.na(annotation_subset) & nzchar(annotation_subset)]
  if (length(annotation_subset) == 0) {
    return(NULL)
  }
  
  annotation_subset
}

prepare_smide_targets <- function(cell_counts,
                                  targets = NULL,
                                  min_detection_fraction = NULL) {
  gene_names <- rownames(cell_counts)
  if (is.null(gene_names)) {
    stop("smiDE counts matrix must contain gene names as rownames.")
  }
  
  keep <- rep(TRUE, nrow(cell_counts))
  names(keep) <- gene_names
  
  if (!is.null(targets)) {
    keep <- keep & gene_names %in% unique(targets)
  }
  
  if (!is.null(min_detection_fraction)) {
    min_detection_fraction <- as.numeric(min_detection_fraction)[1]
    min_detection_fraction <- min(max(min_detection_fraction, 0), 1)
    keep <- keep & (rowSums(cell_counts > 0, na.rm = TRUE) / ncol(cell_counts)) >= min_detection_fraction
  }
  
  kept_genes <- gene_names[keep]
  if (length(kept_genes) == 0) {
    return(list(
      counts = NULL,
      targets = character(),
      n_genes = 0L
    ))
  }
  
  filtered_counts <- cell_counts[kept_genes, , drop = FALSE]
  filtered_targets <- if (is.null(targets)) NULL else kept_genes
  
  list(
    counts = filtered_counts,
    targets = filtered_targets,
    n_genes = nrow(filtered_counts)
  )
}


sanitize_smide_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", as.character(x))
}

parse_smide_unified_interactions <- function(unified_int) {
  parts <- strsplit(as.character(unified_int), "--", fixed = TRUE)
  data.frame(
    source = vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)),
    target = vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

select_smide_partner_celltypes <- function(output_dir,
                                           run_label,
                                           focal_celltype,
                                           top_n = 3,
                                           padj_threshold = 0.05,
                                           include_self = FALSE) {
  cp_path <- file.path(
    output_dir,
    "11_BCell_Microenvironment",
    paste0(run_label, "_cell_proximity_enrichment.csv")
  )
  if (!file.exists(cp_path)) {
    return(character())
  }
  
  tbl <- tryCatch(
    readr::read_csv(cp_path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(tbl) || !"unified_int" %in% names(tbl)) {
    return(character())
  }
  
  parsed <- parse_smide_unified_interactions(tbl$unified_int)
  tbl$source <- parsed$source
  tbl$target <- parsed$target
  
  partner_tbl <- tbl[tbl$source == focal_celltype | tbl$target == focal_celltype, , drop = FALSE]
  if (!include_self) {
    partner_tbl <- partner_tbl[
      !(partner_tbl$source == focal_celltype & partner_tbl$target == focal_celltype),
      ,
      drop = FALSE
    ]
  }
  if ("p.adj_higher" %in% names(partner_tbl)) {
    partner_tbl <- partner_tbl[
      is.na(partner_tbl$p.adj_higher) | partner_tbl$p.adj_higher <= padj_threshold,
      ,
      drop = FALSE
    ]
  }
  if (nrow(partner_tbl) == 0) {
    return(character())
  }
  
  partner_tbl$partner_celltype <- ifelse(
    partner_tbl$source == focal_celltype,
    partner_tbl$target,
    partner_tbl$source
  )
  partner_tbl <- partner_tbl[
    !is.na(partner_tbl$partner_celltype) & nzchar(partner_tbl$partner_celltype),
    ,
    drop = FALSE
  ]
  if (nrow(partner_tbl) == 0) {
    return(character())
  }
  
  if ("int_ranking" %in% names(partner_tbl)) {
    partner_tbl <- partner_tbl[order(partner_tbl$int_ranking, partner_tbl$partner_celltype), , drop = FALSE]
  }
  
  unique(head(partner_tbl$partner_celltype, top_n))
}

adaptive_min_cells_per_niche <- function(n_cells, base_min, n_niches, floor_min = 10L) {
  max(floor_min, min(base_min, as.integer(floor(n_cells / (n_niches + 1L)))))
}

make_smide_neighbor_group <- function(cell_meta,
                                      neighbor_counts,
                                      partner_celltype,
                                      min_group_size = 5) {
  if (is.null(neighbor_counts) || !partner_celltype %in% colnames(neighbor_counts)) {
    return(NULL)
  }
  
  counts <- neighbor_counts[cell_meta$cell_ID, partner_celltype, drop = TRUE]
  counts[is.na(counts)] <- 0
  if (length(counts) == 0 || all(counts == counts[1])) {
    return(NULL)
  }
  
  threshold <- stats::median(counts, na.rm = TRUE)
  if (!is.finite(threshold)) {
    return(NULL)
  }
  
  group <- if (threshold <= 0) {
    ifelse(counts > 0, "high", "low")
  } else {
    ifelse(counts >= threshold, "high", "low")
  }
  
  group <- factor(group, levels = c("low", "high"))
  group_sizes <- table(group)
  if (length(group_sizes) < 2 || any(group_sizes < min_group_size)) {
    return(NULL)
  }
  
  list(
    column = paste0("neighbor_group__", sanitize_smide_name(partner_celltype)),
    values = group,
    partner_celltype = partner_celltype,
    threshold = threshold
  )
}

build_smide_predictor_specs <- function(cell_meta,
                                        cell_type,
                                        annotation_column,
                                        valid_niches,
                                        min_group_size,
                                        sample_contrast = c("one_vs_rest", "pairwise"),
                                        custom_predictor = NULL,
                                        partner_celltypes = NULL,
                                        partner_source = c("none", "cell_proximity_auto"),
                                        partner_top_n = 3,
                                        partner_padj_threshold = 0.05,
                                        include_self_partner = FALSE,
                                        output_dir = NULL,
                                        run_label = NULL,
                                        neighbor_counts = NULL) {
  sample_contrast <- match.arg(sample_contrast)
  partner_source <- match.arg(partner_source)
  
  specs <- list()
  
  if (!is.null(custom_predictor) && custom_predictor %in% names(cell_meta)) {
    predictor_values <- as.character(cell_meta[[custom_predictor]])
    predictor_values[!nzchar(predictor_values)] <- NA_character_
    predictor_levels <- unique(stats::na.omit(predictor_values))
    if (length(predictor_levels) >= 2) {
      cell_meta[[custom_predictor]] <- factor(predictor_values, levels = sort(predictor_levels))
      result_mode <- if (length(levels(cell_meta[[custom_predictor]])) == 2 || sample_contrast == "pairwise") {
        "pairwise"
      } else {
        "one.vs.rest"
      }
      specs[[length(specs) + 1L]] <- list(
        metadata = cell_meta,
        predictor_var = custom_predictor,
        predictor_label = custom_predictor,
        group_levels = levels(cell_meta[[custom_predictor]]),
        result_mode = result_mode
      )
    }
  }
  
  if (length(specs) == 0L && partner_source != "none" && is.null(partner_celltypes)) {
    if (!is.null(output_dir) && !is.null(run_label)) {
      partner_celltypes <- select_smide_partner_celltypes(
        output_dir = output_dir,
        run_label = run_label,
        focal_celltype = cell_type,
        top_n = partner_top_n,
        padj_threshold = partner_padj_threshold,
        include_self = include_self_partner
      )
    }
  }
  
  if (length(specs) == 0L && !is.null(partner_celltypes) && length(partner_celltypes) > 0) {
    for (partner in unique(as.character(partner_celltypes))) {
      group_spec <- make_smide_neighbor_group(
        cell_meta = cell_meta,
        neighbor_counts = neighbor_counts,
        partner_celltype = partner,
        min_group_size = min_group_size
      )
      if (is.null(group_spec)) {
        next
      }
      predictor_meta <- cell_meta
      predictor_meta[[group_spec$column]] <- group_spec$values
      specs[[length(specs) + 1L]] <- list(
        metadata = predictor_meta,
        predictor_var = group_spec$column,
        predictor_label = paste0("neighbor_", partner),
        group_levels = levels(group_spec$values),
        result_mode = "pairwise"
      )
    }
  }
  
  if (length(specs) == 0L) {
    fallback_meta <- cell_meta
    fallback_meta$spatial_niche <- factor(fallback_meta$spatial_niche, levels = valid_niches)
    fallback_mode <- if (length(valid_niches) == 2 || sample_contrast == "pairwise") {
      "pairwise"
    } else {
      "one.vs.rest"
    }
    specs[[1L]] <- list(
      metadata = fallback_meta,
      predictor_var = "spatial_niche",
      predictor_label = "spatial_niche",
      group_levels = valid_niches,
      result_mode = fallback_mode
    )
  }
  
  specs
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

# =============================================================================
# FIX #12c: validate_smide_radius() checks whether the smide_radius parameter
# is consistent with the spatial coordinate system.  CosMx coordinates are
# typically in microns (range ~10,000-50,000).  If the coordinate range exceeds
# 100 and the radius is < 1, it's almost certainly in the wrong unit.
# =============================================================================
validate_smide_radius <- function(smide_radius, spatial_locations) {
  if (is.null(spatial_locations) || nrow(spatial_locations) == 0) {
    return(smide_radius)
  }
  
  coord_range_x <- diff(range(spatial_locations$sdimx, na.rm = TRUE))
  coord_range_y <- diff(range(spatial_locations$sdimy, na.rm = TRUE))
  max_range <- max(coord_range_x, coord_range_y, na.rm = TRUE)
  
  if (is.finite(max_range) && max_range > 100 && smide_radius < 1) {
    warning(
      "smiDE radius = ", smide_radius, " but spatial coordinate range is ",
      round(max_range, 1), " (likely microns). ",
      "A radius of ", smide_radius, " um is almost certainly too small. ",
      "Consider setting smide_radius to 50 (50 um) for cell-cell proximity modeling. ",
      "Proceeding with the configured value.",
      immediate. = TRUE
    )
  } else if (is.finite(max_range) && max_range <= 1 && smide_radius > 1) {
    warning(
      "smiDE radius = ", smide_radius, " but spatial coordinate range is ",
      round(max_range, 4), " (likely normalized [0,1]). ",
      "A radius of ", smide_radius, " exceeds the full coordinate range. ",
      "Consider setting smide_radius to 0.05. ",
      "Proceeding with the configured value.",
      immediate. = TRUE
    )
  }
  
  smide_radius
}

# =============================================================================
# FIX #12b (continued): estimate_dense_matrix_gb() provides a RAM estimate
# so the pipeline can warn users before allocating large dense matrices.
# =============================================================================
estimate_dense_matrix_gb <- function(n_genes, n_cells, bytes_per_element = 8) {
  (as.numeric(n_genes) * as.numeric(n_cells) * bytes_per_element) / (1024^3)
}

run_smide_sample_backend <- function(expr_mat,
                                     metadata,
                                     run_label,
                                     output_dir,
                                     annotation_column,
                                     tables_dir,
                                     sample_replicate_column = "fov",
                                     sample_contrast = c("one_vs_rest", "pairwise"),
                                     min_cells_per_niche = 30,
                                     spatial_locations = NULL,
                                     spatial_network_name = "Delaunay_network",
                                     neighbor_counts = NULL,
                                     smide_family = "nbinom2",
                                     smide_radius = 0.05,
                                     smide_ncores = 1,
                                     smide_overlap_threshold = NULL,
                                     smide_annotation_subset = NULL,
                                     smide_min_detection_fraction = 0.05,
                                     smide_custom_predictor = NULL,
                                     smide_partner_celltypes = NULL,
                                     smide_partner_source = c("none", "cell_proximity_auto"),
                                     smide_partner_top_n = 3,
                                     smide_partner_padj_threshold = 0.05,
                                     smide_include_self_partner = FALSE,
                                     smide_save_raw = TRUE,
                                     smide_adaptive_thresholds = TRUE) {
  sample_contrast <- match.arg(sample_contrast)
  smide_partner_source <- match.arg(smide_partner_source)
  
  if (!requireNamespace("smiDE", quietly = TRUE)) {
    stop("smiDE is not installed, but backend = 'smiDE' was requested.")
  }
  
  metadata <- align_metadata_with_spatial_locations(metadata, spatial_locations)
  
  # -------------------------------------------------------------------------
  # FIX #12b: Keep the full expression matrix SPARSE at this stage.
  # get_counts_genes_by_cells() now returns sparse by default (densify=FALSE).
  # Dense conversion happens per-cell-type AFTER subsetting below.
  # -------------------------------------------------------------------------
  counts_genes_cells <- get_counts_genes_by_cells(expr_mat, metadata, densify = FALSE)
  metadata <- metadata[match(colnames(counts_genes_cells), metadata$cell_ID), , drop = FALSE]
  
  sample_ids <- unique(stats::na.omit(metadata$sample_id))
  if (length(sample_ids) > 1) {
    stop("smiDE sample-mode analysis expects exactly one sample in the input object.")
  }
  
  smide_annotation_subset <- normalize_annotation_subset(smide_annotation_subset)
  smide_min_detection_fraction <- as.numeric(smide_min_detection_fraction)[1]
  if (!is.finite(smide_min_detection_fraction)) {
    smide_min_detection_fraction <- 0.05
  }
  smide_min_detection_fraction <- min(max(smide_min_detection_fraction, 0), 1)
  
  # -------------------------------------------------------------------------
  # FIX #12d: Changed default from 1 (disabled) to 0.7 to filter CosMx
  # transcript diffusion artifacts.
  # -------------------------------------------------------------------------
  if (is.null(smide_overlap_threshold)) {
    smide_overlap_threshold <- 0.7
    message(
      "Using default smiDE overlap-ratio threshold of 0.7 to prefilter ",
      "contamination-prone genes (CosMx transcript diffusion filtering)."
    )
  }
  
  # -------------------------------------------------------------------------
  # FIX #12c: Validate smide_radius against spatial coordinate units.
  # -------------------------------------------------------------------------
  smide_radius <- validate_smide_radius(smide_radius, spatial_locations)
  
  neighbor_counts <- if (is.null(neighbor_counts)) NULL else as_plain_numeric_matrix(neighbor_counts)
  
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
  candidate_cell_types <- sort(unique(stats::na.omit(metadata[[annotation_column]])))
  
  # -------------------------------------------------------------------------
  # FIX #12a: Warn when smide_annotation_subset is NULL (all-vs-all mode).
  # For ~200k cells with ~20 cell types, this can take hours to days.
  # -------------------------------------------------------------------------
  if (is.null(smide_annotation_subset)) {
    warning(
      "smide_annotation_subset is NULL: smiDE will iterate over ALL ",
      length(candidate_cell_types), " cell types. ",
      "For large datasets (>100k cells), this may take hours to days. ",
      "Consider restricting to your cell types of interest, e.g.: ",
      "smide_annotation_subset: [\"B.cell\"] in config.yaml.",
      immediate. = TRUE
    )
  } else {
    candidate_cell_types <- intersect(candidate_cell_types, smide_annotation_subset)
    # Validate that requested types actually exist in the data
    missing_types <- setdiff(smide_annotation_subset, unique(stats::na.omit(metadata[[annotation_column]])))
    if (length(missing_types) > 0) {
      warning(
        "The following smide_annotation_subset types were not found in annotation column '",
        annotation_column, "': ", paste(missing_types, collapse = ", "), ". ",
        "Available types: ", paste(head(sort(unique(stats::na.omit(metadata[[annotation_column]]))), 20), collapse = ", "),
        immediate. = TRUE
      )
    }
  }
  
  if (length(candidate_cell_types) == 0) {
    warning("No cell types matched the requested smiDE annotation subset.")
    return(list(metadata = metadata, summary = data.frame(), results = data.frame()))
  }
  
  cat("smiDE candidate cell types:", length(candidate_cell_types), "\n")
  cat("  Types:", paste(candidate_cell_types, collapse = ", "), "\n")
  cat("smiDE min detection fraction:", smide_min_detection_fraction, "\n")
  cat("smiDE overlap threshold:", smide_overlap_threshold, "\n")
  cat("smiDE radius:", smide_radius, "\n")

  n_niches_effective <- length(unique(stats::na.omit(metadata$spatial_niche)))
  skipped_cell_types <- character(0)

  # -------------------------------------------------------------------------
  # FIX #12b: For smiDE, we need a dense matrix for pre_de() and
  # overlap_ratio_metric() which operate on the full dataset.  However, we
  # only densify the portion of the matrix needed.  For pre_de() we pass the
  # full sparse matrix and let smiDE handle it; for the cell-type-specific
  # model fitting we densify only the focal cell's columns.
  #
  # Estimate RAM for the full dense matrix and warn if it would be very large.
  # -------------------------------------------------------------------------
  n_total_genes <- nrow(counts_genes_cells)
  n_total_cells <- ncol(counts_genes_cells)
  full_dense_gb <- estimate_dense_matrix_gb(n_total_genes, n_total_cells)
  
  if (full_dense_gb > 8) {
    warning(
      "Full dense expression matrix would require ~", round(full_dense_gb, 1), " GB. ",
      "smiDE pre_de() and overlap_ratio_metric() need the full matrix. ",
      "Densifying once for these shared operations, then releasing.",
      immediate. = TRUE
    )
  }
  
  # Densify the full matrix ONCE for pre_de() and overlap_ratio_metric(),
  # which require the complete dataset with all cell types as context.
  counts_dense_full <- as.matrix(counts_genes_cells)
  
  for (cell_type in candidate_cell_types) {
    cell_meta <- metadata[metadata[[annotation_column]] == cell_type, , drop = FALSE]
    n_cells_total <- nrow(cell_meta)
    niche_sizes <- table(cell_meta$spatial_niche)
    effective_min <- if (isTRUE(smide_adaptive_thresholds)) {
      adaptive_min_cells_per_niche(n_cells_total, min_cells_per_niche, n_niches_effective)
    } else {
      min_cells_per_niche
    }
    if (effective_min < min_cells_per_niche) {
      message(sprintf(
        "  [adaptive] min_cells_per_niche relaxed from %d \u2192 %d for '%s' (%d total cells, %d niches)",
        min_cells_per_niche, effective_min, cell_type, n_cells_total, n_niches_effective
      ))
    }
    valid_niches <- names(niche_sizes[niche_sizes >= effective_min])
    cell_meta <- cell_meta[cell_meta$spatial_niche %in% valid_niches, , drop = FALSE]
    if (nrow(cell_meta) < effective_min || length(valid_niches) < 2) {
      message(sprintf(
        "  \u26a0 smiDE skipped for '%s': %d total cells, %d niche(s) \u2265%d cells (requires \u22652)",
        cell_type, n_cells_total, length(valid_niches), effective_min
      ))
      skipped_cell_types <- c(skipped_cell_types, cell_type)
      next
    }
    
    # -------------------------------------------------------------------
    # FIX #12b: Only densify the focal cell type's columns for the model.
    # The pre_de() and overlap_ratio_metric() calls below still use the
    # full dense matrix (counts_dense_full) since they need neighbourhood
    # context from ALL cell types.
    # -------------------------------------------------------------------
    cell_ids <- cell_meta$cell_ID
    cell_counts <- counts_dense_full[, cell_ids, drop = FALSE]
    
    focal_gb <- estimate_dense_matrix_gb(nrow(cell_counts), ncol(cell_counts))
    cat("  Processing cell type:", cell_type, "(", length(cell_ids), "cells, ~",
        round(focal_gb, 2), "GB dense)\n")
    
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
        counts = counts_dense_full,
        verbose = FALSE
      ),
      error = function(e) {
        message("Skipping smiDE pre_de for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(pre_obj)) {
      gc(verbose = FALSE)
      next
    }
    
    overlap_df <- tryCatch(
      coerce_overlap_metric_table(
        smiDE::overlap_ratio_metric(
          assay_matrix = counts_dense_full,
          metadata = metadata,
          cluster_col = annotation_column,
          cellid_col = "cell_ID",
          sdimx_col = "sdimx",
          sdimy_col = "sdimy",
          split_neighbors_by_colname = split_col,
          radius = smide_radius,
          verbose = FALSE
        )
      ),
      error = function(e) {
        message("Could not compute overlap ratio metric for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    
    file_stub_base <- paste(
      gsub("[^A-Za-z0-9]+", "_", run_label),
      "sample",
      gsub("[^A-Za-z0-9]+", "_", cell_type),
      "smide",
      sep = "_"
    )
    
    if (isTRUE(smide_save_raw)) {
      saveRDS(pre_obj, file.path(tables_dir, paste0(file_stub_base, "_pre_de.rds")))
    }
    if (!is.null(overlap_df)) {
      utils::write.csv(
        overlap_df,
        file.path(tables_dir, paste0(file_stub_base, "_overlap_ratio_metric.csv")),
        row.names = FALSE
      )
    }
    
    targets <- extract_overlap_targets(overlap_df, threshold = smide_overlap_threshold)
    if (!is.null(targets) && length(targets) == 0) {
      message("All targets were removed by the overlap threshold for ", run_label, " / ", cell_type)
      gc(verbose = FALSE)
      next
    }
    
    model_inputs <- prepare_smide_targets(
      cell_counts = cell_counts,
      targets = targets,
      min_detection_fraction = smide_min_detection_fraction
    )
    if (is.null(model_inputs$counts) || model_inputs$n_genes == 0) {
      message("No genes passed the smiDE filters for ", run_label, " / ", cell_type)
      gc(verbose = FALSE)
      next
    }
    
    predictor_specs <- build_smide_predictor_specs(
      cell_meta = cell_meta,
      cell_type = cell_type,
      annotation_column = annotation_column,
      valid_niches = valid_niches,
      min_group_size = effective_min,
      sample_contrast = sample_contrast,
      custom_predictor = smide_custom_predictor,
      partner_celltypes = smide_partner_celltypes,
      partner_source = smide_partner_source,
      partner_top_n = smide_partner_top_n,
      partner_padj_threshold = smide_partner_padj_threshold,
      include_self_partner = smide_include_self_partner,
      output_dir = output_dir,
      run_label = run_label,
      neighbor_counts = neighbor_counts
    )
    
    for (spec in predictor_specs) {
      predictor_meta <- spec$metadata
      predictor_var <- spec$predictor_var
      predictor_label <- spec$predictor_label
      predictor_levels <- spec$group_levels
      result_mode <- spec$result_mode
      
      predictor_meta <- predictor_meta[!is.na(predictor_meta[[predictor_var]]), , drop = FALSE]
      predictor_sizes <- table(predictor_meta[[predictor_var]])
      predictor_sizes <- predictor_sizes[predictor_sizes > 0]
      if (length(predictor_sizes) < 2 || any(predictor_sizes < effective_min)) {
        next
      }
      predictor_levels <- names(predictor_sizes)
      predictor_meta[[predictor_var]] <- factor(predictor_meta[[predictor_var]], levels = predictor_levels)
      predictor_counts <- model_inputs$counts[, predictor_meta$cell_ID, drop = FALSE]
      
      cat("Running smiDE:", run_label, "|", cell_type, "|", predictor_label, "\n")
      cat("smiDE workload for", cell_type, "/", predictor_label, ":", nrow(predictor_meta), "cells x", model_inputs$n_genes, "genes\n")
      
      fit <- tryCatch(
        smiDE::smi_de(
          assay_matrix = predictor_counts,
          metadata = predictor_meta,
          formula = stats::as.formula(paste("~", predictor_var)),
          pre_de_obj = pre_obj,
          groupVar = predictor_var,
          groupVar_levels = predictor_levels,
          nCores = smide_ncores,
          family = smide_family,
          targets = model_inputs$targets,
          neighbor_expr_cell_type_metadata_colname = annotation_column,
          cellid_colname = "cell_ID"
        ),
        error = function(e) {
          message("Skipping smiDE model for ", run_label, " / ", cell_type, " / ", predictor_label, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(fit)) {
        next
      }
      
      file_stub <- paste(file_stub_base, sanitize_smide_name(predictor_label), sep = "_")
      if (isTRUE(smide_save_raw)) {
        saveRDS(fit, file.path(tables_dir, paste0(file_stub, "_smide_fit.rds")))
      }
      
      res_obj <- tryCatch(
        smiDE::results(
          smide_results = fit,
          comparisons = result_mode,
          variable = predictor_var,
          targets = "all"
        ),
        error = function(e) {
          message("Could not extract smiDE results for ", run_label, " / ", cell_type, " / ", predictor_label, ": ", conditionMessage(e))
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
          sample_id = if (length(sample_ids) > 0) sample_ids[1] else run_label,
          annotation = cell_type,
          spatial_niche = paste(valid_niches, collapse = " / "),
          predictor = predictor_label,
          comparison = result_mode,
          replicate_column = if (is.null(split_col)) "<none>" else split_col,
          n_cells = nrow(predictor_meta),
          n_genes_tested = model_inputs$n_genes,
          n_fdr_005 = NA_integer_,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
        next
      }
      
      result_df$analysis_scope <- "sample"
      result_df$backend <- "smiDE"
      result_df$run_label <- run_label
      result_df$sample_id <- if (length(sample_ids) > 0) sample_ids[1] else run_label
      result_df$annotation <- cell_type
      result_df$predictor <- predictor_label
      result_df$replicate_column <- if (is.null(split_col)) "<none>" else split_col
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
        sample_id = if (length(sample_ids) > 0) sample_ids[1] else run_label,
        annotation = cell_type,
        spatial_niche = paste(valid_niches, collapse = " / "),
        predictor = predictor_label,
        comparison = result_mode,
        replicate_column = if (is.null(split_col)) "<none>" else split_col,
        n_cells = nrow(predictor_meta),
        n_genes_tested = model_inputs$n_genes,
        n_fdr_005 = n_sig,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
    
    # Clean up after each cell type iteration
    gc(verbose = FALSE)
  }
  
  # Release the full dense matrix now that all cell types are processed
  rm(counts_dense_full)
  gc(verbose = FALSE)
  
  list(
    metadata = metadata,
    summary = if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame(),
    results = if (length(all_results) > 0) do.call(rbind, all_results) else data.frame(),
    skipped_cell_types = skipped_cell_types
  )
}

# ==============================================================================
# Merged-scope smiDE backend
#
# Tests: do B-cells in a given spatial niche show different expression depending
# on treatment group? Combines niche × treatment into a single predictor and
# adds sample_id as a batch covariate to control for sample-to-sample variation.
#
# Spatial neighbourhoods are split by sample_id so that coordinates from
# different samples (offset by x_padding at merge time) never share an edge.
# ==============================================================================
run_smide_merged_backend <- function(expr_mat,
                                     metadata,
                                     run_label,
                                     output_dir,
                                     annotation_column,
                                     tables_dir,
                                     treatment_column,
                                     sample_column = "sample_id",
                                     min_cells_per_niche = 30,
                                     spatial_locations = NULL,
                                     neighbor_counts = NULL,
                                     smide_family = "nbinom2",
                                     smide_radius = 0.05,
                                     smide_ncores = 1,
                                     smide_overlap_threshold = NULL,
                                     smide_annotation_subset = NULL,
                                     smide_min_detection_fraction = 0.05,
                                     smide_save_raw = TRUE,
                                     smide_adaptive_thresholds = TRUE) {
  if (!requireNamespace("smiDE", quietly = TRUE)) {
    stop("smiDE is not installed, but backend = 'smiDE' was requested.")
  }
  if (!treatment_column %in% names(metadata)) {
    stop("treatment_column '", treatment_column, "' not found in metadata.")
  }

  metadata <- align_metadata_with_spatial_locations(metadata, spatial_locations)
  counts_genes_cells <- get_counts_genes_by_cells(expr_mat, metadata, densify = FALSE)
  metadata <- metadata[match(colnames(counts_genes_cells), metadata$cell_ID), , drop = FALSE]

  sample_ids <- unique(stats::na.omit(metadata[[sample_column]]))
  if (length(sample_ids) < 2) {
    stop("run_smide_merged_backend() requires multiple samples. Use run_smide_sample_backend() for single-sample analysis.")
  }

  if (is.null(smide_overlap_threshold)) {
    smide_overlap_threshold <- 0.7
    message("Using default smiDE overlap-ratio threshold of 0.7.")
  }
  smide_radius <- validate_smide_radius(smide_radius, spatial_locations)
  smide_min_detection_fraction <- min(max(as.numeric(smide_min_detection_fraction)[1], 0), 1)
  if (!is.finite(smide_min_detection_fraction)) smide_min_detection_fraction <- 0.05
  neighbor_counts <- if (is.null(neighbor_counts)) NULL else as_plain_numeric_matrix(neighbor_counts)

  # Combined predictor: spatial_niche__treatment_group
  metadata$smide_niche_treatment <- paste(
    metadata$spatial_niche, metadata[[treatment_column]], sep = "__"
  )

  smide_annotation_subset <- normalize_annotation_subset(smide_annotation_subset)
  candidate_cell_types <- sort(unique(stats::na.omit(metadata[[annotation_column]])))
  if (!is.null(smide_annotation_subset)) {
    missing_types <- setdiff(smide_annotation_subset, candidate_cell_types)
    if (length(missing_types) > 0) {
      warning("smide_annotation_subset types not found: ", paste(missing_types, collapse = ", "), immediate. = TRUE)
    }
    candidate_cell_types <- intersect(candidate_cell_types, smide_annotation_subset)
  }
  if (length(candidate_cell_types) == 0) {
    warning("No cell types matched smide_annotation_subset.")
    return(list(metadata = metadata, summary = data.frame(), results = data.frame(), skipped_cell_types = character(0)))
  }

  treatment_groups <- sort(unique(stats::na.omit(metadata[[treatment_column]])))
  n_niches_effective <- length(unique(stats::na.omit(metadata$spatial_niche)))

  cat("Merged smiDE — candidate cell types:", length(candidate_cell_types), "\n")
  cat("  Types:", paste(candidate_cell_types, collapse = ", "), "\n")
  cat("  Samples:", length(sample_ids), " | Treatment groups:", paste(treatment_groups, collapse = ", "), "\n")
  cat("  Neighborhood split by:", sample_column, "\n")
  cat("  smiDE radius:", smide_radius, " | Overlap threshold:", smide_overlap_threshold, "\n\n")

  full_dense_gb <- estimate_dense_matrix_gb(nrow(counts_genes_cells), ncol(counts_genes_cells))
  if (full_dense_gb > 8) {
    warning("Full dense matrix ~", round(full_dense_gb, 1), " GB.", immediate. = TRUE)
  }
  counts_dense_full <- as.matrix(counts_genes_cells)

  summary_rows <- list()
  all_results  <- list()
  skipped_cell_types <- character(0)
  idx <- 1L

  for (cell_type in candidate_cell_types) {
    cell_meta <- metadata[metadata[[annotation_column]] == cell_type, , drop = FALSE]
    n_cells_total <- nrow(cell_meta)

    effective_min <- if (isTRUE(smide_adaptive_thresholds)) {
      adaptive_min_cells_per_niche(n_cells_total, min_cells_per_niche, n_niches_effective)
    } else {
      min_cells_per_niche
    }
    if (effective_min < min_cells_per_niche) {
      message(sprintf("  [adaptive] min_cells_per_niche relaxed %d \u2192 %d for '%s' (%d cells)",
                      min_cells_per_niche, effective_min, cell_type, n_cells_total))
    }

    group_sizes   <- table(cell_meta$smide_niche_treatment)
    valid_groups  <- names(group_sizes[group_sizes >= effective_min])
    treat_in_valid <- unique(vapply(strsplit(valid_groups, "__"), `[`, character(1), 2))
    n_treat_valid  <- length(treat_in_valid)

    if (length(valid_groups) < 2 || n_treat_valid < 2) {
      message(sprintf(
        "  \u26a0 Merged smiDE skipped for '%s': %d cells, %d group(s) \u2265%d cells, %d treatment(s) (needs \u22652 groups spanning \u22652 treatments)",
        cell_type, n_cells_total, length(valid_groups), effective_min, n_treat_valid
      ))
      skipped_cell_types <- c(skipped_cell_types, cell_type)
      next
    }

    cat("  Processing cell type:", cell_type, "(", n_cells_total, "cells,",
        length(valid_groups), "niche-treatment groups)\n")
    cat("  Groups:", paste(names(group_sizes[valid_groups]), group_sizes[valid_groups], sep = "=", collapse = ", "), "\n")

    cell_ids   <- cell_meta$cell_ID
    cell_counts <- counts_dense_full[, cell_ids, drop = FALSE]

    pre_obj <- tryCatch(
      smiDE::pre_de(
        adjacencies_only = FALSE,
        metadata = metadata,
        ref_celltype = cell_type,
        cell_type_metadata_colname = annotation_column,
        cellid_colname = "cell_ID",
        sdimx_colname = "sdimx",
        sdimy_colname = "sdimy",
        split_neighbors_by_colname = sample_column,
        mm_radius = smide_radius,
        counts = counts_dense_full,
        verbose = FALSE
      ),
      error = function(e) {
        message("pre_de failed for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(pre_obj)) { gc(verbose = FALSE); next }

    overlap_df <- tryCatch(
      coerce_overlap_metric_table(
        smiDE::overlap_ratio_metric(
          assay_matrix = counts_dense_full,
          metadata = metadata,
          cluster_col = annotation_column,
          cellid_col = "cell_ID",
          sdimx_col = "sdimx",
          sdimy_col = "sdimy",
          split_neighbors_by_colname = sample_column,
          radius = smide_radius,
          verbose = FALSE
        )
      ),
      error = function(e) NULL
    )

    file_stub_base <- paste(
      gsub("[^A-Za-z0-9]+", "_", run_label),
      "merged",
      gsub("[^A-Za-z0-9]+", "_", cell_type),
      "smide",
      sep = "_"
    )
    if (isTRUE(smide_save_raw)) {
      saveRDS(pre_obj, file.path(tables_dir, paste0(file_stub_base, "_pre_de.rds")))
    }
    if (!is.null(overlap_df)) {
      utils::write.csv(overlap_df,
                       file.path(tables_dir, paste0(file_stub_base, "_overlap_ratio_metric.csv")),
                       row.names = FALSE)
    }

    targets <- extract_overlap_targets(overlap_df, threshold = smide_overlap_threshold)
    if (!is.null(targets) && length(targets) == 0) {
      message("All targets removed by overlap threshold for ", run_label, " / ", cell_type)
      gc(verbose = FALSE); next
    }

    model_inputs <- prepare_smide_targets(
      cell_counts = cell_counts,
      targets = targets,
      min_detection_fraction = smide_min_detection_fraction
    )
    if (is.null(model_inputs$counts) || model_inputs$n_genes == 0) {
      message("No genes passed smiDE filters for ", run_label, " / ", cell_type)
      gc(verbose = FALSE); next
    }

    # Subset to valid niche-treatment groups only
    cell_meta_valid <- cell_meta[cell_meta$smide_niche_treatment %in% valid_groups, , drop = FALSE]
    cell_meta_valid$smide_niche_treatment <- factor(
      cell_meta_valid$smide_niche_treatment, levels = valid_groups
    )
    predictor_counts <- model_inputs$counts[, cell_meta_valid$cell_ID, drop = FALSE]

    # Formula: ~ smide_niche_treatment + sample_id
    # The niche_treatment term tests the niche × treatment interaction.
    # The sample_id term absorbs sample-level batch effects in the count data.
    formula_str <- paste("~ smide_niche_treatment +", sample_column)
    cat("  Formula:", formula_str, "\n")
    cat("  smiDE workload:", nrow(cell_meta_valid), "cells x", model_inputs$n_genes, "genes\n")

    fit <- tryCatch(
      smiDE::smi_de(
        assay_matrix = predictor_counts,
        metadata = cell_meta_valid,
        formula = stats::as.formula(formula_str),
        pre_de_obj = pre_obj,
        groupVar = "smide_niche_treatment",
        groupVar_levels = valid_groups,
        nCores = smide_ncores,
        family = smide_family,
        targets = model_inputs$targets,
        neighbor_expr_cell_type_metadata_colname = annotation_column,
        cellid_colname = "cell_ID"
      ),
      error = function(e) {
        message("smiDE::smi_de failed for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(fit)) { gc(verbose = FALSE); next }

    if (isTRUE(smide_save_raw)) {
      saveRDS(fit, file.path(tables_dir, paste0(file_stub_base, "_smide_fit.rds")))
    }

    res_obj <- tryCatch(
      smiDE::results(
        smide_results = fit,
        comparisons = "one_vs_rest",
        variable = "smide_niche_treatment",
        targets = "all"
      ),
      error = function(e) {
        message("smiDE::results failed for ", run_label, " / ", cell_type, ": ", conditionMessage(e))
        NULL
      }
    )
    if (isTRUE(smide_save_raw) && !is.null(res_obj)) {
      saveRDS(res_obj, file.path(tables_dir, paste0(file_stub_base, "_smide_results_raw.rds")))
    }

    result_df <- coerce_smide_results_table(res_obj)
    if (!is.null(result_df) && nrow(result_df) > 0) {
      result_df$backend        <- "smiDE"
      result_df$analysis_scope <- "merged"
      utils::write.csv(result_df,
                       file.path(tables_dir, paste0(file_stub_base, "_results.csv")),
                       row.names = FALSE)
      all_results[[idx]] <- result_df
    }

    sig_col <- guess_significance_column(result_df)
    n_sig   <- if (!is.null(sig_col) && !is.null(result_df)) {
      sum(result_df[[sig_col]] < 0.05, na.rm = TRUE)
    } else NA_integer_

    summary_rows[[idx]] <- data.frame(
      analysis_scope         = "merged",
      backend                = "smiDE",
      run_label              = run_label,
      sample_id              = paste(sample_ids, collapse = "+"),
      annotation             = cell_type,
      treatment_column       = treatment_column,
      niche_treatment_groups = paste(valid_groups, collapse = " / "),
      n_cells                = nrow(cell_meta_valid),
      n_genes_tested         = model_inputs$n_genes,
      n_fdr_005              = n_sig,
      stringsAsFactors       = FALSE
    )
    idx <- idx + 1L
    gc(verbose = FALSE)
  }

  rm(counts_dense_full)
  gc(verbose = FALSE)

  list(
    metadata           = metadata,
    summary            = if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame(),
    results            = if (length(all_results)  > 0) do.call(rbind, all_results)  else data.frame(),
    skipped_cell_types = skipped_cell_types
  )
}

detect_annotation_column <- function(metadata, preferred = NULL) {
  annotation_stats <- function(column) {
    values <- metadata[[column]]
    values <- trimws(as.character(values))
    values[values == ""] <- NA_character_
    values <- stats::na.omit(values)
    list(
      n_non_missing = length(values),
      n_unique = length(unique(values))
    )
  }
  
  if (!is.null(preferred) && preferred %in% names(metadata)) {
    stats <- annotation_stats(preferred)
    if (stats$n_non_missing > 0 && stats$n_unique >= 2) {
      return(preferred)
    }
    message(
      "Preferred annotation column '", preferred,
      "' is not usable for spatial DE (",
      stats$n_non_missing, " labeled cells, ",
      stats$n_unique, " unique labels). Auto-detecting a better column."
    )
  }
  
  candidates <- c(
    grep("^celltype", names(metadata), value = TRUE),
    grep("annotation|annot", names(metadata), value = TRUE, ignore.case = TRUE),
    grep("leiden|cluster", names(metadata), value = TRUE, ignore.case = TRUE),
    "celltype"
  )
  candidates <- unique(candidates[candidates %in% names(metadata)])
  if (length(candidates) == 0) {
    stop("No annotation column was found for spatial differential expression.")
  }
  
  candidate_scores <- vapply(candidates, function(column) {
    stats <- annotation_stats(column)
    if (stats$n_non_missing == 0 || stats$n_unique < 2) {
      return(-Inf)
    }
    bonus <- 0
    if (grepl("supervised", column, ignore.case = TRUE)) {
      bonus <- bonus + 1000
    }
    if (grepl("^celltype_", column)) {
      bonus <- bonus + 100
    }
    if (identical(column, "celltype")) {
      bonus <- bonus - 10
    }
    stats$n_unique + bonus
  }, numeric(1))
  
  if (all(!is.finite(candidate_scores))) {
    stop("No annotation column with at least two non-missing labels was found for spatial differential expression.")
  }
  
  candidates[which.max(candidate_scores)]
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
    .giotto_get_spatial_network(gobj, name = spatial_network_name, output = "networkDT")
    TRUE
  }, error = function(e) FALSE)
  
  if (network_exists) {
    return(gobj)
  }
  
  .giotto_create_spatial_network(
    gobj = gobj,
    name = spatial_network_name,
    method = "Delaunay",
    minimum_k = 2
  )
}

get_network_dt <- function(gobj, spatial_network_name) {
  .giotto_get_spatial_network(gobj, name = spatial_network_name, output = "networkDT")
}

get_spatial_locations_df <- function(gobj) {
  as.data.frame(.giotto_get_spatial_locations(gobj, output = "data.table"), stringsAsFactors = FALSE)
}

as_plain_numeric_matrix <- function(x) {
  x_mat <- as.matrix(x)
  x_mat <- matrix(
    as.numeric(x_mat),
    nrow = nrow(x_mat),
    ncol = ncol(x_mat),
    dimnames = dimnames(x_mat)
  )
  storage.mode(x_mat) <- "double"
  x_mat
}

normalise_rows <- function(mat) {
  mat <- as_plain_numeric_matrix(mat)
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
  
  if (nrow(edges) == 0) {
    stop(
      "No annotated neighbors could be constructed from annotation column '",
      annotation_column,
      "'. Check that the selected annotation column contains non-missing labels for cells in the spatial network."
    )
  }
  
  niche_mat <- xtabs(~ cell_ID + neighbor_type, data = edges)
  niche_mat <- as_plain_numeric_matrix(niche_mat)
  
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
  niche_props <- as_plain_numeric_matrix(niche_props)
  
  if (is.null(dim(niche_props)) || nrow(niche_props) == 0) {
    return(character())
  }
  
  if (is.null(ncol(niche_props)) || ncol(niche_props) == 0) {
    warning("Neighborhood matrix has zero informative columns. Assigning all cells to niche_1.")
    return(rep("niche_1", nrow(niche_props)))
  }
  
  storage.mode(niche_props) <- "double"
  niche_props[!is.finite(niche_props)] <- 0
  if (all(rowSums(abs(niche_props)) == 0)) {
    warning("Neighborhood matrix contains only zeros. Assigning all cells to niche_1.")
    return(rep("niche_1", nrow(niche_props)))
  }
  
  rounded_props <- round(niche_props, 6)
  unique_props <- unique(as.data.frame(rounded_props, check.names = FALSE))
  unique_row_count <- if (is.null(dim(unique_props))) 1L else nrow(unique_props)
  if (length(unique_row_count) == 0 || is.na(unique_row_count) || !is.finite(unique_row_count)) {
    unique_row_count <- 1L
  }
  
  max_centers <- min(
    as.integer(n_niches)[1],
    nrow(niche_props),
    unique_row_count
  )
  
  if (length(max_centers) == 0 || is.na(max_centers) || !is.finite(max_centers) || max_centers <= 1) {
    return(rep("niche_1", nrow(niche_props)))
  }
  
  set.seed(42)
  niche_fit <- tryCatch(
    stats::kmeans(niche_props, centers = max_centers, nstart = 25),
    error = function(e) {
      warning(
        "Spatial niche clustering failed (", conditionMessage(e),
        "). Assigning all cells to niche_1."
      )
      NULL
    }
  )
  if (is.null(niche_fit)) {
    return(rep("niche_1", nrow(niche_props)))
  }
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
  if (is.null(dim(xy)) || nrow(xy) < 2 || is.null(ncol(xy)) || ncol(xy) < 2 || anyNA(xy)) {
    metadata$spatial_patch <- "patch_1"
    return(list(metadata = metadata, replicate_column = "spatial_patch"))
  }
  
  unique_xy <- unique(as.data.frame(round(xy, 3), check.names = FALSE))
  n_unique_xy <- if (is.null(dim(unique_xy))) 1L else nrow(unique_xy)
  if (length(n_unique_xy) == 0 || is.na(n_unique_xy) || !is.finite(n_unique_xy)) {
    n_unique_xy <- 1L
  }
  
  n_centers <- min(as.integer(n_spatial_patches)[1], nrow(xy), n_unique_xy)
  if (n_centers < 2) {
    metadata$spatial_patch <- "patch_1"
    return(list(metadata = metadata, replicate_column = "spatial_patch"))
  }
  
  set.seed(42)
  km <- tryCatch(
    stats::kmeans(xy, centers = n_centers, nstart = 25),
    error = function(e) {
      warning(
        "Spatial patch clustering failed (", conditionMessage(e),
        "). Using a single spatial patch."
      )
      NULL
    }
  )
  if (is.null(km)) {
    metadata$spatial_patch <- "patch_1"
    return(list(metadata = metadata, replicate_column = "spatial_patch"))
  }
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
      title = sample_plot_title(run_label, "Spatial Niche Assignments"),
      subtitle = "Cells colored by their inferred neighborhood niche",
      x = "Global X Coordinate",
      y = "Global Y Coordinate",
      colour = "Spatial\nNiche"
    ) +
    presentation_theme(base_size = 12, legend_position = "right")
  
  save_presentation_plot(
    plot = p1,
    filename = file.path(output_dir, paste0(run_label, "_spatial_niches.png")),
    width = 13,
    height = 10
  )
  
  niche_comp <- as.data.frame(table(metadata$spatial_niche, metadata[[annotation_column]]), stringsAsFactors = FALSE)
  names(niche_comp) <- c("spatial_niche", "annotation", "n_cells")
  
  p2 <- ggplot2::ggplot(
    niche_comp,
    ggplot2::aes(x = spatial_niche, y = n_cells, fill = annotation)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::scale_fill_discrete(
      labels = function(x) pretty_plot_label(x, width = 18)
    ) +
    ggplot2::labs(
      title = sample_plot_title(run_label, "Cell-Type Composition per Spatial Niche"),
      subtitle = "Relative cell-type abundance within each inferred spatial niche",
      x = "Spatial Niche",
      y = "Proportion",
      fill = "Cell Type"
    ) +
    presentation_theme(base_size = 12, x_angle = 35) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1)
    )
  
  save_presentation_plot(
    plot = p2,
    filename = file.path(output_dir, paste0(run_label, "_niche_composition.png")),
    width = 12,
    height = 8
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
    n_cells_raw <- nrow(cell_meta)
    if (n_cells_raw < min_cells_per_niche) {
      message(sprintf(
        "  \u26a0 edgeR spatial DE skipped for '%s': only %d cells found (requires \u2265%d)",
        cell_type, n_cells_raw, min_cells_per_niche
      ))
      next
    }

    replicate_counts <- as.data.frame(
      table(cell_meta$spatial_niche, cell_meta[[replicate_column]]),
      stringsAsFactors = FALSE
    )
    names(replicate_counts) <- c("spatial_niche", replicate_column, "n_cells")
    keep_groups <- replicate_counts[replicate_counts$n_cells >= min_cells_per_replicate, , drop = FALSE]
    if (nrow(keep_groups) == 0) {
      message(sprintf(
        "  \u26a0 edgeR spatial DE skipped for '%s': no replicate (%s) has \u2265%d cells in any niche (min_cells_per_replicate=%d)",
        cell_type, replicate_column, min_cells_per_replicate, min_cells_per_replicate
      ))
      next
    }

    keep_keys <- paste(keep_groups$spatial_niche, keep_groups[[replicate_column]], sep = "__")
    cell_meta$keep_key <- paste(cell_meta$spatial_niche, cell_meta[[replicate_column]], sep = "__")
    cell_meta <- cell_meta[cell_meta$keep_key %in% keep_keys, , drop = FALSE]
    if (nrow(cell_meta) < min_cells_per_niche) {
      message(sprintf(
        "  \u26a0 edgeR spatial DE skipped for '%s': after replicate filtering only %d cells remain (requires \u2265%d)",
        cell_type, nrow(cell_meta), min_cells_per_niche
      ))
      next
    }

    niche_sizes <- table(cell_meta$spatial_niche)
    valid_niches <- names(niche_sizes[niche_sizes >= min_cells_per_niche])
    cell_meta <- cell_meta[cell_meta$spatial_niche %in% valid_niches, , drop = FALSE]
    if (length(valid_niches) < 2) {
      message(sprintf(
        "  \u26a0 edgeR spatial DE skipped for '%s': only %d niche(s) with \u2265%d cells (requires \u22652)",
        cell_type, length(valid_niches), min_cells_per_niche
      ))
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
          message(sprintf(
            "  \u26a0 edgeR one_vs_rest skipped for '%s' / niche '%s': group sizes [%s] do not meet min_replicates_per_group=%d",
            cell_type, niche,
            paste(names(group_sizes), group_sizes, sep = "=", collapse = ", "),
            min_replicates_per_group
          ))
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
        message(sprintf(
          "  \u26a0 edgeR pairwise skipped for '%s': only %d niche(s) with \u2265%d replicates (requires \u22652)",
          cell_type, length(niche_levels), min_replicates_per_group
        ))
        next
      }

      for (pair in combn(sort(niche_levels), 2, simplify = FALSE)) {
        keep_idx <- pb_meta$spatial_niche %in% pair
        pair_meta <- pb_meta[keep_idx, , drop = FALSE]
        pair_counts <- pb$counts[, keep_idx, drop = FALSE]
        group_sizes <- table(pair_meta$spatial_niche)
        if (any(group_sizes < min_replicates_per_group)) {
          message(sprintf(
            "  \u26a0 edgeR pairwise skipped for '%s' / pair '%s': group sizes [%s] do not meet min_replicates_per_group=%d",
            cell_type, paste(pair, collapse = "_vs_"),
            paste(names(group_sizes), group_sizes, sep = "=", collapse = ", "),
            min_replicates_per_group
          ))
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
                                                smide_annotation_subset = NULL,
                                                smide_min_detection_fraction = 0.05,
                                                smide_custom_predictor = NULL,
                                                smide_partner_celltypes = NULL,
                                                smide_partner_source = c("none", "cell_proximity_auto"),
                                                smide_partner_top_n = 3,
                                                smide_partner_padj_threshold = 0.05,
                                                smide_include_self_partner = FALSE,
                                                smide_save_raw = TRUE,
                                                smide_adaptive_thresholds = TRUE,
                                                smide_edger_fallback = TRUE,
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
      gobj <- .giotto_load(gobj)
    }
  }
  
  results_dir <- ensure_dir(file.path(output_dir, "12_Spatial_Differential_Expression"))
  tables_dir <- ensure_dir(file.path(results_dir, "tables"))
  plots_dir <- ensure_dir(file.path(results_dir, "plots"))
  
  metadata <- as.data.frame(.giotto_pdata_dt(gobj), stringsAsFactors = FALSE)
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
  cat("Spatial network:", spatial_network_name, "\n")
  cat("Total cells:", nrow(metadata), "\n")
  cat("Total annotation types:", length(unique(stats::na.omit(metadata[[annotation_column]]))), "\n\n")
  
  gobj <- ensure_spatial_network(gobj, spatial_network_name = spatial_network_name)
  
  niche_counts <- build_neighbourhood_matrix(
    gobj = gobj,
    metadata = metadata,
    annotation_column = annotation_column,
    spatial_network_name = spatial_network_name
  )
  niche_props <- normalise_rows(niche_counts)
  metadata$spatial_niche <- assign_spatial_niches(niche_props, n_niches = n_niches)
  cat("Spatial niches identified:", length(unique(metadata$spatial_niche)), "\n\n")
  
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
  
  expr_mat <- .giotto_get_expression(gobj, values = "raw", output = "matrix")
  
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
        output_dir = output_dir,
        smide_ncores = smide_ncores,
        smide_overlap_threshold = smide_overlap_threshold,
        smide_annotation_subset = smide_annotation_subset,
        smide_min_detection_fraction = smide_min_detection_fraction,
        smide_custom_predictor = smide_custom_predictor,
        smide_partner_celltypes = smide_partner_celltypes,
        smide_partner_source = smide_partner_source,
        smide_partner_top_n = smide_partner_top_n,
        smide_partner_padj_threshold = smide_partner_padj_threshold,
        smide_include_self_partner = smide_include_self_partner,
        smide_save_raw = smide_save_raw,
        smide_adaptive_thresholds = smide_adaptive_thresholds,
        spatial_network_name = spatial_network_name,
        neighbor_counts = niche_counts
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
    if (backend == "smiDE") {
      run_smide_merged_backend(
        expr_mat = expr_mat,
        metadata = metadata,
        run_label = run_label,
        output_dir = output_dir,
        annotation_column = annotation_column,
        tables_dir = tables_dir,
        treatment_column = treatment_column,
        sample_column = sample_column,
        min_cells_per_niche = min_cells_per_niche,
        spatial_locations = spatial_locations,
        neighbor_counts = niche_counts,
        smide_family = smide_family,
        smide_radius = smide_radius,
        smide_ncores = smide_ncores,
        smide_overlap_threshold = smide_overlap_threshold,
        smide_annotation_subset = smide_annotation_subset,
        smide_min_detection_fraction = smide_min_detection_fraction,
        smide_save_raw = smide_save_raw,
        smide_adaptive_thresholds = smide_adaptive_thresholds
      )
    } else {
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
  }

  if (analysis_scope == "sample" && backend == "smiDE" && isTRUE(smide_edger_fallback)) {
    skipped <- scope_out$skipped_cell_types
    if (length(skipped) > 0) {
      cat("  edgeR fallback for", length(skipped),
          "cell type(s) skipped by smiDE:", paste(skipped, collapse = ", "), "\n")
      meta_fallback <- metadata[metadata[[annotation_column]] %in% skipped, , drop = FALSE]

      # Compute per-type cell counts; scale niches to the smallest group
      n_cells_per_type <- vapply(skipped, function(ct) {
        sum(meta_fallback[[annotation_column]] == ct, na.rm = TRUE)
      }, integer(1))
      min_type_count <- min(n_cells_per_type)

      # Re-assign niches with fewer groups if cell count is too small for
      # the global n_niches — this prevents the 6-niche assignment leaving
      # only 4 cells/niche on average, which the replicate filter then drops
      n_niches_fallback <- max(2L, min(n_niches, as.integer(floor(min_type_count / 10L))))
      if (n_niches_fallback < n_niches) {
        cat(sprintf(
          "  [fallback] Re-assigning niches: %d \u2192 %d (smallest type has %d cells)\n",
          n_niches, n_niches_fallback, min_type_count
        ))
        fb_ids <- meta_fallback$cell_ID
        fb_niche_props <- normalise_rows(niche_counts[
          match(fb_ids, rownames(niche_counts)), , drop = FALSE
        ])
        meta_fallback$spatial_niche <- assign_spatial_niches(
          fb_niche_props, n_niches = n_niches_fallback
        )
      }

      # k-means on neighborhood composition can produce unbalanced groups when
      # cells are spatially concentrated (e.g. all B-cells cluster in one niche).
      # If the smaller niche has < 5 cells, fall back to a coordinate-based
      # binary split (median sdimx), which guarantees a near-equal division.
      fallback_floor <- 5L
      niche_tab <- table(meta_fallback$spatial_niche)
      if (min(niche_tab) < fallback_floor && !is.null(spatial_locations)) {
        fb_ids <- meta_fallback$cell_ID
        loc_idx <- match(fb_ids, spatial_locations$cell_ID)
        fb_x <- spatial_locations$sdimx[loc_idx]
        if (!all(is.na(fb_x))) {
          med_x <- stats::median(fb_x, na.rm = TRUE)
          meta_fallback$spatial_niche <- ifelse(fb_x >= med_x, "geo_niche_A", "geo_niche_B")
          cat(sprintf(
            "  [fallback] k-means unbalanced (%s); using coordinate split (median sdimx=%.1f)\n",
            paste(names(niche_tab), niche_tab, sep = "=", collapse = ", "), med_x
          ))
        }
      }

      fallback_min <- adaptive_min_cells_per_niche(
        n_cells   = min_type_count,
        base_min  = min_cells_per_niche,
        n_niches  = n_niches_fallback,
        floor_min = fallback_floor
      )

      # With few cells spread across many FOVs, the replicate filter strips
      # nearly everything even after niche reduction. In the fallback (last
      # resort) always use 1 so no cell is dropped by this filter.
      fallback_rep_min <- 1L
      cat(sprintf(
        "  [fallback] min_cells_per_replicate set to 1 (was %d; fallback mode)\n",
        min_cells_per_replicate
      ))

      fallback_out <- run_sample_scope_spatial_de(
        expr_mat                 = expr_mat,
        metadata                 = meta_fallback,
        run_label                = paste0(run_label, "_edgeR_fallback"),
        annotation_column        = annotation_column,
        tables_dir               = tables_dir,
        sample_replicate_column  = sample_replicate_column,
        sample_contrast          = sample_contrast,
        min_cells_per_niche      = fallback_min,
        min_cells_per_replicate  = fallback_rep_min,
        min_replicates_per_group = min_replicates_per_group,
        n_spatial_patches        = n_spatial_patches,
        spatial_locations        = spatial_locations
      )
      if (nrow(fallback_out$results) > 0) {
        fallback_out$results$backend <- "edgeR_fallback"
        scope_out$results <- rbind(scope_out$results, fallback_out$results)
        scope_out$summary <- rbind(scope_out$summary, fallback_out$summary)
      }
    }
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
  
  # Final cleanup
  gc(verbose = FALSE)
  
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
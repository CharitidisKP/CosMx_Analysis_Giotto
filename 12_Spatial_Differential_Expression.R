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

# Extract cell polygons as a tidy data.frame (geom, part, x, y, hole, cell_ID).
# Returns NULL when polygon data are not available (e.g. centroid-only
# Giotto objects). Mirrors the helpers in 07_Annotation.R / 10_CCI_Analysis.R.
.extract_polygon_df <- function(gobj) {
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

  # Defensive: drop any SystemControl probes that survived upstream QC.
  # 02_Quality_Control.R strips them but legacy checkpoints may retain them.
  sc_mask <- grepl("^SystemControl", gene_names, ignore.case = TRUE)
  if (any(sc_mask & keep)) {
    n_sc <- sum(sc_mask & keep)
    message(sprintf(
      "  ⚠ smiDE: dropping %d SystemControl probe(s) that survived QC",
      n_sc
    ))
  }
  keep <- keep & !sc_mask
  
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

# Diff input targets vs fitted output → print one summary line and persist
# the failed-gene vector next to the fit log. smiDE::smi_de() prints a
# "Fitting model to target X." line per gene (captured to disk); the only
# surfaced signal here is a single success/warning line.
.report_smide_failed_genes <- function(submitted,
                                       results_df,
                                       fit_log,
                                       context = "") {
  submitted <- unique(as.character(submitted))
  if (length(submitted) == 0) {
    return(invisible(NULL))
  }

  fitted <- character(0)
  if (!is.null(results_df) && is.data.frame(results_df)) {
    target_col <- intersect(c("target", "gene"), colnames(results_df))
    if (length(target_col) > 0) {
      fitted <- unique(as.character(results_df[[target_col[1]]]))
    }
  }

  failed <- setdiff(submitted, fitted)
  ctx <- if (nzchar(context)) paste0(" [", context, "]") else ""

  if (length(failed) == 0) {
    message(sprintf("  ✓ smiDE: all %d gene(s) fitted%s", length(submitted), ctx))
  } else {
    preview <- paste(utils::head(failed, 20), collapse = ", ")
    tail_suffix <- if (length(failed) > 20) ", …" else ""
    message(sprintf(
      "  ⚠ smiDE: %d/%d gene(s) failed to fit%s: %s%s",
      length(failed), length(submitted), ctx, preview, tail_suffix
    ))
    if (nzchar(fit_log)) {
      failed_path <- sub("\\.log$", "_failed_genes.txt", fit_log)
      try(writeLines(failed, failed_path), silent = TRUE)
    }
  }
  invisible(list(submitted = submitted, fitted = fitted, failed = failed))
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

# Guess the fold-change column. Returns NULL if none found.
.guess_fc_column <- function(result_df) {
  if (is.null(result_df) || nrow(result_df) == 0) return(NULL)
  candidates <- grep(
    "^log2fc$|^logfc$|^log2foldchange$|foldchange|estimate|^beta$|logFC",
    names(result_df), value = TRUE, ignore.case = TRUE
  )
  candidates <- unique(candidates[candidates %in% names(result_df)])
  if (length(candidates) == 0) NULL else candidates[1]
}

# Guess lower / upper confidence-interval columns. Returns a 2-char vector
# (lower, upper) or NULL if either end cannot be resolved. Recognises the
# column naming patterns used by smiDE::results() (`*.lower`, `*.upper`,
# `asymp.LCL`, `asymp.UCL`, `CI.L`, `CI.U`, etc.) and edgeR.
.guess_ci_columns <- function(result_df) {
  if (is.null(result_df) || nrow(result_df) == 0) return(NULL)
  nms <- names(result_df)
  lower_pats <- c("lower.?cl$", "asymp\\.lcl$", "ci\\.?l$", "ci_?low$",
                  "\\.lower$", "^lower$", "lcl$", "lwr$", "lo$", "2\\.5")
  upper_pats <- c("upper.?cl$", "asymp\\.ucl$", "ci\\.?u$", "ci_?up$",
                  "\\.upper$", "^upper$", "ucl$", "upr$", "hi$", "97\\.5")
  lo <- unique(unlist(lapply(lower_pats,
                             function(p) grep(p, nms, value = TRUE, ignore.case = TRUE))))
  up <- unique(unlist(lapply(upper_pats,
                             function(p) grep(p, nms, value = TRUE, ignore.case = TRUE))))
  if (length(lo) == 0 || length(up) == 0) return(NULL)
  c(lo[1], up[1])
}

# Guess the gene-identifier column. Returns NULL if no obvious candidate.
.guess_gene_column <- function(result_df) {
  if (is.null(result_df) || nrow(result_df) == 0) return(NULL)
  candidates <- grep(
    "^gene$|^target$|^feature$|^symbol$|^gene_symbol$|^gene_name$|^genes$",
    names(result_df), value = TRUE, ignore.case = TRUE
  )
  candidates <- unique(candidates[candidates %in% names(result_df)])
  if (length(candidates) == 0) NULL else candidates[1]
}

# Emit smiDE-paper-style DE visualisations (volcano + top-hits bar plot) plus
# a lightweight FDR-distribution diagnostic. Publication plots go into a
# `plots/` subdirectory; the FDR histogram stays in `diagnostics/`.
# Graceful no-op when expected columns are missing.
.save_spatial_de_plots <- function(result_df, output_dir, file_stub,
                                   fdr_threshold = 0.05,
                                   lfc_threshold = log2(1.5),
                                   top_n_label = 10,
                                   top_n_bar = 20,
                                   top_n_forest = 15,
                                   top_n_overlay = 4,
                                   expr_mat = NULL,
                                   polygon_df = NULL) {
  tryCatch({
    if (is.null(result_df) || nrow(result_df) == 0) return(invisible(NULL))
    sig_col  <- guess_significance_column(result_df)
    fc_col   <- .guess_fc_column(result_df)
    gene_col <- .guess_gene_column(result_df)

    plots_dir <- file.path(output_dir, "plots")
    diag_dir  <- file.path(output_dir, "diagnostics")
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(diag_dir,  recursive = TRUE, showWarnings = FALSE)

    # --- 1. FDR / p-value histogram (diagnostic) ---
    if (!is.null(sig_col)) {
      sig_vals <- as.numeric(result_df[[sig_col]])
      sig_df <- data.frame(value = sig_vals[!is.na(sig_vals)])
      if (nrow(sig_df) > 0) {
        p_fdr <- ggplot2::ggplot(sig_df, ggplot2::aes(x = value)) +
          ggplot2::geom_histogram(bins = 40, fill = "#4C72B0",
                                  colour = "white", linewidth = 0.2) +
          ggplot2::geom_vline(xintercept = fdr_threshold, linetype = "dashed",
                              colour = "#C44E52") +
          ggplot2::labs(
            title    = paste0(file_stub, " - ", sig_col, " distribution"),
            subtitle = sprintf("n = %d genes; %d with %s < %g",
                               nrow(sig_df),
                               sum(sig_df$value < fdr_threshold, na.rm = TRUE),
                               sig_col, fdr_threshold),
            x = sig_col, y = "Genes"
          ) +
          presentation_theme(base_size = 11)
        save_presentation_plot(
          plot     = p_fdr,
          filename = file.path(diag_dir, paste0(file_stub, "_fdr_hist.png")),
          width    = 8, height = 5, dpi = 300
        )
      }
    }

    if (is.null(sig_col) || is.null(fc_col)) return(invisible(NULL))

    gene_vec <- if (!is.null(gene_col)) {
      as.character(result_df[[gene_col]])
    } else if (!is.null(rownames(result_df)) &&
               !identical(rownames(result_df), as.character(seq_len(nrow(result_df))))) {
      rownames(result_df)
    } else {
      as.character(seq_len(nrow(result_df)))
    }

    de_df <- data.frame(
      gene = gene_vec,
      fc   = as.numeric(result_df[[fc_col]]),
      sig  = as.numeric(result_df[[sig_col]]),
      stringsAsFactors = FALSE
    )
    de_df <- de_df[is.finite(de_df$fc) & is.finite(de_df$sig) & de_df$sig > 0,
                   , drop = FALSE]
    if (nrow(de_df) < 1) return(invisible(NULL))

    # --- 2. Volcano plot (smiDE-aligned publication figure) ---
    if (nrow(de_df) > 5) {
      vol_df <- de_df
      vol_df$neglog10 <- -log10(vol_df$sig)
      vol_df$status <- ifelse(vol_df$sig >= fdr_threshold, "ns",
                       ifelse(vol_df$fc >=  lfc_threshold, "up",
                       ifelse(vol_df$fc <= -lfc_threshold, "down", "ns")))
      vol_df$status <- factor(vol_df$status, levels = c("up", "down", "ns"))

      label_candidates <- vol_df[vol_df$status %in% c("up", "down"), , drop = FALSE]
      label_candidates <- label_candidates[order(label_candidates$sig), , drop = FALSE]
      label_df <- utils::head(label_candidates, top_n_label)

      p_vol <- ggplot2::ggplot(vol_df,
          ggplot2::aes(x = fc, y = neglog10, colour = status)) +
        ggplot2::geom_point(size = 1.3, alpha = 0.75, stroke = 0) +
        ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                            linetype = "dashed", colour = "grey60",
                            linewidth = 0.3) +
        ggplot2::geom_hline(yintercept = -log10(fdr_threshold),
                            linetype = "dashed", colour = "grey60",
                            linewidth = 0.3) +
        ggplot2::scale_colour_manual(
          values = c(up = "#C44E52", down = "#4C72B0", ns = "grey75"),
          breaks = c("up", "down", "ns"),
          labels = c(
            sprintf("Up (%s >= %.2f)",  fc_col, lfc_threshold),
            sprintf("Down (%s <= -%.2f)", fc_col, lfc_threshold),
            "Not significant"
          ),
          name = NULL, drop = FALSE, na.value = "grey75"
        ) +
        ggplot2::labs(
          title    = paste0(file_stub, " volcano"),
          subtitle = sprintf("%s < %g and |%s| >= %.2f flagged",
                             sig_col, fdr_threshold, fc_col, lfc_threshold),
          x = fc_col,
          y = paste0("-log10(", sig_col, ")")
        ) +
        presentation_theme(base_size = 11)

      if (nrow(label_df) > 0) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p_vol <- p_vol + ggrepel::geom_text_repel(
            data = label_df,
            mapping = ggplot2::aes(label = gene),
            size = 3, max.overlaps = Inf, box.padding = 0.4,
            segment.colour = "grey50", min.segment.length = 0,
            show.legend = FALSE, inherit.aes = TRUE
          )
        } else {
          p_vol <- p_vol + ggplot2::geom_text(
            data = label_df,
            mapping = ggplot2::aes(label = gene),
            size = 3, vjust = -0.5, show.legend = FALSE,
            inherit.aes = TRUE
          )
        }
      }
      save_presentation_plot(
        plot     = p_vol,
        filename = file.path(plots_dir, paste0(file_stub, "_volcano.png")),
        width    = 10, height = 8, dpi = 300
      )
    }

    # --- 3. Top-hits ranked bar plot (smiDE-aligned effect-size figure) ---
    hits_df <- de_df[de_df$sig < fdr_threshold, , drop = FALSE]
    if (nrow(hits_df) > 0) {
      hits_df <- hits_df[order(abs(hits_df$fc), decreasing = TRUE), , drop = FALSE]
      hits_df <- utils::head(hits_df, top_n_bar)
      hits_df$direction <- ifelse(hits_df$fc >= 0, "Up", "Down")
      hits_df$gene <- factor(hits_df$gene,
                             levels = hits_df$gene[order(hits_df$fc)])
      p_hits <- ggplot2::ggplot(hits_df,
          ggplot2::aes(x = fc, y = gene, fill = direction)) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_vline(xintercept = 0, colour = "grey30",
                            linewidth = 0.3) +
        ggplot2::scale_fill_manual(
          values = c(Up = "#C44E52", Down = "#4C72B0"),
          name = "Direction"
        ) +
        ggplot2::labs(
          title    = paste0(file_stub, " - top ", nrow(hits_df), " hits"),
          subtitle = sprintf("Ranked by |%s|; %s < %g",
                             fc_col, sig_col, fdr_threshold),
          x = fc_col, y = NULL
        ) +
        presentation_theme(base_size = 11)
      save_presentation_plot(
        plot     = p_hits,
        filename = file.path(plots_dir, paste0(file_stub, "_top_hits.png")),
        width    = 9,
        height   = max(5, 0.28 * nrow(hits_df) + 2),
        dpi      = 300
      )
    }

    # --- 4. Forest plot of top-N effect sizes with 95% CI ---
    # smiDE Fig 1i / 4f-g style: per-gene estimate +/- CI, one row per top gene.
    ci_cols <- .guess_ci_columns(result_df)
    if (!is.null(gene_col) && !is.null(ci_cols)) {
      forest_df <- data.frame(
        gene     = as.character(result_df[[gene_col]]),
        estimate = as.numeric(result_df[[fc_col]]),
        lower    = as.numeric(result_df[[ci_cols[1]]]),
        upper    = as.numeric(result_df[[ci_cols[2]]]),
        sig      = as.numeric(result_df[[sig_col]]),
        stringsAsFactors = FALSE
      )
      forest_df <- forest_df[is.finite(forest_df$estimate) &
                               is.finite(forest_df$lower) &
                               is.finite(forest_df$upper) &
                               is.finite(forest_df$sig),
                             , drop = FALSE]
      sig_forest <- forest_df[forest_df$sig < fdr_threshold, , drop = FALSE]
      if (nrow(sig_forest) > 0) {
        sig_forest <- sig_forest[order(abs(sig_forest$estimate), decreasing = TRUE),
                                 , drop = FALSE]
        sig_forest <- utils::head(sig_forest, top_n_forest)
        sig_forest$direction <- ifelse(sig_forest$estimate >= 0, "Up", "Down")
        sig_forest$gene <- factor(sig_forest$gene,
                                  levels = sig_forest$gene[order(sig_forest$estimate)])
        p_forest <- ggplot2::ggplot(sig_forest,
            ggplot2::aes(x = estimate, y = gene, colour = direction)) +
          ggplot2::geom_vline(xintercept = 0, colour = "grey40",
                              linewidth = 0.3, linetype = "dashed") +
          ggplot2::geom_errorbarh(
            ggplot2::aes(xmin = lower, xmax = upper),
            height = 0.2, linewidth = 0.5
          ) +
          ggplot2::geom_point(size = 2.2) +
          ggplot2::scale_colour_manual(
            values = c(Up = "#C44E52", Down = "#4C72B0"),
            name = "Direction"
          ) +
          ggplot2::labs(
            title    = paste0(file_stub, " - top ", nrow(sig_forest),
                              " hits (effect +/- 95% CI)"),
            subtitle = sprintf("Lower/upper from %s / %s; %s < %g",
                               ci_cols[1], ci_cols[2], sig_col, fdr_threshold),
            x = fc_col, y = NULL
          ) +
          presentation_theme(base_size = 11)
        save_presentation_plot(
          plot     = p_forest,
          filename = file.path(plots_dir, paste0(file_stub, "_forest.png")),
          width    = 9,
          height   = max(5, 0.30 * nrow(sig_forest) + 2),
          dpi      = 300
        )
      }
    }

    # --- 5. Spatial overlay of top-N DE genes (smiDE Fig 5g / 6f-g style) ---
    # Optional — only emitted when the caller supplies an expression matrix
    # and a polygon data.frame. Delegates to `.save_spatial_de_gene_overlay()`
    # which builds a per-gene polygon heatmap and patchworks the panels.
    if (!is.null(expr_mat) && !is.null(polygon_df) && !is.null(gene_col) &&
        top_n_overlay > 0) {
      overlay_hits <- de_df[de_df$sig < fdr_threshold, , drop = FALSE]
      if (nrow(overlay_hits) > 0) {
        overlay_hits <- overlay_hits[order(overlay_hits$sig), , drop = FALSE]
        up_hits   <- utils::head(overlay_hits[overlay_hits$fc > 0, , drop = FALSE],
                                 ceiling(top_n_overlay / 2))
        down_hits <- utils::head(overlay_hits[overlay_hits$fc < 0, , drop = FALSE],
                                 floor(top_n_overlay / 2))
        overlay_genes <- unique(c(up_hits$gene, down_hits$gene))
        # Backfill if one direction was empty
        if (length(overlay_genes) < top_n_overlay) {
          fill <- utils::head(overlay_hits$gene,
                              top_n_overlay - length(overlay_genes))
          overlay_genes <- unique(c(overlay_genes, fill))
        }
        if (length(overlay_genes) > 0) {
          .save_spatial_de_gene_overlay(
            genes       = overlay_genes,
            expr_mat    = expr_mat,
            polygon_df  = polygon_df,
            plots_dir   = plots_dir,
            file_stub   = file_stub
          )
        }
      }
    }
  }, error = function(e) {
    message("  Spatial DE plots skipped: ", conditionMessage(e))
  })
  invisible(NULL)
}

# Spatial overlay of gene expression on cell polygons — one panel per gene,
# patchworked into a single PNG. Renders all cells (no subsetting); the fill
# gradient is log1p of the raw counts found in `expr_mat`. Gracefully skipped
# if patchwork is unavailable or polygon/expression rows cannot be resolved.
.save_spatial_de_gene_overlay <- function(genes,
                                          expr_mat,
                                          polygon_df,
                                          plots_dir,
                                          file_stub,
                                          panel_ncol = NULL) {
  tryCatch({
    if (length(genes) == 0) return(invisible(NULL))
    if (!is.data.frame(polygon_df) ||
        !all(c("x", "y", "cell_ID") %in% names(polygon_df))) return(invisible(NULL))
    genes <- unique(as.character(genes))
    genes <- genes[genes %in% rownames(expr_mat)]
    if (length(genes) == 0) return(invisible(NULL))

    if (!"poly_group" %in% names(polygon_df)) {
      group_parts <- if (all(c("geom", "part") %in% names(polygon_df))) {
        paste(polygon_df$geom, polygon_df$part, sep = "_")
      } else {
        as.character(polygon_df$cell_ID)
      }
      polygon_df$poly_group <- group_parts
    }

    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

    panels <- lapply(genes, function(g) {
      vals <- as.numeric(expr_mat[g, , drop = TRUE])
      names(vals) <- colnames(expr_mat)
      expr_df <- data.frame(cell_ID = names(vals),
                            expr    = log1p(vals),
                            stringsAsFactors = FALSE)
      plot_df <- merge(polygon_df, expr_df, by = "cell_ID", all.x = TRUE)
      plot_df <- plot_df[!is.na(plot_df$expr), , drop = FALSE]
      if (nrow(plot_df) == 0) return(NULL)

      ggplot2::ggplot(plot_df,
          ggplot2::aes(x = x, y = y, group = poly_group, fill = expr)) +
        ggplot2::geom_polygon(colour = NA) +
        ggplot2::scale_fill_viridis_c(option = "magma", name = "log1p(count)") +
        ggplot2::coord_equal() +
        ggplot2::labs(title = g, x = NULL, y = NULL) +
        presentation_theme(base_size = 10) +
        ggplot2::theme(
          axis.text  = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          legend.key.height = grid::unit(0.4, "cm")
        )
    })
    panels <- panels[!vapply(panels, is.null, logical(1))]
    if (length(panels) == 0) return(invisible(NULL))

    n_panel <- length(panels)
    ncol_panel <- if (!is.null(panel_ncol)) panel_ncol else min(2, n_panel)
    nrow_panel <- ceiling(n_panel / ncol_panel)

    combined <- if (requireNamespace("patchwork", quietly = TRUE) && n_panel > 1) {
      Reduce(`+`, panels) + patchwork::plot_layout(ncol = ncol_panel)
    } else {
      panels[[1]]
    }

    save_presentation_plot(
      plot     = combined,
      filename = file.path(plots_dir, paste0(file_stub, "_top_gene_overlay.png")),
      width    = max(10, 6 * ncol_panel),
      height   = max(6, 5 * nrow_panel),
      dpi      = 300
    )
  }, error = function(e) {
    message("  Spatial DE gene overlay skipped: ", conditionMessage(e))
  })
  invisible(NULL)
}

# Run-level overview: counts of FDR < 0.05 genes per cell type and contrast.
# Replaces the old niche re-plot output that duplicated Script 09's niche maps.
.save_spatial_de_summary_plot <- function(summary_df, plots_dir, run_label) {
  tryCatch({
    if (is.null(summary_df) || nrow(summary_df) == 0) return(invisible(NULL))
    if (!"n_fdr_005" %in% names(summary_df)) return(invisible(NULL))

    df <- summary_df
    df$n_fdr_005 <- suppressWarnings(as.integer(df$n_fdr_005))
    df$n_fdr_005[is.na(df$n_fdr_005)] <- 0L
    if (all(df$n_fdr_005 == 0)) return(invisible(NULL))

    label_col <- intersect(c("annotation", "cell_type"), names(df))[1]
    if (is.na(label_col)) label_col <- names(df)[1]
    df$label <- as.character(df[[label_col]])

    facet_candidates <- intersect(
      c("predictor", "comparison", "spatial_niche", "niche_treatment_groups"),
      names(df)
    )
    facet_var <- if (length(facet_candidates) > 0) facet_candidates[1] else NA_character_

    ord <- order(df$n_fdr_005)
    df$label <- factor(df$label, levels = unique(df$label[ord]))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = n_fdr_005, y = label)) +
      ggplot2::geom_col(fill = "#4C72B0", width = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = n_fdr_005),
                         hjust = -0.2, size = 3, colour = "grey30") +
      ggplot2::labs(
        title    = sample_plot_title(run_label,
                                     "Spatial DE: significant genes per cell type"),
        subtitle = "Gene counts with FDR < 0.05 (each bar = one contrast)",
        x = "Significant genes (FDR < 0.05)",
        y = NULL
      ) +
      presentation_theme(base_size = 11) +
      ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 24, b = 10, l = 10))

    if (!is.na(facet_var)) {
      p <- p + ggplot2::facet_wrap(
        stats::as.formula(paste0("~", facet_var)),
        ncol = 1, scales = "free_y"
      )
    }

    n_bars <- nrow(df)
    save_presentation_plot(
      plot     = p,
      filename = file.path(plots_dir, paste0(run_label, "_spatial_de_summary.png")),
      width    = 11,
      height   = max(6, 0.35 * n_bars + 2),
      dpi      = 300
    )
  }, error = function(e) {
    message("  Spatial DE summary plot skipped: ", conditionMessage(e))
  })
  invisible(NULL)
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
                                     smide_adaptive_thresholds = TRUE,
                                     polygon_df = NULL) {
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
      "smide_annotation_subset: [\"B cell\"] in config.yaml.",
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
  skipped_records    <- list()

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
      skipped_records[[length(skipped_records) + 1L]] <- data.frame(
        cell_type = cell_type,
        n_cells   = n_cells_total,
        n_niches  = length(valid_niches),
        reason    = sprintf("below niche threshold (need \u22652 niches with \u2265%d cells)", effective_min),
        stage     = "niche_gate",
        stringsAsFactors = FALSE
      )
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
      skipped_cell_types <- c(skipped_cell_types, cell_type)
      skipped_records[[length(skipped_records) + 1L]] <- data.frame(
        cell_type = cell_type,
        n_cells   = nrow(cell_meta),
        n_niches  = length(valid_niches),
        reason    = "smiDE::pre_de() failed",
        stage     = "pre_de",
        stringsAsFactors = FALSE
      )
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
      skipped_cell_types <- c(skipped_cell_types, cell_type)
      skipped_records[[length(skipped_records) + 1L]] <- data.frame(
        cell_type = cell_type,
        n_cells   = nrow(cell_meta),
        n_niches  = length(valid_niches),
        reason    = sprintf("overlap threshold %.2f removed all targets", smide_overlap_threshold),
        stage     = "overlap_filter",
        stringsAsFactors = FALSE
      )
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
      skipped_cell_types <- c(skipped_cell_types, cell_type)
      skipped_records[[length(skipped_records) + 1L]] <- data.frame(
        cell_type = cell_type,
        n_cells   = nrow(cell_meta),
        n_niches  = length(valid_niches),
        reason    = sprintf("no genes passed min_detection_fraction=%.2f", smide_min_detection_fraction),
        stage     = "detection_filter",
        stringsAsFactors = FALSE
      )
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

      file_stub <- paste(file_stub_base, sanitize_smide_name(predictor_label), sep = "_")
      fit_log <- file.path(tables_dir, paste0(file_stub, "_smide_fit.log"))

      smide_targets_submitted <- model_inputs$targets
      fit <- tryCatch(
        withCallingHandlers(
          {
            capture.output(
              res <- smiDE::smi_de(
                assay_matrix = predictor_counts,
                metadata = predictor_meta,
                formula = stats::as.formula(paste("~", predictor_var)),
                pre_de_obj = pre_obj,
                groupVar = predictor_var,
                groupVar_levels = predictor_levels,
                nCores = smide_ncores,
                family = smide_family,
                targets = smide_targets_submitted,
                neighbor_expr_cell_type_metadata_colname = annotation_column,
                cellid_colname = "cell_ID"
              ),
              file = fit_log, append = FALSE, type = "output"
            )
            res
          },
          message = function(m) invokeRestart("muffleMessage")
        ),
        error = function(e) {
          message("Skipping smiDE model for ", run_label, " / ", cell_type, " / ", predictor_label, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(fit)) {
        next
      }

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

      # Summarise any genes that failed to fit (input vs output diff).
      .report_smide_failed_genes(
        submitted = smide_targets_submitted,
        results_df = result_df,
        fit_log = fit_log,
        context = paste(run_label, cell_type, predictor_label, sep = " / ")
      )

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
      .save_spatial_de_plots(result_df,
                             dirname(tables_dir),
                             paste0(file_stub, "_spatial_de"),
                             expr_mat = expr_mat,
                             polygon_df = polygon_df)

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

  if (length(skipped_records) > 0) {
    skipped_df <- do.call(rbind, skipped_records)
    skip_csv <- file.path(tables_dir, paste0(run_label, "_smide_skipped_summary.csv"))
    utils::write.csv(skipped_df, skip_csv, row.names = FALSE)
    cat(sprintf("\n[smiDE] %d/%d candidate cell type(s) skipped \u2014 see %s\n",
                length(skipped_records), length(candidate_cell_types), skip_csv))
  } else {
    cat(sprintf("\n[smiDE] 0/%d candidate cell types skipped\n",
                length(candidate_cell_types)))
  }

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
                                     smide_adaptive_thresholds = TRUE,
                                     polygon_df = NULL) {
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
  skipped_records    <- list()
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
      skipped_records[[length(skipped_records) + 1L]] <- data.frame(
        cell_type = cell_type,
        n_cells   = n_cells_total,
        n_niches  = length(valid_groups),
        reason    = sprintf("needs \u22652 niche-treatment groups with \u2265%d cells spanning \u22652 treatments (got %d groups, %d treatments)",
                            effective_min, length(valid_groups), n_treat_valid),
        stage     = "group_gate",
        stringsAsFactors = FALSE
      )
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
    if (is.null(pre_obj)) {
      skipped_cell_types <- c(skipped_cell_types, cell_type)
      skipped_records[[length(skipped_records) + 1L]] <- data.frame(
        cell_type = cell_type, n_cells = n_cells_total, n_niches = length(valid_groups),
        reason = "smiDE::pre_de() failed", stage = "pre_de",
        stringsAsFactors = FALSE
      )
      gc(verbose = FALSE); next
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
      skipped_cell_types <- c(skipped_cell_types, cell_type)
      skipped_records[[length(skipped_records) + 1L]] <- data.frame(
        cell_type = cell_type, n_cells = n_cells_total, n_niches = length(valid_groups),
        reason = sprintf("overlap threshold %.2f removed all targets", smide_overlap_threshold),
        stage = "overlap_filter",
        stringsAsFactors = FALSE
      )
      gc(verbose = FALSE); next
    }

    model_inputs <- prepare_smide_targets(
      cell_counts = cell_counts,
      targets = targets,
      min_detection_fraction = smide_min_detection_fraction
    )
    if (is.null(model_inputs$counts) || model_inputs$n_genes == 0) {
      message("No genes passed smiDE filters for ", run_label, " / ", cell_type)
      skipped_cell_types <- c(skipped_cell_types, cell_type)
      skipped_records[[length(skipped_records) + 1L]] <- data.frame(
        cell_type = cell_type, n_cells = n_cells_total, n_niches = length(valid_groups),
        reason = sprintf("no genes passed min_detection_fraction=%.2f", smide_min_detection_fraction),
        stage = "detection_filter",
        stringsAsFactors = FALSE
      )
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

    fit_log <- file.path(tables_dir, paste0(file_stub_base, "_smide_fit.log"))
    smide_targets_submitted <- model_inputs$targets

    fit <- tryCatch(
      withCallingHandlers(
        {
          capture.output(
            res <- smiDE::smi_de(
              assay_matrix = predictor_counts,
              metadata = cell_meta_valid,
              formula = stats::as.formula(formula_str),
              pre_de_obj = pre_obj,
              groupVar = "smide_niche_treatment",
              groupVar_levels = valid_groups,
              nCores = smide_ncores,
              family = smide_family,
              targets = smide_targets_submitted,
              neighbor_expr_cell_type_metadata_colname = annotation_column,
              cellid_colname = "cell_ID"
            ),
            file = fit_log, append = FALSE, type = "output"
          )
          res
        },
        message = function(m) invokeRestart("muffleMessage")
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

    .report_smide_failed_genes(
      submitted = smide_targets_submitted,
      results_df = result_df,
      fit_log = fit_log,
      context = paste(run_label, cell_type, sep = " / ")
    )

    if (!is.null(result_df) && nrow(result_df) > 0) {
      result_df$backend        <- "smiDE"
      result_df$analysis_scope <- "merged"
      utils::write.csv(result_df,
                       file.path(tables_dir, paste0(file_stub_base, "_results.csv")),
                       row.names = FALSE)
      .save_spatial_de_plots(result_df,
                             dirname(tables_dir),
                             paste0(file_stub_base, "_results"),
                             expr_mat = expr_mat,
                             polygon_df = polygon_df)
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

  if (length(skipped_records) > 0) {
    skipped_df <- do.call(rbind, skipped_records)
    skip_csv <- file.path(tables_dir, paste0(run_label, "_smide_merged_skipped_summary.csv"))
    utils::write.csv(skipped_df, skip_csv, row.names = FALSE)
    cat(sprintf("\n[smiDE merged] %d/%d candidate cell type(s) skipped \u2014 see %s\n",
                length(skipped_records), length(candidate_cell_types), skip_csv))
  } else {
    cat(sprintf("\n[smiDE merged] 0/%d candidate cell types skipped\n",
                length(candidate_cell_types)))
  }

  list(
    metadata           = metadata,
    summary            = if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame(),
    results            = if (length(all_results)  > 0) do.call(rbind, all_results)  else data.frame(),
    skipped_cell_types = skipped_cell_types
  )
}

# ==============================================================================
# Per-comparison + multi-stratifier smiDE wrapper.
#
# For each (stratifier ∈ stratify_by) × (comparison ∈ comparisons):
#   1. Subset metadata + counts + spatial_locations + neighbor_counts to cells
#      that match either group_a or group_b of the comparison.
#   2. Add a `smide_de_group` column ("A" or "B") to the subset metadata.
#   3. Call run_smide_merged_backend() with annotation_column set to the
#      stratifier (celltype or leiden cluster column) and treatment_column
#      set to "smide_de_group" — so smiDE's contrast variable is the per-
#      comparison group rather than the global treatment column.
#
# Outputs land in:
#   <output_dir>/12_Spatial_Differential_Expression/<stratifier>/<comparison>/
#
# Runtime caps:
#   smide_min_cells_per_stratum: skip strata with fewer cells than this.
#   smide_max_strata_per_comparison: cap the number of strata processed per
#     comparison (alphabetical order, with a clear log warning when exceeded).
#
# Descriptive comparisons (CART vs Control) are skipped — smiDE needs at
# least 2 biological replicates per side and Control has only 1. Step 13
# handles those contrasts as descriptive log2FC tables.
# ==============================================================================
run_smide_per_comparison_backend <- function(expr_mat,
                                             metadata,
                                             run_label,
                                             output_dir,
                                             stratifier_columns,
                                             comparisons,
                                             tables_dir,
                                             sample_column = "sample_id",
                                             min_cells_per_niche = 30,
                                             spatial_locations = NULL,
                                             neighbor_counts = NULL,
                                             smide_family = "nbinom2",
                                             smide_radius = 0.05,
                                             smide_ncores = 1,
                                             smide_overlap_threshold = NULL,
                                             smide_min_detection_fraction = 0.05,
                                             smide_save_raw = TRUE,
                                             smide_adaptive_thresholds = TRUE,
                                             smide_min_cells_per_stratum = 200,
                                             smide_max_strata_per_comparison = 30,
                                             polygon_df = NULL) {
  if (length(stratifier_columns) == 0L) {
    stop("run_smide_per_comparison_backend(): stratifier_columns is empty.")
  }
  if (length(comparisons) == 0L) {
    stop("run_smide_per_comparison_backend(): no comparisons supplied.")
  }
  stratifier_columns <- stratifier_columns[stratifier_columns %in% names(metadata)]
  if (length(stratifier_columns) == 0L) {
    stop("None of the requested stratifier_columns are present in cell metadata.")
  }

  per_comp_summary  <- list()
  per_comp_results  <- list()

  for (st_label in names(stratifier_columns)) {
    st_col <- stratifier_columns[[st_label]]
    cat("\n========================================\n")
    cat("smiDE | Stratifier: ", st_label, " (column: ", st_col, ")\n", sep = "")
    cat("========================================\n")

    for (cp in comparisons) {
      cp_label <- cp$label %||% "comparison"
      design_kind <- cp$design %||% "unpaired"
      if (design_kind == "descriptive") {
        cat("  - ", cp_label, ": descriptive comparison; skipping smiDE ",
            "(handled as log2FC-only by step 13)\n", sep = "")
        next
      }

      mask_a <- .smide_match_meta_filter(metadata, cp$group_a)
      mask_b <- .smide_match_meta_filter(metadata, cp$group_b)
      keep   <- mask_a | mask_b
      if (sum(keep) < 2L * min_cells_per_niche) {
        cat("  - ", cp_label, ": fewer than ",
            2L * min_cells_per_niche, " cells across both groups; skipping\n",
            sep = "")
        next
      }

      sub_meta <- metadata[keep, , drop = FALSE]
      sub_meta$smide_de_group <- ifelse(mask_a[keep], "A", "B")
      stratum_vec <- as.character(sub_meta[[st_col]])
      stratum_sizes <- table(stratum_vec[!is.na(stratum_vec) & nzchar(stratum_vec)])
      eligible_strata <- names(stratum_sizes[stratum_sizes >= smide_min_cells_per_stratum])

      if (length(eligible_strata) == 0L) {
        cat("  - ", cp_label, ": no ", st_label, " strata reach ",
            smide_min_cells_per_stratum, "-cell floor; skipping\n", sep = "")
        next
      }

      # Apply max_strata cap (alphabetical for determinism).
      eligible_strata <- sort(eligible_strata)
      if (length(eligible_strata) > smide_max_strata_per_comparison) {
        cat("  - ", cp_label, ": ", length(eligible_strata),
            " ", st_label, " strata exceed cap of ",
            smide_max_strata_per_comparison,
            "; running first ", smide_max_strata_per_comparison,
            " (alphabetical)\n", sep = "")
        eligible_strata <- head(eligible_strata,
                                smide_max_strata_per_comparison)
      }

      cat("  ▶ ", cp_label, " (", st_label, "): A=", sum(mask_a[keep]),
          " B=", sum(mask_b[keep]), " cells over ",
          length(eligible_strata), " stratum/strata\n", sep = "")

      cp_dir <- file.path(output_dir, "12_Spatial_Differential_Expression",
                          st_label, .smide_safe_name(cp_label))
      dir.create(cp_dir, recursive = TRUE, showWarnings = FALSE)
      cp_tables <- ensure_dir(file.path(cp_dir, "tables"))

      # Subset the auxiliary inputs to the comparison cells.
      cell_ids_keep <- sub_meta$cell_ID
      sub_expr <- expr_mat[, intersect(colnames(expr_mat), cell_ids_keep),
                           drop = FALSE]
      sub_spatial_locs <- if (!is.null(spatial_locations)) {
        spatial_locations[match(cell_ids_keep,
                                spatial_locations$cell_ID), , drop = FALSE]
      } else NULL
      sub_neighbor <- if (!is.null(neighbor_counts)) {
        neighbor_counts[intersect(rownames(neighbor_counts), cell_ids_keep),
                        , drop = FALSE]
      } else NULL

      iter_label <- paste(run_label, .smide_safe_name(cp_label),
                          st_label, sep = "_")

      iter_result <- tryCatch(
        run_smide_merged_backend(
          expr_mat = sub_expr,
          metadata = sub_meta,
          run_label = iter_label,
          output_dir = cp_dir,
          annotation_column = st_col,
          tables_dir = cp_tables,
          treatment_column = "smide_de_group",
          sample_column = sample_column,
          min_cells_per_niche = min_cells_per_niche,
          spatial_locations = sub_spatial_locs,
          neighbor_counts = sub_neighbor,
          smide_family = smide_family,
          smide_radius = smide_radius,
          smide_ncores = smide_ncores,
          smide_overlap_threshold = smide_overlap_threshold,
          smide_annotation_subset = eligible_strata,
          smide_min_detection_fraction = smide_min_detection_fraction,
          smide_save_raw = smide_save_raw,
          smide_adaptive_thresholds = smide_adaptive_thresholds,
          polygon_df = polygon_df
        ),
        error = function(e) {
          message("smiDE failed for ", iter_label, ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(iter_result)) next

      if (!is.null(iter_result$summary) && nrow(iter_result$summary) > 0) {
        iter_result$summary$comparison       <- cp_label
        iter_result$summary$stratifier_label <- st_label
        per_comp_summary[[length(per_comp_summary) + 1L]] <- iter_result$summary
      }
      if (!is.null(iter_result$results) && nrow(iter_result$results) > 0) {
        iter_result$results$comparison       <- cp_label
        iter_result$results$stratifier_label <- st_label
        per_comp_results[[length(per_comp_results) + 1L]] <- iter_result$results
      }
    }
  }

  combined_summary <- if (length(per_comp_summary) > 0) {
    do.call(rbind, per_comp_summary)
  } else data.frame()
  combined_results <- if (length(per_comp_results) > 0) {
    do.call(rbind, per_comp_results)
  } else data.frame()

  if (nrow(combined_summary) > 0) {
    utils::write.csv(
      combined_summary,
      file.path(tables_dir, paste0(run_label, "_smide_per_comparison_summary.csv")),
      row.names = FALSE
    )
  }
  if (nrow(combined_results) > 0) {
    utils::write.csv(
      combined_results,
      file.path(tables_dir, paste0(run_label, "_smide_per_comparison_results.csv")),
      row.names = FALSE
    )
  }

  list(
    metadata = metadata,
    summary  = combined_summary,
    results  = combined_results,
    skipped_cell_types = character(0)
  )
}

# Filename-safe slugify.
.smide_safe_name <- function(x) {
  s <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  s <- gsub("^_+|_+$", "", s)
  if (!nzchar(s)) "x" else s
}

# Mirror of step 13's pathway match-meta-filter, kept local to avoid
# cross-script sourcing.
.smide_match_meta_filter <- function(meta, filter_spec) {
  if (is.null(filter_spec) || length(filter_spec) == 0L) {
    return(rep(FALSE, nrow(meta)))
  }
  hits <- rep(TRUE, nrow(meta))
  for (key in names(filter_spec)) {
    want <- filter_spec[[key]]
    if (is.null(want) || length(want) == 0L) next
    if (!key %in% names(meta)) {
      cat("  ⚠ filter key '", key,
          "' not in cell metadata; treating as no match\n", sep = "")
      return(rep(FALSE, nrow(meta)))
    }
    have <- as.character(meta[[key]])
    hits <- hits & (have %in% as.character(want))
  }
  hits
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
                                        spatial_locations = NULL,
                                        polygon_df = NULL) {
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
        .save_spatial_de_plots(de_tbl,
                               dirname(tables_dir),
                               paste0(file_stub, "_spatial_de"),
                               expr_mat = expr_mat,
                               polygon_df = polygon_df)
        
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
        .save_spatial_de_plots(de_tbl,
                               dirname(tables_dir),
                               paste0(file_stub, "_spatial_de"),
                               expr_mat = expr_mat,
                               polygon_df = polygon_df)
        
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
                                        min_samples_per_group = 2,
                                        polygon_df = NULL) {
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
        .save_spatial_de_plots(de_tbl,
                               dirname(tables_dir),
                               paste0(file_stub, "_spatial_de"),
                               expr_mat = expr_mat,
                               polygon_df = polygon_df)
        
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
                                                # Per-comparison + multi-stratifier extension
                                                # (see run_smide_per_comparison_backend()):
                                                comparisons = NULL,
                                                stratifier_columns = NULL,
                                                smide_min_cells_per_stratum = 200,
                                                smide_max_strata_per_comparison = 30,
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

  # ── Metadata-column validation (T1.7): fail fast with clear messages ─────
  .validate_de_column <- function(col, label, require_multi_level = TRUE) {
    if (is.null(col) || !nzchar(col)) return(invisible(NULL))
    if (!col %in% names(metadata)) {
      stop(sprintf(
        "%s '%s' is not present in cell metadata. Available columns: %s",
        label, col, paste(names(metadata), collapse = ", ")
      ))
    }
    values <- metadata[[col]]
    n_na <- sum(is.na(values) | !nzchar(as.character(values)))
    n_levels <- length(unique(stats::na.omit(values)))
    if (require_multi_level && n_levels < 2) {
      stop(sprintf(
        "%s '%s' has %d non-NA unique value(s); need \u22652.",
        label, col, n_levels
      ))
    }
    if (n_na > 0 && n_na / length(values) > 0.1) {
      warning(sprintf(
        "%s '%s' has %d/%d (%.1f%%) NA/empty entries.",
        label, col, n_na, length(values), 100 * n_na / length(values)
      ), immediate. = TRUE)
    }
  }
  .validate_de_column(annotation_column, "annotation_column")
  if (analysis_scope == "merged") {
    .validate_de_column(treatment_column, "treatment_column")
    .validate_de_column(patient_column,   "patient_column", require_multi_level = FALSE)
  }

  # Safety net: legacy checkpoints annotated before the 07 normalization may
  # contain dot/underscore-separated labels (e.g. "B.cell" instead of
  # "B cell"). smiDE matches labels by exact string, so normalize in-memory
  # before validating. 07_Annotation.R now normalizes at source.
  ann_col_vec <- metadata[[annotation_column]]
  if (!is.null(ann_col_vec) &&
      any(grepl("[._]", as.character(ann_col_vec)), na.rm = TRUE)) {
    message("  ℹ Normalizing legacy dot/underscore labels in '",
            annotation_column, "' for smiDE compatibility")
    norm_vec <- gsub("\\s+", " ",
                     gsub("[._]+", " ", as.character(ann_col_vec)))
    norm_vec <- trimws(norm_vec)
    metadata[[annotation_column]] <- norm_vec
  }

  if (!is.null(smide_annotation_subset) && length(smide_annotation_subset) > 0) {
    ann_values <- unique(stats::na.omit(metadata[[annotation_column]]))
    missing_labels <- setdiff(smide_annotation_subset, ann_values)
    if (length(missing_labels) == length(smide_annotation_subset)) {
      stop(sprintf(
        "None of the smide_annotation_subset labels (%s) are present in '%s'. Available labels: %s",
        paste(smide_annotation_subset, collapse = ", "),
        annotation_column,
        paste(ann_values, collapse = ", ")
      ))
    } else if (length(missing_labels) > 0) {
      warning(sprintf(
        "smide_annotation_subset labels not in '%s': %s",
        annotation_column, paste(missing_labels, collapse = ", ")
      ), immediate. = TRUE)
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

  # Extract polygon geometry once so backends can emit smiDE-paper-style
  # spatial overlays of top DE genes. Silent fallback to NULL (no overlay)
  # when polygon data are unavailable on centroid-only objects.
  polygon_df <- tryCatch(.extract_polygon_df(gobj), error = function(e) NULL)
  if (is.null(polygon_df)) {
    cat("  Polygon geometry not available; spatial gene overlays will be skipped.\n")
  }

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
        neighbor_counts = niche_counts,
        polygon_df = polygon_df
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
        spatial_locations = spatial_locations,
        polygon_df = polygon_df
      )
    }
  } else {
    if (backend == "smiDE") {
      use_per_comparison <- !is.null(comparisons) && length(comparisons) > 0L &&
                            !is.null(stratifier_columns) &&
                            length(stratifier_columns) > 0L
      if (use_per_comparison) {
        cat("Using per-comparison + multi-stratifier smiDE wrapper:\n")
        cat("  Stratifiers : ", paste(names(stratifier_columns),
                                       collapse = ", "), "\n", sep = "")
        cat("  Comparisons : ",
            paste(vapply(comparisons, function(cp) cp$label %||% "?",
                         character(1)), collapse = ", "), "\n", sep = "")
        run_smide_per_comparison_backend(
          expr_mat = expr_mat,
          metadata = metadata,
          run_label = run_label,
          output_dir = output_dir,
          stratifier_columns = stratifier_columns,
          comparisons = comparisons,
          tables_dir = tables_dir,
          sample_column = sample_column,
          min_cells_per_niche = min_cells_per_niche,
          spatial_locations = spatial_locations,
          neighbor_counts = niche_counts,
          smide_family = smide_family,
          smide_radius = smide_radius,
          smide_ncores = smide_ncores,
          smide_overlap_threshold = smide_overlap_threshold,
          smide_min_detection_fraction = smide_min_detection_fraction,
          smide_save_raw = smide_save_raw,
          smide_adaptive_thresholds = smide_adaptive_thresholds,
          smide_min_cells_per_stratum = smide_min_cells_per_stratum,
          smide_max_strata_per_comparison = smide_max_strata_per_comparison,
          polygon_df = polygon_df
        )
      } else {
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
          smide_adaptive_thresholds = smide_adaptive_thresholds,
          polygon_df = polygon_df
        )
      }
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
        min_samples_per_group = min_samples_per_group,
        polygon_df = polygon_df
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
        spatial_locations        = spatial_locations,
        polygon_df               = polygon_df
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
    .save_spatial_de_summary_plot(summary_df, plots_dir, run_label)
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
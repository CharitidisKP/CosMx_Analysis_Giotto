#!/usr/bin/env Rscript
# ==============================================================================
# 11_B_Cell_Analysis.R
# Focused cell-type spatial interaction summaries
#
# By default this targets B cells and related subtypes for lupus nephritis
# analysis, but the bcell_regex parameter can be set to any cell type pattern
# via config.yaml (interaction.focus_celltype_regex).
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

arrange_feature_plots <- file.path(current_script_dir(), "Helper_Scripts",
                                    "Arrange_Feature_plots.R")
if (!exists("optimal_grid_dims") && file.exists(arrange_feature_plots)) {
  source(arrange_feature_plots)
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

.giotto_cell_proximity_enrichment <- function(gobj,
                                              spatial_network_name,
                                              cluster_column,
                                              number_of_simulations = 250) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("cellProximityEnrichment", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("cellProximityEnrichment", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("cellProximityEnrichment", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("cellProximityEnrichment", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("cellProximityEnrichment", mode = "function")
  }
  
  accessor(
    gobject = gobj,
    spatial_network_name = spatial_network_name,
    cluster_column = cluster_column,
    number_of_simulations = number_of_simulations
  )
}

.giotto_cell_proximity_heatmap <- function(gobj, cp_score, ...) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("cellProximityHeatmap", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("cellProximityHeatmap", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("cellProximityHeatmap", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("cellProximityHeatmap", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("cellProximityHeatmap", mode = "function")
  }
  
  accessor(gobject = gobj, CPscore = cp_score, ...)
}

.giotto_cell_proximity_network <- function(gobj, cp_score, ...) {
  accessor <- NULL
  
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("cellProximityNetwork", envir = asNamespace("Giotto"), inherits = FALSE)) {
    accessor <- get("cellProximityNetwork", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("cellProximityNetwork", envir = asNamespace("GiottoClass"), inherits = FALSE)) {
    accessor <- get("cellProximityNetwork", envir = asNamespace("GiottoClass"))
  } else {
    accessor <- get("cellProximityNetwork", mode = "function")
  }
  
  accessor(gobject = gobj, CPscore = cp_score, ...)
}

detect_annotation_column <- function(metadata, preferred = NULL) {
  if (!is.null(preferred) && preferred %in% names(metadata)) {
    return(preferred)
  }
  
  candidates <- c(
    "celltype",
    grep("^celltype", names(metadata), value = TRUE),
    grep("annotation|annot", names(metadata), value = TRUE, ignore.case = TRUE)
  )
  candidates <- unique(candidates[candidates %in% names(metadata)])
  
  if (length(candidates) == 0) {
    stop("No annotation column was found for focused cell-type analysis.")
  }
  
  candidates[1]
}

.matches_bcell_label <- function(labels, bcell_regex) {
  grepl(bcell_regex, labels, ignore.case = TRUE)
}

.parse_unified_interactions <- function(unified_int) {
  parts <- strsplit(as.character(unified_int), "--", fixed = TRUE)
  data.frame(
    source = vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)),
    target = vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

.filter_interactions_for_bcells <- function(tbl, bcell_regex) {
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(tbl)
  }
  
  if (all(c("source", "target") %in% names(tbl))) {
    return(dplyr::filter(
      tbl,
      .matches_bcell_label(source, bcell_regex) |
        .matches_bcell_label(target, bcell_regex)
    ))
  }
  
  bcell_cols <- grep("cell.*type|celltype|cluster", names(tbl), value = TRUE, ignore.case = TRUE)
  if (length(bcell_cols) >= 2) {
    return(dplyr::filter(
      tbl,
      .matches_bcell_label(rlang::.data[[bcell_cols[1]]], bcell_regex) |
        .matches_bcell_label(rlang::.data[[bcell_cols[2]]], bcell_regex)
    ))
  }
  
  if ("unified_int" %in% names(tbl)) {
    parsed <- .parse_unified_interactions(tbl$unified_int)
    keep <- .matches_bcell_label(parsed$source, bcell_regex) |
      .matches_bcell_label(parsed$target, bcell_regex)
    return(tbl[keep, , drop = FALSE])
  }
  
  tbl[0, , drop = FALSE]
}

.comparison_involves_bcell <- function(comparison_name, bcell_regex) {
  parts <- strsplit(as.character(comparison_name), "_from_", fixed = TRUE)
  labels <- unique(unlist(parts, use.names = FALSE))
  any(.matches_bcell_label(labels, bcell_regex))
}

.build_proximity_heatmap_plot <- function(enrichment_table, sample_id) {
  parsed <- .parse_unified_interactions(enrichment_table$unified_int)
  metric <- enrichment_table$enrichm
  labels <- sort(unique(c(parsed$source, parsed$target)))
  mat <- matrix(
    0,
    nrow = length(labels),
    ncol = length(labels),
    dimnames = list(labels, labels)
  )
  
  for (i in seq_len(nrow(parsed))) {
    src <- parsed$source[i]
    tgt <- parsed$target[i]
    if (!is.na(src) && !is.na(tgt)) {
      mat[src, tgt] <- metric[i]
      mat[tgt, src] <- metric[i]
    }
  }
  
  if (nrow(mat) > 2) {
    ord <- stats::hclust(stats::dist(mat))$order
    mat <- mat[ord, ord, drop = FALSE]
  }
  
  heatmap_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(heatmap_df) <- c("source", "target", "enrichment")
  heatmap_df$source <- factor(
    pretty_plot_label(heatmap_df$source, width = 18),
    levels = pretty_plot_label(rownames(mat), width = 18)
  )
  heatmap_df$target <- factor(
    pretty_plot_label(heatmap_df$target, width = 18),
    levels = pretty_plot_label(colnames(mat), width = 18)
  )
  
  ggplot2::ggplot(
    heatmap_df,
    ggplot2::aes(x = source, y = target, fill = enrichment)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(
      low = "#3B82F6",
      mid = "white",
      high = "#EF4444",
      midpoint = 0,
      name = "Spatial\nEnrichment"
    ) +
    ggplot2::labs(
      title = sample_plot_title(sample_id, "Spatial Proximity Enrichment Heatmap"),
      subtitle = "Red indicates enriched neighborhoods; blue indicates depleted neighborhoods",
      x = "Interacting Cell Type",
      y = "Interacting Cell Type"
    ) +
    presentation_theme(base_size = 12, x_angle = 50) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 9,
        angle = 50,
        hjust = 1,
        vjust = 1,
        lineheight = 0.9
      ),
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.9),
      legend.key.height = grid::unit(0.7, "cm")
    )
}

.build_proximity_network_plot <- function(enrichment_table, sample_id, max_edges = 70) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required to build the custom proximity network plot.")
  }
  
  parsed <- .parse_unified_interactions(enrichment_table$unified_int)
  plot_df <- cbind(enrichment_table, parsed, stringsAsFactors = FALSE)
  plot_df <- plot_df[!is.na(plot_df$source) & !is.na(plot_df$target), , drop = FALSE]
  plot_df$edge_direction <- ifelse(plot_df$enrichm >= 0, "Enriched proximity", "Depleted proximity")
  plot_df$edge_weight <- abs(plot_df$PI_value)
  plot_df$edge_significance <- pmin(plot_df$p.adj_higher, plot_df$p.adj_lower, na.rm = TRUE)
  
  sig_edges <- plot_df[plot_df$edge_significance < 0.05, , drop = FALSE]
  if (nrow(sig_edges) < max_edges) {
    top_edges <- plot_df[order(plot_df$edge_weight, decreasing = TRUE), , drop = FALSE]
    top_edges <- top_edges[!duplicated(top_edges$unified_int), , drop = FALSE]
    keep_n <- min(max_edges, nrow(top_edges))
    plot_df <- top_edges[seq_len(keep_n), , drop = FALSE]
  } else {
    plot_df <- sig_edges[order(sig_edges$edge_weight, decreasing = TRUE), , drop = FALSE]
  }
  
  graph_obj <- igraph::graph_from_data_frame(
    d = plot_df[, c("source", "target")],
    directed = FALSE
  )
  set.seed(42)
  layout_mat <- igraph::layout_with_kk(graph_obj)
  node_names <- igraph::V(graph_obj)$name
  node_df <- data.frame(
    node = node_names,
    x = layout_mat[, 1],
    y = layout_mat[, 2],
    label = pretty_plot_label(node_names, width = 18),
    stringsAsFactors = FALSE
  )
  
  # Build edge data frame with explicit source/target coordinate columns.
  # Two separate lookups are used (instead of two merge() calls) because
  # back-to-back merge() on the same coordinate columns (x, y) causes R to
  # append .x/.y suffixes on the second pass, leaving no plain "x"/"y" columns
  # for ggplot to find — the original source of the "object 'x' not found" crash.
  edge_df <- data.frame(
    source = plot_df$source,
    target = plot_df$target,
    edge_direction = plot_df$edge_direction,
    edge_weight = plot_df$edge_weight,
    stringsAsFactors = FALSE
  )
  # Explicit index-matched lookup avoids merge() suffix collisions entirely.
  src_idx <- match(edge_df$source, node_df$node)
  tgt_idx <- match(edge_df$target, node_df$node)
  edge_df$x    <- node_df$x[src_idx]
  edge_df$y    <- node_df$y[src_idx]
  edge_df$xend <- node_df$x[tgt_idx]
  edge_df$yend <- node_df$y[tgt_idx]
  edge_df$linewidth <- scales::rescale(edge_df$edge_weight, to = c(0.3, 1.8))
  
  palette <- setNames(grDevices::hcl.colors(length(node_names), palette = "Dynamic"), node_names)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        color = edge_direction,
        linewidth = linewidth
      ),
      alpha = 0.55,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Enriched proximity" = "#EF4444",
        "Depleted proximity" = "#3B82F6"
      ),
      name = "Interaction type"
    ) +
    ggplot2::scale_linewidth_identity() +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(x = x, y = y, fill = node),
      shape = 21,
      size = 4.5,
      color = "white",
      stroke = 0.5,
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      labels = pretty_plot_label(names(palette), width = 20),
      name = "Cell Type"
    ) +
    {
      if (requireNamespace("ggrepel", quietly = TRUE)) {
        ggrepel::geom_text_repel(
          data = node_df,
          ggplot2::aes(x = x, y = y, label = label),
          size = 3.6,
          box.padding = 0.35,
          point.padding = 0.15,
          segment.alpha = 0.35,
          max.overlaps = Inf
        )
      } else {
        ggplot2::geom_text(
          data = node_df,
          ggplot2::aes(x = x, y = y, label = label),
          size = 3.3,
          check_overlap = TRUE,
          vjust = -0.8
        )
      }
    } +
    ggplot2::labs(
      title = sample_plot_title(sample_id, "Spatial Proximity Network"),
      subtitle = "Top neighborhood enrichments shown; edge color marks enrichment vs depletion",
      x = NULL,
      y = NULL
    ) +
    presentation_theme(base_size = 12, legend_position = "right") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.box = "vertical"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(ncol = 2, override.aes = list(size = 4.5)),
      color = ggplot2::guide_legend(override.aes = list(linewidth = 1.2))
    )
  
  p
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

# Polygon coordinate extractor — returns a data.frame with columns
# geom, part, x, y, hole, cell_ID (one row per polygon vertex), or NULL
# when no polygon data are available on the Giotto object.
# NOTE: a sibling copy lives in 07_Annotation.R — keep the two in sync.
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


.giotto_get_expression <- function(gobj, values = "normalized",
                                    output = "matrix") {
  fn <- NULL
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("getExpression", envir = asNamespace("Giotto"), inherits = FALSE)) {
    fn <- get("getExpression", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getExpression", envir = asNamespace("GiottoClass"),
                    inherits = FALSE)) {
    fn <- get("getExpression", envir = asNamespace("GiottoClass"))
  } else {
    fn <- get("getExpression", mode = "function")
  }
  fn(gobject = gobj, values = values, output = output)
}


.giotto_get_spatial_locations <- function(gobj, output = "data.table") {
  fn <- NULL
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("getSpatialLocations", envir = asNamespace("Giotto"),
             inherits = FALSE)) {
    fn <- get("getSpatialLocations", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("getSpatialLocations", envir = asNamespace("GiottoClass"),
                    inherits = FALSE)) {
    fn <- get("getSpatialLocations", envir = asNamespace("GiottoClass"))
  } else {
    fn <- get("getSpatialLocations", mode = "function")
  }
  fn(gobject = gobj, output = output)
}


# Replace underscores with spaces for human-readable plot labels.
.pretty <- function(x) gsub("_", " ", as.character(x), fixed = TRUE)


# Per-gene spatial expression plots with the focus population (B cells by
# default) outlined on top. Two outputs per surviving marker gene:
#   * all-FOV image        -> results_dir/B_annotated_gene_expression/<sample>_<GENE>.png
#   * one image per FOV    -> results_dir/B_annotated_gene_expression/<GENE>/<sample>_<GENE>_FOV_<n>.png
#     (only for FOVs that actually contain at least one B cell)
# Skips gracefully if polygons, the annotation column, or the marker
# lists are missing.
plot_bcell_spatial_and_markers <- function(gobj,
                                           sample_id,
                                           results_dir,
                                           annotation_column,
                                           bcell_regex,
                                           bcell_markers     = character(),
                                           subtype_markers   = character(),
                                           highlight_label   = "B cells",
                                           highlight_colour  = "mediumblue") {
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (!annotation_column %in% names(meta)) {
    message("  \u26A0 annotation column '", annotation_column,
            "' missing; skipping B-cell spatial plots")
    return(invisible(NULL))
  }
  is_bcell <- grepl(bcell_regex, as.character(meta[[annotation_column]]),
                    ignore.case = TRUE)
  if (!any(is_bcell)) {
    message("  \u26A0 No cells match bcell_regex; skipping B-cell spatial plots")
    return(invisible(NULL))
  }

  candidates <- unique(c(bcell_markers, subtype_markers))
  candidates <- candidates[nzchar(candidates)]
  if (!length(candidates)) return(invisible(NULL))

  expr <- tryCatch(.giotto_get_expression(gobj, values = "normalized",
                                          output = "matrix"),
                   error = function(e) NULL)
  if (is.null(expr)) {
    message("  \u26A0 Expression matrix unavailable; per-gene plots skipped")
    return(invisible(NULL))
  }
  genes   <- intersect(candidates, rownames(expr))
  missing <- setdiff(candidates, rownames(expr))
  if (length(missing)) {
    message("  \u2139 Marker genes not on panel: ",
            paste(missing, collapse = ", "))
  }
  if (!length(genes)) return(invisible(NULL))

  poly_df <- .extract_polygon_df(gobj)
  if (is.null(poly_df)) {
    message("  \u26A0 No cell polygons on object; B-cell spatial plots skipped")
    return(invisible(NULL))
  }
  poly_df$poly_group <- paste(poly_df$cell_ID, poly_df$geom, poly_df$part,
                              sep = "_")
  bcell_ids <- meta$cell_ID[is_bcell]
  poly_df$.is_bcell <- poly_df$cell_ID %in% bcell_ids
  if ("fov" %in% names(meta)) {
    poly_df$fov <- unname(stats::setNames(meta$fov, meta$cell_ID)[poly_df$cell_ID])
  } else {
    poly_df$fov <- NA_integer_
  }
  fov_with_bcells <- sort(unique(meta$fov[is_bcell & !is.na(meta$fov)]))

  viz_root <- ensure_dir(file.path(results_dir, "B_annotated_gene_expression"))

  make_plot <- function(gene, df, fov_tag = NULL) {
    df$expr <- unname(expr[gene, ][df$cell_ID])
    bcell_df <- df[df$.is_bcell, , drop = FALSE]
    title_txt <- paste0(sample_id, " - ", gene, " expression",
                        if (!is.null(fov_tag)) paste0(" (FOV ", fov_tag, ")") else "")
    subtitle_txt <- paste0(.pretty(highlight_label), " highlighted from ",
                           .pretty(annotation_column), " annotation")
    ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = df,
        mapping = ggplot2::aes(x = x, y = y, group = poly_group, fill = expr),
        colour = "grey30", linewidth = 0.08
      ) +
      ggplot2::geom_polygon(
        data = bcell_df,
        mapping = ggplot2::aes(x = x, y = y, group = poly_group),
        fill = NA, colour = highlight_colour, linewidth = 0.15
      ) +
      ggplot2::scale_fill_gradient(low = "lightgrey", high = "red",
                                    name = .pretty(gene)) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = title_txt, subtitle = subtitle_txt,
                    x = NULL, y = NULL) +
      presentation_theme(base_size = 11, legend_position = "right") +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }

  for (g in genes) {
    save_presentation_plot(
      plot     = make_plot(g, poly_df),
      filename = file.path(viz_root, paste0(sample_id, "_", g, ".png")),
      width = 14, height = 10, dpi = 600
    )
    if (length(fov_with_bcells)) {
      gene_dir <- ensure_dir(file.path(viz_root, g))
      for (fv in fov_with_bcells) {
        df_fv <- poly_df[!is.na(poly_df$fov) & poly_df$fov == fv, , drop = FALSE]
        if (!nrow(df_fv)) next
        save_presentation_plot(
          plot     = make_plot(g, df_fv, fov_tag = fv),
          filename = file.path(gene_dir,
                               paste0(sample_id, "_", g, "_FOV_", fv, ".png")),
          width = 10, height = 10, dpi = 600
        )
      }
    }
  }
  cat("  \u2713 B-cell gene overlay plots saved (", length(genes),
      " gene(s), ", length(fov_with_bcells), " B-cell FOV(s))\n", sep = "")
  invisible(genes)
}


# Per-B-cell-cluster neighbourhood polygon plots coloured by distance to
# nearest B cell. B-cells are first grouped by single-linkage spatial
# clustering (h = bcell_cluster_eps); k-NN context is then found around
# each cluster's centroid with k scaled by cluster size so dense B-cell
# patches get a proportionally larger window.
# Output: <results_dir>/neighbours/cluster_<NNN>_<n>cells/knn_k<K>.png
plot_bcell_neighbourhoods <- function(gobj,
                                      sample_id,
                                      results_dir,
                                      annotation_column,
                                      bcell_regex,
                                      highlight_label    = "B cells",
                                      highlight_colour   = "mediumblue",
                                      bcell_cluster_eps  = 80,
                                      k_base             = 20,
                                      k_increment        = 10,
                                      max_clusters_plot  = 40) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    message("  \u26A0 FNN not installed; B-cell neighbourhood plots skipped")
    return(invisible(NULL))
  }
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (!annotation_column %in% names(meta)) {
    message("  \u26A0 annotation column '", annotation_column,
            "' missing; neighbourhood plots skipped")
    return(invisible(NULL))
  }
  is_bcell <- grepl(bcell_regex, as.character(meta[[annotation_column]]),
                    ignore.case = TRUE)
  if (!any(is_bcell)) return(invisible(NULL))

  sl <- tryCatch(as.data.frame(.giotto_get_spatial_locations(gobj,
                                                               output = "data.table")),
                 error = function(e) NULL)
  if (is.null(sl) || !all(c("sdimx", "sdimy", "cell_ID") %in% names(sl))) {
    message("  \u26A0 spatial locations unavailable; neighbourhood plots skipped")
    return(invisible(NULL))
  }
  cells <- merge(meta[, c("cell_ID", annotation_column), drop = FALSE],
                 sl[, c("cell_ID", "sdimx", "sdimy")],
                 by = "cell_ID", sort = FALSE)
  cells$is_bcell <- grepl(bcell_regex, as.character(cells[[annotation_column]]),
                          ignore.case = TRUE)
  bcell_ids <- cells$cell_ID[cells$is_bcell]

  poly_df <- .extract_polygon_df(gobj)
  if (is.null(poly_df)) {
    message("  \u26A0 No cell polygons; neighbourhood plots skipped")
    return(invisible(NULL))
  }
  poly_df$poly_group <- paste(poly_df$cell_ID, poly_df$geom, poly_df$part,
                              sep = "_")
  poly_df$is_bcell <- poly_df$cell_ID %in% bcell_ids

  bcells_df <- cells[cells$is_bcell, c("cell_ID", "sdimx", "sdimy"), drop = FALSE]
  if (nrow(bcells_df) == 1) {
    bcells_df$cluster_id <- 1L
  } else {
    d <- stats::dist(bcells_df[, c("sdimx", "sdimy")])
    bcells_df$cluster_id <- stats::cutree(
      stats::hclust(d, method = "single"), h = bcell_cluster_eps
    )
  }
  clust_summary <- stats::aggregate(cbind(sdimx, sdimy) ~ cluster_id,
                                     data = bcells_df, FUN = mean)
  clust_summary$n_bcells <- as.integer(
    table(bcells_df$cluster_id)[as.character(clust_summary$cluster_id)]
  )
  clust_summary <- clust_summary[order(-clust_summary$n_bcells), , drop = FALSE]
  clust_summary <- head(clust_summary, max_clusters_plot)

  xy_all <- as.matrix(cells[, c("sdimx", "sdimy")])
  rownames(xy_all) <- cells$cell_ID

  nbr_dir <- ensure_dir(file.path(results_dir, "neighbours"))
  subtitle_common <- paste0(.pretty(highlight_label), " highlighted from ",
                            .pretty(annotation_column), " annotation")

  for (i in seq_len(nrow(clust_summary))) {
    row <- clust_summary[i, , drop = FALSE]
    member_ids <- bcells_df$cell_ID[bcells_df$cluster_id == row$cluster_id]
    k <- k_base + k_increment * (row$n_bcells - 1) + row$n_bcells
    centroid <- c(row$sdimx, row$sdimy)
    nn <- FNN::get.knnx(data = xy_all, query = matrix(centroid, nrow = 1), k = k)
    nbr_ids <- unique(c(member_ids, rownames(xy_all)[nn$nn.index[1, ]]))

    xy_nbrs <- xy_all[nbr_ids, , drop = FALSE]
    xy_members <- xy_all[member_ids, , drop = FALSE]
    dist_nn <- FNN::get.knnx(data = xy_members, query = xy_nbrs, k = 1)
    dist_map <- stats::setNames(dist_nn$nn.dist[, 1], nbr_ids)

    df <- poly_df[poly_df$cell_ID %in% nbr_ids, , drop = FALSE]
    df$distance <- unname(dist_map[df$cell_ID])

    clust_dir <- ensure_dir(file.path(
      nbr_dir, sprintf("cluster_%03d_%dcells", i, row$n_bcells)
    ))
    out_png <- file.path(clust_dir, paste0("knn_k", k, ".png"))

    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = df,
        mapping = ggplot2::aes(x = x, y = y, group = poly_group, fill = distance),
        colour = "grey30", linewidth = 0.08
      ) +
      ggplot2::geom_polygon(
        data = df[df$is_bcell, , drop = FALSE],
        mapping = ggplot2::aes(x = x, y = y, group = poly_group),
        fill = NA, colour = highlight_colour, linewidth = 0.15
      ) +
      viridis::scale_fill_viridis(
        option = "magma", direction = -1,
        name = "Distance to\nnearest B cell"
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        title = paste0(sample_id, " - B cell cluster ", i, " neighbourhood"),
        subtitle = paste0(row$n_bcells, " B cell(s), k = ", k, " \u2014 ",
                          subtitle_common),
        x = NULL, y = NULL
      ) +
      presentation_theme(base_size = 11, legend_position = "right") +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
    save_presentation_plot(p, filename = out_png,
                           width = 10, height = 10, dpi = 600)
  }
  cat("  \u2713 B-cell neighbourhood plots saved (", nrow(clust_summary),
      " cluster(s))\n", sep = "")
  invisible(nrow(clust_summary))
}


# Spatial-niche + distance-to-B-cell plots. Niches are defined by k-NN
# cell-type composition followed by k-means (same recipe as
# 12_Spatial_Differential_Expression.R::assign_spatial_niches).
# Output: <results_dir>/niches/<sample>_niches_annotated.png
#         <results_dir>/niches/<sample>_distance_to_bcells.png
#         <results_dir>/niches/per_fov/<sample>_..._FOV_<n>.png (FOVs with B cells)
plot_bcell_niches <- function(gobj,
                              sample_id,
                              results_dir,
                              annotation_column,
                              bcell_regex,
                              highlight_label  = "B cells",
                              highlight_colour = "mediumblue",
                              n_niches         = 6,
                              niche_knn_k      = 30,
                              um_per_px        = 0.12028) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    message("  \u26A0 FNN not installed; niche plots skipped")
    return(invisible(NULL))
  }
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (!annotation_column %in% names(meta)) {
    message("  \u26A0 annotation column missing; niche plots skipped")
    return(invisible(NULL))
  }
  is_bcell <- grepl(bcell_regex, as.character(meta[[annotation_column]]),
                    ignore.case = TRUE)
  if (!any(is_bcell)) return(invisible(NULL))

  sl <- tryCatch(as.data.frame(.giotto_get_spatial_locations(gobj,
                                                               output = "data.table")),
                 error = function(e) NULL)
  if (is.null(sl) || !all(c("sdimx", "sdimy", "cell_ID") %in% names(sl))) {
    message("  \u26A0 spatial locations unavailable; niche plots skipped")
    return(invisible(NULL))
  }
  meta_cols <- c("cell_ID", annotation_column,
                 intersect("fov", names(meta)))
  cells <- merge(meta[, meta_cols, drop = FALSE],
                 sl[, c("cell_ID", "sdimx", "sdimy")],
                 by = "cell_ID", sort = FALSE)
  cells$cell_type <- as.character(cells[[annotation_column]])
  cells$is_bcell <- grepl(bcell_regex, cells$cell_type, ignore.case = TRUE)
  bcell_ids <- cells$cell_ID[cells$is_bcell]

  poly_df <- .extract_polygon_df(gobj)
  if (is.null(poly_df)) {
    message("  \u26A0 No cell polygons; niche plots skipped")
    return(invisible(NULL))
  }
  poly_df$poly_group <- paste(poly_df$cell_ID, poly_df$geom, poly_df$part,
                              sep = "_")
  poly_df$is_bcell <- poly_df$cell_ID %in% bcell_ids
  if ("fov" %in% names(meta)) {
    poly_df$fov <- unname(stats::setNames(meta$fov, meta$cell_ID)[poly_df$cell_ID])
  } else {
    poly_df$fov <- NA_integer_
  }

  xy_all <- as.matrix(cells[, c("sdimx", "sdimy")])
  rownames(xy_all) <- cells$cell_ID
  ct_levels <- sort(unique(cells$cell_type))
  nn_all <- FNN::get.knn(xy_all, k = niche_knn_k)
  ct_fac <- factor(cells$cell_type, levels = ct_levels)
  nbr_ct <- matrix(as.integer(ct_fac)[nn_all$nn.index], ncol = niche_knn_k)
  niche_counts <- matrix(
    0L, nrow = nrow(cells), ncol = length(ct_levels),
    dimnames = list(cells$cell_ID, ct_levels)
  )
  for (j in seq_len(niche_knn_k)) {
    niche_counts[cbind(seq_len(nrow(cells)), nbr_ct[, j])] <-
      niche_counts[cbind(seq_len(nrow(cells)), nbr_ct[, j])] + 1L
  }
  niche_props <- niche_counts / rowSums(niche_counts)

  n_unique <- nrow(unique(niche_props))
  set.seed(42)
  km <- stats::kmeans(niche_props,
                      centers = min(n_niches, n_unique),
                      nstart = 25, iter.max = 50)
  cells$niche <- paste0("niche_", km$cluster)
  poly_df$niche <- cells$niche[match(poly_df$cell_ID, cells$cell_ID)]

  niche_palette <- stats::setNames(
    grDevices::hcl.colors(length(unique(cells$niche)), palette = "Set 2"),
    sort(unique(cells$niche))
  )

  niche_dir <- ensure_dir(file.path(results_dir, "niches"))
  subtitle_common <- paste0(.pretty(highlight_label), " highlighted from ",
                            .pretty(annotation_column), " annotation")

  make_niche_plot <- function(df_poly, title_suffix = "") {
    ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data = df_poly,
        mapping = ggplot2::aes(x = x, y = y, group = poly_group, fill = niche),
        colour = "grey30", linewidth = 0.08
      ) +
      ggplot2::geom_polygon(
        data = df_poly[df_poly$is_bcell, , drop = FALSE],
        mapping = ggplot2::aes(x = x, y = y, group = poly_group),
        fill = NA, colour = highlight_colour, linewidth = 0.15
      ) +
      ggplot2::scale_fill_manual(values = niche_palette, name = "Niche",
                                  labels = .pretty) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        title = paste0(sample_id, " - spatial niches", title_suffix),
        subtitle = subtitle_common,
        x = NULL, y = NULL
      ) +
      presentation_theme(base_size = 11, legend_position = "right") +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }

  bcell_xy <- xy_all[bcell_ids, , drop = FALSE]
  dist_nn  <- FNN::get.knnx(data = bcell_xy, query = xy_all, k = 1)
  cells$distance_to_bcell_um <- dist_nn$nn.dist[, 1] * um_per_px
  cap <- stats::quantile(cells$distance_to_bcell_um, 0.9, na.rm = TRUE)
  cells$distance_capped <- pmin(cells$distance_to_bcell_um, cap)

  make_distance_plot <- function(df_cells, title_suffix = "") {
    ggplot2::ggplot(
        df_cells,
        ggplot2::aes(x = sdimx, y = sdimy, colour = distance_capped)) +
      ggplot2::geom_point(size = 0.25, alpha = 0.7) +
      viridis::scale_color_viridis(
        name = paste0("Distance to\nnearest B cell\n(\u00b5m, capped at\n",
                      round(cap), ")")
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        title = paste0(sample_id, " - distance to nearest B cell",
                       title_suffix),
        subtitle = paste0(.pretty(highlight_label), " defined by ",
                          .pretty(annotation_column),
                          ", cap = 90th percentile"),
        x = NULL, y = NULL
      ) +
      presentation_theme(base_size = 11, legend_position = "right") +
      ggplot2::theme(
        axis.title = ggplot2::element_blank(),
        axis.text  = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank()
      )
  }

  save_presentation_plot(
    make_niche_plot(poly_df),
    filename = file.path(niche_dir,
                          paste0(sample_id, "_niches_annotated.png")),
    width = 14, height = 10, dpi = 600
  )
  save_presentation_plot(
    make_distance_plot(cells),
    filename = file.path(niche_dir,
                          paste0(sample_id, "_distance_to_bcells.png")),
    width = 14, height = 10, dpi = 600
  )

  # Per-FOV niche + distance plots (only FOVs that contain B cells)
  n_fov_written <- 0L
  if ("fov" %in% names(cells)) {
    per_fov_dir <- ensure_dir(file.path(niche_dir, "per_fov"))
    fov_with_bcells <- sort(unique(cells$fov[cells$is_bcell & !is.na(cells$fov)]))
    for (fv in fov_with_bcells) {
      df_fv_poly <- poly_df[!is.na(poly_df$fov) & poly_df$fov == fv, , drop = FALSE]
      df_fv_cells <- cells[!is.na(cells$fov) & cells$fov == fv, , drop = FALSE]
      if (!nrow(df_fv_poly) || !nrow(df_fv_cells)) next

      save_presentation_plot(
        make_niche_plot(df_fv_poly, title_suffix = paste0(" (FOV ", fv, ")")),
        filename = file.path(per_fov_dir,
                              paste0(sample_id, "_niches_annotated_FOV_",
                                      fv, ".png")),
        width = 10, height = 10, dpi = 600
      )
      save_presentation_plot(
        make_distance_plot(df_fv_cells, title_suffix = paste0(" (FOV ", fv, ")")),
        filename = file.path(per_fov_dir,
                              paste0(sample_id, "_distance_to_bcells_FOV_",
                                      fv, ".png")),
        width = 10, height = 10, dpi = 600
      )
      n_fov_written <- n_fov_written + 1L
    }
  }
  cat("  \u2713 B-cell niche plots saved (", n_fov_written,
      " per-FOV pairs)\n", sep = "")
  invisible(cells$niche)
}


# -----------------------------------------------------------------------------
# B-cell subclustering (reuses scripts 04 + 05)
# -----------------------------------------------------------------------------

plot_bcell_subtype_umap <- function(gobj_bcell,
                                    sample_id,
                                    subtype_markers,
                                    out_dir,
                                    cluster_column = "leiden_bcell_subcluster") {
  umap_df <- .prepare_dim_plot_data(gobj_bcell, "umap", cluster_column)
  umap_df[[cluster_column]] <- factor(umap_df[[cluster_column]])

  leiden_plot <- ggplot2::ggplot(
      umap_df,
      ggplot2::aes(x = dim1, y = dim2, colour = .data[[cluster_column]])
    ) +
    ggplot2::geom_point(size = 0.6, alpha = 0.8) +
    ggplot2::labs(
      title  = sample_plot_title(sample_id, "B-cell subclusters"),
      x      = "UMAP 1",
      y      = "UMAP 2",
      colour = "Subcluster"
    ) +
    presentation_theme(base_size = 12)

  save_presentation_plot(
    plot     = leiden_plot,
    filename = file.path(out_dir, paste0(sample_id, "_bcell_umap_leiden.png")),
    width    = 12,
    height   = 10,
    dpi      = 300
  )

  expr <- tryCatch(
    .giotto_get_expression(gobj_bcell, values = "normalized", output = "matrix"),
    error = function(e) NULL
  )
  if (is.null(expr) || !length(subtype_markers)) {
    return(invisible(NULL))
  }
  genes <- intersect(subtype_markers, rownames(expr))
  missing <- setdiff(subtype_markers, rownames(expr))
  if (length(missing)) {
    message("  \u2139 Subtype markers not on panel: ",
            paste(missing, collapse = ", "))
  }
  if (!length(genes)) return(invisible(NULL))

  dim_base <- umap_df[, c("cell_ID", "dim1", "dim2")]
  gene_plots <- lapply(genes, function(g) {
    df <- dim_base
    df$expr <- unname(expr[g, df$cell_ID])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = dim1, y = dim2, colour = expr)) +
      ggplot2::geom_point(size = 0.6, alpha = 0.85) +
      ggplot2::scale_colour_gradient(
        low = "lightgrey", high = "red", name = .pretty(g)
      ) +
      ggplot2::labs(
        title = sample_plot_title(
          sample_id, paste0(g, " (B-cell subclusters)")
        ),
        x = "UMAP 1", y = "UMAP 2"
      ) +
      presentation_theme(base_size = 11)
    save_presentation_plot(
      plot     = p,
      filename = file.path(out_dir,
                           paste0(sample_id, "_bcell_umap_", g, ".png")),
      width    = 10,
      height   = 8,
      dpi      = 300
    )
    p
  })

  if (requireNamespace("patchwork", quietly = TRUE) && length(gene_plots) > 1) {
    grid  <- optimal_grid_dims(length(gene_plots))
    panel <- patchwork::wrap_plots(gene_plots,
                                   ncol = grid$ncol, nrow = grid$nrow)
    save_presentation_plot(
      plot     = panel,
      filename = file.path(
        out_dir,
        paste0(sample_id, "_bcell_umap_subtype_markers_panel.png")
      ),
      width    = grid$ncol * 5,
      height   = grid$nrow * 4,
      dpi      = 300
    )
  }
  invisible(genes)
}


run_bcell_subclustering <- function(gobj,
                                    sample_id,
                                    results_dir,
                                    annotation_column,
                                    bcell_regex,
                                    subtype_markers     = character(),
                                    min_cells           = 50,
                                    fallback_min_cells  = 20,
                                    n_hvgs              = 250,
                                    n_pcs               = 20,
                                    umap_n_neighbors    = 15,
                                    umap_min_dist       = 0.3,
                                    k_nn                = 10,
                                    leiden_resolution   = 0.4,
                                    leiden_n_iterations = 200,
                                    resolution_sweep    = NULL,
                                    scripts_dir         = NULL,
                                    python_path         = NULL,
                                    save_object         = TRUE,
                                    seed                = 42) {
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (!annotation_column %in% names(meta)) {
    message("B-cell subclustering skipped: annotation column '",
            annotation_column, "' not present.")
    return(invisible(NULL))
  }
  is_bcell <- grepl(bcell_regex, as.character(meta[[annotation_column]]),
                    ignore.case = TRUE)
  bcell_ids <- meta$cell_ID[is_bcell]
  n_bcells <- length(bcell_ids)
  if (n_bcells < fallback_min_cells) {
    message("B-cell subclustering skipped: only ",
            n_bcells, " B-annotated cells (< fallback_min_cells = ",
            fallback_min_cells, ").")
    return(invisible(NULL))
  }
  if (n_bcells < min_cells) {
    warning("B-cell subclustering proceeding with only ", n_bcells,
            " B-annotated cells (< min_cells = ", min_cells,
            "). Subclustering quality may be unreliable; consider raising ",
            "fallback_min_cells or inspecting the results carefully.",
            call. = FALSE)
  }

  subset_fn <- tryCatch(
    get("subsetGiotto", envir = asNamespace("Giotto")),
    error = function(e) get("subsetGiotto", envir = asNamespace("GiottoClass"))
  )
  gobj_bcell <- subset_fn(gobject = gobj, cell_ids = bcell_ids)

  subcluster_dir <- ensure_dir(file.path(results_dir, "subcluster"))
  sid <- paste0(sample_id, "_bcell")

  gobj_bcell <- dimensionality_reduction(
    gobj             = gobj_bcell,
    sample_id        = sid,
    output_dir       = subcluster_dir,
    n_hvgs           = n_hvgs,
    n_pcs            = n_pcs,
    umap_n_neighbors = umap_n_neighbors,
    umap_min_dist    = umap_min_dist,
    spatial_hvg      = FALSE,
    seed             = seed
  )

  gobj_bcell <- perform_clustering(
    gobj                = gobj_bcell,
    sample_id           = sid,
    output_dir          = subcluster_dir,
    k_nn                = k_nn,
    resolution          = leiden_resolution,
    dimensions_to_use   = seq_len(n_pcs),
    scripts_dir         = scripts_dir,
    python_path         = python_path,
    inspect_snn         = FALSE,
    leiden_n_iterations = leiden_n_iterations,
    resolution_sweep    = resolution_sweep,
    seed                = seed
  )

  meta_bcell <- as.data.frame(.giotto_pdata_dt(gobj_bcell))
  if ("leiden_clust" %in% names(meta_bcell)) {
    gobj_bcell <- addCellMetadata(
      gobject        = gobj_bcell,
      new_metadata   = data.frame(
        cell_ID                 = meta_bcell$cell_ID,
        leiden_bcell_subcluster = meta_bcell$leiden_clust,
        stringsAsFactors        = FALSE
      ),
      by_column      = TRUE,
      column_cell_ID = "cell_ID"
    )
  }

  tryCatch({
    plot_bcell_subtype_umap(
      gobj_bcell      = gobj_bcell,
      sample_id       = sample_id,
      subtype_markers = subtype_markers,
      out_dir         = subcluster_dir,
      cluster_column  = "leiden_bcell_subcluster"
    )
  }, error = function(e) {
    message("B-cell subtype UMAP plots skipped: ", conditionMessage(e))
  })

  if (save_object) {
    save_giotto_checkpoint(
      gobj           = gobj_bcell,
      checkpoint_dir = file.path(dirname(results_dir),
                                 "Giotto_Object_BCell_Subcluster"),
      metadata       = list(
        stage             = "bcell_subcluster",
        annotation_column = annotation_column,
        bcell_regex       = bcell_regex,
        n_bcells          = length(bcell_ids)
      )
    )
  }

  invisible(gobj_bcell)
}


run_bcell_microenvironment_analysis <- function(gobj,
                                                sample_id,
                                                output_dir,
                                                annotation_column = NULL,
                                                bcell_regex = "^B\\.cell$",
                                                spatial_network_name = "Delaunay_network",
                                                number_of_simulations = 250,
                                                max_network_edges = 70,
                                                bcell_markers = character(),
                                                subtype_markers = character(),
                                                bcell_subcluster_enabled            = TRUE,
                                                bcell_subcluster_min_cells          = 50,
                                                bcell_subcluster_fallback_min_cells = 20,
                                                bcell_subcluster_n_hvgs             = 250,
                                                bcell_subcluster_n_pcs              = 20,
                                                bcell_subcluster_umap_n_neighbors   = 15,
                                                bcell_subcluster_umap_min_dist      = 0.3,
                                                bcell_subcluster_k_nn               = 10,
                                                bcell_subcluster_resolution         = 0.4,
                                                bcell_subcluster_resolution_sweep   = NULL,
                                                scripts_dir = NULL,
                                                python_path = NULL,
                                                seed = 42,
                                                save_object = FALSE) {
  cat("\n========================================\n")
  cat("STEP 11: Focused Cell-Type Microenvironment\n")
  cat("Sample:", sample_id, "\n")
  cat("Focus regex:", bcell_regex, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    manifest_path <- file.path(gobj, "manifest.json")
    if (file.exists(manifest_path)) {
      gobj <- load_giotto_checkpoint(gobj)
    } else {
      gobj <- .giotto_load(gobj)
    }
  }
  
  results_dir <- ensure_dir(file.path(output_dir, "11_BCell_Microenvironment"))
  metadata <- tibble::as_tibble(.giotto_pdata_dt(gobj))
  annotation_column <- detect_annotation_column(metadata, annotation_column)
  
  abundance_table <- metadata
  abundance_table$annotation <- metadata[[annotation_column]]
  abundance_table <- dplyr::summarise(
    dplyr::group_by(abundance_table, annotation),
    n_cells = dplyr::n(),
    .groups = "drop"
  )
  abundance_table <- dplyr::arrange(abundance_table, dplyr::desc(n_cells))
  abundance_table <- dplyr::mutate(
    abundance_table,
    sample_id = sample_id,
    is_focus_celltype = grepl(bcell_regex, annotation, ignore.case = TRUE)
  )
  
  readr::write_csv(
    abundance_table,
    file.path(results_dir, paste0(sample_id, "_celltype_abundance.csv"))
  )
  
  gobj <- ensure_spatial_network(gobj, spatial_network_name = spatial_network_name)
  
  cp_scores <- .giotto_cell_proximity_enrichment(
    gobj = gobj,
    spatial_network_name = spatial_network_name,
    cluster_column = annotation_column,
    number_of_simulations = number_of_simulations
  )
  
  readr::write_csv(
    tibble::as_tibble(cp_scores$raw_sim_table),
    file.path(results_dir, paste0(sample_id, "_cell_proximity_raw.csv"))
  )
  
  enrichment_table <- tibble::as_tibble(cp_scores$enrichm_res)
  bcell_table <- .filter_interactions_for_bcells(enrichment_table, bcell_regex = bcell_regex)
  
  readr::write_csv(
    enrichment_table,
    file.path(results_dir, paste0(sample_id, "_cell_proximity_enrichment.csv"))
  )
  readr::write_csv(
    bcell_table,
    file.path(results_dir, paste0(sample_id, "_bcell_interactions.csv"))
  )
  
  # Bring forward focused-celltype CCI outputs from step 10 when available.
  liana_path <- file.path(
    output_dir,
    "10_CCI_Analysis",
    "liana",
    paste0(sample_id, "_liana_aggregate.csv")
  )
  if (file.exists(liana_path)) {
    tryCatch({
      liana_tbl <- readr::read_csv(liana_path, show_col_types = FALSE)
      if (all(c("source", "target") %in% names(liana_tbl))) {
        bcell_liana <- .filter_interactions_for_bcells(liana_tbl, bcell_regex = bcell_regex)
        readr::write_csv(
          bcell_liana,
          file.path(results_dir, paste0(sample_id, "_bcell_liana_interactions.csv"))
        )
      }
    }, error = function(e) {
      message("Focused cell-type LIANA summary skipped: ", conditionMessage(e))
    })
  }
  
  nichenet_root <- file.path(output_dir, "10_CCI_Analysis", "nichenet")
  if (dir.exists(nichenet_root)) {
    tryCatch({
      comparison_dirs <- list.dirs(nichenet_root, full.names = TRUE, recursive = TRUE)
      comparison_dirs <- comparison_dirs[comparison_dirs != nichenet_root]
      comparison_dirs <- comparison_dirs[
        vapply(
          basename(comparison_dirs),
          .comparison_involves_bcell,
          logical(1),
          bcell_regex = bcell_regex
        )
      ]
      nichenet_tables <- lapply(comparison_dirs, function(dir_path) {
        csv_path <- file.path(dir_path, paste0(sample_id, "_ligand_activities.csv"))
        if (!file.exists(csv_path)) {
          return(NULL)
        }
        tbl <- readr::read_csv(csv_path, show_col_types = FALSE)
        tbl$comparison <- basename(dir_path)
        tbl
      })
      nichenet_tables <- Filter(Negate(is.null), nichenet_tables)
      if (length(nichenet_tables) > 0) {
        bcell_nichenet <- dplyr::bind_rows(nichenet_tables)
        readr::write_csv(
          bcell_nichenet,
          file.path(results_dir, paste0(sample_id, "_bcell_nichenet_ligand_activities.csv"))
        )
      }
    }, error = function(e) {
      message("Focused cell-type NicheNet summary skipped: ", conditionMessage(e))
    })
  }
  
  tryCatch({
    heatmap_plot <- .build_proximity_heatmap_plot(enrichment_table, sample_id = sample_id)
    save_presentation_plot(
      plot = heatmap_plot,
      filename = file.path(results_dir, paste0(sample_id, "_cell_proximity_heatmap.png")),
      width = 12,
      height = 10
    )
  }, error = function(e) {
    message("cellProximityHeatmap() skipped: ", conditionMessage(e))
  })
  
  tryCatch({
    network_plot <- .build_proximity_network_plot(
      enrichment_table,
      sample_id = sample_id,
      max_edges = max_network_edges
    )
    save_presentation_plot(
      plot = network_plot,
      filename = file.path(results_dir, paste0(sample_id, "_cell_proximity_network.png")),
      width = 13,
      height = 10
    )
  }, error = function(e) {
    message("cellProximityNetwork() skipped: ", conditionMessage(e))
    # Diagnostic: show available columns in enrichment_table to aid future debugging.
    message(
      "  [debug] enrichment_table columns: ",
      paste(names(enrichment_table), collapse = ", ")
    )
  })
  
  tryCatch({
    plot_bcell_spatial_and_markers(
      gobj              = gobj,
      sample_id         = sample_id,
      results_dir       = results_dir,
      annotation_column = annotation_column,
      bcell_regex       = bcell_regex,
      bcell_markers     = bcell_markers,
      subtype_markers   = subtype_markers
    )
  }, error = function(e) {
    message("B-cell spatial/marker plots skipped: ", conditionMessage(e))
  })

  tryCatch({
    plot_bcell_neighbourhoods(
      gobj              = gobj,
      sample_id         = sample_id,
      results_dir       = results_dir,
      annotation_column = annotation_column,
      bcell_regex       = bcell_regex
    )
  }, error = function(e) {
    message("B-cell neighbourhood plots skipped: ", conditionMessage(e))
  })

  tryCatch({
    plot_bcell_niches(
      gobj              = gobj,
      sample_id         = sample_id,
      results_dir       = results_dir,
      annotation_column = annotation_column,
      bcell_regex       = bcell_regex
    )
  }, error = function(e) {
    message("B-cell niche plots skipped: ", conditionMessage(e))
  })

  if (isTRUE(bcell_subcluster_enabled)) {
    tryCatch({
      run_bcell_subclustering(
        gobj                = gobj,
        sample_id           = sample_id,
        results_dir         = results_dir,
        annotation_column   = annotation_column,
        bcell_regex         = bcell_regex,
        subtype_markers     = subtype_markers,
        min_cells           = bcell_subcluster_min_cells,
        fallback_min_cells  = bcell_subcluster_fallback_min_cells,
        n_hvgs              = bcell_subcluster_n_hvgs,
        n_pcs               = bcell_subcluster_n_pcs,
        umap_n_neighbors    = bcell_subcluster_umap_n_neighbors,
        umap_min_dist       = bcell_subcluster_umap_min_dist,
        k_nn                = bcell_subcluster_k_nn,
        leiden_resolution   = bcell_subcluster_resolution,
        resolution_sweep    = bcell_subcluster_resolution_sweep,
        scripts_dir         = scripts_dir,
        python_path         = python_path,
        save_object         = save_object,
        seed                = seed
      )
    }, error = function(e) {
      message("B-cell subclustering skipped: ", conditionMessage(e))
    })
  }

  if (save_object) {
    save_giotto_checkpoint(
      gobj = gobj,
      checkpoint_dir = file.path(output_dir, "Giotto_Object_BCell_Analysis"),
      metadata = list(stage = "bcell_microenvironment", annotation_column = annotation_column)
    )
  }
  
  # FIX #1: Attach summary results as attributes so they travel with the
  # object, then return the (possibly modified) gobj for pipeline checkpoint
  # consistency.  Previously this returned a list, causing the pipeline to
  # lose any object modifications (e.g. spatial network creation).
  attr(gobj, "bcell_analysis") <- list(
    annotation_column = annotation_column,
    abundance_table = abundance_table,
    enrichment_table = enrichment_table,
    bcell_table = bcell_table
  )
  invisible(gobj)
}
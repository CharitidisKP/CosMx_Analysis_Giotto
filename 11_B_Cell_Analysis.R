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

# Giotto-native outlined polygon renderer.
# NOTE: a sibling copy lives in 07_Annotation.R — keep the two in sync.
.spat_in_situ_outlined <- function(gobj,
                                    fill_col,
                                    fill_as_factor,
                                    colour_map   = NULL,
                                    gradient     = c("lightgrey", "red"),
                                    legend_title = "Cell Type",
                                    title_txt) {
  fn <- NULL
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("spatInSituPlotPoints", envir = asNamespace("Giotto"),
             inherits = FALSE)) {
    fn <- get("spatInSituPlotPoints", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoVisuals", quietly = TRUE) &&
             exists("spatInSituPlotPoints", envir = asNamespace("GiottoVisuals"),
                    inherits = FALSE)) {
    fn <- get("spatInSituPlotPoints", envir = asNamespace("GiottoVisuals"))
  } else if (exists("spatInSituPlotPoints", mode = "function")) {
    fn <- get("spatInSituPlotPoints", mode = "function")
  }
  if (is.null(fn)) return(NULL)

  # Giotto's spatInSituPlotPoints() arg names drift between releases — probe
  # the function's formals and only pass arguments it actually accepts.
  accepted <- names(formals(fn))
  pick <- function(val, candidates) {
    hit <- intersect(candidates, accepted)
    if (!length(hit)) return(NULL)
    stats::setNames(list(val), hit[1])
  }

  args <- list(
    gobject                = gobj,
    show_polygon           = TRUE,
    polygon_feat_type      = "cell",
    polygon_fill           = fill_col,
    polygon_fill_as_factor = fill_as_factor,
    polygon_alpha          = 0.75,
    show_image             = FALSE,
    return_plot            = TRUE,
    save_plot              = FALSE
  )
  args <- c(args, pick("grey20",
    c("polygon_line_color", "polygon_color",
      "polygon_stroke_color", "polygon_border_color")))
  args <- c(args, pick(0.15,
    c("polygon_line_size", "polygon_stroke_size",
      "polygon_border_size", "polygon_size")))
  if (fill_as_factor && !is.null(colour_map)) {
    args <- c(args, pick(colour_map,
      c("polygon_fill_code", "polygon_fill_colors",
        "polygon_fill_values")))
  } else if (!fill_as_factor) {
    args <- c(args, pick(gradient,
      c("polygon_fill_gradient", "polygon_fill_gradient_colors",
        "polygon_gradient")))
  }
  args <- args[names(args) %in% accepted | names(args) == "gobject"]

  p <- tryCatch(do.call(fn, args), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  if (!inherits(p, "ggplot") && !is.null(p$ggobj)) p <- p$ggobj
  if (!inherits(p, "ggplot")) return(NULL)

  p +
    ggplot2::labs(title = title_txt, x = NULL, y = NULL) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(title = legend_title, ncol = 1,
                                    override.aes = list(size = 4)),
      color = ggplot2::guide_legend(title = legend_title, ncol = 1,
                                    override.aes = list(size = 4))
    ) +
    presentation_theme(base_size = 11, legend_position = "right") +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}


.giotto_add_cell_metadata <- function(gobj, new_metadata,
                                       by_column = TRUE,
                                       column_cell_ID = "cell_ID") {
  fn <- NULL
  if (requireNamespace("Giotto", quietly = TRUE) &&
      exists("addCellMetadata", envir = asNamespace("Giotto"), inherits = FALSE)) {
    fn <- get("addCellMetadata", envir = asNamespace("Giotto"))
  } else if (requireNamespace("GiottoClass", quietly = TRUE) &&
             exists("addCellMetadata", envir = asNamespace("GiottoClass"),
                    inherits = FALSE)) {
    fn <- get("addCellMetadata", envir = asNamespace("GiottoClass"))
  } else {
    fn <- get("addCellMetadata", mode = "function")
  }
  fn(gobject = gobj, new_metadata = new_metadata,
     by_column = by_column, column_cell_ID = column_cell_ID)
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


# Spatial polygon plots for the focus population (B cells by default):
#  (A) B-cell vs other highlight,
#  (B) per-gene expression polygon plots for marker genes present on the panel,
#  (C) combined multi-panel figure via patchwork.
# All outputs land in results_dir. Skips gracefully if polygons or the
# annotation column are missing.
plot_bcell_spatial_and_markers <- function(gobj,
                                           sample_id,
                                           results_dir,
                                           annotation_column,
                                           bcell_regex,
                                           bcell_markers    = character(),
                                           subtype_markers  = character()) {
  meta <- as.data.frame(.giotto_pdata_dt(gobj))
  if (!annotation_column %in% names(meta)) {
    message("  \u26A0 annotation column '", annotation_column,
            "' missing; skipping B-cell spatial plots")
    return(invisible(NULL))
  }
  ct <- as.character(meta[[annotation_column]])
  is_bcell <- grepl(bcell_regex, ct, ignore.case = TRUE)
  if (!any(is_bcell)) {
    message("  \u26A0 No cells match bcell_regex; skipping B-cell spatial plots")
    return(invisible(NULL))
  }

  # (A) Highlight plot — two-level factor
  flag_df <- data.frame(
    cell_ID = meta$cell_ID,
    .bcell_flag = factor(ifelse(is_bcell, "B cell", "Other"),
                         levels = c("B cell", "Other")),
    stringsAsFactors = FALSE
  )
  gobj_h <- tryCatch(
    .giotto_add_cell_metadata(gobj, new_metadata = flag_df,
                              by_column = TRUE, column_cell_ID = "cell_ID"),
    error = function(e) { message("  \u26A0 addCellMetadata failed: ",
                                  conditionMessage(e)); NULL })
  if (!is.null(gobj_h)) {
    title_h <- paste0(sample_id, " \u2014 B cells highlighted (",
                       annotation_column, ")")
    p_h <- .spat_in_situ_outlined(
      gobj          = gobj_h,
      fill_col      = ".bcell_flag",
      fill_as_factor = TRUE,
      colour_map    = c("B cell" = "#E41A1C", "Other" = "#DDDDDD"),
      legend_title  = "Population",
      title_txt     = title_h
    )
    if (!is.null(p_h)) {
      save_presentation_plot(
        plot     = p_h,
        filename = file.path(results_dir,
                              paste0(sample_id, "_bcells_highlight_spatial.png")),
        width = 20, height = 10, dpi = 300
      )
      cat("  \u2713 B-cell highlight spatial saved\n")
    } else {
      message("  \u26A0 spatInSituPlotPoints unavailable; B-cell highlight skipped")
    }
  }

  # (B) Per-gene polygon plots
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
  keep    <- intersect(candidates, rownames(expr))
  missing <- setdiff(candidates, rownames(expr))
  if (length(missing)) {
    message("  \u2139 Marker genes not on panel: ",
            paste(missing, collapse = ", "))
  }
  if (!length(keep)) return(invisible(NULL))

  plots <- list()
  for (g in keep) {
    title_g <- paste0(sample_id, " \u2014 ", g, " expression")
    p_g <- .spat_in_situ_outlined(
      gobj           = gobj,
      fill_col       = g,
      fill_as_factor = FALSE,
      legend_title   = g,
      title_txt      = title_g
    )
    if (!is.null(p_g)) {
      save_presentation_plot(
        plot     = p_g,
        filename = file.path(results_dir,
                              paste0(sample_id, "_bcell_marker_", g, ".png")),
        width = 14, height = 10, dpi = 300
      )
      plots[[g]] <- p_g
    }
  }
  if (length(plots)) {
    cat("  \u2713 Per-gene polygon plots saved (", length(plots),
        " gene(s))\n", sep = "")
  }

  # (C) Combined multi-panel via patchwork
  if (length(plots) >= 2 && requireNamespace("patchwork", quietly = TRUE)) {
    dims <- tryCatch(optimal_grid_dims(length(plots)),
                     error = function(e) {
                       ncol <- ceiling(sqrt(length(plots)))
                       list(ncol = ncol, nrow = ceiling(length(plots) / ncol))
                     })
    combo <- tryCatch(
      patchwork::wrap_plots(plots, ncol = dims$ncol, nrow = dims$nrow),
      error = function(e) NULL
    )
    if (!is.null(combo)) {
      save_presentation_plot(
        plot     = combo,
        filename = file.path(results_dir,
                              paste0(sample_id, "_bcell_markers_panel.png")),
        width  = dims$ncol * 6,
        height = dims$nrow * 5,
        dpi    = 300
      )
      cat("  \u2713 Combined B-cell marker panel saved\n")
    }
  }
  invisible(plots)
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
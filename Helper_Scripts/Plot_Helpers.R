# Plot_Helpers.R
# Foundation helpers for pipeline visualisations.
# Sourced from 00_Setup.R. Consumed by every pipeline script.
#
# Provides:
#   PER_FOV_LINEWIDTH              constant for per-FOV polygon outlines
#   celltype_palette()             deterministic top-level celltype palette
#   cluster_palette()              deterministic palette for numeric/text clusters
#   plot_cells_polygon()           wrapper over spatInSituPlotPoints
#   plot_cells_polygon_highlight() dual-layer highlight polygon ggplot
#   add_centroid_labels()          white-background centroid labels for UMAP/TSNE/PCA
#   clean_plot_text()              scrubs underscores, dots, em dashes from plot text

## ---- Constants -------------------------------------------------------------

# Per-FOV polygon outline width.
# Per-sample / composite tissue plots retain their existing linewidth.
PER_FOV_LINEWIDTH <- 0.25


## ---- Colour palettes -------------------------------------------------------

# Pool ordering is fixed so that, for a given sorted celltype roster, the
# colour-to-celltype assignment is identical across runs and machines.
# Colours that sit close in hue to "mediumspringgreen" (the B-cell pin) are
# omitted from the pool to avoid accidental visual collisions.
.universal_pool <- c(
  "#E69F00",  # Okabe-Ito orange
  "#56B4E9",  # Okabe-Ito sky blue
  "#F0E442",  # Okabe-Ito yellow
  "#0072B2",  # Okabe-Ito blue
  "#D55E00",  # Okabe-Ito vermilion
  "#CC79A7",  # Okabe-Ito reddish purple
  "#999999",  # Okabe-Ito neutral grey
  "#8DD3C7",  # Set3 teal
  "#BEBADA",  # Set3 lavender
  "#FB8072",  # Set3 salmon
  "#80B1D3",  # Set3 steel
  "#FDB462",  # Set3 tan
  "#B3DE69",  # Set3 sage (distinct from mediumspringgreen in hue)
  "#FCCDE5",  # Set3 pink
  "#D9D9D9",  # Set3 light grey
  "#BC80BD",  # Set3 mauve
  "#FFED6F"   # Set3 pale yellow
)

#' Universal, deterministic celltype palette.
#'
#' Given a discovered set of top-level celltypes, returns a named colour
#' vector with \code{pins} applied verbatim and remaining celltypes assigned
#' deterministically (alphabetical order against the fixed Okabe-Ito + Set3
#' pool). Overflow beyond the pool is filled via \code{colorRampPalette}.
#'
#' @param celltypes character vector of top-level celltypes (post-normalisation)
#' @param pins named character vector of celltype -> colour overrides
#' @return named character vector, one entry per unique input celltype
celltype_palette <- function(celltypes,
                             pins = c("B cell" = "mediumspringgreen")) {

  ct <- unique(stats::na.omit(as.character(celltypes)))
  if (length(ct) == 0) return(stats::setNames(character(0), character(0)))

  pins <- if (is.null(pins)) character(0) else pins
  pinned   <- intersect(ct, names(pins))
  unpinned <- sort(setdiff(ct, pinned))

  pool <- .universal_pool
  if (length(unpinned) > length(pool)) {
    pool <- grDevices::colorRampPalette(pool)(length(unpinned))
  }

  out <- character(length(ct))
  names(out) <- ct
  if (length(pinned)   > 0) out[pinned]   <- unname(pins[pinned])
  if (length(unpinned) > 0) out[unpinned] <- pool[seq_along(unpinned)]
  out
}

#' Deterministic palette for un-named / numeric / text-only clusters.
#'
#' Thin wrapper over \code{generate_cluster_colors()} that sorts the input
#' first so cluster-to-colour assignments are stable across scripts in a run.
#'
#' @param cluster_levels character/numeric vector of cluster labels
#' @param palette RColorBrewer palette name forwarded to generate_cluster_colors
#' @return named character vector
cluster_palette <- function(cluster_levels, palette = "Paired") {
  lev <- sort(unique(as.character(cluster_levels)))
  if (length(lev) == 0) return(stats::setNames(character(0), character(0)))
  if (exists("generate_cluster_colors", mode = "function", inherits = TRUE)) {
    return(generate_cluster_colors(cluster_levels = lev, palette = palette))
  }
  # Fallback if Cluster_Visualisations.R has not been sourced yet.
  if (requireNamespace("RColorBrewer", quietly = TRUE) &&
      palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    cols  <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(max_n, palette)
    )(length(lev))
  } else {
    cols <- scales::hue_pal()(length(lev))
  }
  stats::setNames(cols, lev)
}


## ---- Spatial polygon plotting ---------------------------------------------

#' Standard single-layer polygon plot of cells in tissue space.
#'
#' Thin wrapper over \code{Giotto::spatInSituPlotPoints()} that enforces the
#' per-FOV linewidth rule and accepts an optional celltype/cluster palette.
#'
#' @param gobject         Giotto object
#' @param fill_column     metadata column used for polygon fill
#' @param fill_as_factor  treat the fill column as discrete
#' @param context         "fov" applies PER_FOV_LINEWIDTH; "sample" leaves
#'                        linewidth at the user-supplied value (default 0.5)
#' @param palette         optional named colour vector (celltype or cluster)
#' @param polygon_alpha   fill alpha (default 0.85)
#' @param polygon_line_color outline colour (default "grey20")
#' @param ...             forwarded to spatInSituPlotPoints
plot_cells_polygon <- function(gobject,
                               fill_column,
                               fill_as_factor      = TRUE,
                               context             = c("fov", "sample"),
                               palette             = NULL,
                               polygon_alpha       = 0.85,
                               polygon_line_color  = "grey20",
                               polygon_line_size   = NULL,
                               ...) {
  context <- match.arg(context)
  if (is.null(polygon_line_size)) {
    polygon_line_size <- if (context == "fov") PER_FOV_LINEWIDTH else 0.5
  }

  call_args <- list(
    gobject               = gobject,
    show_polygon          = TRUE,
    polygon_feat_type     = "cell",
    polygon_fill          = fill_column,
    polygon_fill_as_factor = fill_as_factor,
    polygon_alpha         = polygon_alpha,
    polygon_line_color    = polygon_line_color,
    polygon_line_size     = polygon_line_size,
    show_image            = FALSE
  )
  if (!is.null(palette)) {
    call_args$polygon_fill_code <- palette
  }
  extra <- list(...)
  call_args[names(extra)] <- extra

  do.call(Giotto::spatInSituPlotPoints, call_args)
}

#' Dual-layer highlight polygon plot.
#'
#' Draws every cell in \code{background_colour}, then overlays the cells that
#' match \code{highlight_celltype} in \code{highlight_colour}. Used for
#' B-cell / subtype spotlight panels. Consumes pre-built polygon vertex frames
#' (x, y, poly_group, <celltype_col>) so it works in both Giotto and plain
#' ggplot contexts.
#'
#' @param polygon_df        data frame with columns x, y, poly_group, and the
#'                          column named by \code{celltype_col}
#' @param celltype_col      column in polygon_df carrying celltype labels
#' @param highlight_celltype character vector of celltypes to highlight
#' @param background_colour fill for non-highlighted cells
#' @param highlight_colour  fill for highlighted cells
#' @param context           "fov" applies PER_FOV_LINEWIDTH, else 0.15
plot_cells_polygon_highlight <- function(polygon_df,
                                         celltype_col,
                                         highlight_celltype,
                                         background_colour = "grey90",
                                         highlight_colour  = "mediumspringgreen",
                                         context           = c("fov", "sample"),
                                         line_colour       = "grey30") {
  context <- match.arg(context)
  lw <- if (context == "fov") PER_FOV_LINEWIDTH else 0.15

  if (!all(c("x", "y", "poly_group", celltype_col) %in% names(polygon_df))) {
    stop("polygon_df must contain columns: x, y, poly_group, and ", celltype_col)
  }

  is_hit <- as.character(polygon_df[[celltype_col]]) %in% as.character(highlight_celltype)
  bg <- polygon_df[!is_hit, , drop = FALSE]
  fg <- polygon_df[is_hit, , drop = FALSE]

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data   = bg,
      mapping = ggplot2::aes(x = x, y = y, group = poly_group),
      fill   = background_colour,
      colour = line_colour,
      linewidth = lw
    ) +
    ggplot2::geom_polygon(
      data   = fg,
      mapping = ggplot2::aes(x = x, y = y, group = poly_group),
      fill   = highlight_colour,
      colour = line_colour,
      linewidth = lw
    ) +
    ggplot2::coord_equal()
}


## ---- Centroid labels for UMAP / TSNE / PCA --------------------------------

#' Add white-background centroid labels to a dimensionality-reduction plot.
#'
#' Computes median-based centroids per group, then returns a list of ggplot
#' layers you append with \code{+}. Uses \code{ggrepel::geom_label_repel} so
#' labels avoid overlap; the label box is filled white for contrast.
#'
#' @param plot_data data frame with at least group_col, x_col, y_col
#' @param group_col name of the grouping column (cluster or celltype)
#' @param x_col x-coordinate column name
#' @param y_col y-coordinate column name
#' @param label_size text size
#' @param box_padding ggrepel padding
add_centroid_labels <- function(plot_data,
                                group_col,
                                x_col       = "dim1",
                                y_col       = "dim2",
                                label_size  = 3,
                                box_padding = 0.4) {

  if (!all(c(group_col, x_col, y_col) %in% names(plot_data))) {
    stop("plot_data must contain columns: ", group_col, ", ", x_col, ", ", y_col)
  }

  df <- plot_data[, c(group_col, x_col, y_col), drop = FALSE]
  df <- df[stats::complete.cases(df), , drop = FALSE]

  centroids <- stats::aggregate(
    df[, c(x_col, y_col)],
    by  = list(group = df[[group_col]]),
    FUN = stats::median
  )
  names(centroids) <- c(group_col, x_col, y_col)

  list(
    ggrepel::geom_label_repel(
      data          = centroids,
      mapping       = ggplot2::aes(
        x     = .data[[x_col]],
        y     = .data[[y_col]],
        label = .data[[group_col]]
      ),
      fill          = "white",
      colour        = "black",
      size          = label_size,
      fontface      = "bold",
      label.size    = 0.2,
      label.padding = grid::unit(0.2, "lines"),
      box.padding   = box_padding,
      min.segment.length = 0,
      segment.colour = "grey40",
      inherit.aes   = FALSE
    )
  )
}


## ---- Plot text cleaner -----------------------------------------------------

#' Replace underscores and dots with spaces and strip em dashes.
#'
#' Vectorised and safe on NA / non-character input. Used anywhere a pipeline
#' plot needs user-facing text derived from a column name, celltype label,
#' gene symbol, or similar.
#'
#' @param x input vector
#' @return character vector, same length as input; NA preserved
clean_plot_text <- function(x) {
  if (is.null(x) || length(x) == 0) return(x)
  out <- as.character(x)
  out[is.na(x)] <- NA_character_
  not_na <- !is.na(out)
  # Replace underscores and dots with a space, collapse any em dashes to hyphen.
  out[not_na] <- gsub("[_.]+", " ", out[not_na], perl = TRUE)
  out[not_na] <- gsub("—", "-", out[not_na], fixed = TRUE)  # em dash
  out[not_na] <- gsub("–", "-", out[not_na], fixed = TRUE)  # en dash
  out[not_na] <- gsub("\\s+", " ", out[not_na], perl = TRUE)
  out[not_na] <- trimws(out[not_na])
  out
}

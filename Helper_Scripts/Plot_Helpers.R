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


## ---- Sentence-case label utilities ----------------------------------------

# Tokens preserved verbatim (treated as abbreviations / proper nouns).
.DISPLAY_PRESERVE_TOKENS <- c(
  "UMAP", "PCA", "tSNE", "t-SNE", "FOV", "FOVs", "L-R", "LR",
  "QC", "RNA", "DNA", "HCA", "CART", "CCI", "DE", "GSEA", "ORA",
  "T0", "T12", "S1", "S2", "BCR", "TCR", "MHC", "ID", "IDs", "FDR",
  "log2FC", "log2", "CD4", "CD8", "B", "T", "NK", "GO", "KEGG",
  "MSigDB", "GO:BP", "PROGENy", "II", "III", "IV"
)

#' Sentence-case a string while preserving abbreviations.
#'
#' First word capitalised, all other words lowercased, except all-caps tokens
#' (>=2 chars, all-letters) and tokens in .DISPLAY_PRESERVE_TOKENS, which
#' are kept verbatim. Words after a period are also capitalised. Single-letter
#' tokens like "B" or "T" are preserved as-is when they're flanked by
#' another label-token ("B cell", "T cell").
#'
#' @param text character vector (length >=1); NA-safe
#' @return character vector, same length
.title_case_sentence <- function(text) {
  if (is.null(text) || length(text) == 0) return(text)
  out <- as.character(text)
  out[is.na(text)] <- NA_character_

  cap_word <- function(w) {
    if (!nzchar(w)) return(w)
    paste0(toupper(substring(w, 1, 1)), tolower(substring(w, 2)))
  }

  process_one <- function(s) {
    if (is.na(s) || !nzchar(s)) return(s)

    # Split on whitespace into tokens, preserve original spacing on rejoin.
    tokens <- strsplit(s, "(?<=\\s)|(?=\\s)", perl = TRUE)[[1]]
    word_idx <- which(!grepl("^\\s+$", tokens))
    if (length(word_idx) == 0) return(s)

    # Track whether the previous non-space token ended with a period
    # (or this is the very first word) so we capitalise sentence starts.
    cap_next <- TRUE
    for (i in word_idx) {
      tok <- tokens[i]

      # Preserve verbatim if it matches a known token (case-sensitive) OR
      # the alphabetic part is all uppercase (>=2 chars): e.g. "(UMAP)".
      alpha <- gsub("[^A-Za-z]", "", tok)
      preserve <- tok %in% .DISPLAY_PRESERVE_TOKENS ||
                  alpha %in% .DISPLAY_PRESERVE_TOKENS ||
                  (nchar(alpha) >= 2 && alpha == toupper(alpha) &&
                   alpha != tolower(alpha))

      if (!preserve) {
        if (cap_next) {
          tokens[i] <- cap_word(tok)
        } else {
          tokens[i] <- tolower(tok)
        }
      }

      # Sentence-end detection: token ends with `.`, `?`, or `!` (not part
      # of an abbreviation we just preserved).
      cap_next <- !preserve && grepl("[.!?]\\s*$", tok)
      # Always reset to TRUE for the very first word after a sentence end.
    }
    paste(tokens, collapse = "")
  }

  vapply(out, process_one, character(1), USE.NAMES = FALSE)
}

#' Display label: clean punctuation + sentence case.
#'
#' Standard transform applied to every user-facing legend title, axis title,
#' and plot title authored by pipeline code. Combines clean_plot_text() with
#' .title_case_sentence() so the same input ("leiden_clust") becomes the
#' same output ("Leiden cluster") regardless of where it originated.
#'
#' @param text input string or vector
#' @return character vector, same length as input
.display_label <- function(text) {
  .title_case_sentence(clean_plot_text(text))
}


## ---- Per-FOV plot looper --------------------------------------------------

#' Emit one PNG per FOV that contains at least one focus-celltype cell.
#'
#' Generic looper used by scripts 07/08/09/10/11. Iterates the FOVs that
#' contain >=1 cell matching target_celltype_regex (or every FOV if regex is
#' NULL), filters poly_df to that FOV, calls plot_fn(poly_df_fov, fov_id) to
#' build a ggplot, then saves to <out_dir>/<fname_base>_FOV_<fov>.png. Both
#' filename and figure title carry sample_id + FOV.
#'
#' Skips FOVs with too few cells (default <10) to avoid noisy near-empty plots.
#'
#' @param poly_df    polygon vertex frame from .extract_polygon_df(), must have
#'                   columns x, y, cell_ID, geom, part (or equivalent group)
#'                   AND a "fov" column (joined upstream from cell metadata)
#' @param meta       cell metadata data.frame (must include cell_ID, fov, and
#'                   the celltype column referenced by target_celltype_regex)
#' @param celltype_col name of the celltype column in meta (e.g. "celltype")
#' @param target_celltype_regex regex; if NULL, every FOV is included
#' @param plot_fn    function(poly_df_fov, fov_id) -> ggplot
#' @param out_dir    directory; will be created
#' @param fname_base filename prefix; FOV id is appended as "_FOV_<n>.png"
#' @param sample_id  used in the filename prefix and the figure title
#' @param title_prefix human-readable title prefix (e.g. "Annotation"); the
#'                   sample id and FOV id are appended automatically
#' @param min_cells  skip FOVs with fewer than this many target cells (default 5)
#' @param width,height,dpi forwarded to save_presentation_plot
#' @return character vector of the file paths written (invisibly)
.emit_per_fov_plots <- function(poly_df, meta, celltype_col,
                                target_celltype_regex,
                                plot_fn, out_dir, fname_base, sample_id,
                                title_prefix = NULL,
                                min_cells = 5,
                                width = 10, height = 10, dpi = 600) {

  if (is.null(poly_df) || !nrow(poly_df) || !"fov" %in% names(poly_df)) {
    return(invisible(character(0)))
  }
  if (is.null(meta) || !"fov" %in% names(meta)) {
    return(invisible(character(0)))
  }

  is_target <- if (is.null(target_celltype_regex) ||
                   !nzchar(target_celltype_regex) ||
                   !celltype_col %in% names(meta)) {
    rep(TRUE, nrow(meta))
  } else {
    grepl(target_celltype_regex, as.character(meta[[celltype_col]]),
          ignore.case = TRUE, perl = TRUE)
  }

  meta_target <- meta[is_target & !is.na(meta$fov), , drop = FALSE]
  if (!nrow(meta_target)) return(invisible(character(0)))

  fov_counts <- table(meta_target$fov)
  fovs <- names(fov_counts)[fov_counts >= min_cells]
  if (!length(fovs)) return(invisible(character(0)))

  ensure_dir(out_dir)
  written <- character(0)

  for (fv in fovs) {
    df_fv <- poly_df[!is.na(poly_df$fov) & poly_df$fov == fv, , drop = FALSE]
    if (!nrow(df_fv)) next

    p <- tryCatch(plot_fn(df_fv, fv),
                  error = function(e) {
                    message("  [per-FOV] FOV ", fv, " plot failed: ",
                            conditionMessage(e))
                    NULL
                  })
    if (is.null(p)) next

    title_text <- paste0(
      sample_id, " - ",
      if (!is.null(title_prefix) && nzchar(title_prefix))
        paste0(title_prefix, " - ") else "",
      "FOV ", fv
    )
    p <- p + ggplot2::labs(title = title_text)

    fname <- file.path(out_dir,
                       paste0(fname_base, "_FOV_", fv, ".png"))
    save_presentation_plot(p, fname,
                           width = width, height = height, dpi = dpi)
    written <- c(written, fname)
  }

  invisible(written)
}


## ---- Shared proximity heatmap (ComplexHeatmap-based) -----------------------

#' Render a celltype-by-celltype proximity heatmap using ComplexHeatmap.
#'
#' Shared by 09_Spatial_Network and 11_B_Cell_Analysis so both heatmaps look
#' identical in style. 09 passes the full enrichment matrix; 11 passes the
#' same matrix but flags the focus celltype to be reordered/highlighted.
#'
#' Falls back to pheatmap if ComplexHeatmap is unavailable.
#'
#' @param mat            square numeric matrix (celltype x celltype enrichment Z)
#' @param title          plot title string
#' @param subtitle       optional subtitle (NULL = none)
#' @param focus_label    if non-NULL, this row/col is moved to position 1 and
#'                       highlighted with a red border annotation
#' @param filename       output PNG path
#' @param width,height,dpi save dimensions
#' @param colour_breaks  numeric vector of break points (default symmetric)
.plot_proximity_heatmap_complex <- function(mat,
                                            title,
                                            subtitle = NULL,
                                            focus_label = NULL,
                                            filename,
                                            width = 12,
                                            height = 10,
                                            dpi = 300,
                                            colour_breaks = NULL) {

  if (is.null(mat) || !nrow(mat) || !ncol(mat)) return(invisible(NULL))

  # Reorder to put focus row/col first if requested.
  if (!is.null(focus_label) && nzchar(focus_label)) {
    matched <- which(rownames(mat) == focus_label)
    if (length(matched) == 1L) {
      ord <- c(matched, setdiff(seq_len(nrow(mat)), matched))
      mat <- mat[ord, ord, drop = FALSE]
    }
  }

  display_rownames <- .display_label(rownames(mat))
  display_colnames <- .display_label(colnames(mat))

  use_complex <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
                 requireNamespace("circlize", quietly = TRUE)

  if (use_complex) {
    if (is.null(colour_breaks)) {
      mx <- max(abs(mat), na.rm = TRUE)
      mx <- if (!is.finite(mx) || mx == 0) 1 else mx
      colour_breaks <- c(-mx, 0, mx)
    }
    col_fn <- circlize::colorRamp2(colour_breaks,
                                   c("#2166AC", "white", "#B2182B"))

    # Build optional row/col annotation that highlights the focus celltype.
    top_anno <- NULL
    left_anno <- NULL
    if (!is.null(focus_label) && nzchar(focus_label) &&
        focus_label %in% rownames(mat)) {
      hl <- rep("other", nrow(mat))
      hl[1] <- "focus"
      top_anno <- ComplexHeatmap::HeatmapAnnotation(
        focus = hl,
        col = list(focus = c(focus = "#B2182B", other = "grey90")),
        show_legend = FALSE,
        annotation_name_gp = grid::gpar(fontsize = 0)
      )
      left_anno <- ComplexHeatmap::rowAnnotation(
        focus = hl,
        col = list(focus = c(focus = "#B2182B", other = "grey90")),
        show_legend = FALSE,
        annotation_name_gp = grid::gpar(fontsize = 0)
      )
    }

    rownames(mat) <- display_rownames
    colnames(mat) <- display_colnames

    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = "Enrichment",
      col = col_fn,
      cluster_rows = is.null(focus_label),
      cluster_columns = is.null(focus_label),
      rect_gp = grid::gpar(col = "white", lwd = 0.5),
      row_names_side = "left",
      column_names_side = "bottom",
      column_names_rot = 45,
      row_names_gp = grid::gpar(fontsize = 10),
      column_names_gp = grid::gpar(fontsize = 10),
      column_title = title,
      column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
      row_title = if (!is.null(subtitle)) subtitle else NULL,
      row_title_gp = grid::gpar(fontsize = 11, col = "grey20"),
      heatmap_legend_param = list(
        title = "Enrichment",
        legend_direction = "vertical"
      ),
      top_annotation = top_anno,
      left_annotation = left_anno
    )

    grDevices::png(filename = filename, width = width, height = height,
                   units = "in", res = dpi, type = "cairo", bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    ComplexHeatmap::draw(
      ht,
      padding = grid::unit(c(20, 30, 20, 30), "mm"),
      heatmap_legend_side = "right"
    )
    return(invisible(filename))
  }

  # ---- pheatmap fallback ---------------------------------------------------
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("  [proximity-heatmap] neither ComplexHeatmap nor pheatmap available; skipping")
    return(invisible(NULL))
  }
  rownames(mat) <- display_rownames
  colnames(mat) <- display_colnames
  grDevices::png(filename = filename, width = width, height = height,
                 units = "in", res = dpi, type = "cairo", bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  pheatmap::pheatmap(
    mat,
    main = title,
    cellwidth = 18,
    cellheight = 18,
    fontsize_row = 10,
    fontsize_col = 10,
    cluster_rows = is.null(focus_label),
    cluster_cols = is.null(focus_label),
    color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  )
  invisible(filename)
}

#!/usr/bin/env Rscript
# ============================================================================== 
# 09_Spatial_Network.R
# Spatial network construction, neighbourhood enrichment, and spot-level
# niche deconvolution — all self-contained within Giotto.
#
# Column conventions in cellProximityEnrichment output:
#   unified_int  = cell-type pair label string  e.g. "TypeA--TypeB"
#   type_int     = "homo" (self-self) or "hetero" (cross-type) flag
#
# Outputs (under output_dir/09_Spatial_Network/):
#   ├── networks/       Spatial network edge CSVs
#   ├── proximity/      Enrichment table + heatmap + vis plots
#   ├── deconvolution/  PAGE / rank enrichment scores + spatial plots
#   └── spatial_plots/  (reserved for additional spatial plots)
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


# ============================================================================== 
# Internal helper: proximity vis plot (bypasses cellProximityVisPlot Giotto bug) 
# ============================================================================== 

#' @keywords internal
.plot_cell_proximity <- function(gobj,
                                 interaction_name,
                                 celltype_col,
                                 primary_network,
                                 sample_id,
                                 prox_folder,
                                 file_label = NULL) {

  ann_net <- annotateSpatialNetwork(
    gobject              = gobj,
    spatial_network_name = primary_network,
    cluster_column       = celltype_col
  )

  spat_locs <- getSpatialLocations(gobj, output = "data.table")
  meta      <- pDataDT(gobj)

  # Edge endpoints (centroid-to-centroid) - unchanged.
  int_edges <- ann_net[ann_net$unified_int == interaction_name, ]
  edge_dt   <- int_edges[, c("sdimx_begin", "sdimy_begin",
                             "sdimx_end",   "sdimy_end")]

  ct_pair <- strsplit(interaction_name, "--")[[1]]
  in_pair <- meta$cell_ID[as.character(meta[[celltype_col]]) %in% ct_pair]

  # Polygon geometry for cell bodies (replaces geom_point so the plot
  # matches the pipeline-wide polygon convention). Falls back silently
  # to the previous point rendering if polygons are unavailable.
  poly_df <- tryCatch(.extract_polygon_df(gobj), error = function(e) NULL)

  ct_colours <- setNames(scales::hue_pal()(length(ct_pair)), ct_pair)

  if (!is.null(poly_df) && nrow(poly_df) > 0) {
    poly_df$poly_group <- paste(poly_df$cell_ID, poly_df$geom,
                                 poly_df$part, sep = "_")
    poly_df$cell_type  <- as.character(
      meta[[celltype_col]][match(poly_df$cell_ID, meta$cell_ID)]
    )
    poly_df$in_pair <- poly_df$cell_ID %in% in_pair

    bg_poly <- poly_df[!poly_df$in_pair, , drop = FALSE]
    fg_poly <- poly_df[ poly_df$in_pair, , drop = FALSE]

    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(
        data    = bg_poly,
        mapping = ggplot2::aes(x = x, y = y, group = poly_group),
        fill    = "grey90", colour = "grey75", linewidth = 0.06
      ) +
      ggplot2::geom_segment(
        data    = edge_dt,
        mapping = ggplot2::aes(x = sdimx_begin, y = sdimy_begin,
                               xend = sdimx_end, yend = sdimy_end),
        colour  = "grey40", linewidth = 0.2, alpha = 0.6
      ) +
      ggplot2::geom_polygon(
        data    = fg_poly,
        mapping = ggplot2::aes(x = x, y = y, group = poly_group,
                               fill = cell_type),
        colour  = "grey30", linewidth = 0.15, alpha = 0.9
      ) +
      ggplot2::scale_fill_manual(
        values = ct_colours,
        labels = function(x) pretty_plot_label(x),
        name   = "Cell Type"
      )
  } else {
    cell_dt <- merge(
      spat_locs[, c("cell_ID", "sdimx", "sdimy")],
      meta[,     c("cell_ID", celltype_col), with = FALSE],
      by = "cell_ID"
    )
    names(cell_dt)[names(cell_dt) == celltype_col] <- "cell_type"
    cell_dt$in_pair <- cell_dt$cell_type %in% ct_pair
    selected_dt <- cell_dt[cell_dt$in_pair, ]
    other_dt    <- cell_dt[!cell_dt$in_pair, ]
    p <- ggplot2::ggplot() +
      ggplot2::geom_point(
        data = other_dt, mapping = ggplot2::aes(x = sdimx, y = sdimy),
        colour = "grey85", size = 0.25, alpha = 0.35
      ) +
      ggplot2::geom_segment(
        data = edge_dt,
        mapping = ggplot2::aes(x = sdimx_begin, y = sdimy_begin,
                               xend = sdimx_end, yend = sdimy_end),
        colour = "grey40", linewidth = 0.2, alpha = 0.5
      ) +
      ggplot2::geom_point(
        data = selected_dt,
        mapping = ggplot2::aes(x = sdimx, y = sdimy, colour = cell_type),
        size = 1.4, alpha = 0.9
      ) +
      ggplot2::scale_colour_manual(
        values = ct_colours,
        labels = function(x) pretty_plot_label(x),
        name   = "Cell Type"
      )
  }

  # Axis / margin theme prevents axis text being cut off on composite figures.
  p <- p +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = sample_plot_title(sample_id,
        paste0("Cell Proximity: ", pretty_plot_label(interaction_name))),
      x = "Global X Coordinate", y = "Global Y Coordinate"
    ) +
    presentation_theme(base_size = 12, legend_position = "right") +
    ggplot2::theme(
      legend.key.height = grid::unit(0.45, "cm"),
      axis.title.x      = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y      = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      plot.margin       = ggplot2::margin(t = 10, r = 20, b = 20, l = 20)
    )

  safe_name <- gsub("[^A-Za-z0-9]", "_", interaction_name)
  label_tag <- if (!is.null(file_label) && nzchar(file_label)) file_label else sample_id
  save_path <- file.path(prox_folder,
                         paste0(label_tag, "_proximity_visplot_",
                                safe_name, ".png"))
  save_presentation_plot(p, save_path, width = 14, height = 10, dpi = 150)
  cat("\u2713 Proximity vis plot saved:", basename(save_path), "\n")
  invisible(p)
}


# ==============================================================================
# Internal helper: niche-deconvolution polygon plot (replaces spatPlot2D)
# ==============================================================================
#
# Uses .extract_polygon_df() + ggplot2 instead of Giotto's spatPlot2D so that
# the deconvolution score layer rendres with the same polygon convention as
# the rest of the pipeline, gains a named legend, and gets explicit plot /
# axis margins so axis text is never clipped. Falls back to spatPlot2D for
# samples without polygon data.

.plot_deconv_polygons <- function(gobj,
                                  score_col,
                                  enrichment_name,
                                  sample_id,
                                  deconv_folder,
                                  legend_name = NULL,
                                  file_label  = NULL) {

  poly_df <- tryCatch(.extract_polygon_df(gobj), error = function(e) NULL)
  if (is.null(poly_df) || nrow(poly_df) == 0) return(NULL)
  poly_df$poly_group <- paste(poly_df$cell_ID, poly_df$geom,
                               poly_df$part, sep = "_")

  enr <- tryCatch(
    getSpatialEnrichment(gobj, name = enrichment_name, output = "data.table"),
    error = function(e) NULL
  )
  if (is.null(enr) || !(score_col %in% names(enr))) return(NULL)

  poly_df$score <- as.numeric(
    enr[[score_col]][match(poly_df$cell_ID, enr$cell_ID)]
  )

  legend_title <- if (!is.null(legend_name) && nzchar(legend_name)) {
    legend_name
  } else {
    paste(pretty_plot_label(score_col), "score")
  }

  p <- ggplot2::ggplot(poly_df,
                       ggplot2::aes(x = x, y = y,
                                    group = poly_group, fill = score)) +
    ggplot2::geom_polygon(colour = "grey40", linewidth = 0.08) +
    ggplot2::scale_fill_gradient(low = "lightgrey", high = "red",
                                 name = legend_title, na.value = "grey85") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = sample_plot_title(sample_id,
        paste0(enrichment_name, ": ", pretty_plot_label(score_col))),
      x = "Global X Coordinate", y = "Global Y Coordinate"
    ) +
    presentation_theme(base_size = 11, legend_position = "right") +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      plot.margin  = ggplot2::margin(t = 10, r = 20, b = 20, l = 20),
      legend.key.height = grid::unit(0.45, "cm")
    )

  safe_ct <- gsub("[^A-Za-z0-9]", "_", score_col)
  label_tag <- if (!is.null(file_label) && nzchar(file_label)) file_label else sample_id
  save_path <- file.path(deconv_folder,
                         paste0(label_tag, "_", enrichment_name, "_",
                                safe_ct, ".png"))
  save_presentation_plot(p, save_path, width = 12, height = 10, dpi = 200)
  invisible(save_path)
}


# ==============================================================================
# Internal helper: proximity heatmap via pheatmap (no ComplexHeatmap needed)
# ==============================================================================

#' @keywords internal
.plot_proximity_heatmap <- function(prox_csv_path,
                                    sample_id,
                                    prox_folder) {
  
  if (!file.exists(prox_csv_path))
    stop("Proximity enrichment CSV not found: ", prox_csv_path)
  
  enr      <- data.table::fread(prox_csv_path, data.table = FALSE)
  pairs    <- strsplit(enr$unified_int, "--")
  enr$from <- sapply(pairs, `[`, 1)
  enr$to   <- sapply(pairs, `[`, 2)
  
  all_types <- sort(unique(c(enr$from, enr$to)))
  n_types   <- length(all_types)
  
  # Enrichment matrix (symmetrised) — initialise with 0
  mat <- matrix(0, nrow = n_types, ncol = n_types,
                dimnames = list(all_types, all_types))
  
  # Significance mask
  sig_mat <- matrix(FALSE, nrow = n_types, ncol = n_types,
                    dimnames = list(all_types, all_types))
  
  for (i in seq_len(nrow(enr))) {
    r  <- enr$from[i]
    co <- enr$to[i]
    v  <- enr$enrichm[i]
    
    is_sig <- (!is.na(enr$p.adj_higher[i]) & enr$p.adj_higher[i] < 0.05) |
      (!is.na(enr$p.adj_lower[i])  & enr$p.adj_lower[i]  < 0.05)
    
    mat[r, co]     <- v
    mat[co, r]     <- v       # symmetrise
    sig_mat[r, co] <- is_sig
    sig_mat[co, r] <- is_sig
  }
  
  # Compute colour scale from FULL matrix before masking
  max_abs <- max(abs(mat[is.finite(mat)]), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  # Grey out non-significant cells AFTER computing max_abs
  mat_masked           <- mat
  mat_masked[!sig_mat] <- NA
  
  row_order <- seq_len(n_types)
  col_order <- seq_len(n_types)
  if (n_types >= 2) {
    row_order <- tryCatch(
      stats::hclust(stats::dist(mat))$order,
      error = function(e) seq_len(n_types)
    )
    col_order <- tryCatch(
      stats::hclust(stats::dist(t(mat)))$order,
      error = function(e) seq_len(n_types)
    )
  }
  
  mat_masked <- mat_masked[row_order, col_order, drop = FALSE]
  
  breaks <- seq(-max_abs, max_abs, length.out = 101)
  colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)
  
  # Row/col annotation: immune vs structural compartment.
  # Labels match the smiDE-native scheme produced by 07_Annotation.R's
  # .normalize_label() (dots/underscores → spaces).
  immune_types_raw <- c("B.cell", "CD4.T.cell", "CD8.T.cell", "NK.cell",
                        "NKT.cell",
                        "MNP.macrophage", "MNP.monocyte",
                        "MNP.a.classical.monocyte.derived",
                        "MNP.b.non.classical.monocyte.derived",
                        "MNP.c.dendritic.cell", "MNP.d.Tissue.macrophage",
                        "Plasmacytoid.DC", "Plasmacytoid.dendritic.cell",
                        "Mast.cell", "Neutrophil")
  immune_types <- trimws(gsub("\\s+", " ",
                              gsub("[._]+", " ", immune_types_raw)))
  # Accept both schemes for legacy checkpoints that predate normalization.
  immune_types <- unique(c(immune_types, immune_types_raw))
  
  anno_df <- data.frame(
    Compartment = ifelse(all_types %in% immune_types, "Immune", "Structural"),
    row.names   = all_types
  )
  anno_df <- anno_df[rownames(mat_masked), , drop = FALSE]
  anno_colours <- list(
    Compartment = c(Immune     = "#E41A1C",
                    Structural = "#377EB8")
  )
  
  save_path <- file.path(prox_folder,
                         paste0(sample_id, "_proximity_heatmap.png"))
  
  grDevices::png(save_path, width = 16, height = 14, units = "in", res = 150)
  pheatmap::pheatmap(
    mat_masked,
    color             = colors,
    breaks            = breaks,
    na_col            = "grey92",
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = TRUE,
    show_colnames     = TRUE,
    fontsize          = 7.5,
    fontsize_row      = 7.5,
    fontsize_col      = 7.5,
    angle_col         = 45,
    border_color      = NA,
    annotation_row    = anno_df,
    annotation_col    = anno_df,
    annotation_colors = anno_colours,
    main              = paste0(
      sample_id,
      "  \u2014  Cell-type neighbourhood enrichment\n",
      "grey = not significant (FDR \u2265 0.05)  |  ",
      "blue = spatial exclusion  |  red = co-localisation"
    ),
    legend_breaks = c(-max_abs, -max_abs / 2, 0, max_abs / 2, max_abs),
    legend_labels = c(
      sprintf("%.2f  (strong exclusion)",   -max_abs),
      sprintf("%.2f",                       -max_abs / 2),
      "0",
      sprintf("%.2f",                        max_abs / 2),
      sprintf("%.2f  (strong co-loc)",       max_abs)
    )
  )
  grDevices::dev.off()
  
  cat("\u2713 Proximity heatmap saved:", basename(save_path), "\n")
  invisible(save_path)
}


# ============================================================================== 
# Main function 
# ============================================================================== 

#' Build Spatial Networks, Run Neighbourhood Enrichment, and Niche Deconvolution
#'
#' @param gobj              Annotated Giotto object, or path to saved object
#' @param sample_id         Sample identifier string
#' @param output_dir        Root output directory for this sample
#' @param celltype_col      Cell metadata column containing cell type labels.
#'                          Auto-detects first "celltype_.*_supervised$" column;
#'                          falls back to "leiden_clust".
#' @param marker_list       Named list: list(CellType = c("Gene1","Gene2"), ...)
#'                          If NULL the built-in HCA Kidney set is used.
#' @param network_methods   Subset of c("Delaunay","kNN"). Default: both.
#' @param knn_k             k for kNN spatial network (default: 6)
#' @param n_simulations     Permutations for cellProximityEnrichment (default: 1000)
#' @param enrich_methods    Subset of c("PAGE","rank"). Default: both.
#' @param expression_values Expression slot to use (default: "normalized")
#' @param create_plots      Save output plots? (default: TRUE)
#' @param random_seed       Permutation seed (default: 42)
#' @param min_page_genes    Min genes per PAGE signature (default: 5)
#'
#' @return Giotto object with spatial networks and enrichment scores attached

build_spatial_network <- function(gobj,
                                  sample_id,
                                  output_dir,
                                  celltype_col           = NULL,
                                  marker_list            = NULL,
                                  network_methods        = c("Delaunay", "kNN"),
                                  knn_k                  = 6,
                                  n_simulations          = 1000,
                                  enrich_methods         = c("PAGE", "rank"),
                                  expression_values      = "normalized",
                                  create_plots           = TRUE,
                                  random_seed            = 42,
                                  min_page_genes         = 5,
                                  max_delaunay_distance  = 50,
                                  sample_row             = NULL,
                                  sample_sheet_path      = NULL) {
  
  set.seed(random_seed)
  
  cat("\n========================================\n")
  cat("STEP 09: Spatial Network Analysis\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    manifest_path <- file.path(gobj, "manifest.json")
    if (file.exists(manifest_path)) {
      gobj <- load_giotto_checkpoint(gobj)
    } else {
      cat("Loading Giotto object from:", gobj, "\n")
      gobj <- loadGiotto(gobj)
    }
    cat("\u2713 Loaded\n\n")
  }
  
  base_folder   <- file.path(output_dir, "09_Spatial_Network")
  net_folder    <- file.path(base_folder, "networks")
  prox_folder   <- file.path(base_folder, "proximity")
  deconv_folder <- file.path(base_folder, "deconvolution")
  plot_folder   <- file.path(base_folder, "spatial_plots")
  
  for (d in c(net_folder, prox_folder, deconv_folder, plot_folder))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  
  meta_cols <- names(pDataDT(gobj))
  if (is.null(celltype_col)) {
    sup_cols     <- grep("celltype_.*_supervised$", meta_cols, value = TRUE)
    celltype_col <- if (length(sup_cols) > 0) sup_cols[1] else "leiden_clust"
    cat("Cell type column auto-selected:", celltype_col, "\n")
  } else {
    if (!celltype_col %in% meta_cols)
      stop("celltype_col '", celltype_col, "' not found in cell metadata.")
  }
  
  n_types <- length(unique(pDataDT(gobj)[[celltype_col]]))
  cat("Cell types:", n_types, "\n\n")
  
  if (is.null(marker_list)) {
    cat("No marker_list supplied — using built-in HCA Kidney markers.\n\n")
    marker_list <- list(
      Proximal.tubule             = c("SLC3A1","SLC7A9","CUBN","LRP2","UMOD","SLC34A1","GATM","ALDOB","PCK1","SLC17A3"),
      Distinct.proximal.tubule.1  = c("SLC17A3","SLC22A8","SLC22A6","ALDOB","GATM","SLC7A13","PRODH2","SLC16A9","SLC4A4","GPX3"),
      Distinct.proximal.tubule.2  = c("SLC34A1","SLC13A3","SLC5A12","SLC17A1","MIOX","FBP1","G6PC","SLC2A2","PDZK1","TMEM174"),
      Connecting.tubule           = c("CALB1","AQP2","SCNN1A","SCNN1B","SCNN1G","TRPV5","SLC8A1","ATP6V1B1","SLC26A7","KRT8"),
      Principal.cell              = c("AQP2","AQP3","HSD11B2","FXYD4","SCNN1A","SCNN1B","SCNN1G","HSPB1","ATP1B1","CLDN8"),
      Type.A.intercalated.cell    = c("ATP6V1B1","SLC4A1","FOXI1","ATP6V0D2","ATP6V1C2","SLC26A7","KIT","INSRR","CA2","DMBT1"),
      Type.B.intercalated.cell    = c("SLC26A4","CLDN4","AQP6","SLC4A9","HMGA1","CA2","LINC01125","ADRB2","TMEM213","CLNK"),
      Thick.ascending.limb        = c("SLC12A1","UMOD","CLDN16","CLDN19","CASR","KCNJ1","PVALB","SLC9A2","PTGER3","FGF13"),
      Myofibroblast               = c("ACTA2","PDGFRB","COL1A1","COL1A2","COL3A1","TAGLN","MYH11","CNN1","POSTN","NOTCH3"),
      Fibroblast                  = c("PDGFRA","DCN","LUM","COL3A1","COL1A1","COL1A2","MFAP5","C3","CFD","SFRP2"),
      Glomerular.endothelium      = c("PECAM1","EHD3","PLAT","TMEM204","EMCN","PODXL2","HYAL2","CLDN5","KDR","ESAM"),
      Peritubular.cap.endothelium = c("PECAM1","FLT1","KDR","EMCN","ERG","ENG","VWF","PTPRB","ADGRL4","CLEC14A"),
      Ascending.vasa.recta        = c("PLVAP","ENG","PECAM1","IGFBP3","SELE","SELP","VCAM1","ICAM1","JAG1","HEY1"),
      Descending.vasa.recta       = c("PECAM1","IGFBP5","EFNB2","DLL4","NOTCH4","KRT18","CXCL12","GJA5","GJA4","HEY2"),
      Podocyte                    = c("NPHS1","NPHS2","WT1","PODXL","SYNPO","PTPRO","MAFB","TCF21","CLIC5","FGF1"),
      Pelvic.epithelium           = c("UPK1B","UPK2","UPK3A","UPK1A","KRT13","KRT5","PSCA","CD44","KRT7","CLDN3"),
      Epithelial.progenitor       = c("CD24","PROM1","SOX9","KRT19","ALDH1A1","EPCAM","VIM","KRT8","KRT18","CDH1"),
      CD4.T.cell                  = c("CD4","IL7R","CCR7","SELL","TCF7","LEF1","CD3D","CD3E","CD3G","CD5"),
      CD8.T.cell                  = c("CD8A","CD8B","GZMB","PRF1","NKG7","GNLY","GZMA","GZMH","IFNG","FASLG"),
      NK.cell                     = c("NCAM1","NKG7","GNLY","KLRD1","KLRB1","KLRC1","KLRK1","NCR1","NCR3","XCL1"),
      B.cell                      = c("MS4A1","CD19","CD79A","PAX5","CD79B","BANK1","BLK","IGHM","FCRL1","FCER2"),
      MNP.macrophage              = c("CD68","MRC1","C1QA","C1QB","C1QC","FOLR2","LYVE1","TIMD4","MARCO","F13A1"),
      MNP.monocyte                = c("CD14","LYZ","S100A8","FCN1","S100A9","CD163","CD36","VCAN","CLEC7A","ITGAM"),
      Plasmacytoid.DC             = c("LILRA4","CLEC4C","IL3RA","IRF7","TCF4","SELL","PTPRS","CUX2","SERPINF1","LRRC26"),
      Mast.cell                   = c("TPSAB1","CPA3","KIT","HDC","FCER1A","MS4A2","HPGDS","GATA2","PTGDS","SRGAP3"),
      Neutrophil                  = c("FCGR3B","CXCR2","CSF3R","FFAR2","S100A8","S100A9","G0S2","CXCR1","PROK2","KCNJ15"),
      Proliferating.PT            = c("MKI67","TOP2A","PCNA","CDK1","CCNB1","CCNB2","CCNA2","BUB1","AURKB","PLK1")
    )
  }
  
  avail_genes <- rownames(
    getExpression(gobj, values = expression_values, output = "matrix")
  )
  
  marker_list_filtered <- lapply(marker_list,
                                 function(genes) intersect(genes, avail_genes))
  marker_list_page <- marker_list_filtered[
    sapply(marker_list_filtered, length) >= min_page_genes
  ]
  marker_list_rank <- marker_list_filtered[
    sapply(marker_list_filtered, length) >= 2
  ]
  
  cat(sprintf(
    "Marker list: %d cell types with >= %d genes (PAGE), %d with >= 2 genes (rank)\n\n",
    length(marker_list_page), min_page_genes, length(marker_list_rank)
  ))
  
  if (length(marker_list_page) == 0 && "PAGE" %in% enrich_methods) {
    cat("\u26A0 No cell types passed the PAGE gene threshold.\n")
    cat("  Supply a custom marker_list with more genes per type.\n\n")
  }
  
  
  # ============================================================================ 
  # SECTION 1 — SPATIAL NETWORK CONSTRUCTION 
  # ============================================================================ 
  
  cat("================================================================================\n")
  cat("SECTION 1: Spatial Network Construction\n")
  cat("================================================================================\n\n")
  
  network_names <- character(0)
  
  if ("Delaunay" %in% network_methods) {
    cat(sprintf("Building Delaunay spatial network (max edge length: %g um)...\n",
                max_delaunay_distance))
    tryCatch({
      # Try passing maximum_distance_delaunay (supported in Giotto >= 4.0);
      # if the argument is rejected fall back to an uncapped build and
      # filter edges post-hoc.
      built <- tryCatch({
        createSpatialNetwork(gobject                   = gobj,
                             method                    = "Delaunay",
                             minimum_k                 = 2,
                             maximum_distance_delaunay = max_delaunay_distance,
                             name                      = "Delaunay_network")
      }, error = function(e) {
        cat("  note: maximum_distance_delaunay not supported by this Giotto version; falling back to post-hoc edge filtering\n")
        NULL
      })
      if (is.null(built)) {
        gobj <- createSpatialNetwork(gobject = gobj, method = "Delaunay",
                                     minimum_k = 2, name = "Delaunay_network")
      } else {
        gobj <- built
      }
      network_names <- c(network_names, "Delaunay_network")
      del_net <- getSpatialNetwork(gobj, name = "Delaunay_network",
                                   output = "networkDT")
      if ("distance" %in% names(del_net)) {
        edge_len <- as.numeric(del_net$distance)
      } else {
        edge_len <- sqrt((del_net$sdimx_end - del_net$sdimx_begin)^2 +
                         (del_net$sdimy_end - del_net$sdimy_begin)^2)
      }
      n_before <- nrow(del_net)
      over_cap <- edge_len > max_delaunay_distance
      if (any(over_cap, na.rm = TRUE)) {
        cat(sprintf("  Pruning %d/%d edges exceeding %g um cap (post-hoc)\n",
                    sum(over_cap, na.rm = TRUE), n_before, max_delaunay_distance))
        del_net <- del_net[!over_cap, , drop = FALSE]
      }
      cat(sprintf("  Edges kept: %d  |  Unique nodes: %d  |  max length: %.1f um\n",
                  nrow(del_net),
                  length(unique(c(del_net$from, del_net$to))),
                  if (nrow(del_net)) max(edge_len[!over_cap], na.rm = TRUE) else 0))
      write.csv(del_net,
                file.path(net_folder, paste0(sample_id, "_delaunay_network.csv")),
                row.names = FALSE)
      cat("  \u2713 Saved: Delaunay network edges\n\n")
      cat("\u2713 Delaunay network created\n")
    }, error = function(e) cat("\u26A0 Delaunay failed:", conditionMessage(e), "\n\n"))
  }
  
  if ("kNN" %in% network_methods) {
    cat(sprintf("Building kNN spatial network (k = %d)...\n", knn_k))
    tryCatch({
      gobj <- createSpatialNetwork(gobject = gobj, method = "kNN",
                                   k = knn_k,
                                   name = paste0("kNN_network_k", knn_k))
      knn_name <- paste0("kNN_network_k", knn_k)
      network_names <- c(network_names, knn_name)
      cat("\u2713 kNN network created\n")
      knn_net <- getSpatialNetwork(gobj, name = knn_name, output = "networkDT")
      write.csv(knn_net,
                file.path(net_folder,
                          paste0(sample_id, "_knn_k", knn_k, "_network.csv")),
                row.names = FALSE)
      cat("  \u2713 Saved: kNN network edges\n\n")
    }, error = function(e) cat("\u26A0 kNN failed:", conditionMessage(e), "\n\n"))
  }
  
  if (length(network_names) == 0)
    stop("No spatial networks were successfully created. Cannot continue.")
  
  primary_network <- if ("Delaunay_network" %in% network_names) {
    "Delaunay_network"
  } else { network_names[1] }
  
  cat("Primary network for downstream steps:", primary_network, "\n\n")


  # ============================================================================
  # SECTION 1b — NETWORK DIAGNOSTICS (edge-length + degree distributions)
  # ============================================================================

  tryCatch({
    diag_net <- getSpatialNetwork(gobj, name = primary_network,
                                  output = "networkDT")
    if (is.null(diag_net) || nrow(diag_net) == 0) {
      cat("\u26A0 No edges in ", primary_network, "; skipping diagnostics\n\n")
    } else {
      # Edge length
      if ("distance" %in% names(diag_net)) {
        elen <- as.numeric(diag_net$distance)
      } else {
        elen <- sqrt((diag_net$sdimx_end - diag_net$sdimx_begin)^2 +
                     (diag_net$sdimy_end - diag_net$sdimy_begin)^2)
      }
      elen_df <- data.frame(distance = elen)

      # Overall edge-length histogram
      p_elen <- ggplot2::ggplot(elen_df, ggplot2::aes(x = distance)) +
        ggplot2::geom_histogram(bins = 60, fill = "#4C72B0",
                                colour = "white", linewidth = 0.2) +
        ggplot2::labs(
          title    = sample_plot_title(sample_id,
                        paste0("Edge-length distribution - ", primary_network)),
          subtitle = sprintf("n = %d edges; median = %.1f; max = %.1f",
                             nrow(elen_df), stats::median(elen),
                             max(elen, na.rm = TRUE)),
          x = "Edge length (um)", y = "Edges"
        ) +
        presentation_theme(base_size = 12)
      save_presentation_plot(
        plot     = p_elen,
        filename = file.path(net_folder,
                             paste0(sample_id, "_edge_length_hist.png")),
        width    = 9, height = 6, dpi = 300
      )

      # Per-FOV edge-length histograms (subfolder)
      meta_fov <- as.data.frame(pDataDT(gobj))
      if ("fov" %in% names(meta_fov)) {
        fov_lookup <- setNames(meta_fov$fov, meta_fov$cell_ID)
        diag_net$fov_from <- fov_lookup[diag_net$from]
        diag_net$fov_to   <- fov_lookup[diag_net$to]
        diag_net$fov_same <- diag_net$fov_from == diag_net$fov_to
        within_fov <- diag_net[!is.na(diag_net$fov_same) & diag_net$fov_same, ,
                               drop = FALSE]
        if ("distance" %in% names(within_fov)) {
          within_fov$.elen <- as.numeric(within_fov$distance)
        } else {
          within_fov$.elen <- sqrt(
            (within_fov$sdimx_end - within_fov$sdimx_begin)^2 +
            (within_fov$sdimy_end - within_fov$sdimy_begin)^2
          )
        }
        fov_dir <- file.path(net_folder, "edge_length_per_fov")
        dir.create(fov_dir, recursive = TRUE, showWarnings = FALSE)
        fovs <- sort(unique(within_fov$fov_from))
        n_written <- 0L
        for (fv in fovs) {
          sub_df <- within_fov[within_fov$fov_from == fv, , drop = FALSE]
          if (nrow(sub_df) < 20) next
          p_fv <- ggplot2::ggplot(sub_df,
              ggplot2::aes(x = .elen)) +
            ggplot2::geom_histogram(bins = 50, fill = "#4C72B0",
                                    colour = "white", linewidth = 0.2) +
            ggplot2::labs(
              title    = sample_plot_title(sample_id,
                           paste0("Edge-length - FOV ", fv)),
              subtitle = sprintf("n = %d within-FOV edges; median = %.1f",
                                 nrow(sub_df), stats::median(sub_df$.elen)),
              x = "Edge length (um)", y = "Edges"
            ) +
            presentation_theme(base_size = 11)
          save_presentation_plot(
            plot     = p_fv,
            filename = file.path(fov_dir,
              paste0(sample_id, "_edge_length_fov_", fv, ".png")),
            width    = 7, height = 5, dpi = 200
          )
          n_written <- n_written + 1L
        }
        cat(sprintf("  \u2713 Per-FOV edge-length histograms: %d saved under %s/\n",
                    n_written, basename(fov_dir)))
      }

      # Degree distribution per cell type
      if (celltype_col %in% names(meta_fov)) {
        deg_tbl <- c(
          table(as.character(diag_net$from)),
          table(as.character(diag_net$to))
        )
        deg_df <- data.frame(
          cell_ID = names(deg_tbl),
          degree  = as.integer(deg_tbl),
          stringsAsFactors = FALSE
        )
        deg_df <- aggregate(degree ~ cell_ID, data = deg_df, FUN = sum)
        deg_df <- merge(deg_df,
                        meta_fov[, c("cell_ID", celltype_col)],
                        by = "cell_ID")
        names(deg_df)[3] <- "celltype"
        deg_df$celltype <- factor(deg_df$celltype)
        p_deg <- ggplot2::ggplot(deg_df,
            ggplot2::aes(x = celltype, y = degree, fill = celltype)) +
          ggplot2::geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4) +
          ggplot2::labs(
            title = sample_plot_title(sample_id,
                      paste0("Node degree by cell type - ", primary_network)),
            x = NULL, y = "Degree (edges per cell)"
          ) +
          presentation_theme(base_size = 11) +
          ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            legend.position = "none"
          )
        save_presentation_plot(
          plot     = p_deg,
          filename = file.path(net_folder,
                               paste0(sample_id, "_degree_by_celltype.png")),
          width    = max(10, 0.4 * length(levels(deg_df$celltype)) + 4),
          height   = 7, dpi = 300
        )
        cat("  \u2713 Degree-by-celltype plot saved\n")
      }
      cat("\n")
    }
  }, error = function(e) {
    cat("\u26A0 Network diagnostics failed:", conditionMessage(e), "\n\n")
  })


  # ============================================================================
  # SECTION 2 — NEIGHBOURHOOD ENRICHMENT (cellProximityEnrichment)
  # ============================================================================
  
  cat("================================================================================\n")
  cat("SECTION 2: Cell Type Neighbourhood Enrichment\n")
  cat("================================================================================\n\n")
  
  cat(sprintf("Running cellProximityEnrichment (%d permutations)...\n",
              n_simulations))
  cat("  Network:    ", primary_network, "\n")
  cat("  Cell types: ", celltype_col, "\n\n")
  
  prox_results <- tryCatch({
    cellProximityEnrichment(
      gobject               = gobj,
      cluster_column        = celltype_col,
      spatial_network_name  = primary_network,
      adjust_method         = "fdr",
      number_of_simulations = n_simulations
    )
  }, error = function(e) {
    cat("\u26A0 cellProximityEnrichment failed:", conditionMessage(e), "\n")
    NULL
  })
  
  prox_csv_path <- file.path(prox_folder,
                             paste0(sample_id, "_proximity_enrichment.csv"))
  
  if (!is.null(prox_results)) {
    
    # Save enrichment table
    if (!is.null(prox_results$enrichm_res)) {
      enr_out       <- as.data.frame(prox_results$enrichm_res)
      priority_cols <- c("unified_int", "type_int", "enrichm",
                         "p.adj_higher", "p.adj_lower", "PI_value",
                         "original", "simulations",
                         "p_higher_orig", "p_lower_orig", "int_ranking")
      col_order <- c(priority_cols[priority_cols %in% names(enr_out)],
                     setdiff(names(enr_out), priority_cols))
      write.csv(enr_out[, col_order], prox_csv_path, row.names = FALSE)
      cat("\u2713 Enrichment table saved\n")
    }
    
    if (!is.null(prox_results$result)) {
      write.csv(prox_results$result,
                file.path(prox_folder,
                          paste0(sample_id, "_proximity_enrichment_full.csv")),
                row.names = FALSE)
    }
    
    if (create_plots) {
      
      # Proximity heatmap via pheatmap
      tryCatch({
        .plot_proximity_heatmap(
          prox_csv_path = prox_csv_path,
          sample_id     = sample_id,
          prox_folder   = prox_folder
        )
      }, error = function(e) {
        cat("\u26A0 Proximity heatmap failed:", conditionMessage(e), "\n")
      })
      
      # Vis plot — direct ggplot2 (bypasses cellProximityVisPlot Giotto bug)
      tryCatch({
        if (!file.exists(prox_csv_path))
          stop("Proximity enrichment CSV not found")
        
        enr_csv <- data.table::fread(prox_csv_path, data.table = FALSE)
        if (!is.character(enr_csv$unified_int))
          enr_csv$unified_int <- as.character(enr_csv$unified_int)
        
        hetero_csv <- enr_csv[!is.na(enr_csv$type_int) &
                                enr_csv$type_int == "hetero", ]
        if (nrow(hetero_csv) == 0) stop("No hetero interactions in saved CSV")
        
        sig_h   <- hetero_csv[!is.na(hetero_csv$p.adj_higher) &
                                hetero_csv$p.adj_higher < 0.05, ]
        top_int <- if (nrow(sig_h) > 0) {
          sig_h$unified_int[which.max(sig_h$enrichm)]
        } else {
          hetero_csv$unified_int[which.max(hetero_csv$enrichm)]
        }
        
        cat("  Plotting top hetero interaction:", top_int, "\n")
        .plot_cell_proximity(
          gobj             = gobj,
          interaction_name = top_int,
          celltype_col     = celltype_col,
          primary_network  = primary_network,
          sample_id        = sample_id,
          prox_folder      = prox_folder
        )

        # Composite CART slides: emit per-sub-biopsy proximity vis plots
        # into prox_folder/subsamples/. Requires sample_row + sheet path
        # from the orchestrator; silently no-ops for non-composite samples.
        sub_rows <- tryCatch(
          discover_composite_subsamples(sample_row, sample_sheet_path),
          error = function(e) NULL
        )
        if (!is.null(sub_rows) && nrow(sub_rows) > 0) {
          meta_all <- as.data.frame(pDataDT(gobj))
          if ("fov" %in% names(meta_all)) {
            sub_dir <- file.path(prox_folder, "subsamples")
            dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
            cat("  Composite sample - rendering ", nrow(sub_rows),
                " per-sub-biopsy proximity variant(s)\n", sep = "")
            for (k in seq_len(nrow(sub_rows))) {
              sub_r  <- sub_rows[k, , drop = FALSE]
              sub_id <- as.character(sub_r$sample_id)
              fmin   <- as.integer(sub_r$fov_min)
              fmax   <- as.integer(sub_r$fov_max)
              if (anyNA(c(fmin, fmax))) next
              cell_ids <- meta_all$cell_ID[
                !is.na(meta_all$fov) & meta_all$fov >= fmin & meta_all$fov <= fmax
              ]
              if (length(cell_ids) == 0) next
              sub_gobj <- tryCatch(
                subsetGiotto(gobj, cell_ids = cell_ids),
                error = function(e) NULL
              )
              if (is.null(sub_gobj)) next
              tryCatch(
                .plot_cell_proximity(
                  gobj             = sub_gobj,
                  interaction_name = top_int,
                  celltype_col     = celltype_col,
                  primary_network  = primary_network,
                  sample_id        = sub_id,
                  prox_folder      = sub_dir,
                  file_label       = sub_id
                ),
                error = function(e) cat("    \u26A0 ", sub_id,
                  ": proximity sub-plot failed: ",
                  conditionMessage(e), "\n", sep = "")
              )
            }
          }
        }
      }, error = function(e) {
        cat("\u26A0 Proximity vis plot failed:", conditionMessage(e), "\n")
      })
    }
    
    # Console summary
    cat("\n=== Neighbourhood Enrichment Summary ===\n")
    enr_s  <- as.data.frame(prox_results$enrichm_res)
    hetero <- enr_s[!is.na(enr_s$type_int) & enr_s$type_int == "hetero", ]
    homo   <- enr_s[!is.na(enr_s$type_int) & enr_s$type_int == "homo", ]
    cat(sprintf("  Total interactions tested:               %d\n", nrow(enr_s)))
    cat(sprintf("    Homotypic  (self-self):                %d\n", nrow(homo)))
    cat(sprintf("    Heterotypic (cross-type):              %d\n", nrow(hetero)))
    
    if (file.exists(prox_csv_path)) {
      enr_csv    <- data.table::fread(prox_csv_path, data.table = FALSE)
      hetero_csv <- enr_csv[!is.na(enr_csv$type_int) &
                              enr_csv$type_int == "hetero", ]
      
      if ("p.adj_higher" %in% names(hetero_csv)) {
        sig_pos <- sum(hetero_csv$enrichm > 0 &
                         hetero_csv$p.adj_higher < 0.05, na.rm = TRUE)
        sig_neg <- sum(hetero_csv$enrichm < 0 &
                         hetero_csv$p.adj_lower  < 0.05, na.rm = TRUE)
        cat(sprintf(
          "  Sig. co-localised hetero pairs (p.adj_higher < 0.05): %d\n", sig_pos))
        cat(sprintf(
          "  Sig. exclusion    hetero pairs (p.adj_lower  < 0.05): %d\n", sig_neg))
        
        top_c <- hetero_csv[hetero_csv$enrichm > 0 &
                              !is.na(hetero_csv$p.adj_higher) &
                              hetero_csv$p.adj_higher < 0.05, ]
        top_c <- top_c[order(-top_c$enrichm), ]
        cat("\n  Top 5 co-localised hetero pairs:\n")
        for (i in seq_len(min(5, nrow(top_c))))
          cat(sprintf("    %s  enrichm=%.3f  p.adj_higher=%.2e\n",
                      top_c$unified_int[i], top_c$enrichm[i],
                      top_c$p.adj_higher[i]))
        
        top_e <- hetero_csv[hetero_csv$enrichm < 0 &
                              !is.na(hetero_csv$p.adj_lower) &
                              hetero_csv$p.adj_lower < 0.05, ]
        top_e <- top_e[order(top_e$enrichm), ]
        cat("\n  Top 5 exclusion hetero pairs:\n")
        for (i in seq_len(min(5, nrow(top_e))))
          cat(sprintf("    %s  enrichm=%.3f  p.adj_lower=%.2e\n",
                      top_e$unified_int[i], top_e$enrichm[i],
                      top_e$p.adj_lower[i]))
      }
    }
    cat("\n\u2713 Neighbourhood enrichment complete\n\n")
  }
  
  
  # ============================================================================ 
  # SECTION 3 — NICHE DECONVOLUTION 
  # ============================================================================ 
  
  cat("================================================================================\n")
  cat("SECTION 3: Niche Deconvolution (Spot-level Enrichment Scores)\n")
  cat("================================================================================\n\n")
  
  if ("PAGE" %in% enrich_methods) {
    if (length(marker_list_page) == 0) {
      cat("\u26A0 PAGE skipped: no cell types with >=", min_page_genes,
          "genes after filtering.\n\n")
    } else {
      cat(sprintf("Running PAGE enrichment (%d cell types)...\n",
                  length(marker_list_page)))
      tryCatch({
        gobj <- runPAGEEnrich(
          gobject           = gobj,
          sign_matrix       = makeSignMatrixPAGE(
            sign_names = names(marker_list_page),
            sign_list  = marker_list_page
          ),
          expression_values = expression_values,
          reverse_log_scale = FALSE,
          output_enrichment = "zscore",
          name              = "PAGE_enrichment"
        )
        cat("\u2713 PAGE enrichment complete\n")
        
        page_scores <- getSpatialEnrichment(gobj, name = "PAGE_enrichment",
                                            output = "data.table")
        score_cols  <- setdiff(names(page_scores), "cell_ID")
        cat(sprintf("  \u2713 Scores saved (%d cell types scored)\n",
                    length(score_cols)))
        write.csv(page_scores,
                  file.path(deconv_folder,
                            paste0(sample_id, "_PAGE_enrichment_scores.csv")),
                  row.names = FALSE);
        
        if (create_plots && length(score_cols) > 0) {
          cat(sprintf("  Saving spatial plots for all %d scored cell types...\n",
                      length(score_cols)))
          for (ct in score_cols) {
            tryCatch({
              saved <- .plot_deconv_polygons(
                gobj            = gobj,
                score_col       = ct,
                enrichment_name = "PAGE_enrichment",
                sample_id       = sample_id,
                deconv_folder   = deconv_folder,
                legend_name     = paste(pretty_plot_label(ct), "PAGE score")
              )
              if (is.null(saved)) {
                # Fallback: Giotto spatPlot2D when polygon data is missing.
                spatPlot2D(
                  gobject         = gobj,
                  cell_color      = ct,
                  color_as_factor = FALSE,
                  spat_enr_names  = "PAGE_enrichment",
                  point_size      = 0.4,
                  show_image      = FALSE,
                  save_plot       = TRUE,
                  save_param      = list(
                    save_name   = paste0(sample_id, "_PAGE_",
                                         gsub("[^A-Za-z0-9]", "_", ct)),
                    save_dir    = deconv_folder,
                    base_width  = 10,
                    base_height = 8
                  )
                )
              }
            }, error = function(e) cat("  ⚠ PAGE plot failed for", ct, ":", conditionMessage(e), "\n"))
          }
          cat("  \u2713 Spatial enrichment plots saved\n")
        }
      }, error = function(e) {
        cat("\u26A0 PAGE enrichment failed:", conditionMessage(e), "\n\n")
      })
    }
  }
  
  if ("rank" %in% enrich_methods) {
    cat(sprintf("Running rank enrichment (%d cell types)...\n", n_types))
    tryCatch({
      expr_mat  <- getExpression(gobj, values = expression_values,
                                 output = "matrix")
      rank_sign <- makeSignMatrixRank(
        sc_matrix      = expr_mat,
        sc_cluster_ids = pDataDT(gobj)[[celltype_col]],
        ties_method    = "max"
      )
      gobj <- runRankEnrich(gobject = gobj, sign_matrix = rank_sign,
                            expression_values = expression_values,
                            name = "rank_enrichment")
      cat("\u2713 Rank enrichment complete\n")
      
      rank_scores  <- getSpatialEnrichment(gobj, name = "rank_enrichment",
                                           output = "data.table")
      score_cols_r <- setdiff(names(rank_scores), "cell_ID")
      write.csv(rank_scores,
                file.path(deconv_folder,
                          paste0(sample_id, "_rank_enrichment_scores.csv")),
                row.names = FALSE)
      cat(sprintf("  \u2713 Scores saved (%d cell types scored)\n\n",
                  length(score_cols_r)))
    }, error = function(e) {
      cat("\u26A0 Rank enrichment failed:", conditionMessage(e), "\n\n")
    })
  }
  
  if (create_plots) {
    page_present <- !is.null(
      tryCatch(getSpatialEnrichment(gobj, name = "PAGE_enrichment"),
               error = function(e) NULL)
    )
    if (page_present) {
      tryCatch({
        page_scores <- getSpatialEnrichment(gobj, name = "PAGE_enrichment",
                                            output = "data.table")
        score_cols  <- setdiff(names(page_scores), "cell_ID")
        if (length(score_cols) >= 2) {
          plotMetaDataCellsHeatmap(
            gobject        = gobj,
            metadata_cols  = celltype_col,
            value_cols     = score_cols,
            spat_enr_names = "PAGE_enrichment",
            save_plot      = TRUE,
            save_param     = list(
              save_name   = paste0(sample_id, "_PAGE_heatmap_by_celltype"),
              save_dir    = deconv_folder,
              base_width  = 14,
              base_height = 10
            )
          )
          cat("\u2713 PAGE score heatmap by cell type saved\n")
        } else {
          cat("\u26A0 PAGE heatmap skipped: fewer than 2 cell types scored\n")
        }
      }, error = function(e) {
        cat("\u26A0 PAGE heatmap failed:", conditionMessage(e), "\n")
      })
    }
  }
  
  
  # ============================================================================ 
  # SECTION 4 — SAVE & SUMMARISE 
  # ============================================================================ 
  
  cat("================================================================================\n")
  cat("SECTION 4: Save & Summary\n")
  cat("================================================================================\n\n")
  
  tryCatch({
    if (exists("save_giotto_checkpoint")) {
      checkpoint_info <- save_giotto_checkpoint(
        gobj = gobj,
        checkpoint_dir = file.path(output_dir, "Giotto_Object_Spatial"),
        overwrite = TRUE,
        metadata = list(
          sample_id = sample_id,
          pipeline_step = "09_spatial_network"
        )
      )
      if (!identical(checkpoint_info$save_method, "none")) {
        cat("\u2713 Saved: Giotto_Object_Spatial (", checkpoint_info$save_method, ")\n", sep = "")
      } else {
        cat("\u26A0 Save skipped: no serializer succeeded\n")
      }
      if (!is.null(checkpoint_info$error_message) && nzchar(checkpoint_info$error_message)) {
        cat("  saveGiotto fallback:", checkpoint_info$error_message, "\n")
      }
    } else {
      saveGiotto(
        gobject = gobj,
        dir = output_dir,
        foldername = "Giotto_Object_Spatial",
        overwrite = TRUE
      )
      cat("\u2713 Saved: Giotto_Object_Spatial\n")
    }
  }, error = function(e) cat("\u26A0 Save failed:", conditionMessage(e), "\n"))
  
  cat("\n=== Step 09 Summary ===\n")
  cat("  Spatial networks built:  ", length(network_names), "\n")
  cat("  Primary network:         ", primary_network, "\n")
  cat("  Cell types analysed:     ", n_types, "\n")
  cat("  Enrichment methods run:  ",
      paste(enrich_methods[enrich_methods %in% c("PAGE", "rank")],
            collapse = ", "), "\n")
  cat("  Output folder:           ", base_folder, "\n")
  cat("\n\u2713 STEP 09 complete for", sample_id, "\n\n")
  
  invisible(gobj)
}


# Run if sourced directly --------------------------------------------------- 
if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 3) {
    bootstrap_script <- file.path(current_script_dir(), "Helper_Scripts", "Script_Bootstrap.R")
    if (file.exists(bootstrap_script)) {
      source(bootstrap_script, local = .GlobalEnv)
      bootstrap_pipeline_environment(current_script_dir(), load_pipeline_utils = TRUE, verbose = FALSE)
    }
    
    # Composite-awareness: pick up the sample sheet path from env var so a
    # direct Rscript invocation can emit per-sub-biopsy proximity variants.
    sheet_env <- Sys.getenv("COSMX_SAMPLE_SHEET", unset = "")
    sample_sheet_path_arg <- if (nzchar(sheet_env) && file.exists(sheet_env)) sheet_env else NULL
    sample_row_arg <- NULL
    if (!is.null(sample_sheet_path_arg) &&
        exists("safe_read_sheet", mode = "function", inherits = TRUE)) {
      sheet_df <- tryCatch(safe_read_sheet(sample_sheet_path_arg),
                           error = function(e) NULL)
      if (!is.null(sheet_df) && "sample_id" %in% names(sheet_df)) {
        hit <- which(as.character(sheet_df$sample_id) == args[1])
        if (length(hit) > 0) sample_row_arg <- sheet_df[hit[1], , drop = FALSE]
      }
    }

    build_spatial_network(
      gobj              = args[2],
      sample_id         = args[1],
      output_dir        = args[3],
      sample_row        = sample_row_arg,
      sample_sheet_path = sample_sheet_path_arg
    )
  } else {
    stop("Usage: Rscript 09_Spatial_Network.R <sample_id> <input_path> <output_dir>")
  }
}

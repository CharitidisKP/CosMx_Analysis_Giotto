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
if (!exists("save_giotto_checkpoint") && file.exists(pipeline_utils)) {
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
                                 prox_folder) {

  ann_net <- annotateSpatialNetwork(
    gobject              = gobj,
    spatial_network_name = primary_network,
    cluster_column       = celltype_col
  )

  spat_locs <- getSpatialLocations(gobj, output = "data.table")
  meta      <- pDataDT(gobj)

  cell_dt <- merge(
    spat_locs[, c("cell_ID", "sdimx", "sdimy")],
    meta[,     c("cell_ID", celltype_col), with = FALSE],
    by = "cell_ID"
  )
  names(cell_dt)[names(cell_dt) == celltype_col] <- "cell_type"

  int_edges    <- ann_net[ann_net$unified_int == interaction_name, ]
  ct_pair      <- strsplit(interaction_name, "--")[[1]]
  cell_dt$in_pair <- cell_dt$cell_type %in% ct_pair

  edge_dt     <- int_edges[, c("sdimx_begin", "sdimy_begin",
                               "sdimx_end",   "sdimy_end")]
  selected_dt <- cell_dt[cell_dt$in_pair, ]
  other_dt    <- cell_dt[!cell_dt$in_pair, ]

  ct_colours <- setNames(
    scales::hue_pal()(length(ct_pair)),
    ct_pair
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data    = other_dt,
      mapping = ggplot2::aes(x = sdimx, y = sdimy),
      colour  = "grey85", size = 0.25, alpha = 0.35
    ) +
    ggplot2::geom_segment(
      data    = edge_dt,
      mapping = ggplot2::aes(x    = sdimx_begin, y    = sdimy_begin,
                             xend = sdimx_end,   yend = sdimy_end),
      colour = "grey40", linewidth = 0.2, alpha = 0.5
    ) +
    ggplot2::geom_point(
      data    = selected_dt,
      mapping = ggplot2::aes(x = sdimx, y = sdimy, colour = cell_type),
      size    = 1.4, alpha = 0.9
    ) +
    ggplot2::scale_colour_manual(values = ct_colours, name = "Cell type") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title    = paste0("Cell proximity: ", interaction_name),
      subtitle = sprintf("%d interacting edges  |  network: %s",
                         nrow(int_edges), primary_network),
      x = "x (global pixels)", y = "y (global pixels)"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(plot.title      = ggplot2::element_text(face = "bold"),
                   legend.position = "right")

  safe_name <- gsub("[^A-Za-z0-9]", "_", interaction_name)
  save_path <- file.path(prox_folder,
                         paste0(sample_id, "_proximity_visplot_",
                                safe_name, ".png"))
  ggplot2::ggsave(save_path, plot = p, width = 14, height = 10, dpi = 150)
  cat("\u2713 Proximity vis plot saved:", basename(save_path), "\n")
  invisible(p)
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

  breaks <- seq(-max_abs, max_abs, length.out = 101)
  colors <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

  # Row/col annotation: immune vs structural compartment
  immune_types <- c("B.cell", "CD4.T.cell", "CD8.T.cell", "NK.cell",
                    "NKT.cell",
                    "MNP.macrophage", "MNP.monocyte",
                    "MNP.a.classical.monocyte.derived",
                    "MNP.b.non.classical.monocyte.derived",
                    "MNP.c.dendritic.cell", "MNP.d.Tissue.macrophage",
                    "Plasmacytoid.DC", "Plasmacytoid.dendritic.cell",
                    "Mast.cell", "Neutrophil")

  anno_df <- data.frame(
    Compartment = ifelse(all_types %in% immune_types, "Immune", "Structural"),
    row.names   = all_types
  )
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
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
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
                                  celltype_col      = NULL,
                                  marker_list       = NULL,
                                  network_methods   = c("Delaunay", "kNN"),
                                  knn_k             = 6,
                                  n_simulations     = 1000,
                                  enrich_methods    = c("PAGE", "rank"),
                                  expression_values = "normalized",
                                  create_plots      = TRUE,
                                  random_seed       = 42,
                                  min_page_genes    = 5) {

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
    cat("Building Delaunay spatial network...\n")
    tryCatch({
      gobj <- createSpatialNetwork(gobject = gobj, method = "Delaunay",
                                   minimum_k = 2, name = "Delaunay_network")
      network_names <- c(network_names, "Delaunay_network")
      cat("\u2713 Delaunay network created\n")
      del_net <- getSpatialNetwork(gobj, name = "Delaunay_network",
                                   output = "networkDT")
      cat(sprintf("  Edges: %d  |  Unique nodes: %d\n",
                  nrow(del_net),
                  length(unique(c(del_net$from, del_net$to)))))
      write.csv(del_net,
                file.path(net_folder, paste0(sample_id, "_delaunay_network.csv")),
                row.names = FALSE)
      cat("  \u2713 Saved: Delaunay network edges\n\n")
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
            }, error = function(e) NULL)
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
    saveGiotto(gobject = gobj, dir = output_dir,
               foldername = "Giotto_Object_Spatial", overwrite = TRUE)
    cat("\u2713 Saved: Giotto_Object_Spatial\n")
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
    build_spatial_network(
      gobj       = args[2],
      sample_id  = args[1],
      output_dir = args[3]
    )
  }
}
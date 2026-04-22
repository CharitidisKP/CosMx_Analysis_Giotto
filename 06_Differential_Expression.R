# Cluster marker analysis -----------------------------------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 06_Differential_Expression.R  (pipeline alias: 06_markers)
#
# NOTE ON NAMING
#   The historical file name is "06_Differential_Expression.R" but the
#   analysis performed here is *cluster marker discovery* (one-vs-all per
#   Leiden cluster), NOT cross-sample / cross-treatment differential
#   expression. True inter-sample spatial DE is implemented in
#   12_Spatial_Differential_Expression.R.
#
#   The pipeline dispatches this script via the alias "06_markers"
#   (see CosMx_pipeline.R::canonicalize_step_id). The exported function is
#   marker_analysis().
# ==============================================================================

#' Marker Gene Analysis
#'
#' Per-cluster marker discovery (one-vs-all) using scran. This is NOT
#' cross-sample differential expression — see script 12 for that.
#'
#' @param gobj Giotto object or path
#' @param sample_id Sample identifier
#' @param output_dir Output directory
#' @param cluster_column Clustering column name
#' @param top_n Top N markers per cluster
#' @return Giotto object with marker results

# Custom cluster-marker heatmap.
# Replaces Giotto::plotMetaDataHeatmap() with a diverging blue/white/red z-score
# tile heatmap. Genes are row-z-scored across clusters and clipped to [-2, 2]
# so outlier clusters do not dominate the colour scale. Gene order on the y-axis
# is preserved from the input vector (so callers control grouping — main markers
# before subtype markers for the B-cell variant). Panel-absent genes drop
# silently via an intersect with the expression matrix.
.cluster_marker_heatmap <- function(gobj, sample_id, genes, cluster_column,
                                    results_folder,
                                    title_suffix    = "Cluster marker heatmap",
                                    filename_suffix = "marker_heatmap") {
  expr <- getExpression(gobj, values = "normalized", output = "matrix")
  meta <- as.data.frame(pDataDT(gobj))

  genes_in   <- unique(as.character(genes))
  genes_in   <- genes_in[nzchar(genes_in)]
  genes_hit  <- intersect(genes_in, rownames(expr))
  genes_skip <- setdiff(genes_in, genes_hit)
  if (length(genes_skip) > 0) {
    cat(sprintf("  ℹ %d gene(s) absent from expression matrix: %s\n",
                length(genes_skip), paste(genes_skip, collapse = ", ")))
  }
  if (length(genes_hit) == 0) {
    cat("  ℹ no genes available for heatmap; skipping\n")
    return(invisible(NULL))
  }

  cluster_levels <- sort(unique(as.character(meta[[cluster_column]])))
  n_clusters <- length(cluster_levels)

  mean_mat <- matrix(NA_real_,
                     nrow = length(genes_hit),
                     ncol = n_clusters,
                     dimnames = list(genes_hit, cluster_levels))
  for (cl in cluster_levels) {
    keep <- meta[[cluster_column]] == cl
    cell_ids <- meta$cell_ID[keep]
    sub <- expr[genes_hit, cell_ids, drop = FALSE]
    mean_mat[, cl] <- Matrix::rowMeans(sub)
  }

  z_mat <- t(scale(t(mean_mat), center = TRUE, scale = TRUE))
  z_mat[!is.finite(z_mat)] <- 0

  df <- data.frame(
    gene    = factor(rep(genes_hit, times = n_clusters), levels = genes_hit),
    cluster = factor(rep(cluster_levels, each = length(genes_hit)),
                     levels = cluster_levels),
    z       = as.vector(z_mat),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df,
        ggplot2::aes(x = cluster, y = gene, fill = z)) +
    ggplot2::geom_tile(colour = "grey90", linewidth = 0.1) +
    ggplot2::scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-2, 2),
      oob      = scales::squish,
      name     = "z-score"
    ) +
    ggplot2::scale_y_discrete(limits = rev(genes_hit)) +
    ggplot2::labs(
      title = sample_plot_title(sample_id, title_suffix),
      x     = "Cluster",
      y     = "Gene"
    ) +
    presentation_theme(base_size = 11) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 8)
    )

  out_path <- file.path(results_folder,
                        paste0(sample_id, "_", filename_suffix, ".png"))
  save_presentation_plot(
    plot     = p,
    filename = out_path,
    width    = max(6, 0.30 * n_clusters + 3),
    height   = max(8, 0.22 * length(genes_hit) + 3),
    dpi      = 300
  )

  cat("✓ Heatmap saved:", basename(out_path), "\n\n")
  invisible(out_path)
}


marker_analysis <- function(gobj,
                            sample_id,
                            output_dir,
                            cluster_column  = "leiden_clust",
                            top_n           = 25,
                            bcell_markers   = character(),
                            subtype_markers = character()) {
  
  cat("\n========================================\n")
  cat("STEP 06: Cluster Marker Analysis (one-vs-all)\n")
  cat("Sample:", sample_id, "\n")
  cat("Note:  NOT cross-sample DE; see step 12 for that.\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }
  
  results_folder <- file.path(output_dir, "06_Markers")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  cat("Finding marker genes...\n")
  cat("  Cluster column:", cluster_column, "\n")
  cat("  Top N per cluster:", top_n, "\n\n")
  
  # Find markers using scran
  markers_scran <- findMarkers_one_vs_all(
    gobject = gobj,
    method = "scran",
    expression_values = "normalized",
    cluster_column = cluster_column,
    min_feats = 3
  )
  
  cat("✓ Markers found\n\n")
  
  # Extract top markers per cluster
  top_markers <- markers_scran %>%
    as_tibble() %>%
    group_by(cluster) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  # Save all markers
  write_csv(markers_scran, 
            file.path(results_folder, paste0(sample_id, "_all_markers.csv")))
  
  # Save top markers
  write_csv(top_markers,
            file.path(results_folder, paste0(sample_id, "_top_markers.csv")))
  
  cat("✓ Marker tables saved\n\n")
  
  # Create heatmap of top markers (custom ggplot tile heatmap; replaces Giotto's
  # plotMetaDataHeatmap() for visual consistency with the rest of the pipeline).
  cat("Creating marker heatmap...\n")

  top_genes <- top_markers %>%
    group_by(cluster) %>%
    slice_head(n = 10) %>%
    pull(feats)
  top_genes <- unique(top_genes)

  .cluster_marker_heatmap(
    gobj            = gobj,
    sample_id       = sample_id,
    genes           = top_genes,
    cluster_column  = cluster_column,
    results_folder  = results_folder,
    title_suffix    = "Cluster marker heatmap",
    filename_suffix = "marker_heatmap"
  )

  # B-cell-specific heatmap: main B-cell markers first, then subtype markers.
  # Skipped when no markers are provided via config / CLI env vars.
  b_all <- unique(c(as.character(bcell_markers), as.character(subtype_markers)))
  b_all <- b_all[nzchar(b_all)]
  if (length(b_all) > 0) {
    cat("Creating B-cell gene heatmap...\n")
    .cluster_marker_heatmap(
      gobj            = gobj,
      sample_id       = sample_id,
      genes           = b_all,
      cluster_column  = cluster_column,
      results_folder  = results_folder,
      title_suffix    = "B-cell gene heatmap",
      filename_suffix = "bcell_heatmap"
    )
  }
  
  # Violin plots for top markers
  cat("Creating violin plots for top 5 markers per cluster...\n")
  
  top5_per_cluster <- top_markers %>%
    group_by(cluster) %>%
    slice_head(n = 5) %>%
    pull(feats) %>%
    unique()
  
  p_vio <- violinPlot(
    gobject        = gobj,
    feats          = head(top5_per_cluster, 20),  # Limit to avoid overcrowding
    cluster_column = cluster_column,
    strip_text     = 8,
    strip_position = "right",
    show_plot      = FALSE,
    return_plot    = TRUE
  ) + ggplot2::labs(title = sample_plot_title(sample_id, "Top marker violins"))

  save_presentation_plot(
    plot     = p_vio,
    filename = file.path(results_folder,
                         paste0(sample_id, "_marker_violins.png")),
    width    = 16,
    height   = 20,
    dpi      = 300
  )

  cat("✓ Violin plots saved\n\n")

  # B-cell-genes violin variant (panel-present subset; main + subtype markers).
  if (length(b_all) > 0) {
    expr_rows <- rownames(getExpression(gobj, values = "normalized",
                                        output = "matrix"))
    b_present <- intersect(b_all, expr_rows)
    if (length(b_present) > 0) {
      cat("Creating B-cell gene violin plots...\n")
      p_vio_b <- violinPlot(
        gobject        = gobj,
        feats          = b_present,
        cluster_column = cluster_column,
        strip_text     = 8,
        strip_position = "right",
        show_plot      = FALSE,
        return_plot    = TRUE
      ) + ggplot2::labs(title = sample_plot_title(sample_id,
                                                  "B-cell gene violins"))
      save_presentation_plot(
        plot     = p_vio_b,
        filename = file.path(results_folder,
                             paste0(sample_id, "_bcell_violins.png")),
        width    = 14,
        height   = max(8, 1.0 * length(b_present) + 2),
        dpi      = 300
      )
      cat("✓ B-cell violin plots saved\n\n")
    } else {
      cat("  ℹ no B-cell markers present on panel; skipping B-cell violins\n\n")
    }
  }

  # Per-cluster volcano plots -------------------------------------------------
  cat("Creating per-cluster volcano plots...\n")
  tryCatch({
    m_df <- as.data.frame(markers_scran)
    lfc_col <- intersect(c("summary.logFC", "summary_logFC", "logFC",
                           "median.logFC.other", "mean.logFC.other"),
                         names(m_df))[1]
    sig_col <- intersect(c("FDR", "padj", "p.adj", "p.value", "pvalue"),
                         names(m_df))[1]
    if (is.na(lfc_col) || is.na(sig_col)) {
      cat("  note: could not locate logFC / significance columns; skipping volcano\n")
      cat("  available columns:", paste(names(m_df), collapse = ", "), "\n")
    } else {
      volcano_dir <- file.path(results_folder, "volcano_plots")
      dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
      m_df$.lfc <- as.numeric(m_df[[lfc_col]])
      m_df$.sig <- as.numeric(m_df[[sig_col]])
      m_df$.nlog10 <- -log10(pmax(m_df$.sig, 1e-300))
      m_df$.significant <- m_df$.sig < 0.05 & abs(m_df$.lfc) > 0.25
      n_written <- 0L
      for (cl in sort(unique(m_df$cluster))) {
        df_cl <- m_df[m_df$cluster == cl, , drop = FALSE]
        if (nrow(df_cl) < 5) next
        top_labels <- df_cl[order(df_cl$.sig)[seq_len(min(10, nrow(df_cl)))], ]
        p <- ggplot2::ggplot(df_cl,
            ggplot2::aes(x = .lfc, y = .nlog10, colour = .significant)) +
          ggplot2::geom_point(alpha = 0.7, size = 1.2) +
          ggplot2::scale_colour_manual(
            values = c(`TRUE` = "#C44E52", `FALSE` = "grey70"),
            guide = "none"
          ) +
          ggrepel::geom_text_repel(
            data = top_labels,
            mapping = ggplot2::aes(label = feats),
            size = 3, max.overlaps = 12, colour = "black"
          ) +
          ggplot2::geom_hline(yintercept = -log10(0.05),
                              linetype = "dashed", colour = "grey40") +
          ggplot2::geom_vline(xintercept = c(-0.25, 0.25),
                              linetype = "dashed", colour = "grey40") +
          ggplot2::labs(
            title    = sample_plot_title(sample_id,
                                         paste0("Volcano \u2014 cluster ", cl)),
            subtitle = sprintf("%s vs. %s; %d genes",
                               lfc_col, sig_col, nrow(df_cl)),
            x = paste0("log fold-change (", lfc_col, ")"),
            y = paste0("-log10(", sig_col, ")")
          ) +
          presentation_theme(base_size = 11)
        save_presentation_plot(
          plot     = p,
          filename = file.path(volcano_dir,
                               paste0(sample_id, "_volcano_cluster_", cl, ".png")),
          width    = 8, height = 6, dpi = 300
        )
        n_written <- n_written + 1L
      }
      cat(sprintf("  \u2713 %d volcano plots saved under %s/\n",
                  n_written, basename(volcano_dir)))
    }
  }, error = function(e) {
    cat("\u26A0 Volcano plots failed:", conditionMessage(e), "\n")
  })

  # Top-markers dotplot -------------------------------------------------------
  # Axes: clusters on x, genes on y so "many markers" grows height and
  # clusters pack tight on x. The Giotto dotPlot() branch is dropped because
  # it offers no axis-swap hook; the manual ggplot is now the only path.
  cat("Creating marker dotplot...\n")
  tryCatch({
    dp_genes <- unique(top_genes)
    expr     <- getExpression(gobj, values = "normalized", output = "matrix")
    meta     <- as.data.frame(pDataDT(gobj))
    dp_genes <- intersect(dp_genes, rownames(expr))
    if (length(dp_genes) == 0) {
      cat("  \u26A0 No top-marker genes available; skipping dotplot\n")
    } else {
      cluster_levels <- sort(unique(as.character(meta[[cluster_column]])))
      rows <- list()
      for (g in dp_genes) {
        e <- expr[g, meta$cell_ID]
        for (cl in cluster_levels) {
          keep <- meta[[cluster_column]] == cl
          rows[[length(rows) + 1L]] <- data.frame(
            gene = g, cluster = as.character(cl),
            mean_expr = mean(e[keep]),
            pct_expr  = 100 * mean(e[keep] > 0),
            stringsAsFactors = FALSE
          )
        }
      }
      dot_df         <- do.call(rbind, rows)
      dot_df$gene    <- factor(dot_df$gene,    levels = dp_genes)
      dot_df$cluster <- factor(dot_df$cluster, levels = cluster_levels)
      p_dot <- ggplot2::ggplot(dot_df,
          ggplot2::aes(x = cluster, y = gene, size = pct_expr, colour = mean_expr)) +
        ggplot2::geom_point() +
        ggplot2::scale_colour_gradient(low = "lightgrey", high = "red",
                                       name = "mean expr") +
        ggplot2::scale_size_continuous(range = c(0, 6), name = "% cells") +
        ggplot2::scale_y_discrete(limits = rev(dp_genes)) +
        ggplot2::labs(
          title = sample_plot_title(sample_id, "Top marker dotplot"),
          x = "Cluster", y = "Gene"
        ) +
        presentation_theme(base_size = 11)
      save_presentation_plot(
        plot     = p_dot,
        filename = file.path(results_folder,
                             paste0(sample_id, "_marker_dotplot.png")),
        width    = max(6, 0.35 * length(cluster_levels) + 2),
        height   = max(8, 0.30 * length(dp_genes) + 3),
        dpi      = 300
      )
      cat("  \u2713 Dotplot saved\n")
    }
  }, error = function(e) {
    cat("\u26A0 Dotplot failed:", conditionMessage(e), "\n")
  })

  cat("=== Marker Analysis Summary ===\n")
  cat("Total marker genes:", nrow(markers_scran), "\n")
  cat("Top markers per cluster:", top_n, "\n")
  cat("Clusters analyzed:", length(unique(markers_scran$cluster)), "\n")
  
  cat("\n✓ Marker analysis complete for", sample_id, "\n\n")
  
  return(gobj)
}

# Run if sourced directly
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
    
    sample_id  <- args[1]
    input_path <- args[2]
    output_dir <- args[3]

    # B-cell marker lists can be passed via env vars so the CLI entry point
    # works the same way whether it's invoked by the orchestrator (which
    # reads Parameters/config.yaml) or by hand via Rscript. Unset / empty →
    # character(); the B-cell heatmap and violin blocks skip gracefully.
    .parse_gene_env <- function(name) {
      raw <- Sys.getenv(name, unset = "")
      if (!nzchar(raw)) return(character())
      parts <- unlist(strsplit(raw, "[,;]\\s*"))
      parts[nzchar(parts)]
    }
    bcell_markers   <- .parse_gene_env("COSMX_BCELL_MARKERS")
    subtype_markers <- .parse_gene_env("COSMX_BCELL_SUBTYPE_MARKERS")

    gobj <- marker_analysis(
      gobj            = input_path,
      sample_id       = sample_id,
      output_dir      = output_dir,
      bcell_markers   = bcell_markers,
      subtype_markers = subtype_markers
    )

    saveGiotto(gobj, dir = output_dir, foldername = "Giotto_Object_DEG_Markers", overwrite = TRUE)
  } else {
    stop("Usage: Rscript 06_Differential_Expression.R <sample_id> <input_path> <output_dir>")
  }
}

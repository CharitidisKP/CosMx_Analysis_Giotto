# Filter cells and genes based on QC metrics ------------------------------

#!/usr/bin/env Rscript
# ==============================================================================
# 02_quality_control.R
# Quality control and filtering of CosMx data
# ==============================================================================

#' Quality Control for CosMx Sample
#'
#' @param gobj Giotto object (or path to saved object)
#' @param sample_id Sample identifier
#' @param output_dir Output directory for this sample
#' @param gene_min_cells Minimum cells expressing a gene
#' @param cell_min_genes Minimum genes per cell
#' @param cell_max_genes Maximum genes per cell (NULL = no max)
#' @param min_count Minimum total counts per cell
#' @param max_mito_pct Maximum mitochondrial percentage (NULL = no filter)
#' @return Filtered Giotto object

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
if ((!exists("presentation_theme") || !exists("sample_plot_title") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
  source(pipeline_utils)
}

quality_control <- function(gobj,
                            sample_id,
                            output_dir,
                            gene_min_cells = 10,
                            cell_min_genes = 50,
                            cell_max_genes = NULL,
                            min_count = 100,
                            max_mito_pct = NULL,
                            sample_row = NULL,
                            sample_sheet_path = NULL) {
  
  cat("\n========================================\n")
  cat("STEP 02: Quality Control\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  # Load Giotto object if path provided
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("✓ Loaded\n\n")
  }

  # Create output directory
  results_folder <- file.path(output_dir, "02_QC")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

  # Drop SystemControl probes before any stats or filtering. These are QC
  # probes, not real genes - keeping them pollutes downstream InSituType gene
  # lists, dim reduction, and feature plots. Negative/FalseCode probes are
  # retained: InSituType uses them for per-cell background estimation.
  all_feat_ids_pre_ctrl <- fDataDT(gobj)$feat_ID
  ctrl_mask <- grepl("^SystemControl", all_feat_ids_pre_ctrl, ignore.case = TRUE)
  n_controls_removed <- sum(ctrl_mask)
  if (n_controls_removed > 0) {
    cat("Removing", n_controls_removed, "SystemControl probes (Negative/FalseCode retained for InSituType background)...\n")
    gobj <- subsetGiotto(
      gobject  = gobj,
      feat_ids = all_feat_ids_pre_ctrl[!ctrl_mask]
    )
    cat("✓ SystemControl probes removed\n\n")
  }

  # Get initial stats (post-SystemControl removal)
  n_cells_initial <- length(gobj@cell_ID$cell)
  n_genes_initial <- length(gobj@feat_ID$rna)
  
  cat("=== Initial Statistics ===\n")
  cat("Cells:", n_cells_initial, "\n")
  cat("Genes:", n_genes_initial, "\n\n")
  
  # Add statistics if not already present
  cat("Calculating cell/gene statistics...\n")
  
  gobj <- addStatistics(
    gobject = gobj,
    expression_values = "raw"
  )
  
  cat("✓ Statistics calculated\n\n")
  
  # Calculate mitochondrial percentage if MT genes present
  cat("Checking for mitochondrial genes...\n")
  
  all_feat_ids <- fDataDT(gobj)$feat_ID
  mt_genes <- grep("^MT-|^Mt-|^mt-", all_feat_ids, value = TRUE)

  if (length(mt_genes) > 0) {
    cat("Found", length(mt_genes), "mitochondrial genes\n")
    cat("  Examples:", paste(head(mt_genes, 5), collapse = ", "), "\n")

    gobj <- addFeatsPerc(
      gobject = gobj,
      feats = mt_genes,
      vector_name = "mito_pct"
    )

    cat("✓ Mitochondrial percentage calculated\n\n")
  } else {
    # Louder diagnostic: 0 matches usually indicates a non-human panel or a
    # renamed feature prefix (e.g. "mt-Nd1" vs "MT-ND1"), not actual absence
    # of mito genes. Show the first 5 feature IDs so the user can diagnose.
    warning(
      sprintf(
        "No mitochondrial genes matched '^MT-|^Mt-|^mt-' in %d features. Sample '%s' will proceed WITHOUT mito-pct filtering. First 5 feature IDs: %s",
        length(all_feat_ids),
        sample_id,
        paste(head(all_feat_ids, 5), collapse = ", ")
      ),
      immediate. = TRUE, call. = FALSE
    )
    cat("  (Adjust the MT regex in 02_Quality_Control.R if your panel uses a different convention.)\n\n")
  }
  
  # Get cell metadata for QC plots
  metadata <- pDataDT(gobj) %>%
    as_tibble() %>%
    mutate(total_expr_plot = pmax(total_expr, 1))
  
  cat("=== QC Metrics Summary ===\n")
  cat("Total counts per cell:\n")
  cat("  Min:", min(metadata$total_expr, na.rm = TRUE), "\n")
  cat("  Median:", median(metadata$total_expr, na.rm = TRUE), "\n")
  cat("  Max:", max(metadata$total_expr, na.rm = TRUE), "\n")
  cat("\nGenes per cell:\n")
  cat("  Min:", min(metadata$nr_feats, na.rm = TRUE), "\n")
  cat("  Median:", median(metadata$nr_feats, na.rm = TRUE), "\n")
  cat("  Max:", max(metadata$nr_feats, na.rm = TRUE), "\n")
  
  if ("mito_pct" %in% names(metadata)) {
    cat("\nMitochondrial %:\n")
    cat("  Median:", round(median(metadata$mito_pct, na.rm = TRUE), 2), "%\n")
    cat("  95th percentile:", round(quantile(metadata$mito_pct, 0.95, na.rm = TRUE), 2), "%\n")
  }
  cat("\n")
  
  # Create QC plots before filtering
  cat("Creating pre-filtering QC plots...\n")
  
  gene_median_pre  <- stats::median(metadata$nr_feats,        na.rm = TRUE)
  count_median_pre <- stats::median(metadata$total_expr_plot, na.rm = TRUE)

  # Pre-filter distribution plots: black dotted threshold + red dotted median,
  # same line style, different colour. No mean line.
  p1 <- ggplot(metadata, aes(x = nr_feats)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    geom_vline(xintercept = cell_min_genes,  color = "black", linetype = "dotted", linewidth = 0.9) +
    geom_vline(xintercept = gene_median_pre, color = "red",   linetype = "dotted", linewidth = 0.9) +
    labs(
      title    = "Detected genes per cell",
      subtitle = NULL,
      x        = "Detected genes",
      y        = "Number of cells"
    ) +
    presentation_theme(base_size = 12)

  if (!is.null(cell_max_genes)) {
    p1 <- p1 + geom_vline(xintercept = cell_max_genes, color = "black", linetype = "dotted", linewidth = 0.9)
  }

  p2 <- ggplot(metadata, aes(x = total_expr_plot)) +
    geom_histogram(bins = 50, fill = "darkgreen", color = "white") +
    geom_vline(xintercept = min_count,        color = "black", linetype = "dotted", linewidth = 0.9) +
    geom_vline(xintercept = count_median_pre, color = "red",   linetype = "dotted", linewidth = 0.9) +
    scale_x_log10() +
    labs(
      title    = "Total counts per cell",
      subtitle = NULL,
      x        = log10_axis_label("Total counts per cell"),
      y        = "Number of cells"
    ) +
    presentation_theme(base_size = 12)

  if ("mito_pct" %in% names(metadata) && !is.null(max_mito_pct)) {
    mito_median_pre <- stats::median(metadata$mito_pct, na.rm = TRUE)

    p3 <- ggplot(metadata, aes(x = mito_pct)) +
      geom_histogram(bins = 50, fill = "coral", color = "white") +
      geom_vline(xintercept = max_mito_pct,    color = "black", linetype = "dotted", linewidth = 0.9) +
      geom_vline(xintercept = mito_median_pre, color = "red",   linetype = "dotted", linewidth = 0.9) +
      labs(
        title    = "Mitochondrial fraction",
        subtitle = NULL,
        x        = "Mitochondrial reads (%)",
        y        = "Number of cells"
      ) +
      presentation_theme(base_size = 12)

    combined <- (p1 | p2 | p3) +
      plot_annotation(
        title = sample_plot_title(sample_id, "Quality control metrics before filtering"),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
      )
  } else {
    combined <- (p1 | p2) +
      plot_annotation(
        title = sample_plot_title(sample_id, "Quality control metrics before filtering"),
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
      )
  }
  
  save_presentation_plot(
    plot = combined,
    filename = file.path(results_folder, paste0(sample_id, "_qc_pre_filtering.png")),
    width = 18,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  
  cat("✓ Pre-filtering plots saved\n\n")

  # Per-FOV diagnostic heat strip: median genes/cell + cell count per FOV.
  # Quick "is one FOV bad?" view before filtering.
  if ("fov" %in% names(metadata)) {
    fov_summary <- metadata %>%
      dplyr::filter(!is.na(fov)) %>%
      dplyr::group_by(fov) %>%
      dplyr::summarise(
        n_cells       = dplyr::n(),
        median_genes  = stats::median(nr_feats, na.rm = TRUE),
        median_counts = stats::median(total_expr, na.rm = TRUE),
        .groups       = "drop"
      ) %>%
      tidyr::pivot_longer(
        cols      = c("n_cells", "median_genes", "median_counts"),
        names_to  = "metric",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        metric = dplyr::recode(metric,
                               n_cells       = "Cells per FOV",
                               median_genes  = "Median genes per cell",
                               median_counts = "Median counts per cell"),
        # Per-metric z-score so all three rows share a 0-centred fill scale.
        value_z = (value - stats::ave(value, metric, FUN = function(v) mean(v, na.rm = TRUE))) /
                  pmax(stats::ave(value, metric, FUN = function(v) stats::sd(v, na.rm = TRUE)), 1e-6)
      )

    p_fov <- ggplot(fov_summary, aes(x = factor(fov), y = metric, fill = value_z)) +
      geom_tile(colour = "white", linewidth = 0.2) +
      geom_text(aes(label = round(value, 1)), size = 2.6, colour = "grey15") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "Z-score") +
      labs(
        title    = sample_plot_title(sample_id, "QC metrics per FOV"),
        subtitle = NULL,
        x        = "FOV",
        y        = NULL
      ) +
      presentation_theme(base_size = 11, x_angle = 45)

    save_presentation_plot(
      plot     = p_fov,
      filename = file.path(results_folder, paste0(sample_id, "_qc_per_fov_strip.png")),
      width    = max(10, 0.25 * length(unique(fov_summary$fov)) + 4),
      height   = 6,
      dpi      = 300
    )
    cat("  Wrote per-FOV QC strip\n")
  }

  # Spatial QC plots - polygon rendering (cells drawn as real polygons, not
  # points). Each metric keeps its own PNG, and all metrics are additionally
  # combined into one enlarged side-by-side figure. For composite samples,
  # per-sub-biopsy top-bottom stacked variants are written into a subfolder.
  cat("Creating spatial QC plots (polygon rendering)...\n")

  has_mito <- "mito_pct" %in% names(metadata)
  spatial_metrics <- list(
    list(column = "nr_feats",   label = "Detected genes per cell",    slug = "genes_per_cell"),
    list(column = "total_expr", label = "Total counts per cell",      slug = "total_counts")
  )
  if (has_mito) {
    spatial_metrics[[length(spatial_metrics) + 1]] <-
      list(column = "mito_pct", label = "Mitochondrial fraction (%)", slug = "mito_pct")
  }

  # Render one polygon panel per metric (per-sample context: linewidth stays
  # at the existing project default rather than PER_FOV_LINEWIDTH).
  .render_qc_polygon_panel <- function(gobject_local, metric) {
    tryCatch(
      plot_cells_polygon(
        gobject        = gobject_local,
        fill_column    = metric$column,
        fill_as_factor = FALSE,
        context        = "sample",
        polygon_alpha  = 0.9,
        save_plot      = FALSE,
        return_plot    = TRUE,
        show_plot      = FALSE
      ) + labs(title = metric$label),
      error = function(e) {
        cat("  ⚠ ", metric$column, ": ", conditionMessage(e), "\n", sep = "")
        NULL
      }
    )
  }

  # Per-metric individual PNGs + one combined side-by-side PNG for the full sample.
  .save_qc_spatial_group <- function(gobject_local,
                                     out_dir,
                                     file_prefix,
                                     stack_direction = c("side", "stack"),
                                     title_text = NULL) {
    stack_direction <- match.arg(stack_direction)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    panels <- lapply(spatial_metrics, function(m) .render_qc_polygon_panel(gobject_local, m))
    names(panels) <- vapply(spatial_metrics, function(m) m$slug, character(1))
    ok <- !vapply(panels, is.null, logical(1))
    if (!any(ok)) return(invisible(FALSE))

    # Individual files (one per metric). Kept alongside the combined file so
    # both CART composite runs and regular samples have single-metric views.
    for (i in which(ok)) {
      m <- spatial_metrics[[i]]
      save_presentation_plot(
        plot     = panels[[i]],
        filename = file.path(out_dir, paste0(file_prefix, "_spatial_", m$slug, ".png")),
        width    = 14,
        height   = 11,
        dpi      = 600,
        bg       = "white"
      )
    }

    # Combined figure.
    combined <- if (stack_direction == "side") {
      Reduce(`|`, panels[ok])
    } else {
      Reduce(`/`, panels[ok])
    }
    if (!is.null(title_text)) {
      combined <- combined + plot_annotation(
        title = title_text,
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
      )
    }
    combined_width  <- if (stack_direction == "side") 14 * sum(ok) else 16
    combined_height <- if (stack_direction == "side") 11          else 12 * sum(ok)
    combined_name   <- if (stack_direction == "side") {
      paste0(file_prefix, "_spatial_qc_combined.png")
    } else {
      paste0(file_prefix, "_spatial_qc_stacked.png")
    }
    save_presentation_plot(
      plot     = combined,
      filename = file.path(out_dir, combined_name),
      width    = combined_width,
      height   = combined_height,
      dpi      = 600,
      bg       = "white"
    )
    invisible(TRUE)
  }

  tryCatch({
    # Main sample: individual PNGs + side-by-side combined figure.
    .save_qc_spatial_group(
      gobject_local   = gobj,
      out_dir         = results_folder,
      file_prefix     = sample_id,
      stack_direction = "side",
      title_text      = sample_plot_title(sample_id, "Spatial quality control")
    )

    # Composite samples: also emit per-sub-biopsy top-bottom stacked variants.
    sub_rows <- discover_composite_subsamples(sample_row, sample_sheet_path)
    if (!is.null(sub_rows) && nrow(sub_rows) > 0) {
      cat("Composite sample detected - rendering", nrow(sub_rows),
          "per-sub-biopsy stacked variant(s)...\n")
      meta_all <- as.data.frame(pDataDT(gobj))
      if (!"fov" %in% names(meta_all)) {
        cat("  ⚠ No 'fov' column on Giotto object; skipping sub-biopsy split\n")
      } else {
        sub_dir <- file.path(results_folder, "subsamples")
        for (k in seq_len(nrow(sub_rows))) {
          sub_r  <- sub_rows[k, , drop = FALSE]
          sub_id <- as.character(sub_r$sample_id)
          fmin   <- as.integer(sub_r$fov_min)
          fmax   <- as.integer(sub_r$fov_max)
          if (anyNA(c(fmin, fmax))) {
            cat("  ⚠ ", sub_id, ": fov_min/fov_max missing, skipped\n", sep = "")
            next
          }
          cell_ids <- meta_all$cell_ID[
            !is.na(meta_all$fov) & meta_all$fov >= fmin & meta_all$fov <= fmax
          ]
          if (length(cell_ids) == 0) {
            cat("  ⚠ ", sub_id, ": no cells in FOV ", fmin, "-", fmax,
                ", skipped\n", sep = "")
            next
          }
          sub_gobj <- tryCatch(
            subsetGiotto(gobj, cell_ids = cell_ids),
            error = function(e) {
              cat("  ⚠ ", sub_id, ": subsetGiotto failed: ",
                  conditionMessage(e), "\n", sep = "")
              NULL
            }
          )
          if (is.null(sub_gobj)) next
          cat("  → ", sub_id, " (FOV ", fmin, "-", fmax,
              ", ", length(cell_ids), " cells)\n", sep = "")
          .save_qc_spatial_group(
            gobject_local   = sub_gobj,
            out_dir         = sub_dir,
            file_prefix     = sub_id,
            stack_direction = "stack",
            title_text      = sample_plot_title(sub_id, "Spatial quality control")
          )
        }
      }
    }

    cat("✓ Spatial plots saved\n\n")
  }, error = function(e) {
    cat("⚠ Spatial plotting warning:", conditionMessage(e), "\n\n")
  })
  
  # Filter genes
  cat("Filtering genes...\n")
  cat("  Minimum cells per gene:", gene_min_cells, "\n")
  
  gobj <- filterGiotto(
    gobject = gobj,
    feat_type = "rna",
    expression_threshold = 1,
    feat_det_in_min_cells = gene_min_cells,
    min_det_feats_per_cell = 0  # Don't filter cells yet
  )
  
  n_genes_after <- length(gobj@feat_ID$rna)
  cat("✓ Genes after filtering:", n_genes_after, 
      "(removed", n_genes_initial - n_genes_after, ")\n\n")
  
  # Filter cells
  cat("Filtering cells...\n")
  cat("  Minimum genes per cell:", cell_min_genes, "\n")
  if (!is.null(cell_max_genes)) {
    cat("  Maximum genes per cell:", cell_max_genes, "\n")
  }
  cat("  Minimum total counts:", min_count, "\n")
  if (!is.null(max_mito_pct)) {
    cat("  Maximum mitochondrial %:", max_mito_pct, "\n")
  }
  
  # Build cell filter criteria
  cell_metadata <- pDataDT(gobj)
  
  cells_to_keep <- cell_metadata$cell_ID[
    cell_metadata$nr_feats >= cell_min_genes &
      cell_metadata$total_expr >= min_count
  ]
  
  if (!is.null(cell_max_genes)) {
    cells_to_keep <- intersect(
      cells_to_keep,
      cell_metadata$cell_ID[cell_metadata$nr_feats <= cell_max_genes]
    )
  }
  
  if (!is.null(max_mito_pct) && "mito_pct" %in% names(cell_metadata)) {
    cells_to_keep <- intersect(
      cells_to_keep,
      cell_metadata$cell_ID[cell_metadata$mito_pct <= max_mito_pct]
    )
  }
  
  if (length(cells_to_keep) == 0) {
    stop("All cells were removed by QC thresholds. Relax the QC parameters for this sample.")
  }
  
  # Apply cell filter
  gobj <- subsetGiotto(
    gobject = gobj,
    cell_ids = cells_to_keep
  )
  
  n_cells_after <- length(gobj@cell_ID$cell)
  cat("✓ Cells after filtering:", n_cells_after,
      "(removed", n_cells_initial - n_cells_after, ")\n\n")
  
  # Post-filtering plots
  cat("Creating post-filtering QC plots...\n")
  
  metadata_post <- pDataDT(gobj) %>%
    as_tibble() %>%
    mutate(total_expr_plot = pmax(total_expr, 1))
  
  gene_median_post  <- stats::median(metadata_post$nr_feats,        na.rm = TRUE)
  count_median_post <- stats::median(metadata_post$total_expr_plot, na.rm = TRUE)

  # Post-filter distribution plots: single red dotted line at the median. The
  # threshold is no longer drawn because all remaining cells are above it.
  p1_post <- ggplot(metadata_post, aes(x = nr_feats)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white") +
    geom_vline(xintercept = gene_median_post, color = "red", linetype = "dotted", linewidth = 0.9) +
    labs(
      title    = "Detected genes per cell",
      subtitle = NULL,
      x        = "Detected genes",
      y        = "Number of cells"
    ) +
    presentation_theme(base_size = 12)

  p2_post <- ggplot(metadata_post, aes(x = total_expr_plot)) +
    geom_histogram(bins = 50, fill = "darkgreen", color = "white") +
    geom_vline(xintercept = count_median_post, color = "red", linetype = "dotted", linewidth = 0.9) +
    scale_x_log10() +
    labs(
      title    = "Total counts per cell",
      subtitle = NULL,
      x        = log10_axis_label("Total counts per cell"),
      y        = "Number of cells"
    ) +
    presentation_theme(base_size = 12)

  combined_post <- (p1_post | p2_post) +
    plot_annotation(
      title = sample_plot_title(sample_id, "Quality control metrics after filtering"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
    )
  
  save_presentation_plot(
    plot = combined_post,
    filename = file.path(results_folder, paste0(sample_id, "_qc_post_filtering.png")),
    width = 14,
    height = 6,
    dpi = 300,
    bg = "white"
  )

  # Before-vs-after density overlay: shows how filtering reshapes the
  # distribution. Both metrics on log10 scale to pull out the long tail.
  density_df <- dplyr::bind_rows(
    metadata      %>% dplyr::transmute(stage = "Before", nr_feats, total_expr_plot),
    metadata_post %>% dplyr::transmute(stage = "After",  nr_feats, total_expr_plot)
  ) %>%
    dplyr::mutate(stage = factor(stage, levels = c("Before", "After")))

  p_dens_genes <- ggplot(density_df, aes(x = nr_feats, fill = stage, colour = stage)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    scale_fill_manual(values = c(Before = "#777777", After = "#2166AC"),
                      name = "Filter stage") +
    scale_colour_manual(values = c(Before = "#444444", After = "#1A4F8A"),
                        guide = "none") +
    labs(title = "Detected genes per cell", subtitle = NULL,
         x = "Detected genes", y = "Density") +
    presentation_theme(base_size = 11)

  p_dens_counts <- ggplot(density_df, aes(x = total_expr_plot, fill = stage, colour = stage)) +
    geom_density(alpha = 0.35, linewidth = 0.7) +
    scale_x_log10() +
    scale_fill_manual(values = c(Before = "#777777", After = "#2166AC"),
                      name = "Filter stage") +
    scale_colour_manual(values = c(Before = "#444444", After = "#1A4F8A"),
                        guide = "none") +
    labs(title = "Total counts per cell", subtitle = NULL,
         x = log10_axis_label("Total counts per cell"), y = "Density") +
    presentation_theme(base_size = 11)

  density_combined <- (p_dens_genes | p_dens_counts) +
    plot_annotation(
      title = sample_plot_title(sample_id, "QC effect on metric distributions"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
    ) +
    plot_layout(guides = "collect")

  save_presentation_plot(
    plot = density_combined,
    filename = file.path(results_folder, paste0(sample_id, "_qc_before_after_density.png")),
    width = 14,
    height = 5,
    dpi = 300,
    bg = "white"
  )

  cat("✓ Post-filtering plots saved\n\n")
  
  # Save QC summary
  qc_summary <- tibble(
    metric = c("Cells (initial)", "Cells (final)", "Cells (removed)",
               "Genes (initial)", "Genes (final)", "Genes (removed)",
               "Median genes/cell", "Median counts/cell"),
    value = c(n_cells_initial, n_cells_after, n_cells_initial - n_cells_after,
              n_genes_initial, n_genes_after, n_genes_initial - n_genes_after,
              median(metadata_post$nr_feats, na.rm = TRUE),
              median(metadata_post$total_expr, na.rm = TRUE))
  )
  
  write_csv(qc_summary, file.path(results_folder, paste0(sample_id, "_qc_summary.csv")))

  # Audit trail: record the exact thresholds used for this run so downstream
  # analysts can reproduce filtering decisions without re-running the pipeline.
  qc_parameters <- tibble(
    parameter = c("gene_min_cells", "cell_min_genes", "cell_max_genes",
                  "min_count", "max_mito_pct",
                  "n_mt_genes_detected",
                  "n_control_probes_removed",
                  "n_cells_initial", "n_cells_final",
                  "n_genes_initial", "n_genes_final",
                  "timestamp"),
    value = c(as.character(gene_min_cells),
              as.character(cell_min_genes),
              ifelse(is.null(cell_max_genes), "NULL", as.character(cell_max_genes)),
              as.character(min_count),
              ifelse(is.null(max_mito_pct), "NULL", as.character(max_mito_pct)),
              as.character(length(mt_genes)),
              as.character(n_controls_removed),
              as.character(n_cells_initial),
              as.character(n_cells_after),
              as.character(n_genes_initial),
              as.character(n_genes_after),
              format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  )
  write_csv(
    qc_parameters,
    file.path(results_folder, paste0(sample_id, "_qc_parameters.csv"))
  )
  
  cat("=== Final Statistics ===\n")
  print(qc_summary)
  cat("\n")
  
  cat("✓ Quality control complete for", sample_id, "\n\n")
  
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
    
    sample_id <- args[1]
    input_path <- args[2]
    output_dir <- args[3]

    # Resolve optional composite context from an environment variable so the
    # CLI entrypoint can emit per-sub-biopsy stacked spatial plots without
    # changing its argument signature.
    cli_sheet_path <- Sys.getenv("COSMX_SAMPLE_SHEET", unset = "")
    cli_sample_row <- NULL
    if (nzchar(cli_sheet_path) && file.exists(cli_sheet_path) &&
        exists("safe_read_sheet", mode = "function")) {
      sheet_df <- tryCatch(
        as.data.frame(safe_read_sheet(cli_sheet_path), stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (!is.null(sheet_df) && "sample_id" %in% names(sheet_df)) {
        hit <- sheet_df[as.character(sheet_df$sample_id) == sample_id, , drop = FALSE]
        if (nrow(hit) >= 1) cli_sample_row <- hit[1, , drop = FALSE]
      }
    }

    gobj <- quality_control(
      gobj              = input_path,
      sample_id         = sample_id,
      output_dir        = output_dir,
      sample_row        = cli_sample_row,
      sample_sheet_path = if (nzchar(cli_sheet_path)) cli_sheet_path else NULL
    )
    
    # Save filtered object
    saveGiotto(
      gobject = gobj,
      dir = output_dir,
      foldername = "Giotto_Object_Filtered",
      overwrite = TRUE
    )
  } else {
    stop("Usage: Rscript 02_Quality_Control.R <sample_id> <input_path> <output_dir>")
  }
}

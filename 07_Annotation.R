# Cell type annotation using InSituType with multiple reference profiles -------

#!/usr/bin/env Rscript
# 07_Annotation.R

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
if ((!exists("presentation_theme") || !exists("sample_plot_title") ||
     !exists("pretty_plot_label") || !exists("save_presentation_plot")) &&
    file.exists(pipeline_utils)) {
  source(pipeline_utils)
}


# Colour palette ----------------------------------------------------------

.annotation_palette <- function(n) {
  base_cols <- c(
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(12, "Paired"),
    RColorBrewer::brewer.pal(8,  "Dark2"),
    RColorBrewer::brewer.pal(8,  "Set1")
  )
  base_cols <- unique(base_cols)
  if (n <= length(base_cols)) return(base_cols[seq_len(n)])
  grDevices::colorRampPalette(base_cols)(n)
}

#' Build a shared named colour map for a set of cell type labels.
#' The ordering follows n_cells descending (same as the flightpath legend),
#' so both the flightpath and the UMAP receive the exact same colour per type.
#'
#' @param clust_vec  Named or unnamed character vector of per-cell assignments
#' @return Named character vector: names = cell type, values = hex colour
.build_colour_map <- function(clust_vec) {
  ct_counts  <- sort(table(clust_vec), decreasing = TRUE)
  ct_ordered <- names(ct_counts)
  pal        <- .annotation_palette(length(ct_ordered))
  stats::setNames(pal, ct_ordered)
}


# Helper: extract UMAP matrix from Giotto robustly ------------------------

.get_umap_df <- function(gobj) {
  
  umap_mat <- tryCatch(
    getDimReduction(
      gobject          = gobj,
      spat_unit        = "cell",
      feat_type        = "rna",
      reduction        = "cells",
      reduction_method = "umap",
      name             = "umap",
      output           = "matrix"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(umap_mat) && is.matrix(umap_mat) && ncol(umap_mat) >= 2) {
    return(data.frame(
      cell_ID = rownames(umap_mat),
      UMAP_1  = umap_mat[, 1],
      UMAP_2  = umap_mat[, 2],
      stringsAsFactors = FALSE
    ))
  }
  
  cat("    \u26A0 matrix output failed \u2014 trying data.table fallback\n")
  
  umap_dt <- tryCatch(
    getDimReduction(
      gobject          = gobj,
      spat_unit        = "cell",
      feat_type        = "rna",
      reduction        = "cells",
      reduction_method = "umap",
      name             = "umap",
      output           = "data.table"
    ),
    error = function(e) NULL
  )
  
  if (is.null(umap_dt)) {
    cat("    \u26A0 data.table output also failed\n")
    return(NULL)
  }
  
  umap_df    <- as.data.frame(umap_dt)
  non_id     <- colnames(umap_df)[!colnames(umap_df) %in%
                                    c("cell_ID", "sdimx", "sdimy")]
  coord_cols <- non_id[seq_len(min(2, length(non_id)))]
  
  if (length(coord_cols) < 2) {
    cat("    \u26A0 Could not identify 2 coordinate columns\n")
    return(NULL)
  }
  
  id_col <- if ("cell_ID" %in% colnames(umap_df)) {
    umap_df$cell_ID
  } else {
    rownames(umap_df)
  }
  
  data.frame(
    cell_ID = id_col,
    UMAP_1  = umap_df[[coord_cols[1]]],
    UMAP_2  = umap_df[[coord_cols[2]]],
    stringsAsFactors = FALSE
  )
}


# Helper: per-cluster mean confidence from insitutype logliks -------------

# insitutypeML() returns $prob as a plain numeric vector (per-cell max
# posterior for the assigned type only). insitutype() returns a full
# cells x types matrix.
#
# To handle BOTH cases robustly we always derive confidence from $logliks,
# which is a cells x types matrix present in all insitutype result objects.
# Per-cell softmax normalisation converts logliks to posteriors, then we
# take the column mean for each cluster type.
#
# @param insitu_result  Output of insitutypeML() or insitutype()
# @return Named numeric vector: names = cluster type, values = mean posterior

.cluster_mean_confidence <- function(insitu_result) {
  
  logliks <- insitu_result$logliks
  clust   <- insitu_result$clust
  
  if (is.null(logliks) || !is.matrix(logliks))
    stop("insitu_result$logliks is NULL or not a matrix")
  
  # Row-wise softmax (subtract row max for numerical stability)
  log_max  <- apply(logliks, 1, max, na.rm = TRUE)
  exp_liks <- exp(logliks - log_max)
  post_mat <- exp_liks / rowSums(exp_liks, na.rm = TRUE)
  
  ct_types  <- colnames(post_mat)
  mean_conf <- vapply(ct_types, function(ct) {
    cell_idx <- which(clust == ct)
    if (length(cell_idx) == 0) return(NA_real_)
    mean(post_mat[cell_idx, ct], na.rm = TRUE)
  }, numeric(1))
  
  names(mean_conf) <- ct_types
  mean_conf
}


# ── Annotation quality scoring & auto-selection ───────────────────────────
#
# Expected cell-type composition for lupus nephritis / autoimmune kidney
# disease.  Each category carries grepl patterns (case-insensitive) that
# match the cell-type names produced by any of the supported reference
# profiles, plus a biological target proportion and a plausible [min, max]
# range derived from published LN single-cell and spatial studies (e.g.,
# Arazi 2019 Nat Immunol, Monteagudo 2021 Nat Commun, Hinze 2022 ARD).

.expected_composition_ln <- function() {
  list(
    proximal_tubule = list(
      patterns = c("proximal.tubule", "proximal_tubule", "\\bPT\\b",
                   "distinct.proximal", "\\bS1\\b", "\\bS2\\b", "\\bS3\\b"),
      target = 0.35, min = 0.10, max = 0.55
    ),
    loop_henle_TAL = list(
      patterns = c("thick.ascending", "loop.of.henle", "\\bTAL\\b",
                   "henle"),
      target = 0.08, min = 0.02, max = 0.20
    ),
    collecting_duct = list(
      patterns = c("collecting.duct", "principal.cell", "intercalated",
                   "connecting.tubule", "\\bDCT\\b", "\\bCNT\\b"),
      target = 0.07, min = 0.01, max = 0.18
    ),
    endothelial = list(
      patterns = c("endothel", "glomerular.endothel", "peritubular",
                   "vasa.recta", "ascending.vasa", "descending.vasa"),
      target = 0.08, min = 0.02, max = 0.18
    ),
    podocyte = list(
      patterns = c("podocyte", "glomerular.visceral", "\\bPOD\\b"),
      target = 0.03, min = 0.003, max = 0.10
    ),
    fibroblast_stromal = list(
      patterns = c("fibroblast", "myofibroblast", "pericyte",
                   "mesangial", "stromal"),
      target = 0.05, min = 0.005, max = 0.15
    ),
    T_cell = list(
      # elevated in LN
      patterns = c("CD4.T", "CD8.T", "\\bT.cell\\b", "\\bNKT\\b",
                   "T_cell", "regulatory.T"),
      target = 0.09, min = 0.02, max = 0.22
    ),
    B_cell_plasma = list(
      # especially elevated in LN (auto-antibody producers)
      patterns = c("\\bB.cell\\b", "\\bB_cell\\b", "plasma.cell",
                   "plasmablast", "plasmacyt"),
      target = 0.06, min = 0.005, max = 0.18
    ),
    macrophage_MNP = list(
      patterns = c("macrophage", "\\bMNP\\b", "monocyte", "dendritic",
                   "myeloid"),
      target = 0.08, min = 0.02, max = 0.20
    ),
    NK_other_immune = list(
      patterns = c("\\bNK.cell\\b", "natural.killer", "mast.cell",
                   "neutrophil", "plasmacytoid.DC",
                   "plasmacytoid.dendritic"),
      target = 0.02, min = 0.0, max = 0.08
    ),
    pelvic_urothelial = list(
      patterns = c("pelvic", "urothelial", "transitional", "\\bUPK\\b"),
      target = 0.01, min = 0.0, max = 0.06
    )
  )
}


#' Score a single annotation result on confidence and biological composition
#'
#' @param celltype_col   Name of the cell-type metadata column
#' @param score_col      Name of the per-cell confidence score column
#' @param metadata       data.frame / data.table of Giotto cell metadata
#' @param expected_comp  Output of \code{.expected_composition_ln()}
#' @param conf_threshold Threshold for "high-confidence" cells (default 0.8)
#' @return Named list of per-metric and composite scores

score_annotation_quality <- function(celltype_col,
                                     score_col,
                                     metadata,
                                     expected_comp,
                                     conf_threshold = 0.8) {

  scores <- as.numeric(metadata[[score_col]])
  scores <- scores[!is.na(scores)]

  mean_conf      <- mean(scores)
  median_conf    <- median(scores)
  high_conf_frac <- mean(scores >= conf_threshold)

  ct_counts <- table(as.character(metadata[[celltype_col]]))
  ct_props  <- ct_counts / sum(ct_counts)
  ct_names  <- names(ct_props)

  # Map each cell type to expected LN categories via regex
  category_proportions <- vapply(expected_comp, function(cat) {
    matched <- vapply(cat$patterns, function(pat) {
      hits <- grepl(pat, ct_names, ignore.case = TRUE, perl = TRUE)
      sum(ct_props[hits])
    }, numeric(1))
    min(sum(matched), 1.0)
  }, numeric(1))

  # Score each category: full credit inside expected range, linear penalty outside
  comp_scores <- vapply(names(expected_comp), function(cat_name) {
    cat <- expected_comp[[cat_name]]
    obs <- category_proportions[[cat_name]]
    if (obs >= cat$min && obs <= cat$max) {
      1.0 - abs(obs - cat$target) / max(cat$target, 0.01)
    } else {
      gap <- max(cat$min - obs, obs - cat$max, 0)
      max(0.0, 1.0 - gap / max(cat$target, 0.01))
    }
  }, numeric(1))
  composition_score <- mean(comp_scores)

  # Fraction of expected categories with >= 0.5% cells present
  coverage_score <- mean(category_proportions >= 0.005)

  composite_score <- (
    0.35 * mean_conf       +
    0.25 * high_conf_frac  +
    0.30 * composition_score +
    0.10 * coverage_score
  )

  list(
    annotation_column    = celltype_col,
    score_column         = score_col,
    mean_confidence      = round(mean_conf,         4),
    median_confidence    = round(median_conf,       4),
    high_conf_fraction   = round(high_conf_frac,    4),
    composition_score    = round(composition_score, 4),
    coverage_score       = round(coverage_score,    4),
    composite_score      = round(composite_score,   4),
    n_cell_types         = length(unique(as.character(metadata[[celltype_col]]))),
    category_proportions = round(category_proportions, 4)
  )
}


#' Score all supervised annotation columns and let the user pick one
#'
#' When a TTY is attached the user is shown a numbered table and prompted to
#' choose.  In non-interactive / batch mode the column with the highest
#' composite score is selected automatically.
#'
#' Writes \code{annotation_selection.json} and
#' \code{{sample_id}_annotation_scores.csv} to \code{annotation_dir}.
#'
#' @param gobj           Giotto object with annotation metadata already added
#' @param annotation_dir Path to the \code{07_Annotation/} output folder
#' @param sample_id      Sample identifier string
#' @param conf_threshold Threshold passed to \code{score_annotation_quality()}
#' @return List with \code{selected_col} and \code{scores_df}, or NULL on failure

select_best_annotation <- function(gobj,
                                   annotation_dir,
                                   sample_id,
                                   conf_threshold = 0.8) {

  dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)

  metadata      <- as.data.frame(pDataDT(gobj))
  expected_comp <- .expected_composition_ln()
  all_cols      <- names(metadata)

  # Candidate columns: supervised or supervised_refined only (not semi)
  ct_cols <- grep("^celltype_.*_(supervised|supervised_refined)$",
                  all_cols, value = TRUE)

  if (length(ct_cols) == 0) {
    cat("  \u26A0 No supervised annotation columns found; skipping selection.\n")
    return(NULL)
  }

  score_cols <- sub("^celltype_", "score_", ct_cols)
  valid      <- score_cols %in% all_cols
  ct_cols    <- ct_cols[valid]
  score_cols <- score_cols[valid]

  if (length(ct_cols) == 0) {
    cat("  \u26A0 No matching confidence-score columns found; skipping selection.\n")
    return(NULL)
  }

  scores_list <- mapply(function(ct_col, sc_col) {
    tryCatch(
      score_annotation_quality(ct_col, sc_col, metadata,
                               expected_comp, conf_threshold),
      error = function(e) {
        cat("    \u26A0 Could not score '", ct_col, "': ",
            conditionMessage(e), "\n", sep = "")
        NULL
      }
    )
  }, ct_cols, score_cols, SIMPLIFY = FALSE)

  scores_list <- Filter(Negate(is.null), scores_list)

  if (length(scores_list) == 0) {
    cat("  \u26A0 All scoring attempts failed; skipping selection.\n")
    return(NULL)
  }

  scores_df <- do.call(rbind, lapply(scores_list, function(s) {
    data.frame(
      annotation_column  = s$annotation_column,
      score_column       = s$score_column,
      mean_confidence    = s$mean_confidence,
      median_confidence  = s$median_confidence,
      high_conf_fraction = s$high_conf_fraction,
      composition_score  = s$composition_score,
      coverage_score     = s$coverage_score,
      composite_score    = s$composite_score,
      n_cell_types       = s$n_cell_types,
      stringsAsFactors   = FALSE
    )
  }))
  scores_df <- scores_df[order(-scores_df$composite_score), ]

  # ── Display scored table ────────────────────────────────────────────────
  cat(sprintf(
    "\n=== Annotation Selection: %s ===\n", sample_id
  ))
  cat(sprintf(
    "%-3s  %-45s  %9s  %6s  %9s  %11s  %7s\n",
    "#", "Column", "composite", "conf", "high_conf", "composition", "n_types"
  ))
  cat(strrep("-", 95), "\n")
  for (i in seq_len(nrow(scores_df))) {
    r <- scores_df[i, ]
    cat(sprintf(
      "%-3d  %-45s  %9.3f  %6.3f  %9.3f  %11.3f  %7d\n",
      i,
      substr(r$annotation_column, 1, 45),
      r$composite_score,
      r$mean_confidence,
      r$high_conf_fraction,
      r$composition_score,
      r$n_cell_types
    ))
  }
  cat(strrep("-", 95), "\n\n")

  # ── Interactive or auto selection ───────────────────────────────────────
  use_interactive <- interactive() || isatty(stdin())

  if (use_interactive && nrow(scores_df) > 1) {
    choices <- sprintf(
      "%s  [composite=%.3f, conf=%.3f, n_types=%d]",
      scores_df$annotation_column,
      scores_df$composite_score,
      scores_df$mean_confidence,
      scores_df$n_cell_types
    )
    chosen_idx <- menu(choices,
                       title = paste0("Select annotation for ", sample_id,
                                      " (0 = auto-select best):"))
    if (chosen_idx == 0L) {
      chosen_idx <- 1L
      cat("  Auto-selecting row 1 (highest composite score).\n")
    }
  } else {
    chosen_idx <- 1L
    if (nrow(scores_df) > 1) {
      cat("  [Non-interactive] Auto-selecting annotation with highest composite score.\n")
    }
  }

  best     <- scores_df[chosen_idx, ]
  best_col <- best$annotation_column

  cat(sprintf(
    "  \u2713 Selected: %s\n    composite=%.3f  conf=%.3f  high_conf_frac=%.3f  composition=%.3f\n",
    best_col, best$composite_score, best$mean_confidence,
    best$high_conf_fraction, best$composition_score
  ))

  # Write summary CSV
  write.csv(
    scores_df,
    file.path(annotation_dir,
              paste0(sample_id, "_annotation_scores.csv")),
    row.names = FALSE
  )

  # Write machine-readable JSON for the pipeline orchestrator
  selection_json <- list(
    sample_id                  = sample_id,
    selected_annotation_column = best_col,
    selected_score_column      = as.character(best$score_column),
    composite_score            = best$composite_score,
    mean_confidence            = best$mean_confidence,
    high_conf_fraction         = best$high_conf_fraction,
    composition_score          = best$composition_score,
    coverage_score             = best$coverage_score,
    n_cell_types               = best$n_cell_types,
    selection_reason           = paste0(
      "Highest composite score (confidence x0.60 + LN-composition x0.30 ",
      "+ coverage x0.10) across all supervised annotation results."
    ),
    expected_composition_ref   = paste0(
      "Lupus nephritis / autoimmune kidney disease ",
      "(Arazi 2019, Monteagudo 2021, Hinze 2022)"
    ),
    all_scores = lapply(scores_list, function(s) {
      list(
        annotation_column    = s$annotation_column,
        composite_score      = s$composite_score,
        mean_confidence      = s$mean_confidence,
        high_conf_fraction   = s$high_conf_fraction,
        composition_score    = s$composition_score,
        coverage_score       = s$coverage_score,
        n_cell_types         = s$n_cell_types,
        category_proportions = as.list(s$category_proportions)
      )
    }),
    selection_timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )

  jsonlite::write_json(
    selection_json,
    path       = file.path(annotation_dir, "annotation_selection.json"),
    pretty     = TRUE,
    auto_unbox = TRUE
  )
  cat("  \u2713 Selection written to annotation_selection.json\n")

  list(selected_col = best_col, scores_df = scores_df)
}


# Helper: custom UMAP plot ------------------------------------------------

#' @param gobj          Giotto object (must have a umap reduction)
#' @param clust_vec     Named character vector: names = cell_IDs, values = cell types
#' @param colour_map    Named character vector from .build_colour_map()
#' @param profile_name  String used in plot title and file name
#' @param sample_id     String used in plot title and file name
#' @param out_dir       Directory to save the PNG
#' @param ann_type      "supervised", "semi", or "supervised_refined"

plot_giotto_umap <- function(gobj,
                             clust_vec,
                             colour_map,
                             profile_name,
                             sample_id,
                             out_dir,
                             ann_type    = "supervised",
                             point_size  = 0.3,
                             point_alpha = 0.6,
                             label_size  = 3,
                             label_alpha = 0.5,
                             width  = 20,
                             height = 10,
                             dpi    = 300) {
  
  umap_df <- .get_umap_df(gobj)
  
  if (is.null(umap_df)) {
    cat("    \u26A0 Could not extract UMAP coordinates\n")
    return(invisible(NULL))
  }
  
  cat("    UMAP data:", nrow(umap_df), "cells, columns:",
      paste(colnames(umap_df), collapse = ", "), "\n")
  
  ann_df <- tibble::tibble(
    cell_ID  = names(clust_vec),
    CellType = factor(as.character(clust_vec), levels = names(colour_map))
  )
  
  plot_df <- dplyr::inner_join(umap_df, ann_df, by = "cell_ID")
  cat("    After join:", nrow(plot_df), "cells matched\n")
  
  if (nrow(plot_df) == 0) {
    cat("    \u26A0 No overlapping cell IDs between UMAP and annotation\n")
    return(invisible(NULL))
  }
  
  wrapped_levels <- pretty_plot_label(names(colour_map), width = 24)
  wrapped_colour_map <- stats::setNames(unname(colour_map), wrapped_levels)
  plot_df <- plot_df %>%
    dplyr::mutate(
      CellType_label = factor(
        pretty_plot_label(as.character(CellType), width = 24),
        levels = wrapped_levels
      )
    )
  
  centroids <- plot_df %>%
    dplyr::filter(!is.na(CellType_label)) %>%
    dplyr::group_by(CellType_label) %>%
    dplyr::summarise(
      UMAP_1 = stats::median(UMAP_1, na.rm = TRUE),
      UMAP_2 = stats::median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )
  
  title_txt <- sample_plot_title(
    sample_id,
    paste("UMAP Projection -", profile_name, ann_type, "Annotation")
  )
  
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = UMAP_1, y = UMAP_2, colour = CellType_label)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_colour_manual(values = wrapped_colour_map, drop = FALSE) +
    ggrepel::geom_label_repel(
      data         = centroids,
      ggplot2::aes(x = UMAP_1, y = UMAP_2, label = CellType_label),
      colour       = "black",
      fill         = ggplot2::alpha("white", label_alpha),
      size         = label_size,
      label.size   = 0.4,
      label.r      = grid::unit(0.15, "lines"),
      segment.size = 0.2,
      fontface     = "bold",
      show.legend  = FALSE,
      max.overlaps = Inf,
      inherit.aes  = FALSE
    ) +
    ggplot2::labs(
      title = title_txt,
      subtitle = "Cells are colored by annotation, with median cluster labels shown for readability.",
      x = embedding_axis_label("UMAP", 1),
      y = embedding_axis_label("UMAP", 2),
      colour = "Cell Type") +
    ggplot2::guides(colour = ggplot2::guide_legend(
      ncol = 1, override.aes = list(size = 3, alpha = 1))) +
    presentation_theme(base_size = 12, legend_position = "right") +
    ggplot2::theme(
      legend.text = ggplot2::element_text(lineheight = 0.9),
      legend.key.height = grid::unit(0.42, "cm")
    )
  
  fname <- paste0(sample_id, "_umap_", profile_name, "_", ann_type, "_custom.png")
  save_presentation_plot(
    plot = p,
    filename = file.path(out_dir, fname),
    width = width,
    height = height,
    dpi = dpi
  )
  cat("  \u2713 Custom UMAP saved:", fname, "\n")
  invisible(p)
}


# Helper: custom flightpath plot ------------------------------------------

#' @param insitu_result  Output of insitutypeML() / insitutype() / refineClusters()
#' @param colour_map     Named character vector from .build_colour_map()

plot_custom_flightpath <- function(insitu_result,
                                   colour_map,
                                   profile_name,
                                   sample_id,
                                   out_dir,
                                   ann_type    = "supervised",
                                   point_size  = 0.2,
                                   point_alpha = 0.8,
                                   label_size  = 3,
                                   plot_seed   = 1511,
                                   width  = 20,
                                   height = 10,
                                   dpi    = 600) {
  
  if (is.null(insitu_result$logliks) ||
      is.null(insitu_result$profiles) ||
      is.null(insitu_result$clust)) {
    cat("    \u26A0 insitu_result missing logliks / profiles / clust\n")
    return(invisible(NULL))
  }
  
  flight <- tryCatch(
    InSituType::flightpath_layout(
      logliks  = insitu_result$logliks,
      profiles = insitu_result$profiles
    ),
    error = function(e) {
      cat("    \u26A0 flightpath_layout failed:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(flight)) return(invisible(NULL))
  
  cl_levels <- rownames(as.data.frame(flight$clustpos))
  if (is.null(cl_levels) || !length(cl_levels))
    cl_levels <- sort(unique(as.character(insitu_result$clust)))
  
  mean_conf     <- flight$meanconfidence
  mean_conf_map <- if (!is.null(names(mean_conf)) && all(nzchar(names(mean_conf)))) {
    tibble::tibble(Cluster = names(mean_conf), mean_conf = as.numeric(mean_conf))
  } else {
    tibble::tibble(
      Cluster   = cl_levels,
      mean_conf = as.numeric(mean_conf)[
        seq_len(min(length(mean_conf), length(cl_levels)))])
  }
  
  cl_tab <- sort(table(insitu_result$clust), decreasing = TRUE)
  
  Cluster_stats <- tibble::tibble(
    Cluster = factor(names(cl_tab), levels = names(colour_map)),
    n_cells = as.integer(cl_tab)
  ) %>%
    dplyr::left_join(
      mean_conf_map %>%
        dplyr::mutate(Cluster = factor(Cluster, levels = names(colour_map))),
      by = "Cluster"
    ) %>%
    dplyr::mutate(
      Cluster_lab = sprintf(
        "%s\n(n=%s, conf=%.2f)",
        as.character(Cluster),
        scales::comma(n_cells),
        dplyr::if_else(is.na(mean_conf), 0, mean_conf)
      )
    ) %>%
    dplyr::mutate(
      Cluster_lab_pretty = pretty_plot_label(Cluster_lab, width = 24)
    ) %>%
    dplyr::arrange(dplyr::desc(n_cells))
  
  lab_colour_map <- stats::setNames(
    colour_map[as.character(Cluster_stats$Cluster)],
    Cluster_stats$Cluster_lab_pretty
  )
  
  cells_df <- as.data.frame(flight$cellpos) %>%
    dplyr::mutate(
      Cluster = factor(insitu_result$clust, levels = names(colour_map))
    ) %>%
    dplyr::left_join(Cluster_stats, by = "Cluster") %>%
    dplyr::mutate(
      Cluster_lab_pretty = factor(
        pretty_plot_label(Cluster_lab, width = 24),
        levels = Cluster_stats$Cluster_lab_pretty
      )
    )
  
  clustpos_df <- as.data.frame(flight$clustpos) %>%
    tibble::rownames_to_column("Cluster") %>%
    dplyr::mutate(Cluster = factor(Cluster, levels = names(colour_map))) %>%
    dplyr::left_join(Cluster_stats, by = "Cluster") %>%
    dplyr::mutate(Cluster_pretty = pretty_plot_label(Cluster, width = 18))
  
  title_txt <- sample_plot_title(
    sample_id,
    paste("InSituType Flightpath -", profile_name, ann_type)
  )
  set.seed(plot_seed)
  
  p <- ggplot2::ggplot(
    cells_df,
    ggplot2::aes(x = x, y = y, colour = Cluster_lab_pretty)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_colour_manual(values = lab_colour_map, drop = FALSE) +
    ggrepel::geom_label_repel(
      data         = clustpos_df,
      ggplot2::aes(x = x, y = y, label = Cluster_pretty),
      colour = "black", fill = "white", label.size = 0.2,
      fontface = "bold", label.r = grid::unit(0.1, "lines"),
      size = label_size, show.legend = FALSE,
      max.overlaps = Inf, inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title = title_txt,
      subtitle = "Cells are positioned in the InSituType flightpath layout and colored by inferred annotation.",
      colour = "Cluster\n(n cells, confidence)"
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      override.aes = list(size = 4), ncol = 1)) +
    presentation_theme(base_size = 12, legend_position = "right") +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.key.height = grid::unit(0.42, "cm")
    )
  
  fname <- paste0(sample_id, "_flightpath_", profile_name, "_", ann_type, "_custom.png")
  save_presentation_plot(
    plot = p,
    filename = file.path(out_dir, fname),
    width = width,
    height = height,
    dpi = dpi
  )
  cat("  \u2713 Custom flightpath saved:", fname, "\n")
  invisible(p)
}


# Helper: proportions bar chart with percentage labels --------------------

# Shared helper used by supervised, semi-supervised and refined blocks so
# the label/expand logic only lives in one place.
#
# @param summary_df   data.frame with columns: <celltype_col>, n_cells
# @param celltype_col Column name string for the cell type variable
# @param colour_map   Named colour vector
# @param title        Plot title string
# @param out_path     Full file path for the saved PNG

.plot_proportions <- function(summary_df,
                              celltype_col,
                              colour_map,
                              title,
                              out_path,
                              dpi = 300) {
  
  p <- ggplot2::ggplot(
    summary_df,
    ggplot2::aes(
      x    = reorder(!!dplyr::sym(celltype_col), n_cells),
      y    = n_cells,
      fill = !!dplyr::sym(celltype_col))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = colour_map) +
    ggplot2::coord_flip() +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(round(n_cells / sum(n_cells) * 100, 1), "%")
      ),
      hjust    = -0.15,
      size     = 3,
      colour   = "grey20",
      fontface = "bold"
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      title = title,
      x = "Cell Type",
      y = "Number of cells"
    ) +
    presentation_theme(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  
  save_presentation_plot(
    plot = p,
    filename = out_path,
    width = 16,
    height = max(6, nrow(summary_df) * 0.4),
    dpi = dpi
  )
  invisible(p)
}


# Helper: marker gene heatmap ---------------------------------------------

# The heatmap shows z-score-scaled mean log-normalised expression per
# cell type x gene, making it easy to confirm whether each cluster is
# driven by the expected markers (Bruker/NanoString recommended QC).
#
# @param counts_mat      cells x genes matrix (raw counts)
# @param clust_vec       Named character vector: names = cell_IDs,
#                        values = cell type assignments
# @param colour_map      Named colour vector from .build_colour_map()
# @param profile_name    String used in file name / title
# @param sample_id       String used in file name / title
# @param out_dir         Directory to save the PNG
# @param ann_type        "supervised", "semi", or "supervised_refined"
# @param top_n_markers   Top specifically-expressed genes per cluster
# @param min_cells       Minimum cells per cluster to include

.plot_annotation_heatmap <- function(counts_mat,
                                     clust_vec,
                                     colour_map,
                                     profile_name,
                                     sample_id,
                                     out_dir,
                                     ann_type      = "supervised",
                                     top_n_markers = 5,
                                     min_cells     = 5,
                                     width         = 20,
                                     height        = 14,
                                     dpi           = 300) {
  
  shared <- intersect(rownames(counts_mat), names(clust_vec))
  if (length(shared) == 0) {
    cat("    \u26A0 Heatmap: no shared cell IDs between counts_mat and clust_vec\n")
    return(invisible(NULL))
  }
  
  mat   <- counts_mat[shared, , drop = FALSE]
  clust <- clust_vec[shared]
  
  ct_keep <- names(which(table(clust) >= min_cells))
  keep    <- clust %in% ct_keep
  mat     <- mat[keep, , drop = FALSE]
  clust   <- clust[keep]
  
  if (length(unique(clust)) < 2) {
    cat("    \u26A0 Heatmap: fewer than 2 clusters after filtering\n")
    return(invisible(NULL))
  }
  
  # Log-normalise (library-size normalise to median then log1p)
  lib_size <- rowSums(mat)
  med_lib  <- stats::median(lib_size[lib_size > 0])
  norm_mat <- log1p(mat / pmax(lib_size, 1) * med_lib)
  
  # Mean expression per cluster (matrix: genes x cell-types)
  ct_ordered <- names(colour_map)[names(colour_map) %in% unique(clust)]
  mean_expr  <- do.call(cbind, lapply(ct_ordered, function(ct) {
    idx <- which(clust == ct)
    colMeans(norm_mat[idx, , drop = FALSE], na.rm = TRUE)
  }))
  colnames(mean_expr) <- ct_ordered
  
  # Select top_n_markers per cluster by specificity score
  global_mean <- rowMeans(mean_expr, na.rm = TRUE)
  sel_genes   <- unique(unlist(lapply(ct_ordered, function(ct) {
    spec_score <- mean_expr[, ct] - global_mean
    head(names(sort(spec_score, decreasing = TRUE)), top_n_markers)
  })))
  
  if (length(sel_genes) < 2) {
    cat("    \u26A0 Heatmap: not enough marker genes selected\n")
    return(invisible(NULL))
  }
  
  plot_mat <- mean_expr[sel_genes, , drop = FALSE]
  
  # Row-wise z-score, clipped to [-3, 3]
  row_sd  <- pmax(apply(plot_mat, 1, stats::sd, na.rm = TRUE), 1e-6)
  z_mat   <- (plot_mat - rowMeans(plot_mat, na.rm = TRUE)) / row_sd
  z_mat   <- pmin(pmax(z_mat, -3), 3)
  
  ann_col     <- data.frame(CellType = ct_ordered, row.names = ct_ordered)
  ann_colours <- list(CellType = colour_map[ct_ordered])
  
  title_txt <- paste0(sample_id, " \u2014 ", profile_name,
                      " ", ann_type, "\nMarker gene mean expression (z-score)")
  fname <- paste0(sample_id, "_heatmap_", profile_name, "_", ann_type, ".png")
  
  grDevices::png(file.path(out_dir, fname),
                 width = width, height = height, units = "in", res = dpi)
  pheatmap::pheatmap(
    z_mat,
    color             = grDevices::colorRampPalette(
      c("#2166AC", "white", "#B2182B"))(100),
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    annotation_col    = ann_col,
    annotation_colors = ann_colours,
    fontsize_row      = max(5, min(9,  180 / length(sel_genes))),
    fontsize_col      = max(6, min(10, 180 / length(ct_ordered))),
    angle_col         = 45,
    main              = title_txt,
    border_color      = NA,
    legend            = TRUE
  )
  grDevices::dev.off()
  
  cat("  \u2713 Heatmap saved:", fname, "\n")
  invisible(NULL)
}


# Helper: marker gene dot plot --------------------------------------------

# Dot plot combines mean expression (colour) with the fraction of cells
# expressing the gene (dot size), matching the AtoMx SIP recommended
# output format for annotation QC.
#
# @param counts_mat      cells x genes matrix (raw counts)
# @param clust_vec       Named character vector: names = cell_IDs,
#                        values = cell type assignments
# @param colour_map      Named colour vector from .build_colour_map()
# @param top_n_markers   Top specifically-expressed genes per cluster
# @param min_cells       Minimum cells per cluster to include

.plot_annotation_dotplot <- function(counts_mat,
                                     clust_vec,
                                     colour_map,
                                     profile_name,
                                     sample_id,
                                     out_dir,
                                     ann_type      = "supervised",
                                     top_n_markers = 3,
                                     min_cells     = 5,
                                     width         = 20,
                                     height        = 12,
                                     dpi           = 300) {
  
  shared <- intersect(rownames(counts_mat), names(clust_vec))
  if (length(shared) == 0) {
    cat("    \u26A0 Dot plot: no shared cell IDs\n")
    return(invisible(NULL))
  }
  
  mat   <- counts_mat[shared, , drop = FALSE]
  clust <- clust_vec[shared]
  
  ct_keep <- names(which(table(clust) >= min_cells))
  keep    <- clust %in% ct_keep
  mat     <- mat[keep, , drop = FALSE]
  clust   <- clust[keep]
  
  if (length(unique(clust)) < 2) {
    cat("    \u26A0 Dot plot: fewer than 2 clusters after filtering\n")
    return(invisible(NULL))
  }
  
  # Log-normalise
  lib_size <- rowSums(mat)
  med_lib  <- stats::median(lib_size[lib_size > 0])
  norm_mat <- log1p(mat / pmax(lib_size, 1) * med_lib)
  
  ct_ordered <- names(colour_map)[names(colour_map) %in% unique(clust)]
  mean_expr  <- do.call(cbind, lapply(ct_ordered, function(ct) {
    idx <- which(clust == ct)
    colMeans(norm_mat[idx, , drop = FALSE], na.rm = TRUE)
  }))
  colnames(mean_expr) <- ct_ordered
  
  global_mean <- rowMeans(mean_expr, na.rm = TRUE)
  sel_genes   <- unique(unlist(lapply(ct_ordered, function(ct) {
    spec_score <- mean_expr[, ct] - global_mean
    head(names(sort(spec_score, decreasing = TRUE)), top_n_markers)
  })))
  
  if (length(sel_genes) < 2) {
    cat("    \u26A0 Dot plot: not enough marker genes\n")
    return(invisible(NULL))
  }
  
  # Long data frame: mean expression + pct expressing per cluster x gene
  dot_df <- do.call(rbind, lapply(ct_ordered, function(ct) {
    idx     <- which(clust == ct)
    sub_mat <- norm_mat[idx, sel_genes, drop = FALSE]
    data.frame(
      CellType   = ct,
      Gene       = sel_genes,
      MeanExpr   = colMeans(sub_mat, na.rm = TRUE),
      PctExpress = colMeans(sub_mat > 0, na.rm = TRUE) * 100,
      stringsAsFactors = FALSE
    )
  }))
  
  dot_df$CellType <- factor(dot_df$CellType, levels = rev(ct_ordered))
  dot_df$Gene     <- factor(dot_df$Gene,     levels = sel_genes)
  
  title_txt <- sample_plot_title(
    sample_id,
    paste(profile_name, ann_type, "Marker Gene Dot Plot")
  )
  
  p <- ggplot2::ggplot(
    dot_df,
    ggplot2::aes(x = Gene, y = CellType,
                 size  = PctExpress,
                 color = MeanExpr)) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(
      name   = "% Expressing",
      range  = c(0.5, 8),
      limits = c(0, 100)
    ) +
    ggplot2::scale_color_gradientn(
      name    = "Mean log-norm\nexpression",
      colours = c("#2166AC", "#F7F7F7", "#B2182B")
    ) +
    ggplot2::labs(
      title = title_txt,
      subtitle = "Point size reflects the fraction of cells expressing each marker within a cell type.",
      x = "Marker gene",
      y = "Cell Type"
    ) +
    presentation_theme(base_size = 11, x_angle = 45) +
    ggplot2::theme(legend.position = "right")
  
  fname <- paste0(sample_id, "_dotplot_", profile_name, "_", ann_type, ".png")
  save_presentation_plot(
    plot = p,
    filename = file.path(out_dir, fname),
    width = width,
    height = height,
    dpi = dpi
  )
  cat("  \u2713 Dot plot saved:", fname, "\n")
  invisible(p)
}


# Helper: cluster refinement ----------------------------------------------

# refineClusters() is a post-hoc editing tool that operates on the logliks
# matrix. It does NOT accept a per-cell clust vector or align_genes.
#
# Confidence is derived from $logliks via per-cell softmax normalisation
# (see .cluster_mean_confidence). Cluster types whose mean posterior <
# conf_threshold are passed as to_delete; refineClusters() then re-assigns
# those cells to their next-best type via the loglik matrix.
#
# @param gobj              Annotated Giotto object
# @param insitu_result     insitutypeML() result for the profile
# @param counts_mat        cells x genes count matrix (same as annotation)
# @param bg_per_cell       Per-cell background vector
# @param cell_order        Character vector of cell IDs (rownames of counts_mat)
# @param colour_map        Named colour vector from .build_colour_map()
# @param profile_name      Profile name string
# @param sample_id         Sample ID string
# @param results_folder    Path to 07_Annotation/<profile_name>/
# @param conf_threshold    Cluster types with mean posterior < this are deleted
#                          and their cells reassigned (default 0.8)
# @param create_plots      Whether to produce visualisation outputs

refine_annotation <- function(gobj,
                              insitu_result,
                              counts_mat,
                              bg_per_cell,
                              cell_order,
                              colour_map,
                              profile_name,
                              sample_id,
                              results_folder,
                              conf_threshold = 0.8,
                              create_plots   = TRUE) {
  
  cat("\n--- REFINEMENT (mean cluster conf \u2265", conf_threshold, ") ---\n")
  
  refined_folder <- file.path(results_folder, "refined")
  dir.create(refined_folder, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Per-cluster mean confidence -----------------------------------------
  mean_conf_per_type <- tryCatch(
    .cluster_mean_confidence(insitu_result),
    error = function(e) {
      cat("  \u26A0 Could not compute cluster confidence:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(mean_conf_per_type)) return(invisible(gobj))
  
  conf_df <- data.frame(
    cell_type = names(mean_conf_per_type),
    mean_conf = round(mean_conf_per_type, 3),
    n_cells   = as.integer(table(insitu_result$clust)[names(mean_conf_per_type)]),
    stringsAsFactors = FALSE
  )
  conf_df <- conf_df[order(conf_df$mean_conf), ]
  
  cat("  Per-cluster mean confidence (from logliks softmax):\n")
  print(conf_df, row.names = FALSE)
  cat("\n")
  
  readr::write_csv(
    conf_df,
    file.path(refined_folder,
              paste0(sample_id, "_cluster_confidence_scores.csv"))
  )
  
  # 2. Identify types to delete --------------------------------------------
  to_delete <- names(mean_conf_per_type)[
    !is.na(mean_conf_per_type) & mean_conf_per_type < conf_threshold
  ]
  
  if (length(to_delete) == 0) {
    cat("  No cluster types below threshold (", conf_threshold,
        ") \u2014 all types retained, nothing to refine\n")
    return(invisible(gobj))
  }
  
  remaining <- setdiff(colnames(insitu_result$logliks), to_delete)
  if (length(remaining) == 0) {
    cat("  \u26A0 Threshold", conf_threshold,
        "would delete ALL cluster types \u2014 aborting\n")
    cat("  Consider lowering conf_threshold\n")
    return(invisible(gobj))
  }
  
  cat("  Cluster types to delete (mean conf <", conf_threshold, "):\n")
  cat("   ", paste(to_delete, collapse = ", "), "\n")
  cat("  Retained types:", length(remaining), "\n\n")

  # 3. Run refineClusters() ------------------------------------------------
  # refineClusters() requires a matrix with ≥2 columns; a single retained
  # type collapses the logliks subset to a vector and causes a dimension error.
  if (length(remaining) < 2) {
    message("  \u26A0 Skipping refineClusters(): only ", length(remaining),
            " cell type retained after confidence filtering",
            " \u2014 refinement requires \u22652 types")
    return(invisible(gobj))
  }

  cat("  Running InSituType::refineClusters()...\n")

  insitu_refined <- tryCatch({
    InSituType::refineClusters(
      logliks   = insitu_result$logliks,
      to_delete = to_delete,
      counts    = counts_mat,
      neg       = bg_per_cell
    )
  }, error = function(e) {
    cat("  \u26A0 refineClusters() failed:", conditionMessage(e), "\n")
    NULL
  })
  
  if (is.null(insitu_refined)) return(invisible(gobj))
  
  cat("  \u2713 refineClusters() complete\n")
  
  refined_clust    <- insitu_refined$clust
  refined_prob_raw <- insitu_refined$prob
  
  refined_score <- if (is.matrix(refined_prob_raw)) {
    vapply(seq_along(refined_clust), function(i) {
      ct <- refined_clust[i]
      if (is.na(ct) || !ct %in% colnames(refined_prob_raw)) return(NA_real_)
      refined_prob_raw[i, ct]
    }, numeric(1))
  } else {
    as.numeric(refined_prob_raw)
  }
  
  cat("  Cell types after refinement:",
      length(unique(stats::na.omit(refined_clust))), "\n\n")
  
  # 4. Annotation data frame + Giotto metadata -----------------------------
  celltype_col_ref <- paste0("celltype_", profile_name, "_supervised_refined")
  score_col_ref    <- paste0("score_",    profile_name, "_supervised_refined")
  
  annot_ref <- data.frame(
    cell_ID       = cell_order,
    temp_celltype = as.character(refined_clust),
    temp_score    = refined_score,
    stringsAsFactors = FALSE
  )
  names(annot_ref)[2] <- celltype_col_ref
  names(annot_ref)[3] <- score_col_ref
  
  gobj <- addCellMetadata(
    gobject        = gobj,
    new_metadata   = annot_ref,
    by_column      = TRUE,
    column_cell_ID = "cell_ID"
  )
  
  readr::write_csv(
    annot_ref,
    file.path(refined_folder,
              paste0(sample_id, "_supervised_refined_celltypes.csv"))
  )
  
  ref_summary <- dplyr::arrange(
    dplyr::summarise(
      dplyr::group_by(annot_ref, !!dplyr::sym(celltype_col_ref)),
      n_cells    = dplyr::n(),
      mean_score = mean(!!dplyr::sym(score_col_ref), na.rm = TRUE),
      .groups    = "drop"
    ),
    dplyr::desc(n_cells)
  )
  
  readr::write_csv(
    ref_summary,
    file.path(refined_folder,
              paste0(sample_id, "_supervised_refined_summary.csv"))
  )
  
  cat("=== Refined Summary ===\n")
  print(ref_summary)
  cat("\n")
  
  # 5. Colour map ----------------------------------------------------------
  ref_types      <- names(sort(table(stats::na.omit(refined_clust)),
                               decreasing = TRUE))
  colour_map_ref <- colour_map[ref_types]
  missing_types  <- ref_types[is.na(colour_map_ref)]
  
  if (length(missing_types) > 0) {
    n_extra    <- length(missing_types)
    extra_cols <- .annotation_palette(length(colour_map) + n_extra)[
      (length(colour_map) + 1):(length(colour_map) + n_extra)
    ]
    colour_map_ref[missing_types] <- extra_cols
  }
  colour_map_ref <- colour_map_ref[!is.na(colour_map_ref)]
  
  # 6. Visualisations ------------------------------------------------------
  if (create_plots) {
    cat("  Creating refined visualizations...\n")
    
    # Refined — flightpath
    tryCatch({
      plot_custom_flightpath(
        insitu_result = insitu_refined,
        colour_map    = colour_map_ref,
        profile_name  = paste0(profile_name, "_refined"),
        sample_id     = sample_id,
        out_dir       = refined_folder,
        ann_type      = "supervised_refined"
      )
    }, error = function(e) {
      cat("  \u26A0 Flightpath (refined) failed:", conditionMessage(e), "\n")
    })
    
    # Refined — UMAP
    tryCatch({
      clust_named_ref <- stats::setNames(
        as.character(refined_clust), cell_order)
      plot_giotto_umap(
        gobj         = gobj,
        clust_vec    = clust_named_ref,
        colour_map   = colour_map_ref,
        profile_name = paste0(profile_name, "_refined"),
        sample_id    = sample_id,
        out_dir      = refined_folder,
        ann_type     = "supervised_refined"
      )
    }, error = function(e) {
      cat("  \u26A0 UMAP (refined) failed:", conditionMessage(e), "\n")
    })
    
    # Refined — spatial
    tryCatch({
      suppressWarnings({
        spatPlot2D(gobj,
                   cell_color = celltype_col_ref,
                   point_size = 0.5,
                   show_image = FALSE,
                   save_plot  = TRUE,
                   save_param = list(
                     save_name  = paste0(sample_id, "_spatial_",
                                         profile_name, "_supervised_refined"),
                     save_dir   = refined_folder,
                     base_width = 20, base_height = 10))
      })
      cat("  \u2713 Spatial (refined)\n")
    }, error = function(e) {
      cat("  \u26A0 Spatial (refined) failed:", conditionMessage(e), "\n")
    })
    
    # Refined — proportions
    tryCatch({
      .plot_proportions(
        summary_df   = ref_summary,
        celltype_col = celltype_col_ref,
        colour_map   = colour_map_ref,
        title        = paste0(sample_id, " - ", profile_name,
                              " (refined, conf \u2265 ", conf_threshold, ")"),
        out_path     = file.path(refined_folder,
                                 paste0(sample_id,
                                        "_proportions_supervised_refined.png"))
      )
      cat("  \u2713 Proportions (refined)\n")
    }, error = function(e) {
      cat("  \u26A0 Proportions (refined) failed:", conditionMessage(e), "\n")
    })
    
    # Refined — heatmap
    tryCatch({
      .plot_annotation_heatmap(
        counts_mat   = counts_mat,
        clust_vec    = stats::setNames(
          as.character(refined_clust), cell_order),
        colour_map   = colour_map_ref,
        profile_name = paste0(profile_name, "_refined"),
        sample_id    = sample_id,
        out_dir      = refined_folder,
        ann_type     = "supervised_refined"
      )
    }, error = function(e) {
      cat("  \u26A0 Heatmap (refined) failed:", conditionMessage(e), "\n")
    })
    
    # Refined — dot plot
    tryCatch({
      .plot_annotation_dotplot(
        counts_mat   = counts_mat,
        clust_vec    = stats::setNames(
          as.character(refined_clust), cell_order),
        colour_map   = colour_map_ref,
        profile_name = paste0(profile_name, "_refined"),
        sample_id    = sample_id,
        out_dir      = refined_folder,
        ann_type     = "supervised_refined"
      )
    }, error = function(e) {
      cat("  \u26A0 Dot plot (refined) failed:", conditionMessage(e), "\n")
    })
    
    cat("  \u2713 Refined visualizations complete\n\n")
  }
  
  cat("\u2713 Refinement complete for", profile_name, "\n\n")
  invisible(gobj)
}


# load_reference_profile --------------------------------------------------

load_reference_profile <- function(profile_config) {
  
  if (profile_config$type == "url") {
    cat("  Downloading from URL...\n")
    ref_data <- read_profile_csv(url = profile_config$source)
    
  } else if (profile_config$type == "hca") {
    cat("  Downloading HCA profile via SpatialDecon...\n")
    ref_data <- SpatialDecon::download_profile_matrix(
      species    = profile_config$species,
      age_group  = profile_config$age_group,
      matrixname = profile_config$matrixname
    )
    ref_data <- as.matrix(ref_data)
    
  } else if (profile_config$type == "file") {
    cat("  Loading from local file...\n")
    ref_data <- read_profile_csv(url = profile_config$source)
    
  } else {
    stop("Unknown profile type: ", profile_config$type)
  }
  
  return(ref_data)
}


# annotate_cells — main pipeline function ---------------------------------

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

.profile_set_has_immune_reference <- function(profiles) {
  if (is.null(profiles) || length(profiles) == 0) {
    return(FALSE)
  }
  
  profile_text <- vapply(profiles, function(profile) {
    paste(
      profile$name %||% "",
      profile$type %||% "",
      profile$source %||% "",
      profile$matrixname %||% "",
      collapse = " "
    )
  }, character(1))
  
  any(grepl("immune|pbmc|lymph|lupus|bcell|tcell|myeloid", profile_text, ignore.case = TRUE))
}

annotate_cells <- function(gobj,
                           sample_id,
                           output_dir,
                           profiles         = NULL,
                           default_profile  = NULL,
                           align_genes      = TRUE,
                           n_clusts_semi    = 3,
                           n_starts         = 10,
                           cohort_column    = "leiden_clust",
                           min_gene_overlap = 100,
                           create_plots     = TRUE,
                           conf_threshold   = NULL,
                           save_object      = TRUE) {
  
  cat("\n========================================\n")
  cat("STEP 07: Cell Type Annotation\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    cat("Loading Giotto object from:", gobj, "\n")
    gobj <- loadGiotto(gobj)
    cat("\u2713 Loaded\n\n")
  }
  
  results_folder <- file.path(output_dir, "07_Annotation")
  dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
  
  # Count matrix + background ---------------------------------------------
  cat("Preparing data for InSituType...\n")
  
  counts_raw <- getExpression(
    gobject   = gobj,
    feat_type = "rna",
    spat_unit = "cell",
    values    = "raw",
    output    = "matrix"
  )
  
  counts_mat <- t(as.matrix(counts_raw))
  cat("  Count matrix:", nrow(counts_mat), "cells x", ncol(counts_mat), "genes\n")
  
  all_feats <- rownames(counts_raw)
  neg_genes <- grep("^Negative", all_feats, value = TRUE, ignore.case = TRUE)
  
  if (length(neg_genes) > 0) {
    cat("  Using", length(neg_genes), "negative control genes for background\n")
    neg_counts_mat <- counts_raw[neg_genes, , drop = FALSE]
    bg_per_cell    <- colMeans(as.matrix(neg_counts_mat))
    cat("  Background - Mean:", round(mean(bg_per_cell), 3),
        "Median:", round(median(bg_per_cell), 3), "\n")
  } else {
    cat("  No negative probes, using 1% of mean expression\n")
    bg_per_cell <- Matrix::rowMeans(counts_mat) * 0.01
  }
  
  names(bg_per_cell) <- rownames(counts_mat)
  cat("\u2713 Data prepared\n\n")
  
  cell_order <- rownames(counts_mat)
  
  # Default profile index -------------------------------------------------
  if (!is.null(default_profile)) {
    default_idx <- which(sapply(profiles, function(p) p$name == default_profile))
    if (length(default_idx) == 0) {
      cat("\u26A0 Default profile '", default_profile, "' not found, using first\n")
      default_idx <- 1
    }
  } else {
    default_idx <- 1
  }
  
  cat("Default profile:", profiles[[default_idx]]$name, "\n\n")
  
  if (!.profile_set_has_immune_reference(profiles)) {
    cat(
      "⚠ The configured annotation profiles look kidney-centric and may coarsely label infiltrating immune cells.\n",
      "  For lupus nephritis work, consider adding a PBMC, immune atlas, or custom LN immune reference.\n\n",
      sep = ""
    )
  }
  
  # Cohort vector aligned to cell_order -----------------------------------
  metadata_dt <- as.data.frame(pDataDT(gobj))
  
  if (cohort_column %in% names(metadata_dt)) {
    
    if (!"cell_ID" %in% names(metadata_dt))
      stop("pDataDT() has no 'cell_ID' column")
    
    rownames(metadata_dt) <- metadata_dt$cell_ID
    shared_cells          <- intersect(cell_order, rownames(metadata_dt))
    missing_cells         <- setdiff(cell_order,  rownames(metadata_dt))
    
    if (length(missing_cells) > 0)
      cat("\u26A0", length(missing_cells),
          "cells in counts_mat have no metadata row\n")
    
    cohort_vec               <- rep(NA_character_, length(cell_order))
    names(cohort_vec)        <- cell_order
    cohort_vec[shared_cells] <- as.character(
      metadata_dt[shared_cells, cohort_column])
    
    cat("Cohort column '", cohort_column,
        "' aligned to counts_mat \u2014 will be used for semi-supervised\n")
    cat("  Cohort levels:", length(unique(stats::na.omit(cohort_vec))), "\n\n")
    
  } else {
    cohort_vec <- NULL
    cat("\u26A0 Cohort column '", cohort_column,
        "' not found \u2014 semi-supervised will run without cohort vector\n\n")
  }
  
  insitu_results_store <- list()
  
  # Loop over profiles ----------------------------------------------------
  for (i in seq_along(profiles)) {
    
    profile_config <- profiles[[i]]
    profile_name   <- profile_config$name
    is_default     <- (i == default_idx)
    
    cat("================================================================================\n")
    cat("Running annotation with:", profile_name)
    if (is_default) cat(" [DEFAULT]")
    cat("\n================================================================================\n\n")
    
    profile_folder <- file.path(results_folder, profile_name)
    dir.create(profile_folder, recursive = TRUE, showWarnings = FALSE)
    
    tryCatch({
      
      cat("Loading reference profile...\n")
      ref_profiles <- load_reference_profile(profile_config)
      
      cat("\u2713 Profile loaded\n")
      cat("  Cell types:", ncol(ref_profiles), "\n")
      cat("  Genes:", nrow(ref_profiles), "\n")
      
      common_genes <- intersect(rownames(ref_profiles), colnames(counts_mat))
      cat("  Overlapping genes:", length(common_genes), "\n")
      
      if (length(common_genes) < min_gene_overlap) {
        cat("\u2717 Insufficient overlap (", length(common_genes), " < ",
            min_gene_overlap, ")\n  Skipping", profile_name, "\n\n")
        next
      }
      cat("  Sufficient overlap - proceeding\n\n")
      
      # BLOCK A: SUPERVISED ------------------------------------------------
      cat("--- SUPERVISED ---\n")
      cat("Running InSituType (supervised)...\n")
      
      insitu_supervised <- InSituType::insitutypeML(
        x                  = counts_mat,
        neg                = bg_per_cell,
        reference_profiles = ref_profiles,
        align_genes        = align_genes
      )
      
      cat("\u2713 Supervised complete\n")
      cat("  Cell types found:", length(unique(insitu_supervised$clust)), "\n\n")
      
      colour_map_sup   <- .build_colour_map(insitu_supervised$clust)
      celltype_col_sup <- paste0("celltype_", profile_name, "_supervised")
      score_col_sup    <- paste0("score_",    profile_name, "_supervised")
      
      annot_sup <- data.frame(
        cell_ID       = cell_order,
        temp_celltype = as.character(insitu_supervised$clust),
        temp_score    = insitu_supervised$prob,
        stringsAsFactors = FALSE
      )
      names(annot_sup)[2] <- celltype_col_sup
      names(annot_sup)[3] <- score_col_sup
      
      if (is_default) {
        annot_sup$celltype       <- annot_sup[[celltype_col_sup]]
        annot_sup$celltype_score <- annot_sup[[score_col_sup]]
      }
      
      cat("Adding supervised annotations to Giotto object...\n")
      cat("  Dimensions:", nrow(annot_sup), "x", ncol(annot_sup), "\n")
      cat("  Columns:", paste(names(annot_sup), collapse = ", "), "\n")
      
      gobj <- addCellMetadata(
        gobject        = gobj,
        new_metadata   = annot_sup,
        by_column      = TRUE,
        column_cell_ID = "cell_ID"
      )
      
      readr::write_csv(
        annot_sup,
        file.path(profile_folder, paste0(sample_id, "_supervised_celltypes.csv"))
      )
      cat("\u2713 Supervised annotations added and saved\n\n")
      
      sup_summary <- dplyr::arrange(
        dplyr::summarise(
          dplyr::group_by(annot_sup, !!dplyr::sym(celltype_col_sup)),
          n_cells    = dplyr::n(),
          mean_score = mean(!!dplyr::sym(score_col_sup), na.rm = TRUE),
          .groups    = "drop"
        ),
        dplyr::desc(n_cells)
      )
      
      readr::write_csv(
        sup_summary,
        file.path(profile_folder, paste0(sample_id, "_supervised_summary.csv"))
      )
      cat("=== Supervised Summary ===\n")
      print(sup_summary)
      cat("\n")
      
      insitu_results_store[[profile_name]] <- list(
        supervised     = insitu_supervised,
        colour_map     = colour_map_sup,
        counts_mat     = counts_mat,
        bg_per_cell    = bg_per_cell,
        cell_order     = cell_order,
        profile_folder = profile_folder
      )
      
      # BLOCK B: SEMI-SUPERVISED -------------------------------------------
      insitu_semi       <- NULL
      annot_semi        <- NULL
      semi_summary      <- NULL
      colour_map_semi   <- NULL
      celltype_col_semi <- paste0("celltype_", profile_name, "_semi")
      score_col_semi    <- paste0("score_",    profile_name, "_semi")
      
      if (n_clusts_semi > 0) {
        
        cat("--- SEMI-SUPERVISED ---\n")
        cat("Running InSituType (semi-supervised)...\n")
        cat("  Novel clusters:", n_clusts_semi, "\n")
        cat("  Starts:", n_starts, "\n")
        if (!is.null(cohort_vec)) {
          cat("  Cohort column:", cohort_column, "\n")
        } else {
          cat("  Cohort: none\n")
        }
        cat("\n")
        
        insitu_semi <- tryCatch({
          InSituType::insitutype(
            x                  = counts_mat,
            neg                = bg_per_cell,
            reference_profiles = ref_profiles,
            n_clusts           = n_clusts_semi,
            cohort             = cohort_vec,
            align_genes        = align_genes,
            n_starts           = n_starts
          )
        }, error = function(e) {
          err_msg <- conditionMessage(e)
          n_cells <- nrow(counts_mat)
          message(
            "\u26A0 Semi-supervised annotation failed for profile '", profile_name, "'",
            " (sample: ", sample_id, ", cells: ", n_cells, ").\n",
            "  Error: ", err_msg, "\n",
            if (grepl("anchor", err_msg, ignore.case = TRUE)) {
              paste0(
                "  Likely cause: sample is too small (", n_cells, " cells) or",
                " confidence scores are too uniform to select anchor cells.\n",
                "  Tip: try lowering the anchor confidence threshold,",
                " increasing n_clusts_semi, or using refinement = FALSE.\n"
              )
            } else "",
            "  Continuing with supervised-only annotations."
          )

          # Fallback: retry with refinement = FALSE to bypass anchor selection
          message("  \u21B3 Retrying semi-supervised with refinement = FALSE ...")
          fallback_result <- tryCatch({
            InSituType::insitutype(
              x                  = counts_mat,
              neg                = bg_per_cell,
              reference_profiles = ref_profiles,
              n_clusts           = n_clusts_semi,
              cohort             = cohort_vec,
              align_genes        = align_genes,
              n_starts           = n_starts,
              refinement         = FALSE
            )
          }, error = function(e2) {
            message("  \u21B3 Retry also failed: ", conditionMessage(e2))
            message("  \u21B3 Using supervised-only annotations.")
            NULL
          })
          fallback_result
        })
        
        if (!is.null(insitu_semi)) {
          cat("\u2713 Semi-supervised complete\n")
          cat("  Cell types found:", length(unique(insitu_semi$clust)), "\n\n")
          
          colour_map_semi <- .build_colour_map(insitu_semi$clust)
          
          annot_semi <- data.frame(
            cell_ID       = cell_order,
            temp_celltype = as.character(insitu_semi$clust),
            temp_score    = insitu_semi$prob,
            stringsAsFactors = FALSE
          )
          names(annot_semi)[2] <- celltype_col_semi
          names(annot_semi)[3] <- score_col_semi
          
          gobj <- addCellMetadata(
            gobject        = gobj,
            new_metadata   = annot_semi,
            by_column      = TRUE,
            column_cell_ID = "cell_ID"
          )
          
          readr::write_csv(
            annot_semi,
            file.path(profile_folder, paste0(sample_id, "_semi_celltypes.csv"))
          )
          cat("\u2713 Semi-supervised annotations added and saved\n\n")
          
          semi_summary <- dplyr::arrange(
            dplyr::summarise(
              dplyr::group_by(annot_semi, !!dplyr::sym(celltype_col_semi)),
              n_cells    = dplyr::n(),
              mean_score = mean(!!dplyr::sym(score_col_semi), na.rm = TRUE),
              .groups    = "drop"
            ),
            dplyr::desc(n_cells)
          )
          
          readr::write_csv(
            semi_summary,
            file.path(profile_folder, paste0(sample_id, "_semi_summary.csv"))
          )
          cat("=== Semi-Supervised Summary ===\n")
          print(semi_summary)
          cat("\n")
          
          insitu_results_store[[profile_name]]$semi            <- insitu_semi
          insitu_results_store[[profile_name]]$colour_map_semi <- colour_map_semi
        }
        
      } else {
        cat("Semi-supervised skipped (n_clusts_semi = 0)\n\n")
      }
      
      # BLOCK C: VISUALISATIONS --------------------------------------------
      if (create_plots) {
        cat("Creating visualizations...\n")
        
        # Supervised — flightpath
        tryCatch({
          plot_custom_flightpath(
            insitu_result = insitu_supervised,
            colour_map    = colour_map_sup,
            profile_name  = profile_name,
            sample_id     = sample_id,
            out_dir       = profile_folder,
            ann_type      = "supervised"
          )
        }, error = function(e) {
          cat("  \u26A0 Flightpath (supervised) failed:", conditionMessage(e), "\n")
        })
        
        # Supervised — UMAP
        tryCatch({
          clust_named_sup <- stats::setNames(
            as.character(insitu_supervised$clust), cell_order)
          plot_giotto_umap(
            gobj         = gobj,
            clust_vec    = clust_named_sup,
            colour_map   = colour_map_sup,
            profile_name = profile_name,
            sample_id    = sample_id,
            out_dir      = profile_folder,
            ann_type     = "supervised"
          )
        }, error = function(e) {
          cat("  \u26A0 UMAP (supervised) failed:", conditionMessage(e), "\n")
        })
        
        # Supervised — spatial
        tryCatch({
          suppressWarnings({
            spatPlot2D(gobj,
                       cell_color = celltype_col_sup,
                       point_size = 0.5,
                       show_image = FALSE,
                       save_plot  = TRUE,
                       save_param = list(
                         save_name  = paste0(sample_id, "_spatial_",
                                             profile_name, "_supervised"),
                         save_dir   = profile_folder,
                         base_width = 20, base_height = 10))
          })
          cat("  \u2713 Spatial (supervised)\n")
        }, error = function(e) {
          cat("  \u26A0 Spatial (supervised) failed:", conditionMessage(e), "\n")
        })
        
        # Supervised — proportions
        tryCatch({
          .plot_proportions(
            summary_df   = sup_summary,
            celltype_col = celltype_col_sup,
            colour_map   = colour_map_sup,
            title        = paste0(sample_id, " - ", profile_name,
                                  " (supervised)"),
            out_path     = file.path(profile_folder,
                                     paste0(sample_id,
                                            "_proportions_supervised.png"))
          )
          cat("  \u2713 Proportions (supervised)\n")
        }, error = function(e) {
          cat("  \u26A0 Proportions (supervised) failed:", conditionMessage(e), "\n")
        })
        
        # Supervised — heatmap
        tryCatch({
          .plot_annotation_heatmap(
            counts_mat   = counts_mat,
            clust_vec    = stats::setNames(
              as.character(insitu_supervised$clust), cell_order),
            colour_map   = colour_map_sup,
            profile_name = profile_name,
            sample_id    = sample_id,
            out_dir      = profile_folder,
            ann_type     = "supervised"
          )
        }, error = function(e) {
          cat("  \u26A0 Heatmap (supervised) failed:", conditionMessage(e), "\n")
        })
        
        # Supervised — dot plot
        tryCatch({
          .plot_annotation_dotplot(
            counts_mat   = counts_mat,
            clust_vec    = stats::setNames(
              as.character(insitu_supervised$clust), cell_order),
            colour_map   = colour_map_sup,
            profile_name = profile_name,
            sample_id    = sample_id,
            out_dir      = profile_folder,
            ann_type     = "supervised"
          )
        }, error = function(e) {
          cat("  \u26A0 Dot plot (supervised) failed:", conditionMessage(e), "\n")
        })
        
        # Semi-supervised plots (only if semi succeeded)
        if (!is.null(insitu_semi) && !is.null(annot_semi)) {
          
          # Semi — flightpath
          tryCatch({
            plot_custom_flightpath(
              insitu_result = insitu_semi,
              colour_map    = colour_map_semi,
              profile_name  = profile_name,
              sample_id     = sample_id,
              out_dir       = profile_folder,
              ann_type      = "semi"
            )
          }, error = function(e) {
            cat("  \u26A0 Flightpath (semi) failed:", conditionMessage(e), "\n")
          })
          
          # Semi — UMAP
          tryCatch({
            clust_named_semi <- stats::setNames(
              as.character(insitu_semi$clust), cell_order)
            plot_giotto_umap(
              gobj         = gobj,
              clust_vec    = clust_named_semi,
              colour_map   = colour_map_semi,
              profile_name = profile_name,
              sample_id    = sample_id,
              out_dir      = profile_folder,
              ann_type     = "semi"
            )
          }, error = function(e) {
            cat("  \u26A0 UMAP (semi) failed:", conditionMessage(e), "\n")
          })
          
          # Semi — spatial
          tryCatch({
            suppressWarnings({
              spatPlot2D(gobj,
                         cell_color = celltype_col_semi,
                         point_size = 0.5,
                         show_image = FALSE,
                         save_plot  = TRUE,
                         save_param = list(
                           save_name  = paste0(sample_id, "_spatial_",
                                               profile_name, "_semi"),
                           save_dir   = profile_folder,
                           base_width = 20, base_height = 10))
            })
            cat("  \u2713 Spatial (semi)\n")
          }, error = function(e) {
            cat("  \u26A0 Spatial (semi) failed:", conditionMessage(e), "\n")
          })
          
          # Semi — proportions
          if (!is.null(semi_summary)) {
            tryCatch({
              .plot_proportions(
                summary_df   = semi_summary,
                celltype_col = celltype_col_semi,
                colour_map   = colour_map_semi,
                title        = paste0(sample_id, " - ", profile_name,
                                      " (semi-supervised)"),
                out_path     = file.path(profile_folder,
                                         paste0(sample_id,
                                                "_proportions_semi.png"))
              )
              cat("  \u2713 Proportions (semi)\n")
            }, error = function(e) {
              cat("  \u26A0 Proportions (semi) failed:", conditionMessage(e), "\n")
            })
          }
          
          # Semi — heatmap
          tryCatch({
            .plot_annotation_heatmap(
              counts_mat   = counts_mat,
              clust_vec    = stats::setNames(
                as.character(insitu_semi$clust), cell_order),
              colour_map   = colour_map_semi,
              profile_name = profile_name,
              sample_id    = sample_id,
              out_dir      = profile_folder,
              ann_type     = "semi"
            )
          }, error = function(e) {
            cat("  \u26A0 Heatmap (semi) failed:", conditionMessage(e), "\n")
          })
          
          # Semi — dot plot
          tryCatch({
            .plot_annotation_dotplot(
              counts_mat   = counts_mat,
              clust_vec    = stats::setNames(
                as.character(insitu_semi$clust), cell_order),
              colour_map   = colour_map_semi,
              profile_name = profile_name,
              sample_id    = sample_id,
              out_dir      = profile_folder,
              ann_type     = "semi"
            )
          }, error = function(e) {
            cat("  \u26A0 Dot plot (semi) failed:", conditionMessage(e), "\n")
          })
          
        }  # end semi-supervised plots block
        
        cat("\u2713 Visualizations complete\n\n")
      }  # end if (create_plots)
      
      # BLOCK D: OPTIONAL INLINE REFINEMENT --------------------------------
      if (!is.null(conf_threshold) && conf_threshold > 0) {
        gobj <- refine_annotation(
          gobj           = gobj,
          insitu_result  = insitu_supervised,
          counts_mat     = counts_mat,
          bg_per_cell    = bg_per_cell,
          cell_order     = cell_order,
          colour_map     = colour_map_sup,
          profile_name   = profile_name,
          sample_id      = sample_id,
          results_folder = profile_folder,
          conf_threshold = conf_threshold,
          create_plots   = create_plots
        )
      }
      
      cat("\u2713 Complete:", profile_name, "\n\n")
      
    }, error = function(e) {
      cat("\u2717 Error with", profile_name, ":\n  ", conditionMessage(e), "\n\n")
    })
  }  # end for loop over profiles
  
  cat("\u2713 All annotations complete for", sample_id, "\n\n")

  # ── Auto-select best annotation ──────────────────────────────────────────
  cat("Selecting best annotation based on confidence and LN composition...\n")
  annotation_selection <- tryCatch(
    select_best_annotation(
      gobj           = gobj,
      annotation_dir = file.path(output_dir, "07_Annotation"),
      sample_id      = sample_id,
      conf_threshold = if (!is.null(conf_threshold) && conf_threshold > 0)
                         conf_threshold else 0.8
    ),
    error = function(e) {
      cat("  \u26A0 Annotation selection failed:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (!is.null(annotation_selection)) {
    best_col   <- annotation_selection$selected_col
    best_score <- sub("^celltype_", "score_", best_col)
    meta_now   <- as.data.frame(pDataDT(gobj))

    if (best_col %in% names(meta_now) && best_score %in% names(meta_now)) {
      update_df <- data.frame(
        cell_ID        = meta_now$cell_ID,
        celltype       = as.character(meta_now[[best_col]]),
        celltype_score = as.numeric(meta_now[[best_score]]),
        stringsAsFactors = FALSE
      )
      gobj <- addCellMetadata(
        gobject        = gobj,
        new_metadata   = update_df,
        by_column      = TRUE,
        column_cell_ID = "cell_ID"
      )
      cat("  \u2713 Generic 'celltype' column updated to selected annotation\n\n")
    }
  }

  # Save annotated Giotto object ------------------------------------------
  if (save_object) {
    cat("Saving annotated Giotto object to:",
        file.path(output_dir, "Giotto_Object_Annotated"), "\n")
    tryCatch({
      saveGiotto(gobj,
                 dir        = output_dir,
                 foldername = "Giotto_Object_Annotated",
                 overwrite  = TRUE)
      cat("\u2713 Giotto object saved\n\n")
    }, error = function(e) {
      cat("\u26A0 saveGiotto() failed:", conditionMessage(e), "\n\n")
    })
  }
  
  attr(gobj, "insitu_results") <- insitu_results_store
  return(gobj)
}


# Command-line interface --------------------------------------------------

if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 4) {
    script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
    script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
    bootstrap_script <- file.path(script_dir, "Helper_Scripts", "Script_Bootstrap.R")
    if (file.exists(bootstrap_script)) {
      source(bootstrap_script, local = .GlobalEnv)
      bootstrap_pipeline_environment(script_dir, load_pipeline_utils = FALSE, verbose = FALSE)
    }
    
    sample_id   <- args[1]
    input_path  <- args[2]
    output_dir  <- args[3]
    config_file <- if (length(args) >= 4) args[4] else NULL
    
    if (!is.null(config_file) && file.exists(config_file)) {
      config <- yaml::read_yaml(config_file)
      
      helper_path <- file.path(config$paths$scripts_dir,
                               "Helper_Scripts/Helper_Functions.R")
      if (file.exists(helper_path)) source(helper_path)
      
      profiles         <- config$parameters$annotation$profiles
      default_profile  <- config$parameters$annotation$default_profile
      align_genes      <- config$parameters$annotation$align_genes
      n_starts         <- config$parameters$annotation$n_starts
      n_clusts_semi    <- config$parameters$annotation$n_clusts_semi  %||% 3L
      cohort_column    <- config$parameters$annotation$cohort_column
      min_gene_overlap <- config$parameters$annotation$min_gene_overlap %||% 100
      conf_threshold   <- config$parameters$annotation$conf_threshold  %||% NULL
    } else {
      stop("Config required")
    }
    
    gobj <- annotate_cells(
      gobj             = input_path,
      sample_id        = sample_id,
      output_dir       = output_dir,
      profiles         = profiles,
      default_profile  = default_profile,
      align_genes      = align_genes,
      n_starts         = n_starts,
      n_clusts_semi    = n_clusts_semi,
      cohort_column    = cohort_column,
      min_gene_overlap = min_gene_overlap,
      conf_threshold   = conf_threshold,
      save_object      = TRUE
    )
  } else {
    stop("Usage: Rscript 07_Annotation.R <sample_id> <input_path> <output_dir> <config_file>")
  }
}

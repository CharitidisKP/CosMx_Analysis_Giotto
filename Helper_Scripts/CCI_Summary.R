cci_summary_get_expression <- function(gobj, values = "normalized", output = "matrix") {
  if (exists(".giotto_get_expression", mode = "function", inherits = TRUE)) {
    return(get(".giotto_get_expression", mode = "function", inherits = TRUE)(gobj, values = values, output = output))
  }
  
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

cci_summary_muffle_plot_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      warning_message <- conditionMessage(w)
      if (grepl("aes_string\\(", warning_message, fixed = FALSE) ||
          grepl("Viewport has zero dimension", warning_message, ignore.case = TRUE) ||
          grepl("getDistinctColors", warning_message, ignore.case = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

cci_summary_pick_numeric_column <- function(df, preferred = character()) {
  numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
  if (length(numeric_cols) == 0) {
    return(NULL)
  }
  
  preferred_match <- intersect(preferred, numeric_cols)
  if (length(preferred_match) > 0) {
    return(preferred_match[1])
  }
  
  numeric_cols[1]
}

cci_summary_read_csv <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- suppressMessages(
    readr::read_csv(
      path,
      show_col_types = FALSE,
      name_repair = "unique"
    )
  )
  
  if (ncol(df) > 0) {
    first_col <- names(df)[1]
    if (identical(first_col, "...1")) {
      first_values <- df[[1]]
      looks_like_row_index <- is.numeric(first_values) &&
        length(first_values) == nrow(df) &&
        all(first_values == seq_len(nrow(df)))
      if (looks_like_row_index || !("module" %in% names(df) && "gene" %in% names(df))) {
        df <- df[, -1, drop = FALSE]
      }
    }
  }
  
  tibble::as_tibble(df, .name_repair = "unique")
}

create_cci_summary <- function(gobj,
                               sample_id,
                               output_dir,
                               cci_results = NULL,
                               top_n = 15,
                               top_svg_plots = 6) {
  results_dir <- file.path(output_dir, "10_CCI_Analysis")
  summary_dir <- file.path(results_dir, "Summary")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("Creating CCI summary outputs...\n")
  
  summary_rows <- list()
  add_summary_row <- function(section, metric, value) {
    summary_rows[[length(summary_rows) + 1]] <<- data.frame(
      section = section,
      metric = metric,
      value = as.character(value),
      stringsAsFactors = FALSE
    )
  }
  
  section_status <- attr(cci_results, "section_status")
  section_messages <- attr(cci_results, "section_messages")
  if (!is.null(section_status)) {
    status_df <- data.frame(
      section = names(section_status),
      status = unname(section_status),
      message = vapply(
        names(section_status),
        function(section_name) {
          msg <- if (is.null(section_messages)) NULL else section_messages[[section_name]]
          if (is.null(msg) || (length(msg) == 1 && is.na(msg))) "" else as.character(msg)
        },
        character(1)
      ),
      stringsAsFactors = FALSE
    )
    readr::write_csv(
      status_df,
      file.path(summary_dir, paste0(sample_id, "_cci_section_status.csv"))
    )
  }
  
  insitucor_rds <- file.path(results_dir, "insitucor", paste0(sample_id, "_insitucor.rds"))
  insitucor_modules_csv <- file.path(results_dir, "insitucor", paste0(sample_id, "_insitucor_modules.csv"))
  insitucor_res <- NULL
  if (!is.null(cci_results) && !is.null(cci_results$insitucor)) {
    insitucor_res <- cci_results$insitucor
  } else if (file.exists(insitucor_rds)) {
    insitucor_res <- tryCatch(readRDS(insitucor_rds), error = function(e) NULL)
  }
  insitucor_modules <- cci_summary_read_csv(insitucor_modules_csv)
  if (is.null(insitucor_modules) && !is.null(insitucor_res$modules)) {
    insitucor_modules <- tibble::as_tibble(insitucor_res$modules)
  }
  
  if (!is.null(insitucor_modules) && nrow(insitucor_modules) > 0 &&
      all(c("module", "gene", "weight") %in% colnames(insitucor_modules))) {
    module_stats <- insitucor_modules %>%
      dplyr::group_by(module) %>%
      dplyr::summarise(
        n_genes = dplyr::n(),
        mean_weight = mean(weight, na.rm = TRUE),
        max_weight = max(weight, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(n_genes), module)
    
    readr::write_csv(
      module_stats,
      file.path(summary_dir, paste0(sample_id, "_insitucor_module_stats.csv"))
    )
    
    add_summary_row("insitucor", "modules_detected", nrow(module_stats))
    add_summary_row("insitucor", "genes_in_modules", nrow(insitucor_modules))
    add_summary_row("insitucor", "largest_module_size", max(module_stats$n_genes, na.rm = TRUE))
    
    top_modules <- module_stats %>%
      dplyr::slice_head(n = min(top_n, nrow(module_stats)))
    
    p_modules <- ggplot2::ggplot(
      top_modules,
      ggplot2::aes(x = stats::reorder(module, n_genes), y = n_genes, fill = mean_weight)
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_viridis_c(option = "C", name = "Mean weight") +
      ggplot2::labs(
        title = "Top InSituCor modules by gene count",
        x = "Module",
        y = "Number of genes"
      ) +
      ggplot2::theme_minimal(base_size = 11)
    
    ggplot2::ggsave(
      filename = file.path(summary_dir, paste0(sample_id, "_insitucor_module_sizes.png")),
      plot = p_modules,
      width = 10,
      height = 6,
      dpi = 150
    )
    
    if (!is.null(insitucor_res$celltypeinvolvement)) {
      involvement <- as.matrix(insitucor_res$celltypeinvolvement)
      if (length(dim(involvement)) == 2 && nrow(involvement) > 0 && ncol(involvement) > 0) {
        top_module_names <- top_modules$module
        if (!is.null(colnames(involvement)) && all(top_module_names %in% colnames(involvement))) {
          involvement <- involvement[, top_module_names, drop = FALSE]
        } else if (!is.null(rownames(involvement)) && all(top_module_names %in% rownames(involvement))) {
          involvement <- involvement[top_module_names, , drop = FALSE]
          involvement <- t(involvement)
        }
        
        keep_celltypes <- order(apply(involvement, 1, max, na.rm = TRUE), decreasing = TRUE)
        keep_celltypes <- keep_celltypes[seq_len(min(length(keep_celltypes), top_n))]
        involvement_plot <- involvement[keep_celltypes, , drop = FALSE]
        
        pheatmap::pheatmap(
          involvement_plot,
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          main = "InSituCor cell type involvement",
          filename = file.path(summary_dir, paste0(sample_id, "_insitucor_celltype_involvement.png")),
          width = 9,
          height = 7
        )
      }
    }
  }
  
  nnsvg_csv <- file.path(results_dir, "svg", paste0(sample_id, "_nnSVG_results.csv"))
  nnsvg_df <- NULL
  if (!is.null(cci_results) && !is.null(cci_results$nnsvg)) {
    nnsvg_df <- tibble::as_tibble(cci_results$nnsvg, rownames = "gene")
  } else if (file.exists(nnsvg_csv)) {
    nnsvg_df <- read.csv(nnsvg_csv, row.names = 1, check.names = FALSE)
    nnsvg_df$gene <- rownames(nnsvg_df)
    nnsvg_df <- tibble::as_tibble(nnsvg_df)
  }
  
  if (!is.null(nnsvg_df) && nrow(nnsvg_df) > 0) {
    if (!"gene" %in% colnames(nnsvg_df)) {
      nnsvg_df$gene <- rownames(as.data.frame(nnsvg_df))
    }
    order_col <- if ("rank" %in% colnames(nnsvg_df)) "rank" else cci_summary_pick_numeric_column(nnsvg_df)
    if (!is.null(order_col)) {
      nnsvg_df <- nnsvg_df[order(nnsvg_df[[order_col]], decreasing = identical(order_col, "rank") == FALSE), , drop = FALSE]
      if (identical(order_col, "rank")) {
        nnsvg_df <- nnsvg_df[order(nnsvg_df$rank), , drop = FALSE]
      }
    }
    
    top_svg_df <- nnsvg_df %>% dplyr::slice_head(n = min(top_n, nrow(nnsvg_df)))
    metric_col <- cci_summary_pick_numeric_column(
      top_svg_df,
      preferred = c("LR_stat", "prop_sv", "loglik", "variance", "rank")
    )
    if (identical(metric_col, "rank")) {
      top_svg_df$plot_value <- rev(seq_len(nrow(top_svg_df)))
      y_label <- "Top-rank ordering"
    } else if (!is.null(metric_col)) {
      top_svg_df$plot_value <- top_svg_df[[metric_col]]
      y_label <- metric_col
    } else {
      top_svg_df$plot_value <- rev(seq_len(nrow(top_svg_df)))
      y_label <- "Top-rank ordering"
    }
    
    add_summary_row("nnsvg", "genes_retained", nrow(nnsvg_df))
    add_summary_row("nnsvg", "top_gene", top_svg_df$gene[1])
    add_summary_row("nnsvg", "top_n_visualized", nrow(top_svg_df))
    
    readr::write_csv(
      top_svg_df,
      file.path(summary_dir, paste0(sample_id, "_nnsvg_top_genes.csv"))
    )
    
    p_svg <- ggplot2::ggplot(
      top_svg_df,
      ggplot2::aes(x = stats::reorder(gene, plot_value), y = plot_value)
    ) +
      ggplot2::geom_col(fill = "#2C7FB8") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Top nnSVG genes",
        x = "Gene",
        y = y_label
      ) +
      ggplot2::theme_minimal(base_size = 11)
    
    ggplot2::ggsave(
      filename = file.path(summary_dir, paste0(sample_id, "_nnsvg_top_genes.png")),
      plot = p_svg,
      width = 10,
      height = 6,
      dpi = 150
    )
    
    top_svg_genes <- head(top_svg_df$gene, top_svg_plots)
    if (!is.null(gobj) && length(top_svg_genes) > 0) {
      for (gene in top_svg_genes) {
        tryCatch({
          cci_summary_muffle_plot_warnings(
            spatFeatPlot2D(
              gobject = gobj,
              expression_values = "normalized",
              feats = gene,
              point_size = 0.5,
              show_image = FALSE,
              save_plot = TRUE,
              save_param = list(
                save_name = paste0(sample_id, "_nnsvg_spatial_", gene),
                save_dir = summary_dir,
                base_width = 12,
                base_height = 10
              )
            )
          )
        }, error = function(e) {
          NULL
        })
      }
    }
  }
  
  misty_csv <- file.path(results_dir, "misty", paste0(sample_id, "_misty_improvements.csv"))
  misty_df <- cci_summary_read_csv(misty_csv)
  if (!is.null(misty_df) && nrow(misty_df) > 0) {
    metric_col <- cci_summary_pick_numeric_column(
      misty_df,
      preferred = c("improvement", "gain.R2", "gain_R2", "multi.R2", "multi_R2")
    )
    target_col <- intersect(c("target", "Target", "target.marker"), colnames(misty_df))
    if (!is.null(metric_col) && length(target_col) > 0) {
      top_misty <- misty_df %>%
        dplyr::arrange(dplyr::desc(.data[[metric_col]])) %>%
        dplyr::slice_head(n = min(top_n, nrow(misty_df)))
      
      add_summary_row("misty", "targets_modelled", nrow(misty_df))
      add_summary_row("misty", "top_target", top_misty[[target_col[1]]][1])
      
      p_misty <- ggplot2::ggplot(
        top_misty,
        ggplot2::aes(
          x = stats::reorder(.data[[target_col[1]]], .data[[metric_col]]),
          y = .data[[metric_col]]
        )
      ) +
        ggplot2::geom_col(fill = "#1B9E77") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Top MISTy targets by improvement",
          x = "Target",
          y = metric_col
        ) +
        ggplot2::theme_minimal(base_size = 11)
      
      ggplot2::ggsave(
        filename = file.path(summary_dir, paste0(sample_id, "_misty_top_targets.png")),
        plot = p_misty,
        width = 10,
        height = 6,
        dpi = 150
      )
    }
  }
  
  liana_csv <- file.path(results_dir, "liana", paste0(sample_id, "_liana_aggregate.csv"))
  liana_df <- cci_summary_read_csv(liana_csv)
  if (!is.null(liana_df) && nrow(liana_df) > 0) {
    sender_col <- intersect(c("source", "source_group", "source_groups"), colnames(liana_df))
    receiver_col <- intersect(c("target", "target_group", "target_groups"), colnames(liana_df))
    ligand_col <- intersect(c("ligand", "ligand_complex"), colnames(liana_df))
    receptor_col <- intersect(c("receptor", "receptor_complex"), colnames(liana_df))
    metric_col <- cci_summary_pick_numeric_column(
      liana_df,
      preferred = c("aggregate_rank", "magnitude_rank", "specificity_rank", "lrscore")
    )
    if (!is.null(metric_col)) {
      if (identical(metric_col, "aggregate_rank")) {
        liana_df <- liana_df %>% dplyr::arrange(.data[[metric_col]])
      } else {
        liana_df <- liana_df %>% dplyr::arrange(dplyr::desc(.data[[metric_col]]))
      }
      top_liana <- liana_df %>% dplyr::slice_head(n = min(top_n, nrow(liana_df)))
      if (length(ligand_col) > 0 && length(receptor_col) > 0) {
        top_liana$interaction_label <- paste(
          top_liana[[ligand_col[1]]],
          top_liana[[receptor_col[1]]],
          sep = " -> "
        )
      } else {
        top_liana$interaction_label <- seq_len(nrow(top_liana))
      }
      
      add_summary_row("liana", "interactions_ranked", nrow(liana_df))
      
      p_liana <- ggplot2::ggplot(
        top_liana,
        ggplot2::aes(
          x = stats::reorder(interaction_label, .data[[metric_col]]),
          y = .data[[metric_col]]
        )
      ) +
        ggplot2::geom_col(fill = "#A6611A") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Top LIANA interactions",
          x = "Interaction",
          y = metric_col
        ) +
        ggplot2::theme_minimal(base_size = 11)
      
      ggplot2::ggsave(
        filename = file.path(summary_dir, paste0(sample_id, "_liana_top_interactions.png")),
        plot = p_liana,
        width = 10,
        height = 6,
        dpi = 150
      )
    }
  }
  
  nichenet_csv <- list.files(
    file.path(results_dir, "nichenet"),
    pattern = paste0("^", sample_id, "_ligand_activities\\.csv$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(nichenet_csv) > 0) {
    nichenet_df <- cci_summary_read_csv(nichenet_csv[1])
    if (!is.null(nichenet_df) && nrow(nichenet_df) > 0 &&
        all(c("test_ligand", "pearson") %in% colnames(nichenet_df))) {
      top_nichenet <- nichenet_df %>%
        dplyr::arrange(dplyr::desc(pearson)) %>%
        dplyr::slice_head(n = min(top_n, nrow(nichenet_df)))
      
      add_summary_row("nichenet", "ligands_ranked", nrow(nichenet_df))
      add_summary_row("nichenet", "top_ligand", top_nichenet$test_ligand[1])
      
      p_nichenet <- ggplot2::ggplot(
        top_nichenet,
        ggplot2::aes(x = stats::reorder(test_ligand, pearson), y = pearson)
      ) +
        ggplot2::geom_col(fill = "#8C6BB1") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Top NicheNet ligands",
          x = "Ligand",
          y = "Pearson activity score"
        ) +
        ggplot2::theme_minimal(base_size = 11)
      
      ggplot2::ggsave(
        filename = file.path(summary_dir, paste0(sample_id, "_nichenet_top_ligands.png")),
        plot = p_nichenet,
        width = 10,
        height = 6,
        dpi = 150
      )
    }
  }
  
  summary_df <- if (length(summary_rows) > 0) {
    dplyr::bind_rows(summary_rows)
  } else {
    data.frame(
      section = character(),
      metric = character(),
      value = character(),
      stringsAsFactors = FALSE
    )
  }
  
  readr::write_csv(
    summary_df,
    file.path(summary_dir, paste0(sample_id, "_cci_summary_statistics.csv"))
  )
  
  cat("✓ CCI summary outputs saved to:", summary_dir, "\n\n")
  
  invisible(list(
    summary_dir = summary_dir,
    statistics = summary_df
  ))
}

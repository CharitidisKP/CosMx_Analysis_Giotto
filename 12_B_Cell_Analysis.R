#!/usr/bin/env Rscript
# ==============================================================================
# 12_B_Cell_Analysis.R
# B-cell focused spatial interaction summaries
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
    stop("No annotation column was found for B-cell analysis.")
  }
  
  candidates[1]
}

ensure_spatial_network <- function(gobj, spatial_network_name = "Delaunay_network") {
  network_exists <- tryCatch({
    get_spatialNetwork(
      gobject = gobj,
      name = spatial_network_name,
      output = "networkDT"
    )
    TRUE
  }, error = function(e) FALSE)
  
  if (network_exists) {
    return(gobj)
  }
  
  createSpatialNetwork(
    gobject = gobj,
    name = spatial_network_name,
    method = "Delaunay",
    minimum_k = 2
  )
}

run_bcell_microenvironment_analysis <- function(gobj,
                                                sample_id,
                                                output_dir,
                                                annotation_column = NULL,
                                                bcell_regex = "B|Plasma|Plasmablast",
                                                spatial_network_name = "Delaunay_network",
                                                number_of_simulations = 250,
                                                save_object = FALSE) {
  cat("\n========================================\n")
  cat("STEP 12: B-Cell Microenvironment\n")
  cat("Sample:", sample_id, "\n")
  cat("========================================\n\n")
  
  if (is.character(gobj)) {
    manifest_path <- file.path(gobj, "manifest.json")
    if (file.exists(manifest_path)) {
      gobj <- load_giotto_checkpoint(gobj)
    } else {
      gobj <- loadGiotto(gobj)
    }
  }
  
  results_dir <- ensure_dir(file.path(output_dir, "12_BCell_Microenvironment"))
  metadata <- pDataDT(gobj) %>% as_tibble()
  annotation_column <- detect_annotation_column(metadata, annotation_column)
  
  abundance_table <- metadata %>%
    count(.data[[annotation_column]], sort = TRUE, name = "n_cells")
  names(abundance_table)[1] <- "annotation"
  abundance_table <- abundance_table %>%
    mutate(
      sample_id = sample_id,
      is_bcell = grepl(bcell_regex, annotation, ignore.case = TRUE)
    )
  
  readr::write_csv(
    abundance_table,
    file.path(results_dir, paste0(sample_id, "_celltype_abundance.csv"))
  )
  
  gobj <- ensure_spatial_network(gobj, spatial_network_name = spatial_network_name)
  
  cp_scores <- cellProximityEnrichment(
    gobject = gobj,
    spatial_network_name = spatial_network_name,
    cluster_column = annotation_column,
    number_of_simulations = number_of_simulations
  )
  
  readr::write_csv(
    as_tibble(cp_scores$raw_sim_table),
    file.path(results_dir, paste0(sample_id, "_cell_proximity_raw.csv"))
  )
  
  enrichment_table <- as_tibble(cp_scores$enrichm_res)
  bcell_cols <- grep("cell.*type|celltype|cluster", names(enrichment_table), value = TRUE, ignore.case = TRUE)
  if (length(bcell_cols) >= 2) {
    bcell_table <- enrichment_table %>%
      filter(
        grepl(bcell_regex, .data[[bcell_cols[1]]], ignore.case = TRUE) |
          grepl(bcell_regex, .data[[bcell_cols[2]]], ignore.case = TRUE)
      )
  } else {
    bcell_table <- enrichment_table
  }
  
  readr::write_csv(
    enrichment_table,
    file.path(results_dir, paste0(sample_id, "_cell_proximity_enrichment.csv"))
  )
  readr::write_csv(
    bcell_table,
    file.path(results_dir, paste0(sample_id, "_bcell_interactions.csv"))
  )
  
  # Bring forward B-cell focused CCI outputs from step 10 when available.
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
        bcell_liana <- liana_tbl %>%
          dplyr::filter(
            grepl(bcell_regex, source, ignore.case = TRUE) |
              grepl(bcell_regex, target, ignore.case = TRUE)
          )
        readr::write_csv(
          bcell_liana,
          file.path(results_dir, paste0(sample_id, "_bcell_liana_interactions.csv"))
        )
      }
    }, error = function(e) {
      message("B-cell LIANA summary skipped: ", conditionMessage(e))
    })
  }
  
  nichenet_root <- file.path(output_dir, "10_CCI_Analysis", "nichenet")
  if (dir.exists(nichenet_root)) {
    tryCatch({
      comparison_dirs <- list.dirs(nichenet_root, full.names = TRUE, recursive = FALSE)
      comparison_dirs <- comparison_dirs[
        grepl(bcell_regex, basename(comparison_dirs), ignore.case = TRUE)
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
      message("B-cell NicheNet summary skipped: ", conditionMessage(e))
    })
  }
  
  tryCatch({
    heatmap_plot <- cellProximityHeatmap(gobject = gobj, CPscore = cp_scores)
    ggsave(
      filename = file.path(results_dir, paste0(sample_id, "_cell_proximity_heatmap.png")),
      plot = heatmap_plot,
      width = 12,
      height = 10,
      dpi = 300,
      bg = "white"
    )
  }, error = function(e) {
    message("cellProximityHeatmap() skipped: ", conditionMessage(e))
  })
  
  tryCatch({
    net_plot <- cellProximityNetwork(gobject = gobj, CPscore = cp_scores)
    ggsave(
      filename = file.path(results_dir, paste0(sample_id, "_cell_proximity_network.png")),
      plot = net_plot,
      width = 12,
      height = 10,
      dpi = 300,
      bg = "white"
    )
  }, error = function(e) {
    message("cellProximityNetwork() skipped: ", conditionMessage(e))
  })
  
  if (save_object) {
    save_giotto_checkpoint(
      gobj = gobj,
      checkpoint_dir = file.path(output_dir, "Giotto_Object_BCell_Analysis"),
      metadata = list(stage = "bcell_microenvironment", annotation_column = annotation_column)
    )
  }
  
  invisible(
    list(
      annotation_column = annotation_column,
      abundance_table = abundance_table,
      enrichment_table = enrichment_table,
      bcell_table = bcell_table
    )
  )
}

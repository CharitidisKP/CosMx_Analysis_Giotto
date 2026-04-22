## Choose columns function ##
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

## Composite sub-sample discovery ##
# Look up the per-sub-biopsy split rows for a composite sample in the sample
# sheet. Returns a data.frame (rows with sample_id, fov_min, fov_max, ...),
# NULL if the current sample is not a composite, the sheet is unavailable, or
# no split rows exist.
#
# Shared by 02_Quality_Control.R and 08_Visualisation.R; depends on
# safe_read_sheet() from Pipeline_Utils.R, which every pipeline script sources.
discover_composite_subsamples <- function(sample_row, sample_sheet_path) {
  if (is.null(sample_row) || is.null(sample_sheet_path)) return(NULL)
  if (!file.exists(sample_sheet_path)) return(NULL)
  if (!exists("safe_read_sheet", mode = "function", inherits = TRUE)) return(NULL)

  role <- tryCatch(as.character(sample_row$split_role)[1],
                   error = function(e) NA_character_)
  if (is.null(role) || is.na(role) || !nzchar(role) || role != "composite") {
    return(NULL)
  }

  slide_num <- tryCatch(as.character(sample_row$slide_num)[1],
                        error = function(e) NA_character_)

  sheet <- tryCatch(safe_read_sheet(sample_sheet_path),
                    error = function(e) NULL)
  if (is.null(sheet)) return(NULL)
  sheet <- as.data.frame(sheet, stringsAsFactors = FALSE)

  matches_slide <- if (!is.na(slide_num) && "slide_num" %in% names(sheet)) {
    as.character(sheet$slide_num) == slide_num
  } else rep(TRUE, nrow(sheet))

  keep <- matches_slide &
    !is.na(sheet$split_role) &
    as.character(sheet$split_role) == "split"

  subs <- sheet[keep, , drop = FALSE]
  if (nrow(subs) == 0) return(NULL)
  subs <- subs[!is.na(subs$fov_min) & !is.na(subs$fov_max), , drop = FALSE]
  if (nrow(subs) == 0) return(NULL)
  subs
}

## XY coordinate thingies ##
add_xy_to_qc <- function(qc, gobj,
                         spat_unit = "cell",
                         x_candidates = c("x", "sdimx", "x_slide_mm", "CenterX_global_px", "CenterX_local_px"),
                         y_candidates = c("y", "sdimy", "y_slide_mm", "CenterY_global_px", "CenterY_local_px")) {
  
  # 1) If qc already has x/y candidates, use those
  x_col <- pick_col(qc, x_candidates)
  y_col <- pick_col(qc, y_candidates)
  
  if (!is.na(x_col) && !is.na(y_col)) {
    return(qc %>%
             mutate(
               x = suppressWarnings(as.numeric(.data[[x_col]])),
               y = suppressWarnings(as.numeric(.data[[y_col]]))
             ))
  }
  
  # 2) Otherwise pull from Giotto spatLocs raw coordinates
  coords <- gobj@spatial_locs[[spat_unit]]$raw@coordinates
  coords_df <- as.data.frame(coords)
  
  # Try to pick named columns; otherwise fallback to first two
  x2 <- pick_col(coords_df, x_candidates)
  y2 <- pick_col(coords_df, y_candidates)
  
  if (is.na(x2) || is.na(y2)) {
    x2 <- names(coords_df)[1]
    y2 <- names(coords_df)[2]
  }
  
  coords_df2 <- tibble(
    cell_ID = gobj@cell_ID$cell,
    x = suppressWarnings(as.numeric(coords_df[[x2]])),
    y = suppressWarnings(as.numeric(coords_df[[y2]]))
  )
  
  qc %>% left_join(coords_df2, by = "cell_ID")
}


## Area thingies ##
add_area_to_qc <- function(qc, area_candidates = c("Area.um2", "Area", "CellArea", "area")) {
  area_col <- pick_col(qc, area_candidates)
  
  if (is.na(area_col)) {
    qc %>% mutate(area_um2 = NA_real_)
  } else {
    qc %>% mutate(area_um2 = suppressWarnings(as.numeric(.data[[area_col]])))
  }
}

## Background thingies ##
add_bg_to_qc <- function(qc,
                         neg_candidates = c("propNegative", "prop_neg", "NegProbFrac"),
                         neg_count_candidates = c("nCount_negprobes", "negprobe_counts"),
                         false_count_candidates = c("nCount_falsecode", "falsecode_counts"),
                         total_count_col = "Total_Counts") {
  
  neg_col   <- pick_col(qc, neg_candidates)
  negc_col  <- pick_col(qc, neg_count_candidates)
  fcc_col   <- pick_col(qc, false_count_candidates)
  
  qc <- qc %>%
    mutate(
      bg_neg = if (!is.na(neg_col)) suppressWarnings(as.numeric(.data[[neg_col]])) else NA_real_,
      negprobe_cnt = if (!is.na(negc_col)) suppressWarnings(as.numeric(.data[[negc_col]])) else NA_real_,
      falsecode_cnt = if (!is.na(fcc_col)) suppressWarnings(as.numeric(.data[[fcc_col]])) else NA_real_,
      bg_neg = coalesce(bg_neg, negprobe_cnt / .data[[total_count_col]]),
      bg_falsecode = falsecode_cnt / .data[[total_count_col]]
    )
  
  qc
}

# Load reference profiles from Nanostring ---------------------------------
read_profile_csv <- function(url) {
  as.matrix(read.csv(url, row.names = 1, check.names = FALSE))
}

cleanup_memory <- function(remove = NULL,
                           envir = parent.frame(),
                           verbose = TRUE,
                           label = NULL) {
  if (!is.null(remove)) {
    remove <- unique(as.character(remove))
    remove <- remove[nzchar(remove)]
    remove <- remove[vapply(remove, exists, logical(1), envir = envir, inherits = FALSE)]
    if (length(remove) > 0) {
      rm(list = remove, envir = envir)
    }
  }
  
  gc_out <- gc(verbose = FALSE)
  if (isTRUE(verbose)) {
    used_mb <- round(sum(gc_out[, "used"]) / 1024, 1)
    prefix <- if (!is.null(label) && nzchar(label)) paste0(label, ": ") else ""
    cat(prefix, "gc() complete; approx. ", used_mb, " MB tracked in use\n", sep = "")
  }
  invisible(gc_out)
}

sanitize_smide_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", as.character(x))
}

parse_smide_unified_interactions <- function(unified_int) {
  parts <- strsplit(as.character(unified_int), "--", fixed = TRUE)
  data.frame(
    source = vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)),
    target = vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

resolve_smide_gobj <- function(gobj_ref) {
  if (is.null(gobj_ref) || !is.character(gobj_ref)) {
    return(gobj_ref)
  }
  
  manifest_path <- file.path(gobj_ref, "manifest.json")
  if (file.exists(manifest_path) && exists("load_giotto_checkpoint", mode = "function")) {
    return(load_giotto_checkpoint(gobj_ref))
  }
  
  if (exists(".giotto_load", mode = "function")) {
    return(.giotto_load(gobj_ref))
  }
  
  gobj_ref
}

get_smide_neighbor_counts <- function(metadata, annotation_column) {
  neighbor_counts <- getOption("cosmx.smide.neighbor_counts", NULL)
  if (!is.null(neighbor_counts)) {
    return(neighbor_counts)
  }
  
  gobj_ref <- resolve_smide_gobj(getOption("cosmx.smide.gobj", NULL))
  spatial_network_name <- getOption("cosmx.smide.spatial_network_name", NULL)
  if (is.null(gobj_ref) || is.null(spatial_network_name)) {
    return(NULL)
  }
  
  build_neighbourhood_matrix(
    gobj = gobj_ref,
    metadata = metadata,
    annotation_column = annotation_column,
    spatial_network_name = spatial_network_name
  )
}

select_smide_partner_celltypes <- function(output_dir,
                                           run_label,
                                           focal_celltype,
                                           top_n = 3,
                                           padj_threshold = 0.05,
                                           include_self = FALSE) {
  cp_path <- file.path(
    output_dir,
    "11_BCell_Microenvironment",
    paste0(run_label, "_cell_proximity_enrichment.csv")
  )
  if (!file.exists(cp_path)) {
    return(character())
  }
  
  tbl <- tryCatch(
    readr::read_csv(cp_path, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (is.null(tbl) || !"unified_int" %in% names(tbl)) {
    return(character())
  }
  
  parsed <- parse_smide_unified_interactions(tbl$unified_int)
  tbl$source <- parsed$source
  tbl$target <- parsed$target
  
  partner_tbl <- tbl[tbl$source == focal_celltype | tbl$target == focal_celltype, , drop = FALSE]
  if (!include_self) {
    partner_tbl <- partner_tbl[!(partner_tbl$source == focal_celltype & partner_tbl$target == focal_celltype), , drop = FALSE]
  }
  if ("p.adj_higher" %in% names(partner_tbl)) {
    partner_tbl <- partner_tbl[is.na(partner_tbl$p.adj_higher) | partner_tbl$p.adj_higher <= padj_threshold, , drop = FALSE]
  }
  if (nrow(partner_tbl) == 0) {
    return(character())
  }
  
  partner_tbl$partner_celltype <- ifelse(partner_tbl$source == focal_celltype, partner_tbl$target, partner_tbl$source)
  partner_tbl <- partner_tbl[!is.na(partner_tbl$partner_celltype) & nzchar(partner_tbl$partner_celltype), , drop = FALSE]
  if (nrow(partner_tbl) == 0) {
    return(character())
  }
  
  if ("int_ranking" %in% names(partner_tbl)) {
    partner_tbl <- partner_tbl[order(partner_tbl$int_ranking, partner_tbl$partner_celltype), , drop = FALSE]
  }
  
  unique(head(partner_tbl$partner_celltype, top_n))
}

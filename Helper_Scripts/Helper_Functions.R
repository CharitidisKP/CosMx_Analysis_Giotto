## Choose columns function ## 
pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
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

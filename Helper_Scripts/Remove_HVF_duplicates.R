# Function to remove HVF duplicate cells ----------------------------------

remove_hvf_duplicates <- function(gobj, hvgs, output_dir = NULL, verbose = TRUE) {
  
  ## Get normalized expression for HVFs ##
  expr_norm <- getExpression(gobj, values = "normalized")@exprMat
  expr_hvf <- expr_norm[hvgs, , drop = FALSE]
  
  ## Check for duplicates ##
  expr_hvf_dense <- as.matrix(expr_hvf)
  dup_idx <- duplicated(t(expr_hvf_dense))
  
  if (sum(dup_idx) > 0) {
    if (verbose) {
      cat("Found", sum(dup_idx), "cells with duplicate HVF profiles\n")
    }
    
    all_cells <- colnames(expr_norm)
    dup_cells <- all_cells[dup_idx]
    
    ## Diagnostic report ##
    if (!is.null(output_dir)) {
      cell_sigs <- apply(expr_hvf_dense, 2, function(x) paste(round(x, 6), collapse = "||"))
      
      dup_report <- data.frame(
        cell_ID = all_cells,
        signature = cell_sigs,
        is_duplicate = dup_idx,
        total_hvf_counts = colSums(expr_hvf_dense),
        n_hvf_expressed = colSums(expr_hvf_dense > 0),
        stringsAsFactors = FALSE
      )
      
      dup_groups <- dup_report %>%
        group_by(signature) %>%
        filter(n() > 1) %>%
        mutate(group_id = cur_group_id()) %>%
        arrange(group_id, cell_ID)
      
      if (nrow(dup_groups) > 0) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        write_csv(dup_groups, file.path(output_dir, "hvf_duplicate_cells.csv"))
        
        if (verbose) {
          cat("\nDuplicate summary:\n")
          print(dup_groups %>% 
                  group_by(group_id) %>% 
                  summarise(n_cells = n(), 
                            total_hvf = first(total_hvf_counts),
                            n_hvf_expr = first(n_hvf_expressed)))
        }
      }
    }
    
    ## Cells to keep ##
    cells_to_keep <- all_cells[!dup_idx]
    
    if (verbose) cat("\nSubsetting Giotto object...\n")
    
    ## Use subsetGiotto with poly_info = NULL ##
    gobj_sub <- tryCatch({
      subsetGiotto(
        gobject = gobj,
        cell_ids = cells_to_keep,
        poly_info = NULL
      )
    }, error = function(e) {
      if (verbose) cat("subsetGiotto failed, using manual method...\n")
      NULL
    })
    
    ## If subsetGiotto failed, manually subset ##
    if (is.null(gobj_sub)) {
      if (verbose) cat("Manually subsetting object components...\n")
      
      ## Get expression objects ##
      expr_raw_obj <- getExpression(gobj, values = "raw", output = "exprObj")
      expr_norm_obj <- getExpression(gobj, values = "normalized", output = "exprObj")
      
      ## Subset matrices ##
      expr_raw_obj@exprMat <- expr_raw_obj@exprMat[, cells_to_keep, drop = FALSE]
      expr_norm_obj@exprMat <- expr_norm_obj@exprMat[, cells_to_keep, drop = FALSE]
      
      ## Update expression ##
      gobj <- setExpression(gobject = gobj, x = expr_raw_obj, name = "raw")
      gobj <- setExpression(gobject = gobj, x = expr_norm_obj, name = "normalized")
      
      ## Update cell metadata ##
      meta <- pDataDT(gobj)
      meta_sub <- meta[meta$cell_ID %in% cells_to_keep, ]
      gobj@cell_metadata$cell <- meta_sub
      gobj@cell_ID$cell <- cells_to_keep
      
      ## Update spatial locations ##
      if (length(gobj@spatial_locs) > 0) {
        if (verbose) cat("Updating spatial locations...\n")
        spatlocs <- getSpatialLocations(gobj, spat_unit = "cell", output = "data.table")
        spatlocs_sub <- spatlocs[spatlocs$cell_ID %in% cells_to_keep, ]
        gobj <- setSpatialLocations(gobject = gobj, spatlocs = spatlocs_sub, spat_unit = "cell")
      }
      
      gobj_sub <- gobj
    }
    
    ## Handle polygons - simplified approach ##
    if (!is.null(gobj_sub@spatial_info$cell)) {
      if (verbose) cat("Reconstructing polygon data...\n")
      
      poly_obj <- gobj_sub@spatial_info$cell
      all_poly_ids <- poly_obj@unique_ID_cache
      keep_poly_idx <- all_poly_ids %in% cells_to_keep
      
      n_poly_remove <- sum(!keep_poly_idx)
      
      if (n_poly_remove > 0) {
        if (verbose) cat("  Subsetting", n_poly_remove, "polygons...\n")
        
        ## Subset spatVector ##
        sv_new <- poly_obj@spatVector[keep_poly_idx, ]
        
        ## Create a NEW giottoPolygon object by copying the old one ##
        poly_new <- poly_obj  # Start with a copy
        
        ## Update the components ##
        poly_new@spatVector <- sv_new
        poly_new@unique_ID_cache <- all_poly_ids[keep_poly_idx]
        poly_new@spatVectorCentroids <- terra::centroids(sv_new)
        
        ## Update overlap info if it exists ##
        if (length(poly_new@overlaps) > 0) {
          poly_new@overlaps <- list()
        }
        
        ## Assign to object ##
        gobj_sub@spatial_info$cell <- poly_new
        
        if (verbose) cat("Polygons updated\n")
      } else {
        if (verbose) cat("No polygons need removal\n")
      }
    }
    
    if (verbose) {
      cat("\nRemoved", sum(dup_idx), "duplicate cells\n")
      cat("Retained", length(cells_to_keep), "cells\n")
    }
    
    return(list(
      gobject = gobj_sub,
      n_removed = sum(dup_idx),
      removed_cells = dup_cells
    ))
    
  } else {
    if (verbose) cat("No HVF duplicates found\n")
    return(list(
      gobject = gobj,
      n_removed = 0,
      removed_cells = character(0)
    ))
  }
}
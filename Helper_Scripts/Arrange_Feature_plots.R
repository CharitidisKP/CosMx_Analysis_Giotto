## Function to find the best grid divisions for a given number of markers ##
find_optimal_grid_layout <- function(n_markers, allowed_grids = c(4, 6, 9, 16, 25, 36)) {
  
  if (n_markers == 0) return(list(layout = integer(0), total_panels = 0))
  
  ## If markers fit in a single allowed grid, use it ##
  if (n_markers %in% allowed_grids) {
    return(list(
      layout = n_markers,
      total_panels = 1,
      total_markers = n_markers,
      panel_sizes = n_markers
    ))
  }
  
  ## Calculate how many complete panels of each size we can make ##
  best_layout <- NULL
  best_score <- Inf  ## Lower is better
  
  ## Try all combinations, starting with largest grids ##
  for (main_grid in sort(allowed_grids, decreasing = TRUE)) {
    
    n_complete <- floor(n_markers / main_grid)
    remainder <- n_markers %% main_grid
    
    ## Build layout ##
    layout <- rep(main_grid, n_complete)
    
    ## Handle remainder ##
    if (remainder > 0) {
      ## Try to use the largest allowed grid that fits the remainder ##
      remainder_grid <- allowed_grids[allowed_grids >= remainder]
      
      if (length(remainder_grid) > 0) {
        layout <- c(layout, remainder_grid[1])
      } else {
        ## If remainder is smaller than smallest allowed grid ##
        if (remainder <= 3 && n_complete > 0) {
          ## Try to merge with last complete panel if reasonable ##
          new_size <- main_grid + remainder
          if (new_size <= max(allowed_grids) * 1.5) {
            layout[length(layout)] <- new_size
          } else {
            layout <- c(layout, remainder)
          }
        } else {
          layout <- c(layout, remainder)
        }
      }
    }
    
    ## Calculate score (prefer fewer panels with more balanced sizes) ##
    n_panels <- length(layout)
    
    ## Safe variance calculation ##
    if (length(layout) > 1) {
      size_variance <- var(layout)
    } else {
      size_variance <- 0
    }
    
    ## Count problematic single-marker panels ##
    lonely_singles <- sum(layout == 1)
    small_panels <- sum(layout < 4 & layout > 1)
    
    ## Penalize: single markers heavily, small panels moderately ##
    score <- n_panels + (size_variance * 0.1) + (lonely_singles * 100) + (small_panels * 10)
    
    if (score < best_score) {
      best_score <- score
      best_layout <- layout
    }
  }
  
  ## Final check: if we ended up with tiny panels, try to consolidate ##
  if (length(best_layout) > 1) {
    small_idx <- which(best_layout < 4)
    if (length(small_idx) > 1) {
      ## Try to merge small panels ##
      small_sum <- sum(best_layout[small_idx])
      if (small_sum %in% allowed_grids) {
        best_layout <- c(best_layout[-small_idx], small_sum)
      }
    }
  }
  
  list(
    layout = best_layout,
    total_panels = length(best_layout),
    total_markers = sum(best_layout),
    panel_sizes = best_layout
  )
}

## Function to determine grid dimensions for a given number of features ##
optimal_grid_dims <- function(n_features) {
  if (n_features <= 0) return(list(ncol = 1, nrow = 1))
  
  ## Prefer layouts close to square ##
  if (n_features <= 4) {
    return(list(ncol = 2, nrow = 2))
  } else if (n_features <= 6) {
    return(list(ncol = 3, nrow = 2))
  } else if (n_features <= 9) {
    return(list(ncol = 3, nrow = 3))
  } else if (n_features <= 12) {
    return(list(ncol = 4, nrow = 3))
  } else if (n_features <= 16) {
    return(list(ncol = 4, nrow = 4))
  } else if (n_features <= 20) {
    return(list(ncol = 5, nrow = 4))
  } else {
    ## For larger numbers, try to stay squarish ##
    ncol <- ceiling(sqrt(n_features))
    nrow <- ceiling(n_features / ncol)
    return(list(ncol = ncol, nrow = nrow))
  }
}

## Function to split markers into optimally-sized chunks ##
create_marker_panels <- function(markers, allowed_grids = c(4, 6, 9, 16, 25, 36)) {
  
  n_markers <- length(markers)
  
  if (n_markers == 0) {
    warning("No markers provided")
    return(list())
  }
  
  layout_info <- find_optimal_grid_layout(n_markers, allowed_grids)
  
  cat("\n=== Marker Panel Layout ===\n")
  cat("Total markers:", n_markers, "\n")
  cat("Number of panels:", layout_info$total_panels, "\n")
  cat("Panel sizes:", paste(layout_info$panel_sizes, collapse = ", "), "\n\n")
  
  ## Split markers according to layout ##
  panels <- list()
  marker_idx <- 1
  
  for (i in seq_along(layout_info$panel_sizes)) {
    panel_size <- layout_info$panel_sizes[i]
    end_idx <- min(marker_idx + panel_size - 1, n_markers)
    
    panels[[i]] <- list(
      markers = markers[marker_idx:end_idx],
      size = panel_size,
      grid = optimal_grid_dims(panel_size)
    )
    
    marker_idx <- end_idx + 1
  }
  
  panels
}

# ## Test with various marker counts ##
# cat("\n=== Testing Layout Algorithm ===\n\n")
# for (n in 1:30) {
#   layout <- find_optimal_grid_layout(n_markers = n)
#   cat(sprintf("n=%2d: %s (panels: %d)\n", 
#               n, 
#               paste(layout$panel_sizes, collapse = " + "),
#               layout$total_panels))
# }
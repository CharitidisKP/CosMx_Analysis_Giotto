# Enhanced B Cell Marker Visualization Function
create_feature_plot <- function(feature, 
                                umap_df, 
                                expr_matrix,
                                expr_limits = NULL,
                                point_size_low = 0.1,
                                point_size_high = 0.5,
                                alpha_low = 0.3,
                                alpha_high = 0.9) {
  
  ## Get expression values ##
  if (!feature %in% rownames(expr_matrix)) {
    warning(paste("Feature", feature, "not found in expression matrix"))
    return(NULL)
  }
  
  expr_values <- expr_matrix[feature, ]
  
  ## Prepare plotting data ##
  plot_df <- umap_df
  plot_df$expression <- expr_values[plot_df$cell_ID]
  
  ## If no limits provided, calculate from data ##
  if (is.null(expr_limits)) {
    expr_limits <- c(0, max(plot_df$expression, na.rm = TRUE))
  }
  
  ## Separate zero and non-zero expression ##
  plot_df$is_expressed <- plot_df$expression > 0
  
  ## Order cells: non-expressing first, then by expression level ##
  plot_df <- plot_df %>%
    arrange(is_expressed, expression)
  
  ## Create plot ##
  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      data = plot_df %>% filter(!is_expressed),
      color = "gray20",
      size = point_size_low,
      alpha = alpha_low
    ) +
    geom_point(
      data = plot_df %>% filter(is_expressed),
      aes(color = expression),
      size = point_size_high,
      alpha = alpha_high,
      stroke = 0
    ) +
    scale_color_gradientn(
      colors = c("#440154", "springgreen"),
      name = "Expression",
      limits = expr_limits,
      na.value = "gray20"
    ) +
    labs(title = feature, x = "UMAP 1", y = "UMAP 2") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      axis.text = element_text(size = 7),
      axis.title = element_text(size = 9),
      legend.position = "right",
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.25, "cm"),
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 8),
      panel.background = element_rect(fill = "black", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    )
  
  return(p)
}

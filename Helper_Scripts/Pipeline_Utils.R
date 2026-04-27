`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    return(y)
  }
  if (length(x) == 1 && is.atomic(x) && is.na(x)) {
    return(y)
  }
  x
}

resolve_path <- function(path, base_dir = getwd(), mustWork = FALSE) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  expanded <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", expanded)) {
    expanded <- file.path(base_dir, expanded)
  }
  normalizePath(expanded, winslash = "/", mustWork = mustWork)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

parse_cli_csv <- function(value) {
  if (is.null(value) || !nzchar(value)) {
    return(NULL)
  }
  trimws(unlist(strsplit(value, ",", fixed = TRUE)))
}

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

pretty_plot_label <- function(x, width = NULL) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("_", " ", x, fixed = TRUE)
  x <- gsub("\\.", " ", x)
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  if (is.null(width)) {
    return(x)
  }
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

sample_plot_title <- function(sample_id, title) {
  if (is.null(sample_id) || !nzchar(sample_id)) {
    return(title)
  }
  paste(sample_id, title, sep = " - ")
}

# Build a human-friendly plot title from a sample_sheet row. Examples:
#   CART_pt1   T0      -> "CART S1 - Before treatment - <subtitle>"
#   CART_pt2   T12     -> "CART S2 - After treatment - <subtitle>"
#   Conv_pt1   T0      -> "Conventional S1 - Before treatment - <subtitle>"
#   Control_pt1 (any)  -> "Healthy Control - <subtitle>"
#
# sample_row is a one-row data.frame / named list from sample_sheet.csv with
# columns group_id, patient_id, timepoint. Unrecognised groups fall back to
# sample_plot_title().
plot_title_for_sample <- function(sample_row, subtitle,
                                   sample_id_fallback = NULL) {
  pick <- function(field) {
    if (is.null(sample_row)) return(NA_character_)
    val <- tryCatch(as.character(sample_row[[field]])[1],
                    error = function(e) NA_character_)
    if (length(val) == 0 || is.na(val) || !nzchar(val)) NA_character_ else val
  }

  grp <- pick("group_id")
  tp  <- pick("timepoint")
  pid <- pick("patient_id")

  if (is.na(grp)) {
    # No sample-sheet context: fall back to "sample_id - subtitle"
    return(sample_plot_title(sample_id_fallback, subtitle))
  }

  tp_label <- switch(as.character(tp),
                     "T0"  = "Before treatment",
                     "T12" = "After treatment",
                     "")

  grp_label <- switch(as.character(grp),
                      "CART"         = "CART",
                      "Conventional" = "Conventional",
                      "Conv"         = "Conventional",
                      "Control"      = "Healthy Control",
                      as.character(grp))

  # S-tag: CART_pt1 -> S1 ; Conv_pt2 -> S2 ; else blank.
  s_tag <- NA_character_
  if (!is.na(pid)) {
    m <- regmatches(pid, regexec("_pt(\\d+)$", pid))[[1]]
    if (length(m) == 2L) s_tag <- paste0("S", m[2])
  }

  prefix <- if (grp_label == "Healthy Control") {
    "Healthy Control"
  } else if (!is.na(s_tag) && nzchar(tp_label)) {
    paste(grp_label, s_tag, "-", tp_label)
  } else if (!is.na(s_tag)) {
    paste(grp_label, s_tag)
  } else {
    grp_label
  }

  paste(prefix, subtitle, sep = " - ")
}

has_ggtext <- function() {
  requireNamespace("ggtext", quietly = TRUE)
}

element_markdown_safe <- function(...) {
  if (has_ggtext()) {
    return(ggtext::element_markdown(...))
  }
  ggplot2::element_text(...)
}

display_reduction_name <- function(x) {
  value <- as.character(x)
  lower <- tolower(value)
  mapped <- dplyr::case_when(
    lower %in% c("umap") ~ "UMAP",
    lower %in% c("pca") ~ "PCA",
    lower %in% c("tsne", "t-sne") ~ "t-SNE",
    TRUE ~ value
  )
  mapped
}

embedding_axis_label <- function(reduction, dimension) {
  paste(display_reduction_name(reduction), "Dimension", dimension)
}

log10_axis_label <- function(prefix, suffix = "Scale") {
  if (has_ggtext()) {
    return(paste0(prefix, " (log<sub>10</sub> ", suffix, ")"))
  }
  bquote(.(prefix) ~ "(log"[10] ~ .(suffix) * ")")
}

presentation_theme <- function(base_size = 13,
                               legend_position = "right",
                               x_angle = 0,
                               y_angle = 0) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = element_markdown_safe(
        hjust = 0.5,
        face = "bold",
        size = base_size + 2,
        margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = element_markdown_safe(
        hjust = 0.5,
        size = base_size,
        color = "grey20",
        margin = ggplot2::margin(b = 10)
      ),
      axis.title.x = element_markdown_safe(
        face = "bold",
        size = base_size + 1,
        margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
      ),
      axis.title.y = element_markdown_safe(
        face = "bold",
        size = base_size + 1,
        margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
      ),
      axis.text.x = ggplot2::element_text(
        size = base_size - 1,
        color = "grey20",
        angle = x_angle,
        hjust = if (x_angle == 0) 0.5 else 1,
        vjust = if (x_angle == 0) 0.5 else 1
      ),
      axis.text.y = ggplot2::element_text(
        size = base_size - 1,
        color = "grey20",
        angle = y_angle
      ),
      legend.position = legend_position,
      legend.title = element_markdown_safe(
        face = "bold",
        size = base_size,
        margin = ggplot2::margin(b = 4)
      ),
      legend.text = ggplot2::element_text(
        size = base_size - 2,
        lineheight = 0.95
      ),
      legend.spacing.y = grid::unit(0.12, "cm"),
      panel.border = ggplot2::element_rect(
        colour = "grey80",
        fill = NA,
        linewidth = 0.4
      ),
      plot.margin = ggplot2::margin(12, 16, 12, 12)
    )
}

save_presentation_plot <- function(plot,
                                   filename,
                                   width = 10,
                                   height = 8,
                                   dpi = 300,
                                   bg = "white") {
  ext <- tolower(tools::file_ext(filename))

  # PNG output goes through base-R grDevices::png() first - avoids the
  # xfun-version minefield that ggplot2::ggsave() can trip when certain
  # rendering deps are upgraded mid-session. On failure we fall through
  # to ggsave.
  if (ext == "png") {
    ok <- tryCatch({
      grDevices::png(
        filename = filename,
        width = width, height = height, units = "in",
        res = dpi, bg = bg, type = "cairo"
      )
      on.exit(grDevices::dev.off(), add = TRUE)
      print(plot)
      TRUE
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      message("  [save] grDevices::png() failed (", conditionMessage(e),
              "), falling back to ggsave()")
      FALSE
    })
    if (ok) return(invisible(filename))
  }

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = bg,
    limitsize = FALSE
  )
}

write_json_pretty <- function(x, path, auto_unbox = TRUE) {
  jsonlite::write_json(
    x = x,
    path = path,
    pretty = TRUE,
    auto_unbox = auto_unbox,
    null = "null"
  )
}

safe_read_sheet <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    return(readxl::read_excel(path))
  }
  readr::read_csv(path, show_col_types = FALSE)
}

save_giotto_checkpoint <- function(gobj,
                                   checkpoint_dir,
                                   overwrite = TRUE,
                                   metadata = list()) {
  checkpoint_dir <- ensure_dir(checkpoint_dir)
  save_method <- "none"
  error_message <- NULL
  
  giotto_dir <- file.path(checkpoint_dir, "giotto")
  qs_file <- file.path(checkpoint_dir, "object.qs")
  rds_file <- file.path(checkpoint_dir, "object.rds")
  
  unlink(giotto_dir, recursive = TRUE, force = TRUE)
  if (overwrite) {
    unlink(qs_file, force = TRUE)
    unlink(rds_file, force = TRUE)
  }
  
  tryCatch({
    saveGiotto(
      gobject = gobj,
      dir = checkpoint_dir,
      foldername = "giotto",
      overwrite = overwrite
    )
    save_method <- "giotto"
  }, error = function(e) {
    error_message <<- conditionMessage(e)
  })
  
  if (save_method == "none" && requireNamespace("qs2", quietly = TRUE)) {
    tryCatch({
      qs2::qs_save(gobj, qs_file, compress_level = 1)
      save_method <- "qs2"
      error_message <- NULL
    }, error = function(e) {
      error_message <<- conditionMessage(e)
    })
  }
  
  if (save_method == "none") {
    tryCatch({
      saveRDS(gobj, rds_file, compress = "xz")
      save_method <- "rds"
      error_message <- NULL
    }, error = function(e) {
      error_message <<- conditionMessage(e)
    })
  }
  
  if (save_method == "none") {
    stop(
      "save_giotto_checkpoint: all three serialisation backends failed for '",
      checkpoint_dir, "'. Last error: ",
      if (is.null(error_message)) "<none>" else error_message,
      "\n  Tried: saveGiotto(), qs2::qs_save(), saveRDS().",
      "\n  Check disk space, write permissions, and the Giotto object integrity."
    )
  }

  write_json_pretty(
    c(
      list(
        saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        save_method = save_method,
        checkpoint_dir = checkpoint_dir
      ),
      metadata
    ),
    file.path(checkpoint_dir, "manifest.json")
  )

  invisible(
    list(
      checkpoint_dir = checkpoint_dir,
      save_method = save_method,
      error_message = error_message
    )
  )
}

load_giotto_checkpoint <- function(checkpoint_dir) {
  checkpoint_dir <- resolve_path(checkpoint_dir, mustWork = TRUE)

  giotto_dir <- file.path(checkpoint_dir, "giotto")
  qs_file <- file.path(checkpoint_dir, "object.qs")
  rds_file <- file.path(checkpoint_dir, "object.rds")
  manifest_file <- file.path(checkpoint_dir, "manifest.json")

  if (dir.exists(giotto_dir)) {
    return(loadGiotto(giotto_dir))
  }
  if (file.exists(qs_file) && requireNamespace("qs2", quietly = TRUE)) {
    return(qs2::qs_read(qs_file))
  }
  if (file.exists(rds_file)) {
    return(readRDS(rds_file))
  }

  # Manifest without payload signals a prior aborted save - be explicit.
  if (file.exists(manifest_file)) {
    stop(
      "Checkpoint manifest present but no payload (giotto/, object.qs, object.rds) in ",
      checkpoint_dir, ". Rerun the producing step."
    )
  }
  stop("No checkpoint payload found in ", checkpoint_dir)
}

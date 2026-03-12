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
  
  if (save_method == "none" && requireNamespace("qs", quietly = TRUE)) {
    tryCatch({
      qs::qsave(gobj, qs_file, preset = "fast")
      save_method <- "qs"
      error_message <- NULL
    }, error = function(e) {
      error_message <<- conditionMessage(e)
    })
  }
  
  if (save_method == "none") {
    saveRDS(gobj, rds_file, compress = "xz")
    save_method <- "rds"
    error_message <- NULL
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
  
  if (dir.exists(giotto_dir)) {
    return(loadGiotto(giotto_dir))
  }
  if (file.exists(qs_file) && requireNamespace("qs", quietly = TRUE)) {
    return(qs::qread(qs_file))
  }
  if (file.exists(rds_file)) {
    return(readRDS(rds_file))
  }
  
  stop("No checkpoint payload found in ", checkpoint_dir)
}

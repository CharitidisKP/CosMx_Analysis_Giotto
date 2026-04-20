#!/usr/bin/env Rscript

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

repo_dir <- current_script_dir()

default_config_path <- function(base_dir = repo_dir) {
  parameters_config <- file.path(base_dir, "Parameters", "config.yaml")
  legacy_config <- file.path(base_dir, "config.yaml")
  if (file.exists(parameters_config) || !file.exists(legacy_config)) {
    return(parameters_config)
  }
  legacy_config
}

pipeline_utils <- file.path(repo_dir, "Helper_Scripts", "Pipeline_Utils.R")
if (file.exists(pipeline_utils)) {
  source(pipeline_utils, local = TRUE)
}

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

if (!exists("ensure_dir")) {
  ensure_dir <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    normalizePath(path, winslash = "/", mustWork = FALSE)
  }
}

maybe_cleanup_runtime <- function(cfg, label = NULL) {
  if (!isTRUE(cfg$pipeline$memory_cleanup %||% TRUE)) {
    return(invisible(NULL))
  }
  
  cleanup_fn <- get0("cleanup_memory", mode = "function", inherits = TRUE)
  if (is.function(cleanup_fn)) {
    cleanup_fn(label = label %||% "Pipeline", verbose = TRUE)
  } else {
    gc(verbose = FALSE)
  }
}

if (!exists("resolve_path")) {
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
}

if (!exists("parse_cli_csv")) {
  parse_cli_csv <- function(value) {
    if (is.null(value) || !nzchar(value)) {
      return(NULL)
    }
    trimws(unlist(strsplit(value, ",", fixed = TRUE)))
  }
}

canonical_step_ids <- function(step_ids, type = c("sample", "merged")) {
  type <- match.arg(type)
  if (is.null(step_ids) || length(step_ids) == 0) {
    return(step_ids)
  }
  
  normalize_id <- function(x) {
    x <- trimws(tolower(x))
    gsub("[^a-z0-9]+", "_", x)
  }
  
  sample_aliases <- c(
    "01_load_data" = "01_load",
    "02_quality_control" = "02_qc",
    "03_normalize" = "03_norm",
    "03_normalise" = "03_norm",
    "03_normalisation" = "03_norm",
    "03_normalization" = "03_norm",
    "04_reduce" = "04_dimred",
    "04_reduction" = "04_dimred",
    "04_dimensionality_reduction" = "04_dimred",
    "04_dimensionallity_reduction" = "04_dimred",
    "06_de" = "06_markers",
    "06_differential_expression" = "06_markers",
    "06_marker_analysis" = "06_markers",
    "07_annotation" = "07_annotate",
    "08_visualization" = "08_visualize",
    "08_visualisation" = "08_visualize",
    "09_spatial_network" = "09_spatial",
    "10_cci_analysis" = "10_cci",
    "11_b_cell" = "11_bcell",
    "11_b_cell_analysis" = "11_bcell",
    "12_b_cell" = "11_bcell",
    "12_b_cell_analysis" = "11_bcell",
    "12_spatial_de" = "12_spatial_de",
    "12_spatial_differential_expression" = "12_spatial_de",
    "13_spatial_de" = "12_spatial_de",
    "13_spatial_differential_expression" = "12_spatial_de"
  )
  
  merged_aliases <- c(
    "11_batch" = "merge_batch",
    "11_batch_correction" = "merge_batch",
    "merge_batch" = "merge_batch",
    "12_spatial_de" = "12_spatial_de",
    "12_spatial_differential_expression" = "12_spatial_de",
    "13_spatial_de" = "12_spatial_de",
    "13_spatial_differential_expression" = "12_spatial_de"
  )
  
  aliases <- if (type == "sample") sample_aliases else merged_aliases
  normalized <- vapply(step_ids, normalize_id, character(1))
  unname(ifelse(normalized %in% names(aliases), aliases[normalized], normalized))
}

SAMPLE_STEP_ORDER <- c(
  "01_load",
  "02_qc",
  "03_norm",
  "04_dimred",
  "05_cluster",
  "06_markers",
  "07_annotate",
  "08_visualize",
  "09_spatial",
  "10_cci",
  "11_bcell",
  "12_spatial_de"
)

MERGED_STEP_ORDER <- c("merge", "merge_batch", "12_spatial_de")

step_label <- function(step_id) {
  labels <- c(
    "01_load" = "01 Load data",
    "02_qc" = "02 Quality control",
    "03_norm" = "03 Normalization",
    "04_dimred" = "04 Dimensionality reduction",
    "05_cluster" = "05 Clustering",
    "06_markers" = "06 Marker analysis",
    "07_annotate" = "07 Annotation",
    "08_visualize" = "08 Visualisation",
    "09_spatial" = "09 Spatial network analysis",
    "10_cci" = "10 CCI analysis",
    "merge_batch" = "Merge batch correction",
    "11_bcell" = "11 Focused cell-type microenvironment",
    "12_spatial_de" = "12 Spatial differential expression",
    "merge" = "Merge sample objects"
  )
  labels[[step_id]] %||% step_id
}

standardize_colnames <- function(x) {
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  tolower(x)
}

coerce_bool <- function(x, default = TRUE) {
  if (is.null(x) || length(x) == 0) {
    return(default)
  }
  if (is.logical(x)) {
    out <- x[1]
    if (is.na(out)) default else out
  } else if (is.numeric(x)) {
    !is.na(x[1]) && x[1] != 0
  } else {
    value <- tolower(trimws(as.character(x[1])))
    if (!nzchar(value)) {
      default
    } else {
      value %in% c("true", "t", "1", "yes", "y")
    }
  }
}

read_tabular_file <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("readxl is required to read sample sheets in Excel format.")
    }
    return(as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE))
  }
  if (ext %in% c("tsv", "txt")) {
    return(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE))
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

infer_sample_id <- function(directory_name) {
  sample_id <- sub("^.*?Slide\\d+_?", "", directory_name)
  sample_id <- sub("^_+", "", sample_id)
  if (!nzchar(sample_id) || identical(sample_id, directory_name)) {
    sample_id <- directory_name
  }
  sample_id
}

load_config <- function(config_path = default_config_path(repo_dir)) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("yaml is required to read the pipeline configuration.")
  }
  
  config_path <- resolve_path(config_path, base_dir = getwd(), mustWork = TRUE)
  cfg <- yaml::read_yaml(config_path)
  cfg$config_path <- config_path
  cfg$config_dir <- dirname(config_path)
  
  path_defaults <- list(
    project_dir = ".",
    scripts_dir = ".",
    raw_data_dir = "Data/Raw_data",
    output_dir = "Output",
    sample_sheet = "sample_sheet.csv",
    python_path = Sys.getenv("COSMX_PYTHON_PATH", unset = "")
  )
  cfg$paths <- modifyList(path_defaults, cfg$paths %||% list())
  cfg$paths$project_dir <- resolve_path(cfg$paths$project_dir, base_dir = cfg$config_dir, mustWork = FALSE)
  cfg$paths$scripts_dir <- resolve_path(cfg$paths$scripts_dir, base_dir = cfg$config_dir, mustWork = FALSE)
  cfg$paths$raw_data_dir <- resolve_path(cfg$paths$raw_data_dir, base_dir = cfg$config_dir, mustWork = FALSE)
  cfg$paths$output_dir <- resolve_path(cfg$paths$output_dir, base_dir = cfg$config_dir, mustWork = FALSE)
  cfg$paths$sample_sheet <- resolve_path(cfg$paths$sample_sheet, base_dir = cfg$config_dir, mustWork = FALSE)
  cfg$paths$python_path <- path.expand(cfg$paths$python_path %||% "")
  
  # Recover the real Scripts directory for older server layouts where the
  # config file lives in Scripts/Parameters and a relative "." would otherwise
  # resolve to the Parameters directory.
  scripts_dir_candidates <- unique(c(
    cfg$paths$scripts_dir,
    file.path(cfg$config_dir, ".."),
    file.path(cfg$paths$project_dir, "Scripts"),
    cfg$paths$project_dir
  ))
  scripts_dir_candidates <- vapply(
    scripts_dir_candidates,
    resolve_path,
    character(1),
    base_dir = cfg$config_dir,
    mustWork = FALSE
  )
  valid_scripts_dir <- scripts_dir_candidates[
    file.exists(file.path(scripts_dir_candidates, "00_Setup.R"))
  ][1]
  if (!is.na(valid_scripts_dir) && nzchar(valid_scripts_dir)) {
    cfg$paths$scripts_dir <- valid_scripts_dir
  }
  
  cfg$pipeline <- modifyList(
    list(
      mode = "all",
      save_intermediates = TRUE,
      overwrite_existing = FALSE,
      skip_on_error = FALSE
    ),
    cfg$pipeline %||% list()
  )
  
  cfg
}

normalize_sample_table <- function(sample_df, cfg) {
  sample_df <- as.data.frame(sample_df, stringsAsFactors = FALSE, check.names = FALSE)
  names(sample_df) <- standardize_colnames(names(sample_df))
  
  if (!"sample_id" %in% names(sample_df)) {
    if ("sample_id_raw" %in% names(sample_df)) {
      sample_df$sample_id <- sample_df$sample_id_raw
    } else if ("filename" %in% names(sample_df)) {
      sample_df$sample_id <- vapply(sample_df$filename, infer_sample_id, character(1))
    } else if ("directory_name" %in% names(sample_df)) {
      sample_df$sample_id <- vapply(sample_df$directory_name, infer_sample_id, character(1))
    } else {
      stop("Sample table must contain a sample_id or directory_name column.")
    }
  }
  
  if (!"directory_name" %in% names(sample_df)) {
    if ("filename" %in% names(sample_df)) {
      sample_df$directory_name <- sample_df$filename
    } else {
      sample_df$directory_name <- sample_df$sample_id
    }
  }
  
  if (!"subsample_id" %in% names(sample_df)) {
    sample_df$subsample_id <- NA_character_
  }
  if (!"slide_num" %in% names(sample_df)) {
    sample_df$slide_num <- NA_character_
  }
  if (!"include" %in% names(sample_df)) {
    sample_df$include <- TRUE
  }
  if (!"pair_id" %in% names(sample_df)) {
    sample_df$pair_id <- NA_character_
  }
  if (!"group_id" %in% names(sample_df)) {
    sample_df$group_id <- NA_character_
  }
  if (!"patient_id" %in% names(sample_df)) {
    sample_df$patient_id <- NA_character_
  }
  if (!"treatment" %in% names(sample_df)) {
    sample_df$treatment <- NA_character_
  }
  if (!"notes" %in% names(sample_df)) {
    sample_df$notes <- NA_character_
  }
  if (!"fov_min" %in% names(sample_df)) {
    sample_df$fov_min <- NA_character_
  }
  if (!"fov_max" %in% names(sample_df)) {
    sample_df$fov_max <- NA_character_
  }
  
  if (!"data_dir" %in% names(sample_df)) {
    sample_df$data_dir <- file.path(cfg$paths$raw_data_dir, sample_df$directory_name)
  }
  if (!"output_dir" %in% names(sample_df)) {
    sample_df$output_dir <- file.path(cfg$paths$output_dir, paste0("Sample_", sample_df$sample_id))
  }
  
  sample_df$data_dir <- vapply(sample_df$data_dir, resolve_path, character(1), base_dir = cfg$paths$raw_data_dir, mustWork = FALSE)
  sample_df$output_dir <- vapply(sample_df$output_dir, resolve_path, character(1), base_dir = cfg$paths$output_dir, mustWork = FALSE)
  sample_df$include <- vapply(sample_df$include, coerce_bool, logical(1), default = TRUE)
  
  sample_df <- sample_df[sample_df$include, , drop = FALSE]
  
  # Deduplicate composite-slide rows: when multiple rows share the same sample_id
  # (e.g. four CART biopsies on one slide each having a distinct subsample_id),
  # keep only the FIRST row per sample_id and clear fov_min/fov_max so the
  # whole slide is loaded in default (non-split) mode.
  # The full per-subsample rows are preserved in attr(, "split_rows") so that
  # run_pipeline() can restore them when --split is requested.
  dup_ids <- unique(sample_df$sample_id[duplicated(sample_df$sample_id)])
  if (length(dup_ids) > 0) {
    attr(sample_df, "split_rows") <- sample_df[sample_df$sample_id %in% dup_ids, , drop = FALSE]
    sample_df <- sample_df[!duplicated(sample_df$sample_id), , drop = FALSE]
    is_dup <- sample_df$sample_id %in% dup_ids
    sample_df$fov_min[is_dup] <- NA_character_
    sample_df$fov_max[is_dup] <- NA_character_
    sample_df$subsample_id[is_dup] <- NA_character_
  }
  
  rownames(sample_df) <- NULL
  sample_df
}

auto_detect_sample_table <- function(cfg) {
  raw_data_dir <- cfg$paths$raw_data_dir
  if (!dir.exists(raw_data_dir)) {
    stop("Raw data directory does not exist: ", raw_data_dir)
  }
  
  dirs <- list.dirs(raw_data_dir, recursive = FALSE, full.names = FALSE)
  if (length(dirs) == 0) {
    stop("No sample directories found in ", raw_data_dir)
  }
  
  sample_df <- data.frame(
    sample_id = vapply(dirs, infer_sample_id, character(1)),
    directory_name = dirs,
    slide_num = sub(".*Slide(\\d+).*", "\\1", dirs),
    stringsAsFactors = FALSE
  )
  normalize_sample_table(sample_df, cfg)
}

load_sample_table <- function(cfg, sample_sheet_override = NULL) {
  sheet_path <- sample_sheet_override %||% cfg$paths$sample_sheet
  if (!is.null(sheet_path) && file.exists(sheet_path)) {
    return(normalize_sample_table(read_tabular_file(sheet_path), cfg))
  }
  auto_detect_sample_table(cfg)
}

select_samples <- function(samples,
                           sample_ids = NULL,
                           pair_ids = NULL,
                           group_ids = NULL) {
  selected <- samples
  if (!is.null(sample_ids) && length(sample_ids) > 0) {
    selected <- selected[selected$sample_id %in% sample_ids, , drop = FALSE]
  }
  if (!is.null(pair_ids) && length(pair_ids) > 0 && "pair_id" %in% names(selected)) {
    selected <- selected[selected$pair_id %in% pair_ids, , drop = FALSE]
  }
  if (!is.null(group_ids) && length(group_ids) > 0 && "group_id" %in% names(selected)) {
    selected <- selected[selected$group_id %in% group_ids, , drop = FALSE]
  }
  
  rownames(selected) <- NULL
  selected
}

expand_split_samples <- function(samples, cfg) {
  # Replace rows that have a subsample_id with their per-subsample version:
  # sample_id  <- subsample_id
  # output_dir <- Sample_<subsample_id>
  # fov_min / fov_max are kept as-is (already set in the split_rows)
  if (!"subsample_id" %in% names(samples)) {
    return(samples)
  }
  has_sub <- !is.na(samples$subsample_id) & nzchar(as.character(samples$subsample_id))
  if (!any(has_sub)) {
    return(samples)
  }
  samples$sample_id[has_sub] <- as.character(samples$subsample_id[has_sub])
  samples$output_dir[has_sub] <- file.path(
    cfg$paths$output_dir,
    paste0("Sample_", samples$sample_id[has_sub])
  )
  samples$output_dir[has_sub] <- vapply(
    samples$output_dir[has_sub],
    resolve_path,
    character(1),
    base_dir = cfg$paths$output_dir,
    mustWork = FALSE
  )
  rownames(samples) <- NULL
  samples
}

compute_steps <- function(order, selected_steps = NULL, from_step = NULL) {
  steps <- if (is.null(selected_steps) || length(selected_steps) == 0) {
    order
  } else {
    order[order %in% selected_steps]
  }
  
  if (!is.null(from_step) && nzchar(from_step)) {
    idx <- match(from_step, order)
    if (is.na(idx)) {
      stop("Unknown step: ", from_step)
    }
    steps <- steps[match(steps, order) >= idx]
  }
  
  if (length(steps) == 0) {
    stop("No steps selected after applying filters.")
  }
  
  steps
}

checkpoint_dir_for_step <- function(root_dir, step_id) {
  file.path(root_dir, "pipeline_checkpoints", step_id)
}

native_sample_artifact <- function(sample_output_dir, sample_id, step_id) {
  switch(
    step_id,
    "01_load" = file.path(sample_output_dir, "Giotto_Object_Loaded"),
    "02_qc" = file.path(sample_output_dir, "Giotto_Object_Filtered"),
    "03_norm" = file.path(sample_output_dir, "Giotto_Object_Normalized"),
    "04_dimred" = file.path(sample_output_dir, "Giotto_Object_DimReduced"),
    "05_cluster" = file.path(sample_output_dir, "Giotto_Object_Clustered"),
    "06_markers" = file.path(sample_output_dir, "Giotto_Object_DEG_Markers"),
    "07_annotate" = file.path(sample_output_dir, "Giotto_Object_Annotated"),
    "08_visualize" = file.path(sample_output_dir, "Giotto_Object_Annotated"),
    "09_spatial" = file.path(sample_output_dir, "Giotto_Object_Spatial"),
    "10_cci" = file.path(sample_output_dir, "Giotto_Object_Spatial"),
    "11_bcell" = file.path(sample_output_dir, "Giotto_Object_BCell_Analysis"),
    "12_spatial_de" = file.path(sample_output_dir, "Giotto_Object_Spatial_DE"),
    NULL
  )
}

native_merged_artifact <- function(merged_output_dir, step_id) {
  switch(
    step_id,
    "merge" = file.path(merged_output_dir, "Giotto_Object_Merged"),
    "merge_batch" = file.path(merged_output_dir, "Giotto_Object_BatchCorrected"),
    "12_spatial_de" = file.path(merged_output_dir, "Giotto_Object_Spatial_DE"),
    NULL
  )
}

load_saved_object <- function(path) {
  if (is.null(path) || !file.exists(path)) {
    stop("Saved object path does not exist: ", path)
  }
  if (dir.exists(path) && file.exists(file.path(path, "manifest.json"))) {
    return(load_giotto_checkpoint(path))
  }
  if (dir.exists(path)) {
    return(loadGiotto(path))
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(readRDS(path))
  }
  if (grepl("\\.qs$", path, ignore.case = TRUE) && requireNamespace("qs", quietly = TRUE)) {
    return(qs::qread(path))
  }
  stop("Unsupported saved object path: ", path)
}

load_sample_checkpoint <- function(sample_row, step_id) {
  wrapped_path <- checkpoint_dir_for_step(sample_row$output_dir, step_id)
  if (dir.exists(wrapped_path) && file.exists(file.path(wrapped_path, "manifest.json"))) {
    return(load_giotto_checkpoint(wrapped_path))
  }
  
  native_path <- native_sample_artifact(sample_row$output_dir, sample_row$sample_id, step_id)
  if (!is.null(native_path) && file.exists(native_path)) {
    return(load_saved_object(native_path))
  }
  if (identical(step_id, "01_load")) {
    legacy_rds <- file.path(sample_row$output_dir, paste0(sample_row$sample_id, "_cosmx_loaded.rds"))
    if (file.exists(legacy_rds)) {
      return(readRDS(legacy_rds))
    }
  }
  
  stop(
    "No checkpoint found for sample '", sample_row$sample_id,
    "' at step ", step_id, "."
  )
}

load_merged_checkpoint <- function(merged_output_dir, step_id) {
  wrapped_path <- checkpoint_dir_for_step(merged_output_dir, step_id)
  if (dir.exists(wrapped_path) && file.exists(file.path(wrapped_path, "manifest.json"))) {
    return(load_giotto_checkpoint(wrapped_path))
  }
  
  native_path <- native_merged_artifact(merged_output_dir, step_id)
  if (!is.null(native_path) && file.exists(native_path)) {
    return(load_saved_object(native_path))
  }
  
  stop("No merged checkpoint found for step ", step_id, ".")
}

save_step_checkpoint <- function(gobj, root_dir, step_id, metadata = list()) {
  save_giotto_checkpoint(
    gobj = gobj,
    checkpoint_dir = checkpoint_dir_for_step(root_dir, step_id),
    metadata = c(list(step_id = step_id), metadata)
  )
}

# FIX #5: Replaced eval(parse(...)) with safe range/comma parsing.
coerce_dimension_vector <- function(x) {
  if (is.null(x)) {
    return(1:30)
  }
  if (is.numeric(x)) {
    return(as.integer(x))
  }
  # Parse "start:end" range notation safely without eval(parse(...))
  txt <- trimws(as.character(x))
  if (grepl("^\\d+:\\d+$", txt)) {
    bounds <- as.integer(strsplit(txt, ":")[[1]])
    return(seq.int(bounds[1], bounds[2]))
  }
  # Parse comma-separated integers: "1,2,3,10,15"
  if (grepl("^[\\d,\\s]+$", txt, perl = TRUE)) {
    return(as.integer(trimws(strsplit(txt, ",")[[1]])))
  }
  stop("Cannot parse dimension vector: '", txt,
       "'. Use 'start:end' (e.g. '1:30') or comma-separated integers (e.g. '1,2,3').")
}

flatten_marker_genes <- function(marker_cfg) {
  if (is.null(marker_cfg)) {
    return(NULL)
  }
  unique(unlist(marker_cfg, use.names = FALSE))
}

runtime_script_paths <- function() {
  file.path(
    repo_dir,
    c(
      "01_Load_data.R",
      "02_Quality_Control.R",
      "03_Normalisation.R",
      "04_Dimensionallity_Reduction.R",
      "05_Clustering.R",
      "06_Differential_Expression.R",
      "07_Annotation.R",
      "08_Visualisation.R",
      "09_Spatial_Network.R",
      "10_CCI_Analysis.R",
      file.path("Helper_Scripts", "Merge_Batch_Correction.R"),
      "11_B_Cell_Analysis.R",
      "12_Spatial_Differential_Expression.R"
    )
  )
}

ensure_runtime <- function(cfg) {
  runtime_env <- new.env(parent = globalenv())
  
  Sys.setenv(
    COSMX_PROJECT_DIR = cfg$paths$project_dir,
    COSMX_PYTHON_PATH = cfg$paths$python_path,
    COSMX_HELPER_DIR = file.path(cfg$paths$scripts_dir, "Helper_Scripts")
  )
  
  old_disable_cli <- getOption("cosmx.disable_cli")
  options(cosmx.disable_cli = TRUE)
  on.exit(options(cosmx.disable_cli = old_disable_cli), add = TRUE)
  
  source(file.path(repo_dir, "00_Setup.R"), local = runtime_env)
  for (script_path in runtime_script_paths()) {
    source(script_path, local = runtime_env)
  }
  
  runtime_env
}

invoke_sample_step <- function(runtime_env, step_id, gobj, sample_row, cfg) {
  sample_id <- sample_row$sample_id
  output_dir <- sample_row$output_dir
  
  switch(
    step_id,
    ## Changed 01_load after adding the subsampling information to the sample_sheet ##
    "01_load" = {
      gobj <- runtime_env$load_cosmx_sample(
        sample_id = sample_id,
        data_dir  = sample_row$data_dir,
        output_dir = output_dir
      )
      
      # FOV subsetting for composite slides (e.g. 4 biopsies on one CART slide)
      fov_cfg <- cfg$fov_split %||% list()
      if (isTRUE(fov_cfg$enabled %||% FALSE)) {
        fov_min_col <- fov_cfg$fov_min_col %||% "fov_min"
        fov_max_col <- fov_cfg$fov_max_col %||% "fov_max"
        fov_col     <- fov_cfg$fov_column  %||% "fov"
        fov_min_val <- suppressWarnings(as.integer(sample_row[[fov_min_col]]))
        fov_max_val <- suppressWarnings(as.integer(sample_row[[fov_max_col]]))
        
        if (!is.na(fov_min_val) && !is.na(fov_max_val)) {
          message(
            "FOV subsetting: keeping FOV ", fov_min_val, "-", fov_max_val,
            " for sample ", sample_id
          )
          meta <- as.data.frame(
            get("pDataDT", envir = asNamespace("Giotto"))(gobj),
            stringsAsFactors = FALSE
          )
          keep_fovs <- suppressWarnings(as.integer(meta[[fov_col]]))
          keep_ids  <- meta$cell_ID[
            !is.na(keep_fovs) & keep_fovs >= fov_min_val & keep_fovs <= fov_max_val
          ]
          if (length(keep_ids) == 0) {
            stop(
              "FOV subsetting removed all cells for sample '", sample_id,
              "'. Check fov_min/fov_max in sample_sheet.csv ",
              "(fov_min=", fov_min_val, ", fov_max=", fov_max_val, ")."
            )
          }
          gobj <- get("subsetGiotto", envir = asNamespace("Giotto"))(
            gobject  = gobj,
            cell_ids = keep_ids
          )
          message("  Retained ", length(keep_ids), " cells after FOV subsetting.")
        }
      }
      gobj
    },
    "02_qc" = runtime_env$quality_control(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      gene_min_cells = cfg$parameters$qc$gene_min_cells %||% 10,
      cell_min_genes = cfg$parameters$qc$cell_min_genes %||% 50,
      cell_max_genes = cfg$parameters$qc$cell_max_genes %||% NULL,
      min_count = cfg$parameters$qc$min_count %||% 100,
      max_mito_pct = cfg$parameters$qc$max_mito_pct %||% NULL
    ),
    "03_norm" = runtime_env$normalize_expression(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      scalefactor = cfg$parameters$normalization$scalefactor %||% 6000,
      log_transform = cfg$parameters$normalization$log_transform %||% TRUE
    ),
    "04_dimred" = runtime_env$dimensionality_reduction(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      n_hvgs = cfg$parameters$dim_reduction$n_hvgs %||% 500,
      n_pcs = cfg$parameters$dim_reduction$n_pcs %||% 30,
      umap_n_neighbors = cfg$parameters$dim_reduction$umap_n_neighbors %||% 30,
      umap_min_dist = cfg$parameters$dim_reduction$umap_min_dist %||% 0.3,
      spatial_hvg = cfg$parameters$dim_reduction$spatial_hvg %||% FALSE
    ),
    "05_cluster" = runtime_env$perform_clustering(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      k_nn = cfg$parameters$clustering$k_nn %||% 15,
      resolution = cfg$parameters$clustering$resolution %||% 0.3,
      dimensions_to_use = coerce_dimension_vector(cfg$parameters$clustering$dimensions_to_use %||% "1:30"),
      leiden_n_iterations = cfg$parameters$clustering$n_iterations %||% 200,
      resolution_sweep = cfg$parameters$clustering$resolution_sweep %||% NULL,
      scripts_dir = cfg$paths$scripts_dir
    ),
    "06_markers" = runtime_env$marker_analysis(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      cluster_column = cfg$parameters$markers$cluster_column %||% "leiden_clust",
      top_n = cfg$parameters$markers$top_n %||% 25
    ),
    "07_annotate" = runtime_env$annotate_cells(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      profiles = cfg$parameters$annotation$profiles,
      default_profile = cfg$parameters$annotation$default_profile %||% NULL,
      align_genes = cfg$parameters$annotation$align_genes %||% TRUE,
      n_starts = cfg$parameters$annotation$n_starts %||% 10,
      n_clusts_semi = cfg$parameters$annotation$n_clusts_semi %||% 3,
      cohort_column = cfg$parameters$annotation$cohort_column %||% "leiden_clust",
      min_gene_overlap = cfg$parameters$annotation$min_gene_overlap %||% 100,
      create_plots = cfg$parameters$annotation$create_plots %||% TRUE,
      conf_threshold = cfg$parameters$annotation$conf_threshold %||% NULL,
      save_object = TRUE
    ),
    "08_visualize" = runtime_env$create_visualizations(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      celltype_columns = cfg$parameters$visualization$celltype_columns %||% NULL,
      cluster_column = cfg$parameters$visualization$cluster_column %||% "leiden_clust",
      marker_genes = flatten_marker_genes(cfg$parameters$visualization$marker_genes),
      max_cells_preview = cfg$parameters$visualization$max_cells_preview %||% NULL,
      preview_seed = cfg$parameters$visualization$preview_seed %||% 1
    ),
    "09_spatial" = runtime_env$build_spatial_network(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      celltype_col = cfg$interaction$annotation_column %||% NULL,
      n_simulations = cfg$interaction$number_of_simulations %||% 250
    ),
    # FIX #2: run_cci_analysis() does NOT mutate the Giotto object.  All CCI
    # outputs (InSituCor modules, LIANA tables, NicheNet ligand activities,
    # MISTy models, nnSVG SVGs) are written to disk under
    # output_dir/10_CCI_Analysis/.  We intentionally return the unchanged
    # gobj so the pipeline checkpoint chain is unaffected.
    "10_cci" = {
      cci_cfg <- cfg$cci %||% list()
      # Resolve focus cell-type regex from config (mirrors step 11 at L812-825)
      # so LIANA plots (including the dedicated B-cell CCI network) centre on
      # the same focal cell type(s) used by downstream B-cell analysis.
      focus_regex <- cfg$interaction$focus_celltype_regex %||%
                     cfg$interaction$bcell_regex %||%
                     "^B\\.cell$"
      runtime_env$run_cci_analysis(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = output_dir,
        celltype_col = cfg$interaction$annotation_column %||% NULL,
        focus_celltype = focus_regex,
        sender_celltypes = cci_cfg$sender_celltypes %||% NULL,
        receiver_celltype = cci_cfg$receiver_celltype %||% NULL,
        target_genes = cci_cfg$target_genes %||% NULL,
        target_genes_by_receiver = cci_cfg$target_genes_by_receiver %||% NULL,
        nichenet_network_dir = cci_cfg$nichenet_network_dir %||% NULL,
        nichenet_mode = cci_cfg$nichenet_mode %||% "single",
        nichenet_spatial_filter = cci_cfg$nichenet_spatial_filter %||% FALSE,
        nichenet_proximity_enrichment_path = cci_cfg$nichenet_proximity_enrichment_path %||% NULL,
        nichenet_spatial_padj_threshold = cci_cfg$nichenet_spatial_padj_threshold %||% 0.05,
        nichenet_min_cells_per_celltype = cci_cfg$nichenet_min_cells_per_celltype %||% 5,
        nichenet_include_self_pairs = cci_cfg$nichenet_include_self_pairs %||% FALSE,
        run_sections = unlist(cci_cfg$run_sections %||% c(
          insitucor = TRUE,
          liana = TRUE,
          nichenet = TRUE,
          misty = TRUE,
          nnsvg = TRUE
        )),
        cleanup_between_sections = cci_cfg$cleanup_between_sections %||% TRUE
      )
      gobj  # CCI is file-based; gobj passes through unchanged
    },
    # FIX #1 + FIX #10: run_bcell_microenvironment_analysis() now returns the
    # (possibly modified) gobj directly, so the pipeline checkpoint chain
    # carries it forward.  focus_celltype_regex allows switching the focused
    # cell type via config without touching code.
    "11_bcell" = {
      focus_regex <- cfg$interaction$focus_celltype_regex %||%
        cfg$interaction$bcell_regex %||%
        "^B\\.cell$"
      runtime_env$run_bcell_microenvironment_analysis(
        gobj = gobj,
        sample_id = sample_id,
        output_dir = output_dir,
        annotation_column = cfg$interaction$annotation_column %||% NULL,
        bcell_regex = focus_regex,
        spatial_network_name = cfg$interaction$spatial_network_name %||% "Delaunay_network",
        number_of_simulations = cfg$interaction$number_of_simulations %||% 250,
        save_object = TRUE
      )
    },
    "12_spatial_de" = {
      spatial_de_cfg <- cfg$spatial_de %||% list()
      runtime_env$run_spatial_differential_expression(
        gobj = gobj,
        run_label = sample_id,
        output_dir = output_dir,
        analysis_scope = "sample",
        backend = spatial_de_cfg$sample_backend %||% "smiDE",
        annotation_column = spatial_de_cfg$annotation_column %||% NULL,
        sample_column = spatial_de_cfg$sample_column %||% "sample_id",
        treatment_column = spatial_de_cfg$treatment_column %||% "treatment",
        patient_column = spatial_de_cfg$patient_column %||% "patient_id",
        spatial_network_name = spatial_de_cfg$spatial_network_name %||% "Delaunay_network",
        n_niches = spatial_de_cfg$n_niches %||% 6,
        min_cells_per_niche = spatial_de_cfg$min_cells_per_niche %||% 30,
        min_cells_per_sample = spatial_de_cfg$min_cells_per_sample %||% 10,
        min_samples_per_group = spatial_de_cfg$min_samples_per_group %||% 2,
        sample_replicate_column = spatial_de_cfg$sample_replicate_column %||% "fov",
        sample_contrast = spatial_de_cfg$sample_contrast %||% "one_vs_rest",
        min_cells_per_replicate = spatial_de_cfg$min_cells_per_replicate %||% 5,
        min_replicates_per_group = spatial_de_cfg$min_replicates_per_group %||% 2,
        n_spatial_patches = spatial_de_cfg$n_spatial_patches %||% 8,
        smide_family = spatial_de_cfg$smide_family %||% "nbinom2",
        smide_radius = spatial_de_cfg$smide_radius %||% 0.05,
        smide_ncores = spatial_de_cfg$smide_ncores %||% 1,
        smide_overlap_threshold = spatial_de_cfg$smide_overlap_threshold %||% 1,
        smide_annotation_subset = spatial_de_cfg$smide_annotation_subset %||% NULL,
        smide_min_detection_fraction = spatial_de_cfg$smide_min_detection_fraction %||% 0.05,
        smide_custom_predictor = spatial_de_cfg$smide_custom_predictor %||% NULL,
        smide_partner_celltypes = spatial_de_cfg$smide_partner_celltypes %||% NULL,
        smide_partner_source = spatial_de_cfg$smide_partner_source %||% "none",
        smide_partner_top_n = spatial_de_cfg$smide_partner_top_n %||% 3,
        smide_partner_padj_threshold = spatial_de_cfg$smide_partner_padj_threshold %||% 0.05,
        smide_include_self_partner = spatial_de_cfg$smide_include_self_partner %||% FALSE,
        smide_save_raw = spatial_de_cfg$smide_save_raw %||% TRUE,
        save_object = TRUE
      )
    },
    stop("Unsupported sample step: ", step_id)
  )
}

run_sample_pipeline <- function(runtime_env,
                                samples,
                                cfg,
                                selected_steps = NULL,
                                from_step = NULL) {
  steps_to_run <- compute_steps(SAMPLE_STEP_ORDER, selected_steps, from_step)
  results <- vector("list", length = nrow(samples))
  
  for (idx in seq_len(nrow(samples))) {
    sample_row <- samples[idx, , drop = FALSE]
    sample_name <- sample_row$sample_id[[1]]
    sample_output_dir <- ensure_dir(sample_row$output_dir[[1]])
    sample_row$output_dir <- sample_output_dir
    
    message("\n=== Sample: ", sample_name, " ===")
    message("Steps: ", paste(vapply(steps_to_run, step_label, character(1)), collapse = " -> "))
    
    start_time <- Sys.time()
    status <- "SUCCESS"
    last_step <- NA_character_
    error_message <- NA_character_
    
    tryCatch({
      current_gobj <- NULL
      if (steps_to_run[1] != SAMPLE_STEP_ORDER[1]) {
        previous_step <- SAMPLE_STEP_ORDER[match(steps_to_run[1], SAMPLE_STEP_ORDER) - 1]
        current_gobj <- load_sample_checkpoint(sample_row, previous_step)
      }
      
      for (step_id in steps_to_run) {
        message("Running ", step_label(step_id), " for ", sample_name)
        current_gobj <- invoke_sample_step(runtime_env, step_id, current_gobj, sample_row, cfg)

        # After annotation: read auto-selection JSON and override annotation
        # columns for all downstream steps of this sample. cfg is a local
        # copy inside run_sample_pipeline(), so other samples are unaffected.
        if (step_id == "07_annotate") {
          selection_file <- file.path(
            sample_output_dir, "07_Annotation", "annotation_selection.json"
          )
          if (file.exists(selection_file)) {
            sel <- tryCatch(
              jsonlite::read_json(selection_file),
              error = function(e) {
                message("  [Auto-annotation] Could not read selection JSON: ",
                        conditionMessage(e))
                NULL
              }
            )
            if (!is.null(sel)) {
              sel_col <- sel$selected_annotation_column
              if (!is.null(sel_col) && nzchar(sel_col)) {
                cfg$interaction$annotation_column <- sel_col
                cfg$spatial_de$annotation_column  <- sel_col
                message(
                  "  [Auto-annotation] Downstream steps will use '", sel_col,
                  "' (composite score: ",
                  round(as.numeric(sel$composite_score), 3), ")"
                )
              }
            }
          }
        }

        if (isTRUE(cfg$pipeline$save_intermediates)) {
          save_step_checkpoint(
            current_gobj,
            root_dir = sample_output_dir,
            step_id = step_id,
            metadata = list(sample_id = sample_name)
          )
        }
        maybe_cleanup_runtime(cfg, label = paste(sample_name, step_id))
        last_step <- step_id
      }
    }, error = function(e) {
      status <<- "FAILED"
      error_message <<- conditionMessage(e)
      if (!isTRUE(cfg$pipeline$skip_on_error)) {
        stop(e)
      }
    })
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    results[[idx]] <- data.frame(
      sample_id = sample_name,
      status = status,
      last_step = last_step,
      output_dir = sample_output_dir,
      time_seconds = elapsed,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results)
}

merged_output_dir <- function(cfg, run_label = NULL) {
  label <- run_label %||% cfg$merged$run_label %||% format(Sys.time(), "%Y%m%d_%H%M%S")
  label <- gsub("[^A-Za-z0-9._-]+", "_", label)
  ensure_dir(file.path(cfg$paths$output_dir, "Merged", label))
}

load_merge_inputs <- function(samples, cfg) {
  input_step <- cfg$merged$input_step %||% "07_annotate"
  objs <- vector("list", length = nrow(samples))
  for (idx in seq_len(nrow(samples))) {
    gobj <- load_sample_checkpoint(samples[idx, , drop = FALSE], input_step)
    objs[[idx]] <- .strip_polygon_centroids(gobj)
  }
  objs
}

# Nullify spatVectorCentroids on all giottoPolygon slots so joinGiottoObjects
# does not encounter a centroid/polygon count mismatch from stale checkpoints.
.strip_polygon_centroids <- function(gobj) {
  tryCatch({
    poly_slot <- gobj@spatial_info
    if (!is.null(poly_slot) && length(poly_slot) > 0) {
      for (nm in names(poly_slot)) {
        poly <- poly_slot[[nm]]
        if (inherits(poly, "giottoPolygon")) {
          poly@spatVectorCentroids <- NULL
          gobj@spatial_info[[nm]] <- poly
        }
      }
    }
    gobj
  }, error = function(e) {
    message("  ⚠ Could not strip polygon centroids from checkpoint (non-fatal): ", conditionMessage(e))
    gobj
  })
}

run_merged_pipeline <- function(runtime_env,
                                samples,
                                cfg,
                                selected_steps = NULL,
                                from_step = NULL,
                                run_label = NULL) {
  steps_to_run <- compute_steps(MERGED_STEP_ORDER, selected_steps, from_step)
  output_dir <- merged_output_dir(cfg, run_label = run_label)
  merged_name <- basename(output_dir)
  
  status <- "SUCCESS"
  last_step <- NA_character_
  error_message <- NA_character_
  start_time <- Sys.time()
  
  tryCatch({
    current_gobj <- NULL
    if (steps_to_run[1] != MERGED_STEP_ORDER[1]) {
      previous_step <- MERGED_STEP_ORDER[match(steps_to_run[1], MERGED_STEP_ORDER) - 1]
      current_gobj <- load_merged_checkpoint(output_dir, previous_step)
    }
    
    for (step_id in steps_to_run) {
      message("Running ", step_label(step_id), " for merged run ", merged_name)
      if (step_id == "merge") {
        current_gobj <- runtime_env$merge_giotto_samples(
          gobject_list = load_merge_inputs(samples, cfg),
          sample_table = samples,
          output_dir = output_dir,
          join_method = cfg$merged$join_method %||% "shift",
          x_padding = cfg$merged$x_padding %||% 1000,
          save_object = TRUE
        )
      } else if (step_id == "merge_batch") {
        current_gobj <- runtime_env$batch_correct_merged_object(
          gobj = current_gobj,
          sample_id = merged_name,
          output_dir = output_dir,
          batch_column = cfg$merged$batch_column %||% "slide_id",
          n_hvgs = cfg$merged$n_hvgs %||% cfg$samples$dim_reduction$n_hvgs %||% 500,
          n_pcs = cfg$merged$n_pcs %||% cfg$samples$dim_reduction$n_pcs %||% 30,
          dimensions_to_use = cfg$merged$dimensions_to_use %||% "1:30",
          umap_n_neighbors = cfg$merged$umap_n_neighbors %||% 30,
          umap_min_dist = cfg$merged$umap_min_dist %||% 0.3,
          k_nn = cfg$merged$k_nn %||% 15,
          resolution = cfg$merged$resolution %||% 0.3,
          create_plots = TRUE,
          save_object = TRUE
        )
      } else if (step_id == "12_spatial_de") {
        spatial_de_cfg <- cfg$spatial_de %||% list()
        current_gobj <- runtime_env$run_spatial_differential_expression(
          gobj = current_gobj,
          run_label = merged_name,
          output_dir = output_dir,
          analysis_scope = "merged",
          backend = spatial_de_cfg$sample_backend %||% "smiDE",
          annotation_column = spatial_de_cfg$annotation_column %||% NULL,
          sample_column = spatial_de_cfg$sample_column %||% "sample_id",
          treatment_column = spatial_de_cfg$treatment_column %||% "treatment",
          patient_column = spatial_de_cfg$patient_column %||% "patient_id",
          spatial_network_name = spatial_de_cfg$spatial_network_name %||% "Delaunay_network",
          n_niches = spatial_de_cfg$n_niches %||% 6,
          min_cells_per_niche = spatial_de_cfg$min_cells_per_niche %||% 30,
          min_cells_per_sample = spatial_de_cfg$min_cells_per_sample %||% 10,
          min_samples_per_group = spatial_de_cfg$min_samples_per_group %||% 2,
          smide_family = spatial_de_cfg$smide_family %||% "nbinom2",
          smide_radius = spatial_de_cfg$smide_radius %||% 0.05,
          smide_ncores = spatial_de_cfg$smide_ncores %||% 1,
          smide_overlap_threshold = spatial_de_cfg$smide_overlap_threshold %||% 1,
          smide_annotation_subset = spatial_de_cfg$smide_annotation_subset %||% NULL,
          smide_min_detection_fraction = spatial_de_cfg$smide_min_detection_fraction %||% 0.05,
          smide_save_raw = spatial_de_cfg$smide_save_raw %||% TRUE,
          save_object = TRUE
        )
      } else {
        stop("Unsupported merged step: ", step_id)
      }
      
      if (isTRUE(cfg$pipeline$save_intermediates)) {
        save_step_checkpoint(
          current_gobj,
          root_dir = output_dir,
          step_id = step_id,
          metadata = list(run_label = merged_name)
        )
      }
      maybe_cleanup_runtime(cfg, label = paste("merged", merged_name, step_id))
      last_step <- step_id
    }
  }, error = function(e) {
    status <<- "FAILED"
    error_message <<- conditionMessage(e)
    if (!isTRUE(cfg$pipeline$skip_on_error)) {
      stop(e)
    }
  })
  
  data.frame(
    run_label = merged_name,
    status = status,
    last_step = last_step,
    output_dir = output_dir,
    time_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

write_results_csv <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, path, row.names = FALSE)
}

parse_cli_args <- function(args) {
  opts <- list(
    config = default_config_path(repo_dir),
    sample_sheet = NULL,
    mode = NULL,
    samples = NULL,
    pairs = NULL,
    groups = NULL,
    sample_steps = NULL,
    merged_steps = NULL,
    from_step = NULL,
    merged_from_step = NULL,
    run_label = NULL,
    dry_run = FALSE,
    split = FALSE
  )
  
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop("Unexpected positional argument: ", arg)
    }
    key <- sub("^--", "", arg)
    if (key == "dry-run") {
      opts$dry_run <- TRUE
      i <- i + 1L
      next
    }
    if (key == "split") {
      opts$split <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) {
      stop("Missing value for argument ", arg)
    }
    value <- args[[i + 1L]]
    i <- i + 2L
    
    if (key %in% c("samples", "pairs", "groups", "sample-steps", "merged-steps")) {
      target_name <- switch(
        key,
        "sample-steps" = "sample_steps",
        "merged-steps" = "merged_steps",
        key
      )
      opts[[target_name]] <- parse_cli_csv(value)
    } else if (key %in% c("sample-sheet", "from-step", "merged-from-step", "run-label")) {
      target_name <- switch(
        key,
        "sample-sheet" = "sample_sheet",
        "from-step" = "from_step",
        "merged-from-step" = "merged_from_step",
        "run-label" = "run_label"
      )
      opts[[target_name]] <- value
    } else if (key %in% c("config", "mode")) {
      opts[[key]] <- value
    } else {
      stop("Unknown argument: ", arg)
    }
  }
  
  opts
}

run_pipeline <- function(cli_opts) {
  cfg <- load_config(cli_opts$config)
  if (is.null(cli_opts$mode) || !nzchar(cli_opts$mode)) {
    cli_opts$mode <- cfg$pipeline$mode %||% "all"
  }
  
  sample_table <- load_sample_table(cfg, sample_sheet_override = cli_opts$sample_sheet)
  selected_samples <- select_samples(
    sample_table,
    sample_ids = cli_opts$samples,
    pair_ids = cli_opts$pairs,
    group_ids = cli_opts$groups
  )
  
  if (nrow(selected_samples) == 0) {
    stop("No samples matched the selected filters.")
  }
  
  # When --split is requested, expand composite-slide rows (those with a
  # subsample_id) into per-subsample rows, each with its own sample_id and
  # output_dir.  The full split_rows were preserved by normalize_sample_table()
  # as an attribute; we recover them here and re-select to honour --samples /
  # --groups / --pairs filters.
  if (isTRUE(cli_opts$split)) {
    split_rows <- attr(sample_table, "split_rows")
    if (!is.null(split_rows) && nrow(split_rows) > 0) {
      # Rows from split_rows whose parent sample_id was selected
      selected_split <- split_rows[split_rows$sample_id %in% selected_samples$sample_id, , drop = FALSE]
      # Rows that had no sub-splitting (no entry in split_rows)
      non_split <- selected_samples[!selected_samples$sample_id %in% split_rows$sample_id, , drop = FALSE]
      selected_samples <- rbind(non_split, selected_split)
      rownames(selected_samples) <- NULL
    }
    selected_samples <- expand_split_samples(selected_samples, cfg)
  }
  
  sample_steps <- compute_steps(
    SAMPLE_STEP_ORDER,
    canonical_step_ids(cli_opts$sample_steps, type = "sample"),
    canonical_step_ids(cli_opts$from_step, type = "sample")
  )
  merged_steps <- compute_steps(
    MERGED_STEP_ORDER,
    canonical_step_ids(cli_opts$merged_steps, type = "merged"),
    canonical_step_ids(cli_opts$merged_from_step, type = "merged")
  )
  
  if (isTRUE(cli_opts$dry_run)) {
    return(list(
      config = cfg$config_path,
      mode = cli_opts$mode,
      sample_count = nrow(selected_samples),
      sample_ids = selected_samples$sample_id,
      sample_steps = sample_steps,
      merged_steps = merged_steps,
      run_label = cli_opts$run_label %||% cfg$merged$run_label %||% NULL
    ))
  }
  
  runtime_env <- ensure_runtime(cfg)
  
  run_dir <- ensure_dir(file.path(
    cfg$paths$output_dir,
    "pipeline_runs",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  ))
  
  sample_results <- NULL
  merged_results <- NULL
  
  if (cli_opts$mode %in% c("all", "separate", "paired")) {
    sample_results <- run_sample_pipeline(
      runtime_env = runtime_env,
      samples = selected_samples,
      cfg = cfg,
      selected_steps = sample_steps,
      from_step = cli_opts$from_step
    )
    write_results_csv(sample_results, file.path(run_dir, "sample_pipeline_results.csv"))
  }
  
  if (cli_opts$mode %in% c("all", "merged")) {
    merged_results <- run_merged_pipeline(
      runtime_env = runtime_env,
      samples = selected_samples,
      cfg = cfg,
      selected_steps = merged_steps,
      from_step = cli_opts$merged_from_step,
      run_label = cli_opts$run_label
    )
    write_results_csv(merged_results, file.path(run_dir, "merged_pipeline_results.csv"))
  }
  
  invisible(list(
    run_dir = run_dir,
    sample_results = sample_results,
    merged_results = merged_results
  ))
}

if (!interactive() && !isTRUE(getOption("cosmx.disable_cli", FALSE))) {
  cli_opts <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  result <- run_pipeline(cli_opts)
  if (isTRUE(cli_opts$dry_run)) {
    print(result)
  }
}
# Manual single-sample pipeline runner for RStudio on the server

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

read_config_if_present <- function(config_path) {
  if (is.null(config_path) || !file.exists(config_path)) {
    return(NULL)
  }
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("yaml is required to read the pipeline config: ", config_path)
  }
  yaml::read_yaml(config_path)
}

flatten_marker_genes <- function(marker_config) {
  if (is.null(marker_config)) {
    return(NULL)
  }
  unique(as.character(unlist(marker_config, recursive = TRUE, use.names = FALSE)))
}

step_ids <- c(
  "01_load",
  "02_qc",
  "03_normalize",
  "04_dimred",
  "05_cluster",
  "06_markers",
  "07_annotation",
  "08_visualization",
  "09_spatial_network",
  "10_cci",
  "11_bcell",
  "12_spatial_de"
)

step_labels <- c(
  "01_load" = "Load Data",
  "02_qc" = "Quality Control",
  "03_normalize" = "Normalization",
  "04_dimred" = "Dimensionality Reduction",
  "05_cluster" = "Clustering",
  "06_markers" = "Marker Analysis",
  "07_annotation" = "Annotation",
  "08_visualization" = "Visualization",
  "09_spatial_network" = "Spatial Network",
  "10_cci" = "CCI Analysis",
  "11_bcell" = "B-Cell Microenvironment",
  "12_spatial_de" = "Spatial Differential Expression"
)

coerce_step_id <- function(step_id) {
  if (is.null(step_id) || !nzchar(step_id)) {
    return(NULL)
  }
  step_id <- trimws(as.character(step_id[[1]]))
  if (step_id %in% step_ids) {
    return(step_id)
  }
  numeric_candidate <- suppressWarnings(as.integer(step_id))
  if (!is.na(numeric_candidate) && numeric_candidate >= 1 && numeric_candidate <= length(step_ids)) {
    return(step_ids[[numeric_candidate]])
  }
  stop(
    "Unknown step identifier: ", step_id, "\n",
    "Valid values are: ", paste(step_ids, collapse = ", ")
  )
}

previous_step_id <- function(step_id) {
  idx <- match(step_id, step_ids)
  if (is.na(idx) || idx <= 1) {
    return(NULL)
  }
  step_ids[[idx - 1]]
}

should_run_step <- function(step_id, start_index, stop_index) {
  idx <- match(step_id, step_ids)
  !is.na(idx) && idx >= start_index && idx <= stop_index
}

checkpoint_step_dir <- function(checkpoint_root, step_id) {
  file.path(checkpoint_root, step_id)
}

append_manual_log <- function(log_file, ..., .sep = "") {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = .sep))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

write_step_status <- function(status_file,
                              step_id,
                              status,
                              note = NULL,
                              checkpoint_dir = NULL) {
  status_df <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    step_id = step_id,
    step_label = unname(step_labels[[step_id]] %||% step_id),
    status = status,
    checkpoint_dir = checkpoint_dir %||% NA_character_,
    note = note %||% NA_character_,
    stringsAsFactors = FALSE
  )
  utils::write.table(
    status_df,
    file = status_file,
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(status_file),
    append = file.exists(status_file),
    qmethod = "double"
  )
}

save_manual_checkpoint <- function(gobj,
                                   checkpoint_root,
                                   step_id,
                                   sample_id,
                                   output_dir,
                                   log_file,
                                   status_file) {
  if (!exists("save_giotto_checkpoint")) {
    stop("save_giotto_checkpoint() is not available. Make sure 00_Setup.R loaded Helper_Scripts/Pipeline_Utils.R")
  }
  checkpoint_dir <- checkpoint_step_dir(checkpoint_root, step_id)
  result <- save_giotto_checkpoint(
    gobj = gobj,
    checkpoint_dir = checkpoint_dir,
    overwrite = TRUE,
    metadata = list(
      sample_id = sample_id,
      output_dir = output_dir,
      step_id = step_id,
      step_label = unname(step_labels[[step_id]] %||% step_id)
    )
  )
  append_manual_log(
    log_file,
    "Checkpoint saved after ", step_id,
    " using method: ", result$save_method %||% "unknown",
    " (", checkpoint_dir, ")"
  )
  write_step_status(
    status_file = status_file,
    step_id = step_id,
    status = "checkpoint_saved",
    checkpoint_dir = checkpoint_dir,
    note = paste("save_method =", result$save_method %||% "unknown")
  )
  invisible(checkpoint_dir)
}

run_gobject_step <- function(step_id,
                             expr,
                             checkpoint_root,
                             sample_id,
                             output_dir,
                             save_checkpoints,
                             log_file,
                             status_file) {
  step_label <- unname(step_labels[[step_id]] %||% step_id)
  append_manual_log(log_file, "Starting ", step_id, " | ", step_label)
  write_step_status(status_file, step_id, "started")
  step_start <- Sys.time()
  gobj <- force(expr)
  elapsed_minutes <- round(as.numeric(difftime(Sys.time(), step_start, units = "mins")), 2)
  append_manual_log(log_file, "Completed ", step_id, " | ", step_label, " in ", elapsed_minutes, " min")
  write_step_status(
    status_file = status_file,
    step_id = step_id,
    status = "completed",
    note = paste("elapsed_min =", elapsed_minutes)
  )
  if (isTRUE(save_checkpoints)) {
    save_manual_checkpoint(
      gobj = gobj,
      checkpoint_root = checkpoint_root,
      step_id = step_id,
      sample_id = sample_id,
      output_dir = output_dir,
      log_file = log_file,
      status_file = status_file
    )
  }
  gobj
}

run_side_effect_step <- function(step_id,
                                 expr,
                                 gobj,
                                 checkpoint_root,
                                 sample_id,
                                 output_dir,
                                 save_checkpoints,
                                 log_file,
                                 status_file) {
  step_label <- unname(step_labels[[step_id]] %||% step_id)
  append_manual_log(log_file, "Starting ", step_id, " | ", step_label)
  write_step_status(status_file, step_id, "started")
  step_start <- Sys.time()
  result <- force(expr)
  elapsed_minutes <- round(as.numeric(difftime(Sys.time(), step_start, units = "mins")), 2)
  append_manual_log(log_file, "Completed ", step_id, " | ", step_label, " in ", elapsed_minutes, " min")
  write_step_status(
    status_file = status_file,
    step_id = step_id,
    status = "completed",
    note = paste("elapsed_min =", elapsed_minutes)
  )
  if (isTRUE(save_checkpoints) && !is.null(gobj)) {
    save_manual_checkpoint(
      gobj = gobj,
      checkpoint_root = checkpoint_root,
      step_id = step_id,
      sample_id = sample_id,
      output_dir = output_dir,
      log_file = log_file,
      status_file = status_file
    )
  }
  result
}

# ------------------------------------------------------------------------------
# Server-specific paths
# ------------------------------------------------------------------------------
analysis_dir <- normalizePath(
  Sys.getenv("COSMX_ANALYSIS_DIR", unset = "/mnt/home/koncha/P_lab/CosMx_analysis"),
  winslash = "/",
  mustWork = TRUE
)

scripts_dir <- if (file.exists(file.path(analysis_dir, "Scripts", "00_Setup.R"))) {
  file.path(analysis_dir, "Scripts")
} else {
  analysis_dir
}
scripts_dir <- normalizePath(scripts_dir, winslash = "/", mustWork = TRUE)

parameters_dir <- if (dir.exists(file.path(analysis_dir, "Parameters"))) {
  file.path(analysis_dir, "Parameters")
} else if (dir.exists(file.path(scripts_dir, "Parameters"))) {
  file.path(scripts_dir, "Parameters")
} else {
  analysis_dir
}
parameters_dir <- normalizePath(parameters_dir, winslash = "/", mustWork = FALSE)

helper_dir <- file.path(scripts_dir, "Helper_Scripts")
config_path <- file.path(parameters_dir, "config.yaml")
config <- read_config_if_present(config_path) %||% list()

Sys.setenv(COSMX_PROJECT_DIR = analysis_dir)
if (dir.exists(helper_dir)) {
  Sys.setenv(COSMX_HELPER_DIR = normalizePath(helper_dir, winslash = "/", mustWork = TRUE))
}
if (!is.null(config$paths$python_path) && nzchar(config$paths$python_path)) {
  Sys.setenv(COSMX_PYTHON_PATH = config$paths$python_path)
  Sys.setenv(RETICULATE_PYTHON = config$paths$python_path)
}

# ------------------------------------------------------------------------------
# Sample-specific settings
# ------------------------------------------------------------------------------
sample_id <- "18H12037"
data_dir <- normalizePath(
  Sys.getenv(
    "COSMX_SAMPLE_DIR",
    unset = file.path(analysis_dir, "Data", "Raw_data", "MCFASCLST1196_Slide418H12037")
  ),
  winslash = "/",
  mustWork = FALSE
)
output_dir <- file.path(analysis_dir, "Output", paste0("Sample_", sample_id))
checkpoint_root <- file.path(output_dir, "00_Checkpoints")
log_dir <- file.path(output_dir, "00_Manual_Run_Logs")

save_intermediate_checkpoints <- TRUE
start_from_step <- "01_load"
stop_after_step <- NULL

run_marker_analysis_step <- TRUE
run_visualization_step <- TRUE
run_cci_step <- TRUE
run_bcell_step <- TRUE
run_spatial_de_step <- TRUE

required_script_paths <- c(
  file.path(scripts_dir, "00_Setup.R"),
  file.path(scripts_dir, "01_Load_data.R"),
  file.path(scripts_dir, "02_Quality_Control.R"),
  file.path(scripts_dir, "03_Normalisation.R"),
  file.path(scripts_dir, "04_Dimensionallity_Reduction.R"),
  file.path(scripts_dir, "05_Clustering.R"),
  file.path(scripts_dir, "06_Differential_Expression.R"),
  file.path(scripts_dir, "07_Annotation.R"),
  file.path(scripts_dir, "08_Visualisation.R"),
  file.path(scripts_dir, "09_Spatial_Network.R"),
  file.path(scripts_dir, "10_CCI_Analysis.R"),
  file.path(scripts_dir, "11_B_Cell_Analysis.R"),
  file.path(scripts_dir, "12_Spatial_Differential_Expression.R")
)
missing_script_paths <- required_script_paths[!file.exists(required_script_paths)]
if (length(missing_script_paths) > 0) {
  stop(
    "The following required scripts are missing on the server:\n - ",
    paste(missing_script_paths, collapse = "\n - ")
  )
}

source(file.path(scripts_dir, "00_Setup.R"))
source(file.path(scripts_dir, "01_Load_data.R"))
source(file.path(scripts_dir, "02_Quality_Control.R"))
source(file.path(scripts_dir, "03_Normalisation.R"))
source(file.path(scripts_dir, "04_Dimensionallity_Reduction.R"))
source(file.path(scripts_dir, "05_Clustering.R"))
source(file.path(scripts_dir, "06_Differential_Expression.R"))
source(file.path(scripts_dir, "07_Annotation.R"))
source(file.path(scripts_dir, "08_Visualisation.R"))
source(file.path(scripts_dir, "09_Spatial_Network.R"))
source(file.path(scripts_dir, "10_CCI_Analysis.R"))
source(file.path(scripts_dir, "11_B_Cell_Analysis.R"))
source(file.path(scripts_dir, "12_Spatial_Differential_Expression.R"))

if (!dir.exists(data_dir)) {
  stop(
    "data_dir does not exist: ", data_dir, "\n",
    "Set COSMX_SAMPLE_DIR or edit data_dir to the absolute raw-data folder for sample ",
    sample_id, "."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_root, recursive = TRUE, showWarnings = FALSE)

run_log_file <- file.path(log_dir, paste0(sample_id, "_manual_run.log"))
step_status_file <- file.path(log_dir, paste0(sample_id, "_step_status.csv"))

start_from_step <- coerce_step_id(start_from_step %||% "01_load")
stop_after_step <- coerce_step_id(stop_after_step)
start_index <- match(start_from_step, step_ids)
stop_index <- if (is.null(stop_after_step)) length(step_ids) else match(stop_after_step, step_ids)

if (stop_index < start_index) {
  stop("stop_after_step must be the same as or later than start_from_step")
}

cat("\n========================================\n")
cat("Manual Single-Sample Test\n")
cat("========================================\n")
cat("Analysis dir: ", analysis_dir, "\n")
cat("Scripts dir:  ", scripts_dir, "\n")
cat("Config:       ", if (file.exists(config_path)) config_path else "<none>", "\n")
cat("Sample:       ", sample_id, "\n")
cat("Data dir:     ", data_dir, "\n")
cat("Output dir:   ", output_dir, "\n")
cat("Checkpoints:  ", checkpoint_root, "\n")
cat("Run log:      ", run_log_file, "\n")
cat("Step status:  ", step_status_file, "\n")
cat("Start step:   ", start_from_step, "\n")
cat("Stop step:    ", if (is.null(stop_after_step)) "<end>" else stop_after_step, "\n")
cat("========================================\n\n")

append_manual_log(run_log_file, paste(rep("=", 72), collapse = ""))
append_manual_log(
  run_log_file,
  "Run started for sample ", sample_id,
  " | start_from_step = ", start_from_step,
  " | stop_after_step = ", if (is.null(stop_after_step)) "<end>" else stop_after_step
)

annotation_cfg <- config$parameters$annotation %||% list()
clustering_cfg <- config$parameters$clustering %||% list()
visual_cfg <- config$parameters$visualization %||% list()
marker_cfg <- config$parameters$markers %||% list()
cci_cfg <- config$cci %||% list()
interaction_cfg <- config$interaction %||% list()
spatial_de_cfg <- config$spatial_de %||% list()

if (is.null(annotation_cfg$profiles) || length(annotation_cfg$profiles) == 0) {
  stop("No annotation profiles were found in the config file. Step 07 cannot run.")
}

if (isTRUE(save_intermediate_checkpoints)) {
  if (!exists("save_giotto_checkpoint") || !exists("load_giotto_checkpoint")) {
    stop(
      "Checkpoint helpers are unavailable. Make sure 00_Setup.R loaded Helper_Scripts/Pipeline_Utils.R"
    )
  }
}

gobj <- NULL
if (start_index > 1) {
  resume_after_step <- previous_step_id(start_from_step)
  resume_checkpoint <- checkpoint_step_dir(checkpoint_root, resume_after_step)
  if (!dir.exists(resume_checkpoint)) {
    stop(
      "Requested start_from_step = ", start_from_step,
      " but the required checkpoint does not exist: ", resume_checkpoint, "\n",
      "Run the pipeline up to ", resume_after_step, " first, or set start_from_step = '01_load'."
    )
  }
  append_manual_log(
    run_log_file,
    "Loading checkpoint from ", resume_checkpoint,
    " before resuming at ", start_from_step
  )
  gobj <- load_giotto_checkpoint(resume_checkpoint)
  write_step_status(
    status_file = step_status_file,
    step_id = start_from_step,
    status = "resume_loaded",
    checkpoint_dir = resume_checkpoint,
    note = paste("loaded checkpoint saved after", resume_after_step)
  )
}

# Step 01
if (should_run_step("01_load", start_index, stop_index)) {
  gobj <- run_gobject_step(
    step_id = "01_load",
    expr = load_cosmx_sample(
      sample_id = sample_id,
      data_dir = data_dir,
      output_dir = output_dir
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 02
if (should_run_step("02_qc", start_index, stop_index)) {
  gobj <- run_gobject_step(
    step_id = "02_qc",
    expr = quality_control(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 03
if (should_run_step("03_normalize", start_index, stop_index)) {
  gobj <- run_gobject_step(
    step_id = "03_normalize",
    expr = normalize_expression(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 04
if (should_run_step("04_dimred", start_index, stop_index)) {
  gobj <- run_gobject_step(
    step_id = "04_dimred",
    expr = dimensionality_reduction(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 05
if (should_run_step("05_cluster", start_index, stop_index)) {
  gobj <- run_gobject_step(
    step_id = "05_cluster",
    expr = perform_clustering(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      leiden_n_iterations = clustering_cfg$n_iterations %||% 200,
      resolution_sweep = clustering_cfg$resolution_sweep %||% NULL,
      scripts_dir = scripts_dir
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 06
if (should_run_step("06_markers", start_index, stop_index) && isTRUE(run_marker_analysis_step)) {
  gobj <- run_gobject_step(
    step_id = "06_markers",
    expr = marker_analysis(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      cluster_column = marker_cfg$cluster_column %||% "leiden_clust",
      top_n = marker_cfg$top_n %||% 25
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 07
if (should_run_step("07_annotation", start_index, stop_index)) {
  gobj <- run_gobject_step(
    step_id = "07_annotation",
    expr = annotate_cells(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      profiles = annotation_cfg$profiles,
      default_profile = annotation_cfg$default_profile %||% annotation_cfg$profiles[[1]]$name,
      align_genes = annotation_cfg$align_genes %||% TRUE,
      n_clusts_semi = annotation_cfg$n_clusts_semi %||% 5,
      n_starts = annotation_cfg$n_starts %||% 10,
      cohort_column = annotation_cfg$cohort_column %||% "leiden_clust",
      min_gene_overlap = annotation_cfg$min_gene_overlap %||% 100,
      create_plots = annotation_cfg$create_plots %||% TRUE,
      conf_threshold = annotation_cfg$conf_threshold %||% NULL,
      save_object = FALSE
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 08
if (should_run_step("08_visualization", start_index, stop_index) && isTRUE(run_visualization_step)) {
  gobj <- run_gobject_step(
    step_id = "08_visualization",
    expr = create_visualizations(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      celltype_columns = visual_cfg$celltype_columns %||% NULL,
      cluster_column = visual_cfg$cluster_column %||% "leiden_clust",
      marker_genes = flatten_marker_genes(visual_cfg$marker_genes),
      max_cells_preview = visual_cfg$max_cells_preview %||% NULL,
      preview_seed = visual_cfg$preview_seed %||% 1
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 09
if (should_run_step("09_spatial_network", start_index, stop_index)) {
  gobj <- run_gobject_step(
    step_id = "09_spatial_network",
    expr = build_spatial_network(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      celltype_col = spatial_de_cfg$annotation_column %||% interaction_cfg$annotation_column %||% NULL
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 10
if (should_run_step("10_cci", start_index, stop_index) && isTRUE(run_cci_step)) {
  cci_results <- run_side_effect_step(
    step_id = "10_cci",
    expr = run_cci_analysis(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      celltype_col = interaction_cfg$annotation_column %||% spatial_de_cfg$annotation_column %||% NULL,
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
        nichenet = FALSE,
        misty = TRUE,
        nnsvg = FALSE
      )),
      cleanup_between_sections = cci_cfg$cleanup_between_sections %||% TRUE
    ),
    gobj = gobj,
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 11
if (should_run_step("11_bcell", start_index, stop_index) && isTRUE(run_bcell_step)) {
  bcell_results <- run_side_effect_step(
    step_id = "11_bcell",
    expr = run_bcell_microenvironment_analysis(
      gobj = gobj,
      sample_id = sample_id,
      output_dir = output_dir,
      annotation_column = interaction_cfg$annotation_column %||% spatial_de_cfg$annotation_column %||% NULL,
      bcell_regex = interaction_cfg$bcell_regex %||% "^B\\.cell$",
      spatial_network_name = interaction_cfg$spatial_network_name %||% "Delaunay_network",
      number_of_simulations = interaction_cfg$number_of_simulations %||% 250,
      save_object = FALSE
    ),
    gobj = gobj,
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

# Step 12
if (should_run_step("12_spatial_de", start_index, stop_index) && isTRUE(run_spatial_de_step)) {
  gobj <- run_gobject_step(
    step_id = "12_spatial_de",
    expr = run_spatial_differential_expression(
      gobj = gobj,
      run_label = sample_id,
      output_dir = output_dir,
      analysis_scope = "sample",
      backend = spatial_de_cfg$sample_backend %||% "smiDE",
      annotation_column = spatial_de_cfg$annotation_column %||% interaction_cfg$annotation_column %||% NULL,
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
      save_object = FALSE
    ),
    checkpoint_root = checkpoint_root,
    sample_id = sample_id,
    output_dir = output_dir,
    save_checkpoints = save_intermediate_checkpoints,
    log_file = run_log_file,
    status_file = step_status_file
  )
}

if (!is.null(gobj)) {
  saveGiotto(
    gobject = gobj,
    dir = output_dir,
    foldername = "Giotto_Object_Manual_Test",
    overwrite = TRUE
  )
}

append_manual_log(run_log_file, "Manual single-sample pipeline complete for ", sample_id)
if (!is.null(stop_after_step)) {
  append_manual_log(run_log_file, "Run stopped after requested step ", stop_after_step)
}

cat("\n✓ Manual single-sample pipeline complete for ", sample_id, "\n", sep = "")

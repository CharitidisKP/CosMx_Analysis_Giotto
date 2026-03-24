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
config <- read_config_if_present(config_path)

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

run_marker_analysis_step <- TRUE
run_visualization_step <- TRUE
run_cci_step <- TRUE
run_bcell_step <- TRUE
run_spatial_de_step <- TRUE

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
source(file.path(scripts_dir, "12_B_Cell_Analysis.R"))
source(file.path(scripts_dir, "13_Spatial_Differential_Expression.R"))

if (!dir.exists(data_dir)) {
  stop(
    "data_dir does not exist: ", data_dir, "\n",
    "Set COSMX_SAMPLE_DIR or edit data_dir to the absolute raw-data folder for sample ",
    sample_id, "."
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========================================\n")
cat("Manual Single-Sample Test\n")
cat("========================================\n")
cat("Analysis dir: ", analysis_dir, "\n")
cat("Scripts dir:  ", scripts_dir, "\n")
cat("Config:       ", if (file.exists(config_path)) config_path else "<none>", "\n")
cat("Sample:       ", sample_id, "\n")
cat("Data dir:     ", data_dir, "\n")
cat("Output dir:   ", output_dir, "\n")
cat("========================================\n\n")

annotation_cfg <- config$parameters$annotation %||% list()
visual_cfg <- config$parameters$visualization %||% list()
marker_cfg <- config$parameters$markers %||% list()
cci_cfg <- config$cci %||% list()
interaction_cfg <- config$interaction %||% list()
spatial_de_cfg <- config$spatial_de %||% list()

if (is.null(annotation_cfg$profiles) || length(annotation_cfg$profiles) == 0) {
  stop("No annotation profiles were found in the config file. Step 07 cannot run.")
}

# Step 01
gobj <- load_cosmx_sample(
  sample_id = sample_id,
  data_dir = data_dir,
  output_dir = output_dir
)

# Step 02
gobj <- quality_control(
  gobj = gobj,
  sample_id = sample_id,
  output_dir = output_dir
)

# Step 03
gobj <- normalize_expression(
  gobj = gobj,
  sample_id = sample_id,
  output_dir = output_dir
)

# Step 04
gobj <- dimensionality_reduction(
  gobj = gobj,
  sample_id = sample_id,
  output_dir = output_dir
)

# Step 05
gobj <- perform_clustering(
  gobj = gobj,
  sample_id = sample_id,
  output_dir = output_dir, 
  inspect_snn = FALSE
)

# Step 06
if (isTRUE(run_marker_analysis_step)) {
  gobj <- marker_analysis(
    gobj = gobj,
    sample_id = sample_id,
    output_dir = output_dir,
    cluster_column = marker_cfg$cluster_column %||% "leiden_clust",
    top_n = marker_cfg$top_n %||% 25
  )
}

# Step 07
gobj <- annotate_cells(
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
  save_object = TRUE
)

# Step 08
if (isTRUE(run_visualization_step)) {
  gobj <- create_visualizations(
    gobj = gobj,
    sample_id = sample_id,
    output_dir = output_dir,
    celltype_columns = visual_cfg$celltype_columns %||% NULL,
    cluster_column = visual_cfg$cluster_column %||% "leiden_clust",
    marker_genes = flatten_marker_genes(visual_cfg$marker_genes)
  )
}

# Step 09
gobj <- build_spatial_network(
  gobj = gobj,
  sample_id = sample_id,
  output_dir = output_dir,
  celltype_col = spatial_de_cfg$annotation_column %||% interaction_cfg$annotation_column %||% NULL
)

# Step 10
if (isTRUE(run_cci_step)) {
  cci_results <- run_cci_analysis(
    gobj = gobj,
    sample_id = sample_id,
    output_dir = output_dir,
    celltype_col = interaction_cfg$annotation_column %||% spatial_de_cfg$annotation_column %||% NULL,
    sender_celltypes = cci_cfg$sender_celltypes %||% NULL,
    receiver_celltype = cci_cfg$receiver_celltype %||% NULL,
    nichenet_network_dir = cci_cfg$nichenet_network_dir %||% NULL,
    run_sections = cci_cfg$run_sections %||% c(
      insitucor = TRUE,
      liana = TRUE,
      nichenet = FALSE,
      misty = TRUE,
      nnsvg = TRUE
    )
  )
}

# Step 12
if (isTRUE(run_bcell_step)) {
  bcell_results <- run_bcell_microenvironment_analysis(
    gobj = gobj,
    sample_id = sample_id,
    output_dir = output_dir,
    annotation_column = interaction_cfg$annotation_column %||% spatial_de_cfg$annotation_column %||% NULL,
    bcell_regex = interaction_cfg$bcell_regex %||% "B|Plasma|Plasmablast",
    spatial_network_name = interaction_cfg$spatial_network_name %||% "Delaunay_network",
    number_of_simulations = interaction_cfg$number_of_simulations %||% 250,
    save_object = FALSE
  )
}

# Step 13
if (isTRUE(run_spatial_de_step)) {
  gobj <- run_spatial_differential_expression(
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
    smide_overlap_threshold = spatial_de_cfg$smide_overlap_threshold %||% NULL,
    smide_save_raw = spatial_de_cfg$smide_save_raw %||% TRUE,
    save_object = FALSE
  )
}

saveGiotto(
  gobject = gobj,
  dir = output_dir,
  foldername = "Giotto_Object_Manual_Test",
  overwrite = TRUE
)

cat("\n✓ Manual single-sample pipeline complete for ", sample_id, "\n", sep = "")

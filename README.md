# CosMx / Giotto Analysis Pipeline

Spatially resolved transcriptomic analysis of NanoString CosMx data using the
[Giotto](https://giottosuite.readthedocs.io) suite.  The pipeline processes
raw CosMx output directories through QC, normalisation, dimensionality
reduction, clustering, annotation, cell-cell interaction analysis, and spatial
differential expression.  All steps are orchestrated by a single R driver
(`Scripts/CosMx_pipeline.R`) launched via the Apptainer wrapper script
(`Scripts/Run_Scripts/Run_Giotto_Pipeline.sh`).

---

## Repository layout

```
CosMx_analysis/
├── Scripts/
│   ├── CosMx_pipeline.R               # Main pipeline driver
│   ├── 00_Setup.R                     # Environment / library setup
│   ├── 01_Load_data.R                 # Step 01 – Data loading & FOV subsetting
│   ├── 02_Quality_Control.R           # Step 02 – QC filtering
│   ├── 03_Normalisation.R             # Step 03 – Normalisation
│   ├── 04_Dimensionallity_Reduction.R # Step 04 – PCA / UMAP
│   ├── 05_Clustering.R                # Step 05 – Leiden clustering
│   ├── 06_Differential_Expression.R   # Step 06 – Marker genes
│   ├── 07_Annotation.R                # Step 07 – Cell-type annotation
│   ├── 08_Visualisation.R             # Step 08 – Plots & feature maps
│   ├── 09_Spatial_Network.R           # Step 09 – Delaunay spatial network
│   ├── 10_CCI_Analysis.R              # Step 10 – Cell–cell interaction (CCI)
│   ├── 11_B_Cell_Analysis.R           # Step 11 – Focused immune microenvironment
│   ├── 12_Spatial_Differential_Expression.R  # Step 12 – Spatial DE
│   ├── Helper_Scripts/
│   │   ├── Merge_Batch_Correction.R   # Merge + Harmony batch correction
│   │   ├── Pipeline_Utils.R           # Shared utilities
│   │   └── ...
│   ├── Parameters/
│   │   ├── config.yaml                # All pipeline parameters
│   │   └── sample_sheet.csv           # Sample manifest
│   └── Run_Scripts/
│       └── Run_Giotto_Pipeline.sh     # Apptainer launch wrapper
├── Data/
│   └── Raw_data/                      # One subdirectory per physical slide
└── Output/                            # All results written here
```

---

## Quick-start

### Prerequisites

| Requirement | Notes |
|---|---|
| Apptainer ≥ 1.2 | Used to run the containerised R environment |
| `giotto_rstudio_latest.sif` | Singularity image under `Environment/Docker_Image/` |
| Raw CosMx directories | One per physical slide under `Data/Raw_data/` |
| Populated `sample_sheet.csv` | See [Sample sheet](#sample-sheet) section |

---

## The sample sheet

`Scripts/Parameters/sample_sheet.csv` defines every logical sample the
pipeline will process.  **One row = one independently processed sample.**
For composite slides containing multiple biopsies (e.g. the CART slide with
four tissue sections), add one row per biopsy and use `fov_min`/`fov_max` to
subset the FOV range after loading.

### Columns

| Column | Required | Description |
|---|---|---|
| `sample_id` | ✅ | Unique identifier used for all output files and checkpoint directories |
| `directory_name` | ✅ | Name of the raw data folder under `raw_data_dir`; multiple rows may share the same `directory_name` for composite slides |
| `slide_num` | ✅ | Physical slide number (integer) |
| `slide_id` | ✅ | **Run/batch identifier for Harmony batch correction** — samples prepared and sequenced together share the same `slide_id` (see [Batch correction](#batch-correction)) |
| `patient_id` | — | Patient identifier for paired/repeated-measures analysis |
| `pair_id` | — | Explicitly links T0/T12 biopsies from the same patient for paired DE |
| `group_id` | — | Treatment group (`CART`, `Control`, `Conventional`) — used for `--groups` filtering |
| `treatment` | — | Treatment label passed through to the merged spatial DE step |
| `timepoint` | — | Biopsy timepoint (`T0`, `T12`) — informational |
| `include` | — | `TRUE`/`FALSE` — set to `FALSE` to exclude a sample without deleting the row |
| `fov_min` | — | Minimum FOV number to retain (inclusive). Leave blank to keep all FOVs |
| `fov_max` | — | Maximum FOV number to retain (inclusive). Leave blank to keep all FOVs |
| `notes` | — | Free-text notes |

### Current sample manifest

| sample_id | directory_name | slide_num | slide_id | group_id | timepoint | fov_min | fov_max |
|---|---|---|---|---|---|---|---|
| CART_T0_S1 | MCFASCIST1250Slide1… | 1 | batch_CART | CART | T0 | 55 | 71 |
| CART_T12_S1 | MCFASCIST1250Slide1… | 1 | batch_CART | CART | T12 | 42 | 54 |
| CART_T0_S2 | MCFASCIST1250Slide1… | 1 | batch_CART | CART | T0 | 29 | 41 |
| CART_T12_S2 | MCFASCIST1250Slide1… | 1 | batch_CART | CART | T12 | 1 | 28 |
| 1281 | MCFASCIST1250Slide2… | 2 | batch_Control | Control | T0 | — | — |
| 19H28111 | MCFASCLST1196Slide31… | 31 | batch_Conv_2 | Conventional | T0 | — | — |
| 17H13349 | MCFASCLST1196_Slide3… | 3 | batch_Conv_1 | Conventional | T0 | — | — |
| 20H24159 | MCFASCLST1196Slide32… | 32 | batch_Conv_2 | Conventional | T12 | — | — |
| 18H12037 | MCFASCLST1196_Slide4… | 4 | batch_Conv_1 | Conventional | T12 | — | — |

> **CART pairing**: `CART_T0_S1` ↔ `CART_T12_S1` (pair_id `pair_CART_T12_S1`) and
> `CART_T0_S2` ↔ `CART_T12_S2` (pair_id `pair_CART_T12_S2`).  Both T0/T12 pairs
> derive from the same physical slide (`MCFASCIST1250Slide11979876960634322`);
> they are separated at load time by their FOV ranges.

---

## Configuration (`config.yaml`)

All pipeline behaviour is controlled by `Scripts/Parameters/config.yaml`.

### `paths`

```yaml
paths:
  project_dir: "../.."          # Root of the project (relative to config file)
  scripts_dir: ".."             # Directory containing the *.R scripts
  raw_data_dir: "../../Data/Raw_data"
  output_dir: "../../Output"
  sample_sheet: "sample_sheet.csv"
  python_path: null             # Override Python path; null = auto-detect
```

### `pipeline`

```yaml
pipeline:
  mode: "all"                   # "all" | "separate" | "merged" | "paired"
  save_intermediates: true      # Save a Giotto checkpoint after each step
  overwrite_existing: false     # Re-run steps even if a checkpoint exists
  skip_on_error: false          # Continue to next sample on failure
  memory_cleanup: true          # Call gc() between steps
```

### `merged`

Controls the multi-sample merge and Harmony batch correction step.

```yaml
merged:
  input_step: "07_annotate"     # Which per-sample checkpoint to load as merge input
  join_method: "shift"          # Coordinate shift method for joinGiottoObjects
  x_padding: 1000               # µm padding between samples in merged space
  batch_column: "slide_id"      # Column used by Harmony (see Batch correction)
  dimensions_to_use: "1:30"
  umap_n_neighbors: 30
  umap_min_dist: 0.3
  k_nn: 15
  resolution: 0.3
  run_label: null               # Custom label for the merged run; null = timestamp
```

### `fov_split`

Enables automatic FOV subsetting at load time for composite slides.

```yaml
fov_split:
  enabled: true
  fov_column: "fov"             # Metadata column containing the FOV number
  fov_min_col: "fov_min"        # sample_sheet column for lower bound
  fov_max_col: "fov_max"        # sample_sheet column for upper bound
```

### `spatial_de`

Controls step 12 (spatial differential expression).

```yaml
spatial_de:
  sample_backend: "smiDE"       # "smiDE" | "edgeR"
  annotation_column: "celltype_HCA_Kidney_supervised"
  smide_radius: 50              # Neighbourhood radius in µm (CosMx coords are µm)
  smide_overlap_threshold: 0.85 # Max overlap ratio to retain a gene (0–1)
  smide_min_detection_fraction: 0.01  # Min fraction of cells expressing a gene
  smide_annotation_subset: ["B.cell"]  # Restrict to these cell types; null = all
  smide_ncores: 4
  n_niches: 6
  min_cells_per_niche: 10
```

> **`smide_radius` note**: CosMx spatial coordinates are in microns.  The
> correct value for cell–cell proximity modelling is **50** (50 µm ≈ one cell
> diameter).  The default `0.05` is a normalised-coordinate value and will
> produce `No genes passed smiDE filters` for CosMx data.

### `parameters`

Per-step analytical parameters (QC thresholds, clustering resolution, etc.).
See the full `config.yaml` for all keys and inline comments.

---

## Batch correction

Harmony is run during the merged pipeline step using a single metadata column
as the batch variable (`merged.batch_column`).  The recommended value is
`"slide_id"`, which encodes the **physical preparation and sequencing run**:

| `slide_id` | Samples | Rationale |
|---|---|---|
| `batch_CART` | CART_T0_S1, CART_T12_S1, CART_T0_S2, CART_T12_S2 | All four biopsies on the same physical slide; same library prep and run |
| `batch_Control` | 1281 | Standalone slide; separate run |
| `batch_Conv_1` | 17H13349, 18H12037 | Slides 3 + 4 prepared and run together |
| `batch_Conv_2` | 19H28111, 20H24159 | Slides 31 + 32 prepared and run together |

**Why not `group_id`?**  Harmony would then correct for CART vs Control vs
Conventional — collapsing the very biological signal you want to study.

**Why not `slide_num`?**  Slides 3 and 4 (and 31 and 32) were run together but
have different slide numbers; treating them as independent batches
over-corrects.

---

## Pipeline steps

| Step ID | Alias(es) | Script | Description |
|---|---|---|---|
| `01_load` | `01_load_data` | `01_Load_data.R` | Load CosMx CSVs → Giotto object; subset by FOV if `fov_min`/`fov_max` set |
| `02_qc` | `02_quality_control` | `02_Quality_Control.R` | Filter low-quality cells and genes |
| `03_norm` | `03_normalisation`, `03_normalization` | `03_Normalisation.R` | Library-size normalisation and log-transform |
| `04_dimred` | `04_dimensionality_reduction` | `04_Dimensionallity_Reduction.R` | HVG selection, PCA, UMAP |
| `05_cluster` | — | `05_Clustering.R` | k-NN graph + Leiden clustering |
| `06_markers` | `06_de`, `06_differential_expression` | `06_Differential_Expression.R` | Per-cluster marker genes |
| `07_annotate` | `07_annotation` | `07_Annotation.R` | Reference-profile cell-type annotation |
| `08_visualize` | `08_visualisation` | `08_Visualisation.R` | Feature plots, UMAP panels |
| `09_spatial` | `09_spatial_network` | `09_Spatial_Network.R` | Delaunay spatial network + cell-proximity enrichment |
| `10_cci` | `10_cci_analysis` | `10_CCI_Analysis.R` | InSituCor, LIANA, MISTy CCI analysis |
| `11_bcell` | `11_b_cell_analysis` | `11_B_Cell_Analysis.R` | Focused immune-cell microenvironment |
| `12_spatial_de` | `12_spatial_differential_expression` | `12_Spatial_Differential_Expression.R` | Spatial DE via smiDE or edgeR pseudobulk |
| **Merged** | | | |
| `merge` | — | `Merge_Batch_Correction.R` | Join all per-sample Giotto objects |
| `merge_batch` | `11_batch_correction` | `Merge_Batch_Correction.R` | Harmony batch correction + merged UMAP/clustering |
| `12_spatial_de` *(merged)* | — | `12_Spatial_Differential_Expression.R` | Cohort-level spatial DE (edgeR pseudobulk) |

---

## Running the pipeline

All examples use the Apptainer wrapper:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh [OPTIONS]
```

### CLI options

| Option | Value | Description |
|---|---|---|
| `--config` | path | Override the config file path |
| `--mode` | `all` \| `separate` \| `merged` \| `paired` | Pipeline mode (default: value in `config.yaml`) |
| `--samples` | `id1,id2,...` | Comma-separated `sample_id` values to process |
| `--groups` | `CART,Control,...` | Filter by `group_id` column in sample sheet |
| `--pairs` | `pair_id,...` | Filter by `pair_id` column in sample sheet |
| `--sample-steps` | `01_load,02_qc,...` | Which per-sample steps to run |
| `--merged-steps` | `merge,merge_batch,...` | Which merged steps to run |
| `--from-step` | `05_cluster` | Start per-sample pipeline from this step (resume) |
| `--merged-from-step` | `merge_batch` | Start merged pipeline from this step (resume) |
| `--run-label` | `my_run` | Custom label for the merged output directory |
| `--dry-run` | *(flag)* | Print what would run without executing |

---

## Use-case examples

### 1 — Run everything for all samples automatically

Processes every `include: TRUE` sample through all 12 steps, then runs the
full merged pipeline (merge → batch correction → merged spatial DE).

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh
```

Or explicitly:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh --mode all
```

---

### 2 — Run one sample through all steps

Process sample `1281` (Control) from loading through spatial DE.

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --samples 1281 \
  --mode separate
```

---

### 3 — Run one sample through selected steps only

Resume `1281` from clustering onwards (steps 05–12), skipping already-complete
loading, QC, normalisation, and dimensionality reduction:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --samples 1281 \
  --from-step 05_cluster \
  --mode separate
```

Or run only specific steps explicitly:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --samples 1281 \
  --sample-steps 07_annotate,08_visualize,09_spatial \
  --mode separate
```

---

### 4 — Rerun spatial DE for all CART sub-samples

The CART slide is split into four logical samples.  To rerun step 12 for all
four simultaneously:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --groups CART \
  --sample-steps 12_spatial_de \
  --mode separate
```

Or by listing individual sample IDs:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --samples CART_T0_S1,CART_T12_S1,CART_T0_S2,CART_T12_S2 \
  --sample-steps 12_spatial_de \
  --mode separate
```

---

### 5 — Run the merged pipeline only (all merge steps)

Assumes per-sample step 07 checkpoints already exist.

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh --mode merged
```

With a custom run label:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --mode merged \
  --run-label cohort_v2
```

---

### 6 — Rerun batch correction and merged spatial DE only (skip re-merging)

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --mode merged \
  --merged-from-step merge_batch
```

---

### 7 — Paired CART analysis (T0 vs T12)

The CART biopsies are linked by `pair_id`:

- `pair_CART_T0` → `CART_T0_S1` + `CART_T0_S2`  
- `pair_CART_T12` → `CART_T12_S1` + `CART_T12_S2`

Within-slide pairings for longitudinal comparison:

- **Pair A**: `CART_T0_S1` ↔ `CART_T12_S1` (FOV 55–71 and 42–54)
- **Pair B**: `CART_T0_S2` ↔ `CART_T12_S2` (FOV 29–41 and 1–28)

To process only the CART paired samples through the full per-sample pipeline:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --groups CART \
  --mode separate
```

Then run a targeted merged analysis for the CART group only with a dedicated
run label:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --groups CART \
  --mode merged \
  --run-label CART_paired
```

The merged spatial DE step (`12_spatial_de` in merged mode) will use the
`treatment` and `patient_id` columns from the sample sheet to test T0 vs T12
within each CART biopsy pair using edgeR pseudobulk.

---

### 8 — Preview any run without executing (`--dry-run`)

Prints the resolved config path, selected samples, and step order without
running anything:

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --groups CART \
  --sample-steps 12_spatial_de \
  --mode separate \
  --dry-run
```

---

## Output structure

```
Output/
├── Sample_CART_T0_S1/
│   ├── 01_Data_Loading/
│   ├── 02_Quality_Control/
│   ├── ...
│   ├── 12_Spatial_Differential_Expression/
│   │   ├── tables/
│   │   └── plots/
│   ├── Giotto_Object_Annotated/     # Step 07 checkpoint
│   └── pipeline_checkpoints/        # Per-step checkpoints
├── Sample_CART_T12_S1/
├── Sample_1281/
├── Merged/
│   └── <run_label>/
│       ├── 10_Merged/
│       │   └── Batch_Correction/
│       ├── 12_Spatial_Differential_Expression/
│       ├── Giotto_Object_Merged/
│       └── Giotto_Object_BatchCorrected/
└── pipeline_log_<timestamp>.log
```

---

## Troubleshooting

### `No samples matched the selected filters`

The `--samples` argument must match `sample_id` values in `sample_sheet.csv`
exactly.  The old composite ID `1979_8769_6063_4320` no longer exists — it has
been replaced by `CART_T0_S1`, `CART_T12_S1`, `CART_T0_S2`, `CART_T12_S2`.
Use `--groups CART` to select all four at once.

### `incomplete final line found on config.yaml`

The `config.yaml` file is missing a trailing newline.  Open it and add a blank
line at the very end, or run:

```bash
echo "" >> Scripts/Parameters/config.yaml
```

### `Batch column 'slide_num' is not present in merged metadata`

The merged Giotto object does not carry a `slide_num` column in its cell
metadata.  Ensure `batch_column` in `config.yaml` is set to `"slide_id"` and
that the `slide_id` column is populated in `sample_sheet.csv`.

### `No genes passed the smiDE filters`

Two common causes:

1. **Wrong radius unit** — `smide_radius` must be `50` (µm) for CosMx data,
   not `0.05`.
2. **Filters too strict for small populations** — for B cells (~100–200 cells),
   set `smide_overlap_threshold: 0.85` and `smide_min_detection_fraction: 0.01`.

### `renv` out-of-sync warning

```
The project is out-of-sync -- use `renv::status()` for details.
```

This is a non-fatal warning from the container's `renv` lockfile.  The pipeline
will still run.  To resolve, enter the container interactively and run
`renv::restore()`.

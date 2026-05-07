# CosMx / Giotto Analysis Pipeline

Spatially resolved transcriptomic analysis of NanoString CosMx data using the
[Giotto](https://giottosuite.readthedocs.io) suite. The pipeline runs raw
CosMx exports through QC, normalisation, dimensionality reduction,
clustering, annotation, cell-cell interaction analysis, spatial differential
expression, and per-sample pathway enrichment, orchestrated by a single R
driver (`CosMx_pipeline.R`) launched through an Apptainer wrapper
(`Run_Scripts/Run_Giotto_Pipeline.sh`).

This README documents the single-sample workflow only. Cross-sample merged
and paired analyses are out of scope until the merged side is fixed
(ETA: Eventually).

---

## Repository layout

```
CosMx_analysis/
  Scripts/
    CosMx_pipeline.R                       - Pipeline driver / orchestrator
    00_Setup.R                             - Library + Python env bootstrap
    01_Load_data.R                         - CosMx CSVs to Giotto object
    02_Quality_Control.R                   - Per-cell / per-gene QC
    03_Normalisation.R                     - Library-size normalisation
    04_Dimensionallity_Reduction.R         - HVGs, PCA, UMAP
    05_Clustering.R                        - k-NN + Leiden
    06_Differential_Expression.R           - Per-cluster marker genes (scran)
    07_Annotation.R                        - InSituType reference annotation
    08_Visualisation.R                     - Cluster + feature plots
    09_Spatial_Network.R                   - Delaunay network + cell-proximity
    10_CCI_Analysis.R                      - InSituCor / LIANA / NicheNet / MISTy / nnSVG
    11_B_Cell_Analysis.R                   - Focused immune microenvironment
    12_Spatial_Differential_Expression.R   - smiDE niche DE (sample scope)
    13_Pathway_Analysis.R                  - MSigDB / Reactome / GSEA (sample scope)
    Helper_Scripts/
      Pipeline_Utils.R                     - Checkpointing, plotting, safe_run
      Plot_Helpers.R                       - presentation_theme, save_presentation_plot
      Auto_Detect_Samples.R                - Build sample_sheet.csv from raw_data_dir
      Split_Raw_CosMx_Slide.R              - Per-FOV physical slide splitter
      ...                                  - Smaller utilities
    Parameters/
      config.yaml                          - All pipeline parameters
      sample_sheet.csv                     - Sample manifest
      Install_CCI_dependencies.R           - On-demand R-pkg installer
    Run_Scripts/
      Run_Giotto_Pipeline.sh               - Apptainer launch wrapper
  Data/Raw_data/                           - One subdirectory per physical slide
  Environment/Docker_Image/                - giotto_rstudio_*.sif
  Output/                                  - All results land here
```

---

## Quick start

### Prerequisites

| Requirement | Notes |
|---|---|
| Apptainer >= 1.2 | Containerised R environment |
| `giotto_rstudio_latest.sif` | Place under `Environment/Docker_Image/` |
| Conda env with `umap`, `igraph`, `leidenalg` | E.g. `giotto_Py_3_11`. Set `COSMX_PYTHON_PATH` to its `python3` |
| Raw CosMx directories | One per physical slide under `Data/Raw_data/` |
| Populated `sample_sheet.csv` | See [Sample sheet](#sample-sheet) |

### Run a single sample end-to-end

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh --samples 1281
```

Each per-sample run lives in its own R / Apptainer process and writes to
`Output/Sample_<id>/`.

---

## Sample sheet

`Scripts/Parameters/sample_sheet.csv` defines every logical sample.

### Columns

| Column | Required | Description |
|---|---|---|
| `sample_id` | yes | Physical slide / directory identifier |
| `subsample_id` | no | Logical sub-sample (e.g. `CART_T0_S1`); used by `--split` |
| `directory_name` | yes | Folder name under `paths.raw_data_dir` |
| `slide_num` | yes | Physical slide number |
| `slide_id` | yes | Slide / batch identifier |
| `split_role` | no | `composite`, `split`, or empty - which view a row belongs to (see below) |
| `patient_id` | no | Patient identifier |
| `pair_id` | no | Links T0 / T12 biopsies from the same patient |
| `group_id` | no | Treatment group (`CART`, `Control`, `Conventional`) |
| `treatment` | no | Treatment label |
| `timepoint` | no | `T0`, `T12` - informational |
| `include` | no | `FALSE` to exclude without deleting the row |
| `fov_min`, `fov_max` | no | Inclusive FOV bounds, applied at load time when `--split` is on |
| `notes` | no | Free-text |

### Split mode (default) vs composite mode

The launcher defaults to `--split`. The pipeline picks rows based on
`split_role`:

- **`--split` (default)** - keeps rows with `split_role` in {`split`, empty}.
  Composite slides are expanded into per-biopsy runs, each restricted to
  its `fov_min`/`fov_max` range and written to `Sample_<subsample_id>/`.
- **`--no-split`** - keeps rows with `split_role` in {`composite`, empty}.
  The CART slide is loaded as one whole-slide sample written to
  `Sample_<sample_id>/`. Useful for composite-view UMAPs only.

Sample selection rules:

- `include = FALSE` is the default off switch.
- `--samples X` always runs `X`, bypassing both `include` and `split_role`.
- `--groups CART,...` filters by `group_id`.

### Current manifest

| sample_id | subsample_id | slide_num | slide_id | group_id | timepoint | fov_min | fov_max |
|---|---|---|---|---|---|---|---|
| 1979_8769_6063_4320 | CART_T0_S1 | 1 | batch_CART | CART | T0 | 55 | 71 |
| 1979_8769_6063_4320 | CART_T12_S1 | 1 | batch_CART | CART | T12 | 42 | 54 |
| 1979_8769_6063_4320 | CART_T0_S2 | 1 | batch_CART | CART | T0 | 29 | 41 |
| 1979_8769_6063_4320 | CART_T12_S2 | 1 | batch_CART | CART | T12 | 1 | 28 |
| 1281 | - | 2 | batch_Control | Control | T0 | - | - |
| 17H13349 | - | 3 | batch_Conv_1 | Conventional | T0 | - | - |
| 18H12037 | - | 4 | batch_Conv_1 | Conventional | T12 | - | - |
| 19H28111 | - | 31 | batch_Conv_2 | Conventional | T0 | - | - |
| 20H24159 | - | 32 | batch_Conv_2 | Conventional | T12 | - | - |

---

## Pipeline steps (per-sample, run in order)

| Step ID | Aliases | Script | Description |
|---|---|---|---|
| `01_load` | `01_load_data` | `01_Load_data.R` | Load CosMx CSVs to Giotto; FOV subset if `--split` |
| `02_qc` | `02_quality_control` | `02_Quality_Control.R` | Filter low-quality cells + genes |
| `03_norm` | `03_normalisation` | `03_Normalisation.R` | Library-size normalise + log-transform |
| `04_dimred` | `04_dimensionality_reduction` | `04_Dimensionallity_Reduction.R` | HVGs, PCA, UMAP |
| `05_cluster` | - | `05_Clustering.R` | k-NN + Leiden |
| `06_markers` | `06_de` | `06_Differential_Expression.R` | Per-cluster marker genes |
| `07_annotate` | `07_annotation` | `07_Annotation.R` | InSituType reference annotation |
| `08_visualize` | `08_visualisation` | `08_Visualisation.R` | Feature + cluster plots |
| `09_spatial` | `09_spatial_network` | `09_Spatial_Network.R` | Delaunay network + cell-proximity enrichment |
| `10_cci` | `10_cci_analysis` | `10_CCI_Analysis.R` | InSituCor / LIANA / NicheNet / MISTy / nnSVG |
| `11_bcell` | `11_b_cell_analysis` | `11_B_Cell_Analysis.R` | Focused immune microenvironment |
| `12_spatial_de` | `12_spatial_differential_expression` | `12_Spatial_Differential_Expression.R` | smiDE niche DE (sample scope) |
| `13_pathway` | `13_pathway_analysis` | `13_Pathway_Analysis.R` | Per-sample GSEA / ORA on cluster + celltype DE |

---

## CLI options

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh [OPTIONS]
```

Default mode is per-sample (`--mode separate`). Cross-sample modes
(`--mode merged`, `--mode paired`, `--pairs`) exist in the codebase but are
not documented here until the merged side is fixed (ETA: Eventually).

| Option | Value | Description |
|---|---|---|
| `--samples` | `id1,id2,...` | Run specific `sample_id` / `subsample_id` rows |
| `--groups` | `CART,Control,...` | Filter by `group_id` |
| `--sample-steps` | `01_load,02_qc,...` | Restrict per-sample steps |
| `--from-step` | `05_cluster` | Resume per-sample pipeline from this step |
| `--split` | flag | Expand composite slides into per-subsample runs (**default**) |
| `--no-split` | flag | Composite-view: keep CART rows as one whole-slide sample |
| `--config` | path | Override config file path |
| `--dry-run` | flag | Print resolved plan without executing |

---

## Configuration (`Parameters/config.yaml`)

Highlights only - see inline comments in `config.yaml` for every key.

### `paths`

```yaml
project_dir: "../.."
scripts_dir: ".."
raw_data_dir: "../../Data/Raw_data"
output_dir: "../../Output"
sample_sheet: "sample_sheet.csv"
python_path: null              # null = use COSMX_PYTHON_PATH env var
```

### `pipeline`

```yaml
save_intermediates: true       # qs::qsave a Giotto checkpoint after each step
overwrite_existing: false      # currently advisory only
skip_on_error: false           # continue to next sample if one fails
memory_cleanup: true           # gc() between steps
```

### `cci` - Step 10 sub-section toggles + nnSVG raster

```yaml
run_sections:
  insitucor: true
  liana: true
  nichenet: true
  misty: true
  nnsvg: true
nnsvg_raster_resolution: auto  # SEraster pixel side length; auto = pick from extent
```

nnSVG flow: tissue-wide pass always tries `SEraster::rasterizeGeneExpression`;
if it fails (auto-pick failed or bins sub-cellular), falls back to a 5k
random cell subsample so the pass always produces results. The B-cell-focused
pass runs per-cell unconditionally.

### `spatial_de` - Step 12 (smiDE)

```yaml
sample_backend: "smiDE"        # "smiDE" | "edgeR"
annotation_column: "celltype_HCA_Kidney_supervised"
smide_radius: 50               # micrometres (CosMx coords are um)
smide_overlap_threshold: 0.8   # 1.0 = disabled
smide_min_detection_fraction: 0.01
smide_predictor_mode: "group"  # "group" | "count" (WTX-guide continuous)
smide_ncores: 4                # only consumed by smi_de(); pre_de/results are single-threaded
n_niches: 6
min_cells_per_niche: 10
```

### `pathway` - Step 13 (per-sample scope)

```yaml
enabled: true
de_engine: "pseudobulk_deseq2"
species: "Homo sapiens"
databases: [hallmark, kegg, reactome, biocarta, go_bp]
gsea: { ... }
```

The `comparisons:` block in `pathway` drives cross-sample contrasts and is
only consumed by the merged scope (out of scope here).

### Other sections

`reproducibility` (RNG seed), `fov_split` (load-time FOV subsetting),
`interaction` (Step 09/11 focal-cell config), `parameters` (per-step QC +
clustering thresholds).

---

## Use-case examples

### Run one sample end-to-end

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh --samples 1281
```

### Resume one sample from a specific step

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --samples 1281 --from-step 07_annotate
```

### Re-run only spatial DE for the CART biopsies

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --groups CART --sample-steps 12_spatial_de
```

The default `--split` is already in effect, so each biopsy runs
independently against its FOV range.

### Composite view (CART as one whole slide)

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh --no-split
```

### Dry-run preview

```bash
./Scripts/Run_Scripts/Run_Giotto_Pipeline.sh \
  --groups CART --sample-steps 12_spatial_de --dry-run
```

---

## Output structure

```
Output/
  Sample_CART_T0_S1/                                  - one dir per (sub)sample
    01_Data_Loading/
    02_Quality_Control/
    ...
    10_CCI_Analysis/
      insitucor/
      liana/
      nichenet/
      misty/
      svg/                                            - nnSVG outputs
    12_Spatial_Differential_Expression/
      tables/
      plots/
    13_Pathway/
    pipeline_checkpoints/                             - per-step Giotto checkpoints
      <step_id>/manifest.json
    Pipeline_log_<sample>_<timestamp>.log
  Sample_CART_T12_S1/
  Sample_1281/
  ...
  pipeline_log_<timestamp>.log                        - global log (tee from launcher)
```

Per-step checkpoints under `pipeline_checkpoints/<step_id>/manifest.json`
let you resume mid-pipeline (`--from-step`).

---

## Troubleshooting

**`No samples matched the selected filters`**

`--samples` matches `sample_id` or `subsample_id` in `sample_sheet.csv`.
For all CART biopsies use `--groups CART`.

**`No genes passed the smiDE filters`**

- `spatial_de.smide_radius` must be `50` (micrometres) for CosMx, not `0.05`.
- For rare cell types (e.g. B cells, ~100-200 cells), set
  `smide_overlap_threshold: 0.8` and `smide_min_detection_fraction: 0.01`.

**nnSVG: `Spots (cells): N` instead of `Spots (pixels): N`**

The tissue-wide pass fell back to a 5k cell subsample because SEraster
rasterisation could not produce usable bins (auto-pick failed or pixels
were sub-cellular). The pass still completes; P14/P15 plots (which require
rasterised celltype counts) are skipped. To force rasterisation, set
`cci.nnsvg_raster_resolution` to a numeric value.

**Step crashes with `Execution halted` and no error message**

The orchestrator sinks `message`/`stop` output to the per-sample log only
(`split = FALSE` at `CosMx_pipeline.R:1046`). Check
`Output/Sample_<id>/Pipeline_log_<sample>_<timestamp>.log` for the actual
stack trace.


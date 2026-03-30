# Parameters and Support Files

This folder holds configuration and one-time environment setup files for the
CosMx Giotto pipeline.

## Files in this folder

- `config.yaml`
  - main pipeline configuration
  - paths are resolved relative to this folder
- `sample_sheet_template.csv`
  - example sample-sheet structure
  - copy and adapt this into your real server-side sample sheet if needed
- `Install_CCI_Dependencies.R`
  - one-time installer for optional step 10 dependencies
  - also includes optional spatial DE extras such as `smiDE`
  - do not run this as part of the analysis pipeline itself
- `RUNTIME_SETUP.md`
  - notes about the expected R/Python/container runtime

## Expected repo layout

- top level: pipeline scripts (`00_Setup.R` to `12_Spatial_Differential_Expression.R`) and `CosMx_pipeline.R`
- helper scripts: reusable utilities under `Helper_Scripts/`, including `Merge_Batch_Correction.R`
- `Helper_Scripts/`: shared helper utilities
- `Shiny_Portal/`: Shiny app files
- `Parameters/`: configuration and dependency setup files
- `test_scripts/`: manual test and launch scripts

## Typical manual workflow

1. Edit `config.yaml` for the server paths you actually use.
2. Install optional CCI dependencies once if you plan to run step 10.
3. If you installed `smiDE`, run `test_scripts/Inspect_smiDE.R` once to inspect the installed API.
4. Open `test_scripts/Manual_Single_Sample_Test.R` in RStudio for one-sample testing.

By default, `spatial_de.sample_backend` is now set to `smiDE` for single-sample
spatially resolved DE. Merged-object spatial DE still uses `edgeR`.

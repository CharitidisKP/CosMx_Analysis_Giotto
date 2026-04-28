#!/bin/bash
# ==============================================================================
# CosMx Pipeline Launcher
# Usage:
#   ./Run_Giotto_Pipeline.sh                   # DEFAULT split view: 4 CART biopsies + Control + 4 Conv (9 samples)
#   ./Run_Giotto_Pipeline.sh --no-split        # composite view: CART_composite + Control + 4 Conv (6 samples)
#   ./Run_Giotto_Pipeline.sh --samples 1281    # one sample (bypasses include=FALSE if set)
#   ./Run_Giotto_Pipeline.sh --samples CART_T0_S1
#                                              # runs a single CART split (default already selects splits)
#   ./Run_Giotto_Pipeline.sh --samples 1281 --sample-steps 01_load,02_qc
#   ./Run_Giotto_Pipeline.sh --pairs pair_CART_S1
#                                              # T0 + T12 for CART patient 1 (paired contrast)
#   ./Run_Giotto_Pipeline.sh --dry-run         # preview without running
#   ./Run_Giotto_Pipeline.sh --overwrite       # regenerate existing checkpoints
#                                              #  + section outputs (default skips them)
#
# Skip-by-default: per-step Giotto checkpoints (steps 01-09, 11) and
# canonical CCI / smiDE section outputs (CSVs) are detected on disk and
# skipped on re-run. Pass --overwrite to force regeneration.
#
# Server pre-flight (run before launching with multiple workers later):
#   nproc; free -g; uptime
# Pick a worker count yourself based on server load.
#
# Default set to --split 2026-04-27: all of the cross-sample comparisons in
# config.yaml::pathway.comparisons (CART/Conv before-vs-after, CART vs Conv
# at each timepoint, CART vs Control at each timepoint) require T0/T12 +
# patient labels at the cell level, which only exist when the CART composite
# is split into its 4 biopsies. --no-split is an alternative escape hatch
# (CART_composite + Control + 4 Conv = 6 samples) for the rare case you want
# a composite-view overview UMAP — not the regular run.
#
# Sample selection rules:
#   - `include` column in sample_sheet.csv is the default on/off switch.
#   - `split_role` column ("composite"/"split"/"") controls which view is
#     active: --split (default) keeps "split" and ""; --no-split keeps
#     "composite" and "".
#   - --samples X always runs X, bypassing include and split_role.
#   - --pairs / --groups layer on top but do NOT resurrect include=FALSE rows.
#
# All arguments are passed directly to CosMx_pipeline.R.
# Output is logged to Output/pipeline_log_<timestamp>.log.
#
# Environment overrides (set before invocation to customise):
#   PROJECT_DIR     Absolute path to the CosMx project root.
#                   Default: $HOME/P_lab/CosMx_analysis
#   SIF             Path to the Apptainer .sif image.
#                   Default: $PROJECT_DIR/Environment/Docker_Image/giotto_rstudio_latest.sif
#   CONFIG          Path to config.yaml.
#                   Default: $PROJECT_DIR/Scripts/Parameters/config.yaml
#   R_LIBS_USER     R library path inside the container. Leave UNSET to let
#                   renv activate from $PROJECT_DIR/.Rprofile (the standard
#                   path — packages live in $PROJECT_DIR/renv/library/).
#                   Only set this if you have a non-renv setup.
#   COSMX_PYTHON_PATH
#                   Python interpreter for reticulate/InSituType.
#                   Default: $HOME/anaconda3/envs/giotto_Py_3_11/bin/python3
#   TMPDIR_HOST     Host-side temp directory to bind into /tmp.
#                   Default: $PROJECT_DIR/tmp
# ==============================================================================

set -euo pipefail

# ---- Host paths (override via environment variables) ------------------------
PROJECT_DIR="${PROJECT_DIR:-$HOME/P_lab/CosMx_analysis}"
SIF="${SIF:-$PROJECT_DIR/Environment/Docker_Image/giotto_rstudio_latest.sif}"
CONFIG="${CONFIG:-$PROJECT_DIR/Scripts/Parameters/config.yaml}"
TMPDIR_HOST="${TMPDIR_HOST:-$PROJECT_DIR/tmp}"
LOG_DIR="$PROJECT_DIR/Output"
# renv (in $PROJECT_DIR/renv) handles the library path via .Rprofile when
# Rscript is started from $PROJECT_DIR. Only honour an externally-set
# R_LIBS_USER (advanced override); do NOT default to one — that bypasses renv.
R_LIBS_USER="${R_LIBS_USER:-}"
COSMX_PYTHON_PATH="${COSMX_PYTHON_PATH:-$HOME/anaconda3/envs/giotto_Py_3_11/bin/python3}"

# ---- Pre-flight checks ------------------------------------------------------
if [[ ! -f "$SIF" ]]; then
  echo "ERROR: Apptainer image not found at: $SIF" >&2
  echo "  Set SIF=/path/to/image.sif or place it at the default location." >&2
  exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: Pipeline config not found at: $CONFIG" >&2
  echo "  Set CONFIG=/path/to/config.yaml or place it at the default location." >&2
  exit 1
fi

# ---- Python interpreter sanity check ----------------------------------------
# COSMX_PYTHON_PATH often gets overridden by a stale export in ~/.bashrc or
# ~/.profile (most common: /usr/bin/python3 which lacks umap/igraph/leidenalg).
# Fail loudly before the pipeline runs for an hour and then can't do UMAP.
if [[ ! -x "$COSMX_PYTHON_PATH" ]]; then
  echo "ERROR: COSMX_PYTHON_PATH does not point to an executable: $COSMX_PYTHON_PATH" >&2
  echo "  Expected: conda env with umap, igraph, leidenalg (e.g. giotto_Py_3_11)." >&2
  echo "  Fix: unset COSMX_PYTHON_PATH, or export the correct path before launch." >&2
  exit 1
fi
if [[ "$COSMX_PYTHON_PATH" == "/usr/bin/python3" ]]; then
  echo "WARNING: COSMX_PYTHON_PATH=/usr/bin/python3 — this is the system python," >&2
  echo "  which usually lacks umap/igraph/leidenalg. Giotto dim-reduction and" >&2
  echo "  clustering will fail. Check 'env | grep -i python' and your shell rc." >&2
fi

mkdir -p "$TMPDIR_HOST" "$LOG_DIR"
[[ -n "$R_LIBS_USER" ]] && mkdir -p "$R_LIBS_USER"

# Default to --split unless the user explicitly opted in/out. Padding with
# spaces around "$*" avoids matching --splitfoo or --no-splittery.
case " $* " in
  *" --split "*|*" --no-split "*) ;;
  *) set -- --split "$@" ;;
esac

# ---- Resolve current user (fall back to id -un inside containerless envs) ---
USER_NAME="${USER:-$(id -un)}"
LOG_NAME="${LOGNAME:-$USER_NAME}"

# ---- Build log filename from date + --samples argument ----------------------
LOG_TAG="$(date +%Y-%m-%d_%H%M%S)"
prev=""
for i in "$@"; do
  if [[ "$prev" == "--samples" ]]; then
    # Sanitise special characters that could break the filename.
    SAMPLE_TAG="$(echo "$i" | sed 's/[^A-Za-z0-9_,-]/_/g')"
    LOG_TAG="${SAMPLE_TAG}_${LOG_TAG}"
    break
  fi
  prev="$i"
done

LOG_FILE="$LOG_DIR/pipeline_log_${LOG_TAG}.log"

echo "================================================"
echo "CosMx Pipeline"
echo "  Project: $PROJECT_DIR"
echo "  Config:  $CONFIG"
echo "  Image:   $(basename "$SIF")"
echo "  User:    $USER_NAME"
echo "  R libs:  ${R_LIBS_USER:-<renv via $PROJECT_DIR/renv>}"
echo "  Python:  $COSMX_PYTHON_PATH"
echo "  TmpDir:  $TMPDIR_HOST"
echo "  Log:     $LOG_FILE"
echo "  Args:    $*"
echo "================================================"
echo ""

# Use APPTAINERENV_* prefixing in addition to --env. Some apptainer builds
# drop --env values under --cleanenv; the APPTAINERENV_ prefix is honored
# unconditionally by every apptainer/singularity version and survives
# --cleanenv. Belt-and-braces.
export APPTAINERENV_USER="$USER_NAME"
export APPTAINERENV_LOGNAME="$LOG_NAME"
export APPTAINERENV_R_LIBS_USER="$R_LIBS_USER"
export APPTAINERENV_R_LIBS=""
export APPTAINERENV_COSMX_PYTHON_PATH="$COSMX_PYTHON_PATH"
export APPTAINERENV_RETICULATE_PYTHON="$COSMX_PYTHON_PATH"
export APPTAINERENV_VROOM_TEMP_PATH="$TMPDIR_HOST"
export APPTAINERENV_TMPDIR="$TMPDIR_HOST"

apptainer exec --cleanenv \
  --pwd "$PROJECT_DIR" \
  --env USER="$USER_NAME" \
  --env LOGNAME="$LOG_NAME" \
  --env R_LIBS_USER="$R_LIBS_USER" \
  --env R_LIBS= \
  --env COSMX_PYTHON_PATH="$COSMX_PYTHON_PATH" \
  --env RETICULATE_PYTHON="$COSMX_PYTHON_PATH" \
  --env VROOM_TEMP_PATH="$TMPDIR_HOST" \
  --env TMPDIR="$TMPDIR_HOST" \
  --bind ~/rs-nss/passwd:/etc/passwd:ro \
  --bind ~/rs-nss/group:/etc/group:ro \
  --bind "$HOME":"$HOME" \
  --bind "$TMPDIR_HOST":/tmp \
  "$SIF" \
  Rscript "$PROJECT_DIR/Scripts/CosMx_pipeline.R" \
    --config "$CONFIG" \
    "$@" \
  2>&1 | tee "$LOG_FILE"

# Clean up temp files after pipeline completes or fails
echo "Cleaning up temporary files..."
rm -rf "$TMPDIR_HOST"/*
echo "✓ Temp files cleared"

echo ""
echo "================================================"
echo "Done. Log saved to: $LOG_FILE"
echo "================================================"

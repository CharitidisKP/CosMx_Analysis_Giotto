#!/bin/bash
# ==============================================================================
# CosMx Pipeline Launcher
# Usage:
#   ./run_pipeline.sh                          # all samples, all steps
#   ./run_pipeline.sh --samples 1281           # one sample
#   ./run_pipeline.sh --samples 1281 --sample-steps 01_load,02_qc
#   ./run_pipeline.sh --dry-run                # preview without running
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
#   R_LIBS_USER     R library path inside the container.
#                   Default: $HOME/Rlibs_giotto
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
R_LIBS_USER="${R_LIBS_USER:-$HOME/Rlibs_giotto}"
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

mkdir -p "$TMPDIR_HOST" "$LOG_DIR" "$R_LIBS_USER"

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
echo "  R libs:  $R_LIBS_USER"
echo "  TmpDir:  $TMPDIR_HOST"
echo "  Log:     $LOG_FILE"
echo "  Args:    $*"
echo "================================================"
echo ""

apptainer exec --cleanenv \
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

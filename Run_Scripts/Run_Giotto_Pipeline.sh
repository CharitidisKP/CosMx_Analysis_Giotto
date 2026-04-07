#!/bin/bash
# ==============================================================================
# CosMx Pipeline Launcher
# Usage:
#   ./run_pipeline.sh                          # all samples, all steps
#   ./run_pipeline.sh --samples 1281           # one sample
#   ./run_pipeline.sh --samples 1281 --sample-steps 01_load,02_qc
#   ./run_pipeline.sh --dry-run                # preview without running
#
# All arguments are passed directly to CosMx_pipeline.R
# Output is logged to Output/pipeline_log_<timestamp>.log
# ==============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/P_lab/CosMx_analysis"
SIF="$PROJECT_DIR/Environment/Docker_Image/giotto_rstudio_latest.sif"
CONFIG="$PROJECT_DIR/Scripts/Parameters/config.yaml"
TMPDIR_HOST="$PROJECT_DIR/tmp"
LOG_DIR="$PROJECT_DIR/Output"

# Ensure tmp and Output dirs exist
mkdir -p "$TMPDIR_HOST" "$LOG_DIR"

# Build log filename from date + any --samples argument
LOG_TAG="$(date +%Y-%m-%d_%H%M%S)"
prev=""
for i in "$@"; do
  if [[ "$prev" == "--samples" ]]; then
    LOG_TAG="${i}_${LOG_TAG}"
    break
  fi
  prev="$i"
done

LOG_FILE="$LOG_DIR/pipeline_log_${LOG_TAG}.log"

echo "================================================"
echo "CosMx Pipeline"
echo "  Config:  $CONFIG"
echo "  Image:   $(basename $SIF)"
echo "  Log:     $LOG_FILE"
echo "  Args:    $@"
echo "================================================"
echo ""

apptainer exec --cleanenv \
  --env USER=koncha \
  --env LOGNAME=koncha \
  --env R_LIBS_USER="$HOME/Rlibs_giotto" \
  --env R_LIBS= \
  --env VROOM_TEMP_PATH="/mnt/home/koncha/P_lab/CosMx_analysis/tmp" \
  --env TMPDIR="/mnt/home/koncha/P_lab/CosMx_analysis/tmp" \
  --bind ~/rs-nss/passwd:/etc/passwd:ro \
  --bind ~/rs-nss/group:/etc/group:ro \
  --bind /mnt/home/koncha:/mnt/home/koncha \
  --bind "$TMPDIR_HOST":/tmp \
  "$SIF" \
  Rscript /mnt/home/koncha/P_lab/CosMx_analysis/Scripts/CosMx_pipeline.R \
    --config /mnt/home/koncha/P_lab/CosMx_analysis/Scripts/Parameters/config.yaml \
    "$@" \
  2>&1 | tee "$LOG_FILE"
  
  2>&1 | tee "$LOG_FILE"

# Clean up temp files after pipeline completes or fails
echo "Cleaning up temporary files..."
rm -rf "$TMPDIR_HOST"/*
echo "✓ Temp files cleared"

echo ""
echo "================================================"
echo "Done. Log saved to: $LOG_FILE"
echo "================================================"
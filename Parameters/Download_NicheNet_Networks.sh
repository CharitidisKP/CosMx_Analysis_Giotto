#!/bin/bash
# ==============================================================================
# Download_NicheNet_Networks.sh
# One-time fetch of NicheNet 2.0 prior networks from Zenodo.
# https://doi.org/10.5281/zenodo.7074291
#
# Usage:
#   bash Parameters/Download_NicheNet_Networks.sh          # default dest
#   NICHENET_DIR=/custom/path bash Parameters/Download_NicheNet_Networks.sh
#
# After download, step 10 auto-finds the files under
#   ${NICHENET_DIR}/ or $HOME/P_lab/CosMx_analysis/Reference/NicheNet_networks/
# (default search path in 10_CCI_Analysis.R::run_nichenet, via
#  .resolve_nichenet_prior_files()).
# ==============================================================================

set -euo pipefail

NICHENET_DIR="${NICHENET_DIR:-$HOME/P_lab/CosMx_analysis/Reference/NicheNet_networks}"
mkdir -p "$NICHENET_DIR"

echo "================================================"
echo "NicheNet 2.0 prior networks"
echo "  Destination: $NICHENET_DIR"
echo "================================================"

# File set that .resolve_nichenet_prior_files() accepts (human-mouse versions
# live under the same zenodo record; we fetch the human-only v2 files).
declare -a FILES=(
  "ligand_target_matrix_nsga2r_final.rds"
  "lr_network_human_21122021.rds"
  "weighted_networks_nsga2r_final.rds"
)

BASE_URL="https://zenodo.org/record/7074291/files"

for f in "${FILES[@]}"; do
  dest="$NICHENET_DIR/$f"
  if [[ -f "$dest" && -s "$dest" ]]; then
    echo "  [skip] $f (exists, $(du -h "$dest" | cut -f1))"
    continue
  fi
  echo "  [get]  $f ..."
  curl --fail --location --progress-bar "$BASE_URL/$f" -o "$dest"
done

echo ""
echo "Downloaded files:"
ls -lh "$NICHENET_DIR"

echo ""
echo "Set this in ~/.zshrc or ~/.bashrc to persist:"
echo "  export COSMX_NICHENET_DIR=\"$NICHENET_DIR\""
echo ""
echo "Or leave it unset — step 10 auto-searches the default path."

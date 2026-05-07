#!/bin/bash
# Run the CosMx pipeline for each in-scope sample in its own Apptainer process.
# Avoids cross-sample memory accumulation that OOM-killed the conventional batch
# on 2026-05-06 (python child grew to 130 GB anon-rss before cgroup killed it).
# Each sample is launched via Run_Giotto_Pipeline.sh --samples <id>, sequentially.
# Continues past per-sample failures and prints a summary at the end.
#
# Usage:
#   ./Run_Samples_Serially.sh                                  # 8 in-scope samples (default)
#   ./Run_Samples_Serially.sh --samples 17H13349,18H12037      # explicit subset
#   ./Run_Samples_Serially.sh --dry-run                        # preview only
#   ./Run_Samples_Serially.sh --tmux                           # detach into a tmux session
#   ./Run_Samples_Serially.sh --tmux my_run --samples 19H28111 # named session, single sample
#   ./Run_Samples_Serially.sh -- --sample-steps 02_qc,03_norm  # forward flags to launcher
#
# Default sample list is derived from Parameters/sample_sheet.csv (treatment in
# CART/Conventional, split_role != composite, sample_id != TEST_5k).
# tmux mode re-execs the script in a detached session and prints the attach command.

set -uo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="${PROJECT_DIR:-$(cd "$SCRIPT_DIR/.." && pwd)}"
LAUNCHER="$SCRIPT_DIR/Run_Giotto_Pipeline.sh"
SHEET="${SHEET:-$PROJECT_DIR/Parameters/sample_sheet.csv}"

if [[ ! -x "$LAUNCHER" ]]; then
  echo "ERROR: launcher not found or not executable: $LAUNCHER" >&2
  exit 2
fi
if [[ ! -f "$SHEET" ]]; then
  echo "ERROR: sample sheet not found: $SHEET" >&2
  exit 2
fi

# ---- Argument parsing -------------------------------------------------------
SAMPLES_OVERRIDE=""
DRY_RUN=0
TMUX_MODE=0
TMUX_SESSION=""
PASSTHROUGH=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --samples)
      SAMPLES_OVERRIDE="${2:-}"
      shift 2
      ;;
    --samples=*)
      SAMPLES_OVERRIDE="${1#--samples=}"
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --tmux)
      TMUX_MODE=1
      shift
      if [[ $# -gt 0 && "$1" != --* ]]; then
        TMUX_SESSION="$1"
        shift
      fi
      ;;
    --)
      shift
      PASSTHROUGH+=("$@")
      break
      ;;
    *)
      PASSTHROUGH+=("$1")
      shift
      ;;
  esac
done

# ---- tmux re-exec -----------------------------------------------------------
# When --tmux is requested and we're not already inside a tmux session, spawn a
# detached session that re-runs this script without --tmux. User attaches with
# the command we print. Detach with Ctrl-b d, reattach later.
if [[ "$TMUX_MODE" -eq 1 && -z "${TMUX:-}" ]]; then
  if ! command -v tmux >/dev/null 2>&1; then
    echo "ERROR: tmux not installed; cannot honour --tmux." >&2
    exit 2
  fi
  if [[ -z "$TMUX_SESSION" ]]; then
    TMUX_SESSION="cosmx_serial_$(date +%Y%m%d_%H%M%S)"
  fi
  CMD=("$0")
  if [[ -n "$SAMPLES_OVERRIDE" ]]; then
    CMD+=(--samples "$SAMPLES_OVERRIDE")
  fi
  if [[ "$DRY_RUN" -eq 1 ]]; then
    CMD+=(--dry-run)
  fi
  if [[ ${#PASSTHROUGH[@]} -gt 0 ]]; then
    CMD+=(-- "${PASSTHROUGH[@]}")
  fi
  tmux new-session -d -s "$TMUX_SESSION" "${CMD[@]}"
  echo "Started tmux session: $TMUX_SESSION"
  echo "Attach with: tmux attach -t $TMUX_SESSION"
  echo "Detach inside the session with: Ctrl-b d"
  exit 0
fi

# ---- Resolve sample list ----------------------------------------------------
if [[ -n "$SAMPLES_OVERRIDE" ]]; then
  SAMPLES_CSV="$SAMPLES_OVERRIDE"
else
  SAMPLES_CSV="$(awk -F',' 'NR>1 && ($9=="CART" || $9=="Conventional") && $12 != "composite" && $1 != "TEST_5k" {print $1}' "$SHEET" | paste -sd ',' -)"
fi

if [[ -z "$SAMPLES_CSV" ]]; then
  echo "ERROR: no samples to run (sheet filter returned empty and no --samples override)." >&2
  exit 2
fi

IFS=',' read -r -a SAMPLES <<< "$SAMPLES_CSV"

echo "Project dir : $PROJECT_DIR"
echo "Launcher    : $LAUNCHER"
echo "Sample sheet: $SHEET"
echo "Samples     : ${#SAMPLES[@]} (${SAMPLES_CSV})"
if [[ ${#PASSTHROUGH[@]} -gt 0 ]]; then
  echo "Pass-through: ${PASSTHROUGH[*]}"
fi
echo

# ---- Per-sample loop --------------------------------------------------------
# Indexed arrays in lockstep with SAMPLES (no `declare -A`; bash 3.2 friendly).
FAILED=()
EXIT_CODES=()
ELAPSED=()
LOG_PATHS=()
RUN_START=$(date +%s)

for i in "${!SAMPLES[@]}"; do
  SID="${SAMPLES[$i]}"
  POS=$((i + 1))
  TOTAL=${#SAMPLES[@]}
  LOG_PATH="$PROJECT_DIR/Output/Sample_${SID}/Pipeline_log_${SID}.log"
  LOG_PATHS+=("$LOG_PATH")
  echo "================================================================"
  echo "[$POS/$TOTAL] Sample: $SID"
  echo "Started     : $(date -Iseconds 2>/dev/null || date)"
  echo "Log (target): $LOG_PATH"
  echo "================================================================"

  CMD=("$LAUNCHER" --samples "$SID")
  if [[ ${#PASSTHROUGH[@]} -gt 0 ]]; then
    CMD+=("${PASSTHROUGH[@]}")
  fi

  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "DRY RUN: ${CMD[*]}"
    EXIT_CODES+=(0)
    ELAPSED+=(0)
    echo
    continue
  fi

  SAMPLE_START=$(date +%s)
  set +u
  "${CMD[@]}"
  RC=$?
  set -u
  SAMPLE_END=$(date +%s)
  EXIT_CODES+=("$RC")
  ELAPSED+=("$((SAMPLE_END - SAMPLE_START))")
  if [[ $RC -ne 0 ]]; then
    FAILED+=("$SID")
    echo ">> Sample $SID FAILED with exit $RC after ${ELAPSED[$i]}s. Continuing." >&2
  else
    echo ">> Sample $SID OK in ${ELAPSED[$i]}s."
  fi
  echo
done

RUN_END=$(date +%s)
TOTAL_ELAPSED=$((RUN_END - RUN_START))

# ---- Summary ----------------------------------------------------------------
echo "================================================================"
echo "Run summary"
echo "Total elapsed: ${TOTAL_ELAPSED}s"
echo "================================================================"
printf "%-20s %-6s %-8s %s\n" "sample_id" "exit" "secs" "log"
printf "%-20s %-6s %-8s %s\n" "---------" "----" "----" "---"
for i in "${!SAMPLES[@]}"; do
  printf "%-20s %-6s %-8s %s\n" "${SAMPLES[$i]}" "${EXIT_CODES[$i]:-?}" "${ELAPSED[$i]:-?}" "${LOG_PATHS[$i]}"
done

if [[ ${#FAILED[@]} -gt 0 ]]; then
  echo
  echo "FAILED (${#FAILED[@]}/${#SAMPLES[@]}): ${FAILED[*]}" >&2
  echo "Re-run failed samples with: $0 --samples $(IFS=,; echo "${FAILED[*]}")" >&2
  exit 1
fi

exit 0

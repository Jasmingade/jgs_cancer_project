#!/usr/bin/env bash
#SBATCH --job-name=03_sum
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --output=02_proteomics/logs/summary/summary_%A_%a.log
#SBATCH --error=02_proteomics/logs/summary/summary_%A_%a.err

set -eo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda deactivate || true

R_ENV="${R_ENV:-/home/people/s184275/r-env}"
SUMMARY_SCRIPT="02_proteomics/pipeline/plotting/05_cox_global_summary.R"
SUMMARY_PLOT_SCRIPT="02_proteomics/pipeline/plotting/06_cox_global_plots.R"
LOGDIR="02_proteomics/logs/summary"

COX_ROOT="${COX_ROOT:-02_proteomics/out/cox_filtered}"
SUMMARY_OUT_ROOT="${SUMMARY_OUT_ROOT:-02_proteomics/out/summary}"
SUMMARY_FDR="${SUMMARY_FDR:-0.05}"
RUN_SUMMARY_EXPORT="${RUN_SUMMARY_EXPORT:-true}"
RUN_SUMMARY_PLOTS="${RUN_SUMMARY_PLOTS:-true}"

mkdir -p "$LOGDIR" "$SUMMARY_OUT_ROOT"

is_true() {
  local flag
  flag="$(echo "${1:-}" | tr '[:upper:]' '[:lower:]')"
  case "$flag" in
    1|true|yes|y|on) return 0 ;;
    *) return 1 ;;
  esac
}

if [[ ! -f "$SUMMARY_SCRIPT" ]]; then
  echo "[ERROR] Summary script not found: $SUMMARY_SCRIPT" >&2
  exit 1
fi
if is_true "$RUN_SUMMARY_PLOTS" && [[ ! -f "$SUMMARY_PLOT_SCRIPT" ]]; then
  echo "[ERROR] Plotting script not found: $SUMMARY_PLOT_SCRIPT" >&2
  exit 1
fi
if [[ ! -d "$COX_ROOT" ]]; then
  echo "[ERROR] Cox root does not exist: $COX_ROOT" >&2
  exit 1
fi

root_name=$(basename "$COX_ROOT")
OUT_DIR="$SUMMARY_OUT_ROOT/$root_name"
mkdir -p "$OUT_DIR"

echo "===== Cox global summary launcher ====="
date
echo "R env: $R_ENV"
echo "Cox root: $COX_ROOT"
echo "Out dir: $OUT_DIR"
echo "FDR threshold: $SUMMARY_FDR"
echo "Run summary export: $RUN_SUMMARY_EXPORT"
echo "Run summary plots : $RUN_SUMMARY_PLOTS"
echo "======================================="

if is_true "$RUN_SUMMARY_EXPORT"; then
  START_TIME=$(date +%s)
  echo "[RUN] Starting summary export at $(date)"
  conda run -p "$R_ENV" --no-capture-output Rscript "$SUMMARY_SCRIPT" \
    "$COX_ROOT" \
    "$OUT_DIR" \
    "$SUMMARY_FDR"
  END_TIME=$(date +%s)
  echo "[DONE] Finished summary export in $(( END_TIME - START_TIME )) seconds"
else
  echo "[SKIP] Summary export disabled (RUN_SUMMARY_EXPORT=$RUN_SUMMARY_EXPORT)"
fi

if is_true "$RUN_SUMMARY_PLOTS"; then
  if [[ ! -f "$OUT_DIR/cox_results_significant.csv" ]]; then
    echo "[ERROR] Significant summary not found in $OUT_DIR. Run export first." >&2
    exit 1
  fi
  START_TIME=$(date +%s)
  echo "[RUN] Starting summary plotting at $(date)"
  conda run -p "$R_ENV" --no-capture-output Rscript "$SUMMARY_PLOT_SCRIPT" \
    "$OUT_DIR" \
    "$SUMMARY_FDR"
  END_TIME=$(date +%s)
  echo "[DONE] Finished summary plotting in $(( END_TIME - START_TIME )) seconds"
else
  echo "[SKIP] Summary plotting disabled (RUN_SUMMARY_PLOTS=$RUN_SUMMARY_PLOTS)"
fi

echo "===== Summary job done ====="
date

#!/usr/bin/env bash
#SBATCH --job-name=02_plots
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=02_proteomics/logs/plots/%x-%j.log
#SBATCH --error=02_proteomics/logs/plots/%x-%j.err

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh

R_ENV="${R_ENV:-/home/people/s184275/r-env}"
PLOT_DIR="02_proteomics/pipeline/plotting"
LOG_DIR="02_proteomics/logs/plots"

FILTERED_ROOT="${FILTERED_ROOT:-02_proteomics/out/preprocessed/filtered}"
SANITY_OUT="${SANITY_OUT:-02_proteomics/out/plots/sanity}"
EXPR_OUT="${EXPR_OUT:-02_proteomics/out/plots/expression}"

RUN_PLOT_SANITY=${RUN_PLOT_SANITY:-true}
RUN_PLOT_EXPRESSION=${RUN_PLOT_EXPRESSION:-false}

mkdir -p "$LOG_DIR" "$SANITY_OUT" "$EXPR_OUT"

echo "===== Proteomics plotting launcher ====="
date
echo "R env: $R_ENV"
echo "Filtered root: $FILTERED_ROOT"
echo "Sanity out: $SANITY_OUT (toggle: $RUN_PLOT_SANITY)"
echo "Expression out: $EXPR_OUT (toggle: $RUN_PLOT_EXPRESSION)"
echo "======================================="

run_plot() {
  local toggle="$1"
  local script="$2"
  shift 2
  if [[ "$toggle" != "true" ]]; then
    echo "[SKIP] $script disabled"
    return 0
  fi
  if [[ ! -d "$FILTERED_ROOT" ]]; then
    echo "[ERROR] Filtered expression root missing: $FILTERED_ROOT" >&2
    return 1
  fi
  echo "[RUN] $script"
  conda run -p "$R_ENV" --no-capture-output Rscript "$PLOT_DIR/$script" "$@"
}

run_plot "$RUN_PLOT_SANITY" "01_plot_sanity.R" "$FILTERED_ROOT" "$SANITY_OUT"
run_plot "$RUN_PLOT_EXPRESSION" "02_plot_expression.R" "$FILTERED_ROOT" "$EXPR_OUT"

echo "===== Proteomics plotting launcher done ====="
date

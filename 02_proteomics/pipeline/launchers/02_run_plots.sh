#!/usr/bin/env bash
#SBATCH --job-name=02_plots
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=36G
#SBATCH --output=02_proteomics/logs/plots/%x-%j.log
#SBATCH --error=02_proteomics/logs/plots/%x-%j.err

#set -eu pipefail
source ~/miniconda3/etc/profile.d/conda.sh

conda deactivate || true
conda deactivate || true

R_ENV="${R_ENV:-/home/people/s184275/r-env}"
PLOT_DIR="02_proteomics/pipeline/plotting"
LOG_DIR="02_proteomics/logs/plots"

FILTERED_ROOT="${FILTERED_ROOT:-02_proteomics/out/preprocessed/filtered}"
COX_ROOT="${COX_ROOT:-02_proteomics/out/ruv/univariate}"
COX_FILTERED_ROOT="${COX_FILTERED_ROOT:-02_proteomics/out/cox_filtered}"
PLOT_CI_WIDTH_MAX="${PLOT_CI_WIDTH_MAX:-10}"
SANITY_OUT="${SANITY_OUT:-02_proteomics/out/plots/sanity}"
EXPR_OUT="${EXPR_OUT:-02_proteomics/out/plots/expression}"
FOREST_OUT="${FOREST_OUT:-02_proteomics/out/plots/forest}"
RUN_PLOT_SANITY=${RUN_PLOT_SANITY:-true}
RUN_FILTER_COX_OUTPUT=${RUN_FILTER_COX_OUTPUT:-false}
RUN_PLOT_EXPRESSION=${RUN_PLOT_EXPRESSION:-false}
RUN_PLOT_FOREST=${RUN_PLOT_FOREST:-false}
FOREST_FDR_THRESHOLD=${FOREST_FDR_THRESHOLD:-0.05}
FOREST_TOP_N=${FOREST_TOP_N:-25}

mkdir -p "$LOG_DIR" "$SANITY_OUT" "$EXPR_OUT"

echo "===== Proteomics plotting launcher ====="
date
echo "R env: $R_ENV"
echo "Filtered root: $FILTERED_ROOT"
echo "Cox root: $COX_ROOT"
echo "Filtered Cox root: $COX_FILTERED_ROOT (toggle: $RUN_FILTER_COX_OUTPUT, max CI width: $PLOT_CI_WIDTH_MAX)"
echo "Sanity out: $SANITY_OUT (toggle: $RUN_PLOT_SANITY)"
echo "Expression out: $EXPR_OUT (toggle: $RUN_PLOT_EXPRESSION)"
echo "Forest out: $FOREST_OUT (toggle: $RUN_PLOT_FOREST, FDR $FOREST_FDR_THRESHOLD, top N $FOREST_TOP_N)"
echo "======================================="

EXPR_INPUT_ROOT="$COX_ROOT"

filter_cox_inputs() {
  if [[ "$RUN_FILTER_COX_OUTPUT" != "true" ]]; then
    echo "[SKIP] CI-width filtering disabled; using $COX_ROOT directly"
    EXPR_INPUT_ROOT="$COX_ROOT"
    return 0
  fi
  if [[ ! -d "$COX_ROOT" ]]; then
    echo "[ERROR] Cox root missing: $COX_ROOT" >&2
    return 1
  fi
  mkdir -p "$COX_FILTERED_ROOT"
  echo "[RUN] 00_filter_cox_ci.R (ci_width <= $PLOT_CI_WIDTH_MAX)"
  conda run -p "$R_ENV" --no-capture-output Rscript \
    "$PLOT_DIR/00_filter_cox_ci.R" \
    "$COX_ROOT" \
    "$COX_FILTERED_ROOT" \
    "$PLOT_CI_WIDTH_MAX"
  EXPR_INPUT_ROOT="$COX_FILTERED_ROOT"
}

run_plot() {
  local toggle="$1"
  local script="$2"
  local input_root="$3"
  local output_dir="$4"
  shift 2
  if [[ "$toggle" != "true" ]]; then
    echo "[SKIP] $script disabled"
    return 0
  fi
  if [[ ! -d "$input_root" ]]; then
    echo "[ERROR] Input root missing: $input_root" >&2
    return 1
  fi
  echo "[RUN] $script"
  conda run -p "$R_ENV" --no-capture-output Rscript "$PLOT_DIR/$script" "$input_root" "$output_dir"
}

filter_cox_inputs || exit 1
run_plot "$RUN_PLOT_SANITY" "01_plot_sanity.R" "$FILTERED_ROOT" "$SANITY_OUT"
run_plot "$RUN_PLOT_EXPRESSION" "02_plot_expression_bee.R" "$EXPR_INPUT_ROOT" "$EXPR_OUT"

if [[ "$RUN_PLOT_FOREST" == "true" ]]; then
  mkdir -p "$FOREST_OUT"
  echo "[RUN] 03_forest_by_cancer.R"
  conda run -p "$R_ENV" --no-capture-output Rscript \
    "$PLOT_DIR/03_forest_by_cancer.R" \
    "$EXPR_INPUT_ROOT" \
    "$FOREST_OUT" \
    "$FOREST_FDR_THRESHOLD" \
    "$FOREST_TOP_N"
fi

echo "===== Proteomics plotting launcher done ====="
date

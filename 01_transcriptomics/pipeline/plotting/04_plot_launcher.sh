#!/usr/bin/env bash
#SBATCH --job-name=05_plots
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=01_transcriptomics/logs/05_plots/%x-%A_%a.log
#SBATCH --error=01_transcriptomics/logs/05_plots/%x-%A_%a.err
#SBATCH --array=0-2

set -euo pipefail

# =====================================================================
# CONFIG
# =====================================================================
R_ENV="/home/people/s184275/r-env"

PLOTDIR="01_transcriptomics/pipeline/plotting"
LOGDIR="01_transcriptomics/logs/05_plots"
OUTDIR="01_transcriptomics/out/05_plots"

mkdir -p "$LOGDIR" "$OUTDIR"

# =====================================================================
# MODEL TOGGLES (SET AT SUBMISSION TIME)
#
# sbatch --export=RUN_PLOT_M0=true,RUN_PLOT_M1=false,...
# =====================================================================
export RUN_PLOT_M0=${RUN_PLOT_M0:-true}
export RUN_PLOT_M1=${RUN_PLOT_M1:-true}
export RUN_PLOT_M2=${RUN_PLOT_M2:-true}
export RUN_PLOT_M3=${RUN_PLOT_M3:-true}
export RUN_PLOT_M4=${RUN_PLOT_M4:-true}
export RUN_PLOT_M5=${RUN_PLOT_M5:-true}

# Default for whether models should use *_full results (for p-hists)
export M1_USE_FULL_RESULTS=${M1_USE_FULL_RESULTS:-true}
export M2_USE_FULL_RESULTS=${M2_USE_FULL_RESULTS:-false}
export M3_USE_FULL_RESULTS=${M3_USE_FULL_RESULTS:-false}

echo "[INFO] Plot toggles:"
echo "  M0: $RUN_PLOT_M0"
echo "  M1: $RUN_PLOT_M1"
echo "  M2: $RUN_PLOT_M2"
echo "  M3: $RUN_PLOT_M3"
echo "  M4: $RUN_PLOT_M4"
echo "  M1_USE_FULL_RESULTS (for p-hists): $M1_USE_FULL_RESULTS"
echo "  M2_USE_FULL_RESULTS (for p-hists): $M2_USE_FULL_RESULTS"
echo "  M3_USE_FULL_RESULTS (for p-hists): $M3_USE_FULL_RESULTS"
echo

# =====================================================================
# ARRAY INDEX → MODEL
#
# 0 → M0 (sanity)
# 1 → M1 (expression only)
# 2 → M2 (mutation only)
# 3 → M3 (expression + mutation)
# 4 → M4 (interaction)
# =====================================================================

TASK=${SLURM_ARRAY_TASK_ID}

echo "[INFO] Running plot task index: $TASK"

case $TASK in

  0)
    if [[ "$RUN_PLOT_M0" == "true" ]]; then
      echo "[RUN] M0 – sanity plots"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        "$PLOTDIR/05a_plot_m0_sanity.R"
    else
      echo "[SKIP] M0 disabled"
    fi
    ;;

  1)
    if [[ "$RUN_PLOT_M1" == "true" ]]; then
      echo "[RUN] M1 – expression-only plots"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        "$PLOTDIR/05b_plot_m1_expression.R"
    else
      echo "[SKIP] M1 disabled"
    fi
    ;;

  2)
    if [[ "$RUN_PLOT_M2" == "true" ]]; then
      echo "[RUN] M2 – mutation-only plots"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        "$PLOTDIR/05c_plot_m2_mutation.R"
    else
      echo "[SKIP] M2 disabled"
    fi
    ;;

  3)
    if [[ "$RUN_PLOT_M3" == "true" ]]; then
      echo "[RUN] M3 – combined expression + mutation plots"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        "$PLOTDIR/05d_plot_m3_combined.R"
    else
      echo "[SKIP] M3 disabled"
    fi
    ;;

  4)
    if [[ "$RUN_PLOT_M4" == "true" ]]; then
      echo "[RUN] M4 – interaction model plots"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        "$PLOTDIR/05e_plot_m4_interactions.R"
    else
      echo "[SKIP] M4 disabled"
    fi
    ;;

  *)
    echo "[ERROR] Invalid array index: $TASK"
    exit 1
    ;;
esac

echo "[DONE] Plot model $TASK completed."

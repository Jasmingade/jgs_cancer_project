#!/bin/bash
#SBATCH --job-name=run_all_plots
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=01_transcriptomics/logs/04_plots/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/04_plots/%x-%A.err

# ============================================================
# 04_run_all_plots.sh
# ------------------------------------------------------------
# Runs all visualization scripts for post-CoxPH analyses:
#   ① Mutation-per-gene survival (baseline)
#   ② Expression + mutation combined CoxPH
#   ③ Isoform × mutation interactions
#   ④ Expression-only baseline (gene / iso_log / iso_frac)
# ============================================================

set -euo pipefail
mkdir -p 01_transcriptomics/logs/04_plots

say() { echo "[$(date '+%H:%M:%S')] $*"; }

say "=== Starting all post-CoxPH plotting tasks ==="
date; echo "PWD: $(pwd)"

R_ENV_PATH="/home/people/s184275/r-env"

# Helper to run Rscript safely inside conda
run_rscript() {
  local script="$1"
  local label="$2"
  if [[ -f "$script" ]]; then
    say "[RUN] $label → $script"
    conda run -p "$R_ENV_PATH" --no-capture-output Rscript "$script"
    say "[DONE] $label completed successfully."
  else
    say "[WARN] $label script not found: $script"
  fi
}

# ============================================================
# ① Mutation-per-gene survival (baseline)
# ============================================================
run_rscript "01_transcriptomics/pipeline/scripts/04b_plot_mutation_baseline.R" \
            "Mutation baseline plots"

# ============================================================
# ② Expression + mutation combined CoxPH
# ============================================================
run_rscript "01_transcriptomics/pipeline/scripts/04c_plot_exp_mutation_combined.R" \
            "Expression+Mutation combined plots"

# ============================================================
# ③ Isoform × mutation interactions
# ============================================================
run_rscript "01_transcriptomics/pipeline/scripts/04d_plot_iso_mutation_interactions.R" \
            "Isoform×Mutation interaction plots"

# ============================================================
# ④ Expression-only baseline (gene / iso_log / iso_frac)
# ============================================================
run_rscript "01_transcriptomics/pipeline/scripts/04a_plot_baseline_expression.R" \
            "Expression-only baseline plots"

say "=== All plotting scripts finished successfully ==="
date

#!/usr/bin/env bash
#SBATCH --job-name=03d_postprocess_coxph
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=01_transcriptomics/logs/03_coxph/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/03_coxph/%x-%A.err
set -euo pipefail

# -------------------------------------------------------------------
# CONFIG
# -------------------------------------------------------------------
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}
SCRIPT="01_transcriptomics/pipeline/scripts/03d_postprocess_and_plot_all_coxph.R"

OUT_BASE="01_transcriptomics/out"
LOG_DIR="01_transcriptomics/logs/03_coxph"
mkdir -p "$LOG_DIR"

# Expected folders
MODELS=(
  "03_univariate_coxph"
  "03a_mutation_univariate_coxph"
  "03b_exp_mutation_univariate_coxph"
  "03c_iso_mut_univariate_coxph"
)
MIN_FILES=5   # minimum .cox_results.csv files required per model to consider it valid

# -------------------------------------------------------------------
# FUNCTIONS
# -------------------------------------------------------------------
ts() { date "+%Y-%m-%d %H:%M:%S"; }
say() { echo -e "[$(ts)] $*"; }

# -------------------------------------------------------------------
# HEADER
# -------------------------------------------------------------------
echo "============================================================"
say "[START] Post-processing and summary plotting of CoxPH results"
echo "------------------------------------------------------------"
echo "  R environment : $R_ENV"
echo "  Script        : $SCRIPT"
echo "  Output base   : $OUT_BASE"
echo "  Log directory : $LOG_DIR"
echo "------------------------------------------------------------"

# -------------------------------------------------------------------
# CHECK MODEL OUTPUTS
# -------------------------------------------------------------------
missing_models=()
total_files=0
for m in "${MODELS[@]}"; do
  dir="$OUT_BASE/$m"
  count=$(find "$dir" -type f -name "*.cox_results.csv" 2>/dev/null | wc -l || echo 0)
  total_files=$((total_files + count))
  if (( count < MIN_FILES )); then
    say "[WARN] $m appears incomplete: only $count result files found (<$MIN_FILES)"
    missing_models+=("$m")
  else
    say "[OK] $m contains $count result files"
  fi
done

echo "------------------------------------------------------------"
say "[INFO] Total result files detected: $total_files"
if (( ${#missing_models[@]} > 0 )); then
  say "[WARN] Missing/incomplete models:"
  for mm in "${missing_models[@]}"; do
    echo "  - $mm"
  done
  echo
  read -p "Continue post-processing anyway? (y/N): " resp
  [[ "$resp" =~ ^[Yy]$ ]] || { say "[ABORT] Please rerun missing analyses before post-processing."; exit 1; }
fi

# -------------------------------------------------------------------
# RUN POST-PROCESSING
# -------------------------------------------------------------------
say "[RUN] Launching R post-processing script..."
conda run -p "$R_ENV" --no-capture-output Rscript "$SCRIPT"

if [[ $? -eq 0 ]]; then
  say "[DONE] Post-processing completed successfully"
else
  say "[ERROR] Post-processing failed"
  exit 1
fi

echo "------------------------------------------------------------"
say "[SUMMARY] Outputs available under:"
echo "  - $OUT_BASE/03d_summary_plots/"
echo "------------------------------------------------------------"
say "[FINISH] $(ts)"
echo "============================================================"

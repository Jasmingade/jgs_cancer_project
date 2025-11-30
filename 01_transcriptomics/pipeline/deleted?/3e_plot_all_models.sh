#!/usr/bin/env bash
#SBATCH --job-name=03_plot_all_models
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=01_transcriptomics/logs/03_plots/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/03_plots/%x-%A.err
set -euo pipefail

# ============================================================
# CONFIGURATION
# ============================================================
BASE_OUT="01_transcriptomics/out"
PLOT_OUT_DIR="${BASE_OUT}/03d_summary_plots"
INDEX_TSV="${PLOT_OUT_DIR}/index_all_models.tsv"
PLOT_SCRIPT="01_transcriptomics/pipeline/scripts/03_plot_all_coxph.R"
CONDA_ENV_PATH="/home/people/s184275/r-env"

mkdir -p "${PLOT_OUT_DIR}" "01_transcriptomics/logs/03_plots"

# ============================================================
# STEP 1. Build the unified index file (all 4 model types)
# ============================================================
echo "[STEP 1] Building unified index file for all model outputs..."
find "${BASE_OUT}" -type f -name "*.cox_results.csv" | \
awk -v OFS='\t' '
{
  path=$0
  n=split(path, parts, "/")

  # detect model directory starting with 03
  model=""
  for (i=1; i<=n; i++) if (parts[i] ~ /^03/) { model=parts[i]; break }

  file=parts[n]
  sub(/\.cox_results\.csv$/, "", file)
  split(file, a, "_")
  cancer=a[2]
  dtype=a[length(a)]

  # normalize iso_log / iso_frac naming
  if (dtype == "log") dtype="iso_log"
  if (dtype == "frac") dtype="iso_frac"

  # assign model-based types
  if (model == "03a_mutation_univariate_coxph")       dtype="mutation_baseline"
  else if (model == "03b_exp_mutation_univariate_coxph") dtype="mutation_combined"
  else if (model == "03c_iso_mut_univariate_coxph")      dtype="mutation_interaction"

  print path, cancer, dtype
}' > "${INDEX_TSV}"

echo "[OK] Index file created: ${INDEX_TSV}"
echo "  → $(wc -l < ${INDEX_TSV}) entries"

# ============================================================
# STEP 2. Generate HR boxplots & summaries
# ============================================================
FIG_NAME="HR_boxplots_all_models.png"
FIG_TITLE="Hazard Ratio Distributions across Models"
WINSOR=0.01
CLIP="0.01,0.99"

args=(
  --index "${INDEX_TSV}"
  --out   "${PLOT_OUT_DIR}/${FIG_NAME}"
  --title "${FIG_TITLE}"
  --winsor "${WINSOR}"
  --clip "${CLIP}"
  --width 16
  --height 9
  --dpi 300
)

echo "[STEP 2] Running plotting script..."
conda run -p "${CONDA_ENV_PATH}" --no-capture-output \
  Rscript "${PLOT_SCRIPT}" "${args[@]}"
echo "[DONE] Plots saved under: ${PLOT_OUT_DIR}/"

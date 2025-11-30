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

# Control behavior
FORCE_REBUILD=${FORCE_REBUILD:-false}     # true = overwrite even if exists
SIGNIFICANCE_MODE=${SIGNIFICANCE_MODE:-"FDR"}  # FDR or p
SIG_THRESHOLD=${SIG_THRESHOLD:-"0.05"}    # threshold (e.g. 0.01 or 0.05)

mkdir -p "${PLOT_OUT_DIR}" "01_transcriptomics/logs/03_plots"

FIG_NAME="HR_boxplots_all_models_sigonly.png"
FIG_TITLE="Significant Hazard Ratio Distributions (FDR < 0.05)"
WINSOR=0.01
CLIP="0.01,0.99"
FIG_PATH="${PLOT_OUT_DIR}/${FIG_NAME}"

# ============================================================
# SKIP IF ALREADY DONE (unless forced)
# ============================================================
if [[ -f "${FIG_PATH}" && "${FORCE_REBUILD}" != "true" ]]; then
  echo "[SKIP] ${FIG_PATH} already exists. Use FORCE_REBUILD=true to overwrite."
  exit 0
fi

# ============================================================
# STEP 1. Build the unified index file (all 4 model types)
# ============================================================
echo "[STEP 1] Building unified index file for all model outputs..."
echo -e "path\tcancer\tdata_type" > "${INDEX_TSV}"

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
  if (model == "03a_mutation_univariate_coxph")          dtype="mutation_baseline"
  else if (model == "03b_exp_mutation_univariate_coxph") dtype="mutation_combined"
  else if (model == "03c_iso_mut_univariate_coxph")      dtype="mutation_interaction"

  print path, cancer, dtype
}' >> "${INDEX_TSV}"

echo "[OK] Index file created: ${INDEX_TSV}"
head -n 10 "${INDEX_TSV}" | sed 's/^/  /'
echo "  → $(($(wc -l < ${INDEX_TSV}) - 1)) entries (excluding header)"

# ============================================================
# STEP 2. Generate HR boxplots & summaries (only significant)
# ============================================================
if [[ "${SIGNIFICANCE_MODE}" == "FDR" ]]; then
  FILTER_EXPR="FDR < ${SIG_THRESHOLD}"
  FIG_TITLE="Significant Hazard Ratio Distributions (FDR < ${SIG_THRESHOLD})"
elif [[ "${SIGNIFICANCE_MODE}" == "p" ]]; then
  FILTER_EXPR="p < ${SIG_THRESHOLD}"
  FIG_TITLE="Significant Hazard Ratio Distributions (p < ${SIG_THRESHOLD})"
else
  echo "[WARN] Unknown SIGNIFICANCE_MODE=${SIGNIFICANCE_MODE}, defaulting to FDR < ${SIG_THRESHOLD}"
  FILTER_EXPR="FDR < ${SIG_THRESHOLD}"
fi

args=(
  --index "${INDEX_TSV}"
  --out   "${FIG_PATH}"
  --title "${FIG_TITLE}"
  --filter "${FILTER_EXPR}"
  --winsor "${WINSOR}"
  --clip "${CLIP}"
  --width 16
  --height 9
  --dpi 300
)

echo "[STEP 2] Running plotting script (filter=${FILTER_EXPR})..."
conda run -p "${CONDA_ENV_PATH}" --no-capture-output \
  Rscript "${PLOT_SCRIPT}" "${args[@]}"

echo "[DONE] Significant HR plots saved in: ${PLOT_OUT_DIR}/"

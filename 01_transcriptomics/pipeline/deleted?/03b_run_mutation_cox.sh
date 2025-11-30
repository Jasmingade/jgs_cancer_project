#!/usr/bin/env bash
#SBATCH --job-name=03b_mutation_cox
#SBATCH --cpus-per-task=6
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --array=0-300
#SBATCH --output=01_transcriptomics/logs/03_mut/%x-%A_%a.log
#SBATCH --error=01_transcriptomics/logs/03_mut/%x-%A_%a.err
set -euo pipefail

# ============================================================
# CONFIGURATION
# ============================================================
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}
CANCERS_YAML=${CANCERS_YAML:-"01_transcriptomics/config/cancers.yaml"}
NORM_DIR=${NORM_DIR:-"01_transcriptomics/out/02_norm"}
MUT_DIR=${MUT_DIR:-"01_transcriptomics/out/02_mutation/01_per_gene_baseline/split_by_cancer_type"}
OUT_DIR=${OUT_DIR:-"01_transcriptomics/out/03_univariate_coxph_mutation"}
SCRIPT=${SCRIPT:-"01_transcriptomics/pipeline/scripts/03b_mutation_univariate_coxph.R"}
COV_YAML=${COV_YAML:-"01_transcriptomics/config/covariates.yaml"}

mkdir -p "$OUT_DIR" "01_transcriptomics/logs/03_mut"

mapfile -t CANCERS < <(sed -n 's/^[[:space:]]*-[[:space:]]*//p' "$CANCERS_YAML")

# Determine cancer index
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  TASK=0
else
  TASK=${SLURM_ARRAY_TASK_ID}
fi
CANCER="${CANCERS[$TASK]}"
echo "[TASK] Running mutation CoxPH for: $CANCER"

# Input paths
MUT_IN="${MUT_DIR}/TCGA_${CANCER}_mutation.normalized.csv"
MANI_IN="${NORM_DIR}/TCGA_${CANCER}_gene.sample_manifest.csv"
OUT_RES="${OUT_DIR}/TCGA_${CANCER}_mutation.cox_results.csv"
OUT_SUM="${OUT_DIR}/TCGA_${CANCER}_mutation.summary.txt"

# Sanity checks
if [[ ! -f "$MUT_IN" ]]; then
  echo "[SKIP] No mutation covariates for $CANCER: $MUT_IN"
  exit 0
fi
if [[ ! -f "$MANI_IN" ]]; then
  echo "[SKIP] No manifest for $CANCER: $MANI_IN"
  exit 0
fi

echo "[RUN] $CANCER"
conda run -p "$R_ENV" --no-capture-output \
  Rscript "$SCRIPT" "$MUT_IN" "$MANI_IN" "$COV_YAML" "$OUT_RES" "$OUT_SUM"
echo "[DONE] $CANCER → $OUT_RES"

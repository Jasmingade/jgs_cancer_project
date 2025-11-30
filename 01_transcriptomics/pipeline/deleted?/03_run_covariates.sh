#!/usr/bin/env bash
#SBATCH --job-name=03_covariates_only
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --output=01_transcriptomics/logs/03_cov/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/03_cov/%x-%A.err
set -euo pipefail

CANCERS_YAML=${CANCERS_YAML:-"01_transcriptomics/config/cancers.yaml"}
MANI_DIR=${MANI_DIR:-"01_transcriptomics/out/02_norm_batch"}   # where *_sample_manifest.csv live
COV_YAML=${COV_YAML:-"01_transcriptomics/config/covariates.yaml"}
OUT_DIR=${OUT_DIR:-"01_transcriptomics/out/03_covariates"}
SCRIPT=${SCRIPT:-"01_transcriptomics/pipeline/scripts/03_covariates_only.R"}
CONDA_ENV_PATH=${CONDA_ENV_PATH:-"/home/people/s184275/r-env"}

mkdir -p "$OUT_DIR" "01_transcriptomics/logs/03_cov"
mapfile -t CANCERS < <(sed -n 's/^[[:space:]]*-[[:space:]]*//p' "$CANCERS_YAML")

for C in "${CANCERS[@]}"; do
  mani="${MANI_DIR}/TCGA_${C}_gene.sample_manifest.csv"  # any data_type manifest works; covariates identical
  [[ -f "$mani" ]] || { echo "[WARN] Missing manifest for $C"; continue; }
  out_tsv="${OUT_DIR}/TCGA_${C}_covariates.tsv"
  out_png="${OUT_DIR}/TCGA_${C}_covariates_forest.png"
  echo "[RUN] $C"
  conda run -p "$CONDA_ENV_PATH" --no-capture-output \
    Rscript "$SCRIPT" "$mani" "$COV_YAML" "$out_tsv" "$out_png"
done
echo "[DONE] Covariates-only models"

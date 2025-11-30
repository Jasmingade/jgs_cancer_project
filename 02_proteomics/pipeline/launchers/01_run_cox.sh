#!/usr/bin/env bash
#SBATCH --job-name=02_cox
#SBATCH --cpus-per-task=6
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --output=02_proteomics/logs/cox/%x-%j.log
#SBATCH --error=02_proteomics/logs/cox/%x-%j.err
#SBATCH --array=0-14

set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda deactivate || true


R_ENV="${R_ENV:-/home/people/s184275/r-env}"
EXPR_ROOT="${EXPR_ROOT:-02_proteomics/out/preprocessed/filtered}"
CLIN_DIR="${CLIN_DIR:-02_proteomics/data/clinical}"
COV_YAML="${COV_YAML:-02_proteomics/config/covariates.yaml}"
BATCH_DIR="${BATCH_DIR:-02_proteomics/data/batch_annotation}"
OUT_BASE="${OUT_BASE:-02_proteomics/out/cox}"

RUN_COX_GENE=${RUN_COX_GENE:-true}
RUN_COX_ISO_LOG=${RUN_COX_ISO_LOG:-true}
RUN_COX_ISO_FRAC=${RUN_COX_ISO_FRAC:-true}

mkdir -p "$OUT_BASE" "02_proteomics/logs/cox"

echo "===== Proteomics Cox launcher ====="
date
echo "R env: $R_ENV"
echo "Expr root: $EXPR_ROOT"
echo "Clinical dir: $CLIN_DIR"
echo "Batch dir: $BATCH_DIR"
echo "Covariates YAML: $COV_YAML"
echo "Out dir: $OUT_BASE"
echo "RUN_COX_GENE: $RUN_COX_GENE"
echo "RUN_COX_ISO_LOG: $RUN_COX_ISO_LOG"
echo "RUN_COX_ISO_FRAC: $RUN_COX_ISO_FRAC"
echo "==================================="

collect_datasets() {
  local root="$EXPR_ROOT/gene"
  if [[ ! -d "$root" ]]; then return; fi
  find "$root" -maxdepth 1 -type f -name '*_gene.csv' -printf '%f\n' 2>/dev/null |
    sed 's/_gene\.csv$//' | sort -u
}

run_dtype_for_dataset() {
  local dtype="$1"
  local toggle="$2"
  local dataset="$3"
  local clin_file="$4"
  if [[ "$toggle" != "true" ]]; then
    return
  fi
  local expr_file="${EXPR_ROOT}/${dtype}/${dataset}_${dtype}.csv"
  if [[ ! -f "$expr_file" ]]; then
    echo "[WARN] Missing ${dtype} matrix for ${dataset} → $expr_file"
    return
  fi
  local out_dir="$OUT_BASE/${dtype}/${dataset}"
  mkdir -p "$out_dir"
  local out_res="$out_dir/${dataset}.${dtype}.cox_results.csv"
  local out_sum="$out_dir/${dataset}.${dtype}.summary.txt"
  echo "[RUN] dtype=$dtype dataset=$dataset"
  conda run -p "$R_ENV" --no-capture-output Rscript \
    02_proteomics/pipeline/cox/01_univariate_coxph.R \
    "$expr_file" \
    "$clin_file" \
    "$COV_YAML" \
    "$BATCH_DIR" \
    "$out_res" \
    "$out_sum"
}

mapfile -t DATASETS < <(collect_datasets)
if [[ ${#DATASETS[@]} -eq 0 ]]; then
  echo "[WARN] No datasets found under $EXPR_ROOT/gene"
  exit 0
fi

tasks=("${DATASETS[@]}")
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  idx=${SLURM_ARRAY_TASK_ID}
  if (( idx >= ${#DATASETS[@]} )); then
    echo "[SKIP] Array index $idx out of range (datasets=${#DATASETS[@]})"
    exit 0
  fi
  tasks=("${DATASETS[$idx]}")
fi

for dataset in "${tasks[@]}"; do
  study_platform="$dataset"
  sample_type=""
  if [[ "$dataset" =~ _normal$ || "$dataset" =~ _tumor$ ]]; then
    sample_type="${dataset##*_}"
    study_platform="${dataset%_*}"
  fi
  study="${study_platform%%_*}"
  clin_file="${CLIN_DIR}/${study_platform}_clinical.csv"
  if [[ ! -f "$clin_file" ]]; then
    fallback="${CLIN_DIR}/${study}_clinical.csv"
    if [[ -f "$fallback" ]]; then
      clin_file="$fallback"
    fi
  fi
  if [[ ! -f "$clin_file" ]]; then
    echo "[WARN] Missing clinical file for $dataset → ${CLIN_DIR}/${study_platform}_clinical.csv"
    continue
  fi

  run_dtype_for_dataset gene "$RUN_COX_GENE" "$dataset" "$clin_file"
  run_dtype_for_dataset iso_log "$RUN_COX_ISO_LOG" "$dataset" "$clin_file"
  run_dtype_for_dataset iso_frac "$RUN_COX_ISO_FRAC" "$dataset" "$clin_file"
done

echo "===== Proteomics Cox launcher done ====="
date

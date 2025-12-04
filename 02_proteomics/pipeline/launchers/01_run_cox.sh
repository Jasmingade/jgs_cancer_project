#!/usr/bin/env bash
#SBATCH --job-name=02_cox
#SBATCH --cpus-per-task=6
#SBATCH --time=02:00:00
#SBATCH --mem=36G
#SBATCH --output=02_proteomics/logs/cox/cox_%A_%a.log
#SBATCH --error=02_proteomics/logs/cox/cox_%A_%a.err
#SBATCH --array=0-14

#set -eu pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda deactivate || true
conda deactivate || true


R_ENV="${R_ENV:-/home/people/s184275/r-env}"
EXPR_ROOT="${EXPR_ROOT:-02_proteomics/out/preprocessed_none/filtered}"
CLIN_DIR="${CLIN_DIR:-02_proteomics/data/clinical}"
COV_YAML="${COV_YAML:-02_proteomics/config/covariates.yaml}"


OUT_BASE="${OUT_BASE:-02_proteomics/out/cox_none_quan}"

RUN_COX_PER_TYPE=${RUN_COX_PER_TYPE:-true}
RUN_COX_UNIV=${RUN_COX_UNIV:-true}
RUN_COX_COXMOS=${RUN_COX_COXMOS:-false}
RUN_COX_GLMNET=${RUN_COX_GLMNET:-false}
RUN_COX_DEEPSURV=${RUN_COX_DEEPSURV:-false}

canonical_dataset_id() {
  local ds="$1"
  if [[ -z "$ds" ]]; then
    echo "$ds"
    return
  fi
  local suffix="${ds##*_}"
  local suffix_lower="${suffix,,}"
  if [[ "$suffix_lower" == "reference" ]]; then
    echo "$ds"
    return
  fi
  if [[ "$suffix_lower" == "normal" || "$suffix_lower" == "tumor" ]]; then
    echo "${ds%_*}"
    return
  fi
  echo "$ds"
}

mkdir -p "$OUT_BASE" "02_proteomics/logs/cox"

echo "===== Proteomics Cox launcher ====="
date
echo "R env: $R_ENV"
echo "Expr root: $EXPR_ROOT"
echo "Clinical dir: $CLIN_DIR"
echo "Covariates YAML: $COV_YAML"
echo "Out dir: $OUT_BASE"
echo "RUN_COX_PER_TYPE: $RUN_COX_PER_TYPE"
echo "RUN_COX_UNIV: $RUN_COX_UNIV"
echo "RUN_COX_COXMOS: $RUN_COX_COXMOS"
echo "RUN_COX_GLMNET: $RUN_COX_GLMNET"
echo "RUN_COX_DEEPSURV: $RUN_COX_DEEPSURV"
echo "==================================="
echo "[INFO] Methods enabled → univariate: $RUN_COX_UNIV, coxmos: $RUN_COX_COXMOS, glmnet: $RUN_COX_GLMNET, deepsurv: $RUN_COX_DEEPSURV"
if [[ "$RUN_COX_PER_TYPE" != "true" && "$RUN_COX_TUMOR_NORMAL" != "true" ]]; then
  echo "[WARN] Both RUN_COX_PER_TYPE and RUN_COX_TUMOR_NORMAL are false → nothing to run"
  exit 0
fi

collect_datasets() {
  local root="$EXPR_ROOT/gene"
  if [[ ! -d "$root" ]]; then return; fi
  find "$root" -maxdepth 1 -type f -name '*_gene.csv' -printf '%f\n' 2>/dev/null |
    sed 's/_gene\.csv$//' |
    { grep -v '_reference$' || true; } |
    sort -u
}

run_dtype_for_dataset() {
  local dtype="$1"
  local mode="$2"   # univariate | coxmos | glmnet | deepsurv
  local toggle="$3"
  local dataset_src="$4"
  local dataset_key="$5"
  local clin_file="$6"
  if [[ "$toggle" != "true" ]]; then
    return
  fi
  local expr_file="${EXPR_ROOT}/${dtype}/${dataset_key}_${dtype}.csv"
  if [[ ! -f "$expr_file" && "$dataset_src" != "$dataset_key" ]]; then
    local legacy="${EXPR_ROOT}/${dtype}/${dataset_src}_${dtype}.csv"
    if [[ -f "$legacy" ]]; then
      expr_file="$legacy"
    fi
  fi
  if [[ ! -f "$expr_file" ]]; then
    echo "[WARN] Missing ${dtype} matrix for ${dataset_key} → ${EXPR_ROOT}/${dtype}/${dataset_key}_${dtype}.csv"
    return
  fi
  local method_dir="$OUT_BASE"
  case "$mode" in
    coxmos)   method_dir="$OUT_BASE/coxmos" ;;
    glmnet)   method_dir="$OUT_BASE/glmnet" ;;
    deepsurv) method_dir="$OUT_BASE/deepsurv" ;;
    *)        method_dir="$OUT_BASE/univariate" ;;
  esac
  local out_dir="$method_dir/${dtype}/${dataset_key}"
  mkdir -p "$out_dir"
  case "$mode" in
    coxmos)
      local out_res="$out_dir/${dataset_key}.${dtype}.coxmos.csv"
      local out_sum="$out_dir/${dataset_key}.${dtype}.coxmos.summary.txt"
      echo "[RUN] coxmos dtype=$dtype dataset=$dataset_key → $out_res"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        02_proteomics/pipeline/cox/02_coxmos_cox.R \
        "$expr_file" \
        "$clin_file" \
        "$COV_YAML" \
        "$out_res" \
        "$out_sum"
      ;;
    glmnet)
      local out_res="$out_dir/${dataset_key}.${dtype}.glmnet.csv"
      local out_sum="$out_dir/${dataset_key}.${dtype}.glmnet.summary.txt"
      echo "[RUN] glmnet dtype=$dtype dataset=$dataset_key → $out_res"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        02_proteomics/pipeline/cox/03_glmnet_cox.R \
        "$expr_file" \
        "$clin_file" \
        "$COV_YAML" \
        "$out_res" \
        "$out_sum"
      ;;
    deepsurv)
      local out_res="$out_dir/${dataset_key}.${dtype}.deepsurv.csv"
      local out_sum="$out_dir/${dataset_key}.${dtype}.deepsurv.summary.txt"
      echo "[RUN] deepsurv dtype=$dtype dataset=$dataset_key → $out_res"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        02_proteomics/pipeline/cox/04_deepsurv_cox.R \
        "$expr_file" \
        "$clin_file" \
        "$COV_YAML" \
        "$out_res" \
        "$out_sum"
      ;;
    *)
      local out_res="$out_dir/${dataset_key}.${dtype}.cox_results.csv"
      local out_sum="$out_dir/${dataset_key}.${dtype}.summary.txt"
      echo "[RUN] univariate dtype=$dtype dataset=$dataset_key → $out_res"
      conda run -p "$R_ENV" --no-capture-output Rscript \
        02_proteomics/pipeline/cox/01_univariate_coxph.R \
        "$expr_file" \
        "$clin_file" \
        "$COV_YAML" \
        "$out_res" \
        "$out_sum"
      ;;
  esac
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

if [[ "$RUN_COX_PER_TYPE" == "true" ]]; then
for dataset in "${tasks[@]}"; do
    dataset_src="$dataset"
    dataset_key="$(canonical_dataset_id "$dataset_src")"
    study_platform="$dataset_key"
    study="${study_platform%%_*}"
    clin_file="${CLIN_DIR}/${study_platform}_clinical.csv"
    if [[ ! -f "$clin_file" && "$dataset_src" != "$dataset_key" ]]; then
      legacy_clin="${CLIN_DIR}/${dataset_src}_clinical.csv"
      if [[ -f "$legacy_clin" ]]; then
        clin_file="$legacy_clin"
      fi
    fi
    if [[ ! -f "$clin_file" ]]; then
      fallback="${CLIN_DIR}/${study}_clinical.csv"
      if [[ -f "$fallback" ]]; then
        clin_file="$fallback"
      fi
    fi
    if [[ ! -f "$clin_file" ]]; then
      echo "[WARN] Missing clinical file for $dataset_key → ${CLIN_DIR}/${study_platform}_clinical.csv"
      continue
    fi

    for dtype in gene iso_log iso_frac; do
      run_dtype_for_dataset "$dtype" "univariate" "$RUN_COX_UNIV" "$dataset_src" "$dataset_key" "$clin_file"
      run_dtype_for_dataset "$dtype" "coxmos" "$RUN_COX_COXMOS" "$dataset_src" "$dataset_key" "$clin_file"
      run_dtype_for_dataset "$dtype" "glmnet" "$RUN_COX_GLMNET" "$dataset_src" "$dataset_key" "$clin_file"
      run_dtype_for_dataset "$dtype" "deepsurv" "$RUN_COX_DEEPSURV" "$dataset_src" "$dataset_key" "$clin_file"
    done
  done
fi

echo "===== Proteomics Cox launcher done ====="
date

#!/usr/bin/env bash
#SBATCH --job-name=02_pre
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=36G
#SBATCH --output=02_proteomics/logs/preprocess/preprocess_%A_%a.log
#SBATCH --error=02_proteomics/logs/preprocess/preprocess_%A_%a.err
#SBATCH --array=1-1

set -eo pipefail

#source ~/miniconda3/bin/activate
source ~/miniconda3/etc/profile.d/conda.sh

conda deactivate || true
conda deactivate || true

# =====================================================================
# CONFIGURATION
# Override via: sbatch --export=ALL,RUN_SPLIT_CLINICAL=true,PROT_CLINICAL_ALL=...
# =====================================================================
R_ENV="${R_ENV:-/home/people/s184275/r-env}"
PREPROCESS_DIR="02_proteomics/pipeline/preprocess"
LOGDIR="02_proteomics/logs/preprocess"
PROT_STUDY_LIST="${PROT_STUDY_LIST:-02_proteomics/data/raw/study_list.txt}"
PROT_ENST_ENSG_MAP="${PROT_ENST_ENSG_MAP:-02_proteomics/data/raw/ENST-ENSG_mapping.csv}"

# Toggle scripts
RUN_SPLIT_CLINICAL=${RUN_SPLIT_CLINICAL:-false}
RUN_LOG_TRANSFORM=${RUN_LOG_TRANSFORM:-true}
RUN_CYCLIC_LOESS=${RUN_CYCLIC_LOESS:-true}
RUN_RUV_BATCH=${RUN_RUV_BATCH:-true}
RUN_FILTER_POS_EXPR=${RUN_FILTER_POS_EXPR:-true}
RUN_PREPROCESS_PLOTS=${RUN_PREPROCESS_PLOTS:-false}

# Splitting clinical
PROT_CLINICAL_ALL=${PROT_CLINICAL_ALL:-02_proteomics/data/raw/clinical_all_merged.csv}
PROT_GENE_DIR=${PROT_GENE_DIR:-02_proteomics/data/gene}
PROT_CLINICAL_DIR=${PROT_CLINICAL_DIR:-02_proteomics/data/clinical}
PROT_CANCER_CONFIG=${PROT_CANCER_CONFIG:-02_proteomics/config/cancers.yaml}

# Log transform
PROT_LOG_IN_DIR=${PROT_LOG_IN_DIR:-02_proteomics/data}
PROT_LOG_OUT_DIR=${PROT_LOG_OUT_DIR:-02_proteomics/out/preprocessed/log_transform}
PROT_LOG_DATATYPES=${PROT_LOG_DATATYPES:-gene,iso_log}

# Normalization
PROT_NORM_DATA_ROOT=${PROT_NORM_DATA_ROOT:-$PROT_LOG_OUT_DIR}
PROT_NORM_OUT_DIR=${PROT_NORM_OUT_DIR:-02_proteomics/out/preprocessed/normalization}
PROT_NORM_DATATYPES=${PROT_NORM_DATATYPES:-gene,iso_log}
PROT_BATCH_ANNOT_DIR=${PROT_BATCH_ANNOT_DIR:-02_proteomics/data/batch_annotation}

# Batch correction
PROT_BATCH_IN_DIR=${PROT_BATCH_IN_DIR:-$PROT_NORM_OUT_DIR}
PROT_BATCH_OUT_DIR=${PROT_BATCH_OUT_DIR:-02_proteomics/out/preprocessed/batch_corrected}
PROT_BATCH_DATATYPES=${PROT_BATCH_DATATYPES:-gene,iso_log}
PROT_RUV_K=${PROT_RUV_K:-1}
PROT_BATCH_METHOD=${PROT_BATCH_METHOD:-none}   # ruv, combat, or none
PROT_BATCH_ANNOT_DIR=${PROT_BATCH_ANNOT_DIR:-02_proteomics/data/batch_annotation}

# Filtering
PROT_FILTER_IN_DIR_WAS_DEFAULT="false"
if [[ -z "${PROT_FILTER_IN_DIR+x}" ]]; then
  PROT_FILTER_IN_DIR="$PROT_BATCH_OUT_DIR"
  PROT_FILTER_IN_DIR_WAS_DEFAULT="true"
fi
PROT_FILTER_OUT_DIR=${PROT_FILTER_OUT_DIR:-02_proteomics/out/preprocessed/filtered}
PROT_FILTER_MIN_PROP=${PROT_FILTER_MIN_PROP:-0.2}
PROT_FILTER_DATATYPES=${PROT_FILTER_DATATYPES:-gene,iso_log}
PROT_FILTER_CLIN_DIR=${PROT_FILTER_CLIN_DIR:-02_proteomics/data/clinical}
PROT_PLOT_RAW_DIR=${PROT_PLOT_RAW_DIR:-02_proteomics/data}
PROT_PLOT_OUT_DIR=${PROT_PLOT_OUT_DIR:-02_proteomics/out/preprocessed/plots/preprocess_summary}

BATCH_CORRECTION_METHOD_NONE="false"
if [[ "${PROT_BATCH_METHOD,,}" == "none" ]]; then
  BATCH_CORRECTION_METHOD_NONE="true"
  if [[ "$PROT_FILTER_IN_DIR_WAS_DEFAULT" == "true" ]]; then
    PROT_FILTER_IN_DIR="$PROT_NORM_OUT_DIR"
  fi
fi

collect_dataset_ids() {
  local root="$1"
  local gene_root="$root/gene"
  if [[ ! -d "$gene_root" ]]; then return; fi
  find "$gene_root" -maxdepth 1 -type f -name '*_gene.csv' -printf '%f\n' 2>/dev/null |
    sed 's/_gene\.csv$//' | sort -u
}

collect_dataset_ids_no_ref() {
  collect_dataset_ids "$1" | { grep -v '_reference$' || true; }
}

mkdir -p "$LOGDIR" "$PROT_CLINICAL_DIR"

echo "===== Proteomics preprocess launcher ====="
date
echo "PWD: $(pwd)"
echo "R env: $R_ENV"
echo "Clinical merged: $PROT_CLINICAL_ALL"
echo "Gene dir: $PROT_GENE_DIR"
echo "Clinical out dir: $PROT_CLINICAL_DIR"
echo "Cancer config: $PROT_CANCER_CONFIG"
echo "RUN_LOG_TRANSFORM: $RUN_LOG_TRANSFORM"
echo "RUN_SPLIT_CLINICAL: $RUN_SPLIT_CLINICAL"
echo "RUN_CYCLIC_LOESS: $RUN_CYCLIC_LOESS"
echo "RUN_RUV_BATCH: $RUN_RUV_BATCH"
echo "RUN_FILTER_POS_EXPR: $RUN_FILTER_POS_EXPR"
echo "RUN_PREPROCESS_PLOTS: $RUN_PREPROCESS_PLOTS"
echo "PROT_LOG_IN_DIR: $PROT_LOG_IN_DIR"
echo "PROT_LOG_OUT_DIR: $PROT_LOG_OUT_DIR"
echo "PROT_LOG_DATATYPES: $PROT_LOG_DATATYPES"
echo "PROT_NORM_DATA_ROOT: $PROT_NORM_DATA_ROOT"
echo "PROT_NORM_OUT_DIR: $PROT_NORM_OUT_DIR"
echo "PROT_NORM_DATATYPES: $PROT_NORM_DATATYPES"
echo "PROT_BATCH_IN_DIR: $PROT_BATCH_IN_DIR"
echo "PROT_BATCH_OUT_DIR: $PROT_BATCH_OUT_DIR"
echo "PROT_BATCH_DATATYPES: $PROT_BATCH_DATATYPES"
echo "PROT_RUV_K: $PROT_RUV_K"
echo "PROT_BATCH_METHOD: $PROT_BATCH_METHOD"
echo "PROT_BATCH_ANNOT_DIR: $PROT_BATCH_ANNOT_DIR"
echo "PROT_FILTER_IN_DIR: $PROT_FILTER_IN_DIR"
echo "PROT_FILTER_OUT_DIR: $PROT_FILTER_OUT_DIR"
echo "PROT_FILTER_MIN_PROP: $PROT_FILTER_MIN_PROP"
echo "PROT_FILTER_DATATYPES: $PROT_FILTER_DATATYPES"
echo "PROT_FILTER_CLIN_DIR: $PROT_FILTER_CLIN_DIR"
echo "PROT_PLOT_RAW_DIR: $PROT_PLOT_RAW_DIR"
echo "PROT_PLOT_OUT_DIR: $PROT_PLOT_OUT_DIR"
echo "PROT_STUDY_LIST: $PROT_STUDY_LIST"
echo "=========================================="
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo "SLURM job: $SLURM_JOB_ID (array index=${SLURM_ARRAY_TASK_ID:-none})"
fi

CURRENT_STUDY=""
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  if [[ ! -f "$PROT_STUDY_LIST" ]]; then
    echo "[ERROR] Study list not found: $PROT_STUDY_LIST" >&2
    exit 1
  fi
  mapfile -t STUDY_IDS < <(grep -v '^\s*#' "$PROT_STUDY_LIST" | grep -v '^\s*$')
  if [[ ${#STUDY_IDS[@]} -eq 0 ]]; then
    echo "[ERROR] Study list is empty: $PROT_STUDY_LIST" >&2
    exit 1
  fi
  idx=${SLURM_ARRAY_TASK_ID}
  if (( idx >= ${#STUDY_IDS[@]} )); then
    echo "[SKIP] Array index $idx exceeds study list entries (${#STUDY_IDS[@]})"
    exit 0
  fi
  CURRENT_STUDY="${STUDY_IDS[$idx]}"
  echo "[INFO] Selected study for this task: $CURRENT_STUDY (index ${SLURM_ARRAY_TASK_ID})"
fi

FILTERED_DATASETS=()
filter_datasets_for_study() {
  local study="$1"
  shift || true
  FILTERED_DATASETS=()
  if [[ -z "$study" ]]; then
    FILTERED_DATASETS=("$@")
    return 0
  fi
  local item
  for item in "$@"; do
    if [[ "$item" == "$study"* ]]; then
      FILTERED_DATASETS+=("$item")
    fi
  done
}

run_split_clinical() {
  local script="$PREPROCESS_DIR/00_split_clinical_by_study.R"
  if [[ "$RUN_SPLIT_CLINICAL" != "true" ]]; then
    echo "[SKIP] Split clinical-by-study disabled"
    return 0
  fi

  if [[ ! -s "$PROT_CLINICAL_ALL" ]]; then
    echo "[ERROR] Missing merged clinical CSV: $PROT_CLINICAL_ALL" >&2
    return 1
  fi

  if [[ ! -d "$PROT_GENE_DIR" ]]; then
    echo "[ERROR] Missing proteomics expression directory: $PROT_GENE_DIR" >&2
    return 1
  fi
  if [[ ! -f "$PROT_CANCER_CONFIG" ]]; then
    echo "[WARN] Cancer config not found: $PROT_CANCER_CONFIG (continuing)" >&2
  fi

  echo "[RUN] Split proteomics clinical metadata by study"
  conda run -p "$R_ENV" --no-capture-output Rscript "$script" \
    "$PROT_CLINICAL_ALL" "$PROT_GENE_DIR" "$PROT_CLINICAL_DIR" "$PROT_CANCER_CONFIG"
}

run_cyclic_loess_norm() {
  local script="$PREPROCESS_DIR/01_cyclic_loess_normalize.R"
  if [[ "$RUN_CYCLIC_LOESS" != "true" ]]; then
    echo "[SKIP] Cyclic loess normalization disabled"
    return 0
  fi

  if [[ ! -d "$PROT_NORM_DATA_ROOT" ]]; then
    echo "[ERROR] Missing proteomics data root: $PROT_NORM_DATA_ROOT" >&2
    return 1
  fi

  mapfile -t DATASETS < <(collect_dataset_ids "$PROT_NORM_DATA_ROOT")
  if [[ ${#DATASETS[@]} -eq 0 ]]; then
    echo "[WARN] No datasets found under $PROT_NORM_DATA_ROOT/gene"
    return 0
  fi

  mkdir -p "$PROT_NORM_OUT_DIR"
  filter_datasets_for_study "$CURRENT_STUDY" "${DATASETS[@]}"
  local tasks=("${FILTERED_DATASETS[@]}")
  if [[ ${#tasks[@]} -eq 0 ]]; then
    if [[ -n "$CURRENT_STUDY" ]]; then
      echo "[SKIP] No normalization datasets found for study $CURRENT_STUDY"
    else
      echo "[SKIP] No normalization datasets available"
    fi
    return 0
  fi

  for dataset in "${tasks[@]}"; do
    echo "[RUN] Normalization for $dataset"
    conda run -p "$R_ENV" --no-capture-output Rscript "$script" \
      "$PROT_NORM_DATA_ROOT" "$PROT_NORM_OUT_DIR" "$PROT_NORM_DATATYPES" "$dataset"
  done
}

run_ruv_batch() {
  local method="${PROT_BATCH_METHOD,,}"
  local script
  case "$method" in
    none)
      echo "[SKIP] Batch correction disabled via PROT_BATCH_METHOD=none"
      return 0
      ;;
    ruv|"")
      script="$PREPROCESS_DIR/03_ruv_batch_correct.R"
      ;;
    combat)
      script="$PREPROCESS_DIR/03_combat_batch_correct.R"
      ;;
    *)
      echo "[ERROR] Unknown PROT_BATCH_METHOD: $PROT_BATCH_METHOD (expected ruv or combat)" >&2
      return 1
      ;;
  esac
  if [[ "$RUN_RUV_BATCH" != "true" ]]; then
    echo "[SKIP] Batch correction disabled"
    return 0
  fi

  mapfile -t DATASETS < <(collect_dataset_ids_no_ref "$PROT_BATCH_IN_DIR")
  if [[ ${#DATASETS[@]} -eq 0 ]]; then
    echo "[WARN] No datasets found for batch correction in $PROT_BATCH_IN_DIR/gene"
    return 0
  fi

  filter_datasets_for_study "$CURRENT_STUDY" "${DATASETS[@]}"
  local targets=("${FILTERED_DATASETS[@]}")
  if [[ ${#targets[@]} -eq 0 ]]; then
    if [[ -n "$CURRENT_STUDY" ]]; then
      echo "[SKIP] No batch-correction datasets found for study $CURRENT_STUDY"
    else
      echo "[SKIP] No batch-correction datasets available"
    fi
    return 0
  fi

  mkdir -p "$PROT_BATCH_OUT_DIR"
  for dataset in "${targets[@]}"; do
    echo "[RUN] ${PROT_BATCH_METHOD^^} batch correction for $dataset"
    if [[ "$method" == "ruv" ]]; then
      conda run -p "$R_ENV" --no-capture-output Rscript "$script" \
        "$PROT_BATCH_IN_DIR" "$PROT_BATCH_OUT_DIR" "$PROT_RUV_K" "$PROT_BATCH_DATATYPES" "$dataset" "$PROT_BATCH_ANNOT_DIR"
    else
      conda run -p "$R_ENV" --no-capture-output Rscript "$script" \
        "$PROT_BATCH_IN_DIR" "$PROT_BATCH_OUT_DIR" "$PROT_BATCH_DATATYPES" "$dataset" "$PROT_BATCH_ANNOT_DIR"
    fi
  done
}

# Log transform
run_log_transform() {
  local script="$PREPROCESS_DIR/00_log_transform.R"
  if [[ "$RUN_LOG_TRANSFORM" != "true" ]]; then
    echo "[SKIP] Log transform disabled"
    return 0
  fi

  if [[ ! -d "$PROT_LOG_IN_DIR" ]]; then
    echo "[ERROR] Missing proteomics input directory: $PROT_LOG_IN_DIR" >&2
    return 1
  fi

  mapfile -t DATASETS < <(collect_dataset_ids "$PROT_LOG_IN_DIR")
  if [[ ${#DATASETS[@]} -eq 0 ]]; then
    echo "[WARN] No datasets found under $PROT_LOG_IN_DIR/gene"
    return 0
  fi

  mkdir -p "$PROT_LOG_OUT_DIR"
  filter_datasets_for_study "$CURRENT_STUDY" "${DATASETS[@]}"
  local tasks=("${FILTERED_DATASETS[@]}")
  if [[ ${#tasks[@]} -eq 0 ]]; then
    if [[ -n "$CURRENT_STUDY" ]]; then
      echo "[SKIP] No log-transform datasets found for study $CURRENT_STUDY"
    else
      echo "[SKIP] No log-transform datasets available"
    fi
    return 0
  fi

  local joined
  joined=$(echo "$PROT_LOG_DATATYPES" | tr -d ' ')

  for ds in "${tasks[@]}"; do
    echo "[RUN] Log transform for $ds"
    conda run -p "$R_ENV" --no-capture-output Rscript "$script" \
      "$PROT_LOG_IN_DIR" \
      "$PROT_LOG_OUT_DIR" \
      "$joined" \
      "$ds"
  done
}

run_filter_positive_expr() {
  local script="$PREPROCESS_DIR/02_filter_positive_expression.R"
  if [[ "$RUN_FILTER_POS_EXPR" != "true" ]]; then
    echo "[SKIP] Feature filtering disabled"
    return 0
  fi

  if [[ ! -d "$PROT_FILTER_IN_DIR" ]]; then
    echo "[ERROR] Filter input dir missing: $PROT_FILTER_IN_DIR" >&2
    return 1
  fi

  mkdir -p "$PROT_FILTER_OUT_DIR"
  local target_prefix=""
  if [[ -n "$CURRENT_STUDY" ]]; then
    target_prefix="$CURRENT_STUDY"
    echo "[RUN] Filtering datasets for study prefix $CURRENT_STUDY"
  else
    echo "[RUN] Filtering all datasets"
  fi

  conda run -p "$R_ENV" --no-capture-output Rscript "$script" \
    "$PROT_FILTER_IN_DIR" \
    "$PROT_FILTER_OUT_DIR" \
    "$PROT_FILTER_MIN_PROP" \
    "$PROT_FILTER_DATATYPES" \
    "$PROT_FILTER_CLIN_DIR" \
    "$target_prefix"
}

script="$PREPROCESS_DIR/04_plots.R"
run_preprocess_plots() {
  local script="$PREPROCESS_DIR/04_plots.R"
  if [[ "$RUN_PREPROCESS_PLOTS" != "true" ]]; then
    echo "[SKIP] Preprocess plotting disabled"
    return 0
  fi
  if [[ "$BATCH_CORRECTION_METHOD_NONE" == "true" ]]; then
    echo "[SKIP] Preprocess plotting requires batch correction outputs; skipped because PROT_BATCH_METHOD=none"
    return 0
  fi
  if [[ ! -f "$script" ]]; then
    echo "[ERROR] Plot script not found: $script" >&2
    return 1
  fi
  if [[ ! -d "$PROT_NORM_OUT_DIR/gene" ]]; then
    echo "[ERROR] Missing normalized gene directory: $PROT_NORM_OUT_DIR/gene" >&2
    return 1
  fi
  if [[ ! -d "$PROT_BATCH_OUT_DIR/gene" ]]; then
    echo "[ERROR] Missing batch-corrected gene directory: $PROT_BATCH_OUT_DIR/gene" >&2
    return 1
  fi
  mkdir -p "$PROT_PLOT_OUT_DIR"
  if [[ -n "$CURRENT_STUDY" ]]; then
    echo "[RUN] Preprocess plotting for study $CURRENT_STUDY"
    conda run -p "$R_ENV" --no-capture-output Rscript "$script" \
      "$PROT_PLOT_RAW_DIR" \
      "$PROT_NORM_OUT_DIR" \
      "$PROT_BATCH_OUT_DIR" \
      "$PROT_BATCH_ANNOT_DIR" \
      "$PROT_PLOT_OUT_DIR" \
      "$CURRENT_STUDY"
  else
    echo "[RUN] Preprocess plotting for all studies"
    conda run -p "$R_ENV" --no-capture-output Rscript "$script" \
      "$PROT_PLOT_RAW_DIR" \
      "$PROT_NORM_OUT_DIR" \
      "$PROT_BATCH_OUT_DIR" \
      "$PROT_BATCH_ANNOT_DIR" \
      "$PROT_PLOT_OUT_DIR"
  fi
}

# Add more preprocess steps here as they are implemented
run_split_clinical
run_log_transform
run_cyclic_loess_norm
run_ruv_batch
run_filter_positive_expr
run_preprocess_plots

echo "===== Proteomics preprocess launcher done ====="
date

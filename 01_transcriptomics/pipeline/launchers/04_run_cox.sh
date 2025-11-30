#!/usr/bin/env bash
#SBATCH --job-name=04_cox_models
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --output=01_transcriptomics/logs/04_cox/%x-%A_%a.log
#SBATCH --error=01_transcriptomics/logs/04_cox/%x-%A_%a.err
#SBATCH --array=0-99

set -euo pipefail

# =====================================================================
# CONFIG
# =====================================================================
R_ENV="/home/people/s184275/r-env"

CFG_YAML="01_transcriptomics/config/cancers.yaml"
TX2GENE="01_transcriptomics/data/raw/tx2gene.csv"

EXPR_DIR="01_transcriptomics/out/02_norm"
MUT_DIR="01_transcriptomics/out/03_mutation"
OUT_BASE="01_transcriptomics/out/04_coxph"

mkdir -p "$OUT_BASE" "01_transcriptomics/logs/04_cox"

# =====================================================================
# MODEL TOGGLES (DEFAULTS)
#
# >>> Set these at job submit time:
# sbatch --export=RUN_MODEL1_EXPR_ONLY=true,RUN_MODEL2_MUT_ONLY=false ...
#
# If unset → defaults below:
# =====================================================================
export RUN_MODEL1_EXPR_ONLY=${RUN_MODEL1_EXPR_ONLY:-false}
export RUN_MODEL2_MUT_ONLY=${RUN_MODEL2_MUT_ONLY:-true}
export RUN_MODEL3_EXPR_PLUS_MUT=${RUN_MODEL3_EXPR_PLUS_MUT:-false}
export RUN_MODEL4_ISO_MUT_INTERACT=${RUN_MODEL4_ISO_MUT_INTERACT:-false}

echo "[INFO] Model selection:"
echo "  Model 1 expr-only:           $RUN_MODEL1_EXPR_ONLY"
echo "  Model 2 mutation-only:       $RUN_MODEL2_MUT_ONLY"
echo "  Model 3 expr + mut:          $RUN_MODEL3_EXPR_PLUS_MUT"
echo "  Model 4 isoform × mutation:  $RUN_MODEL4_ISO_MUT_INTERACT"
echo

# =====================================================================
# PARSE YAML: CANCERS + DATATYPES
# =====================================================================
readarray -t CANCERS < <(
  awk '/^cancers:/{f=1;next} /^[^ ]/{if(f)exit} 
       f {sub("- ",""); print}' "$CFG_YAML" | sed 's/^[ \t]*//;s/[ \t]*$//'
)

readarray -t DATATYPES < <(
  awk '/^datatypes:/{f=1;next} /^[^ ]/{if(f)exit} 
       f {sub("- ",""); print}' "$CFG_YAML" | sed 's/^[ \t]*//;s/[ \t]*$//'
)
NC=${#CANCERS[@]}
ND=${#DATATYPES[@]}
TOTAL=$((NC * ND))

if [[ $SLURM_ARRAY_TASK_ID -ge $TOTAL ]]; then
  echo "[SKIP] Task index out of range."
  exit 0
fi

task=${SLURM_ARRAY_TASK_ID}
cancer_idx=$(( task / ND ))
dtype_idx=$(( task % ND ))

CANCER=${CANCERS[$cancer_idx]}
DTYPE=${DATATYPES[$dtype_idx]}

CANCER="TCGA_${CANCER}"

echo "[TASK] ArrayID=$SLURM_ARRAY_TASK_ID  CANCER=$CANCER  DTYPE=$DTYPE"
echo "---------------------------------------------------------"

# =====================================================================
# INPUT FILES
# =====================================================================
expr_norm="${EXPR_DIR}/${CANCER}_${DTYPE}.normalized.csv"
clin_file="01_transcriptomics/data/clinical/${CANCER}_clinical.csv"
mut_any_file="${MUT_DIR}/${CANCER}/gene/${CANCER}_gene_ensembl_coding_any.csv"

# =====================================================================
# CHECK INPUTS
# =====================================================================
if [[ ! -f "$expr_norm" ]]; then
  echo "[SKIP] Missing expression file: $expr_norm"
  exit 0
fi

if [[ ! -f "$clin_file" ]]; then
  echo "[SKIP] Missing clinical file: $clin_file"
  exit 0
fi

# Mutation file is optional — some cancers have 0 mutations
if [[ ! -f "$mut_any_file" ]]; then
  echo "[WARN] No mutation ANY file for $CANCER — Models 2, 3, 4 will skip."
fi

mkdir -p "$OUT_BASE/$CANCER"

# =====================================================================
# RUN
# =====================================================================
echo "[RUN] CoxPH models for $CANCER / $DTYPE"
echo

conda run -p "$R_ENV" --no-capture-output Rscript \
  01_transcriptomics/pipeline/scripts/04_run_cox_models.R \
  "$CANCER" \
  "$DTYPE" \
  "$TX2GENE" \
  "$OUT_BASE"

echo "[DONE] $CANCER / $DTYPE"

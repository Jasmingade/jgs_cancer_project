#!/usr/bin/env bash
#SBATCH --job-name=03_univariate_coxph
#SBATCH --cpus-per-task=6
#SBATCH --time=06:00:00
#SBATCH --mem=24G
#SBATCH --output=01_transcriptomics/logs/03_coxph/%x-%A.out
#SBATCH --error=01_transcriptomics/logs/03_coxph/%x-%A.err
#SBATCH --array=0-300

set -euo pipefail
mkdir -p 01_transcriptomics/logs/03_coxph 01_transcriptomics/out/03_univariate_coxph

# -------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}
CFG_YAML="01_transcriptomics/config/cancers.yaml"
COV_YAML="01_transcriptomics/config/covariates.yaml"
MUT_BASE="01_transcriptomics/out/02_mutation/02_combined_expression_mutation/split_by_cancer_type"

# Read cancers and datatypes from YAML
readarray -t CANCERS < <(awk '/^cancers:/{f=1;next}/^[^ ]/{if(f)exit}f' "$CFG_YAML" \
  | sed -n 's/^[[:space:]]*-[[:space:]]*//p')
readarray -t DATATYPES < <(awk '/^datatypes:/{f=1;next}/^[^ ]/{if(f)exit}f' "$CFG_YAML" \
  | sed -n 's/^[[:space:]]*-[[:space:]]*//p')

# Compute total tasks
NC=${#CANCERS[@]}
ND=${#DATATYPES[@]}
TOTAL=$((NC * ND))

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "[HINT] Submit as: sbatch --array=0-$((TOTAL-1)) $0"
  exit 0
fi

TASK=${SLURM_ARRAY_TASK_ID}
CTYPE_INDEX=$(( TASK / ND ))
DTYPE_INDEX=$(( TASK % ND ))

CANCER=${CANCERS[$CTYPE_INDEX]}
DTYPE=${DATATYPES[$DTYPE_INDEX]}
echo "------------------------------------------------------------"
echo "[TASK] Cancer=$CANCER | DType=$DTYPE"
echo "------------------------------------------------------------"

# Paths
expr_norm="01_transcriptomics/out/02_norm/TCGA_${CANCER}_${DTYPE}.normalized.csv"
mani="01_transcriptomics/out/02_norm/TCGA_${CANCER}_${DTYPE}.sample_manifest.csv"

if [[ ! -f "$expr_norm" || ! -f "$mani" ]]; then
  echo "[SKIP] Missing expression or manifest file for $CANCER / $DTYPE"
  exit 0
fi

# -------------------------------------------------------------------
# MODES: baseline / combined / interaction
# -------------------------------------------------------------------
for MODE in baseline combined interaction; do
  OUT_DIR="01_transcriptomics/out/03_univariate_coxph/${MODE}"
  mkdir -p "$OUT_DIR"

  out_res="${OUT_DIR}/TCGA_${CANCER}_${DTYPE}.cox_results.csv"
  out_sum="${OUT_DIR}/TCGA_${CANCER}_${DTYPE}.summary.txt"

  # Determine expected mutation matrix (for combined/interaction)
  mut_path="${MUT_BASE}/TCGA_${CANCER}_mutation.normalized.csv"
  if [[ "$MODE" != "baseline" && ! -f "$mut_path" ]]; then
    echo "[SKIP] No mutation data for $CANCER (needed for $MODE mode): $mut_path"
    continue
  fi

  echo "[RUN] Mode=${MODE} | Cancer=${CANCER} | DType=${DTYPE}"
  echo "  Expr : $expr_norm"
  echo "  Mani : $mani"
  echo "  Mode : $MODE"
  echo "  Out  : $out_res"

  { time conda run -p "$R_ENV" --no-capture-output \
      Rscript 01_transcriptomics/pipeline/scripts/03_univariate_coxph.R \
        "$expr_norm" "$mani" "$COV_YAML" "$out_res" "$out_sum" "$MODE"; } 2>&1

  echo "[DONE] ${CANCER}/${DTYPE} (${MODE}) → $out_res"
done

echo "------------------------------------------------------------"
echo "[ALL DONE] $CANCER / $DTYPE completed across available modes"

#!/bin/bash
#SBATCH --job-name=normalize_tcga
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --array=0-110
#SBATCH --output=01_transcriptomics/logs/02_norm/%x-%A_%a.log
#SBATCH --error=01_transcriptomics/logs/02_norm/%x-%A_%a.err

set -euo pipefail
mkdir -p 01_transcriptomics/logs/02_norm

echo "=== Normalizing all TCGA transcriptomics data (array) ==="
date; echo "PWD: $(pwd)"

DATATYPES=("gene" "iso_log" "iso_frac")
mapfile -t CANCERS < <(sed -n 's/^[[:space:]]*-[[:space:]]*//p' 01_transcriptomics/config/cancers.yaml)

ND=${#DATATYPES[@]} # number of data types
NC=${#CANCERS[@]}   # number of cancers
TOTAL=$(( NC * ND ))

TASK=${SLURM_ARRAY_TASK_ID:-0}
if (( TASK >= TOTAL )); then
  echo "Task index $TASK >= TOTAL $TOTAL → nothing to do"; exit 0
fi


if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "[HINT] sbatch --array=0-$((TOTAL-1)) $0"
  exit 0
fi

CIDX=$(( TASK / ND ))
DIDX=$(( TASK % ND ))
cancer="${CANCERS[$CIDX]}"
dtype="${DATATYPES[$DIDX]}"

expr_in="01_transcriptomics/data/${dtype}/RNA_${cancer}_${dtype}.csv"
clin_in="01_transcriptomics/data/clinical/TCGA_${cancer}_clinical.csv"
out_expr="01_transcriptomics/out/02_norm/TCGA_${cancer}_${dtype}.normalized.csv"
out_mani="01_transcriptomics/out/02_norm/TCGA_${cancer}_${dtype}.sample_manifest.csv"

echo "[Task ${SLURM_ARRAY_TASK_ID}] cancer=${cancer} dtype=${dtype}"
echo "expr_in: $expr_in"
echo "clin_in: $clin_in"

# Threading/env for the R script
export DT_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export RUN_PLOTS="${RUN_PLOTS:-TRUE}"
export MASTER_SUMMARY="01_transcriptomics/out/02_norm/master_summary.csv"
export TPM_CUTOFF="${TPM_CUTOFF:-1.0}"
export FRACTION_CUTOFF="${FRACTION_CUTOFF:-0.10}"

if [[ -f "$expr_in" && -f "$clin_in" ]]; then
  mkdir -p "$(dirname "$out_expr")" "$(dirname "$out_mani")"
  echo "[$(date '+%H:%M:%S')] Normalizing ${cancer} ${dtype}..."
  conda run -p /home/people/s184275/r-env --no-capture-output \
    Rscript 01_transcriptomics/pipeline/scripts/02_norm_log2_tpm.R \
    "$expr_in" "$clin_in" 01_transcriptomics/config/covariates.yaml "$out_expr" "$out_mani"
  echo "Wrote: $out_expr"
else
  echo "Skipping ${cancer} ${dtype} (missing files)"
  ls -lh "$expr_in" "$clin_in" 2>/dev/null || true
fi

date; echo "=== Task done ==="

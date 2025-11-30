#!/bin/bash
#SBATCH --job-name=tcga_coxph
#SBATCH --cpus-per-task=6
#SBATCH --time=06:00:00
#SBATCH --mem=24G
#SBATCH --output=01_transcriptomics/logs/03_coxph/%x-%A.out
#SBATCH --error=01_transcriptomics/logs/03_coxph/%x-%A.err
#SBATCH --array=0-98

set -euo pipefail
mkdir -p 01_transcriptomics/logs/03_coxph 01_transcriptomics/out/03_univariate_coxph

# Define sets
readarray -t CANCERS < <(sed -n 's/^[[:space:]]*-[[:space:]]*//p' 01_transcriptomics/config/cancers.yaml)

# Read cancers and datatypes from YAML (two lists)
readarray -t CANCERS < <(awk '/^cancers:/{f=1;next}/^[^ ]/{if(f)exit}f' 01_transcriptomics/config/cancers.yaml \
                         | sed -n 's/^[[:space:]]*-[[:space:]]*//p')

readarray -t DATATYPES < <(awk '/^datatypes:/{f=1;next}/^[^ ]/{if(f)exit}f' 01_transcriptomics/config/cancers.yaml \
                           | sed -n 's/^[[:space:]]*-[[:space:]]*//p')

# Calculate total tasks
NC=${#CANCERS[@]}
ND=${#DATATYPES[@]}
TOTAL=$((NC * ND))

# If not submitted as array yet, print how to submit
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "Submit as: sbatch --array=0-$((TOTAL-1)) $0"
  exit 0
fi

TASK=${SLURM_ARRAY_TASK_ID}
CTYPE_INDEX=$(( TASK / ND ))
DTYPE_INDEX=$(( TASK % ND ))

CANCER=${CANCERS[$CTYPE_INDEX]}
DTYPE=${DATATYPES[$DTYPE_INDEX]}

echo "Task $TASK -> CANCER=$CANCER, DTYPE=$DTYPE"

expr_norm="01_transcriptomics/out/02_norm_batch/TCGA_${CANCER}_${DTYPE}.normalized.csv"
mani="01_transcriptomics/out/02_norm_batch/TCGA_${CANCER}_${DTYPE}.sample_manifest.csv"
out_res="01_transcriptomics/out/03_univariate_coxph/TCGA_${CANCER}_${DTYPE}.cox_results.csv"
out_sum="01_transcriptomics/out/03_univariate_coxph/TCGA_${CANCER}_${DTYPE}.cox_summary.txt"

if [[ -f "$expr_norm" && -f "$mani" ]]; then
  echo "Running: $CANCER / $DTYPE with $SLURM_CPUS_PER_TASK cores"
  conda run -p /home/people/s184275/r-env --no-capture-output \
    Rscript 01_transcriptomics/pipeline/scripts/03_univariate_coxph.R \
      "$expr_norm" "$mani" 01_transcriptomics/config/covariates.yaml "$out_res" "$out_sum"
else
  echo "Missing inputs for $CANCER $DTYPE; skipping."
fi

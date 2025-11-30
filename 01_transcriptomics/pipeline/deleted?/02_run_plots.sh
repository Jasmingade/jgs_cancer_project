#!/bin/bash
#SBATCH --job-name=02_plots_norm
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --array=0-98
#SBATCH --output=01_transcriptomics/logs/02_plots/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/02_plots/%x-%A.err

set -euo pipefail
mkdir -p 01_transcriptomics/logs/02_plots 01_transcriptomics/out/02_plots

DATATYPES=("gene" "iso_log" "iso_frac")
mapfile -t CANCERS < <(sed -n 's/^[[:space:]]*-[[:space:]]*//p' 01_transcriptomics/config/cancers.yaml)

ND=${#DATATYPES[@]}
NC=${#CANCERS[@]}
TOTAL=$(( NC * ND ))

TASK=${SLURM_ARRAY_TASK_ID:-0}
(( TASK < TOTAL )) || { echo "TASK $TASK >= TOTAL $TOTAL"; exit 0; }

CIDX=$(( TASK / ND ))
DIDX=$(( TASK % ND ))
CANCER="${CANCERS[$CIDX]}"
DTYPE="${DATATYPES[$DIDX]}"

expr_in="01_transcriptomics/out/02_norm_batch/TCGA_${CANCER}_${DTYPE}.normalized.csv"
clin_in="01_transcriptomics/out/02_norm_batch/TCGA_${CANCER}_${DTYPE}.sample_manifest.csv"
cov_yaml="01_transcriptomics/config/covariates.yaml"
plot_dir="01_transcriptomics/out/02_plots"

echo "[Task $TASK] $CANCER $DTYPE (normalized input)"
if [[ -f "$expr_in" && -f "$clin_in" ]]; then
  export INPUT_IS_LOG2=TRUE
  export CANCER="$CANCER"
  export DTYPE="$DTYPE"
  conda run -p /home/people/s184275/r-env --no-capture-output \
    Rscript 01_transcriptomics/pipeline/scripts/02_norm_plots.R \
    "$expr_in" "$clin_in" "$cov_yaml" "$plot_dir"
else
  echo "Missing files -> skip"
  ls -lh "$expr_in" "$clin_in" 2>/dev/null || true
fi

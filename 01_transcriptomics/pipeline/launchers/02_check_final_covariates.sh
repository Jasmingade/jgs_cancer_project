#!/bin/bash
#SBATCH --job-name=check_covariates
#SBATCH --cpus-per-task=2
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --output=01_transcriptomics/logs/02_norm/check_covariates.log

echo "[INFO] Running manifest covariate completeness check after normalization array..."
conda run -p /home/people/s184275/r-env --no-capture-output \
  Rscript 01_transcriptomics/pipeline/scripts/02_check_final_covariates.R
echo "[DONE] Covariate completeness summary saved."
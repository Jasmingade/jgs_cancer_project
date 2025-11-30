#!/bin/bash
#SBATCH --job-name=model3_overlap
#SBATCH --output=01_transcriptomics/logs/overlap_plot/model3_overlap_%A.out
#SBATCH --error=01_transcriptomics/logs/overlap_plot/model3_overlap_%A.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

# ---------------------------------------------------
# Activate personal R environment
# ---------------------------------------------------
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}

echo "[INFO] Activating personal R environment: $R_ENV"
export R_LIBS_USER="$R_ENV"

# Ensure logs directory exists
mkdir -p logs

echo "[INFO] Starting Model 3 overlap plotting script"
echo "[INFO] Running on: $(hostname)"
echo "[INFO] Timestamp: $(date)"

# ---------------------------------------------------
# Run script
# ---------------------------------------------------
        conda run -p "$R_ENV" --no-capture-output \
            Rscript 01_transcriptomics/pipeline/plotting/05f_plot_overlap_sig_features.R

echo "[INFO] Finished at: $(date)"

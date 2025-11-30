#!/bin/bash
#SBATCH --job-name=collect_uni
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --output=01_transcriptomics/logs/04_collect/%x-%j.out
#SBATCH --error=01_transcriptomics/logs/04_collect/%x-%j.err

set -euo pipefail
mkdir -p 01_transcriptomics/logs/04_collect 01_transcriptomics/out/04_univariate_collect

conda run -p /home/people/s184275/r-env --no-capture-output \
  Rscript 01_transcriptomics/pipeline/scripts/04_collect_univariate_results.R

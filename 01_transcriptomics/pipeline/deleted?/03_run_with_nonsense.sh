#!/usr/bin/env bash
#SBATCH --job-name=03_nonsense_pipeline
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --array=0-999      # will be reset dynamically below
#SBATCH --output=01_transcriptomics/logs/03_nonsense/%x-%A_%a.log
#SBATCH --error=01_transcriptomics/logs/03_nonsense/%x-%A_%a.err

set -euo pipefail

# -----------------------------
# User-configurable paths
# -----------------------------
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}

MC3_GZ=${MC3_GZ:-"01_transcriptomics/data/mutation/mc3.v0.2.8.PUBLIC.xena.all_mutation_positions.gz"}
NONSENSE_TSV=${NONSENSE_TSV:-"01_transcriptomics/out/mutation/mc3_nonsense_only.tsv"}

CANCERS_YAML=${CANCERS_YAML:-"01_transcriptomics/config/cancers.yaml"}
DATATYPES_DEFAULT="gene iso_frac iso_log"
read -r -a DATATYPES <<< "${DATATYPES:-$DATATYPES_DEFAULT}"

PLOT_LOG_DIR="01_transcriptomics/logs/03_nonsense"
OUT_MUT_DIR=${OUT_MUT_DIR:-"01_transcriptomics/out/mutation"}
NORM_DIR=${NORM_DIR:-"01_transcriptomics/out/02_norm_batch"}
COX_OUT_DIR=${COX_OUT_DIR:-"01_transcriptomics/out/03_univariate_coxph"}

SCRIPT_02A=${SCRIPT_02A:-"01_transcriptomics/pipeline/scripts/02a_extract_nonsense_mutations.R"}
SCRIPT_02B=${SCRIPT_02B:-"01_transcriptomics/pipeline/scripts/02b_build_nonsense_covariates.R"}
SCRIPT_03=${SCRIPT_03:-"01_transcriptomics/pipeline/scripts/03_univariate_coxph.R"}
COV_YAML=${COV_YAML:-"01_transcriptomics/config/covariates.yaml"}

mkdir -p "$PLOT_LOG_DIR" "$OUT_MUT_DIR" "$COX_OUT_DIR"

# -----------------------------
# Parse cancers from YAML
# -----------------------------
mapfile -t CANCERS < <(sed -n 's/^[[:space:]]*-[[:space:]]*//p' "$CANCERS_YAML")
NC=${#CANCERS[@]}
ND=${#DATATYPES[@]}
TOTAL=$(( NC * ND ))           # tasks > 0 map to cancer×dtype
# We reserve TASK 0 for building the nonsense subset once.
# Update array size dynamically on first task (harmless if already set correctly)
if [[ -n "${SLURM_ARRAY_TASK_MAX:-}" && "$SLURM_ARRAY_TASK_MAX" -gt 1 ]]; then
  : # SLURM already set it
else
  echo "[INFO] Suggested array: 0-$TOTAL"
fi

TASK=${SLURM_ARRAY_TASK_ID:-0}

# -----------------------------
# Helper: build nonsense subset ONCE (with a lock)
# -----------------------------
build_nonsense_once() {
  local lock="${NONSENSE_TSV}.lockdir"
  if [[ -f "$NONSENSE_TSV" ]]; then
    echo "[NONSENSE] Found: $NONSENSE_TSV"
    return 0
  fi
  if mkdir "$lock" 2>/dev/null; then
    echo "[NONSENSE] Building nonsense-only subset → $NONSENSE_TSV"
    conda run -p "$R_ENV" --no-capture-output \
      Rscript "$SCRIPT_02A" "$MC3_GZ" "$NONSENSE_TSV"
    rmdir "$lock" || true
  else
    echo "[NONSENSE] Another task is building it; waiting..."
    # Wait up to ~30 minutes; adjust if needed
    for _ in {1..180}; do
      [[ -f "$NONSENSE_TSV" ]] && break
      sleep 10
    done
    [[ -f "$NONSENSE_TSV" ]] || { echo "[ERROR] Timeout waiting for $NONSENSE_TSV"; exit 1; }
  fi
}

# -----------------------------
# TASK 0: build nonsense subset and exit
# -----------------------------
if (( TASK == 0 )); then
  build_nonsense_once
  echo "[DONE] Task 0 complete."
  exit 0
fi

# -----------------------------
# Map TASK>0 → (cancer, dtype)
# -----------------------------
t=$(( TASK - 1 ))
CIDX=$(( t / ND ))
DIDX=$(( t % ND ))
if (( CIDX >= NC )); then
  echo "[INFO] Task $TASK out of range (NC=$NC, ND=$ND). Exit."
  exit 0
fi
CANCER="${CANCERS[$CIDX]}"
DTYPE="${DATATYPES[$DIDX]}"

echo "[TASK $TASK] Cancer=$CANCER DType=$DTYPE"

# Ensure nonsense subset exists (build if missing)
build_nonsense_once

# -----------------------------
# Build per-cancer nonsense covariates (if missing)
# -----------------------------
COVAR_NONS="${OUT_MUT_DIR}/TCGA_${CANCER}_nonsense_covariates.csv"
if [[ ! -f "$COVAR_NONS" ]]; then
  # Prefer the gene manifest for case_id universe; fall back to the current dtype if needed
  MANI_GENE="${NORM_DIR}/TCGA_${CANCER}_gene.sample_manifest.csv"
  MANI_FALLBACK="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.sample_manifest.csv"

  if [[ -f "$MANI_GENE" ]]; then
    MANI_USE="$MANI_GENE"
  elif [[ -f "$MANI_FALLBACK" ]]; then
    MANI_USE="$MANI_FALLBACK"
  else
    echo "[WARN] No manifest found for $CANCER; looked for:"
    echo "       $MANI_GENE"
    echo "       $MANI_FALLBACK"
    exit 0
  fi

  echo "[COVAR] Building nonsense covariates for $CANCER using $MANI_USE"
  conda run -p "$R_ENV" --no-capture-output \
    Rscript "$SCRIPT_02B" \
      "$NONSENSE_TSV" \
      "$MANI_USE" \
      "$CANCER" \
      "$COVAR_NONS" \
      10 0.05
else
  echo "[COVAR] Exists: $COVAR_NONS"
fi

# Symlink to the filename your 03_univariate_coxph.R auto-merge expects
COVAR_LINK="${OUT_MUT_DIR}/TCGA_${CANCER}_mutation_covariates.csv"
if [[ ! -e "$COVAR_LINK" ]]; then
  ln -s "$(basename "$COVAR_NONS")" "$COVAR_LINK" || cp -f "$COVAR_NONS" "$COVAR_LINK"
  echo "[COVAR] Linked $COVAR_LINK → $(basename "$COVAR_NONS")"
fi

# -----------------------------
# Run 03_univariate_coxph.R for this cancer × dtype
# -----------------------------
EXPR_IN="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.normalized.csv"
MANI_IN="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.sample_manifest.csv"

if [[ ! -f "$EXPR_IN" || ! -f "$MANI_IN" ]]; then
  echo "[SKIP] Missing inputs:"
  ls -lh "$EXPR_IN" "$MANI_IN" 2>/dev/null || true
  exit 0
fi

OUT_RES="${COX_OUT_DIR}/TCGA_${CANCER}_${DTYPE}.cox_results.csv"
OUT_SUM="${COX_OUT_DIR}/TCGA_${CANCER}_${DTYPE}.summary.txt"

echo "[RUN] 03_univariate_coxph.R for $CANCER $DTYPE"
conda run -p "$R_ENV" --no-capture-output \
  Rscript "$SCRIPT_03" \
    "$EXPR_IN" \
    "$MANI_IN" \
    "$COV_YAML" \
    "$OUT_RES" \
    "$OUT_SUM"

echo "[DONE] $CANCER $DTYPE → $OUT_RES"

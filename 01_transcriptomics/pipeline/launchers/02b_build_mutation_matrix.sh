#!/usr/bin/env bash
#SBATCH --job-name=02_mut_matrix
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --mem=36G
#SBATCH --output=01_transcriptomics/logs/02_mutation/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/02_mutation/%x-%A.err
set -euo pipefail

# ============================================================
# CONFIGURATION
# ============================================================
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}
MC3_GZ=${MC3_GZ:-"01_transcriptomics/data/mutation/mc3.v0.2.8.PUBLIC.xena.all_mutation_positions.gz"}

OUT_DIR=${OUT_DIR:-"01_transcriptomics/out/02_mutation"}
GLOBAL_CSV="${OUT_DIR}/mutation_covariates.global.csv"
NORM_DIR=${NORM_DIR:-"01_transcriptomics/out/02_norm"}
LOG_DIR=${LOG_DIR:-"01_transcriptomics/logs/02_mutation"}

SCRIPT_02A=${SCRIPT_02A:-"01_transcriptomics/pipeline/scripts/02a_all_mutations.R"}
SCRIPT_02B=${SCRIPT_02B:-"01_transcriptomics/pipeline/scripts/02b_build_mutation_covariates.R"}
SCRIPT_02C=${SCRIPT_02C:-"01_transcriptomics/pipeline/scripts/02c_split_mutation_types.R"}
SCRIPT_02D=${SCRIPT_02D:-"01_transcriptomics/pipeline/scripts/02d_summarize_mutations_by_cancer.R"}

ALL_MUT_TSV="${OUT_DIR}/mc3_all_mutations.tsv"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$NORM_DIR"

# -------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------
MIN_PREV=${MIN_PREV:-0.05}       # Default mutation prevalence cutoff
ID_TYPE=${ID_TYPE:-ensembl}      # Default gene ID type (ensembl|hgnc)

# Environmental variable
ID_TYPE=${ID_TYPE:-ensembl}

{ time conda run -p "$R_ENV" --no-capture-output \
    Rscript "$SCRIPT_02B" \
      "$ALL_MUT_TSV" \
      "$OUT_DIR" \
      "$MIN_PREV" \
      "$ID_TYPE"; } 2>&1

# Color setup
C_RESET="\033[0m"
C_BLUE="\033[1;34m"
C_GREEN="\033[1;32m"
C_YELLOW="\033[1;33m"
C_RED="\033[1;31m"

ts() { date "+%Y-%m-%d %H:%M:%S"; }

echo -e "${C_BLUE}============================================================${C_RESET}"
echo -e "${C_BLUE}[START] Mutation Matrix Build Job${C_RESET}"
date
echo "------------------------------------------------------------"
echo "  Input MC3 file : $MC3_GZ"
echo "  Output folder  : $OUT_DIR"
echo "  R environment  : $R_ENV"
echo "  Output TSV     : $ALL_MUT_TSV"
echo "  Output CSV     : $GLOBAL_CSV"
echo "  Gene ID mode   : $ID_TYPE"
echo "------------------------------------------------------------"

# ============================================================
# STEP A: Extract clean mutation table
# ============================================================
echo -e "${C_YELLOW}[STEP A] $(ts) - Extract clean mutation table${C_RESET}"
if [[ -f "$ALL_MUT_TSV" ]]; then
  echo "[SKIP] Existing: $ALL_MUT_TSV"
else
  echo "[RUN] Rscript $SCRIPT_02A"
  { time conda run -p "$R_ENV" --no-capture-output \
      Rscript "$SCRIPT_02A" "$MC3_GZ" "$ALL_MUT_TSV"; } 2>&1
  echo "[OK] $(ts) - Wrote $ALL_MUT_TSV"
fi
echo

# ============================================================
# STEP B: Build per-gene mutation matrices (3 analyses)
# ============================================================
echo -e "${C_YELLOW}[STEP B] $(ts) - Build GLOBAL mutation covariates (all cases)${C_RESET}"

ID_TYPE=${ID_TYPE:-ensembl}   # or set hgnc for HGNC naming
echo "[CFG] ID_TYPE = $ID_TYPE"

if [[ -f "$GLOBAL_CSV" ]]; then
  echo "[SKIP] Existing: $GLOBAL_CSV"
else
  echo "[RUN] Rscript $SCRIPT_02B"
  { time conda run -p "$R_ENV" --no-capture-output \
      Rscript "$SCRIPT_02B" \
        "$ALL_MUT_TSV" \
        "$OUT_DIR" \
        "$MIN_PREV" \
        "$ID_TYPE"; } 2>&1
  echo "[OK] $(ts) - Wrote per-gene mutation matrices"
fi
echo
# ============================================================
# STEP C: Split per-cancer matrices (baseline/combined/interaction)
# ============================================================
for ANALYSIS in 01_per_gene_baseline 02_combined_expression_mutation 03_interaction_models; do
  case "$ANALYSIS" in
    01_per_gene_baseline) MATRIX="$OUT_DIR/$ANALYSIS/mutation_covariates.by_gene.csv" ;;
    02_combined_expression_mutation) MATRIX="$OUT_DIR/$ANALYSIS/mutation_covariates.combined.csv" ;;
    03_interaction_models) MATRIX="$OUT_DIR/$ANALYSIS/mutation_covariates.interaction.csv" ;;
  esac
  OUT_SUBDIR="$OUT_DIR/${ANALYSIS}/split_by_cancer_type"
  
  if [[ -d "$OUT_SUBDIR" && -n "$(find "$OUT_SUBDIR" -name '*.normalized.csv' -type f 2>/dev/null)" ]]; then
    echo "[SKIP] Existing split files for ${ANALYSIS}"
    continue
  fi

  if [[ -f "$MATRIX" ]]; then
    echo "[RUN] Splitting mutation matrix for ${ANALYSIS}"
    { time conda run -p "$R_ENV" --no-capture-output \
        Rscript "$SCRIPT_02C" "$MATRIX" "$NORM_DIR" "$OUT_DIR/$ANALYSIS"; } 2>&1
    echo "[OK] $(ts) - Split matrices written under $OUT_DIR/$ANALYSIS/split_by_cancer_type"
  else
    echo "[WARN] Missing expected matrix: $MATRIX"
  fi
done


# ============================================================
# STEP D: Summarize per-cancer mutation prevalence
# ============================================================
echo -e "${C_YELLOW}[STEP D] $(ts) - Summarize mutation prevalence by cancer type${C_RESET}"
OUT_SUM_LONG="$OUT_DIR/mutation_gene_prevalence_by_cancer.long.csv"
if [[ -f "$OUT_SUM_LONG" ]]; then
  echo "[SKIP] Existing: $OUT_SUM_LONG"
else
  { time conda run -p "$R_ENV" --no-capture-output \
      Rscript "$SCRIPT_02D" "$OUT_DIR/01_per_gene_baseline/mutation_covariates.by_gene.csv" "$NORM_DIR" "$OUT_DIR"; } 2>&1
  echo "[OK] $(ts) - Wrote mutation prevalence summaries"
fi
echo

# ============================================================
# STEP E: Overview of outputs
# ============================================================
echo -e "${C_BLUE}------------------------------------------------------------${C_RESET}"
echo "[SUMMARY] Mutation outputs available in:"
echo "  - Baseline:     $OUT_DIR/01_per_gene_baseline/"
echo "  - Combined:     $OUT_DIR/02_combined_expression_mutation/"
echo "  - Interaction:  $OUT_DIR/03_interaction_models/"
echo
echo "[INFO] Example files:"
find "$OUT_DIR" -type f -name "mutation_covariates*.csv" | sed 's/^/  - /'
echo
echo "[DONE] $(ts) - Mutation data ready for downstream CoxPH analyses"
echo -e "${C_BLUE}============================================================${C_RESET}"

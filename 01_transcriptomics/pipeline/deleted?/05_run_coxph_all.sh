#!/usr/bin/env bash
#SBATCH --job-name=05_cox_mut
#SBATCH --cpus-per-task=2
#SBATCH --time=06:00:00
#SBATCH --mem=12G
#SBATCH --array=0-999
#SBATCH --output=01_transcriptomics/logs/05_cox_mut/%x-%A_%a.log
#SBATCH --error=01_transcriptomics/logs/05_cox_mut/%x-%A_%a.err
set -euo pipefail

# -------------------- CONFIG --------------------
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}

CANCERS_YAML=${CANCERS_YAML:-"01_transcriptomics/config/cancers.yaml"}
COV_YAML=${COV_YAML:-"01_transcriptomics/config/covariates.yaml"}

# Inputs produced earlier
NORM_DIR=${NORM_DIR:-"01_transcriptomics/out/02_norm_batch"}
MUT_GLOBAL=${MUT_GLOBAL:-"01_transcriptomics/out/04_mutation/mutation_covariates.global.csv"}

# Outputs
COX_DIR=${COX_DIR:-"01_transcriptomics/out/05_univariate_coxph_mut"}
LOG_DIR=${LOG_DIR:-"01_transcriptomics/logs/05_cox_mut"}
mkdir -p "$COX_DIR" "$LOG_DIR"

# Scripts
SCRIPT_COX=${SCRIPT_COX:-"01_transcriptomics/pipeline/scripts/03_univariate_coxph.R"}

# Which mutation “data types” you declared in cancers.yaml
# Map datatype -> column name in MUT_GLOBAL
declare -A MUT_MAP=(
  [mutation_nonsense]="nonsense"
  [mutation_missense]="missense"
  [mutation_frameshift]="frameshift"
  [mutation_splice]="splice"
)

# -------------------- PARSE YAML --------------------
# robust parse of lists (works with Windows line endings)
readarray -t CANCERS   < <(tr -d '\r' < "$CANCERS_YAML" | awk '/^cancers:/{flag=1;next}/^[^ ]/{if(flag)exit}flag' | sed -n 's/^[[:space:]]*-[[:space:]]*//p')
readarray -t DATATYPES < <(tr -d '\r' < "$CANCERS_YAML" | awk '/^datatypes:/{flag=1;next}/^[^ ]/{if(flag)exit}flag' | sed -n 's/^[[:space:]]*-[[:space:]]*//p')

NC=${#CANCERS[@]}
ND=${#DATATYPES[@]}
TOTAL=$((NC*ND))
TASK="${SLURM_ARRAY_TASK_ID:-0}"

if [[ -z "${SLURM_ARRAY_TASK_MAX:-}" || "${SLURM_ARRAY_TASK_MAX:-0}" -lt $TOTAL ]]; then
  echo "[INFO] Found $NC cancers × $ND datatypes = $TOTAL tasks. Suggested --array=0-$((TOTAL-1))."
fi

# Map array index -> (cancer, datatype)
CIDX=$(( TASK / ND ))
DIDX=$(( TASK % ND ))
if (( CIDX >= NC )); then
  echo "[INFO] Task $TASK out of range (NC=$NC, ND=$ND). Exit."; exit 0
fi
CANCER="${CANCERS[$CIDX]}"
DTYPE="${DATATYPES[$DIDX]}"

echo "------------------------------------------------------------"
echo "[TASK $TASK] Cancer=$CANCER  DTYPE=$DTYPE"
echo "NORM_DIR=$NORM_DIR"
echo "COX_DIR =$COX_DIR"
echo "------------------------------------------------------------"

# -------------------- RESOLVE INPUTS --------------------
# For expression/isoform datatypes: expect prebuilt normalized matrices
if [[ "$DTYPE" != mutation_* ]]; then
  EXPR_IN="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.normalized.csv"
  MANI_IN="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.sample_manifest.csv"
  if [[ ! -f "$EXPR_IN" || ! -f "$MANI_IN" ]]; then
    echo "[SKIP] Missing ${EXPR_IN} or ${MANI_IN}"; exit 0
  fi
else
  # For mutation datatypes: auto-materialize a binary matrix from MUT_GLOBAL (if absent)
  col="${MUT_MAP[$DTYPE]:-}"
  if [[ -z "$col" ]]; then
    echo "[ERROR] No mapping for mutation datatype '$DTYPE' in MUT_MAP."; exit 1
  fi
  if [[ ! -f "$MUT_GLOBAL" ]]; then
    echo "[ERROR] Missing global mutation covariates: $MUT_GLOBAL"; exit 1
  fi
  # Build per-cancer matrix if missing
  EXPR_IN="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.normalized.csv"
  MANI_IN="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.sample_manifest.csv"
  if [[ ! -f "$EXPR_IN" || ! -f "$MANI_IN" ]]; then
    echo "[BUILD] Creating matrix for $DTYPE (col='$col') for $CANCER from $(basename "$MUT_GLOBAL")"
    MANI_GENE="${NORM_DIR}/TCGA_${CANCER}_gene.sample_manifest.csv"
    if [[ ! -f "$MANI_GENE" ]]; then
      for dt in iso_log iso_frac; do
        cand="${NORM_DIR}/TCGA_${CANCER}_${dt}.sample_manifest.csv"
        [[ -f "$cand" ]] && MANI_GENE="$cand" && break
      done
    fi
    if [[ ! -f "$MANI_GENE" ]]; then
      echo "[SKIP] No manifest found for $CANCER in $NORM_DIR"; exit 0
    fi
    # Use awk to join by case_id and output 0/1 column
    # (Assumes headers present; robust to column order)
    tmp_cases=$(mktemp)
    tmp_cols=$(mktemp)
    awk -F, 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} {print $h["case_id"]}' "$MANI_GENE" > "$tmp_cases"
    # Extract header to find mutation column index
    IFS=, read -r -a hdr < <(head -1 "$MUT_GLOBAL")
    idx=-1
    for i in "${!hdr[@]}"; do [[ "${hdr[$i]}" == "$col" ]] && idx=$((i+1)) && break; done
    if (( idx < 0 )); then echo "[ERROR] Column '$col' not found in $(basename "$MUT_GLOBAL")"; exit 1; fi
    # Build case_id,value CSV
    {
      echo "case_id,${DTYPE}"
      # Create a lookup of case_id -> value from MUT_GLOBAL
      awk -F, -v ci="$idx" '
        NR==1{
          for(i=1;i<=NF;i++)h[$i]=i
          cid=h["case_id"]
          next
        }
        { printf "%s,%s\n", $cid, ($ci==""?0:$ci) }
      ' "$MUT_GLOBAL" \
      | sort -t, -k1,1 \
      | join -t, -1 1 -2 1 -o 1.1,2.2 - <(sort -u "$tmp_cases") \
      | awk -F, 'BEGIN{OFS=","} { if(NF==1) print $1,0; else print }'
    } > "$EXPR_IN"
    cp -f "$MANI_GENE" "$MANI_IN"
    rm -f "$tmp_cases" "$tmp_cols"
    echo "[OK] Wrote ${EXPR_IN} and ${MANI_IN}"
  else
    echo "[FOUND] ${EXPR_IN}"
  fi
fi

# -------------------- RUN COXPH --------------------
OUT_RES="${COX_DIR}/TCGA_${CANCER}_${DTYPE}.cox_results.csv"
OUT_SUM="${COX_DIR}/TCGA_${CANCER}_${DTYPE}.summary.txt"

echo "[RUN] CoxPH: ${CANCER} × ${DTYPE}"
conda run -p "$R_ENV" --no-capture-output \
  Rscript "$SCRIPT_COX" \
    "$EXPR_IN" \
    "$MANI_IN" \
    "$COV_YAML" \
    "$OUT_RES" \
    "$OUT_SUM"

echo "[DONE] ${CANCER} ${DTYPE} → $(basename "$OUT_RES")"

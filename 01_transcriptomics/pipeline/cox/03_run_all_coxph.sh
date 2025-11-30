#!/usr/bin/env bash
#SBATCH --job-name=03_run_all_coxph
#SBATCH --cpus-per-task=6
#SBATCH --time=06:00:00
#SBATCH --output=01_transcriptomics/logs/03_coxph/global/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/03_coxph/global/%x-%A.err
#SBATCH --mem=24G
#SBATCH --array=0-132
set -euo pipefail

# ============================================================
# CONFIGURATION
# ============================================================
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}
CFG_YAML="01_transcriptomics/config/cancers.yaml"
COV_YAML="01_transcriptomics/config/covariates.yaml"

BASE_OUT="01_transcriptomics/out"
BASE_LOG="01_transcriptomics/logs/03_coxph"
NORM_DIR="${BASE_OUT}/02_norm"
MUT_DIR="${BASE_OUT}/03_mutation"

SCRIPT_EXPR="01_transcriptomics/pipeline/cox/03a_univariate_coxph.R"
SCRIPT_MUT="01_transcriptomics/pipeline/cox/03b_mutation_univariate_coxph.R"
SCRIPT_EXP_MUT="01_transcriptomics/pipeline/cox/03c_exp_mutation_univariate_coxph.R"
SCRIPT_ISO_MUT="01_transcriptomics/pipeline/cox/03d_iso_mut_univariate_coxph.R"

# If FORCE_RERUN=true → always rerun.
# If FORCE_RERUN=false → skip when output exists.
FORCE_RERUN=${FORCE_RERUN:-false}

# Model-specific log dirs
LOG_EXPR="${BASE_LOG}/03a_expr"
LOG_MUT="${BASE_LOG}/03b_mut"
LOG_EXP_MUT="${BASE_LOG}/03c_exp_mut"
LOG_ISO_MUT="${BASE_LOG}/03d_iso_mut"

mkdir -p "$BASE_LOG" \
         "$LOG_EXPR" "$LOG_MUT" "$LOG_EXP_MUT" "$LOG_ISO_MUT" \
         "${BASE_OUT}/03a_univariate_coxph" \
         "${BASE_OUT}/03b_mutation_univariate_coxph" \
         "${BASE_OUT}/03c_exp_mutation_univariate_coxph" \
         "${BASE_OUT}/03d_iso_mut_univariate_coxph"

# ============================================================
# EXECUTION CONTROL
# ============================================================
RUN_01=${RUN_01:-false}  # Model 2 ------ Mutation-per-gene (per mutation group)
RUN_02=${RUN_02:-false}  # Model 3 ------ Model Expression + Mutation combined
RUN_03=${RUN_03:-true}   # Model 4 ------ Isoform × Mutation interaction
RUN_04=${RUN_04:-false}  # Model 1 ------ Expression-only

# Mutation groups
MUT_GROUPS=(
  "truncating_LOF"
  "missense_or_inframe"
  "rna"
  "splice"
)

# ============================================================
# PARSE YAML CONFIG
# ============================================================
readarray -t CANCERS < <(
  awk '/^cancers:/{f=1;next}/^[^ ]/{if(f)exit}f' "$CFG_YAML" \
  | sed -n 's/^[[:space:]]*-[[:space:]]*//p'
)

readarray -t DATATYPES < <(
  awk '/^datatypes:/{f=1;next}/^[^ ]/{if(f)exit}f' "$CFG_YAML" \
  | sed -n 's/^[[:space:]]*-[[:space:]]*//p'
)
NC=${#CANCERS[@]}
ND=${#DATATYPES[@]}
TOTAL=$((NC * ND))

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "[HINT] sbatch --array=0-$((TOTAL-1)) $0"
    exit 0
fi

TASK=${SLURM_ARRAY_TASK_ID}
CTYPE_INDEX=$(( TASK / ND ))
DTYPE_INDEX=$(( TASK % ND ))
CANCER=${CANCERS[$CTYPE_INDEX]}
DTYPE=${DATATYPES[$DTYPE_INDEX]}

echo "[TASK] CANCER=$CANCER | DTYPE=$DTYPE"
echo "------------------------------------------------------------"

# ============================================================
# PATHS
# ============================================================
expr_norm="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.normalized.csv"
mani="${NORM_DIR}/TCGA_${CANCER}_${DTYPE}.sample_manifest.csv"

if [[ ! -f "$expr_norm" || ! -f "$mani" ]]; then
    echo "[SKIP] Missing inputs for $CANCER / $DTYPE"
    exit 0
fi

# Helper for FORCE_RERUN logic
should_rerun() {
  local out_csv="$1"
  if [[ "$FORCE_RERUN" == "true" ]]; then
    return 0  # always rerun
  fi
  if [[ -f "$out_csv" ]]; then
    return 1  # skip
  fi
  return 0      # no file → run
}

# ============================================================
# ① Expression-only baseline (03a)
# ============================================================
if [[ "${RUN_04}" == "true" ]]; then
    expr_out="${BASE_OUT}/03a_univariate_coxph/${DTYPE}"
    mkdir -p "$expr_out"

    out_csv="${expr_out}/TCGA_${CANCER}_${DTYPE}.cox_results.csv"
    out_sum="${expr_out}/TCGA_${CANCER}_${DTYPE}.summary.txt"
    log_file="${LOG_EXPR}/TCGA_${CANCER}_${DTYPE}.log"

    if ! should_rerun "$out_csv"; then
        echo "[SKIP] 03a Expression-only already done: $out_csv"
    else
        echo "[RUN] 03a Expression-only (CANCER=$CANCER, DTYPE=$DTYPE)"
        conda run -p "$R_ENV" --no-capture-output \
            Rscript "$SCRIPT_EXPR" \
                "$expr_norm" "$mani" "$COV_YAML" "$out_csv" "$out_sum" \
            2>&1 | tee "$log_file"
    fi
else
    echo "[SKIP] RUN_04=false → skipping 03a Expression-only"
fi
echo

# ============================================================
# ② Mutation-only CoxPH PER MUTATION GROUP (03b)
#    NOTE: run once per cancer → only when DTYPE == gene
# ============================================================
if [[ "${RUN_01}" == "true" ]]; then
    if [[ "$DTYPE" != "gene" ]]; then
      echo "[SKIP] Mutation-only: run only for DTYPE=gene (current=$DTYPE)"
    else
      out_dir="${BASE_OUT}/03b_mutation_univariate_coxph/TCGA_${CANCER}"
      mkdir -p "$out_dir"

      for grp in "${MUT_GROUPS[@]}"; do
          mut_in="${MUT_DIR}/TCGA_${CANCER}/gene/TCGA_${CANCER}_gene_ensembl_${grp}.csv"
          out_csv="${out_dir}/TCGA_${CANCER}_mutation_${grp}.cox_results.csv"
          out_sum="${out_dir}/TCGA_${CANCER}_mutation_${grp}.summary.txt"
          log_file="${LOG_MUT}/TCGA_${CANCER}_${grp}.log"

          if [[ ! -f "$mut_in" ]]; then
              echo "[SKIP] No mutation matrix for $CANCER group=$grp"
              continue
          fi

          if ! should_rerun "$out_csv"; then
              echo "[SKIP] Mutation-only (grp=$grp) already done: $out_csv"
              continue
          fi

          echo "[RUN] 03b Mutation-only (CANCER=$CANCER, grp=$grp)"
          conda run -p "$R_ENV" --no-capture-output \
            Rscript "$SCRIPT_MUT" "$mut_in" "$mani" "$COV_YAML" "$out_csv" "$out_sum" \
            2>&1 | tee "$log_file"
      done
    fi
else
    echo "[SKIP] RUN_01=false → skipping 03b mutation-only"
fi
echo

# ============================================================
# ③ Expression + Mutation (combined, 03c)
# ============================================================
if [[ "${RUN_02}" == "true" ]]; then
    out_dir2="${BASE_OUT}/03c_exp_mutation_univariate_coxph/TCGA_${CANCER}"
    mkdir -p "$out_dir2"

    for grp in "${MUT_GROUPS[@]}"; do
        mut_in="${MUT_DIR}/TCGA_${CANCER}/gene/TCGA_${CANCER}_gene_ensembl_${grp}.csv"
        out_csv="${out_dir2}/TCGA_${CANCER}_${DTYPE}_${grp}.cox_results.csv"
        out_sum="${out_dir2}/TCGA_${CANCER}_${DTYPE}_${grp}.summary.txt"
        log_file="${LOG_EXP_MUT}/TCGA_${CANCER}_${DTYPE}_${grp}.log"

        if [[ ! -f "$mut_in" ]]; then
            echo "[SKIP] No mutation matrix for combined model group=$grp"
            continue
        fi

        if ! should_rerun "$out_csv"; then
            echo "[SKIP] 03c Expression+Mutation (grp=$grp) already done: $out_csv"
            continue
        fi

        echo "[RUN] 03c Expression + mutation (CANCER=$CANCER, DTYPE=$DTYPE, grp=$grp)"
        conda run -p "$R_ENV" --no-capture-output \
          Rscript "$SCRIPT_EXP_MUT" \
              "$expr_norm" "$mut_in" "$mani" "$COV_YAML" "$out_csv" "$out_sum" \
          2>&1 | tee "$log_file"
    done
else
    echo "[SKIP] RUN_02=false → skipping 03c Expression+Mutation"
fi
echo

# ============================================================
# ④ Isoform × Mutation Interaction (03d)
#    (likely only meaningful for iso_log / iso_frac)
# ============================================================
if [[ "${RUN_03}" == "true" ]]; then
      tx2gene="01_transcriptomics/data/raw/tx2gene.csv"
      out_dir3="${BASE_OUT}/03d_iso_mut_univariate_coxph/TCGA_${CANCER}"
      mkdir -p "$out_dir3"

      for grp in "${MUT_GROUPS[@]}"; do
          mut_in="${MUT_DIR}/TCGA_${CANCER}/gene/TCGA_${CANCER}_gene_ensembl_${grp}.csv"
          out_csv="${out_dir3}/TCGA_${CANCER}_${DTYPE}_${grp}.cox_results.csv"
          out_sum="${out_dir3}/TCGA_${CANCER}_${DTYPE}_${grp}.summary.txt"
          log_file="${LOG_ISO_MUT}/TCGA_${CANCER}_${DTYPE}_${grp}.log"

          if [[ ! -f "$mut_in" ]]; then
              echo "[SKIP] No mutation matrix for interaction model group=$grp"
              continue
          fi

          if ! should_rerun "$out_csv"; then
              echo "[SKIP] 03d Isoform×mutation (grp=$grp) already done: $out_csv"
              continue
          fi

          echo "[RUN] 03d Isoform × mutation (CANCER=$CANCER, DTYPE=$DTYPE, grp=$grp)"
          conda run -p "$R_ENV" --no-capture-output \
            Rscript "$SCRIPT_ISO_MUT" \
                "$expr_norm" "$mut_in" "$mani" "$COV_YAML" "$tx2gene" "$out_csv" "$out_sum" \
            2>&1 | tee "$log_file"
      done
    fi
else
    echo "[SKIP] RUN_03=false → skipping 03d Isoform×Mutation"
fi


echo "[DONE] Completed all requested models for $CANCER / $DTYPE"

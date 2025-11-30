#!/usr/bin/env bash
#SBATCH --job-name=05_run_summaries
#SBATCH --cpus-per-task=6
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --output=01_transcriptomics/logs/05_summaries/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/05_summaries/%x-%A.err

set -euo pipefail

# ============================================================
# CONFIG
# ============================================================
R_ENV=${R_ENV:-"/home/people/s184275/r-env"}
CFG_YAML="01_transcriptomics/config/cancers.yaml"
COV_YAML="01_transcriptomics/config/covariates.yaml"

ROOT_OUT="01_transcriptomics/out"
NORM_DIR="${ROOT_OUT}/02_norm"
MUT_DIR="${ROOT_OUT}/03_mutation"
OUT_03E="${ROOT_OUT}/03e_gene_expr_mut"
GENE_LIST_DIR="${ROOT_OUT}/03e_gene_lists"
LOGDIR="01_transcriptomics/logs/05_summaries"

mkdir -p "$LOGDIR" "$OUT_03E"

# Mutation groups (keep in sync with other launchers)
MUT_GROUPS=(
  "truncating_LOF"
  "missense_or_inframe"
  "rna"
  "splice"
)

# Toggles
RUN_PREP_GENELISTS=${RUN_PREP_GENELISTS:-true}
RUN_GENE_EXPR_MUT=${RUN_GENE_EXPR_MUT:-true}
RUN_COLLECT_GLOBAL=${RUN_COLLECT_GLOBAL:-true}
RUN_BUILD_MODALITY=${RUN_BUILD_MODALITY:-true}
RUN_MODALITY_SUMMARY=${RUN_MODALITY_SUMMARY:-true}
FORCE_RERUN=${FORCE_RERUN:-false}

# ============================================================
# HELPERS
# ============================================================
run_r() {
  local script="$1"
  local log_file="$2"
  echo "------------------------------------------------------------"
  echo "[RUN] $script"
  conda run -p "$R_ENV" --no-capture-output Rscript "$script" 2>&1 | tee "$log_file"
}

should_rerun() {
  local out_csv="$1"
  if [[ "$FORCE_RERUN" == "true" ]]; then
    return 0
  fi
  [[ ! -f "$out_csv" ]]
}

# ------------------------------------------------------------
# Parse cancers from config
# ------------------------------------------------------------
readarray -t CANCERS < <(
  awk '/^cancers:/{f=1;next}/^[^ ]/{if(f)exit}f' "$CFG_YAML" \
  | sed -n 's/^[[:space:]]*-[[:space:]]*//p'
)

echo "[INFO] Cancers loaded from config: ${CANCERS[*]}"
echo "[INFO] Mutation groups: ${MUT_GROUPS[*]}"
echo "[INFO] FORCE_RERUN=$FORCE_RERUN"

# ============================================================
# 1) Prepare gene lists (significant overlap 03a ∩ 03b)
# ============================================================
if [[ "$RUN_PREP_GENELISTS" == "true" ]]; then
  run_r "01_transcriptomics/pipeline/cox/03e_prepare_gene_lists.R" \
        "${LOGDIR}/03e_prepare_gene_lists.log"
else
  echo "[SKIP] RUN_PREP_GENELISTS=false"
fi

# ============================================================
# 2) Gene-matched Expr + Mut Cox (03e_gene_expr_mut.R)
#     per cancer × mutation group
# ============================================================
if [[ "$RUN_GENE_EXPR_MUT" == "true" ]]; then
  for cancer in "${CANCERS[@]}"; do
    expr_in="${NORM_DIR}/TCGA_${cancer}_gene.normalized.csv"
    mani_in="${NORM_DIR}/TCGA_${cancer}_gene.sample_manifest.csv"

    if [[ ! -f "$expr_in" || ! -f "$mani_in" ]]; then
      echo "[SKIP] Missing expression/manifest for ${cancer} → skip all mut groups"
      continue
    fi

    for mg in "${MUT_GROUPS[@]}"; do
      mut_in="${MUT_DIR}/TCGA_${cancer}/gene/TCGA_${cancer}_gene_ensembl_${mg}.csv"
      if [[ ! -f "$mut_in" ]]; then
        echo "[SKIP] No mutation matrix for ${cancer} / ${mg}"
        continue
      fi

      out_dir="${OUT_03E}/TCGA_${cancer}"
      mkdir -p "$out_dir"

      out_csv="${out_dir}/TCGA_${cancer}_gene_${mg}.cox_results.csv"
      out_sum="${out_dir}/TCGA_${cancer}_gene_${mg}.summary.txt"
      gene_list="${GENE_LIST_DIR}/TCGA_${cancer}_gene_${mg}_gene_list.txt"

      if ! should_rerun "$out_csv"; then
        echo "[SKIP] 03e already exists (use FORCE_RERUN=true to overwrite): $out_csv"
        continue
      fi

      log_file="${LOGDIR}/03e_gene_expr_mut_${cancer}_${mg}.log"

      if [[ -f "$gene_list" ]]; then
        conda run -p "$R_ENV" --no-capture-output Rscript \
          "01_transcriptomics/pipeline/cox/03e_gene_expr_mut.R" \
          "$expr_in" "$mut_in" "$mani_in" "$COV_YAML" "$out_csv" "$out_sum" "$gene_list" \
          2>&1 | tee "$log_file"
      else
        conda run -p "$R_ENV" --no-capture-output Rscript \
          "01_transcriptomics/pipeline/cox/03e_gene_expr_mut.R" \
          "$expr_in" "$mut_in" "$mani_in" "$COV_YAML" "$out_csv" "$out_sum" \
          2>&1 | tee "$log_file"
      fi
    done
  done
else
  echo "[SKIP] RUN_GENE_EXPR_MUT=false"
fi

# ============================================================
# 3) Collect global summaries
# ============================================================
if [[ "$RUN_COLLECT_GLOBAL" == "true" ]]; then
  run_r "01_transcriptomics/pipeline/cox/03e_collect_global_summaries.R" \
        "${LOGDIR}/03e_collect_global_summaries.log"
else
  echo "[SKIP] RUN_COLLECT_GLOBAL=false"
fi

# ============================================================
# 4) Build gene-modality long table (05a)
# ============================================================
if [[ "$RUN_BUILD_MODALITY" == "true" ]]; then
  run_r "01_transcriptomics/pipeline/cox/05a_build_gene_modality_table.R" \
        "${LOGDIR}/05a_build_gene_modality_table.log"
else
  echo "[SKIP] RUN_BUILD_MODALITY=false"
fi

# ============================================================
# 5) Modality contribution summary (06)
# ============================================================
if [[ "$RUN_MODALITY_SUMMARY" == "true" ]]; then
  run_r "01_transcriptomics/pipeline/cox/06_modality_contribution_summary.R" \
        "${LOGDIR}/06_modality_contribution_summary.log"
else
  echo "[SKIP] RUN_MODALITY_SUMMARY=false"
fi

echo "------------------------------------------------------------"
echo "[DONE] 05_run_summaries complete."

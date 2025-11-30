#!/usr/bin/env bash
#SBATCH --job-name=03_plot_cox_boxes
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --output=01_transcriptomics/logs/03_plots/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/03_plots/%x-%A.err
set -euo pipefail

RESULTS_DIR=${RESULTS_DIR:-"01_transcriptomics/out/03_univariate_coxph"}
PLOT_OUT_DIR=${PLOT_OUT_DIR:-"01_transcriptomics/out/03_plots_mut"}
INDEX_OUT_DIR=${INDEX_OUT_DIR:-"01_transcriptomics/out/03_plots_mut/indexes"}
PLOT_SCRIPT=${PLOT_SCRIPT:-"01_transcriptomics/pipeline/scripts/03_plot_coxph.R"}
CANCERS_YAML=${CANCERS_YAML:-"01_transcriptomics/config/cancers.yaml"}

DATATYPES_DEFAULT="gene iso_frac iso_log mutation"
read -r -a DATATYPES <<< "${DATATYPES:-$DATATYPES_DEFAULT}"

FIG_NAME=${FIG_NAME:-"cox_hr_boxplots.png"}
FIG_TITLE=${FIG_TITLE:-"Hazard Ratio distributions by cancer"}
FILTER_EXPR=${FILTER_EXPR:-""}      # e.g., "FDR < 0.05"
LOGY=${LOGY:-"true"}                 # true|false
WINSOR=${WINSOR:-"0"}                # e.g., 0.01
CLIP=${CLIP:-"0.01,0.99"}            # "" to disable
ORDER=${ORDER:-"alphabet"}
WIDTH_IN=${WIDTH_IN:-"16"}
HEIGHT_IN=${HEIGHT_IN:-"9"}
DPI=${DPI:-"300"}
JITTER_N=${JITTER_N:-"0"}
JITTER_WIDTH=${JITTER_WIDTH:-"0.25"}

CONDA_ENV_PATH=${CONDA_ENV_PATH:-"/home/people/s184275/r-env"}

mkdir -p "${PLOT_OUT_DIR}" "${INDEX_OUT_DIR}" "01_transcriptomics/logs/03_plots"
[[ -f "${PLOT_SCRIPT}" ]] || { echo "[ERROR] Missing ${PLOT_SCRIPT}" >&2; exit 1; }
[[ -d "${RESULTS_DIR}" ]] || { echo "[ERROR] Missing ${RESULTS_DIR}" >&2; exit 1; }
[[ -f "${CANCERS_YAML}" ]] || { echo "[ERROR] Missing ${CANCERS_YAML}" >&2; exit 1; }

mapfile -t CANCERS < <(sed -n 's/^[[:space:]]*-[[:space:]]*//p' "${CANCERS_YAML}")
(( ${#CANCERS[@]} > 0 )) || { echo "[ERROR] No cancers parsed." >&2; exit 1; }

INDEX_TSV="${INDEX_OUT_DIR}/cox_index_$(date +%Y%m%d_%H%M%S).tsv"
printf "path\tcancer\tdata_type\n" > "${INDEX_TSV}"
present=0; missing=0
for cancer in "${CANCERS[@]}"; do
  for dtype in "${DATATYPES[@]}"; do
    fp="${RESULTS_DIR}/TCGA_${cancer}_${dtype}.cox_results.csv"
    if [[ -f "${fp}" ]]; then
      printf "%s\t%s\t%s\n" "${fp}" "${cancer}" "${dtype}" >> "${INDEX_TSV}"
      ((present++)) || true
    else
      echo "[WARN] Missing ${fp}" >&2
      ((missing++)) || true
    fi
  done
done
(( present > 0 )) || { echo "[ERROR] No result files found (.cox_results.csv)." >&2; exit 1; }
echo "[INFO] Index: ${INDEX_TSV} (${present} files; missing ${missing})"

args=(
  --index "${INDEX_TSV}"
  --out   "${PLOT_OUT_DIR}/${FIG_NAME}"
  --title "${FIG_TITLE}"
  --order "${ORDER}"
  --width  "${WIDTH_IN}"
  --height "${HEIGHT_IN}"
  --dpi    "${DPI}"
  --jitter_n "${JITTER_N}"
  --jitter_width "${JITTER_WIDTH}"
)
[[ -n "${FILTER_EXPR}" ]] && args+=( --filter "${FILTER_EXPR}" )
[[ "${LOGY}" =~ ^([Tt][Rr][Uu][Ee]|true)$ ]] && args+=( --logy )
[[ -n "${WINSOR}" && "${WINSOR}" != "0" ]] && args+=( --winsor "${WINSOR}" )
[[ -n "${CLIP}" ]] && args+=( --clip "${CLIP}" )

echo "[RUN] ${PLOT_SCRIPT} ${args[*]}"
conda run -p "${CONDA_ENV_PATH}" --no-capture-output \
  Rscript "${PLOT_SCRIPT}" "${args[@]}"
echo "[DONE] Figure at: ${PLOT_OUT_DIR}/${FIG_NAME}"

#!/usr/bin/env bash
#SBATCH --job-name=03_plot_coxph
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=01_transcriptomics/logs/03_plots/%x-%A.log
#SBATCH --error=01_transcriptomics/logs/03_plots/%x-%A.err
set -euo pipefail

RESULTS_BASE=${RESULTS_BASE:-"01_transcriptomics/out/03_univariate_coxph"}
PLOT_OUT_DIR=${PLOT_OUT_DIR:-"01_transcriptomics/out/03_plots_mut"}
INDEX_DIR=${INDEX_DIR:-"${PLOT_OUT_DIR}/indexes"}
SCRIPT_R=${SCRIPT_R:-"01_transcriptomics/pipeline/scripts/03_plot_coxph.R"}
CANCERS_YAML=${CANCERS_YAML:-"01_transcriptomics/config/cancers.yaml"}

mkdir -p "$PLOT_OUT_DIR" "$INDEX_DIR"

# --- Parse cancer list ---
mapfile -t CANCERS < <(sed -n 's/^[[:space:]]*-[[:space:]]*//p' "$CANCERS_YAML")
(( ${#CANCERS[@]} > 0 )) || { echo "[ERROR] No cancers found in YAML"; exit 1; }

# --- Define modes ---
MODES=(baseline combined interaction)

# --- Run per-mode plots ---
for MODE in "${MODES[@]}"; do
  RES_DIR="${RESULTS_BASE}/${MODE}"
  OUT_IMG="${PLOT_OUT_DIR}/cox_boxplot_${MODE}.png"
  INDEX_TSV="${INDEX_DIR}/cox_index_${MODE}.tsv"
  mkdir -p "$(dirname "$OUT_IMG")"

  echo "[MODE] ${MODE}  -> ${RES_DIR}"
  printf "path\tcancer\tdata_type\n" > "$INDEX_TSV"

  for f in ${RES_DIR}/TCGA_*_*.cox_results.csv; do
    [[ -f "$f" ]] || continue
    fname=$(basename "$f")
    cancer=$(echo "$fname" | cut -d'_' -f2)
    dtype=$(echo "$fname" | sed -E "s/^TCGA_${cancer}_//; s/\.cox_results\.csv//")
    printf "%s\t%s\t%s\n" "$f" "$cancer" "$dtype" >> "$INDEX_TSV"
  done

  echo "[RUN] Building boxplot for ${MODE}"
  conda run -p /home/people/s184275/r-env --no-capture-output \
    Rscript "$SCRIPT_R" \
      --mode "$MODE" \
      --index "$INDEX_TSV" \
      --out_dir "$PLOT_OUT_DIR"
done

# --- Build comparison plot & summary ---
echo "[RUN] Building comparison plot across modes"
conda run -p /home/people/s184275/r-env --no-capture-output \
  Rscript "$SCRIPT_R" \
    --mode "compare" \
    --in_dir "$PLOT_OUT_DIR"

echo "[DONE] All plots written to: $PLOT_OUT_DIR"

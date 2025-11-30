#!/usr/bin/env Rscript
# ============================================================
# 04b_plot_mutation_baseline.R
# ------------------------------------------------------------
# Summarizes univariate mutation-level CoxPH survival models:
#   - HR violin + boxplots per cancer
#   - Counts of significant mutated genes
#   - P-value histograms across cohorts
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(scales)
})

say <- function(...) message(sprintf(...))
outdir <- "01_transcriptomics/out/04_visuals/mutation_baseline"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

FDR_thresholds <- c(0.05, 0.1, 0.2)
PVAL_thresholds <- c(0.001, 0.01, 0.05)

# ============================================================
# Load all mutation-based CoxPH results
# ============================================================
say("[LOAD] Collecting mutation CoxPH results...")
res_files <- list.files("01_transcriptomics/out/03b_mutation_univariate_coxph",
                        pattern = "cox_results\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(res_files)) stop("No mutation CoxPH results found.")

res_all <- rbindlist(lapply(res_files, function(f) {
  dt <- fread(f, showProgress = FALSE)
  setnames(dt, tolower(names(dt)))
  dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", basename(f))]
  dt[, data_type := "mutation"]
  dt
}), fill = TRUE)
say("[INFO] Combined %d models from %d cancers", nrow(res_all), uniqueN(res_all$cancer))

# ============================================================
# HR plots
# ============================================================
res_sig <- res_all[FDR < 0.05 & is.finite(hr) & hr > 0]
if (!nrow(res_sig)) {
  say("[WARN] No significant mutation-level features found.")
} else {
  p_violin <- ggplot(res_sig, aes(x = cancer, y = hr)) +
    geom_violin(fill = "#D55E00", alpha = 0.5) +
    geom_boxplot(width = 0.15, color = "black", outlier.shape = NA) +
    geom_hline(yintercept = 1, color = "#c90028", linetype = "dashed") +
    scale_y_log10(labels = function(x) paste0(x, "×")) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mutation Hazard Ratios per Cancer",
         y = "Hazard Ratio (log10 scale)", x = "Cancer Type")
  ggsave(file.path(outdir, "Mutation_HR_violin_box_FDR05.pdf"), p_violin, width = 10, height = 7)
}

# ============================================================
# Significant mutation counts
# ============================================================
sig_summary <- rbindlist(lapply(FDR_thresholds, function(th) {
  res_all[, .(n_signif = sum(FDR < th, na.rm = TRUE),
              threshold = paste0("FDR<", th)), by = cancer]
}))
p_bar <- ggplot(sig_summary, aes(x = reorder(cancer, n_signif), y = n_signif + 1, fill = threshold)) +
  geom_col(color = "black") +
  scale_y_log10() +
  scale_fill_brewer(palette = "Reds") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Significant Mutated Genes per Cancer",
       y = "Count (log10 scale)", x = "Cancer Type")
ggsave(file.path(outdir, "Mutation_significant_genes_barplot.pdf"), p_bar, width = 10, height = 6)

# ============================================================
# P-value histograms
# ============================================================
pvals <- res_all[!is.na(p) & is.finite(p)]
p_hist <- ggplot(pvals, aes(x = p)) +
  geom_histogram(bins = 40, fill = "#D55E00", color = "black") +
  facet_wrap(~ cancer, scales = "free_y") +
  theme_bw(base_size = 9) +
  labs(title = "P-value Distributions (Mutation Baseline Models)",
       x = "p-value", y = "Count")
ggsave(file.path(outdir, "Mutation_Pvalue_histograms.pdf"), p_hist, width = 14, height = 8)

say("[DONE] Mutation baseline plots saved successfully.")

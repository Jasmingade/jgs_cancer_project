#!/usr/bin/env Rscript
# ============================================================
# 04d_plot_iso_mutation_interactions.R
# ------------------------------------------------------------
# Visualizes isoform × mutation interaction models:
#   - HR distribution for interaction terms
#   - Count of significant interactions per cancer
#   - Volcano-like interaction plots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(ggrepel)
})

say <- function(...) message(sprintf(...))
outdir <- "01_transcriptomics/out/04_visuals/iso_mut_interactions"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

FDR_thresholds <- c(0.05, 0.1)

# ============================================================
# Load all interaction CoxPH results
# ============================================================
say("[LOAD] Collecting Isoform×Mutation interaction results...")
res_files <- list.files("01_transcriptomics/out/03d_iso_mut_univariate_coxph",
                        pattern = "cox_results\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(res_files)) stop("No interaction CoxPH result files found.")

res_all <- rbindlist(lapply(res_files, function(f) {
  dt <- fread(f)
  dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", basename(f))]
  dt
}), fill = TRUE)

say("[INFO] Loaded %d interaction models across %d cancers", nrow(res_all), uniqueN(res_all$cancer))

# ============================================================
# Volcano plot per cancer
# ============================================================
for (c in unique(res_all$cancer)) {
  dt <- res_all[cancer == c & is.finite(interact_HR)]
  if (!nrow(dt)) next
  p_volcano <- ggplot(dt, aes(x = logHR_interaction, y = -log10(p_interaction))) +
    geom_point(aes(color = FDR_interaction < 0.05), alpha = 0.7) +
    scale_color_manual(values = c("TRUE" = "#E69F00", "FALSE" = "grey60")) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    theme_bw(base_size = 12) +
    labs(title = sprintf("%s — Isoform×Mutation Interaction Effects", c),
         x = "log(HR interaction term)", y = "-log10(p)", color = "FDR<0.05")
  ggsave(file.path(outdir, sprintf("%s_iso_mut_interaction_volcano.pdf", c)),
         p_volcano, width = 8, height = 6)
}

# ============================================================
# Interaction count summary
# ============================================================
sig_summary <- res_all[, .(n_signif = sum(FDR_interaction < 0.05, na.rm = TRUE)), by = cancer]
p_bar <- ggplot(sig_summary, aes(x = reorder(cancer, n_signif), y = n_signif)) +
  geom_col(fill = "#E69F00", color = "black") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Significant Isoform×Mutation Interactions per Cancer",
       y = "Count of significant interactions", x = "Cancer Type")
ggsave(file.path(outdir, "Isoform_Mutation_interaction_barplot.pdf"), p_bar, width = 10, height = 6)
say("[DONE] Isoform×Mutation interaction plots saved.")

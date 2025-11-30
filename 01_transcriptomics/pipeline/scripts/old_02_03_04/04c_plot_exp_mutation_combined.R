#!/usr/bin/env Rscript
# ============================================================
# 04c_plot_exp_mutation_combined.R
# ------------------------------------------------------------
# Visualizes combined CoxPH models (Expression + Mutation):
#   - Compare HR distributions of expression vs mutation
#   - Joint significance heatmaps
#   - Multi-threshold barplots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(scales)
})

say <- function(...) message(sprintf(...))
outdir <- "01_transcriptomics/out/04_visuals/exp_mutation_combined"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

FDR_thresholds <- c(0.05, 0.1)

# ============================================================
# Load all combined CoxPH results
# ============================================================
say("[LOAD] Collecting Expression+Mutation CoxPH results...")
res_files <- list.files("01_transcriptomics/out/03c_exp_mutation_univariate_coxph",
                        pattern = "cox_results\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(res_files)) stop("No combined CoxPH result files found.")

res_all <- rbindlist(lapply(res_files, function(f) {
  dt <- fread(f)
  dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", basename(f))]
  dt[, model_type := "exp_mut"]
  dt
}), fill = TRUE)

say("[INFO] Loaded %d combined models across %d cancers", nrow(res_all), uniqueN(res_all$cancer))

# ============================================================
# Separate HRs for expression and mutation covariates
# ============================================================
hr_long <- melt(res_all[, .(cancer, expr_HR, mut_HR)], id.vars = "cancer",
                variable.name = "covariate", value.name = "HR")

p_violin <- ggplot(hr_long, aes(x = cancer, y = HR, fill = covariate)) +
  geom_violin(alpha = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  geom_hline(yintercept = 1, color = "#c90028", linetype = "dashed") +
  scale_y_log10(breaks = c(0.5, 1, 2, 4), labels = c("0.5×", "1×", "2×", "4×")) +
  scale_fill_manual(values = c(expr_HR = "#0072B2", mut_HR = "#D55E00"),
                    labels = c("Expression", "Mutation")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Combined CoxPH: Expression vs Mutation Effects",
       y = "Hazard Ratio (log10 scale)", x = "Cancer Type")

ggsave(file.path(outdir, "Exp_Mut_combined_HR_violin_box.pdf"), p_violin, width = 12, height = 8)

# ============================================================
# Joint significance heatmap
# ============================================================
res_all[, sig_expr := expr_FDR < 0.05]
res_all[, sig_mut := mut_FDR < 0.05]

heat_dt <- res_all[, .N, by = .(cancer, sig_expr, sig_mut)]
heat_dt[, status := paste0(ifelse(sig_expr, "Expr", ""), ifelse(sig_mut, "+Mut", ""))]
heat_dt[, status := fcase(
  status == "", "None",
  status == "Expr", "Expr_only",
  status == "Mut", "Mut_only",
  status == "Expr+Mut", "Both"
)]

p_heat <- ggplot(heat_dt, aes(x = cancer, y = status, fill = N)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#2166AC", trans = "log10") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Joint Significance of Expression and Mutation",
       fill = "Feature Count (log10)", x = "Cancer", y = "Significance Category")

ggsave(file.path(outdir, "Exp_Mut_joint_significance_heatmap.pdf"), p_heat, width = 12, height = 6)
say("[DONE] Combined model plots saved.")

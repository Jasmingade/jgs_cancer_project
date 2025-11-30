#!/usr/bin/env Rscript
# ============================================================
# 03d_postprocess_and_plot_all_coxph.R
# ------------------------------------------------------------
# Combines results from:
#   - 03_univariate_coxph.R
#   - 03a_mutation_univariate_coxph.R
#   - 03b_exp_mutation_univariate_coxph.R
#   - 03c_iso_mut_univariate_coxph.R
# Generates:
#   ① Merged comparison table (per cancer × model)
#   ② HR boxplots for each model type
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

say <- function(...) message(sprintf(...))
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a

# -----------------------------
# Input/output directories
# -----------------------------
base_out <- "01_transcriptomics/out"
out_dir  <- file.path(base_out, "03d_summary_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dirs <- list(
  expr       = file.path(base_out, "03_univariate_coxph"),
  mutation   = file.path(base_out, "03a_mutation_univariate_coxph"),
  combined   = file.path(base_out, "03b_exp_mutation_univariate_coxph"),
  interaction= file.path(base_out, "03c_iso_mut_univariate_coxph")
)

# -----------------------------
# Helper: load all Cox results
# -----------------------------
load_results <- function(path, model) {
  files <- list.files(path, pattern = "\\.cox_results\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) return(data.table())
  say("[%s] Found %d result files", model, length(files))
  rbindlist(lapply(files, function(f) {
    dt <- fread(f, showProgress = FALSE)
    # Extract cancer from filename (TCGA_<CANCER>_*.csv)
    cancer <- str_match(basename(f), "TCGA_([A-Z0-9]+)_")[,2] %||% "UNKNOWN"
    dt[, `:=`(model = model, cancer = cancer)]
    dt
  }), fill = TRUE)
}


# -----------------------------
# Load and merge all models
# -----------------------------
res_list <- list(
  expression  = load_results(dirs$expr, "Expression-only"),
  mutation    = load_results(dirs$mutation, "Mutation-only"),
  combined    = load_results(dirs$combined, "Expr+Mutation"),
  interaction = load_results(dirs$interaction, "Expr×Mutation")
)

merged <- rbindlist(res_list, fill = TRUE)
if (nrow(merged) == 0) {
  stop("No results found across any model folders.")
}


# Ensure numeric columns (guard against string "NA"/"Inf")
num_cols <- c("HR", "p", "FDR", "logHR")
for (cn in num_cols) {
  if (cn %in% names(merged)) {
    merged[[cn]] <- suppressWarnings(as.numeric(merged[[cn]]))
  }
}

# Drop rows with invalid or non-finite HR
merged <- merged[!is.na(HR) & is.finite(HR) & HR > 0]
if (nrow(merged) == 0) stop("No valid HR values found across all models.")

# ------------------------------------------------------------
# Robust numeric cleaning of HR column
# ------------------------------------------------------------
# Coerce HR to numeric and remove invalid/overflow values
merged$HR <- suppressWarnings(as.numeric(merged$HR))
merged <- merged[!is.na(HR) & is.finite(HR) & HR > 0 & HR < 1e5]  # filter out absurd outliers
if (nrow(merged) == 0) stop("No valid HR values found across all models.")

# Coerce other numeric columns
for (cn in c("p", "FDR", "logHR")) {
  if (cn %in% names(merged)) {
    merged[[cn]] <- suppressWarnings(as.numeric(merged[[cn]]))
  }
}

# QC print
say("[QC] HR summary after cleaning: min=%.3f, max=%.3f, n=%d",
    min(merged$HR, na.rm=TRUE), max(merged$HR, na.rm=TRUE), nrow(merged))


# Compute log10 safely
merged[, HR_log10 := log10(HR)]

# Cleanup and formatting
merged[, cancer := factor(cancer, levels = sort(unique(cancer)))]
merged[, model := factor(model, levels = c("Expression-only","Mutation-only","Expr+Mutation","Expr×Mutation"))]

say("[QC] HR summary after cleaning: min=%.3f, max=%.3f, n=%d",
    min(merged$HR, na.rm=TRUE), max(merged$HR, na.rm=TRUE), nrow(merged))

# -----------------------------
# Summary per model/cancer
# -----------------------------
summary_dt <- merged[, .(
  n_features = .N,
  median_HR = median(HR, na.rm = TRUE),
  median_logHR = median(logHR, na.rm = TRUE),
  sig_p = sum(p < 0.05, na.rm = TRUE),
  sig_FDR = sum(FDR < 0.05, na.rm = TRUE)
), by = .(cancer, model)]

fwrite(merged, file.path(out_dir, "merged_all_models_long.csv"))
fwrite(summary_dt, file.path(out_dir, "summary_counts_by_model.csv"))
say("[OK] Wrote merged comparison tables to: %s", out_dir)

# -----------------------------
# Boxplots of HR per model
# -----------------------------
plot_hr_box <- function(dt, model_name, out_file) {
  if (nrow(dt) == 0) return(NULL)
  p <- ggplot(dt, aes(x = cancer, y = HR, fill = cancer)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
    scale_y_log10() +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle=60, hjust=1, size=7),
          legend.position="none") +
    labs(title = paste0("Hazard Ratios: ", model_name),
         subtitle = "Log10 scale; each box = HR distribution per cancer",
         x = "Cancer type", y = "Hazard Ratio (HR, log10)") +
    geom_hline(yintercept = 1, linetype="dashed", color="red")
  ggsave(out_file, p, width = 9, height = 4)
  say("[PLOT] %s → %s", model_name, out_file)
}

# Generate one boxplot per model
for (m in unique(merged$model)) {
  sub <- merged[model == m]
  plot_hr_box(sub, m, file.path(out_dir, paste0("HR_boxplot_", gsub("[^A-Za-z0-9]+","_",m), ".pdf")))
}

# -----------------------------
# Combined multi-model boxplot
# -----------------------------
p_all <- ggplot(merged, aes(x = model, y = HR, fill = model)) +
  geom_boxplot(outlier.size=0.5, alpha=0.8) +
  scale_y_log10() +
  theme_bw(base_size = 11) +
  labs(title = "Overall HR distributions across model types",
       x = "Model type", y = "Hazard Ratio (log10)") +
  geom_hline(yintercept = 1, linetype="dashed", color="red") +
  theme(legend.position="none")
ggsave(file.path(out_dir, "HR_boxplot_all_models.pdf"), p_all, width = 7, height = 5)

say("[DONE] HR boxplots written to %s", out_dir)

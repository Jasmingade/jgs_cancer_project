#!/usr/bin/env Rscript
# ============================================================
# 04a_plot_baseline_expression.R
# ------------------------------------------------------------
# Generates visual summaries for univariate expression-only CoxPH models:
#   - HR boxplots (per data type and cancer)
#   - Significant feature barplots
#   - P-value histogram matrix
#   - Optional covariate sanity forest plot (if available)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(ggbeeswarm)
  library(tidyr)
  library(dplyr)
})

say <- function(...) message(sprintf(...))
outdir <- "01_transcriptomics/out/04_visuals/expression_baseline"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load all CoxPH results
# ============================================================
say("[LOAD] Collecting expression-only CoxPH results...")

# --- Find all result files ---
data_types <- c("gene", "iso_log", "iso_frac")
res_files <- unlist(lapply(data_types, function(dt)
  list.files(file.path("01_transcriptomics/out/03a_univariate_coxph", dt),
             pattern = "cox_results\\.csv$", full.names = TRUE, recursive = TRUE)))

# --- Check presence ---
if (length(res_files) == 0)
  stop("No CoxPH result files found. Did you run 03a_univariate_coxph.R?")

say("[INFO] Found %d result files across %d data types", length(res_files), length(data_types))

# --- Quick integrity check ---
say("[DEBUG] Checking HR column integrity per file...")
for (f in res_files) {
  dt <- fread(f, nrows = 10)
  if (!"HR" %in% names(dt) && !"hr" %in% names(dt)) {
    say("[WARN] File %s missing HR column.", f)
    next
  }

  colname <- if ("HR" %in% names(dt)) "HR" else "hr"
  vals <- suppressWarnings(as.numeric(dt[[colname]]))
  if (any(is.na(vals))) {
    say("[WARN] File %s contains non-numeric HR entries.", f)
  } else if (any(vals > 1000 | vals < 0.001)) {
    say("[WARN] File %s has extreme HR range: %.3e–%.3e", f, min(vals), max(vals))
  }
}

# ============================================================
# Combine, clean, and normalize HR values
# ============================================================
res_all <- rbindlist(lapply(res_files, function(f) {
  # --- Read file safely ---
  dt <- tryCatch(fread(f, strip.white = TRUE), error = function(e) {
    say("[ERROR] Failed to read file: %s", f)
    return(NULL)
  })
  if (is.null(dt) || nrow(dt) == 0) {
    say("[WARN] File %s is empty. Skipping.", f)
    return(NULL)
  }

  # --- Normalize column names ---
  setnames(dt, tolower(names(dt)))
  if ("hr" %in% names(dt))    setnames(dt, "hr", "HR")
  if ("fdr" %in% names(dt))   setnames(dt, "fdr", "FDR")
  if ("loghr" %in% names(dt)) setnames(dt, "loghr", "logHR")
  if ("p" %in% names(dt)) {
    dt[, p := suppressWarnings(as.numeric(as.character(p)))]
    }
    say("[DEBUG] %s: mean(p) = %.3g (non-NA %d/%d)",
    basename(f),
    mean(dt$p, na.rm = TRUE),
    sum(!is.na(dt$p)), nrow(dt))

  # --- Ensure HR column exists ---
  if (!"HR" %in% names(dt)) {
    if ("logHR" %in% names(dt)) {
      say("[INFO] %s: Reconstructing HR from logHR", basename(f))
      dt[, HR := exp(as.numeric(logHR))]
    } else if ("beta" %in% names(dt)) {
      say("[INFO] %s: Reconstructing HR from beta", basename(f))
      dt[, HR := exp(as.numeric(beta))]
    } else {
      say("[WARN] %s: Missing HR/logHR/beta column, skipping file.", basename(f))
      return(NULL)
    }
  }

  # --- Coerce HR safely ---
  dt[, HR := suppressWarnings(as.numeric(HR))]

  # --- Compute descriptive stats ---
  med_hr <- median(dt$HR, na.rm = TRUE)
  rng_hr <- range(dt$HR, na.rm = TRUE)

  # --- Handle only extreme outliers per feature ---
  extreme_idx <- which(dt$HR > 1e6 | dt$HR < 1e-6 | !is.finite(dt$HR))
  if (length(extreme_idx) > 0) {
    say("[WARN] %s: %d extreme HR values (%.2e–%.2e). Replacing with NA instead of log-transforming.",
        basename(f), length(extreme_idx), rng_hr[1], rng_hr[2])
    dt[extreme_idx, HR := NA_real_]
  }

  # --- Recalculate summary after cleanup ---
  rng_hr <- range(dt$HR, na.rm = TRUE)
  med_hr <- median(dt$HR, na.rm = TRUE)

  # --- Check overall HR range sanity ---
  if (all(dt$HR > 0.05 & dt$HR < 50, na.rm = TRUE)) {
    say("[INFO] %s: HR range (%.2f–%.2f) looks fine. Keeping as-is.",
        basename(f), rng_hr[1], rng_hr[2])
  } else if (med_hr < 0) {
    say("[INFO] %s: Detected log(HR) scale. Exponentiating to get HR.", basename(f))
    dt[, HR := exp(HR)]
  } else {
    say("[INFO] %s: HR distribution (%.2f–%.2f) within plausible limits. Keeping as-is.",
        basename(f), rng_hr[1], rng_hr[2])
  }

  say("[DEBUG] %s: HR median %.2f, range %.2e–%.2e after cleanup.",
      basename(f), median(dt$HR, na.rm = TRUE),
      min(dt$HR, na.rm = TRUE), max(dt$HR, na.rm = TRUE))

  # --- Filter only after HR is stable ---
  dt <- dt[is.finite(HR) & HR > 0]

  # --- Annotate metadata ---
  dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", basename(f))]
  dt[, data_type := sub("^TCGA_[A-Z0-9]+_(.*?)\\.cox_results\\.csv$", "\\1", basename(f))]

  return(dt)
}), fill = TRUE)

say("[INFO] Combined table: %d rows", nrow(res_all))





# Define color palette for data types
type_colors <- c(
    gene     = "#0072B2",  # blue
    iso_log  = "#009E73",  # green
    iso_frac = "#E69F00"   # orange
  )
  
# ============================================================
# HR boxplots per cancer (grouped by data type)
# ============================================================
say("[PLOT] Generating grouped HR boxplots (significant features only)")

# --- Filter significant features only ---
res_sig <- res_all[FDR < 0.05 & is.finite(HR) & HR > 0]
if (nrow(res_sig) == 0) {
  say("[WARN] No significant features found (FDR < 0.05). Skipping boxplot.")
} else {
  say("[INFO] Plotting %d significant features.", nrow(res_sig))
  
  # 🔍 Check whether HR looks log-transformed
  say("[DEBUG] HR summary before plotting:")
  print(summary(res_sig$HR))
  say("[DEBUG] HR range: %.3f–%.3f", min(res_sig$HR, na.rm = TRUE), max(res_sig$HR, na.rm = TRUE))
  
  if (median(res_sig$HR, na.rm = TRUE) < 0.5 || median(res_sig$HR, na.rm = TRUE) > 5) {
    say("[WARN] HR values seem off — check if already log-transformed.")
  }

  p_box <- ggplot(res_sig, aes(x = cancer, y = HR, color = data_type, group = interaction(cancer, data_type))) +
    geom_boxplot(
      fill = "white",               # white interior
      outlier.size = 0.6,
      position = position_dodge(width = 0.9, preserve = "single"),
      width = 0.8,                  # wider boxes
      linewidth = 0.7
    ) +
    scale_fill_identity() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#c90028ff", linewidth = 0.6) +
    scale_y_log10(
        breaks = c(0.25, 0.5, 1, 2, 4, 8),
        labels = function(x) paste0(x, "×"),
        limits = c(0.3, 8)
    ) +
    scale_color_manual(values = type_colors, name = "Data Type") +
    theme_bw(base_size = 12) +
    theme(
        panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 9)
        ) +
    labs(
      title = "Significant Hazard Ratios per Cancer (FDR < 0.05)",
      y = "Hazard Ratio (log10 scale)",
      x = "Cancer Type"
    )

  out_file <- file.path(outdir, "HR_boxplots_expression_significant_grouped.pdf")
  ggsave(out_file, p_box, width = 15, height = 8)
  say("[DONE] Grouped HR boxplot (significant only) saved: %s", out_file)
}



# ============================================================
# Violin + Boxplot plots per cancer (grouped by data type)
# ============================================================
say("[PLOT] Generating violin + box overlay HR plot (significant features only)")

res_sig <- res_all[FDR < 0.05 & is.finite(HR) & HR > 0]
if (nrow(res_sig) == 0) {
  say("[WARN] No significant features found (FDR < 0.05). Skipping violin plot.")
} else {
  say("[INFO] Plotting %d significant features.", nrow(res_sig))
  
  # 🔍 Check whether HR looks log-transformed
  say("[DEBUG] HR summary before plotting:")
  print(summary(res_sig$HR))
  say("[DEBUG] HR range: %.3f–%.3f", min(res_sig$HR, na.rm = TRUE), max(res_sig$HR, na.rm = TRUE))
  
  if (median(res_sig$HR, na.rm = TRUE) < 0.5 || median(res_sig$HR, na.rm = TRUE) > 5) {
    say("[WARN] HR values seem off — check if already log-transformed.")
  }

  res_sig <- res_sig %>%
    complete(cancer, data_type)

  dodge_width <- 0.7

  p_violin_box <- ggplot(res_sig, aes(x = cancer, y = HR, fill = data_type, group = interaction(cancer, data_type))) +
    # --- Violin layer ---
    geom_violin(
      position = position_dodge(width = dodge_width, preserve = "single"),
      width = 1.8,
      alpha = 0.4,
      #scale = "width",
      color = NA,             # remove black border (cleaner)
      trim = TRUE,
      na.rm = TRUE
    ) +
    
    # --- Boxplot layer ---
    geom_boxplot(
      aes(group = interaction(cancer, data_type)),  # ensure grouping for dodge
      width = 0.15, 
      position = position_dodge(width = dodge_width, preserve = "single"),
      outlier.shape = NA,
      linewidth = 0.3,
      alpha = 0.6,
      color = "black",
      na.rm = TRUE
    ) +

    # --- Jitter overlay ---
    geom_quasirandom(aes(group = interaction(cancer, data_type)), 
                    dodge.width = dodge_width,
                    shape = 21,
                    color = "black",
                    size = 0.2, 
                    alpha = 0.4,
                    stroke = 0.2) +

    # --- Horizontal reference line ---
    geom_hline(yintercept = 1, linetype = "dashed", color = "#c90028", linewidth = 0.6) +

    scale_y_log10(
      breaks = c(0.25, 0.5, 1, 2, 4, 8),
      labels = function(x) paste0(x, "×"),
      limits = c(0.3, 8)
    ) +

    scale_fill_manual(values = type_colors, name = "Data Type") +

    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
      legend.position = "right",
      legend.key.size = unit(0.5, "cm")
    ) +
    labs(
      title = "Significant Hazard Ratios per Cancer (FDR < 0.05)",
      subtitle = "Violin + box + points overlay; dashed line = HR = 1",
      x = "Cancer Type",
      y = "Hazard Ratio (log10 scale)"
    )

  out_file <- file.path(outdir, "HR_violin_box_expression_refined.pdf")
  ggsave(out_file, p_violin_box, width = 15, height = 8)
  say("[DONE] Refined grouped HR violin-box plot saved: %s", out_file)
}

print(table(res_sig$cancer, res_sig$data_type))

say("[DEBUG] Checking which groups were plotted in the violin layer...")

say("[DEBUG] Inspecting plotted groups...")
violins <- ggplot_build(p_violin_box)$data[[1]] %>% dplyr::count(PANEL, group)
boxes   <- ggplot_build(p_violin_box)$data[[2]] %>% dplyr::count(PANEL, group)
say("[DEBUG] Violin groups drawn: %d", nrow(violins))
say("[DEBUG] Boxplot groups drawn: %d", nrow(boxes))

# ============================================================
# Significant feature barplot (FDR < 0.05)
# ============================================================
say("[PLOT] Generating significant feature barplot")

sig_summary <- res_all[, .(
  n_signif = sum(FDR < 0.05, na.rm = TRUE)
), by = .(cancer, data_type)]

p_bar <- ggplot(sig_summary, aes(x = reorder(cancer, n_signif), y = n_signif + 1,
                                 color = data_type, fill = "white")) +
  geom_col(
    position = position_dodge(width = 0.9),
    width = 0.7,
    linewidth = 0.8,
    show.legend = TRUE
  ) +

  scale_y_log10(labels = scales::comma) +
  scale_color_manual(values = type_colors, name = "Data Type") +
  scale_fill_manual(values = c("white" = "white"), guide = "none") +

  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Significant Features per Cancer (FDR < 0.05)",
    subtitle = "White bars with colored outlines; y-axis in log10 scale",
    x = "Cancer Type",
    y = "Count of Significant Features (log10)",
    color = "Data Type",
    fill = "white"
  )


ggsave(file.path(outdir, "Significant_features_barplot.pdf"), p_bar, width = 15, height = 6)
say("[DONE] Significant feature barplot saved.")

# ============================================================
# P-value histogram matrix
# ============================================================
say("[PLOT] Generating p-value histogram matrix")
pvals <- res_all[!is.na(p) & is.finite(p) & p >= 0 & p <= 1]
say("[DEBUG] P-value histogram: %d rows with valid p", nrow(pvals))
print(head(pvals))

p_hist <- ggplot(res_all, aes(x = p)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  facet_grid(data_type ~ cancer, scales = "free_y") +
  theme_bw(base_size = 9) +
  theme(
    strip.text.x = element_text(angle = 0, hjust = 0),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    title = "P-value Distributions by Cancer and Data Type",
    x = "P-value",
    y = "Count"
  )

ggsave(file.path(outdir, "Pvalue_histogram_matrix.pdf"), p_hist, width = 16, height = 8)
say("[DONE] P-value histogram matrix saved.")


for (dt in unique(res_all$data_type)) {
  say("[PLOT] P-value histogram for %s", dt)
  sub <- res_all[data_type == dt & !is.na(p) & p >= 0 & p <= 1]
  
  p_hist <- ggplot(sub, aes(x = p)) +
    geom_histogram(bins = 40, fill = "steelblue", color = "white") +
    facet_wrap(~ cancer, scales = "free_y", ncol = 6) +
    theme_bw(base_size = 9) +
    theme(
      strip.text.x = element_text(angle = 45, hjust = 0),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6)
    ) +
    labs(
      title = sprintf("P-value Distributions — %s", dt),
      x = "P-value",
      y = "Count"
    )

  out_file <- file.path(outdir, sprintf("Pvalue_histogram_%s.pdf", dt))
  ggsave(out_file, p_hist, width = 10, height = 8)
  say("[DONE] P-value histogram saved: %s", out_file)
}

# ============================================================
# Quantitative summary of significant features
# ============================================================
say("[SUMMARY] Computing per-cancer significance summary...")

if (exists("res_all")) {

  summary_dt <- res_all %>%
    filter(is.finite(HR), HR > 0) %>%
    group_by(cancer, data_type) %>%
    summarise(
      total_features = n(),
      sig_features   = sum(FDR < 0.05, na.rm = TRUE),
      frac_signif    = sig_features / total_features,
      median_HR      = median(HR[FDR < 0.05], na.rm = TRUE),
      IQR_HR         = IQR(HR[FDR < 0.05], na.rm = TRUE)
    ) %>%
    arrange(desc(frac_signif))

  # --- Save quantitative summary ---
  out_summary <- file.path(outdir, "significance_summary_by_cancer.csv")
  fwrite(summary_dt, out_summary)
  say("[DONE] Summary table saved: %s", out_summary)

  # --- Optional barplot for quick overview ---
  p_summary <- ggplot(summary_dt, aes(x = reorder(cancer, -frac_signif),
                                      y = frac_signif, fill = data_type)) +
    geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = type_colors, name = "Data Type") +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 9)
    ) +
    labs(
      title = "Fraction of Significant Features per Cancer Type",
      x = "Cancer Type",
      y = "Fraction (FDR < 0.05)"
    )

  out_plot <- file.path(outdir, "fraction_significant_features_barplot.pdf")
  ggsave(out_plot, p_summary, width = 12, height = 6)
  say("[DONE] Fraction significant barplot saved: %s", out_plot)
}


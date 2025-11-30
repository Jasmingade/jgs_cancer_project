#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# Short log helper
say <- function(...) message(sprintf(...))

# ============================================================
# Load precomputed results
# ============================================================
manifest_counts_path <- "01_transcriptomics/out/02_norm/manifest_covariate_counts.csv"
covariate_dir <- "01_transcriptomics/out/03_univariate_coxph"

say("[LOAD] Reading sanity check and manifest data...")
mani_cov <- fread(manifest_counts_path)

all_cov_files <- list.files(
  path = covariate_dir,
  pattern = "_covariate_sanity\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
cov_all <- rbindlist(lapply(all_cov_files, fread), fill = TRUE)


# ============================================================
# UNIFIED FOREST PLOT: Age & Stage sanity check (gene-level only)
# ============================================================
say("[PLOT] Generating unified age/stage forest plot (gene-level only)")

manifest_counts_path <- "01_transcriptomics/out/02_norm/manifest_covariate_counts.csv"
if (!file.exists(manifest_counts_path)) {
  die("Manifest covariate counts file not found: %s", manifest_counts_path)
}
mani_cov <- fread(manifest_counts_path)
mani_cov[, `:=`(
  pct_age = round(100 * n_age / n_total, 1),
  pct_stage = round(100 * n_stage / n_total, 1)
)]

all_cov_files <- list.files(
  path = "01_transcriptomics/out/03_univariate_coxph",
  pattern = "_covariate_sanity\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(all_cov_files)) {
  say("[WARN] No covariate sanity CSVs found.")
  quit(status = 0)
}

say("[INFO] Found %d covariate sanity results", length(all_cov_files))

cov_all <- rbindlist(lapply(all_cov_files, function(f) {
  dt <- fread(f)
  cancer <- sub("^.*TCGA_([A-Z0-9]+)_.*$", "\\1", basename(f))
  data_type <- if (grepl("/gene/", f)) "gene"
  else if (grepl("/iso_log/", f)) "iso_log"
  else if (grepl("/iso_frac/", f)) "iso_frac"
  else "unknown"

  mani_path <- file.path("01_transcriptomics/out/02_norm",
                         sprintf("TCGA_%s_%s.sample_manifest.csv", cancer, data_type))
  n_samples <- if (file.exists(mani_path)) nrow(fread(mani_path)) else NA_integer_
  dt[, `:=`(cancer = cancer, data_type = data_type, n_samples = n_samples)]
  dt
}), fill = TRUE)

cov_all <- cov_all[data_type == "gene"]
if (nrow(cov_all) == 0) {
  say("[WARN] No results for selected data type ('gene'). Skipping plot.")
  quit(status = 0)
}

cov_all <- merge(cov_all, mani_cov, by = "cancer", all.x = TRUE)

# Filter and flag
MIN_STAGE_PCT_FOR_MODEL <- 20
MIN_STAGE_N_FOR_MODEL   <- 20
LOW_STAGE_PCT_FLAG      <- 25
say("[FILTER] Keeping 'stage' fits with >= %d%% coverage AND >= %d staged samples",
    MIN_STAGE_PCT_FOR_MODEL, MIN_STAGE_N_FOR_MODEL)

valid_stage <- cov_all[
  covariate == "stage" &
  !is.na(pct_stage) & !is.na(n_stage) &
  pct_stage >= MIN_STAGE_PCT_FOR_MODEL &
  n_stage   >= MIN_STAGE_N_FOR_MODEL
]
invalid_stage <- cov_all[covariate == "stage" &
                         (is.na(pct_stage) | pct_stage < MIN_STAGE_PCT_FOR_MODEL |
                          is.na(n_stage)   | n_stage   < MIN_STAGE_N_FOR_MODEL)]
if (nrow(invalid_stage) > 0) {
  say("[INFO] Skipped %d stage fits (<%d%% coverage or <%d staged samples): %s",
      nrow(invalid_stage), MIN_STAGE_PCT_FOR_MODEL, MIN_STAGE_N_FOR_MODEL,
      paste(unique(invalid_stage$cancer), collapse = ", "))
}

cov_all <- rbindlist(list(cov_all[covariate == "age"], valid_stage), fill = TRUE)
cov_all[, flag_low_n := (covariate == "stage" & pct_stage < LOW_STAGE_PCT_FLAG)]

# Classify significance (p-value)
cov_all <- cov_all[is.finite(HR) & HR > 0]
cov_all[, sig_group := ifelse(p < 0.05, "Significant", "Not significant")]
cov_all[, cancer_label := sprintf("%s (n=%s)", cancer, ifelse(is.na(n_total), "?", n_total))]

cov_all <- unique(cov_all, by = c("cancer", "covariate"))
cov_all <- cov_all[order(n_total, decreasing = TRUE)]
cov_all[, cancer_label := factor(cancer_label, levels = unique(cancer_label))]

cov_all[, HR := pmin(pmax(HR, 0.25), 4)]
cov_all[, HR_lo := pmax(HR_lo, 0.25)]
cov_all[, HR_hi := pmin(HR_hi, 4)]

sig_colors <- c("Significant" = "#1B9E77", "Not significant" = "#7f7f7f")
shape_map <- c("FALSE" = 16, "TRUE" = 1)



facet_dt <- data.table(
  covariate = c("age", "stage"),
  HR_min = c(0.8, 0.8),
  HR_max = c(1.4, 5)
)

cov_all <- merge(cov_all, facet_dt, by = "covariate", all.x = TRUE)
cov_all <- cov_all[HR >= HR_min & HR <= HR_max]

mean_hr <- cov_all[sig_group == "Significant",
                   .(mean_HR = exp(mean(log(HR), na.rm = TRUE))),
                   by = covariate]

mean_hr <- merge(mean_hr, facet_dt, by = "covariate", all.x = TRUE)

# Debug
print("DEBUG:")
print(cov_all[, .(min_HR = min(HR), max_HR = max(HR)), by = covariate])

# Plot
p <- ggplot(cov_all, aes(y = cancer_label, x = HR, color = sig_group, shape = flag_low_n)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = HR_lo, xmax = HR_hi), width = 0.25) +

  # Mean HR
  geom_vline(
    data = mean_hr,
    aes(xintercept = mean_HR),
    color = "#c90028ff", linetype = "dashed", linewidth = 0.6
  ) +
    geom_vline(
      xintercept = 1, linetype = "solid", color = "black"
    ) +

  # Add text annotation for mean HR line
  geom_label_repel(
    data = mean_hr,
    aes(x = mean_HR, y = 10, label = sprintf("Mean HR = %.2f", mean_HR)),
    vjust = 0, hjust = 0, color = "#c90028ff", size = 2.5, fill = "#ffffffff",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~covariate, nrow = 1, scales = "free_x", strip.position = "top",
    labeller = as_labeller(c(
      age = "Age (per year increase)",
      stage = "Tumor stage (I-IV)"
    ))
  ) +

  # Log10 x-axis with clean breaks
  scale_x_log10(
    breaks = function(x) scales::log_breaks(n = 5)(x),
    labels = function(x) sprintf("%.1f×", x)
  ) +

  # Custom colors and shapes
  scale_color_manual(values = sig_colors, name = "Significance (p-value < 0.05)") +
  scale_shape_manual(values = shape_map, name = "Low Coverage (<25%)") +

  # --- Polished theme ---
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.y = element_text(size = 8, color = "gray15"),
    axis.text.x = element_text(size = 9, color = "gray20"),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(size = 13, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 10, hjust = 0, color = "gray25"),
    plot.margin = margin(5, 15, 5, 5)
  ) +
  labs(
    title = "Clinical Covariates (Age & Stage)",
    subtitle = sprintf(
      "Faceted by covariate; solid black = mean HR among significant cancers; dashed red = HR = 1\nOpen circles = low stage coverage (<%d%%)",
      LOW_STAGE_PCT_FLAG),
    x = "Hazard Ratio (log10 scale)",
    y = "Cancer Type"
  ) +
  geom_text(
    data = cov_all[covariate == "age"], 
    aes(label = paste0(pct_age, "%")), 
    hjust = -0.3, size = 2.8, color = "black") +
  geom_text(
    data = cov_all[covariate == "stage"], 
    aes(label = paste0(pct_stage, "%")), 
    hjust = -0.3, size = 2.8, color = "black"
  ) 

forest_out_dir <- "01_transcriptomics/out/03_univariate_coxph/forest_plots"
dir.create(forest_out_dir, recursive = TRUE, showWarnings = FALSE)
out_png <- file.path(forest_out_dir, "HR_forestplot_age_stage.png")
ggsave(out_png, p, width = 12, height = 6)
say("[PLOT] Unified age/stage forest plot saved: %s", out_png)
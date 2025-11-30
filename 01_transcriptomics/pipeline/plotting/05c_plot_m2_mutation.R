#!/usr/bin/env Rscript
# 05c_plot_m2_mutation.R
# Model 2 – Mutation-only CoxPH (03b)
# -----------------------------------
# Uses:
#   - Significant outputs:  TCGA_<CANCER>_mutation_<GROUP>.cox_results.csv
#     (for all HR-based plots and counts)
#   - Optional FULL outputs: TCGA_<CANCER>_mutation_<GROUP>.cox_results_full.csv
#     (ONLY for p-value histograms, when M2_USE_FULL_RESULTS=true)
#
# Env toggle:
#   M2_USE_FULL_RESULTS = "true"/"1"/"yes" → load *_full.csv additionally for p-hist
#   otherwise                              → skip p-hist (everything else still works)
#
# Produces:
#   - m2_hr_boxplot_mutation_log2_zoom.png
#   - m2_hr_jitter_mutation_log2_zoom.png
#   - m2_significant_mutation_barplot.png
#   - mutation_sig_all_results.csv
#   - (optional) mutation_full_all_results_for_phist.csv
#   - (optional) m2_pvalue_histograms.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(dplyr)
})

say     <- function(...) message(sprintf(...))
to_bool <- function(x) tolower(x) %in% c("1","true","yes")

# ---------------------------------------------------------------
# PATHS & TOGGLES
# ---------------------------------------------------------------
root <- "01_transcriptomics/out/03b_mutation_univariate_coxph"
plot_dir <- "01_transcriptomics/out/05_plots/model2_mut"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

USE_FULL <- to_bool(Sys.getenv("M2_USE_FULL_RESULTS", "false"))
say("[INFO] M2_USE_FULL_RESULTS (for p-hist) = %s", USE_FULL)
say("[INFO] Loading mutation Cox results (SIGNIFICANT files) from: %s", root)
MIN_SIG_PER_GROUP <- as.integer(Sys.getenv("M2_MIN_SIG_PER_GROUP", "2"))
if (is.na(MIN_SIG_PER_GROUP) || MIN_SIG_PER_GROUP < 1) {
  MIN_SIG_PER_GROUP <- 1
}
ylim_lo <- as.numeric(Sys.getenv("M2_YLIM_LO", "-1.5"))
ylim_hi <- as.numeric(Sys.getenv("M2_YLIM_HI", "10"))
if (!is.finite(ylim_lo) || !is.finite(ylim_hi) || ylim_lo >= ylim_hi) {
  stop("Invalid y-axis limits defined via M2_YLIM_LO / M2_YLIM_HI")
}
say("[INFO] Minimum significant features per cancer/group = %d", MIN_SIG_PER_GROUP)
say("[INFO] log2(HR) plotting window = [%.2f, %.2f]", ylim_lo, ylim_hi)

# ---------------------------------------------------------------
# 0) Always load SIGNIFICANT ONLY results for main plots
#     (TCGA_<CANCER>_mutation_<GROUP>.cox_results.csv)
# ---------------------------------------------------------------
files_sig <- list.files(
  root,
  pattern = "cox_results\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
# This pattern does NOT match *_full.csv, so it's safe.

if (length(files_sig) == 0)
  stop("No significant cox_results.csv files found in 03b_mutation_univariate_coxph")

say("[INFO] Found %d sig-only mutation result files", length(files_sig))

# ---------------------------------------------------------------
# Loader for mutation CoxPH results
# ---------------------------------------------------------------
load_mut_file <- function(f) {
  dt <- fread(f, showProgress = FALSE)

  needed <- c("feature", "beta", "HR", "p")
  if (!all(needed %in% names(dt))) {
    warning(sprintf("[WARN] Skipping %s — missing required columns", f))
    return(NULL)
  }

  dt[, p  := as.numeric(p)]
  dt[, HR := as.numeric(HR)]

  # Add FDR if missing
  if (!"FDR" %in% names(dt)) {
    dt[, FDR := p.adjust(p, "BH")]
  } else {
    dt[, FDR := as.numeric(FDR)]
  }

  fname <- basename(f)

  # Only create cancer/mut_group if not already there
  if (!"cancer" %in% names(dt)) {
    dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", fname)]
  }
  if (!"mut_group" %in% names(dt)) {
    dt[, mut_group :=
          sub("^TCGA_[A-Z0-9]+_mutation_(.*)\\.cox_results.*$", "\\1", fname)]
  }

  # Keep valid rows
  dt <- dt[
    is.finite(HR) & HR > 0 &
    is.finite(p)  & p >= 0 & p <= 1
  ]

  dt
}

# ---------------------------------------------------------------
# 1) Load sig-only files for main plots
# ---------------------------------------------------------------
res_list_sig <- lapply(files_sig, load_mut_file)
res_list_sig <- Filter(Negate(is.null), res_list_sig)
res_all_sig  <- rbindlist(res_list_sig, fill = TRUE)

say("[INFO] Loaded %d total mutation rows (significant files combined)",
    nrow(res_all_sig))

# Export combined SIG results
fwrite(res_all_sig, file.path(plot_dir, "mutation_sig_all_results.csv"))

# Distinct mutation groups (from sig results)
mut_levels_sig <- sort(unique(res_all_sig$mut_group))
mut_colors_sig <- scales::hue_pal()(length(mut_levels_sig))
names(mut_colors_sig) <- mut_levels_sig

# ---------------------------------------------------------------
# 2) Filter significant rows (FDR<0.05) & log2 transform
# ---------------------------------------------------------------
res_sig <- res_all_sig[FDR < 0.05 & is.finite(HR) & HR > 0]

if (nrow(res_sig) == 0)
  stop("No significant mutation features found (FDR<0.05). Cannot plot.")

say("[INFO] Significant rows: %d", nrow(res_sig))

group_counts <- res_sig[, .(n_signif = .N), by = .(cancer, mut_group)]
if (MIN_SIG_PER_GROUP > 1) {
  keep_pairs <- group_counts[n_signif >= MIN_SIG_PER_GROUP,
                             .(cancer, mut_group)]
  res_sig <- res_sig[keep_pairs, on = .(cancer, mut_group)]
  if (nrow(res_sig) == 0) {
    stop(sprintf(
      "No cancer/mutation groups have at least %d significant features.",
      MIN_SIG_PER_GROUP
    ))
  }
  say("[INFO] After requiring >=%d sig features per cancer/group: %d rows kept",
      MIN_SIG_PER_GROUP, nrow(res_sig))
} else {
  say("[INFO] Minimum per cancer/group set to 1 → keeping all significant rows")
}

res_sig[, logHR := log2(HR)]

say("[INFO] logHR summary (all significant):")
print(summary(res_sig$logHR))

# For plotting, cap extreme values so boxes are not blown up
res_sig[, logHR_plot := pmin(pmax(logHR, ylim_lo), ylim_hi)]

# Pre-compute counts per cancer/mutation combo for labeling
box_counts <- res_sig[, .(n_signif = .N), by = .(cancer, mut_group)]
box_counts[, label := sprintf("n_sig=%d", n_signif)]
box_counts[, label_y := ylim_hi - 0.15]

# ---------------------------------------------------------------
# 3) BOXPLOT — log2(HR), colored by mut_group
# ---------------------------------------------------------------
p_box <- ggplot(
  res_sig,
  aes(x = cancer,
      y = logHR_plot,
      color = mut_group,
      group = interaction(cancer, mut_group))
) +
  geom_boxplot(
    fill = "white",
    width = 0.7,
    outlier.size = 0.4,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#c90028") +
  geom_text(
    data = box_counts,
    aes(x = cancer, y = 2, label = label),
    position = position_dodge(width = 0.8),
    size = 3,
    show.legend = FALSE,
    color = "#000000",
    angle = 90
  ) +
  scale_y_continuous(
    limits = c(ylim_lo, ylim_hi),
    breaks = seq(ylim_lo, ylim_hi, by = 1),
    labels = function(x) sprintf("HR=%.2f×", 2^x)
  ) +
  scale_color_manual(values = mut_colors_sig, name = "Mutation Group") +
  labs(
    title = "Significant Hazard Ratios per Cancer (Mutation-only CoxPH)",
    subtitle = "log2(HR) scale; values outside [-1.5, 10] are capped for visibility \n Filtered for ≤ 2 significant features per box",
    x = "Cancer Type",
    y = "Hazard Ratio (log2 scale)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "right",
    panel.grid.major.y = element_line(color = "grey85")
  )

ggsave(
  file.path(plot_dir, "m2_hr_boxplot_mutation_log2_zoom.png"),
  p_box, width = 16, height = 8, dpi = 300
)

say("[DONE] Saved mutation log2(HR) boxplot")

# ---------------------------------------------------------------
# 4) INTERPRETABLE MUTATION HR PLOT (Jitter only)
# ---------------------------------------------------------------
p_jitter <- ggplot(
  res_sig,
  aes(x = cancer,
      y = logHR,
      color = mut_group)
) +
  geom_jitter(
    position = position_jitter(width = 0.25, height = 0),
    size = 1.0,
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "#c90028",
    linewidth = 0.7
  ) +
  scale_y_continuous(
    limits = c(ylim_lo, ylim_hi),
    breaks = seq(ylim_lo, ylim_hi, by = 1),
    labels = function(x) sprintf("HR=%.2f×", 2^x)
  ) +
  scale_color_manual(values = mut_colors_sig, name = "Mutation Group") +
  labs(
    title    = "Significant Hazard Ratios (Mutation-only CoxPH)",
    subtitle = "log2(HR) scale",
    x = "Cancer Type",
    y = "Hazard Ratio (log2 scale)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major.y = element_line(color = "grey85"),
    legend.position = "right"
  )

ggsave(
  file.path(plot_dir, "m2_hr_jitter_mutation_log2_zoom.png"),
  p_jitter, width = 16, height = 8, dpi = 300
)

say("[DONE] Saved mutation jitter HR plot")

# ---------------------------------------------------------------
# 5) BARPLOT — Significant mutation feature counts
# ---------------------------------------------------------------
sig_counts <- res_all_sig[FDR < 0.05,
                          .(n_signif = .N),
                          by = .(cancer, mut_group)]

p_counts <- ggplot(sig_counts,
                   aes(x = cancer, y = n_signif + 1, fill = mut_group)) +
  geom_col(position = "dodge") +
  scale_y_log10() +
  scale_fill_manual(values = mut_colors_sig) +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Significant Mutation Features per Cancer",
    y = "Count (log10 scale)"
  )

ggsave(file.path(plot_dir, "m2_significant_mutation_barplot.png"),
       p_counts, width = 18, height = 10, dpi = 300)

# ---------------------------------------------------------------
# 6) P-VALUE HISTOGRAMS — ONLY if full results are requested
# ---------------------------------------------------------------
if (USE_FULL) {
  say("[INFO] Loading FULL *_full.csv files for p-value histograms...")
  files_full <- list.files(
    root,
    pattern = "cox_results_full\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(files_full) == 0) {
    say("[WARN] No *_full.csv files found; skipping p-value histograms.")
  } else {
    res_list_full <- lapply(files_full, load_mut_file)
    res_list_full <- Filter(Negate(is.null), res_list_full)
    res_all_full  <- rbindlist(res_list_full, fill = TRUE)
    say("[INFO] Loaded %d rows from full mutation result files",
        nrow(res_all_full))

    fwrite(res_all_full,
           file.path(plot_dir, "mutation_full_all_results_for_phist.csv"))

    pvals <- res_all_full[p >= 0 & p <= 1]

    # Colors based on all groups present in FULL results
    mut_levels_full <- sort(unique(pvals$mut_group))
    mut_colors_full <- scales::hue_pal()(length(mut_levels_full))
    names(mut_colors_full) <- mut_levels_full

    p_hist <- ggplot(pvals, aes(x = p, fill = mut_group)) +
      geom_histogram(bins = 40, position = "identity", alpha = 0.8) +
      facet_grid(mut_group ~ cancer, scales = "free_y") +
      scale_fill_manual(values = mut_colors_full, name = "Mutation group") +
      theme_bw(10) +
      labs(
        title = "P-value Distributions (Mutation-only, FULL models)",
        x = "p-value",
        y = "Count"
      )

    ggsave(file.path(plot_dir, "m2_pvalue_histograms.png"),
           p_hist, width = 22, height = 12, dpi = 300)
    say("[DONE] Saved p-value histograms → m2_pvalue_histograms.png")
  }
} else {
  say("[INFO] Skipping p-value histograms (set M2_USE_FULL_RESULTS=true in your plot launcher to include them).")
}

say(sprintf("[DONE] Mutation model (M2) plotting complete → %s", plot_dir))

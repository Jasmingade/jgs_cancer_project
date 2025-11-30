#!/usr/bin/env Rscript
# 05b_plot_m1_expression.R
# Model 1 – Expression-only CoxPH (03a)
# --------------------------------------
# Uses:
#   - Significant outputs:  TCGA_<CANCER>_<DTYPE>.cox_results.csv
#     (for all HR-based plots and diagnostics)
#   - Optional FULL outputs: TCGA_<CANCER>_<DTYPE>.cox_results_full.csv
#     (ONLY for p-value histograms, when M1_USE_FULL_RESULTS=true)
#
# Env toggle:
#   M1_USE_FULL_RESULTS = "true"/"1"/"yes" → load *_full.csv additionally for p-hist
#   otherwise                              → skip p-hist (everything else still works)
#
# Produces:
#   - m1_hr_boxplot_clean.png
#   - m1_hr_boxplot_faceted.png
#   - m1_logHR_distribution_stats.csv
#   - m1_sample_event_counts.csv
#   - m1_logHR_stats_with_events.csv
#   - m1_iso_frac_filter_flags.csv
#   - m1_diagnostic_iso_frac.png
#   - (optional) m1_pvalue_histograms.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(dplyr)
  library(rdist)
  library(grid)
  library(ggtext)
})

say     <- function(...) message(sprintf(...))
to_bool <- function(x) tolower(x) %in% c("1","true","yes")

root     <- "01_transcriptomics/out/03a_univariate_coxph"
plot_dir <- "01_transcriptomics/out/05_plots/model1_expr"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

USE_FULL <- to_bool(Sys.getenv("M1_USE_FULL_RESULTS", "false"))
say("[INFO] M1_USE_FULL_RESULTS (for p-hist) = %s", USE_FULL)
say("[INFO] Loading expression Cox results (SIGNIFICANT files) from: %s", root)

# ------------------------------------------------------------
# 0) Always load SIGNIFICANT ONLY results for main plots
#     (TCGA_<CANCER>_<DTYPE>.cox_results.csv)
# ------------------------------------------------------------
files_sig <- list.files(
  root,
  pattern = "cox_results\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
# This pattern does NOT match *_full.csv, so it's safe.

if (length(files_sig) == 0) {
  stop("No significant cox_results.csv files found in 03a_univariate_coxph")
}

say("[INFO] Found %d sig-only Cox result files", length(files_sig))

# Loader for 03a outputs (works for both sig and full, but here we use it on sig)
load_expr_file <- function(f){
  dt <- fread(f)
  needed <- c("feature","beta","HR","p")
  if (!all(needed %in% names(dt))) {
    warning(sprintf("[WARN] Skipping %s – missing required columns", f))
    return(NULL)
  }

  dt[, p  := as.numeric(p)]
  dt[, HR := as.numeric(HR)]

  # FDR
  if (!"FDR" %in% names(dt)) {
    dt[, FDR := p.adjust(p, method = "BH")]
  } else {
    dt[, FDR := as.numeric(FDR)]
  }

  fname <- basename(f)

  # Only create cancer/data_type if not already in file
  if (!"cancer" %in% names(dt)) {
    dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", fname)]
  }
  if (!"data_type" %in% names(dt)) {
    dt[, data_type :=
         fifelse(grepl("_gene", fname),     "gene",
         fifelse(grepl("_iso_log", fname),  "iso_log",
         fifelse(grepl("_iso_frac", fname), "iso_frac", "unknown")))]
  }

  # Keep valid rows
  dt <- dt[
    is.finite(HR) & HR > 0 &
    is.finite(p)  & p >= 0 & p <= 1
  ]

  dt
}

res_list_sig <- lapply(files_sig, load_expr_file)
res_list_sig <- Filter(Negate(is.null), res_list_sig)

res_all_sig <- rbindlist(res_list_sig, fill = TRUE)
say("[INFO] Loaded %d total rows (significant files combined)", nrow(res_all_sig))

# Save the combined significant results used in this script
fwrite(res_all_sig, file.path(plot_dir, "expression_sig_all_results.csv"))

type_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00"
)

# ------------------------------------------------------------
# 1) HR in log space (for all significant rows)
# ------------------------------------------------------------
res_sig <- res_all_sig[FDR < 0.05 & is.finite(HR) & HR > 0]

if (nrow(res_sig) == 0) {
  stop("No significant features (FDR<0.05) found for plotting.")
}

trim_lo <- 0.2
trim_hi <- 10

res_sig[, HR_trim := pmax(pmin(HR, trim_hi), trim_lo)]
res_sig[, logHR   := log2(HR_trim)]

# ------------------------------------------------------------
# 1b) Per cancer × datatype distribution summary
#     (log2(HR) stats + sample/event counts) – used for filtering
# ------------------------------------------------------------
dist_stats <- res_sig[
  ,
  .(
    n_sig        = .N,
    median_logHR = median(logHR, na.rm = TRUE),
    q1_logHR     = quantile(logHR, 0.25, na.rm = TRUE),
    q3_logHR     = quantile(logHR, 0.75, na.rm = TRUE),
    IQR_logHR    = IQR(logHR, na.rm = TRUE),
    min_logHR    = min(logHR, na.rm = TRUE),
    max_logHR    = max(logHR, na.rm = TRUE)
  ),
  by = .(cancer, data_type)
]

fwrite(dist_stats,
       file.path(plot_dir, "m1_logHR_distribution_stats.csv"))

# --- sample + event counts per cancer × data_type ---
norm_dir <- "01_transcriptomics/out/02_norm"

mani_files <- list.files(
  norm_dir,
  pattern = "sample_manifest\\.csv$",
  full.names = TRUE
)

mani_list <- lapply(mani_files, function(f) {
  mani <- fread(f)
  cancer    <- sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", basename(f))
  data_type <- sub("^TCGA_[A-Z0-9]+_(.*?)\\.sample_manifest\\.csv$", "\\1", basename(f))

  data.table(
    cancer    = cancer,
    data_type = data_type,
    n_samples = nrow(mani),
    n_events  = sum(mani$OS_event == 1, na.rm = TRUE)
  )
})

mani_stats <- rbindlist(mani_list, fill = TRUE)
fwrite(mani_stats,
       file.path(plot_dir, "m1_sample_event_counts.csv"))

# --- merge distribution stats with event counts ---
dist_stats_ev <- merge(
  dist_stats,
  mani_stats,
  by = c("cancer", "data_type"),
  all.x = TRUE
)

fwrite(dist_stats_ev,
       file.path(plot_dir, "m1_logHR_stats_with_events.csv"))

# ------------------------------------------------------------
# 1c) Define iso_frac filter (used for plotting)
#      keep_iso_frac = (n_events >= 40) & (IQR_logHR <= 2)
# ------------------------------------------------------------
iso_filter <- dist_stats_ev[
  data_type == "iso_frac",
  .(
    cancer,
    data_type,
    n_events,
    IQR_logHR,
    keep_iso_frac = !is.na(n_events) & n_events >= 10 &
                    !is.na(IQR_logHR) & IQR_logHR <= 4
  )
]

fwrite(iso_filter,
       file.path(plot_dir, "m1_iso_frac_filter_flags.csv"))

bad_iso_cancers  <- iso_filter[is.na(keep_iso_frac) | keep_iso_frac == FALSE, cancer]
good_iso_cancers <- iso_filter[keep_iso_frac == TRUE, cancer]

say("[INFO] iso_frac cancers kept for plotting (n_events>=10 & IQR<=4): %s",
    if (length(good_iso_cancers)) paste(sort(unique(good_iso_cancers)), collapse = ", ") else "<none>")
say("[INFO] iso_frac cancers DROPPED from plotting: %s",
    if (length(bad_iso_cancers)) paste(sort(unique(bad_iso_cancers)), collapse = ", ") else "<none>")

# Filter res_sig for plotting:
res_plot <- copy(res_sig)
if (length(bad_iso_cancers) > 0) {
  res_plot <- res_plot[
    !(data_type == "iso_frac" & cancer %in% bad_iso_cancers)
  ]
}
if (nrow(res_plot) == 0) {
  stop("After iso_frac filtering there are no rows left to plot.")
}

# y-range floor based on filtered data
y_floor <- min(res_plot$logHR, na.rm = TRUE) - 0.5

# ------------------------------------------------------------
# 1d) Per cancer × data_type stats for boxplot annotations
#      (on the filtered data)
# ------------------------------------------------------------
sig_counts <- res_plot[, .(n_signif = .N), by = .(cancer, data_type)]

fmt_or_q <- function(x) ifelse(is.na(x), "?", x)
first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_integer_ else x[1]
}

# attach n_samples/n_events from mani_stats
combo_dt <- unique(res_plot[, .(cancer, data_type)])
mani_stats_annot <- merge(combo_dt, mani_stats,
                          by = c("cancer", "data_type"), all.x = TRUE)

annot_base <- merge(sig_counts, mani_stats_annot,
                    by = c("cancer", "data_type"), all.x = TRUE)

annot_heights <- res_plot[
  , .(y_pos = max(logHR, na.rm = TRUE) + 0.3),
  by = .(cancer, data_type)
]

annotation_dt <- merge(annot_base, annot_heights,
                       by = c("cancer", "data_type"), all.x = TRUE)

annotation_dt[, label := sprintf("n=%s | events=%s | sig=%s",
                                 fmt_or_q(n_samples),
                                 fmt_or_q(n_events),
                                 fmt_or_q(n_signif))]
visible_boxes <- unique(res_plot[, .(cancer, data_type)])
annotation_dt_valid <- annotation_dt[
  visible_boxes,
  on = .(cancer, data_type),
  nomatch = 0L   # drop combos that are not drawn
][is.finite(y_pos)]
annotation_dt_valid[, signif_label := fmt_or_q(n_signif)]

y_cap <- if (nrow(annotation_dt_valid)) {
  max(annotation_dt_valid$y_pos, na.rm = TRUE) + 0.2
} else {
  max(res_plot$logHR, na.rm = TRUE) + 0.5
}

# axis labels per cancer (aggregated over datatypes)
cancer_axis_labels_dt <- mani_stats[
  ,
  .(
    n_samples = first_non_na(n_samples),
    n_events  = first_non_na(n_events)
  ),
  by = cancer
]

cancer_axis_labels <- setNames(
  sprintf(
    "%s\nn=%s | events=%s",
    cancer_axis_labels_dt$cancer,
    fmt_or_q(cancer_axis_labels_dt$n_samples),
    fmt_or_q(cancer_axis_labels_dt$n_events)
  ),
  cancer_axis_labels_dt$cancer
)

# ------------------------------------------------------------
# 2) HR boxplot (filtered, trimmed log2(HR))
# ------------------------------------------------------------
p <- ggplot(
  res_plot,
  aes(x = cancer, y = logHR,
      color = data_type,
      group = interaction(cancer, data_type))
) +
  geom_boxplot(
    fill = "white",
    outlier.shape = NA,
    width = 0.65,
    position = position_dodge(width = 0.8)
  ) +
  geom_hline(
    yintercept = 0, linetype = "dashed",
    color = "#c90028", linewidth = 0.6
  ) +
  geom_text(
    data = annotation_dt_valid,
    aes(
      x = cancer,
      y = y_floor,
      label = sprintf("n_sig = %s", signif_label),
      group = data_type
    ),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0.2,
    hjust = 0,
    position = position_dodge(width = 0.8),
    show.legend = FALSE,
    color = "#000000",
    angle = 90
  ) +
  scale_y_continuous(
    breaks = seq(-3, 3, 1),
    labels = function(x) sprintf("HR = %.2f×", 2^x)
  ) +
  scale_x_discrete(
    labels = function(x) {
      lbl <- cancer_axis_labels[x]
      lbl[is.na(lbl)] <- x[is.na(lbl)]
      lbl
    }
  ) +
  scale_color_manual(values = type_colors, name = "Data Type") +
  labs(
    title = "Significant Hazard Ratios per Cancer (FDR < 0.05)",
    subtitle = "Expression CoxPH — log2(HR), trimmed to HR 0.2–10\niso_frac cancers filtered by n_events ≥ 40 and IQR(log2(HR)) ≤ 4",
    x = "Cancer Type",
    y = "Hazard Ratio (log2 scale)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major.y = element_line(color = "grey85"),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(y_floor, y_cap))

outfile <- file.path(plot_dir, "m1_hr_boxplot_clean.png")
ggsave(outfile, p, width = 16, height = 8, dpi = 300)
say("[DONE] Saved clean M1 HR boxplot → %s", outfile)

# ------------------------------------------------------------
# 3) Optional faceted view (one panel per data_type)
# ------------------------------------------------------------
p_faceted <- ggplot(
  res_sig,
  aes(x = cancer, y = logHR,
      color = data_type)
) +
  geom_boxplot(
    fill = "white",
    outlier.size = 0.2,
    width = 0.6
  ) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "#c90028", linewidth = 0.5) +
  scale_y_continuous(
    breaks = seq(-3, 3, 1),
    labels = function(x) sprintf("HR = %.2f×", 2^x)
  ) +
  scale_color_manual(values = type_colors, guide = "none") +
  labs(
    title = "Significant Hazard Ratios per Cancer (Faceted by Data Type)",
    subtitle = "Expression CoxPH — log2(HR), trimmed to HR 0.2–10",
    x = "Cancer Type",
    y = "Hazard Ratio (log2 scale)"
  ) +
  geom_text(
    data = annotation_dt_valid,
    aes(x = cancer, y = y_pos, label = signif_label),
    inherit.aes = FALSE,
    size = 2.3,
    vjust = -0.2,
    color = "#000000"
  ) +
  scale_x_discrete(
    labels = function(x) {
      lbl <- cancer_axis_labels[x]
      lbl[is.na(lbl)] <- x[is.na(lbl)]
      lbl
    }
  ) +
  facet_wrap(~data_type, ncol = 1, scales = "fixed") +
  theme_bw(11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.y = element_line(color = "grey85")
  )

ggsave(file.path(plot_dir, "m1_hr_boxplot_faceted.png"),
       p_faceted, width = 14, height = 12, dpi = 300)

# ============================================================
# 4) Iso_frac diagnostic: IQR(log2(HR)) vs #events (UNfiltered)
# ============================================================

iso_diag <- dist_stats_ev[data_type == "iso_frac" & !is.na(IQR_logHR)]

if (nrow(iso_diag) > 0) {
  q <- ggplot(
        iso_diag,
        aes(x = n_events, y = IQR_logHR, label = cancer)
      ) +
    geom_point(color = "#E69F00", size = 2) +
    geom_text(nudge_y = 0.05, size = 2) +
    theme_bw(base_size = 12) +
    labs(
      title = "Iso_frac variability vs number of events",
      subtitle = "IQR of log2(HR) for significant iso_frac features (M1) – used for filtering",
      x = "Number of OS events per cancer",
      y = "IQR of log2(HR)"
    )

  ggsave(file.path(plot_dir, "m1_diagnostic_iso_frac.png"),
         q, width = 8, height = 6, dpi = 300)
  say("[DONE] Saved iso_frac diagnostic plot → m1_diagnostic_iso_frac.png")
} else {
  say("[WARN] No iso_frac stats available for diagnostic plot (no significant iso_frac or missing IQR).")
}

# ============================================================
# 5) P-value histograms — ONLY if full results are requested
# ============================================================
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
    res_list_full <- lapply(files_full, load_expr_file)
    res_list_full <- Filter(Negate(is.null), res_list_full)
    res_all_full  <- rbindlist(res_list_full, fill = TRUE)
    say("[INFO] Loaded %d rows from full result files", nrow(res_all_full))

    fwrite(res_all_full,
           file.path(plot_dir, "expression_full_all_results_for_phist.csv"))

    pvals <- res_all_full[p >= 0 & p <= 1]

    p3 <- ggplot(pvals, aes(x = p, fill = data_type)) +
      geom_histogram(bins = 40, position = "identity", alpha = 0.8) +
      facet_grid(data_type ~ cancer, scales = "free_y") +
      scale_fill_manual(values = type_colors, name = "Data type") +
      theme_bw(10) +
      labs(
        title = "P-value Distributions (Expression data)",
        x = "p-value",
        y = "Count"
      )

    ggsave(file.path(plot_dir, "m1_pvalue_histograms.png"),
           p3, width = 22, height = 12, dpi = 300)
    say("[DONE] Saved p-value histograms → m1_pvalue_histograms.png")
  }
} else {
  say("[INFO] Skipping p-value histograms (set M1_USE_FULL_RESULTS=true in '01_transcriptomics/pipeline/plotting/04_plot_launcher.sh' to include them).")
}

say(sprintf("[DONE] Expression model (M1) plotting complete → %s", plot_dir))

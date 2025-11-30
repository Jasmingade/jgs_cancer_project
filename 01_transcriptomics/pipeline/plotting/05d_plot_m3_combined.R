#!/usr/bin/env Rscript
# 05d_plot_m3_combined.R
# ------------------------------------------------------------
# Plotting for combined Expression + Mutation CoxPH (Model 3)
#
# Uses:
#   - SIGNIFICANT outputs:
#       01_transcriptomics/out/03c_exp_mutation_univariate_coxph/
#         TCGA_<CANCER>_<DTYPE>_<MUTGROUP>.cox_results.csv
#     for:
#       * Expression-term boxplots (log2(HR_expr))
#       * Mutation-term jitter plots (log2(HR_mut))
#       * Expression / mutation significant counts
#       * ALL-modality grouped boxplot
#         (iso_log / iso_frac optionally collapsed to ENSG)
#
#   - Optional FULL outputs:
#       TCGA_<CANCER>_<DTYPE>_<MUTGROUP>.cox_results_full.csv
#     ONLY for:
#       * P-value histograms (when M3_USE_FULL_RESULTS = true)
#
# Env toggle:
#   M3_USE_FULL_RESULTS = "true"/"1"/"yes" → load *_full.csv additionally for
#                                            p-value histograms
#   otherwise                               → skip p-value histograms
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(dplyr)
})

say         <- function(...) message(sprintf(...))
to_bool     <- function(x) tolower(x) %in% c("1","true","yes")
strip_version <- function(x) sub("\\.\\d+$", "", x)

# ------------------------------------------------------------
# Optional: isoform (ENST) -> gene (ENSG) mapping for collapsing
# ------------------------------------------------------------
tx2gene_file <- "01_transcriptomics/data/raw/tx2gene.csv"
tx2gene <- NULL

if (file.exists(tx2gene_file)) {
  tx2gene <- fread(tx2gene_file)
  names(tx2gene) <- tolower(names(tx2gene))

  tx_col <- intersect(
    c("tx_id", "transcript_id", "ensembl_transcript_id"),
    names(tx2gene)
  )[1]
  gene_col <- intersect(
    c("gene_id", "ensembl_gene_id", "gene"),
    names(tx2gene)
  )[1]

  if (is.na(tx_col) || is.na(gene_col)) {
    warning("[WARN] tx2gene file found but no suitable tx/gene columns — isoforms will NOT be collapsed.")
    tx2gene <- NULL
  } else {
    tx2gene <- unique(tx2gene[, .(
      tx_id   = strip_version(get(tx_col)),
      gene_id = strip_version(get(gene_col))
    )])
    say("[INFO] Loaded tx2gene mapping: %d transcript→gene rows", nrow(tx2gene))
  }
} else {
  warning(sprintf("[WARN] tx2gene mapping not found at %s — isoforms will NOT be collapsed.", tx2gene_file))
}

# ------------------------------------------------------------
# Paths + env toggle
# ------------------------------------------------------------
root     <- "01_transcriptomics/out/03c_exp_mutation_univariate_coxph"
plot_dir <- "01_transcriptomics/out/05_plots/model3_combined"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

USE_FULL <- to_bool(Sys.getenv("M3_USE_FULL_RESULTS", "false"))
say("[INFO] M3_USE_FULL_RESULTS (for p-hist) = %s", USE_FULL)

MIN_SIG_PER_GROUP <- as.integer(Sys.getenv("M3_MIN_SIG_PER_GROUP", "3"))
if (is.na(MIN_SIG_PER_GROUP) || MIN_SIG_PER_GROUP < 1) {
  MIN_SIG_PER_GROUP <- 1
}
say("[INFO] Minimum significant features per cancer/group = %d", MIN_SIG_PER_GROUP)

say("[INFO] Loading combined (expr + mut) Cox results (SIGNIFICANT files) from: %s", root)

# ------------------------------------------------------------
# 0) Always load SIGNIFICANT ONLY results for main plots
# ------------------------------------------------------------
files_sig <- list.files(
  root,
  pattern    = "cox_results\\.csv$",
  full.names = TRUE,
  recursive  = TRUE
)

if (length(files_sig) == 0)
  stop("No significant combined cox_results.csv files found in 03c_exp_mutation_univariate_coxph")

say("[INFO] Found %d sig-only combined result files", length(files_sig))

# ------------------------------------------------------------
# Loader for combined Expression + Mutation results
# ------------------------------------------------------------
load_combined_file <- function(f) {
  dt <- fread(f, showProgress = FALSE)

  needed <- c("feature",
              "beta_expr","HR_expr","p_expr",
              "beta_mut", "HR_mut", "p_mut")
  if (!all(needed %in% names(dt))) {
    warning(sprintf("[WARN] Skipping %s — missing required columns", f))
    return(NULL)
  }

  # Numeric conversions
  dt[, p_expr  := as.numeric(p_expr)]
  dt[, p_mut   := as.numeric(p_mut)]
  dt[, HR_expr := as.numeric(HR_expr)]
  dt[, HR_mut  := as.numeric(HR_mut)]

  # FDR computation if missing
  if (!"FDR_expr" %in% names(dt)) {
    dt[, FDR_expr := p.adjust(p_expr, "BH")]
  } else {
    dt[, FDR_expr := as.numeric(FDR_expr)]
  }

  if (!"FDR_mut" %in% names(dt)) {
    dt[, FDR_mut := p.adjust(p_mut,  "BH")]
  } else {
    dt[, FDR_mut := as.numeric(FDR_mut)]
  }

  fname <- basename(f)
  base  <- sub("\\.cox_results.*$", "", fname)

  # expected: TCGA_<CANCER>_(gene|iso_log|iso_frac)_(MUT_GROUP)
  m  <- regexec("^TCGA_([A-Z0-9]+)_(gene|iso_log|iso_frac)_(.+)$", base)
  mm <- regmatches(base, m)[[1]]

  if (length(mm) == 0) {
    warning(sprintf("[WARN] Unexpected filename pattern (cannot parse cancer/DTYPE/mut_group): %s", fname))
    return(NULL)
  }

  cancer    <- mm[2]
  data_type <- mm[3]
  mut_group <- mm[4]

  dt[, cancer    := cancer]
  dt[, data_type := data_type]
  dt[, mut_group := mut_group]

  # Keep clean rows
  dt <- dt[
    is.finite(HR_expr) & HR_expr > 0 &
    is.finite(HR_mut)  & HR_mut  > 0 &
    is.finite(p_expr)  & p_expr >= 0 & p_expr <= 1 &
    is.finite(p_mut)   & p_mut  >= 0 & p_mut  <= 1
  ]

  dt
}

# ------------------------------------------------------------
# Load all SIG-ONLY combined results
# ------------------------------------------------------------
res_list_sig <- lapply(files_sig, load_combined_file)
res_list_sig <- Filter(Negate(is.null), res_list_sig)

if (length(res_list_sig) == 0)
  stop("No valid combined result tables loaded from sig-only files")

res_all_sig <- rbindlist(res_list_sig, fill = TRUE)
say("[INFO] Loaded %d rows (combined model, sig-only files)", nrow(res_all_sig))

# Save combined sig-only results used in this script
fwrite(res_all_sig, file.path(plot_dir, "combined_sig_all_results.csv"))

# ------------------------------------------------------------
# Color palettes
# ------------------------------------------------------------
expr_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00")

mut_groups <- sort(unique(res_all_sig$mut_group))
# Palette for plots where color = mut_group
mut_colors <- setNames(
  scales::hue_pal()(length(mut_groups)),
  mut_groups)

# ============================================================
# 1. Expression-term boxplot (filtered log2(HR_expr))
# ============================================================
res_expr_sig <- res_all_sig[
  FDR_expr < 0.05 &
  HR_expr > 0 & is.finite(HR_expr)
]

if (nrow(res_expr_sig) > 0) {
  say("[INFO] Expression: %d significant terms (combined model)", nrow(res_expr_sig))

  # keep only "reasonable" HRs (edit thresholds as you like)
  hr_lo_expr <- 0.25   # HR >= 0.25×
  hr_hi_expr <- 8      # HR <= 8×

  res_expr_sig <- res_expr_sig[
    HR_expr >= hr_lo_expr & HR_expr <= hr_hi_expr
  ]

  res_expr_sig[, logHR_expr := log2(HR_expr)]

  expr_box_counts <- res_expr_sig[, .(n_signif = .N), by = .(cancer, data_type)]
  expr_box_counts[, label := sprintf("n_sig=%d", n_signif)]
  expr_box_counts[, label_y := log2(hr_hi_expr) - 0.15]

  p_expr <- ggplot(
    res_expr_sig,
    aes(x = cancer,
        y = logHR_expr,
        color = data_type,
        group = interaction(cancer, data_type))
    ) +
    geom_boxplot(
      fill = "white",
      width = 0.6,
      outlier.shape = NA,
      position = position_dodge(width = 0.75)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "#c90028", linewidth = 0.6) +
    geom_text(
      data = expr_box_counts,
      aes(x = cancer, y = label_y, label = label,
          group = interaction(cancer, data_type)),
      position = position_dodge(width = 0.75),
      size = 3,
      show.legend = FALSE,
      color = "#000000",
      angle = 90
    ) +
    scale_y_continuous(
      breaks = seq(-3, 3, 1),
      labels = function(x) sprintf("HR=%.2f×", 2^x)
    ) +
    scale_color_manual(values = expr_colors, name = "Expression Type") +
    labs(
      title = "Expression Effects in Combined Model (M3)",
      subtitle = "Significant expr terms (FDR_expr<0.05), log2(HR) filtered",
      y = "Expression log2(HR)",
      x = "Cancer Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey85")
    )

  ggsave(file.path(plot_dir, "m3_expr_log2HR_boxplot.png"),
         p_expr, width = 16, height = 8, dpi = 300)
}

# ============================================================
# 2. Mutation-term boxplot (filtered log2(HR_mut))
# ============================================================
res_mut_sig <- res_all_sig[
  FDR_mut < 0.05 &
  HR_mut > 0 & is.finite(HR_mut)
]

if (nrow(res_mut_sig) > 0) {
  say("[INFO] Mutation: %d significant terms (combined model)", nrow(res_mut_sig))

  hr_lo_mut <- 0.25   # HR >= 0.25× <-no
  hr_hi_mut <- 10     # HR <= 10×

  res_mut_sig <- res_mut_sig[
    HR_mut >= hr_lo_mut & HR_mut <= hr_hi_mut
  ]

  res_mut_sig[, logHR_mut := log2(HR_mut)]

  mut_box_counts <- res_mut_sig[, .(n_signif = .N), by = .(cancer, mut_group)]
  mut_box_counts[, label := sprintf("n_sig=%d", n_signif)]
  mut_box_counts[, label_y := log2(hr_hi_mut) - 0.15]

  p_mut <- ggplot(
    res_mut_sig,
    aes(x = cancer,
        y = logHR_mut,
        color = mut_group,
        group = interaction(cancer, mut_group))
    ) +
    geom_boxplot(
      fill = "white",
      width = 0.6,
      outlier.size = 0.6,
      position = position_dodge(width = 0.75)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#c90028") +
    geom_text(
      data = mut_box_counts,
      aes(x = cancer, y = label_y, label = label,
          group = interaction(cancer, mut_group)),
      position = position_dodge(width = 0.75),
      size = 3,
      show.legend = FALSE,
      color = "#000000",
      angle = 90
    ) +
    scale_y_continuous(
      limits = c(log2(hr_lo_mut), log2(hr_hi_mut)),
      breaks = seq(log2(hr_lo_mut), log2(hr_hi_mut), by = 1),
      labels = function(x) sprintf("HR=%.2f×", 2^x)
    ) +
    scale_color_manual(values = mut_colors, name = "Mutation Type") +
    labs(
      title = "Mutation Effects in Combined Model (M3)",
      subtitle = "Significant mutation terms (FDR_mut<0.05), log2(HR) filtered",
      y = "Mutation log2(HR)",
      x = "Cancer Type"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey85")
    )

  ggsave(file.path(plot_dir, "m3_mut_log2HR_boxplot.png"),
         p_mut, width = 16, height = 8, dpi = 300)
}

# ============================================================
# 3. Significant counts per cancer
# ============================================================
expr_counts <- res_all_sig[
  , .(n_signif = sum(FDR_expr < 0.05, na.rm = TRUE)),
  by = .(cancer, data_type)
]

p_expr_count <- ggplot(expr_counts,
                       aes(x = cancer, y = n_signif + 1, fill = data_type)) +
  geom_col(position = "dodge") +
  scale_y_log10() +
  scale_fill_manual(values = expr_colors) +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Significant Expression Features (Combined M3)",
       y = "log10(count+1)")

ggsave(file.path(plot_dir, "m3_expr_significant_barplot.png"),
       p_expr_count, width = 18, height = 8, dpi = 300)

mut_counts <- res_all_sig[
  , .(n_signif = sum(FDR_mut < 0.05, na.rm = TRUE)),
  by = .(cancer, mut_group)
]

p_mut_count <- ggplot(mut_counts,
                      aes(x = cancer, y = n_signif + 1, fill = mut_group)) +
  geom_col(position = "dodge") +
  scale_y_log10() +
  scale_fill_manual(values = setNames(scales::hue_pal()(length(mut_groups)), mut_groups)) +
  theme_bw(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Significant Mutation Features (Combined M3)",
       y = "log10(count+1)")

ggsave(file.path(plot_dir, "m3_mut_significant_barplot.png"),
       p_mut_count, width = 18, height = 8, dpi = 300)

# ============================================================
# 4. P-value histograms — ONLY if full results are requested
# ============================================================
if (USE_FULL) {
  say("[INFO] Loading FULL *_full.csv files for p-value histograms...")
  files_full <- list.files(
    root,
    pattern    = "cox_results_full\\.csv$",
    full.names = TRUE,
    recursive  = TRUE
  )

  if (length(files_full) == 0) {
    say("[WARN] No *_full.csv files found; skipping p-value histograms.")
  } else {
    res_list_full <- lapply(files_full, load_combined_file)
    res_list_full <- Filter(Negate(is.null), res_list_full)
    res_all_full  <- rbindlist(res_list_full, fill = TRUE)
    say("[INFO] Loaded %d rows from full combined result files", nrow(res_all_full))

    fwrite(res_all_full,
           file.path(plot_dir, "combined_full_all_results_for_phist.csv"))

    p_long <- rbind(
      res_all_full[, .(cancer, term = "expr", p = p_expr)],
      res_all_full[, .(cancer, term = "mut",  p = p_mut)]
    )

    p_long <- p_long[is.finite(p) & p >= 0 & p <= 1]

    p_hist <- ggplot(p_long, aes(x = p)) +
      geom_histogram(bins = 40, fill = "grey80") +
      facet_grid(term ~ cancer, scales = "free_y") +
      theme_bw(9) +
      labs(title = "P-value Distributions (Combined Model M3, FULL results)",
           x = "p-value")

    ggsave(file.path(plot_dir, "m3_pvalue_histograms_expr_vs_mut.png"),
           p_hist, width = 22, height = 10, dpi = 300)
    say("[DONE] Saved p-value histograms → m3_pvalue_histograms_expr_vs_mut.png")
  }
} else {
  say("[INFO] Skipping p-value histograms (set M3_USE_FULL_RESULTS=true in plot launcher to include them).")
}

# ============================================================
# 5. ALL modalities combined grouped boxplot
#     (isoforms optionally collapsed to ENSG using tx2gene)
# ============================================================
say("[INFO] Building ALL-modality combined boxplot...")

expr_sig <- res_all_sig[
  FDR_expr < 0.05 & HR_expr > 0 & is.finite(HR_expr)
]
group_counts <- expr_sig[, .(n_signif = .N), by = .(cancer, data_type)]
if (MIN_SIG_PER_GROUP > 1) {
  keep_pairs <- group_counts[n_signif >= MIN_SIG_PER_GROUP,
                             .(cancer, data_type)]
  expr_sig <- expr_sig[keep_pairs, on = .(cancer, data_type)]
  if (nrow(expr_sig) == 0) {
    stop(sprintf(
      "No cancer/mutation groups have at least %d significant features.",
      MIN_SIG_PER_GROUP
    ))
  }
  say("[INFO] After requiring >=%d sig features per cancer/group for expression data: %d rows kept",
      MIN_SIG_PER_GROUP, nrow(expr_sig))
} else {
  say("[INFO] Minimum per cancer/group set to 1 → keeping all significant rows")
}


mut_sig <- res_all_sig[
  FDR_mut < 0.05 & HR_mut > 0 & is.finite(HR_mut)
]
group_counts <- mut_sig[, .(n_signif = .N), by = .(cancer, mut_group)]
if (MIN_SIG_PER_GROUP > 1) {
  keep_pairs <- group_counts[n_signif >= MIN_SIG_PER_GROUP,
                             .(cancer, mut_group)]
  mut_sig <- mut_sig[keep_pairs, on = .(cancer, mut_group)]
  if (nrow(mut_sig) == 0) {
    stop(sprintf(
      "No cancer/mutation groups have at least %d significant features.",
      MIN_SIG_PER_GROUP
    ))
  }
  say("[INFO] After requiring >=%d sig features per cancer/group for mutation data: %d rows kept",
      MIN_SIG_PER_GROUP, nrow(mut_sig))
} else {
  say("[INFO] Minimum per cancer/group set to 1 → keeping all significant rows")
}


if (nrow(expr_sig) > 0 || nrow(mut_sig) > 0) {

  # sensible HR ranges (edit if you like)
  hr_lo_expr <- 0.25   # 0.25×
  hr_hi_expr <- 8      # 8×
  hr_lo_mut  <- 0.25   # 0.25×
  hr_hi_mut  <- 16     # 16×

  # ---- Expression: filter + log2(HR_expr) ----
  expr_sig <- expr_sig[
    HR_expr >= hr_lo_expr & HR_expr <= hr_hi_expr
  ]
  expr_sig[, logHR := log2(HR_expr)]
  expr_sig[, modality := data_type]

  # (collapse isoforms as you already do)
  if (!is.null(tx2gene)) {
    expr_sig[, feature_clean := strip_version(feature)]
    expr_sig[, gene_id := feature_clean]
    iso_idx <- expr_sig$data_type %in% c("iso_log", "iso_frac")
    if (any(iso_idx)) {
      expr_sig[iso_idx, gene_id :=
                tx2gene[match(feature_clean, tx_id), gene_id]]
    }
    expr_sig <- expr_sig[!is.na(gene_id)]
    expr_collapsed <- expr_sig[
      , .(logHR = median(logHR, na.rm = TRUE)),
      by = .(cancer, modality, gene_id)
    ]
  } else {
    expr_collapsed <- expr_sig[, .(cancer, modality, logHR)]
  }

  # ---- Mutation: filter + log2(HR_mut) ----
  mut_sig <- mut_sig[
    HR_mut >= hr_lo_mut & HR_mut <= hr_hi_mut
  ]
  mut_sig[, logHR := log2(HR_mut)]
  mut_sig[, modality := mut_group]

  # -------------------------
  # Combine & plot
  # -------------------------
  combined <- rbindlist(list(
    expr_collapsed[, .(cancer, modality, logHR)],
    mut_sig[,        .(cancer, modality, logHR)]
  ), fill = TRUE)

  mut_groups <- sort(unique(res_all_sig$mut_group))  # put this before using it

  combined[, modality := factor(
    modality,
    levels = c("gene", "iso_log", "iso_frac", mut_groups)
  )]

  # Y-axis limits in *log2(HR)* space
  log_lo <- log2(hr_lo_expr) # -2
  log_hi <- log2(hr_hi_mut)  # 4

  say("[DEBUG] Combined logHR summary:")
  print(summary(combined$logHR))
  say("[DEBUG] Combined counts per modality:")
  print(table(combined$modality))

  combined_counts <- combined[, .(n_signif = .N), by = .(cancer, modality)]
  combined_counts[, label := sprintf("n_sig=%d", n_signif)]
  combined_counts[, label_y := log_hi - 0.15]

  all_colors <- c(expr_colors, mut_colors)

  p_all <- ggplot(
    combined,
    aes(x = cancer, y = logHR, color = modality,
        group = interaction(cancer, modality))
  ) +
    geom_boxplot(
      fill = "white",
      width = 0.55,
      outlier.shape = NA,
      position = position_dodge(width = 0.75)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#c90028") +
    geom_text(
      data = combined_counts,
      aes(x = cancer, y = label_y, label = label,
          group = interaction(cancer, modality)),
      position = position_dodge(width = 0.75),
      size = 4,
      show.legend = FALSE,
      color = "#000000",
      angle = 90
    ) +
    scale_y_continuous(
      limits = c(log_lo, log_hi),
      breaks = seq(log_lo, log_hi, by = 1),
      labels = function(x) sprintf("HR=%.2f×", 2^x)
    ) +
    scale_color_manual(values = all_colors, name = "Modality") +
    labs(
      title    = "log2(HR): ALL Modalities Across All Cancers (M3)",
      subtitle = "Expression (isoforms collapsed to ENSG when mapping available) + Mutation; HR trimmed in HR space",
      x        = "Cancer Type",
      y        = "log2(HR)"
    ) +
    theme_bw(12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "grey85")
    )

  ggsave(
    file.path(plot_dir, "m3_combined_log2HR_boxplot_ALL_modalities_ENSG.png"),
    p_all, width = 25, height = 12, dpi = 300
  )
}

say(sprintf("[DONE] Combined Model 3 plotting complete → %s", plot_dir))

#!/usr/bin/env Rscript
# 05e_plot_m4_interactions.R
# Model 4 – Isoform × Mutation interaction (03d)
# ----------------------------------------------
# Loads per-cancer, per-datatype, per-mut_group interaction CoxPH results from:
#   01_transcriptomics/out/03d_iso_mut_univariate_coxph
#
# Assumes each .csv contains at least:
#   feature, HR, p
# plus optional FDR. If FDR missing → computed from p.
#
# Produces:
#   - Boxplot of significant interaction HRs per cancer × data_type
#   - Heatmap of top features (feature × cancer) log2(HR)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
})

say <- function(...) message(sprintf(...))
strip_version <- function(x) sub("\\.\\d+$", "", x)

root    <- "01_transcriptomics/out/03d_iso_mut_univariate_coxph"
out_dir <- "01_transcriptomics/out/05_plots/model4_interactions"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

say("[INFO] Loading interaction Cox results from: %s", root)

# Prefer *_full outputs if present
files_full <- list.files(
  root,
  pattern = "cox_results_full\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(files_full) > 0) {
  files <- files_full
  say("[INFO] Using FULL interaction results (_full.csv): %d files", length(files))
} else {
  files <- list.files(
    root,
    pattern = "cox_results\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  say("[INFO] No *_full files found, using cox_results.csv: %d files", length(files))
}

if (length(files) == 0) stop("No interaction CoxPH results found!")

# ------------------------------------------------------------
# Loader for interaction results
# ------------------------------------------------------------
load_int_file <- function(f) {
  dt <- fread(f, showProgress = FALSE)

  # Map from HR_int / p_int if necessary
  if ("HR_int" %in% names(dt) && !"HR" %in% names(dt)) {
    dt[, HR := HR_int]
  }
  if ("p_int" %in% names(dt) && !"p" %in% names(dt)) {
    dt[, p := p_int]
  }

  needed <- c("feature", "HR", "p")
  if (!all(needed %in% names(dt))) {
    warning(sprintf("[WARN] Skipping %s — missing required columns", f))
    return(NULL)
  }

  dt[, HR := as.numeric(HR)]
  dt[, p  := as.numeric(p)]

  if (!"FDR" %in% names(dt)) {
    dt[, FDR := p.adjust(p, "BH")]
  } else {
    dt[, FDR := as.numeric(FDR)]
  }

  fname <- basename(f)
  base  <- sub("\\.cox_results.*$", "", fname)

  # expected: TCGA_<CANCER>_(gene|iso_log|iso_frac)_<MUTGROUP>
  m  <- regexec("^TCGA_([A-Z0-9]+)_(gene|iso_log|iso_frac)_(.+)$", base)
  mm <- regmatches(base, m)[[1]]

  if (length(mm) == 0) {
    warning(sprintf("[WARN] Cannot parse cancer/data_type/mut_group from %s", fname))
    return(NULL)
  }

  cancer    <- mm[2]
  data_type <- mm[3]
  mut_group <- mm[4]

  dt[, cancer    := cancer]
  dt[, data_type := data_type]
  dt[, mut_group := mut_group]

  dt <- dt[
    is.finite(HR) & HR > 0 &
    is.finite(p)  & p >= 0 & p <= 1
  ]

  dt
}

res_list <- lapply(files, load_int_file)
res_list <- Filter(Negate(is.null), res_list)

if (length(res_list) == 0)
  stop("No valid interaction tables loaded")

int_all <- rbindlist(res_list, fill = TRUE)
say("[INFO] Loaded %d rows (interaction model)", nrow(int_all))

# Strip Ensembl versions for convenience (ENSGxxxx.x → ENSGxxxx)
int_all[, feature_clean := strip_version(feature)]

fwrite(int_all, file.path(out_dir, "m4_interaction_all_results.csv"))

# ------------------------------------------------------------
# Color palettes
# ------------------------------------------------------------
expr_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00"
)

mut_groups <- sort(unique(int_all$mut_group))
mut_colors <- scales::hue_pal()(length(mut_groups))

# ------------------------------------------------------------
# 1) Boxplot – mutation groups, faceted by expression data_type
# ------------------------------------------------------------
int_sig <- int_all[FDR < 0.05 & is.finite(HR) & HR > 0]

if (nrow(int_sig) == 0) {
  say("[INFO] No significant interactions (FDR<0.05) – boxplot skipped.")
} else {
  int_sig[, logHR := log2(HR)]

  ## nice ordering for factors
  int_sig[, cancer    := factor(cancer, levels = sort(unique(cancer)))]
  int_sig[, data_type := factor(data_type, levels = c("gene","iso_frac","iso_log"))]
  int_sig[, mut_group := factor(mut_group,
                                levels = c("missense_or_inframe","rna","splice","truncating_LOF"))]

  ## per-box n_sig labels (smaller and a bit above the boxes)
  box_counts_mut <- int_sig[
    , .(n_signif = .N), by = .(cancer, mut_group, data_type)
  ]
  box_counts_mut[, label   := sprintf("n=%d", n_signif)]
  box_counts_mut[, label_y := 2.5]

  p_box <- ggplot(
    int_sig,
    aes(x = cancer, y = logHR,
        color = mut_group,
        group = interaction(cancer, mut_group))
  ) +
    geom_boxplot(
      fill = "white",
      width = 0.65,
      outlier.shape = NA,
      position = position_dodge(width = 0.75),
      linewidth = 0.3
    ) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color   = "#c90028",
               linewidth = 0.4) +
    geom_text(
      data = box_counts_mut,
      aes(x = cancer, y = label_y, label = label,
          group = interaction(cancer, mut_group)),
      position   = position_dodge(width = 0.75),
      size       = 2.3,
      show.legend = FALSE,
      color      = "black",
      angle      = 90,
      vjust      = 0
    ) +
    facet_wrap(~ data_type, nrow = 1) +
    scale_y_continuous(
      limits = c(-3, 3),
      breaks = -3:3,
      labels = function(x) sprintf("HR=%.2f×", 2^x)
    ) +
    scale_color_manual(
      values = setNames(mut_colors, mut_groups),
      name   = "Mutation group"
    ) +
    labs(
      title    = "Isoform × Mutation Interaction Effects Across Cancers",
      subtitle = "Significant interaction terms (FDR < 0.05), log2(HR)",
      x        = "Cancer",
      y        = "Interaction HR (log2 scale)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x        = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      legend.position    = "right",
      strip.background   = element_rect(fill = "grey92", colour = NA),
      strip.text         = element_text(face = "bold")
    )

  out_png1 <- file.path(out_dir, "m4_mut_per_datatype.png")
  ggsave(out_png1, p_box, width = 20, height = 10, dpi = 300)
  say("[DONE] Saved interaction HR boxplot → %s", out_png1)
}


# ------------------------------------------------------------
# 2) Boxplot – expression data types, faceted by mutation group
# ------------------------------------------------------------
if (nrow(int_sig) == 0) {
  say("[INFO] No significant interactions (FDR<0.05) – second boxplot skipped.")
} else {
  box_counts_dt <- int_sig[
    , .(n_signif = .N), by = .(cancer, data_type, mut_group)
  ]
  box_counts_dt[, label   := sprintf("n=%d", n_signif)]
  box_counts_dt[, label_y := 2.5]

  p_mut <- ggplot(
    int_sig,
    aes(x = cancer, y = logHR,
        color = data_type,
        group = interaction(cancer, data_type))
  ) +
    geom_boxplot(
      fill = "white",
      width = 0.65,
      outlier.shape = NA,
      position = position_dodge(width = 0.75),
      linewidth = 0.3
    ) +
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color   = "#c90028",
               linewidth = 0.4) +
    geom_text(
      data = box_counts_dt,
      aes(x = cancer, y = label_y, label = label,
          group = interaction(cancer, data_type)),
      position   = position_dodge(width = 0.75),
      size       = 2.3,
      show.legend = FALSE,
      color      = "black",
      angle      = 90,
      vjust      = 0
    ) +
    facet_wrap(~ mut_group, nrow = 2) +
    scale_y_continuous(
      limits = c(-3, 3),
      breaks = -3:3,
      labels = function(x) sprintf("HR=%.2f×", 2^x)
    ) +
    scale_color_manual(
      values = expr_colors,
      name   = "Expression type"
    ) +
    labs(
      title    = "Isoform × Mutation Interaction Effects Across Cancers",
      subtitle = "Significant interaction terms (FDR < 0.05), log2(HR)",
      x        = "Cancer",
      y        = "Interaction HR (log2 scale)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x        = element_text(angle = 60, hjust = 1, vjust = 1, size = 7),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      legend.position    = "right",
      strip.background   = element_rect(fill = "grey92", colour = NA),
      strip.text         = element_text(face = "bold")
    )

  out_png2 <- file.path(out_dir, "m4_datatype_per_mut.png")
  ggsave(out_png2, p_mut, width = 12, height = 6, dpi = 300)
  say("[DONE] Saved interaction HR boxplot → %s", out_png2)
}



# ------------------------------------------------------------
# 2) Heatmap of top interaction features (feature × cancer)
# ------------------------------------------------------------
if (nrow(int_sig) == 0) {
  say("[INFO] No significant interactions for heatmap.")
} else {
  # define effect strength and pick top features (using feature_clean)
  int_sig[, logHR := log2(HR)]
  int_sig[, effect_strength := abs(logHR)]

  top_n <- 40L

  # number of distinct (version-stripped) features
  n_feats <- uniqueN(int_sig$feature_clean)
  top_n_use <- min(top_n, n_feats)

  # top features by effect_strength
  top_features <- int_sig[
    order(-effect_strength),
    head(unique(feature_clean), top_n_use)
  ]

  hm_dt <- int_sig[feature_clean %in% top_features,
                   .(logHR = median(logHR, na.rm = TRUE)),
                   by = .(feature_clean, cancer)]

  # cast to matrix: feature_clean × cancer
  hm_wide <- dcast(hm_dt, feature_clean ~ cancer, value.var = "logHR")
  mat <- as.matrix(hm_wide[, -1, with = FALSE])
  rownames(mat) <- hm_wide$feature_clean

  # ---- CLEAN MATRIX FOR CLUSTERING ----
  # drop rows/cols that are all NA or single non-NA (no distance info)
  keep_row <- rowSums(!is.na(mat)) > 1
  keep_col <- colSums(!is.na(mat)) > 1
  mat <- mat[keep_row, keep_col, drop = FALSE]

  # if anything left, replace remaining NA with 0 (no effect)
  if (length(mat)) {
    mat[is.na(mat)] <- 0
  } else {
    stop("Heatmap matrix is empty after filtering rows/cols with only NA.")
  }

  # colour function
  max_abs <- max(abs(mat), na.rm = TRUE)
  max_abs <- max(max_abs, 0.5)
  col_fun <- colorRamp2(
    c(-max_abs, 0, max_abs),
    c("#377EB8", "white", "#E41A1C")
  )

  out_png2 <- file.path(out_dir, "m4_interaction_heatmap.png")
  png(out_png2, width = 14, height = 10, units = "in", res = 300)
  ht <- Heatmap(
    mat,
    name = "log2(HR)",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(
      title  = "log2(HR)",
      at     = c(-max_abs, 0, max_abs),
      labels = sprintf("%.2f", c(-max_abs, 0, max_abs))
    )
  )
  draw(ht)
  dev.off()
  say("[DONE] Saved interaction heatmap → %s", out_png2)
}

say("[DONE] M4 interaction plotting complete → %s", out_dir)

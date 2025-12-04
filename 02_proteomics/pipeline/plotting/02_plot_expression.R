#!/usr/bin/env Rscript

# Plot significant model outputs across methods (cox, penalized, glmnet).
# Inputs:
#   <cox_root> <out_dir>
# Expects:
#   - Univariate Cox:      *cox_results.csv (with HR)
#   - Penalized (coxmos):  *.penalized.csv (beta)
#   - glmnet Cox:          *.glmnet.csv (beta)
# Writes per-method outputs under <out_dir>/<method>/:
#   - expression_distribution_stats.csv
#   - expression_boxplot.png
#   - expression_median_iqr.png

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
})

say <- function(...) message(sprintf(...))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 02_plot_expression.R <cox_root> <out_dir>")
}
cox_root <- args[[1]]
out_dir  <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_meta <- function(stem) {
  parts <- strsplit(stem, "\\.")[[1]]
  dataset <- parts[1]
  dtype <- if (length(parts) >= 2) parts[2] else "unknown"
  sample_type <- if (grepl("_normal$", dataset, ignore.case = TRUE)) "normal" else
                 if (grepl("_tumor$", dataset, ignore.case = TRUE)) "tumor" else "tumor"
  study <- sub("_(TMT|ITRAQ|LABELFREE).*", "", dataset)
  list(dataset = dataset, dtype = dtype, sample_type = sample_type, study_id = study)
}

type_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00",
  unknown  = "grey50"
)

load_files_for_method <- function(files, method) {
  rbindlist(lapply(files, function(f) {
    dt <- fread(f)
    meta <- parse_meta(tools::file_path_sans_ext(basename(f)))
    if (method == "cox") {
      if (!"HR" %in% names(dt)) return(NULL)
      dt[, HR := as.numeric(HR)]
      dt <- dt[is.finite(HR) & HR > 0]
      if (!nrow(dt)) return(NULL)
      dt[, value := log2(HR)]
      dt[, value_label := "log2(HR)"]
    } else if (method %in% c("penalized", "glmnet")) {
      if (!"beta" %in% names(dt)) return(NULL)
      dt[, beta := as.numeric(beta)]
      dt <- dt[is.finite(beta)]
      if (!nrow(dt)) return(NULL)
      dt[, value := beta]
      dt[, value_label := "beta"]
    } else {
      return(NULL)
    }
    dt[, c("dataset","data_type","sample_type","study_id","method") :=
         list(meta$dataset, meta$dtype, meta$sample_type, meta$study_id, method)]
    dt
  }), fill = TRUE)
}

plot_method <- function(dt, method_name, cancer_map) {
  if (!nrow(dt)) {
    say("[info] %s: no usable rows", method_name)
    return()
  }
  dt[is.na(sample_type), sample_type := "tumor"]
  dt[is.na(data_type), data_type := "unknown"]
  dt <- dt[is.finite(value)]
  if (!nrow(dt)) {
    say("[warn] %s: all values non-finite", method_name)
    return()
  }

  dt[, cancer := {
    if (!is.null(cancer_map) && length(cancer_map)) {
      mapped <- vapply(study_id, function(s) {
        if (!is.null(cancer_map[[s]])) paste0(as.character(cancer_map[[s]]), ":", s) else s
      }, character(1))
      mapped
    } else {
      study_id
    }
  }]

  dist_stats <- dt[
    ,
    .(
      n_sig   = .N,
      med_val = median(value, na.rm = TRUE),
      iqr_val = IQR(value, na.rm = TRUE),
      min_val = min(value, na.rm = TRUE),
      max_val = max(value, na.rm = TRUE),
      value_label = first(value_label)
    ),
    by = .(dataset, study_id, data_type, sample_type, cancer)
  ]

  method_out <- file.path(out_dir, method_name)
  dir.create(method_out, recursive = TRUE, showWarnings = FALSE)
  fwrite(dist_stats, file.path(method_out, "expression_distribution_stats.csv"))
  say("[write] %s/expression_distribution_stats.csv", method_name)

  trim_lo <- quantile(dt$value, 0.01, na.rm = TRUE)
  trim_hi <- quantile(dt$value, 0.99, na.rm = TRUE)
  dt[, value_trim := pmax(pmin(value, trim_hi), trim_lo)]

  dataset_order <- unique(dist_stats[order(med_val)][, dataset])
  dt[, dataset := factor(dataset, levels = dataset_order)]
  dist_stats[, dataset := factor(dataset, levels = dataset_order)]
  dt[, cancer := factor(cancer)]
  cancer_order <- dt[, .(med = median(value_trim, na.rm = TRUE)), by = cancer][order(med)]$cancer
  dt[, cancer := factor(cancer, levels = cancer_order)]

  y_lab <- unique(dt$value_label)
  if (length(y_lab) != 1) y_lab <- method_name

  y_floor <- min(dt$value_trim, na.rm = TRUE) - 0.5
  y_cap <- max(dt$value_trim, na.rm = TRUE) + 0.5

  p_box <- ggplot(
    dt,
    aes(
      x = cancer,
      y = value_trim,
      color = data_type,
      group = interaction(cancer, data_type)
    )
  ) +
    geom_boxplot(
      fill = "white",
      outlier.shape = NA,
      width = 0.65,
      position = position_dodge(width = 0.75)
    ) +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      color = "#c90028", linewidth = 0.6
    ) +
    facet_grid(~ sample_type, scales = "free_y") +
    scale_color_manual(values = type_colors, name = "Data Type") +
    labs(
      title = sprintf("Significant results per cancer (%s)", toupper(method_name)),
      subtitle = sprintf("Trimmed to 1–99%% (%s)", y_lab),
      x = "Cancer Type",
      y = y_lab
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.y = element_line(color = "grey85"),
      legend.position = "right"
    ) +
    coord_cartesian(ylim = c(y_floor, y_cap))

  p_scatter <- ggplot(dist_stats,
                      aes(x = med_val, y = iqr_val,
                          color = data_type, shape = sample_type)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = type_colors) +
    labs(
      title = sprintf("Dataset median vs IQR (%s)", toupper(method_name)),
      x = sprintf("Median %s", y_lab),
      y = sprintf("IQR %s", y_lab)
    ) +
    theme_bw()

  ggsave(file.path(method_out, "expression_boxplot.png"),
         p_box, width = 12, height = 8, dpi = 300)
  say("[write] %s/expression_boxplot.png", method_name)

  ggsave(file.path(method_out, "expression_median_iqr.png"),
         p_scatter, width = 7, height = 5, dpi = 300)
  say("[write] %s/expression_median_iqr.png", method_name)
}

methods <- list(
  list(name = "cox",       pattern = "cox_results\\.csv$", exclude = "_full\\.csv$", loader = function(files) load_files_for_method(files, "cox")),
  list(name = "penalized", pattern = "\\.penalized\\.csv$", exclude = NULL,         loader = function(files) load_files_for_method(files, "penalized")),
  list(name = "glmnet",    pattern = "\\.glmnet\\.csv$",    exclude = NULL,         loader = function(files) load_files_for_method(files, "glmnet"))
)

cancer_cfg <- tryCatch(yaml::read_yaml("02_proteomics/config/cancers.yaml"),
                       error = function(e) NULL)
cancer_map <- if (!is.null(cancer_cfg) && "cancers" %in% names(cancer_cfg)) cancer_cfg$cancers else list()

for (m in methods) {
  files <- list.files(cox_root, pattern = m$pattern, full.names = TRUE, recursive = TRUE)
  if (!length(files)) {
    say("[info] %s: no files matched", m$name)
    next
  }
  if (!is.null(m$exclude)) files <- files[!grepl(m$exclude, files)]
  if (!length(files)) {
    say("[info] %s: files excluded by pattern", m$name)
    next
  }
  dt <- m$loader(files)
  if (!is.null(dt) && nrow(dt)) {
    plot_method(dt, m$name, cancer_map)
  } else {
    say("[info] %s: no usable rows after loading", m$name)
  }
}

message("[done] Wrote plots and stats to ", out_dir)

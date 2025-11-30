#!/usr/bin/env Rscript

# Basic sanity plots for proteomics preprocessing.
# Produces median intensity density plots and missingness heatmaps per PDC study,
# split by tumour/normal sample types.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 01_plot_sanity.R <filtered_expr_root> <out_dir>")
}

expr_root <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_matrix <- function(path) {
  dt <- fread(path)
  if (ncol(dt) < 2) return(NULL)
  feat <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- feat
  mat
}

parse_dataset <- function(file, dtype) {
  stem <- tools::file_path_sans_ext(basename(file))
  stem <- sub(paste0("_", dtype, "$"), "", stem, ignore.case = TRUE)
  tokens <- strsplit(stem, "_")[[1]]
  study <- tokens[1]
  platform <- if (length(tokens) >= 2) tokens[2] else NA_character_
  sample_type <- NA_character_
  last_token <- tokens[length(tokens)]
  if (last_token %in% c("tumor", "normal")) {
    sample_type <- last_token
  }
  data.table(
    dtype = dtype,
    study = study,
    platform = platform,
    sample_type = sample_type,
    file = file
  )
}

collect_runs <- function(dtype) {
  base <- file.path(expr_root, dtype)
  if (!dir.exists(base)) return(data.table())
  files <- list.files(base, pattern = paste0("_", dtype, "\\.csv$"), full.names = TRUE)
  if (!length(files)) return(data.table())
  rbindlist(lapply(files, parse_dataset, dtype = dtype), fill = TRUE)
}

run_table <- rbindlist(list(collect_runs("gene"), collect_runs("iso_log")), fill = TRUE)
if (!nrow(run_table)) stop("No filtered expression files found under ", expr_root)
run_table[is.na(sample_type), sample_type := "all"]

plot_sanity <- function(sub_dt) {
  study <- unique(sub_dt$study)
  sample_type <- unique(sub_dt$sample_type)
  med_stats <- list()
  missing_stats <- list()

  for (dtype in unique(sub_dt$dtype)) {
    rows <- sub_dt[dtype == dtype]
    for (f in rows$file) {
      mat <- load_matrix(f)
      if (is.null(mat)) next
      med <- apply(mat, 2, median, na.rm = TRUE)
      med_stats[[length(med_stats) + 1L]] <- data.table(
        dtype = dtype,
        sample = colnames(mat),
        median = med
      )
      missing <- apply(mat, 2, function(x) mean(is.na(x)))
      missing_stats[[length(missing_stats) + 1L]] <- data.table(
        dtype = dtype,
        sample = colnames(mat),
        missing = missing
      )
    }
  }

  if (!length(med_stats) || !length(missing_stats)) return(NULL)
  med_dt <- rbindlist(med_stats)
  miss_dt <- rbindlist(missing_stats)

  p1 <- ggplot(med_dt, aes(x = median, fill = dtype)) +
    geom_density(alpha = 0.5) +
    theme_bw() +
    labs(title = sprintf("%s %s – Sample medians", study, sample_type),
         x = "Median intensity", y = "Density")

  p2 <- ggplot(miss_dt, aes(x = sample, y = missing, fill = dtype)) +
    geom_col(position = "dodge") +
    theme_bw() + coord_flip() +
    labs(title = sprintf("%s %s – Missing fraction", study, sample_type),
         x = "Sample", y = "Fraction missing")

  list(p1, p2)
}

for (study in unique(run_table$study)) {
  study_types <- unique(run_table[study == !!study, sample_type])
  for (stype in intersect(c("tumor", "normal", "all"), study_types)) {
    sub <- run_table[study == !!study & sample_type == !!stype]
    if (!nrow(sub)) next
    plots <- plot_sanity(sub)
    if (is.null(plots)) next
    out_file <- file.path(out_dir, sprintf("%s_%s_sanity.pdf", study, stype))
    pdf(out_file, width = 12, height = 6)
    grid.arrange(grobs = plots, ncol = 2)
    dev.off()
    message("[plot] wrote ", out_file)
  }
}

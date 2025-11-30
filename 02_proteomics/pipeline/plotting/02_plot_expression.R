#!/usr/bin/env Rscript

# Expression-focused plots for proteomics per study/sample-type.
# Summaries:
#   * Volcano-like scatter (median vs MAD) per run
#   * Boxplots of sample medians grouped by platform

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 02_plot_expression.R <filtered_expr_root> <out_dir>")
}

expr_root <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_matrix <- function(path) {
  dt <- fread(path)
  if (ncol(dt) < 2) return(NULL)
  if (names(dt)[1] == "") setnames(dt, 1, "feature")
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
    dataset = stem,
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

run_table <- collect_runs("gene")
if (!nrow(run_table)) stop("No gene matrices found under ", expr_root)
run_table[is.na(sample_type), sample_type := "all"]
run_table[is.na(platform), platform := "unknown"]

for (study in unique(run_table$study)) {
  study_types <- unique(run_table[study == !!study, sample_type])
  for (stype in intersect(c("tumor", "normal", "all"), study_types)) {
    sub <- run_table[study == !!study & sample_type == !!stype]
    if (!nrow(sub)) next

    med_summary <- list()
    run_summary <- list()

    for (i in seq_len(nrow(sub))) {
      row <- sub[i]
      mat <- load_matrix(row$file)
      if (is.null(mat)) next
      sample_medians <- apply(mat, 2, median, na.rm = TRUE)
      med_summary[[length(med_summary) + 1L]] <- data.table(
        sample = names(sample_medians),
        median = as.numeric(sample_medians),
        platform = row$platform,
        dataset = row$dataset
      )

      feat_med <- apply(mat, 1, median, na.rm = TRUE)
      feat_mad <- apply(mat, 1, mad, na.rm = TRUE)
      run_summary[[length(run_summary) + 1L]] <- data.table(
        feature = rownames(mat),
        median = feat_med,
        mad = feat_mad,
        dataset = row$dataset
      )
    }

    if (!length(med_summary) || !length(run_summary)) next

    med_dt <- rbindlist(med_summary, fill = TRUE)
    run_dt <- rbindlist(run_summary, fill = TRUE)

    p_box <- ggplot(med_dt, aes(x = platform, y = median, fill = platform)) +
      geom_boxplot() +
      theme_bw() +
      labs(title = sprintf("%s %s – Sample median distribution", study, stype),
           x = "Platform", y = "Sample median")

    p_scatter <- ggplot(run_dt, aes(x = median, y = mad, color = dataset)) +
      geom_point(alpha = 0.4) +
      theme_bw() +
      labs(title = sprintf("%s %s – Feature median vs MAD", study, stype),
           x = "Feature median", y = "MAD")

    out_file <- file.path(out_dir, sprintf("%s_%s_expression.pdf", study, stype))
    pdf(out_file, width = 12, height = 6)
    grid.arrange(p_box, p_scatter, ncol = 2)
    dev.off()
    message("[plot] Wrote ", out_file)
  }
}

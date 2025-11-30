#!/usr/bin/env Rscript

# Generate per-study normalization/batch-correction plots showing
# raw -> normalized -> batch-corrected distributions for gene, iso_log, iso_frac.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 || length(args) > 6) {
  stop("Usage: 04_plots.R <raw_root> <norm_root> <batch_root> <batch_annotation_dir> <out_dir> [dataset_id]")
}

raw_root   <- args[[1]]
norm_root  <- args[[2]]
batch_root <- args[[3]]
batch_dir  <- args[[4]]
out_dir    <- args[[5]]
target_dataset <- if (length(args) == 6) args[[6]] else NA_character_

say <- function(fmt, ...) cat(sprintf(paste0("[plot] ", fmt, "\n"), ...))

max_plot_features <- as.integer(Sys.getenv("PREPROCESS_PLOT_MAX_FEATURES", "4000"))

canonicalize_single <- function(val) {
  if (length(val) != 1L) val <- val[1]
  if (is.null(val) || is.na(val) || !nzchar(val)) return(NA_character_)
  clean <- gsub("[^a-zA-Z0-9_ ]", "", val)
  clean <- gsub(" +", "_", trimws(tolower(clean)))
  parts <- unlist(strsplit(clean, "_", fixed = TRUE))
  parts <- parts[parts != ""]
  if (!length(parts)) return(NA_character_)
  if (length(parts) >= 3) parts <- parts[seq_len(3)]
  paste(parts, collapse = "_")
}

canonicalize_batch <- function(vals) {
  if (is.null(vals)) return(rep(NA_character_, length(vals)))
  vapply(vals, canonicalize_single, character(1), USE.NAMES = FALSE)
}

infer_sample_type <- function(dataset) {
  if (grepl("_normal$", dataset, ignore.case = TRUE)) return("normal")
  if (grepl("_tumor$", dataset, ignore.case = TRUE)) return("tumor")
  NA_character_
}

load_case_batches <- function(study_key) {
  if (!nzchar(batch_dir)) return(NULL)
  file <- file.path(batch_dir, paste0(study_key, ".csv"))
  if (!file.exists(file)) {
    pattern <- paste0("^", study_key, "\\.csv$")
    cand <- list.files(batch_dir, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
    if (!length(cand)) return(NULL)
    file <- cand[1]
    say("Using batch annotation file %s for study key %s", basename(file), study_key)
  }
  dt <- tryCatch(fread(file), error = function(e) NULL)
  if (is.null(dt) || !"case_id" %in% names(dt) || !"folder_name" %in% names(dt)) {
    return(NULL)
  }
  dt[, case_id := toupper(trimws(case_id))]
  if ("sample_type" %in% names(dt)) {
    dt[, sample_type := canonicalize_batch(sample_type)]
  } else {
    dt[, sample_type := NA_character_]
  }
  dt[, folder_name := canonicalize_batch(folder_name)]
  dt <- dt[nzchar(case_id) & nzchar(folder_name)]
  if (!nrow(dt)) return(NULL)
  dt <- dt[order(case_id, sample_type, folder_name)]
  dt <- dt[!duplicated(data.table(case_id, sample_type))]
  dt
}

load_matrix <- function(path) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  if (ncol(dt) < 2) return(NULL)
  feat <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- feat
  storage.mode(mat) <- "double"
  rn <- colnames(mat)
  keep <- !grepl("POOLED|QC", rn, ignore.case = TRUE)
  mat[, keep, drop = FALSE]
}

gather_values <- function(mat, dtype_label, stage_label, order_idx, ordered_batches,
                          dataset_name, entries_list) {
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) return(entries_list)
  mat_use <- mat
  if (nrow(mat_use) > max_plot_features) {
    idx <- sample.int(nrow(mat_use), max_plot_features)
    mat_use <- mat_use[idx, , drop = FALSE]
  }
  mat_use <- mat_use[, order_idx, drop = FALSE]
  values <- as.vector(mat_use)
  sample_order <- rep(seq_along(order_idx), each = nrow(mat_use))
  batches <- rep(ordered_batches, each = nrow(mat_use))
  entries_list[[length(entries_list) + 1L]] <- data.table(
    dataset = dataset_name,
    dtype = dtype_label,
    stage = stage_label,
    value = values,
    sample_order = sample_order,
    batch = batches
  )
  entries_list
}

ensure_common_cols <- function(mats, reference_cols) {
  common <- reference_cols
  for (m in mats) {
    if (!is.null(m)) {
      common <- intersect(common, colnames(m))
    }
  }
  if (!length(common)) return(NULL)
  lapply(mats, function(m) {
    if (is.null(m)) return(NULL)
    m[, common, drop = FALSE]
  })
}

raw_gene_dir  <- file.path(raw_root, "gene")
raw_iso_dir   <- file.path(raw_root, "iso_log")
norm_gene_dir <- file.path(norm_root, "gene")
norm_iso_dir  <- file.path(norm_root, "iso_log")
batch_gene_dir <- file.path(batch_root, "gene")
batch_iso_dir  <- file.path(batch_root, "iso_log")

if (!dir.exists(raw_gene_dir)) stop("Missing raw gene directory: ", raw_gene_dir)

datasets <- list.files(raw_gene_dir, pattern = "_gene\\.csv$", full.names = FALSE)
datasets <- sub("_gene\\.csv$", "", datasets)
datasets <- datasets[!grepl("_reference$", datasets, ignore.case = TRUE)]
if (!is.na(target_dataset)) {
  datasets <- datasets[datasets == target_dataset]
}
if (!length(datasets)) stop("No datasets found for plotting.")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (dataset in sort(datasets)) {
  say("Processing dataset %s", dataset)
  study_key <- sub("_(normal|tumor)$", "", dataset, ignore.case = TRUE)
  dataset_type <- infer_sample_type(dataset)

  raw_gene   <- load_matrix(file.path(raw_gene_dir,  paste0(dataset, "_gene.csv")))
  norm_gene  <- load_matrix(file.path(norm_gene_dir, paste0(dataset, "_gene.csv")))
  batch_gene <- load_matrix(file.path(batch_gene_dir, paste0(dataset, "_gene.csv")))
  raw_iso    <- load_matrix(file.path(raw_iso_dir,   paste0(dataset, "_iso_log.csv")))
  norm_iso   <- load_matrix(file.path(norm_iso_dir,  paste0(dataset, "_iso_log.csv")))
  batch_iso  <- load_matrix(file.path(batch_iso_dir, paste0(dataset, "_iso_log.csv")))

  if (is.null(raw_gene) || is.null(norm_gene) || is.null(batch_gene)) {
    say("Skipping %s (missing gene matrices)", dataset)
    next
  }

  sample_ids <- colnames(raw_gene)
  case_batches <- load_case_batches(study_key)
  sample_batches <- rep(NA_character_, length(sample_ids))
  if (!is.null(case_batches)) {
    dt_sub <- case_batches
    if (!is.na(dataset_type) && "sample_type" %in% names(case_batches)) {
      sel <- dt_sub[sample_type == dataset_type]
      if (nrow(sel)) dt_sub <- sel
    }
    idx <- match(toupper(sample_ids), dt_sub$case_id)
    sample_batches <- dt_sub$folder_name[idx]
  }
  sample_batches[is.na(sample_batches)] <- "unknown"
  order_idx <- order(sample_batches, sample_ids)
  ordered_batches <- sample_batches[order_idx]

  gene_mats <- ensure_common_cols(list(raw_gene, norm_gene, batch_gene), sample_ids)
  if (is.null(gene_mats)) {
    say("No overlapping gene samples for %s", dataset)
    next
  }
  graw <- gene_mats[[1]]
  gnorm <- gene_mats[[2]]
  gbatch <- gene_mats[[3]]

  plot_entries <- list()
  plot_entries <- gather_values(graw, "gene", "raw", order_idx, ordered_batches, dataset, plot_entries)
  plot_entries <- gather_values(gnorm, "gene", "normalized", order_idx, ordered_batches, dataset, plot_entries)
  plot_entries <- gather_values(gbatch, "gene", "batch_corrected", order_idx, ordered_batches, dataset, plot_entries)

  iso_mats <- ensure_common_cols(list(raw_iso, norm_iso, batch_iso), sample_ids)
  if (!is.null(iso_mats)) {
    iraw <- iso_mats[[1]]
    inorm <- iso_mats[[2]]
    ibatch <- iso_mats[[3]]
    plot_entries <- gather_values(iraw, "iso_log", "raw", order_idx, ordered_batches, dataset, plot_entries)
    plot_entries <- gather_values(inorm, "iso_log", "normalized", order_idx, ordered_batches, dataset, plot_entries)
    plot_entries <- gather_values(ibatch, "iso_log", "batch_corrected", order_idx, ordered_batches, dataset, plot_entries)

    if (!is.null(iraw)) {
      raw_frac <- plogis(iraw)
      plot_entries <- gather_values(raw_frac, "iso_frac", "raw", order_idx, ordered_batches, dataset, plot_entries)
    }
    if (!is.null(inorm)) {
      norm_frac <- plogis(inorm)
      plot_entries <- gather_values(norm_frac, "iso_frac", "normalized", order_idx, ordered_batches, dataset, plot_entries)
    }
    if (!is.null(ibatch)) {
      batch_frac <- plogis(ibatch)
      plot_entries <- gather_values(batch_frac, "iso_frac", "batch_corrected", order_idx, ordered_batches, dataset, plot_entries)
    }
  }

  if (!length(plot_entries)) {
    say("No plot entries for %s", dataset)
    next
  }

  plot_dt <- rbindlist(plot_entries, fill = TRUE)
  plot_dt[, dtype := factor(dtype, levels = c("gene", "iso_log", "iso_frac"))]
  plot_dt <- plot_dt[!is.na(dtype)]
  plot_dt[, stage := factor(stage, levels = c("raw", "normalized", "batch_corrected"))]
  plot_dt[, batch := factor(batch)]

  p <- ggplot(plot_dt, aes(x = sample_order, y = value, fill = batch)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    facet_grid(dtype ~ stage, scales = "free_y") +
    labs(
      title = sprintf("Raw → normalized → batch-corrected: %s", dataset),
      x = "Sample order (grouped by batch)",
      y = "Expression intensity",
      fill = "Batch"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom"
    )

  out_pdf <- file.path(out_dir, sprintf("%s_preprocess_boxplots.pdf", dataset))
  out_png <- file.path(out_dir, sprintf("%s_preprocess_boxplots.png", dataset))
  ggsave(out_pdf, p, width = 10, height = 7, units = "in")
  ggsave(out_png, p, width = 10, height = 7, units = "in", dpi = 300)
  say("Saved plot for %s → %s", dataset, out_pdf)
}

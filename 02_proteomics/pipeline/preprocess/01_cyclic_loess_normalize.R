#!/usr/bin/env Rscript

# Normalize combined proteomics matrices stored as CSVs:
#   <data_root>/<dtype>/<dataset>_<dtype>.csv
# Outputs are written to <out_root>/<dtype>/ with matching filenames.
# iso_frac is regenerated from normalized iso_log.

suppressPackageStartupMessages({
  library(data.table)
})
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("Package 'limma' is required. Install via install.packages('limma').")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 4) {
  stop("Usage: 01_cyclic_loess_normalize.R <data_root> <out_root> [comma_separated_dtypes] [dataset_id]")
}

data_root <- args[[1]]
out_root <- args[[2]]
dtype_arg <- if (length(args) >= 3) args[[3]] else "gene,iso_log"
target_dataset <- if (length(args) == 4) args[[4]] else NA_character_

data_types <- unique(trimws(strsplit(dtype_arg, ",")[[1]]))
data_types <- data_types[nzchar(data_types)]
if (!length(data_types)) stop("No data types specified.")

if ("iso_frac" %in% data_types) {
  warning("[norm] 'iso_frac' is derived from iso_log; removing from requested types.")
  data_types <- setdiff(data_types, "iso_frac")
}

if (!dir.exists(data_root)) stop("Data root not found: ", data_root)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

say <- function(fmt, ...) cat(sprintf(paste0("[norm] ", fmt, "\n"), ...))

normalize_matrix <- function(mat, label) {
  if (!is.matrix(mat) || !nrow(mat) || !ncol(mat)) {
    warning(sprintf("[norm] %s has empty matrix, skipping", label))
    return(NULL)
  }
  if (ncol(mat) < 2) {
    warning(sprintf("[norm] %s has <2 samples; copying without normalization", label))
    return(mat)
  }
  limma::normalizeCyclicLoess(mat, method = "fast")
}

build_dt <- function(mat) {
  dt <- data.table(feature = rownames(mat))
  dt <- cbind(dt, as.data.table(mat, keep.rownames = FALSE))
  setnames(dt, c("feature", colnames(mat)))
  dt
}

iso_frac_name <- function(base_name) {
  if (grepl("_iso_log(\\.csv)?$", base_name, ignore.case = TRUE)) {
    sub("_iso_log(\\.csv)?$", "_iso_frac.csv", base_name, ignore.case = TRUE)
  } else if (tolower(base_name) == "iso_log.csv") {
    "iso_frac.csv"
  } else {
    sub("\\.csv$", "_iso_frac.csv", base_name, ignore.case = TRUE)
  }
}

for (dtype in intersect(c("gene","iso_log"), data_types)) {
  base_dir <- file.path(data_root, dtype)
  if (!dir.exists(base_dir)) {
    warning(sprintf("[norm] Missing directory for %s: %s", dtype, base_dir))
    next
  }

  files <- list.files(base_dir, pattern = paste0(dtype, "\\.csv$"), full.names = TRUE, recursive = FALSE)
  if (!is.na(target_dataset)) {
    files <- files[sub(paste0("_", dtype, "\\.csv$"), "", basename(files), ignore.case = TRUE) == target_dataset]
  }
  if (!length(files)) {
    warning(sprintf("[norm] No %s CSV files selected under %s", dtype, base_dir))
    next
  }

  say("Scanning %s datasets under %s", dtype, base_dir)

  dtype_out_dir <- file.path(out_root, dtype)
  dir.create(dtype_out_dir, recursive = TRUE, showWarnings = FALSE)

  for (expr_file in sort(files)) {
    dataset_name <- sub(paste0("_", dtype, "\\.csv$"), "", basename(expr_file), ignore.case = TRUE)
    say("[1/3] Normalizing dataset %s (%s)", dataset_name, dtype)
    expr_dt <- fread(expr_file)
    if (ncol(expr_dt) < 2) {
      warning(sprintf("[norm] %s has <2 columns, skipping", expr_file))
      next
    }
    feature_col <- names(expr_dt)[1]
    expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
    rownames(expr_mat) <- expr_dt[[feature_col]]
    storage.mode(expr_mat) <- "double"

    norm_mat <- normalize_matrix(expr_mat, dataset_name)
    if (is.null(norm_mat)) next

    norm_dt <- build_dt(norm_mat)
    out_file <- file.path(dtype_out_dir, basename(expr_file))
    fwrite(norm_dt, out_file)
    say("Wrote %s normalized → %s", dtype, out_file)

    if (dtype == "iso_log") {
      eps <- 1e-6
      frac_mat <- plogis(norm_mat)
      frac_mat <- pmin(pmax(frac_mat, eps), 1 - eps)
      frac_dt <- build_dt(frac_mat)
      frac_dir <- file.path(out_root, "iso_frac")
      dir.create(frac_dir, recursive = TRUE, showWarnings = FALSE)
      frac_out <- file.path(frac_dir, iso_frac_name(basename(expr_file)))
      fwrite(frac_dt, frac_out)
      say("   Derived iso_frac → %s", frac_out)
    }
  }
  say("Completed %s datasets (%d files processed)", dtype, length(files))
}

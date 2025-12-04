#!/usr/bin/env Rscript

# Log2-transform raw proteomics matrices (gene, iso_log) with a +1 pseudocount.
# Usage:
#   00_log_transform.R <input_root> <output_root> [comma_sep_dtypes] [target_dataset]

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 4) {
  stop("Usage: 00_log_transform.R <input_root> <output_root> [comma_sep_dtypes] [target_dataset]")
}

input_root  <- args[[1]]
output_root <- args[[2]]
dtype_arg   <- if (length(args) >= 3) args[[3]] else "gene,iso_log"
target_dataset <- if (length(args) == 4) args[[4]] else NA_character_

say <- function(fmt, ...) cat(sprintf(paste0("[log] ", fmt, "\n"), ...))

canonical_dataset_name <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(x)
  if (grepl("_reference$", x, ignore.case = TRUE)) return(x)
  sub("_(normal|tumor)$", "", x, ignore.case = TRUE)
}

matches_target <- function(file_base, target) {
  if (is.na(target) || !nzchar(target)) return(TRUE)
  canon <- canonical_dataset_name(file_base)
  target %in% c(file_base, canon)
}

dtypes <- unique(trimws(strsplit(dtype_arg, ",")[[1]]))
dtypes <- dtypes[nzchar(dtypes)]
if (!length(dtypes)) stop("[log] No dtypes specified.")

if (!dir.exists(input_root)) stop("[log] Input root not found: ", input_root)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

transform_file <- function(path, dtype) {
  dt <- fread(path)
  if (ncol(dt) < 2) {
    say("Skipping %s (%s) <2 columns", basename(path), dtype)
    return()
  }
  dataset_base <- sub(paste0("_", dtype, "\\.csv$"), "", basename(path), ignore.case = TRUE)
  dataset_canon <- canonical_dataset_name(dataset_base)
  feats <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  storage.mode(mat) <- "double"
  mat <- log2(mat + 1)
  out_dt <- data.table(feature = feats)
  out_dt <- cbind(out_dt, as.data.table(mat))
  out_dir <- file.path(output_root, dtype)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_name <- sprintf("%s_%s.csv", dataset_canon, dtype)
  out_path <- file.path(out_dir, out_name)
  fwrite(out_dt, out_path)
  if (!identical(dataset_base, dataset_canon)) {
    say("Transformed → %s (canonical name from %s)", out_path, dataset_base)
  } else {
    say("Transformed → %s", out_path)
  }
}

for (dtype in dtypes) {
  if (dtype == "iso_frac") {
    say("Skipping iso_frac in log transform")
    next
  }
  in_dir <- file.path(input_root, dtype)
  if (!dir.exists(in_dir)) {
    say("Missing input dir for %s: %s", dtype, in_dir)
    next
  }
  files <- list.files(in_dir, pattern = paste0(dtype, "\\.csv$"), full.names = TRUE, recursive = FALSE)
  if (!is.na(target_dataset)) {
    files <- Filter(function(f) {
      base <- sub(paste0("_", dtype, "\\.csv$"), "", basename(f), ignore.case = TRUE)
      matches_target(base, target_dataset)
    }, files)
  }
  if (!length(files)) {
    say("No %s files to transform", dtype)
    next
  }
  for (f in files) transform_file(f, dtype)
}

say("Done.")

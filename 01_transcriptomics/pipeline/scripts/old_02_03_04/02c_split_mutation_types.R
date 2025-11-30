#!/usr/bin/env Rscript
# ============================================================
# 02c_split_mutation_types.R
# Optional step: export each mutation type as its own covariate matrix
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: 02c_split_mutation_types.R <by_gene_type_csv> <out_dir>")

in_file <- args[[1]]
out_dir <- args[[2]]

mut <- fread(in_file)
stopifnot("case_id" %in% names(mut))

split_dir <- file.path(out_dir, "split_by_type")
dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)

cols <- setdiff(names(mut), "case_id")
types <- unique(sub(".*_", "", cols))
message("[INFO] Found ", length(types), " mutation types")

for (tp in types) {
  cols_tp <- c("case_id", grep(paste0("_", tp, "$"), cols, value = TRUE))
  if (length(cols_tp) < 2) next
  out_file <- file.path(split_dir, sprintf("mutation_covariates.%s.csv", tp))
  fwrite(mut[, ..cols_tp], out_file)
  message("[DONE] Wrote ", tp, ": ", out_file)
}

message("[DONE] Split-by-type matrices written to: ", split_dir)

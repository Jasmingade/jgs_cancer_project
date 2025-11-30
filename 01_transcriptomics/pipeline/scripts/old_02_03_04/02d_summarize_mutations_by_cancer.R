#!/usr/bin/env Rscript
# ============================================================
# 02d_summarize_mutations_by_cancer.R
# Summarizes mutation burden and per-type frequencies by cancer
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: 02d_summarize_mutations_by_cancer.R <mut_tsv> <out_csv>")

mut_tsv <- args[[1]]
out_csv <- args[[2]]

mut <- fread(mut_tsv)
stopifnot(all(c("case_id", "gene", "effect") %in% names(mut)))

summary_dt <- mut[, .(
  n_cases = uniqueN(case_id),
  n_genes = uniqueN(gene),
  mean_mutations_per_case = .N / uniqueN(case_id)
), by = .(effect)]

dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
fwrite(summary_dt, out_csv)
message("[DONE] Mutation summary saved to: ", out_csv)

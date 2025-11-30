#!/usr/bin/env Rscript
# ============================================================
# 02a_all_mutations.R
# Cleans and normalizes mutation data for downstream analyses
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: 02a_all_mutations.R <01_transcriptomics/data/mutation/mc3.v0.2.8.PUBLIC.xena.all_mutation_positions.gz> <out_tsv>")

mc3_gz  <- args[[1]]
out_tsv <- args[[2]]

message("[LOAD] Reading mutation file: ", mc3_gz)
mc3 <- fread(mc3_gz, sep = "\t", showProgress = TRUE)

# Standardize key columns
if ("Hugo_Symbol" %in% names(mc3)) setnames(mc3, "Hugo_Symbol", "gene")
if ("Tumor_Sample_Barcode" %in% names(mc3)) setnames(mc3, "Tumor_Sample_Barcode", "sample")
if ("Variant_Classification" %in% names(mc3)) setnames(mc3, "Variant_Classification", "effect")

need <- c("sample", "gene", "effect")
stopifnot(all(need %in% names(mc3)))

# Normalize IDs
mc3[, case_id := substr(gsub("\\.", "-", sample), 1, 12)]
mc3 <- mc3[, .(case_id, sample, gene, effect)]

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)
fwrite(mc3, out_tsv, sep = "\t")
message(sprintf("[DONE] Wrote cleaned mutation file: %s (rows=%d, genes=%d)",
                out_tsv, nrow(mc3), length(unique(mc3$gene))))

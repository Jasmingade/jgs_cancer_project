#!/usr/bin/env Rscript
# ============================================================
# 02b_build_mutation_covariates.R
# Builds per-gene AND per-gene-type binary mutation matrices (0/1)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: 02b_build_mutation_covariates.R <mut_tsv> <out_dir> [min_prev=0.05]")

mut_tsv  <- args[[1]]
out_dir  <- args[[2]]
min_prev <- ifelse(length(args) >= 3, as.numeric(args[[3]]), 0.05)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("[LOAD] Reading mutation table: ", mut_tsv)
mut <- fread(mut_tsv)
stopifnot(all(c("case_id", "gene", "effect") %in% names(mut)))
message("[INFO] Rows: ", nrow(mut), " | Genes: ", length(unique(mut$gene)))

# Normalize mutation types
canon_effect <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", trimws(x)))
  x <- sub("^missense_mutation$", "missense", x)
  x <- sub("^nonsense_mutation$", "nonsense", x)
  x <- sub("^frame_shift_del$", "frameshift_del", x)
  x <- sub("^frame_shift_ins$", "frameshift_ins", x)
  x <- sub("^in_frame_del$", "inframe_del", x)
  x <- sub("^in_frame_ins$", "inframe_ins", x)
  x <- sub("^splice_site$", "splice", x)
  x <- sub("^silent$", "silent", x)
  x
}
mut[, effect := canon_effect(effect)]

# ============================================================
# (1) PER-GENE baseline binary matrix
# ============================================================
message("[BUILD] Per-gene mutation presence")
gene_case <- unique(mut[, .(case_id, gene)])
wide_gene <- dcast(gene_case, case_id ~ gene, fun.aggregate = length)
for (j in setdiff(names(wide_gene), "case_id")) set(wide_gene, j = j, value = as.integer(wide_gene[[j]] > 0))
prev_gene <- sapply(setdiff(names(wide_gene), "case_id"), function(cn) mean(wide_gene[[cn]] == 1, na.rm = TRUE))
keep_gene <- names(prev_gene)[prev_gene >= min_prev]
wide_gene <- wide_gene[, c("case_id", keep_gene), with = FALSE]
message("[INFO] Retained ", length(keep_gene), " genes ≥ ", 100*min_prev, "% prevalence")

# ============================================================
# (2) PER-GENE-TYPE binary matrix (NEW)
# ============================================================
message("[BUILD] Per-gene+type matrix (e.g., TP53_missense)")
mut[, feature := paste0(gene, "_", effect)]
gene_case_eff <- unique(mut[, .(case_id, feature)])
wide_eff <- dcast(gene_case_eff, case_id ~ feature, fun.aggregate = length)
for (j in setdiff(names(wide_eff), "case_id")) set(wide_eff, j = j, value = as.integer(wide_eff[[j]] > 0))
prev_eff <- sapply(setdiff(names(wide_eff), "case_id"), function(cn) mean(wide_eff[[cn]] == 1, na.rm = TRUE))
keep_eff <- names(prev_eff)[prev_eff >= min_prev]
wide_eff <- wide_eff[, c("case_id", keep_eff), with = FALSE]
message("[INFO] Retained ", length(keep_eff), " gene+type features ≥ ", 100*min_prev, "% prevalence")

# ============================================================
# OUTPUTS
# ============================================================
out_paths <- list(
  baseline = file.path(out_dir, "01_per_gene_baseline/mutation_covariates.by_gene.csv"),
  by_type  = file.path(out_dir, "04_by_gene_and_type/mutation_covariates.by_gene_type.csv")
)
for (p in out_paths) dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)

fwrite(wide_gene, out_paths$baseline)
fwrite(wide_eff, out_paths$by_type)

message("[DONE] Mutation matrices written:")
message(" - Baseline: ", out_paths$baseline)
message(" - By-type : ", out_paths$by_type)

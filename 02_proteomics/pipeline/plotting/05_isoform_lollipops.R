#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------
# CLI + setup
# ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 05_isoform_lollipops.R <cox_long_csv> <cox_sig_csv> <out_dir>")
}

cox_long_path <- args[[1]]
cox_sig_path  <- args[[2]]
summary_dir   <- args[[3]]

dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

say <- function(...) message(sprintf(...))

say("[init] long = %s", cox_long_path)
say("[init] sig  = %s", cox_sig_path)
say("[init] out  = %s", summary_dir)

# Colours (consistent with earlier scripts)
type_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00"
)

# ------------------------------------------------------------------
# Load long + significant Cox tables
# ------------------------------------------------------------------
long_dt <- fread(cox_long_path)
sig_dt  <- fread(cox_sig_path)

# basic sanity
if (!all(c("cancer","data_type","gene_id","gene","feature_id","logHR","direction") %in% names(long_dt))) {
  stop("long_dt is missing required columns (expected at least: cancer, data_type, gene_id, gene, feature_id, logHR, direction)")
}
if (!all(c("cancer","data_type","gene_id","gene","feature_id","logHR","direction") %in% names(sig_dt))) {
  stop("sig_dt is missing required columns (expected at least: cancer, data_type, gene_id, gene, feature_id, logHR, direction)")
}

# ------------------------------------------------------------------
# Isoform-focused summaries
# ------------------------------------------------------------------
iso_summary <- sig_dt[data_type %in% c("iso_log","iso_frac"),
                      .(
                        n_iso_sig      = .N,
                        n_isoforms     = uniqueN(feature_id),
                        n_datatypes    = uniqueN(data_type),
                        has_risk       = any(direction == "risk"),
                        has_protective = any(direction == "protective")
                      ),
                      by = .(cancer, study, gene_id, gene)]
iso_summary[, heterogeneity := has_risk & has_protective]

# Base per-gene isoform summary (reused in C2/C3)
iso_gene_summary_base <- sig_dt[data_type %in% c("iso_log", "iso_frac"),
                                .(
                                  n_isoforms = uniqueN(feature_id),
                                  n_risk     = sum(direction == "risk"),
                                  n_prot     = sum(direction == "protective")
                                ),
                                by = .(cancer, gene_id, gene)]
iso_gene_summary_base[, total := n_risk + n_prot]

# Gene-level Cox results (reused)
if (!"pval" %in% names(long_dt) && "p" %in% names(long_dt)) {
  long_dt[, pval := as.numeric(p)]
}
gene_info_all <- long_dt[
  data_type == "gene",
  .SD[which.min(pval)],
  by = .(cancer, gene_id, gene)
]

# ------------------------------------------------------------------
# DIAGNOSTICS: how many genes per cancer for each plot flavour
# ------------------------------------------------------------------

# C1: heterogeneous genes with >= min_isoforms_mixed
min_isoforms_mixed <- 3

C1_pre <- iso_summary[
  n_isoforms >= min_isoforms_mixed & heterogeneity == TRUE,
  .(n_genes_hetero = uniqueN(gene_id)),
  by = cancer
][order(-n_genes_hetero)]

C1_post <- merge(
  iso_summary[
    n_isoforms >= min_isoforms_mixed & heterogeneity == TRUE,
    .(cancer, gene_id, gene)
  ],
  gene_info_all[, .(cancer, gene_id, gene)],
  by = c("cancer", "gene_id", "gene"),
  all = FALSE
)[, .(n_genes_plottable = uniqueN(gene_id)), by = cancer][order(-n_genes_plottable)]

fwrite(C1_pre,  file.path(summary_dir, "C1_hetero_genes_per_cancer_premerge.csv"))
fwrite(C1_post, file.path(summary_dir, "C1_hetero_genes_per_cancer_postmerge.csv"))
say("[diag] Wrote C1 diagnostics")

# C2: both directions, stricter
min_isoforms <- 5

iso_gene_summary_C2 <- iso_gene_summary_base[
  total > 0 & n_isoforms >= min_isoforms & n_risk > 0 & n_prot > 0
]

C2_pre <- iso_gene_summary_C2[
  , .(n_genes_mixeddir = uniqueN(gene_id)),
  by = cancer
][order(-n_genes_mixeddir)]

C2_post <- merge(
  iso_gene_summary_C2[, .(cancer, gene_id, gene)],
  gene_info_all[, .(cancer, gene_id, gene)],
  by = c("cancer", "gene_id", "gene"),
  all = FALSE
)[, .(n_genes_plottable = uniqueN(gene_id)), by = cancer][order(-n_genes_plottable)]

fwrite(C2_pre,  file.path(summary_dir, "C2_mixed_genes_per_cancer_premerge.csv"))
fwrite(C2_post, file.path(summary_dir, "C2_mixed_genes_per_cancer_postmerge.csv"))
say("[diag] Wrote C2 diagnostics")

# C3: softer criteria, any direction, classified by prop_risk
iso_gene_summary_C3 <- iso_gene_summary_base[
  total > 0 & n_isoforms >= min_isoforms
]
iso_gene_summary_C3[, prop_risk := n_risk / total]
iso_gene_summary_C3[, category := fifelse(
  prop_risk >= 2/3, "Risk-dominated",
  fifelse(prop_risk <= 1/3, "Protective-dominated", "Balanced mixed")
)]

C3_pre <- iso_gene_summary_C3[
  , .(n_genes_eligible = uniqueN(gene_id)),
  by = cancer
][order(-n_genes_eligible)]

C3_pre_cat <- iso_gene_summary_C3[
  , .(n_genes = uniqueN(gene_id)),
  by = .(cancer, category)
][order(cancer, category)]

fwrite(C3_pre,     file.path(summary_dir, "C3_genes_per_cancer_premerge.csv"))
fwrite(C3_pre_cat, file.path(summary_dir, "C3_genes_per_cancer_category_premerge.csv"))
say("[diag] Wrote C3 diagnostics")

say("[done] Isoform lollipop pipeline finished.")

#!/usr/bin/env Rscript

# 01_qc.R
# Input:
#   1) expr_in: matrix CSV (rows = features, cols = case_id)
#   2) clinical_in: CSV with case_id, OS_time, OS_event, + covariates
#   3) thresholds_yaml: YAML with transcriptomics thresholds (optional)
#   4) out_expr_qc: output CSV (QC-passed matrix, features x case_id)
#
# Side outputs (next to out_expr_qc):
#   *_qc_summary.txt, *_pca.png, *_density.png
#
# Example:
#   Rscript 01_qc.R data/iso/BRCA.csv clinical/BRCA.csv 00_config/features.yaml \
#                   01_transcriptomics/01_qc/out/BRCA.iso.qc_passed.csv

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: 01_qc.R <expr_in.csv> <clinical_in.csv> <thresholds.yaml> <out_expr_qc.csv>")
}
expr_in <- args[[1]]
clinical_in <- args[[2]]
yaml_in <- args[[3]]
out_expr <- args[[4]]

# ---------- Load ----------
message("[01_qc] Loading inputs")
X <- fread(expr_in)
X <- as.data.frame(X)
rownames(X) <- X[[1]]
X[[1]] <- NULL

clin <- fread(clinical_in)
clin <- as.data.frame(clin)

cfg <- list(
  transcriptomics = list(
    min_expr = 1.0,                # e.g., TPM > 1
    min_prop_samples = 0.2,        # expressed in >= 20% of samples
    max_zero_prop_sample = 0.9,    # drop sample if >90% zeros
    min_total_expr_quantile = 0.01 # drop bottom 1% library-size samples
  )
)
if (file.exists(yaml_in)) {
  y <- yaml::read_yaml(yaml_in)
  if (!is.null(y$transcriptomics)) cfg$transcriptomics <- modifyList(cfg$transcriptomics, y$transcriptomics)
}

# ---------- Basic checks ----------
stopifnot("case_id" %in% colnames(clin))
if (!all(sapply(X, is.numeric))) {
  X[] <- lapply(X, function(col) as.numeric(as.character(col)))
}

# Intersect samples with clinical
common <- intersect(colnames(X), as.character(clin$case_id))
if (length(common) < 10) stop("Too few overlapping samples between expression and clinical.")
X <- X[, common, drop = FALSE]
clin <- clin[match(common, clin$case_id), , drop = FALSE]

# ---------- Sample QC ----------
message("[01_qc] Sample QC")
total_expr <- colSums(X, na.rm = TRUE)
qcut <- quantile(total_expr, probs = cfg$transcriptomics$min_total_expr_quantile, na.rm = TRUE)
keep_total <- total_expr >= qcut

zero_prop <- colMeans(X == 0 | is.na(X))
keep_zero <- zero_prop <= cfg$transcriptomics$max_zero_prop_sample

keep_samples <- keep_total & keep_zero
X <- X[, keep_samples, drop = FALSE]
clin <- clin[keep_samples, , drop = FALSE]

# ---------- Feature QC ----------
message("[01_qc] Feature QC")
min_expr <- cfg$transcriptomics$min_expr
min_prop <- cfg$transcriptomics$min_prop_samples
keep_feat <- rowMeans(X > min_expr, na.rm = TRUE) >= min_prop
X <- X[keep_feat, , drop = FALSE]

# ---------- Diagnostics ----------
# Density plot (log2(TPM+1) assumption for visualization)
logX <- log2(X + 1)
dens_df <- data.frame(value = as.numeric(as.matrix(logX)))
png(sub("\\.csv$", "_density.png", out_expr), width = 900, height = 600)
ggplot(dens_df, aes(value)) + geom_density() + ggtitle("Density of log2(TPM+1)") + xlab("log2(TPM+1)")
dev.off()

# PCA
centered <- t(scale(t(logX), center = TRUE, scale = TRUE))
pca <- prcomp(t(centered), scale. = FALSE)
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], case_id = colnames(X))
png(sub("\\.csv$", "_pca.png", out_expr), width = 900, height = 700)
ggplot(pca_df, aes(PC1, PC2)) + geom_point() + ggtitle("PCA of samples (log2(TPM+1), centered)")
dev.off()

# ---------- Write outputs ----------
dir.create(dirname(out_expr), showWarnings = FALSE, recursive = TRUE)
fwrite(data.table(feature = rownames(X), X), out_expr)

sink(sub("\\.csv$", "_qc_summary.txt", out_expr))
cat("=== 01_qc summary ===\n")
cat("Input expr:", expr_in, "\n")
cat("Input clinical:", clinical_in, "\n")
cat("Samples kept:", ncol(X), "\n")
cat("Features kept:", nrow(X), "\n")
cat("min_expr:", min_expr, " | min_prop_samples:", min_prop, "\n")
cat("Dropped samples (low total expr or high zeros):", sum(!keep_samples), "\n")
sink()

message("[01_qc] Done. Wrote: ", out_expr)

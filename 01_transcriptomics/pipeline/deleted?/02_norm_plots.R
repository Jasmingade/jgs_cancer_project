#!/usr/bin/env Rscript

# 02_qc_summarize_per_cancer.R
# Build ONE concise PDF per cancer comparing gene / iso_log / iso_frac side-by-side:
# - Violin+box (detected features)
# - Violin+box (total TPM; log10 axis)
# - Violin+box (median log2(TPM+1))
# - PCA per data type (color by stage if available)
#
# Usage:
# Rscript 02_qc_summarize_per_cancer.R <CANCER> [PLOT_ROOT] [NORM_ROOT] [OUT_DIR]
#
# Defaults:
#   PLOT_ROOT = 01_transcriptomics/out/02_plots
#   NORM_ROOT = 01_transcriptomics/out/02_norm_batch
#   OUT_DIR   = 01_transcriptomics/out/02_plots/TCGA_<CANCER>

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(matrixStats)
  library(patchwork)  # arrange plots
})

die <- function(...) { message(sprintf(...)); quit(status=1) }
say <- function(...) message(sprintf(...))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  die("Usage: 02_qc_summarize_per_cancer.R <CANCER> [PLOT_ROOT] [NORM_ROOT] [OUT_DIR]")
}
CANCER    <- args[[1]]
PLOT_ROOT <- if (length(args) >= 2) args[[2]] else "01_transcriptomics/out/02_plots"
NORM_ROOT <- if (length(args) >= 3) args[[3]] else "01_transcriptomics/out/02_norm_batch"
OUT_DIR   <- if (length(args) >= 4) args[[4]] else file.path(PLOT_ROOT, paste0("TCGA_", CANCER))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

DTYPES <- c("gene","iso_log","iso_frac")

qc_path <- function(dtype) file.path(PLOT_ROOT, paste0("TCGA_", CANCER), dtype,
                                     sprintf("TCGA_%s_%s_qc_metrics.csv", CANCER, dtype))
expr_path <- function(dtype) file.path(NORM_ROOT,
                                       sprintf("TCGA_%s_%s.normalized.csv", CANCER, dtype))
mani_path <- function(dtype) file.path(NORM_ROOT,
                                       sprintf("TCGA_%s_%s.sample_manifest.csv", CANCER, dtype))

# ---------- load QC metrics (must exist for the plot set) ----------
qc_list <- list()
avail <- character(0)
for (d in DTYPES) {
  fp <- qc_path(d)
  if (file.exists(fp)) {
    dt <- tryCatch(fread(fp), error=function(e) NULL)
    if (!is.null(dt) && nrow(dt)) {
      dt[, dtype := d]
      qc_list[[d]] <- dt
      avail <- c(avail, d)
    }
  }
}
if (length(avail) == 0) die("No QC metrics found for TCGA_%s under %s", CANCER, PLOT_ROOT)

qc_all <- rbindlist(qc_list, fill=TRUE, use.names=TRUE)

# ---------- small helpers for consistent theming ----------
theme_min <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())

violin_box <- function(df, x, y, title, log10_y = FALSE) {
  p <- ggplot(df, aes_string(x=x, y=y)) +
    geom_violin(trim=TRUE, scale="width") +
    geom_boxplot(width=0.15, outlier.alpha=0.2) +
    labs(title=title, x=NULL, y=NULL) +
    theme_min
  if (log10_y) p <- p + scale_y_log10()
  p
}

# ---------- make distribution panels (one facet per dtype) ----------
p_detect <- violin_box(qc_all, "dtype", "detected_ge1", "Detected features (per sample)")
p_total  <- violin_box(qc_all, "dtype", "total_tpm", "Total TPM (per sample, log10)", log10_y = TRUE)
p_median <- violin_box(qc_all, "dtype", "median_log2", "Median log2(TPM+1) (per sample)")

# ---------- PCA per dtype (recompute from normalized matrices) ----------
pca_plot <- function(dtype) {
  ex <- expr_path(dtype)
  mani <- mani_path(dtype)
  if (!file.exists(ex) || !file.exists(mani)) return(NULL)

  Xdt <- tryCatch(fread(ex), error=function(e) NULL)
  Mdt <- tryCatch(fread(mani), error=function(e) NULL)
  if (is.null(Xdt) || is.null(Mdt) || !nrow(Xdt) || !nrow(Mdt)) return(NULL)

  if (!"feature" %in% names(Xdt)) {
    # assume first column is IDs
    setnames(Xdt, names(Xdt)[1], "feature")
  }
  feats <- Xdt$feature
  Xdt[, feature := NULL]
  for (j in seq_len(ncol(Xdt))) set(Xdt, j=j, value=as.numeric(Xdt[[j]]))
  X <- as.matrix(Xdt); rownames(X) <- feats

  # Align samples
  samp <- intersect(colnames(X), Mdt$case_id)
  if (length(samp) < 10) return(NULL)
  X <- X[, samp, drop=FALSE]
  setkey(Mdt, case_id); Mdt <- Mdt[J(samp)]

  # already log2(TPM+1) → center/scale rows, top variance
  vars <- rowVars(X, na.rm=TRUE); ord <- order(vars, decreasing = TRUE)
  top <- head(ord, min(2000, length(ord)))
  Xp <- t(scale(t(X[top, , drop=FALSE]), center=TRUE, scale=TRUE))
  pc <- tryCatch(prcomp(t(Xp), center=FALSE, scale.=FALSE), error=function(e) NULL)
  if (is.null(pc)) return(NULL)

  pcdt <- data.table(sample = rownames(pc$x), PC1 = pc$x[,1], PC2 = pc$x[,2])
  if ("stage" %in% names(Mdt)) pcdt <- merge(pcdt, Mdt[, .(sample=case_id, stage)], by="sample", all.x=TRUE)

  gg <- ggplot(pcdt, aes(PC1, PC2)) +
    geom_point(aes(color = stage), size=1.8, alpha=0.85, na.rm=TRUE) +
    labs(title = sprintf("PCA — %s", dtype), x="PC1", y="PC2") +
    theme_min + theme(legend.position="none")
  gg
}

p_pca <- lapply(avail, pca_plot)
names(p_pca) <- avail
p_pca <- Filter(Negate(is.null), p_pca)

# If one or two PCAs are missing (e.g., missing normalized matrix), just layout what we have.
pca_row <- NULL
if (length(p_pca) > 0) {
  # keep order gene, iso_log, iso_frac when available
  in_order <- intersect(DTYPES, names(p_pca))
  p_list <- lapply(in_order, function(d) p_pca[[d]])
  # arrange horizontally
  pca_row <- Reduce(`+`, p_list) + plot_layout(ncol = length(p_list))
}

# ---------- assemble final layout ----------
top_row    <- p_detect + p_total + p_median + plot_layout(ncol = 3)
final_plot <- if (!is.null(pca_row)) top_row / pca_row else top_row
final_plot <- final_plot + plot_annotation(
  title = sprintf("TCGA %s — QC Summary (gene vs iso_log vs iso_frac)", CANCER),
  subtitle = "Distributions across samples per data type; PCA colored by stage when available",
  theme = theme(plot.title = element_text(face="bold", size=14))
)

# ---------- write PDF ----------
out_pdf <- file.path(OUT_DIR, sprintf("TCGA_%s_qc_summary_minimal.pdf", CANCER))
ggsave(out_pdf, final_plot, width = 11, height = if (is.null(pca_row)) 4.2 else 7.8, dpi = 150)
say("QC summary written: %s", out_pdf)

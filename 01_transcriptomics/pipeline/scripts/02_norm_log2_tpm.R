#!/usr/bin/env Rscript
# ============================================================
# 02_norm_log2_tpm.R (final version)
# ------------------------------------------------------------
# - Normalizes expression matrices (gene, iso_log, iso_frac)
# - gene / iso_log → log2(TPM + 1)
# - iso_frac → raw clipped [ε, 1−ε] fractions
# - Parallel-friendly, YAML-driven covariates
# - Optional QC plots for iso_frac
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(matrixStats)
  library(preprocessCore)
  library(ggplot2)
})

# -------------------- runtime knobs --------------------
get_threads <- function() {
  env <- Sys.getenv("DT_THREADS", "")
  if (nzchar(env)) return(as.integer(env))
  pmax(1L, parallel::detectCores(logical = TRUE) - 1L)
}
DT_THREADS <- get_threads()
data.table::setDTthreads(DT_THREADS)
Sys.setenv(OMP_NUM_THREADS = 1, MKL_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1)

RUN_PLOTS <- as.logical(Sys.getenv("RUN_PLOTS", "FALSE"))
PLOT_DIR  <- Sys.getenv("PLOT_DIR", "01_transcriptomics/out/02_norm/plots")
datatype_env <- Sys.getenv("DATATYPE", "")

# -------------------- helpers --------------------
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
first12 <- function(x) substr(x, 1, 12)
say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }

# -------------------- args --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
  die("Usage: 02_norm_log2_tpm.R <expr_in.csv> <clinical_in.csv> <covars.yaml> <out_expr.csv> <out_manifest.csv>")

expr_in  <- args[[1]]
clin_in  <- args[[2]]
cov_yaml <- args[[3]]
out_expr <- args[[4]]
out_mani <- args[[5]]

say("=== 02_norm_log2_tpm (threads=%d) ===", DT_THREADS)
say("Expr in:   %s", expr_in)
say("Clinical:  %s", clin_in)
say("Covars YML:%s", cov_yaml)
say("Out expr:  %s", out_expr)
say("Out mani:  %s", out_mani)

# -------------------- read config --------------------
if (!file.exists(cov_yaml)) die("Config YAML not found: %s", cov_yaml)
cfg <- yaml::read_yaml(cov_yaml)
stage_levels   <- cfg$stage_levels %||% c("I","II","III","IV")
baseline_covs  <- cfg$baseline_covariates %||% c("age","sex","stage")
conditional_covs <- cfg$conditional_covariates %||% list()

# -------------------- load expression --------------------
if (!file.exists(expr_in)) die("Expression file not found: %s", expr_in)
Xdt <- fread(expr_in)
if (!("feature" %in% names(Xdt))) setnames(Xdt, names(Xdt)[1], "feature")
feats <- Xdt$feature
Xdt[, feature := NULL]
for (j in seq_len(ncol(Xdt))) set(Xdt, j = j, value = as.numeric(Xdt[[j]]))
X <- as.matrix(Xdt)
rownames(X) <- feats
say("Expr dims (raw): %d features × %d samples", nrow(X), ncol(X))

# -------------------- load clinical --------------------
if (!file.exists(clin_in)) die("Clinical file not found: %s", clin_in)
clin <- fread(clin_in)
need_cols <- c("case_id","OS_time","OS_event")
if (!all(need_cols %in% names(clin)))
  die("Clinical file missing required columns: %s", paste(setdiff(need_cols, names(clin)), collapse=", "))

if ("age" %in% names(clin)) clin[, age := as.numeric(age)]
if ("sex" %in% names(clin)) {
  clin[, sex := tolower(as.character(sex))]
  clin[sex %in% c("f","female"), sex := "female"]
  clin[sex %in% c("m","male"),   sex := "male"]
  clin[, sex := factor(sex, levels=c("female","male"))]
}
if ("Subtype" %in% names(clin)) setnames(clin, "Subtype", "subtype")

# -------------------- align samples --------------------
expr_cases <- colnames(X)
clin_cases <- as.character(clin$case_id)
common <- intersect(expr_cases, clin_cases)
say("Overlap: %d samples (expr %d, clinical %d)", length(common), length(expr_cases), length(clin_cases))

if (length(common) == 0) {
  expr12 <- first12(expr_cases)
  clin12 <- first12(clin_cases)
  common12 <- intersect(expr12, clin12)
  say("Patient-level overlap: %d", length(common12))

  pats <- sort(unique(common12))
  idx_list <- split(seq_along(expr_cases), expr12)
  Xpat <- matrix(NA_real_, nrow = nrow(X), ncol = length(pats),
                 dimnames = list(rownames(X), pats))
  for (i in seq_along(pats)) {
    idx <- idx_list[[pats[i]]]
    if (!is.null(idx)) Xpat[, i] <- rowMedians(X[, idx, drop = FALSE], na.rm = TRUE)
  }
  X <- Xpat
  clin[, case_id12 := first12(case_id)]
  clin <- unique(clin, by="case_id12")
  setkey(clin, case_id12)
  common <- intersect(colnames(X), clin$case_id12)
  if (length(common)==0) die("No overlap after patient collapse.")
  X <- X[, common, drop=FALSE]
  clin <- clin[list(common), on="case_id12"]
  clin[, case_id := case_id12]
}
stopifnot(identical(colnames(X), clin$case_id))
say("[INFO] Final aligned samples: %d", ncol(X))

# -------------------- detect datatype --------------------
datatype <- if (nzchar(datatype_env)) datatype_env else {
  if (grepl("iso[_-]?frac", expr_in, ignore.case=TRUE)) "iso_frac" else
  if (grepl("iso[_-]?log", expr_in, ignore.case=TRUE)) "iso_log" else
  if (grepl("gene", expr_in, ignore.case=TRUE)) "gene" else "unknown"
}
say("Detected datatype: %s", datatype)

# ============================================================
# Normalization
# ============================================================
out_dir <- dirname(out_expr)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

if (datatype %in% c("gene", "iso_log")) {
  say("[NORM] Applying log2(TPM+1) normalization for %s...", datatype)
  Xn <- log2(X + 1)
} else if (datatype == "iso_frac") {
  say("[NORM] Applying raw clipping (ε-adjustment) for iso_frac fractions...")
  eps <- 1e-6
  Xn <- pmin(pmax(X, eps), 1 - eps)
  say("[QC] iso_frac zeros after clip: %.2f%% | ones: %.2f%%", 
      100 * mean(Xn == 0), 100 * mean(Xn == 1))
} else {
  say("[WARN] Unknown data type: %s → skipping normalization.", datatype)
  Xn <- X
}

# ============================================================
# Outputs
# ============================================================
fwrite(data.table(feature = rownames(Xn), Xn), out_expr)
say("[DONE] Normalized expression written: %s", out_expr)

keep_cov <- intersect(names(clin), c("case_id","OS_time","OS_event","age","sex","stage","subtype"))
fwrite(clin[, ..keep_cov], out_mani)
say("[DONE] Sample manifest written: %s", out_mani)

# ============================================================
# Optional QC plot for iso_frac
# ============================================================
if (RUN_PLOTS && datatype == "iso_frac") {
  say("[PLOT] Generating QC plots for iso_frac (raw clipped)...")
  if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive=TRUE, showWarnings=FALSE)

  pdf_out <- file.path(PLOT_DIR, sprintf("QC_iso_frac_raw_%s.pdf", basename(expr_in)))
  pdf(pdf_out, width=10, height=5)
  df_long <- melt(data.table(feature=rownames(Xn), Xn), id.vars="feature", 
                  variable.name="sample", value.name="expr")
  p_hist <- ggplot(df_long, aes(x=expr)) +
    geom_histogram(bins=80, fill="steelblue", color="black") +
    theme_bw() + labs(title="iso_frac distribution (raw clipped)",
                      x="Fraction value", y="Count")
  p_box <- ggplot(df_long, aes(y=expr)) +
    geom_boxplot(fill="gray70", color="black", outlier.shape=NA) +
    theme_bw() + labs(title="Distribution spread (raw clipped)",
                      y="Fraction value", x="")
  print(p_hist)
  print(p_box)
  dev.off()
  say("[DONE] QC PDF saved: %s (%.1f KB)", pdf_out, file.info(pdf_out)$size/1024)
}

say("=== Normalization complete ✓ ===")

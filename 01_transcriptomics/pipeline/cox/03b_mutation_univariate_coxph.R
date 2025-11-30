#!/usr/bin/env Rscript
# ============================================================
# 03b_mutation_univariate_coxph.R
# ------------------------------------------------------------
# Runs univariate Cox models per mutation feature (per cancer,
# per mutation group), using the same structure as baseline
# expression CoxPH (03a_univariate_coxph.R).
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(parallel)
  library(yaml)
})

say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }

# ============================================================
# Arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
    die("Usage: 03b_mutation_univariate_coxph.R <mutation_matrix.csv> <manifest.csv> <covariates.yaml> <out_results.csv> <out_summary.txt>")
}

mut_in  <- args[[1]]   # mutation matrix (feature rows, case columns)
mani_in <- args[[2]]   # manifest with case_id, OS_time, OS_event, covariates
cov_yaml <- args[[3]]  # YAML with baseline + conditional covariates
out_res  <- args[[4]]  # final Cox results table
out_sum  <- args[[5]]  # summary text file

outdir <- dirname(out_res)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

say("=== 03b_mutation_univariate_coxph ===")
say("[INFO] Mutation file: %s", mut_in)
say("[INFO] Manifest: %s", mani_in)

# Extract cancer + mutation group (for metadata)
cancer <- sub("^.*TCGA_([A-Z0-9]+).*", "\\1", basename(mut_in))
mut_group <- sub("^.*ensembl_(.*)\\.csv$", "\\1", basename(mut_in))
if (mut_group == basename(mut_in)) mut_group <- "unknown"

run_id <- sprintf("[03b][%s][%s]", cancer, mut_group)

say("=== 03b_mutation_univariate_coxph ===")
say("%s Cancer=%s | MutationGroup=%s", run_id, cancer, mut_group)

# ============================================================
# Load mutation matrix
# ============================================================
mut <- fread(mut_in)

if (!"feature_id" %in% names(mut)) {
    die("Mutation file must have feature_id as first column.")
}

feature_ids <- mut$feature_id
case_ids <- setdiff(names(mut), "feature_id")

mut_mat <- as.matrix(mut[, ..case_ids])
rownames(mut_mat) <- feature_ids

say("[INFO] Loaded mutation matrix: %d features × %d samples",
    nrow(mut_mat), ncol(mut_mat))

# ============================================================
# Load clinical manifest
# ============================================================
mani <- fread(mani_in)
cfg  <- yaml::read_yaml(cov_yaml)

if (!all(c("case_id","OS_time","OS_event") %in% names(mani)))
    die("Manifest must contain: case_id, OS_time, OS_event")

# ============================================================
# Align cases
# ============================================================
mani <- mani[complete.cases(mani[, .(OS_time, OS_event)])]

common <- intersect(mani$case_id, colnames(mut_mat))
if (length(common) < 20) die("Too few overlapping samples.")

mut_mat <- mut_mat[, common, drop = FALSE]
mani <- mani[match(common, mani$case_id)]
y <- with(mani, Surv(OS_time, OS_event))

say("[INFO] Overlapping cases: %d", length(common))

# ============================================================
# Covariates
# ============================================================
covariates <- cfg$baseline_covariates

if (!is.null(cfg$conditional_covariates)) {
  for (cond in cfg$conditional_covariates) {
    nm <- cond$name
    thr <- cond$coverage_threshold
    if (nm %in% names(mani) && mean(!is.na(mani[[nm]])) >= thr)
      covariates <- c(covariates, nm)
  }
}

covariates <- unique(covariates)
Xcov <- mani[, covariates, with = FALSE]

# remove covariates with no variation
Xcov <- Xcov[, lapply(.SD, function(x)
  if (length(unique(na.omit(x))) > 1) x else NULL
)]

say("[INFO] Covariates used: %s", paste(names(Xcov), collapse=", "))

# ============================================================
# Build Cox model formula
# ============================================================
cov_terms <- paste(names(Xcov), collapse = " + ")

if (nchar(cov_terms) > 0) {
    fml <- as.formula(paste("y ~ mut +", cov_terms))
} else {
    fml <- as.formula("y ~ mut")
}

say("[INFO] Cox formula: %s", deparse(fml))

# ============================================================
# Univariate CoxPH per mutation feature
# ============================================================
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "6"))
say("[INFO] Using %d cores", ncores)

feat_idx <- seq_len(nrow(mut_mat))

res_list <- mclapply(feat_idx, function(i) {
  tryCatch({
    vals <- as.numeric(mut_mat[i, ])
    # Skip features with all NA
    if (all(is.na(vals))) return(data.frame())

    # >>> IMPORTANT FILTER <<<
    # Require at least 5 mutated samples
    if (sum(vals, na.rm = TRUE) < 3)
      return(data.frame())

    df <- data.frame(mut = vals, Xcov, y)

    fit <- suppressWarnings(try(coxph(fml, data=df, ties="efron"), silent=TRUE))
    if (inherits(fit, "try-error") || isFALSE(fit$converged)) return(data.frame())

    s <- summary(fit)
    if (!"mut" %in% rownames(s$coef)) return(data.frame())

    z_mut        <- s$coef["mut","z"]
    wald_mut     <- z_mut^2
    delta_LL_mut <- wald_mut / 2

    data.frame(
      feature      = rownames(mut_mat)[i],
      beta         = s$coef["mut","coef"],
      HR           = exp(s$coef["mut","coef"]),
      z            = z_mut,
      p            = s$coef["mut","Pr(>|z|)"],
      se           = s$coef["mut","se(coef)"],
      delta_LL_mut = delta_LL_mut
    )
  }, error = function(e) data.frame())
}, mc.cores=ncores)

res_list <- Filter(function(x) nrow(x) > 0, res_list)

if (length(res_list) == 0) {
  # No valid fits AFTER filters -> treat as "empty result", not an error
  say("[WARN] No valid CoxPH models fitted for cancer=%s, mut_group=%s (after filters). Writing empty output.",
      cancer, mut_group)

  # empty result table with the expected columns
  res <- data.table(
    feature    = character(),
    beta       = numeric(),
    HR         = numeric(),
    z          = numeric(),
    p          = numeric(),
    se         = numeric(),
    cancer     = character(),
    mut_group  = character(),
    FDR        = numeric()
  )

  # write empty results file so downstream scripts don't break
  data.table::fwrite(res, out_res)

  # write a minimal summary based only on saved output
  sink(out_sum)
  cat("=== Mutation Cox Summary ===\n")
  cat("Cancer:", cancer, "\n")
  cat("Mutation group:", mut_group, "\n")
  cat("Total features:", nrow(res), "\n")
  cat("Valid p-values:", 0, "\n")
  cat("Invalid p-values (removed):", 0, "\n")
  cat("Significant (p<0.05):", 0, "\n")
  cat("Significant (FDR<0.05):", 0, "\n")
  sink()

  say("%s [WARN] Empty result set written → %s", run_id, out_res)
  say("%s Summary saved to %s", run_id, out_sum)

  quit(status = 0)  # SUCCESS → bash set -e is happy
}

# Combine & annotate
res <- rbindlist(res_list, fill = TRUE)

res[, p := as.numeric(p)]
res <- res[is.finite(p)]
res[, FDR := p.adjust(p, "BH")]
res[, cancer := cancer]
res[, mut_group := mut_group]

total_features   <- nrow(res)
valid_p_total    <- sum(!is.na(res$p))
invalid_p_total  <- total_features - valid_p_total

# ---- FULL results ----
res_full <- copy(res)
out_res_full <- sub("\\.cox_results\\.csv$", ".cox_results_full.csv", out_res)
fwrite(res_full, out_res_full)

# ---- SIGNIFICANT subset ----
res_sig <- res_full[FDR < 0.05 & is.finite(HR) & HR > 0]
fwrite(res_sig, out_res)

say("[DONE] Wrote %d full mutation models → %s", nrow(res_full), out_res_full)
say("[DONE] Wrote %d significant mutation models → %s", nrow(res_sig), out_res)

# ---- SUMMARY ----
sink(out_sum)
cat("=== Mutation Cox Summary ===\n")
cat("Cancer:", cancer, "\n")
cat("Mutation group:", mut_group, "\n")
cat("Total features:", total_features, "\n")
cat("Valid p-values:", valid_p_total, "\n")
cat("Invalid p-values (removed):", invalid_p_total, "\n")

cat("Significant (p<0.05):", sum(res_sig$p < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR<0.05):", sum(res_sig$FDR < 0.05, na.rm = TRUE), "\n")
sink()

say("%s Summary saved to %s", run_id, out_sum)

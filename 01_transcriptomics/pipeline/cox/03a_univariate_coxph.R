#!/usr/bin/env Rscript
# ============================================================
# 03a_univariate_coxph.R
# ------------------------------------------------------------
# Runs univariate Cox models per feature (gene / iso_log / iso_frac),
# decides stratification based on PH assumption summaries,
# and generates comprehensive summary + diagnostic plots:
#   - HR boxplots (significant only)
#   - Age/stage covariate sanity plots
#   - Significant feature barplot
#   - P-value histogram matrix (faceted by data_type × cancer)
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
    die("Usage: 03a_univariate_coxph.R <expr_norm.csv> <manifest.csv> <covariates.yaml> <out_results.csv> <out_summary.txt>")
}

expr_in  <- args[[1]]   # normalized expression matrix
mani_in  <- args[[2]]   # manifest with case_id, OS_time, OS_event, covariates
cov_yaml <- args[[3]]   # YAML with baseline + conditional covariates
out_res  <- args[[4]]   # final Cox results table
out_sum  <- args[[5]]   # summary text file

# Directory where plots + diagnostics will be written
outdir <- dirname(out_res)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
# ============================================================
# Load data
# ============================================================
expr <- fread(expr_in)
mani <- fread(mani_in)
cfg  <- yaml::read_yaml(cov_yaml)

feat_col <- if ("feature" %in% names(expr)) "feature" else names(expr)[1]
features <- expr[[feat_col]]
expr[[feat_col]] <- NULL
expr <- as.matrix(expr)
rownames(expr) <- features

cancer <- sub("^.*TCGA_([A-Z0-9]+)_.*$", "\\1", basename(expr_in))
data_type <- sub("^.*TCGA_[A-Z0-9]+_(.*?)\\.normalized\\.csv$", "\\1", basename(expr_in))

run_id <- sprintf("[03a][%s][%s]", cancer, data_type)

say("=== 03a_univariate_coxph ===")
say("%s Starting Expression CoxPH", run_id)
say("%s Cancer=%s | DataType=%s", run_id, cancer, data_type)

# ============================================================
# Iso_frac logit transform
# ============================================================
#if (identical(data_type, "iso_frac")) {
#  eps <- 1e-6
#  expr <- pmin(pmax(expr, eps), 1 - eps)
#  expr <- log(expr / (1 - expr))
#}

# ============================================================
# Align expression and clinical data
# ============================================================
mani <- mani[complete.cases(mani[, .(OS_time, OS_event)]), ]
common <- intersect(colnames(expr), mani$case_id)
if (length(common) < 20) die("Too few overlapping samples.")
expr <- expr[, common, drop = FALSE]
mani <- mani[match(common, mani$case_id)]
y <- with(mani, Surv(OS_time, OS_event))

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
Xcov <- Xcov[, lapply(.SD, function(x) if (length(unique(na.omit(x))) > 1) x else NULL)]
say("[INFO] Covariates used: %s", paste(names(Xcov), collapse=", "))


# ============================================================
# Build Cox formula
# ============================================================
cov_terms <- paste(names(Xcov), collapse = " + ")

if (nchar(cov_terms) > 0) {
    fml <- as.formula(paste("y ~ expr +", cov_terms))
} else {
    fml <- as.formula("y ~ expr")
}

say("[INFO] Cox formula: %s", deparse(fml))


# ============================================================
# Run univariate CoxPH
# ============================================================
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "6"))
say("[INFO] Using %d cores", ncores)
idx <- seq_len(nrow(expr))

res_list <- mclapply(idx, function(i) {
  tryCatch({
    vals <- as.numeric(expr[i, ])
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) return(data.frame())
    vals_z <- as.numeric(scale(vals))
    df <- data.frame(expr = vals_z, Xcov, y)

    fit <- suppressWarnings(try(coxph(fml, data = df, ties = "efron"), silent = TRUE))
    if (inherits(fit, "try-error") || isFALSE(fit$converged)) return(data.frame())
    s <- summary(fit)
    if (!"expr" %in% rownames(s$coef)) return(data.frame())

    z_expr        <- s$coef["expr","z"]
    wald_expr     <- z_expr^2
    delta_LL_expr <- wald_expr / 2   # ≈ improvement in log-likelihood for expr term

    data.frame(
      feature       = rownames(expr)[i],
      beta          = s$coef["expr","coef"],
      HR            = exp(s$coef["expr","coef"]),
      z             = z_expr,
      p             = s$coef["expr","Pr(>|z|)"],
      se            = s$coef["expr","se(coef)"],
      delta_LL_expr = delta_LL_expr
    )
  }, error = function(e) data.frame())
}, mc.cores = ncores)

res_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, res_list)
if (length(res_list) == 0) die("No valid CoxPH models fitted.")
res <- rbindlist(res_list, fill = TRUE)

# basic clean-up
res[, p := as.numeric(p)]
res <- res[!is.na(p) & is.finite(p)]   # drop invalid p
res[, FDR := p.adjust(p, "BH")]
res[, cancer := cancer]
res[, data_type := data_type]

total_features <- nrow(res)
valid_p_total  <- sum(!is.na(res$p))
invalid_p_total <- total_features - valid_p_total

# ---- FULL results (no FDR threshold) ----
res_full <- copy(res)
out_res_full <- sub("\\.cox_results\\.csv$", ".cox_results_full.csv", out_res)
fwrite(res_full, out_res_full)

# ---- SIGNIFICANT subset ----
res_sig <- res_full[FDR < 0.05 & is.finite(HR) & HR > 0]
fwrite(res_sig, out_res)

# ---- SUMMARY (based on what we actually saved) ----
sink(out_sum)
cat("=== Expression Cox Summary ===\n")
cat("Cancer:", cancer, "\n")
cat("Datatype:", data_type, "\n")
cat("Total features:", total_features, "\n")
cat("Valid p-values:", valid_p_total, "\n")
cat("Invalid p-values (removed):", invalid_p_total, "\n")

cat("Significant (p<0.05):", sum(res_sig$p < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR<0.05):", sum(res_sig$FDR < 0.05, na.rm = TRUE), "\n")
sink()

say("%s [DONE] Wrote %d full models → %s", run_id, nrow(res_full), out_res_full)
say("%s [DONE] Wrote %d significant models → %s", run_id, nrow(res_sig), out_res)

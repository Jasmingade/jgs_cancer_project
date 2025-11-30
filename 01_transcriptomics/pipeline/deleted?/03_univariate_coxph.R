#!/usr/bin/env Rscript
# ============================================================
# 03_univariate_coxph.R
# ------------------------------------------------------------
# Runs univariate CoxPH per feature (gene / isoform / iso_frac)
# Supports three modes:
#   ① baseline    -> Survival ~ Expression + Covariates
#   ② combined    -> Survival ~ Expression + Mutation + Covariates
#   ③ interaction -> Survival ~ Expression * Mutation + Covariates
# ------------------------------------------------------------
# Compatible with outputs from 02_norm and 02_mutation/
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(yaml)
  library(parallel)
})

say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a

# ============================================================
# Args
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6)
  die("Usage: 03_univariate_coxph.R <expr_norm.csv> <manifest.csv> <covariates_yaml> <out_results.csv> <out_summary.txt> <mode>")

expr_in <- args[[1]]
mani_in <- args[[2]]
cov_yaml <- args[[3]]
out_res  <- args[[4]]
out_sum  <- args[[5]]
mode     <- tolower(args[[6]])

if (!mode %in% c("baseline", "combined", "interaction"))
  die("Invalid mode. Must be one of: baseline, combined, interaction")

say("=== 03_univariate_coxph ===")
say("Mode: %s", mode)
say("Expression matrix: %s", expr_in)
say("Manifest: %s", mani_in)

# ============================================================
# Load input
# ============================================================
if (!file.exists(expr_in)) die("Expression file missing: %s", expr_in)
if (!file.exists(mani_in)) die("Manifest file missing: %s", mani_in)
if (!file.exists(cov_yaml)) die("YAML file missing: %s", cov_yaml)

expr <- fread(expr_in)
mani <- fread(mani_in)
cfg  <- yaml::read_yaml(cov_yaml)

# Identify feature column (if present)
feat_col <- if ("feature" %in% names(expr)) "feature" else names(expr)[1]
features <- expr[[feat_col]]
expr[[feat_col]] <- NULL
expr <- as.matrix(expr)
rownames(expr) <- features

# Align samples
mani <- mani[complete.cases(mani[, .(OS_time, OS_event)]), ]
common <- intersect(colnames(expr), mani$case_id)
expr <- expr[, common, drop=FALSE]
mani <- mani[match(common, mani$case_id)]
say("Aligned %d samples with survival info", length(common))

y <- with(mani, Surv(OS_time, OS_event))

# ============================================================
# Load mutation data (if needed)
# ============================================================
mut_mat <- NULL
if (mode %in% c("combined", "interaction")) {
  cancer_id <- sub(".*TCGA_([A-Z]+)_.+", "\\1", expr_in)
  mut_path <- sprintf("01_transcriptomics/out/02_mutation/02_combined_expression_mutation/split_by_cancer_type/TCGA_%s_mutation.normalized.csv", cancer_id)
  if (file.exists(mut_path)) {
    say("Loading mutation matrix: %s", mut_path)
    mut_mat <- fread(mut_path)
    mut_mat <- mut_mat[case_id %in% mani$case_id]
    setkey(mut_mat, case_id)
  } else {
    say("[WARN] No mutation matrix found for %s — reverting to baseline mode", cancer_id)
    mode <- "baseline"
  }
}

# ============================================================
# Covariates
# ============================================================
covariates <- cfg$baseline_covariates %||% c("age", "sex", "stage")
if (!is.null(cfg$conditional_covariates)) {
  cond_cov <- sapply(cfg$conditional_covariates, `[[`, "name")
  covariates <- unique(c(covariates, cond_cov))
}
Xcov <- data.frame(row.names = mani$case_id)
for (cov in covariates) if (cov %in% names(mani)) Xcov[[cov]] <- mani[[cov]]
say("Covariates included: %s", paste(names(Xcov), collapse = ", "))

# ============================================================
# Model builder
# ============================================================
make_formula <- function() {
  switch(mode,
         "baseline"    = "y ~ expr + %s",
         "combined"    = "y ~ expr + mutation + %s",
         "interaction" = "y ~ expr * mutation + %s")
}

# ============================================================
# Parallel per-feature Cox fits
# ============================================================
idx <- seq_len(nrow(expr))
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("NCORES", "2")))
say("Using %d cores for %d features", ncores, nrow(expr))

res_list <- mclapply(idx, function(i) {
  vals <- as.numeric(expr[i, ])
  if (all(is.na(vals)) || sd(vals, na.rm=TRUE) == 0) return(NULL)

  df <- data.frame(expr = vals, Xcov)
  if (!is.null(mut_mat) && mode %in% c("combined", "interaction")) {
    gene <- features[i]
    if (gene %in% names(mut_mat)) {
      df$mutation <- mut_mat[match(rownames(Xcov), mut_mat$case_id), get(gene)]
    } else {
      df$mutation <- 0
    }
  }

  rhs <- paste(names(Xcov), collapse = " + ")
  f_str <- sprintf(make_formula(), rhs)
  fit <- try(coxph(as.formula(f_str), data=df, ties="efron"), silent=TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  s <- try(summary(fit), silent=TRUE)
  if (inherits(s, "try-error")) return(NULL)
  co <- s$coef
  if (!"expr" %in% rownames(co)) return(NULL)

  ci <- try(suppressWarnings(confint(fit, 'expr')), silent = TRUE)
  lo <- hi <- NA_real_
  if (!inherits(ci, "try-error") && all(is.finite(ci))) {
    lo <- exp(ci[1]); hi <- exp(ci[2])
  }

  data.frame(
    feature = features[i],
    n       = nrow(model.frame(fit)),
    beta    = co["expr","coef"],
    HR      = co["expr","exp(coef)"],
    z       = co["expr","z"],
    p       = co["expr","Pr(>|z|)"],
    se      = co["expr","se(coef)"],
    HR_lo   = lo,
    HR_hi   = hi,
    cindex  = s$concordance[1]
  )
}, mc.cores = ncores)

# ============================================================
# Combine and post-process
# ============================================================
res_dt <- rbindlist(res_list, fill = TRUE)
if (is.null(res_dt) || nrow(res_dt) == 0) {
  say("[WARN] No valid feature fits — writing empty output.")
  res_dt <- data.table(feature=character(), HR=numeric(), beta=numeric(),
                       z=numeric(), p=numeric(), se=numeric(),
                       HR_lo=numeric(), HR_hi=numeric(), cindex=numeric())
}

res_dt[, `:=`(logHR = log(HR),
              FDR = p.adjust(p, "BH"))]

dir.create(dirname(out_res), recursive = TRUE, showWarnings = FALSE)
fwrite(res_dt, out_res)

# ============================================================
# Summary output
# ============================================================
sink(out_sum)
cat(sprintf("=== CoxPH summary (%s) ===\n", toupper(mode)))
cat("Expression:", expr_in, "\n")
cat("Samples:", ncol(expr), "\n")
cat("Features tested:", nrow(expr), "\n")
cat("Valid models:", nrow(res_dt), "\n")
cat("Median HR:", round(median(res_dt$HR, na.rm=TRUE), 3), "\n")
cat("Significant (p<0.05):", sum(res_dt$p < 0.05, na.rm=TRUE), "\n")
cat("Significant (FDR<0.05):", sum(res_dt$FDR < 0.05, na.rm=TRUE), "\n")
sink()

say("[DONE] Mode=%s → %s", mode, out_res)

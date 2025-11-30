#!/usr/bin/env Rscript
# ============================================================
# 03b_mutation_univariate_coxph.R
# ------------------------------------------------------------
# Runs univariate CoxPH per mutation feature (0/1 per gene)
# adjusted for the same clinical covariates used in omics models.
# Output schema matches 03_univariate_coxph.R, so plotting scripts
# work unchanged for cross‐omics comparison.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(yaml)
  library(parallel)
})

`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }

# ============================================================
# Arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
  die("Usage: 03b_mutation_univariate_coxph.R <mutation_covariates.csv> <manifest.csv> <covariates_yaml> <out_results.csv> <out_summary.txt>")

mut_in <- args[[1]]   # mutation matrix (case_id + genes)
mani_in <- args[[2]]  # manifest (OS_time, OS_event, age, sex, stage, ...)
cov_yaml <- args[[3]] # covariate config
out_res  <- args[[4]] # results CSV
out_sum  <- args[[5]] # summary text

say("=== 03b_mutation_univariate_coxph ===")
say("Mutation matrix : %s", mut_in)
say("Manifest        : %s", mani_in)
say("Covariates YAML : %s", cov_yaml)
say("Output results  : %s", out_res)
say("Output summary  : %s", out_sum)

# ============================================================
# Load data
# ============================================================
if (!file.exists(mut_in))  die("Mutation covariates file missing: %s", mut_in)
if (!file.exists(mani_in)) die("Manifest file missing: %s", mani_in)
if (!file.exists(cov_yaml)) die("YAML file missing: %s", cov_yaml)

mut <- fread(mut_in)
mani <- fread(mani_in)
cfg  <- yaml::read_yaml(cov_yaml)

# required columns
if (!all(c("case_id","OS_time","OS_event", "age", "stage") %in% names(mani)))
  die("Manifest must include case_id, OS_time, OS_event, age, and stage.")
if (!("case_id" %in% names(mut)))
  die("Mutation matrix must include case_id.")

# align samples
mani <- mani[complete.cases(mani[, .(OS_time, OS_event)]), ]
common <- intersect(mani$case_id, mut$case_id)
if (length(common) < 20)
  die("Too few overlapping cases with survival + mutation (n=%d).", length(common))

setkey(mani, case_id)
setkey(mut, case_id)
dat <- mut[mani, on = "case_id"]   # keeps manifest order
dat <- dat[case_id %in% common]

say("Aligned %d cases between mutation and clinical data", length(common))

# ============================================================
# Survival + covariates
# ============================================================
y <- with(dat, Surv(OS_time, OS_event))

coverage_prop <- function(v) {
  x <- tolower(trimws(as.character(v)))
  miss <- is.na(x) | x == "" | x == "-" | x == "na" | x == "n/a" |
          x == "unknown" | x == "unk"
  mean(!miss)
}

covariates <- cfg$baseline_covariates %||% c("age","sex","stage")

# conditional covariates (e.g. subtype)
if (!is.null(cfg$conditional_covariates)) {
  for (cond in cfg$conditional_covariates) {
    cov_name <- cond$name
    ref_var  <- cond$include_if %||% cov_name
    thresh   <- as.numeric(cond$coverage_threshold %||% 1.0)
    ref_col  <- if (paste0(ref_var, "_raw") %in% names(dat)) paste0(ref_var, "_raw") else ref_var
    if (ref_col %in% names(dat)) {
      cov_ok <- coverage_prop(dat[[ref_col]])
      if (cov_ok >= thresh && cov_name %in% names(dat)) {
        covariates <- unique(c(covariates, cov_name))
        say("Including conditional covariate '%s' (%.1f%% coverage)", cov_name, 100*cov_ok)
      }
    }
  }
}

say("Clinical covariates used: %s", paste(covariates, collapse=", "))

# ============================================================
# Select mutation features
# ============================================================
mut_cols <- setdiff(names(mut), "case_id")

# keep features with reasonable prevalence or variance
is_binary <- function(v) all(na.omit(unique(v)) %in% c(0,1,TRUE,FALSE))
keep <- c()
for (cn in mut_cols) {
  v <- dat[[cn]]
  if (all(is.na(v))) next
  if (is_binary(v)) {
    p <- mean(v %in% c(1, TRUE), na.rm = TRUE)
    if (p >= 0.01 && p <= 0.99) keep <- c(keep, cn)
  } else {
    if (sd(v, na.rm = TRUE) > 0) keep <- c(keep, cn)
  }
}
mut_cols <- unique(keep)
if (!length(mut_cols))
  die("No usable mutation features after QC.")
say("Using %d mutation features across %d cases", length(mut_cols), nrow(dat))

# ============================================================
# Model setup
# ============================================================
cov_part <- paste(covariates[covariates %in% names(dat)], collapse = " + ")
base_rhs <- if (nzchar(cov_part)) paste("+", cov_part) else ""
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("NCORES", "2")))
say("Using %d cores", ncores)

# ============================================================
# Per-feature Cox regression
# ============================================================
idx <- seq_along(mut_cols)
res_list <- mclapply(idx, function(i) {
  cn <- mut_cols[i]
  vals <- dat[[cn]]
  if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) return(NULL)

  df <- data.frame(mut = vals)
  for (cv in covariates) if (cv %in% names(dat)) df[[cv]] <- dat[[cv]]

  fit <- try(coxph(as.formula(paste("y ~ mut", base_rhs)),
                   data = df, ties = "efron", na.action = na.omit),
             silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)

  s <- summary(fit)
  ci <- try(suppressWarnings(confint(fit, 'mut')), silent = TRUE)
  lo <- hi <- NA_real_
  if (!inherits(ci, "try-error") && all(is.finite(ci))) {
    lo <- exp(ci[1]); hi <- exp(ci[2])
  }

  data.frame(
    feature = cn,
    n       = nrow(model.frame(fit)),
    beta    = s$coef["mut","coef"],
    HR      = s$coef["mut","exp(coef)"],
    z       = s$coef["mut","z"],
    p       = s$coef["mut","Pr(>|z|)"],
    cindex  = s$concordance[1],
    se      = s$coef["mut","se(coef)"],
    HR_lo   = lo,
    HR_hi   = hi
  )
}, mc.cores = ncores)

# ============================================================
# Combine results and adjust p-values
# ============================================================
res <- rbindlist(res_list, fill = TRUE)
if (is.null(res) || !nrow(res)) {
  say("No valid Cox fits; writing empty results")
  res <- data.table(feature=character(), n=integer(), beta=numeric(),
                    HR=numeric(), z=numeric(), p=numeric(),
                    cindex=numeric(), se=numeric(),
                    HR_lo=numeric(), HR_hi=numeric())
}

res <- res[!is.na(p)]
if (!("HR" %in% names(res))) res[, HR := exp(beta)]
if (!("HR_lo" %in% names(res)) || !("HR_hi" %in% names(res))) {
  res[, `:=`(HR_lo = exp(log(HR) - 1.96*se),
             HR_hi = exp(log(HR) + 1.96*se))]
}
res[, `:=`(logHR = log(HR), FDR = p.adjust(p, "BH"))]

# ============================================================
# Write output
# ============================================================
dir.create(dirname(out_res), recursive = TRUE, showWarnings = FALSE)
fwrite(res, out_res)

sink(out_sum)
cat("=== 03b_mutation_univariate_coxph summary ===\n")
cat("Mutation features tested:", length(mut_cols), "\n")
cat("Features with valid fit:", nrow(res), "\n")
cat("Median HR:", round(median(res$HR, na.rm = TRUE), 3), "\n")
cat("Significant (p<0.05):", sum(res$p < 0.05, na.rm = TRUE), "\n")
cat("Significant (FDR<0.05):", sum(res$FDR < 0.05, na.rm = TRUE), "\n")
sink()

say("[DONE] Mutation Cox results written to %s", out_res)
say("[DONE] Summary written to %s", out_sum)

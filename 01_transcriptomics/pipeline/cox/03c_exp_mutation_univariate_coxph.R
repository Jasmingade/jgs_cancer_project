#!/usr/bin/env Rscript
# ============================================================
# 03c_exp_mutation_univariate_coxph.R
# ------------------------------------------------------------
# Combined Expression + Mutation Cox model:
#
#   Surv(OS_time, OS_event) ~ expr_feature + mut_group + covariates
#
# mut_group = sample-level binary flag (ANY mutated gene in group)
# Works for gene, iso_log, iso_frac.
#
# Called from bash with:
#   Rscript 03c_exp_mutation_univariate_coxph.R \
#       <expr_norm.csv> <mut.csv> <manifest.csv> \
#       <covariates.yaml> <out_results.csv> <out_summary.txt>
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
# Args
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  die(paste(
    "Usage:",
    "03c_exp_mutation_univariate_coxph.R",
    "<expr_norm.csv> <mut.csv> <manifest.csv>",
    "<covariates.yaml> <out_results.csv> <out_summary.txt>"
  ))
}

expr_in  <- args[[1]]
mut_in   <- args[[2]]
mani_in  <- args[[3]]
cov_yaml <- args[[4]]
out_res  <- args[[5]]
out_sum  <- args[[6]]

dir.create(dirname(out_res), recursive = TRUE, showWarnings = FALSE)

say("=== 03c_exp_mutation_univariate_coxph ===")
say("[INFO] expr_in = %s", expr_in)
say("[INFO] mut_in  = %s", mut_in)
say("[INFO] mani_in = %s", mani_in)

# ============================================================
# LOAD EXPRESSION
# ============================================================
expr <- fread(expr_in)
feat_col <- if ("feature" %in% names(expr)) "feature" else names(expr)[1]
features <- expr[[feat_col]]
expr[[feat_col]] <- NULL
expr_mat <- as.matrix(expr)
rownames(expr_mat) <- features

# ============================================================
# LOAD MUTATION → collapse to sample-level indicator
# ============================================================
mut_dt <- fread(mut_in)

if (!"feature_id" %in% names(mut_dt)) {
  die("Mutation file must contain 'feature_id' as first column.")
}

sample_cols <- setdiff(names(mut_dt), "feature_id")

# convert to numeric 0/1
mut_mat <- as.matrix(mut_dt[, ..sample_cols])
mut_mat[is.na(mut_mat)] <- 0
mut_mat <- (mut_mat > 0) * 1

# MUTATION GROUP = sample is mutated in ANY gene of this group
mut_group <- apply(mut_mat, 2, function(x) as.integer(any(x == 1)))

mut_group <- data.table(
  case_id = sample_cols,
  mut     = as.integer(mut_group)
)

say("[INFO] Mut-group: %d mutated / %d samples",
    sum(mut_group$mut), nrow(mut_group))

# ============================================================
# LOAD MANIFEST
# ============================================================
mani <- fread(mani_in)
stopifnot(all(c("case_id", "OS_time", "OS_event") %in% names(mani)))

# ============================================================
# ALIGN CASE IDS
# ============================================================
common <- Reduce(
  intersect,
  list(colnames(expr_mat), mut_group$case_id, mani$case_id)
)

if (length(common) < 20) {
  die(sprintf("Too few overlapping samples: %d", length(common)))
}

common <- sort(common)

expr_mat <- expr_mat[, common, drop = FALSE]
mut_vec  <- mut_group[match(common, mut_group$case_id), mut]
mani     <- mani[match(common, mani$case_id)]

y <- with(mani, Surv(OS_time, OS_event))

events_mut  <- sum(mani$OS_event[mut_vec == 1] == 1, na.rm = TRUE)
events_wild <- sum(mani$OS_event[mut_vec == 0] == 1, na.rm = TRUE)

if (events_mut < 3 || events_wild < 3) {
  die(sprintf("Too few events in one of the mutation groups (mut=%d, wt=%d)",
              events_mut, events_wild))
}
say("[INFO] Aligned samples: %d", length(common))

# ============================================================
# COVARIATES  (from YAML, with coverage check)
# ============================================================
cfg <- yaml::read_yaml(cov_yaml)
covariates <- cfg$baseline_covariates

if (!is.null(cfg$conditional_covariates)) {
  for (cond in cfg$conditional_covariates) {
    nm  <- cond$name
    thr <- cond$coverage_threshold
    if (nm %in% names(mani) && mean(!is.na(mani[[nm]])) >= thr) {
      covariates <- c(covariates, nm)
    }
  }
}

covariates <- unique(covariates)

# subset columns that actually exist
covariates <- covariates[covariates %in% names(mani)]
Xcov <- mani[, covariates, with = FALSE]

# drop covariates with no variability
Xcov <- Xcov[, lapply(.SD, function(x) {
  if (length(unique(na.omit(x))) > 1) x else NULL
})]

cov_used <- names(Xcov)
say("[INFO] Covariates used after filtering: %s",
    if (length(cov_used)) paste(cov_used, collapse = ", ") else "<none>")

# ---- Build formula string safely ----
if (length(cov_used) > 0) {
  cov_terms  <- paste(cov_used, collapse = " + ")
  fml_string <- paste("y ~ expr + mut +", cov_terms)
} else {
  fml_string <- "y ~ expr + mut"
}

# remove any trailing '+'
fml_string <- gsub("\\+\\s*$", "", fml_string)
fml <- as.formula(fml_string)

say("[INFO] Cox formula: %s", fml_string)

# ============================================================
# RUN COXPH PER EXPRESSION FEATURE
# ============================================================
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "6"))
if (is.na(ncores) || ncores < 1) ncores <- 2

say("[INFO] Testing %d expression features", nrow(expr_mat))
say("[INFO] Using %d cores", ncores)

res_list <- mclapply(rownames(expr_mat), function(feat) {
  tryCatch({
    x <- as.numeric(expr_mat[feat, ])

    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(data.frame())

    df <- data.frame(
      expr = as.numeric(scale(x)),
      mut  = mut_vec,
      Xcov,
      y
    )

    fit <- suppressWarnings(try(coxph(fml, data = df, ties = "efron"),
                                silent = TRUE))
    if (inherits(fit, "try-error") || isFALSE(fit$converged)) return(data.frame())

    s <- summary(fit)

    z_expr        <- s$coef["expr","z"]
    wald_expr     <- z_expr^2
    delta_LL_expr <- wald_expr / 2

    z_mut         <- s$coef["mut","z"]
    wald_mut      <- z_mut^2
    delta_LL_mut  <- wald_mut / 2

    data.frame(
      feature       = feat,
      beta_expr     = s$coef["expr","coef"],
      HR_expr       = exp(s$coef["expr","coef"]),
      p_expr        = s$coef["expr","Pr(>|z|)"],

      beta_mut      = s$coef["mut","coef"],
      HR_mut        = exp(s$coef["mut","coef"]),
      p_mut         = s$coef["mut","Pr(>|z|)"],

      delta_LL_expr = delta_LL_expr,
      delta_LL_mut  = delta_LL_mut
    )
  }, error = function(e) data.frame())
}, mc.cores = ncores)

res_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, res_list)

# ----- metadata (for logging & summary) -----
cancer    <- sub("^.*TCGA_([A-Z0-9]+)_.*$", "\\1", basename(expr_in))
data_type <- sub("^.*TCGA_[A-Z0-9]+_(.*?)\\.normalized\\.csv$", "\\1", basename(expr_in))
run_id    <- sprintf("[03c][%s][%s]", cancer, data_type)

say("=== 03c_exp_mutation_univariate_coxph ===")
say("%s expr_in=%s", run_id, expr_in)
say("%s mut_in=%s", run_id, mut_in)
say("%s mani_in=%s", run_id, mani_in)

# ----- branch: no valid models after fitting -----
if (length(res_list) == 0) {
  say("%s [WARN] No valid Cox models fitted (Expression+Mutation). Writing empty output.",
      run_id)

  res <- data.table(
    feature    = character(),
    beta_expr  = numeric(),
    HR_expr    = numeric(),
    p_expr     = numeric(),
    beta_mut   = numeric(),
    HR_mut     = numeric(),
    p_mut      = numeric(),
    FDR_expr   = numeric(),
    FDR_mut    = numeric(),
    cancer     = character(),
    data_type  = character()
  )

  fwrite(res, out_res)

  sink(out_sum)
  cat("=== Expression + Mutation (mut-group) Cox Summary ===\n")
  cat("Cancer:", cancer, "\n")
  cat("Datatype:", data_type, "\n")
  cat("Features tested:", nrow(expr_mat), "\n")
  cat("Valid models:", 0, "\n\n")
  cat("Significant expr term (FDR<0.05):", 0, "\n")
  cat("Significant mut term  (FDR<0.05):", 0, "\n")
  sink()

  say("%s [WARN] Empty Expression+Mutation result written → %s", run_id, out_res)
  say("%s Summary saved to %s", run_id, out_sum)

  quit(status = 0)  # SUCCESS
}

# ----- normal case: we have valid fits -----
res <- rbindlist(res_list)

# FDRs
res[, FDR_expr := p.adjust(p_expr, "BH")]
res[, FDR_mut  := p.adjust(p_mut,  "BH")]
res[, cancer    := cancer]
res[, data_type := data_type]

# metadata
res[, cancer    := cancer]
res[, data_type := data_type]

# ---- FULL results (no FDR or HR filtering) ----
res_full <- copy(res)
out_res_full <- sub("\\.cox_results\\.csv$", ".cox_results_full.csv", out_res)
fwrite(res_full, out_res_full)

# ---- SIGNIFICANT subset ----
res_sig <- res_full[
  (FDR_expr < 0.05 & is.finite(HR_expr) & HR_expr > 0) |
  (FDR_mut  < 0.05 & is.finite(HR_mut)  & HR_mut  > 0)
]
fwrite(res_sig, out_res)

# ---- SUMMARY ----
sink(out_sum)
cat("=== Expression + Mutation (mut-group) Cox Summary ===\n")
cat("Cancer:", cancer, "\n")
cat("Datatype:", data_type, "\n")
cat("Features tested:", nrow(expr_mat), "\n")
cat("Valid models:", nrow(res_full), "\n\n")
cat("Significant expr term (FDR<0.05):", sum(res_sig$FDR_expr < 0.05, na.rm = TRUE), "\n")
cat("Significant mut term  (FDR<0.05):", sum(res_sig$FDR_mut  < 0.05, na.rm = TRUE), "\n")
sink()

say("%s [DONE] Saved FULL Expression+Mutation CoxPH results → %s", run_id, out_res_full)
say("%s [DONE] Saved SIGNIFICANT Expression+Mutation CoxPH results → %s", run_id, out_res)
say("%s Summary saved to %s", run_id, out_sum)

#!/usr/bin/env Rscript
# ============================================================
# 03d_expr_mut_interaction_coxph.R   (generic, all expr modalities)
# ------------------------------------------------------------
# Interaction model per expression feature:
#
#   Surv(OS_time, OS_event) ~ expr * mutation + covariates
#
# expr : feature × sample matrix
#        (first col = feature, remaining cols = case_ids)
#        Can be gene / iso_log / iso_frac normalised matrices.
#
# mut  : feature × sample matrix for one mutation group
#        (first col = feature_id or similar, remaining cols = case_ids)
#        This is collapsed to a *sample-level* 0/1 covariate:
#          mutation = 1 if ANY mutated feature in this group for that sample
#
# mani : clinical + covariates, with at least:
#        case_id, OS_time, OS_event
#
# Outputs:
#   - <out_results>_full.csv  (all valid fits)
#   - <out_results>.csv       (significant interactions only, FDR_int<0.05)
#   Each row contains main expr, main mut, interaction terms + ΔLL per term
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(yaml)
  library(parallel)
})

say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }

# ------------------------------------------------------------
# Args (support 6 or 7 so launcher with tx2gene still works)
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  die(paste(
    "Usage:",
    "03d_expr_mut_univariate_coxph.R",
    "<expr.csv> <mutation.csv> <manifest.csv> <cov_yaml>",
    "[optional_tx2gene.csv]",
    "<out_results.csv> <out_summary.txt>"
  ))
}

expr_in  <- args[[1]]
mut_in   <- args[[2]]
mani_in  <- args[[3]]
cov_yaml <- args[[4]]

if (length(args) == 6) {
  # old style: exactly 6 args
  out_res <- args[[5]]
  out_sum <- args[[6]]
} else {
  # new style: extra tx2gene argument in position 5 (ignored here)
  out_res <- args[[6]]
  out_sum <- args[[7]]
}

# ------------------------------------------------------------
# Meta: cancer, data_type, mut_group, run_id
# ------------------------------------------------------------
cancer    <- sub("^.*TCGA_([A-Z0-9]+)_.*$", "\\1", basename(expr_in))
data_type <- sub("^.*TCGA_[A-Z0-9]+_(.*?)\\.normalized\\.csv$", "\\1", basename(expr_in))

# mutation group from the mutation matrix filename
# examples: ...gene_ensembl_truncating_LOF.csv, ...gene_ensembl_missense_or_inframe.csv
mut_group <- sub("^.*ensembl_(.*)\\.csv$", "\\1", basename(mut_in))
if (identical(mut_group, basename(mut_in))) mut_group <- "unknown"

run_id <- sprintf("[03d][%s][%s][%s]", cancer, data_type, mut_group)

say("=== 03d_expr_mut_interaction_coxph ===")
say("%s Expr:      %s", run_id, expr_in)
say("%s Mutation:  %s", run_id, mut_in)
say("%s Manifest:  %s", run_id, mani_in)
say("%s Cov YAML:  %s", run_id, cov_yaml)
say("%s Out res:   %s", run_id, out_res)
say("%s Out sum:   %s", run_id, out_sum)

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
expr <- fread(expr_in)
mut  <- fread(mut_in)
mani <- fread(mani_in)
cfg  <- yaml::read_yaml(cov_yaml)

# --- expression matrix: feature × sample ---
feat_col_expr <- if ("feature" %in% names(expr)) "feature" else names(expr)[1]
features <- expr[[feat_col_expr]]
expr[[feat_col_expr]] <- NULL
expr_mat <- as.matrix(expr)
rownames(expr_mat) <- features

# --- mutation matrix: feature × sample (same layout) ---
feat_col_mut <- if ("feature" %in% names(mut)) "feature" else names(mut)[1]
mut_features <- mut[[feat_col_mut]]
mut[[feat_col_mut]] <- NULL
mut_mat <- as.matrix(mut)
rownames(mut_mat) <- mut_features

# --- manifest check ---
needed_mani <- c("case_id", "OS_time", "OS_event")
if (!all(needed_mani %in% names(mani))) {
  die("Manifest must contain columns: %s", paste(needed_mani, collapse = ", "))
}
mani <- mani[complete.cases(mani[, .(OS_time, OS_event)]), ]

# ------------------------------------------------------------
# Align sample IDs (case_ids are column names)
# ------------------------------------------------------------
common <- Reduce(
  intersect,
  list(colnames(expr_mat), colnames(mut_mat), mani$case_id)
)

if (length(common) < 20) {
  die(sprintf("%s Too few overlapping samples: %d", run_id, length(common)))
}

common <- sort(common)

expr_mat <- expr_mat[, common, drop = FALSE]
mut_mat  <- mut_mat[,  common, drop = FALSE]
mani     <- mani[match(common, mani$case_id)]

say("%s Aligned samples: %d", run_id, length(common))

# ------------------------------------------------------------
# Build mutation covariate (sample-level 0/1)
#   1 if ANY mutation in this group for that sample
# ------------------------------------------------------------
mut_mat[is.na(mut_mat)] <- 0
mut_mat <- (mut_mat > 0) * 1

mut_vec <- apply(mut_mat, 2, function(x) as.integer(any(x == 1)))
mani$mutation <- mut_vec

events_mut  <- sum(mani$OS_event[mani$mutation == 1] == 1, na.rm = TRUE)
events_wild <- sum(mani$OS_event[mani$mutation == 0] == 1, na.rm = TRUE)
say("%s Mutation covariate (%s): %d mutated / %d wild; events mutated=%d, wild=%d",
    run_id, mut_group,
    sum(mani$mutation == 1), sum(mani$mutation == 0),
    events_mut, events_wild)

if (events_mut < 5 || events_wild < 5) {
  die(sprintf("%s Too few events in one of the mutation groups (mut=%d, wt=%d)",
              run_id, events_mut, events_wild))
}

# ------------------------------------------------------------
# Covariates (baseline + conditional, with coverage checks)
# ------------------------------------------------------------
covariates <- cfg$baseline_covariates

if (!is.null(cfg$conditional_covariates)) {
  for (cond in cfg$conditional_covariates) {
    nm  <- cond$name
    thr <- cond$coverage_threshold
    if (nm %in% names(mani)) {
      cov_coverage <- mean(!is.na(mani[[nm]]))
      if (!is.na(cov_coverage) && cov_coverage >= thr) {
        covariates <- c(covariates, nm)
      }
    }
  }
}

covariates <- unique(covariates)
covariates <- covariates[covariates %in% names(mani)]

Xcov <- mani[, ..covariates]

# drop covariates with no variability
Xcov <- Xcov[, lapply(.SD, function(x) {
  if (length(unique(na.omit(x))) > 1) x else NULL
})]

cov_used <- names(Xcov)
say("%s Covariates used: %s", run_id,
    if (length(cov_used)) paste(cov_used, collapse = ", ") else "<none>")

# ------------------------------------------------------------
# Build Cox formula once
#   Surv(OS_time, OS_event) ~ expr * mutation + covariates
# ------------------------------------------------------------
if (length(cov_used) > 0) {
  cov_terms  <- paste(cov_used, collapse = " + ")
  fml_string <- paste0("Surv(OS_time, OS_event) ~ expr * mutation + ", cov_terms)
} else {
  fml_string <- "Surv(OS_time, OS_event) ~ expr * mutation"
}
fml <- as.formula(fml_string)
say("%s Cox formula: %s", run_id, fml_string)

# ------------------------------------------------------------
# Parallel settings
# ------------------------------------------------------------
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
if (is.na(ncores) || ncores < 1) ncores <- 2
say("%s Using %d cores", run_id, ncores)

# ------------------------------------------------------------
# Fit interaction model per feature
# ------------------------------------------------------------
idx <- seq_len(nrow(expr_mat))

res_list <- mclapply(idx, function(i) {
  tryCatch({
    vals <- as.numeric(expr_mat[i, ])

    # skip features with no variation
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0)
      return(NULL)

    df <- data.frame(
      OS_time  = mani$OS_time,
      OS_event = mani$OS_event,
      expr     = as.numeric(scale(vals)),   # 1 SD scale
      mutation = mani$mutation,
      Xcov
    )

    fit <- suppressWarnings(
      try(coxph(fml, data = df, ties = "efron"), silent = TRUE)
    )
    if (inherits(fit, "try-error"))
      return(NULL)

    s  <- summary(fit)
    co <- s$coef

    # safe extractor for term rows
    get_row <- function(term) {
      rn <- rownames(co)
      if (!is.null(rn) && term %in% rn) co[term, ] else rep(NA_real_, ncol(co))
    }

    expr_row <- get_row("expr")
    mut_row  <- get_row("mutation")
    int_row  <- get_row("expr:mutation")

    # convenience: function for ΔLL from Wald z
    delta_from_row <- function(row) {
      if (is.null(row) || is.na(row["z"])) return(NA_real_)
      z <- as.numeric(row["z"])
      if (is.na(z)) return(NA_real_)
      (z * z) / 2
    }

    data.frame(
      feature   = rownames(expr_mat)[i],

      beta_expr = expr_row["coef"],
      HR_expr   = exp(expr_row["coef"]),
      p_expr    = expr_row["Pr(>|z|)"],
      delta_LL_expr = delta_from_row(expr_row),

      beta_mut  = mut_row["coef"],
      HR_mut    = exp(mut_row["coef"]),
      p_mut     = mut_row["Pr(>|z|)"],
      delta_LL_mut  = delta_from_row(mut_row),

      beta_int  = int_row["coef"],
      HR_int    = exp(int_row["coef"]),
      p_int     = int_row["Pr(>|z|)"],
      delta_LL_int  = delta_from_row(int_row)
    )
  }, error = function(e) {
    NULL
  })
}, mc.cores = ncores)

res_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, res_list)
if (length(res_list) == 0) die(sprintf("%s No valid Cox models fitted.", run_id))

res <- rbindlist(res_list, fill = TRUE)

# ------------------------------------------------------------
# Add FDRs per term + meta (cancer, data_type, mut_group)
# ------------------------------------------------------------
res[, FDR_expr := p.adjust(p_expr, "BH")]
res[, FDR_mut  := p.adjust(p_mut,  "BH")]
res[, FDR_int  := p.adjust(p_int,  "BH")]

res[, cancer    := cancer]
res[, data_type := data_type]
res[, mut_group := mut_group]

# ------------------------------------------------------------
# Save full results
# ------------------------------------------------------------
out_res_full <- sub("\\.cox_results\\.csv$", ".cox_results_full.csv", out_res)
if (identical(out_res_full, out_res)) {
  out_res_full <- sub("\\.csv$", "_full.csv", out_res)
}
dir.create(dirname(out_res_full), recursive = TRUE, showWarnings = FALSE)
fwrite(res, out_res_full)
say("%s Wrote full results → %s", run_id, out_res_full)

# ------------------------------------------------------------
# FILTER SIGNIFICANT FEATURES ONLY (by interaction term)
# ------------------------------------------------------------
res_sig <- res[
  !is.na(FDR_int) & FDR_int < 0.05 &
  is.finite(HR_int) & HR_int > 0
]

dir.create(dirname(out_res), recursive = TRUE, showWarnings = FALSE)
fwrite(res_sig, out_res)
say("%s Wrote significant-only results (FDR_int < 0.05) → %s",
    run_id, out_res)

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------
sink(out_sum)
cat("=== 03d_expr_mut_interaction_coxph summary ===\n")
cat("Cancer:", cancer, "\n")
cat("Datatype (expression):", data_type, "\n")
cat("Mutation group:", mut_group, "\n")
cat("Features tested:", nrow(expr_mat), "\n")
cat("Features with valid fit:", nrow(res), "\n\n")

cat("Significant expr term  (FDR_expr < 0.05):", sum(res$FDR_expr < 0.05, na.rm = TRUE), "\n")
cat("Significant mut term   (FDR_mut  < 0.05):", sum(res$FDR_mut  < 0.05,  na.rm = TRUE), "\n")
cat("Significant expr:mut   (FDR_int  < 0.05):", sum(res$FDR_int  < 0.05,  na.rm = TRUE), "\n")
cat("Significant rows saved (FDR_int  < 0.05):", nrow(res_sig), "\n")
sink()

say("%s [DONE] Completed interaction CoxPH", run_id)

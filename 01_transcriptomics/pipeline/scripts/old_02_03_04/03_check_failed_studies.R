#!/usr/bin/env Rscript
# ============================================================
# check_failed_studies.R
# ------------------------------------------------------------
# Diagnoses why certain cancers failed ("no_valid_fits" or
# "nonfinite HR") in 03_univariate_coxph.R.
#
# It checks for:
#   - Sample overlap between expression & manifest
#   - Number of survival events (OS_event)
#   - Expression variance across features
#   - Flags likely causes of failure
#
# Outputs:
#   01_transcriptomics/out/03_univariate_coxph/failed_study_diagnostics.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

say <- function(...) message(sprintf(...))

# ============================================================
# 1. Load cohort summary
# ============================================================
cohort_path <- "01_transcriptomics/out/03_univariate_coxph/cohort_summary.csv"
if (!file.exists(cohort_path)) stop("Cohort summary not found: ", cohort_path)
cohort <- fread(cohort_path)

# Keep only failed cases
failed <- cohort[failure_reason %in% c("no_valid_fits", "nonfinite HR")]
say("[INFO] Found %d failed studies to diagnose.", nrow(failed))
if (nrow(failed) == 0) quit(save = "no")

# ============================================================
# 2. Iterate over failed studies
# ============================================================
diagnostics <- list()

for (i in seq_len(nrow(failed))) {
  cancer <- failed$cancer[i]
  dt <- failed$data_type[i]

  expr_file <- sprintf("01_transcriptomics/out/02_norm/TCGA_%s_%s.normalized.csv", cancer, dt)
  mani_file <- sprintf("01_transcriptomics/out/02_norm/TCGA_%s_%s.sample_manifest.csv", cancer, dt)

  if (!file.exists(expr_file) || !file.exists(mani_file)) {
    say("[WARN] %s_%s: Missing expression or manifest file — skipping.", cancer, dt)
    next
  }

  expr <- fread(expr_file)
  mani <- fread(mani_file)

  # Identify IDs and overlap
  feat_col <- if ("feature" %in% names(expr)) "feature" else names(expr)[1]
  expr_ids <- setdiff(colnames(expr), feat_col)
  mani_ids <- mani$case_id
  overlap_ids <- intersect(expr_ids, mani_ids)
  overlap_n <- length(overlap_ids)

  # Event count
  events <- if (all(c("OS_time", "OS_event") %in% names(mani))) {
    mani <- mani[complete.cases(mani[, .(OS_time, OS_event)])]
    sum(mani$OS_event == 1, na.rm = TRUE)
  } else NA_integer_

  # Expression variance
  if (length(expr_ids) > 0) {
    expr_mat <- as.matrix(expr[, ..expr_ids])
    sds <- apply(expr_mat, 1, sd, na.rm = TRUE)
    median_sd <- median(sds, na.rm = TRUE)
    pct_zero_var <- mean(sds == 0, na.rm = TRUE) * 100
  } else {
    median_sd <- NA
    pct_zero_var <- NA
  }

  # ============================================================
  # 3. Flag likely cause
  # ============================================================
  if (is.na(events) || events == 0) {
    cause <- "no_events"
  } else if (overlap_n < 20) {
    cause <- "too_few_samples"
  } else if (!is.na(pct_zero_var) && pct_zero_var > 50) {
    cause <- "flat_expression_data"
  } else if (median_sd < 0.01) {
    cause <- "very_low_variance"
  } else if (nrow(expr) == 0) {
    cause <- "empty_expression_matrix"
  }
    else {
        cause <- "unknown"
    }

  # Collect diagnostic row
  diagnostics[[length(diagnostics) + 1]] <- data.table(
    cancer = cancer,
    data_type = dt,
    overlap_n = overlap_n,
    n_samples_expr = length(expr_ids),
    n_samples_mani = length(mani_ids),
    n_events = events,
    median_sd = round(median_sd, 5),
    pct_zero_var = round(pct_zero_var, 2),
    likely_cause = cause
  )

  say("[CHECK] %s_%s → overlap=%d | events=%s | median SD=%.4f | zero-var=%.1f%% | cause=%s",
      cancer, dt, overlap_n, ifelse(is.na(events), "NA", events),
      median_sd, pct_zero_var, cause)
}

# ============================================================
# 4. Save combined diagnostics table
# ============================================================
diag_dt <- rbindlist(diagnostics, fill = TRUE)
out_file <- "01_transcriptomics/out/03_univariate_coxph/failed_study_diagnostics.csv"
fwrite(diag_dt, out_file)

say("[DONE] Diagnostics summary written to: %s", out_file)

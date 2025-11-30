#!/usr/bin/env Rscript
# ============================================================
# 03_check_successful_studies.R
# ------------------------------------------------------------
# Diagnostics for *succeeded* CoxPH runs from 03_univariate_coxph.R.
# Computes metrics comparable to check_failed_studies.R:
#   - sample overlap
#   - number of survival events
#   - expression variation
#
# Goal: detect inconsistencies where one data type succeeded
# while others from same cancer failed (despite shared clinical data).
#
# Output:
#   01_transcriptomics/out/03_univariate_coxph/successful_study_diagnostics.csv
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

# Filter successful studies
succeeded <- cohort[failure_reason == "ok"]
say("[INFO] Found %d successful studies to analyze.", nrow(succeeded))
if (nrow(succeeded) == 0) quit(save = "no")

# ============================================================
# 2. Iterate over succeeded studies
# ============================================================
diagnostics <- list()

for (i in seq_len(nrow(succeeded))) {
  cancer <- succeeded$cancer[i]
  dt <- succeeded$data_type[i]

  expr_file <- sprintf("01_transcriptomics/out/02_norm/TCGA_%s_%s.normalized.csv", cancer, dt)
  mani_file <- sprintf("01_transcriptomics/out/02_norm/TCGA_%s_%s.sample_manifest.csv", cancer, dt)

  if (!file.exists(expr_file) || !file.exists(mani_file)) {
    say("[WARN] %s_%s: Missing expression or manifest file — skipping.", cancer, dt)
    next
  }

  expr <- fread(expr_file)
  mani <- fread(mani_file)

  # Identify overlap
  feat_col <- if ("feature" %in% names(expr)) "feature" else names(expr)[1]
  expr_ids <- setdiff(colnames(expr), feat_col)
  mani_ids <- mani$case_id
  overlap_ids <- intersect(expr_ids, mani_ids)
  overlap_n <- length(overlap_ids)

  # Count events
  events <- if (all(c("OS_time", "OS_event") %in% names(mani))) {
    mani <- mani[complete.cases(mani[, .(OS_time, OS_event)])]
    sum(mani$OS_event == 1, na.rm = TRUE)
  } else NA_integer_

  # Expression variance metrics
  if (length(expr_ids) > 0) {
    expr_mat <- as.matrix(expr[, ..expr_ids])
    sds <- apply(expr_mat, 1, sd, na.rm = TRUE)
    median_sd <- median(sds, na.rm = TRUE)
    pct_zero_var <- mean(sds == 0, na.rm = TRUE) * 100
  } else {
    median_sd <- NA
    pct_zero_var <- NA
  }

  # Add diagnostic row
  diagnostics[[length(diagnostics) + 1]] <- data.table(
    cancer = cancer,
    data_type = dt,
    overlap_n = overlap_n,
    n_samples_expr = length(expr_ids),
    n_samples_mani = length(mani_ids),
    n_events = events,
    median_sd = round(median_sd, 5),
    pct_zero_var = round(pct_zero_var, 2)
  )

  say("[CHECK] %s_%s → overlap=%d | events=%s | median SD=%.4f | zero-var=%.1f%%",
      cancer, dt, overlap_n, ifelse(is.na(events), "NA", events),
      median_sd, pct_zero_var)
}

# ============================================================
# 3. Save results
# ============================================================
diag_dt <- rbindlist(diagnostics, fill = TRUE)
out_file <- "01_transcriptomics/out/03_univariate_coxph/successful_study_diagnostics.csv"
fwrite(diag_dt, out_file)

say("[DONE] Successful study diagnostics written to: %s", out_file)

# ============================================================
# 4. Optional: cross-check for inconsistencies
# ============================================================
say("[INFO] Checking for potential inconsistencies...")

failed_file <- "01_transcriptomics/out/03_univariate_coxph/failed_study_diagnostics.csv"

if (file.exists(failed_file)) {
  failed <- fread(failed_file)

  combined <- merge(
    diag_dt[, .(cancer, data_type, n_events, median_sd, pct_zero_var)],
    failed[, .(cancer, data_type, n_events, median_sd, pct_zero_var, likely_cause)],
    by = "cancer", all = TRUE, suffixes = c("_success", "_fail")
  )

  # --- Safely flag inconsistent event counts ---
  combined[, consistent_events := abs(n_events_success - n_events_fail) <= 2 | is.na(n_events_fail)]

  # ✅ use parentheses to ensure data.table evaluates inside itself
  combined[(!consistent_events), flag := "MISMATCH_EVENT_COUNTS"]

  # --- Add extra diagnostics ---
  combined[, flag := fifelse(is.na(flag), "", flag)]
  combined[, consistency_note := fifelse(
    flag == "MISMATCH_EVENT_COUNTS",
    sprintf("Different #events (success=%s vs fail=%s)", n_events_success, n_events_fail),
    "Consistent"
  )]

  # --- Save combined output ---
  out_consistency <- "01_transcriptomics/out/03_univariate_coxph/study_consistency_overview.csv"
  fwrite(combined, out_consistency)
  say("[DONE] Cross-type consistency overview saved: %s", out_consistency)
} else {
  say("[SKIP] No failed diagnostics file found — skipping consistency check.")
}
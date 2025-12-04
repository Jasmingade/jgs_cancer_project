#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 4) {
  stop("Usage: 00_filter_cox_ci.R <cox_root> <out_root> [ci_width_threshold] [max_abs_log2HR]")
}

cox_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_root <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)

ci_thresh <- if (length(args) >= 3) as.numeric(args[[3]]) else 10
if (!is.finite(ci_thresh) || ci_thresh <= 0) {
  stop("ci_width_threshold must be a positive numeric value.")
}

log2_thresh <- if (length(args) == 4) as.numeric(args[[4]]) else 10
if (!is.finite(log2_thresh) || log2_thresh <= 0) {
  stop("max_abs_log2HR must be a positive numeric value.")
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

say <- function(fmt, ...) message(sprintf(fmt, ...))

files <- list.files(cox_root, pattern = "cox_results\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(files)) {
  say("[filter] No cox_results.csv files found under %s", cox_root)
  quit(status = 0)
}

agg_list <- list()
base_len <- nchar(cox_root)
processed <- 0L
kept_rows <- 0L

for (f in sort(files)) {
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt)) {
    say("[filter] Skipping unreadable file: %s", f)
    next
  }
  if (!all(c("beta", "se", "HR") %in% names(dt))) {
    say("[filter] Skipping %s (missing beta/se/HR)", f)
    next
  }

  # Cast necessary columns to numeric
  dt[, beta := as.numeric(beta)]
  dt[, se   := as.numeric(se)]
  dt[, HR   := as.numeric(HR)]

  dt <- dt[is.finite(beta) & is.finite(se) & is.finite(HR) & HR > 0]
  if (!nrow(dt)) next

  processed <- processed + nrow(dt)

  # Compute CI
  dt[, ci_low   := beta - 1.96 * se]
  dt[, ci_high  := beta + 1.96 * se]
  dt[, ci_width := ci_high - ci_low]

  # Compute log2(HR)
  dt[, log2HR := log2(HR)]

  # Stable = passing both CI and log2(HR) filters
  stable <- dt[
    is.finite(ci_width) &
    ci_width <= ci_thresh &
    is.finite(log2HR) &
    abs(log2HR) <= log2_thresh
  ]

  kept_rows <- kept_rows + nrow(stable)

  rel_path <- substring(normalizePath(f, winslash = "/", mustWork = TRUE), base_len + 2)
  out_file <- file.path(out_root, rel_path)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  fwrite(stable, out_file)

  if (nrow(stable)) {
    agg_list[[length(agg_list) + 1L]] <- cbind(data.table(source_file = rel_path), stable)
  }
}

if (length(agg_list)) {
  combined <- rbindlist(agg_list, fill = TRUE)
  fwrite(combined, file.path(out_root, "significant_results_stable_ci_log2HR.csv"))
  say("[filter] Wrote %d total stable rows to %s",
      nrow(combined),
      file.path(out_root, "significant_results_stable_ci_log2HR.csv"))
} else {
  fwrite(data.table(), file.path(out_root, "significant_results_stable_ci_log2HR.csv"))
  say("[filter] No stable rows found; wrote empty summary.")
}

say("[filter] Processed %d rows across %d files; kept %d stable rows (ci_width <= %.2f & |log2HR| <= %.2f)",
    processed, length(files), kept_rows, ci_thresh, log2_thresh)

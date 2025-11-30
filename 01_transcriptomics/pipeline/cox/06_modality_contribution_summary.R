#!/usr/bin/env Rscript
# ============================================================
# 06_modality_contribution_summary.R
# ------------------------------------------------------------
# Uses *_full outputs from:
#   03a_univariate_coxph (expression-only)
#   03b_mutation_univariate_coxph (mutation-only)
#   03c_exp_mutation_univariate_coxph (expr+mut-group)
#
# For each cancer and modality, computes:
#   - approximate Δ log-likelihood (delta_LL) from Wald chi-square
#   - per-modality share of total delta_LL (contribution fraction)
#
# This is heuristic, not a formal causal estimate, but it gives:
#   “Which modality accounts for more of the survival signal in each cancer?”
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

say <- function(...) message(sprintf(...))
strip_version <- function(x) sub("\\.\\d+$", "", x)

root <- "01_transcriptomics/out"

# ------------------------------------------------------------
# Helper: compute delta_LL from p-value
# ------------------------------------------------------------
compute_delta_LL <- function(p) {
  p <- as.numeric(p)
  p <- p[is.finite(p) & p > 0 & p < 1]
  if (length(p) == 0) return(numeric(0))
  chisq <- qchisq(1 - p, df = 1)
  chisq / 2
}

# ------------------------------------------------------------
# 1) Load 03a FULL (expression-only)
# ------------------------------------------------------------
dir_03a <- file.path(root, "03a_univariate_coxph")

files_03a <- list.files(
  dir_03a,
  pattern = "cox_results_full\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

say("[INFO] Found %d 03a full result files", length(files_03a))

load_03a <- function(f) {
  dt <- fread(f)
  if (!all(c("feature", "beta", "HR", "p", "cancer", "data_type") %in% names(dt))) {
    warning(sprintf("[WARN] Skipping 03a file (missing cols): %s", f))
    return(NULL)
  }
  dt[, source := "03a_expr_only"]
  dt[, term   := "expr"]   # expression term
  dt[, modality := data_type]   # gene / iso_log / iso_frac
  dt[, delta_LL := compute_delta_LL(p)]
  dt
}

res_03a_list <- lapply(files_03a, load_03a)
res_03a_list <- Filter(Negate(is.null), res_03a_list)
res_03a <- if (length(res_03a_list)) rbindlist(res_03a_list, fill = TRUE) else data.table()

# ------------------------------------------------------------
# 2) Load 03b FULL (mutation-only)
# ------------------------------------------------------------
dir_03b <- file.path(root, "03b_mutation_univariate_coxph")

files_03b <- list.files(
  dir_03b,
  pattern = "cox_results_full\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

say("[INFO] Found %d 03b full result files", length(files_03b))

load_03b <- function(f) {
  dt <- fread(f)
  if (!all(c("feature", "beta", "HR", "p", "cancer", "mut_group") %in% names(dt))) {
    warning(sprintf("[WARN] Skipping 03b file (missing cols): %s", f))
    return(NULL)
  }
  dt[, source := "03b_mut_only"]
  dt[, term   := "mut"]
  dt[, modality := paste0("mut_", mut_group)]  # e.g. mut_missense_or_inframe
  dt[, delta_LL := compute_delta_LL(p)]
  dt
}

res_03b_list <- lapply(files_03b, load_03b)
res_03b_list <- Filter(Negate(is.null), res_03b_list)
res_03b <- if (length(res_03b_list)) rbindlist(res_03b_list, fill = TRUE) else data.table()

# ------------------------------------------------------------
# 3) Load 03c FULL (expr + mut-group)
# ------------------------------------------------------------
dir_03c <- file.path(root, "03c_exp_mutation_univariate_coxph")

files_03c <- list.files(
  dir_03c,
  pattern = "cox_results_full\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

say("[INFO] Found %d 03c full result files", length(files_03c))

# parse cancer, data_type, mut_group from filename
parse_03c_filename <- function(f) {
  base <- sub("\\.cox_results_full\\.csv$", "", basename(f))
  m <- regexec("^TCGA_([A-Z0-9]+)_(gene|iso_log|iso_frac)_(.+)$", base)
  mm <- regmatches(base, m)[[1]]
  if (length(mm) == 0) {
    return(list(cancer = NA_character_, data_type = NA_character_, mut_group = NA_character_))
  }
  list(cancer = mm[2], data_type = mm[3], mut_group = mm[4])
}

load_03c <- function(f) {
  dt <- fread(f)
  needed <- c("feature", "beta_expr", "HR_expr", "p_expr",
              "beta_mut",  "HR_mut",  "p_mut")
  if (!all(needed %in% names(dt))) {
    warning(sprintf("[WARN] Skipping 03c file (missing cols): %s", f))
    return(NULL)
  }

  info <- parse_03c_filename(f)
  dt[, cancer    := info$cancer]
  dt[, data_type := info$data_type]
  dt[, mut_group := info$mut_group]

  dt[, source := "03c_expr_mut"]

  # expr term row
  expr_dt <- dt[, .(
    feature, cancer, data_type, mut_group,
    term     = "expr",
    modality = data_type,
    p        = as.numeric(p_expr),
    HR       = as.numeric(HR_expr)
  )]

  # mut-group term row
  mut_dt <- dt[, .(
    feature, cancer, data_type, mut_group,
    term     = "mut",
    modality = paste0("mutgroup_", mut_group),
    p        = as.numeric(p_mut),
    HR       = as.numeric(HR_mut)
  )]

  combo <- rbind(expr_dt, mut_dt, fill = TRUE)
  combo[, delta_LL := compute_delta_LL(p)]

  combo
}

res_03c_list <- lapply(files_03c, load_03c)
res_03c_list <- Filter(Negate(is.null), res_03c_list)
res_03c <- if (length(res_03c_list)) rbindlist(res_03c_list, fill = TRUE) else data.table()

# ------------------------------------------------------------
# 4) Combine all terms and summarise
# ------------------------------------------------------------
all_terms <- rbindlist(list(res_03a, res_03b, res_03c), fill = TRUE)

# keep only finite delta_LL
all_terms <- all_terms[is.finite(delta_LL) & delta_LL > 0]

say("[INFO] Combined table has %d rows", nrow(all_terms))

# Per cancer × modality summary
contrib <- all_terms[, .(
  n_terms       = .N,
  sum_delta_LL  = sum(delta_LL, na.rm = TRUE)
), by = .(cancer, modality)]

# normalise within each cancer
contrib[, total_delta_LL := sum(sum_delta_LL, na.rm = TRUE), by = cancer]
contrib[, frac_contrib := ifelse(total_delta_LL > 0, sum_delta_LL / total_delta_LL, NA_real_)]

# Save
out_file <- file.path(root, "06_modality_contribution_summary.csv")
fwrite(contrib, out_file)

say("[DONE] Saved modality contribution summary → %s", out_file)

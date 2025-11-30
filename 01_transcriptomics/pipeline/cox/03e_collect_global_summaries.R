#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

root <- "01_transcriptomics/out"

# ---------- 03a: Expression-only (per cancer × data_type) ----------
res_files_03a <- list.files(
  file.path(root, "03a_univariate_coxph"),
  pattern = "\\.cox_results\\.csv$", recursive = TRUE, full.names = TRUE
)

if (length(res_files_03a)) {
  summ_list_03a <- lapply(res_files_03a, function(f) {
    dt <- fread(f)

    base <- basename(f)
    cancer <- sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", base)
    data_type <- sub("^TCGA_[A-Z0-9]+_(.*?)\\.cox_results\\.csv$", "\\1", base)

    n_features_saved <- nrow(dt)
    n_FDR_lt_0_05 <- if ("FDR" %in% names(dt)) {
      sum(dt$FDR < 0.05, na.rm = TRUE)
    } else {
      NA_integer_
    }

    data.table(
      cancer           = cancer,
      data_type        = data_type,
      n_features_saved = n_features_saved,
      n_FDR_lt_0_05    = n_FDR_lt_0_05,
      results_file     = f
    )
  })

  summ_03a <- rbindlist(summ_list_03a, use.names = TRUE, fill = TRUE)
  fwrite(summ_03a, file.path(root, "03a_univariate_coxph", "03a_global_summary.csv"))
}

# ---------- 03b: Mutation-only (per cancer × mut_group) ----------
res_files_03b <- list.files(
  file.path(root, "03b_mutation_univariate_coxph"),
  pattern = "\\.cox_results\\.csv$", recursive = TRUE, full.names = TRUE
)

if (length(res_files_03b)) {
  summ_list_03b <- lapply(res_files_03b, function(f) {
    dt <- fread(f)

    base <- basename(f)
    tmp  <- sub("^TCGA_", "", base)
    tmp  <- sub("\\.cox_results\\.csv$", "", tmp)
    cancer    <- sub("^([^_]+)_.*$", "\\1", tmp)
    mut_group <- sub("^[^_]+_mutation_", "", tmp)

    n_features_saved <- nrow(dt)
    n_FDR_lt_0_05 <- if ("FDR" %in% names(dt)) {
      sum(dt$FDR < 0.05, na.rm = TRUE)
    } else {
      NA_integer_
    }

    data.table(
      cancer           = cancer,
      mut_group        = mut_group,
      n_features_saved = n_features_saved,
      n_FDR_lt_0_05    = n_FDR_lt_0_05,
      results_file     = f
    )
  })

  summ_03b <- rbindlist(summ_list_03b, use.names = TRUE, fill = TRUE)
  fwrite(summ_03b, file.path(root, "03b_mutation_univariate_coxph", "03b_global_summary.csv"))
}

# ---------- helper for 03c/03d filename parsing ----------
# Expect results like:
#   03c: TCGA_<CANCER>_<DTYPE>_<MUTGROUP>.cox_results.csv
#   03d: TCGA_<CANCER>_<DTYPE>_<MUTGROUP>.cox_results.csv
parse_cancer_dtype_mut_from_results <- function(f) {
  base <- basename(f)
  base <- sub("^TCGA_", "", base)
  base <- sub("\\.cox_results\\.csv$", "", base)  # now: <CANCER>_<DTYPE>_<MUTGROUP>

  cancer <- sub("^([^_]+)_.*$", "\\1", base)

  data_type <- NA_character_
  mut_group <- NA_character_

  if (grepl("_gene_", base)) {
    data_type <- "gene"
    mut_group <- sub("^[^_]+_gene_", "", base)
  } else if (grepl("_iso_log_", base)) {
    data_type <- "iso_log"
    mut_group <- sub("^[^_]+_iso_log_", "", base)
  } else if (grepl("_iso_frac_", base)) {
    data_type <- "iso_frac"
    mut_group <- sub("^[^_]+_iso_frac_", "", base)
  }

  list(cancer = cancer, data_type = data_type, mut_group = mut_group)
}

# ---------- 03c: Expression + Mutation (per cancer × data_type × mut_group) ----------
res_files_03c <- list.files(
  file.path(root, "03c_exp_mutation_univariate_coxph"),
  pattern = "\\.cox_results\\.csv$", recursive = TRUE, full.names = TRUE
)

if (length(res_files_03c)) {
  summ_list_03c <- lapply(res_files_03c, function(f) {
    dt <- fread(f)
    id <- parse_cancer_dtype_mut_from_results(f)

    n_features_saved <- nrow(dt)

    # These only work if the columns exist in your saved output
    n_expr_FDR_lt_0_05 <- if ("FDR_expr" %in% names(dt)) {
      sum(dt$FDR_expr < 0.05, na.rm = TRUE)
    } else {
      NA_integer_
    }

    n_mut_FDR_lt_0_05 <- if ("FDR_mut" %in% names(dt)) {
      sum(dt$FDR_mut < 0.05, na.rm = TRUE)
    } else {
      NA_integer_
    }

    n_both_FDR_lt_0_05 <- if (all(c("FDR_expr", "FDR_mut") %in% names(dt))) {
      sum(dt$FDR_expr < 0.05 & dt$FDR_mut < 0.05, na.rm = TRUE)
    } else {
      NA_integer_
    }

    data.table(
      cancer              = id$cancer,
      data_type           = id$data_type,
      mut_group           = id$mut_group,
      n_features_saved    = n_features_saved,
      n_expr_FDR_lt_0_05  = n_expr_FDR_lt_0_05,
      n_mut_FDR_lt_0_05   = n_mut_FDR_lt_0_05,
      n_both_FDR_lt_0_05  = n_both_FDR_lt_0_05,
      results_file        = f
    )
  })

  summ_03c <- rbindlist(summ_list_03c, use.names = TRUE, fill = TRUE)
  fwrite(summ_03c, file.path(root, "03c_exp_mutation_univariate_coxph", "03c_global_summary.csv"))
}

# ---------- 03d: Isoform × Mutation interaction ----------
res_files_03d <- list.files(
  file.path(root, "03d_iso_mut_univariate_coxph"),
  pattern = "\\.cox_results\\.csv$", recursive = TRUE, full.names = TRUE
)

if (length(res_files_03d)) {
  summ_list_03d <- lapply(res_files_03d, function(f) {
    dt <- fread(f)
    id <- parse_cancer_dtype_mut_from_results(f)

    n_features_saved <- nrow(dt)

    n_FDR_lt_0_05 <- if ("FDR" %in% names(dt)) {
      sum(dt$FDR < 0.05, na.rm = TRUE)
    } else {
      NA_integer_
    }

    median_HR <- if ("HR" %in% names(dt)) {
      median(dt$HR, na.rm = TRUE)
    } else {
      NA_real_
    }

    data.table(
      cancer           = id$cancer,
      data_type        = id$data_type,
      mut_group        = id$mut_group,
      n_features_saved = n_features_saved,
      n_FDR_lt_0_05    = n_FDR_lt_0_05,
      median_HR        = median_HR,
      results_file     = f
    )
  })

  summ_03d <- rbindlist(summ_list_03d, use.names = TRUE, fill = TRUE)
  fwrite(summ_03d, file.path(root, "03d_iso_mut_univariate_coxph", "03d_global_summary.csv"))
}

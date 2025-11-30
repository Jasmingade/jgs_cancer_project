#!/usr/bin/env Rscript

# Perform RUV-based batch correction per study/platform/sample_type using the
# combined per-study CSV layout:
#   <input_root>/<dtype>/<dataset_id>_<dtype>.csv
# The corrected matrices are written to:
#   <output_root>/<dtype>/<dataset_id>_<dtype>.csv
# Optional iso_frac matrices are regenerated from iso_log.

suppressPackageStartupMessages({
  library(data.table)
})
if (!requireNamespace("ruv", quietly = TRUE)) {
  stop("Package 'ruv' is required. Install via install.packages('ruv').")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 6) {
  stop("Usage: 03_ruv_batch_correct.R <input_root> <output_root> <k> [comma_separated_dtypes] [target_dataset] [batch_annotation_dir]")
}

input_root <- args[[1]]
output_root <- args[[2]]
k <- as.integer(args[[3]])
if (is.na(k) || k < 1) stop("k must be a positive integer.")
dtype_arg <- if (length(args) >= 4) args[[4]] else "gene,iso_log"
target_dataset <- if (length(args) >= 5) args[[5]] else NA_character_
batch_dir <- if (length(args) == 6) args[[6]] else "02_proteomics/data/batch_annotation"

dtypes <- unique(trimws(strsplit(dtype_arg, ",")[[1]]))
dtypes <- dtypes[nzchar(dtypes)]
if (!length(dtypes)) stop("No dtypes specified.")

if (!dir.exists(input_root)) stop("Input root not found: ", input_root)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

say <- function(fmt, ...) cat(sprintf(paste0("[ruv] ", fmt, "\n"), ...))

canonicalize_single <- function(val) {
  if (length(val) != 1L) val <- val[1]
  if (is.null(val) || is.na(val) || !nzchar(val)) return(NA_character_)
  clean <- gsub("[^a-zA-Z0-9_ ]", "", val)
  clean <- gsub(" +", "_", trimws(tolower(clean)))
  parts <- unlist(strsplit(clean, "_", fixed = TRUE))
  parts <- parts[parts != ""]
  if (!length(parts)) return(NA_character_)
  if (length(parts) >= 3) parts <- parts[seq_len(3)]
  paste(parts, collapse = "_")
}

canonicalize_batch <- function(vals) {
  if (is.null(vals)) return(rep(NA_character_, length(vals)))
  vapply(vals, canonicalize_single, character(1), USE.NAMES = FALSE)
}

reference_batch_labels <- function(cols) {
  if (is.null(cols)) return(character())
  vapply(cols, function(col) {
    if (is.null(col) || is.na(col) || !nzchar(col)) return(NA_character_)
    col_clean <- gsub("[^a-zA-Z0-9_ ]", "", col)
    col_clean <- gsub(" +", "_", trimws(tolower(col_clean)))
    suffix <- sub("^[^_]*_", "", col_clean)
    if (!nzchar(suffix)) suffix <- col_clean
    canonicalize_single(suffix)
  }, character(1), USE.NAMES = FALSE)
}

batch_cache <- new.env(parent = emptyenv())

load_case_batches <- function(study_key) {
  if (!nzchar(batch_dir)) return(NULL)
  file <- file.path(batch_dir, paste0(study_key, ".csv"))
  if (!file.exists(file)) {
    pattern <- paste0("^", study_key, "\\.csv$")
    cand <- list.files(batch_dir, pattern = pattern, ignore.case = TRUE, full.names = TRUE)
    if (!length(cand)) return(NULL)
    file <- cand[1]
    say("Using batch annotation file %s for study key %s", basename(file), study_key)
  }
  if (exists(study_key, envir = batch_cache, inherits = FALSE)) {
    return(get(study_key, envir = batch_cache))
  }
  dt <- tryCatch(fread(file), error = function(e) NULL)
  if (is.null(dt) || !"case_id" %in% names(dt) || !"folder_name" %in% names(dt)) {
    return(NULL)
  }
  dt[, case_id := toupper(trimws(case_id))]
  if ("sample_type" %in% names(dt)) {
    dt[, sample_type := canonicalize_batch(sample_type)]
  } else {
    dt[, sample_type := NA_character_]
  }
  dt[, folder_name := canonicalize_batch(folder_name)]
  dt <- dt[nzchar(case_id) & nzchar(folder_name)]
  if (!nrow(dt)) return(NULL)
  dt <- dt[order(case_id, sample_type, folder_name)]
  dt <- dt[!duplicated(data.table(case_id, sample_type))]
  map <- dt[, .(case_id, sample_type, batch_id = folder_name)]
  assign(study_key, map, envir = batch_cache)
  map
}

infer_sample_type <- function(dataset) {
  if (grepl("_normal$", dataset, ignore.case = TRUE)) return("normal")
  if (grepl("_tumor$", dataset, ignore.case = TRUE)) return("tumor")
  NA_character_
}

iso_frac_name <- function(base_name) {
  if (grepl("_iso_log(\\.csv)?$", base_name, ignore.case = TRUE)) {
    sub("_iso_log(\\.csv)?$", "_iso_frac.csv", base_name, ignore.case = TRUE)
  } else if (tolower(base_name) == "iso_log.csv") {
    "iso_frac.csv"
  } else {
    sub("\\.csv$", "_iso_frac.csv", base_name, ignore.case = TRUE)
  }
}

dataset_id <- function(path, dtype) {
  sub(paste0("_", dtype, "\\.csv$"), "", basename(path), ignore.case = TRUE)
}

reference_candidates <- function(file, dtype) {
  base <- sub(paste0("_", dtype, "\\.csv$"), "", basename(file), ignore.case = TRUE)
  dirn <- dirname(file)
  cand <- file.path(dirn, paste0(base, "_reference_", dtype, ".csv"))
  if (grepl("_(normal|tumor)$", base, ignore.case = TRUE)) {
    base2 <- sub("_(normal|tumor)$", "", base, ignore.case = TRUE)
    cand <- c(cand, file.path(dirn, paste0(base2, "_reference_", dtype, ".csv")))
  }
  unique(cand)
}

align_to <- function(mat, target_feats) {
  if (is.null(mat)) return(NULL)
  feat <- rownames(mat)
  keep <- feat %in% target_feats
  mat2 <- mat[keep, , drop = FALSE]
  mat2 <- mat2[match(target_feats, feat[keep]), , drop = FALSE]
  rownames(mat2) <- target_feats
  mat2
}

safe_median <- function(x) if (all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)

process_file <- function(file, dtype) {
  dt <- fread(file)
  if (ncol(dt) < 2) {
    say("Skipping %s (<2 columns)", basename(file))
    return()
  }
  feats <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- feats
  storage.mode(mat) <- "double"
  dataset <- dataset_id(file, dtype)
  dataset_type <- infer_sample_type(dataset)

  say("Starting batch correction prep for %s (%s)", dataset, dtype)
  ref_paths <- reference_candidates(file, dtype)
  existing_refs <- ref_paths[file.exists(ref_paths)]
  ref_file <- if (length(existing_refs)) existing_refs[1] else NA_character_
  has_reference <- !is.na(ref_file)
  ref_mat <- NULL
  ref_cols <- 0L
  if (has_reference) {
    ref_dt <- fread(ref_file)
    if (ncol(ref_dt) < 2) {
      say("Reference %s has <2 columns; ignoring", basename(ref_file))
      has_reference <- FALSE
    } else {
      ref_feats <- ref_dt[[1]]
      ref_mat <- as.matrix(ref_dt[, -1, with = FALSE])
      rownames(ref_mat) <- ref_feats
      storage.mode(ref_mat) <- "double"
      ref_mat <- align_to(ref_mat, rownames(mat))
      colnames(ref_mat) <- paste0("REF_", colnames(ref_mat))
      ref_cols <- ncol(ref_mat)
    }
  }

  full_mat <- if (has_reference) cbind(mat, ref_mat) else mat
  is_reference <- c(rep(FALSE, ncol(mat)), if (has_reference) rep(TRUE, ncol(ref_mat)) else logical())
  sample_names <- colnames(full_mat)
  study_key <- sub("_(normal|tumor|reference)$", "", dataset, ignore.case = TRUE)
  case_batches <- load_case_batches(study_key)
  study_case_ids <- toupper(colnames(mat))
  study_batches <- rep(NA_character_, length(study_case_ids))
  study_cases_matched <- 0L
  if (!is.null(case_batches)) {
    dt_sub <- case_batches
    if (!is.na(dataset_type) && "sample_type" %in% names(case_batches)) {
      sel <- dt_sub[sample_type == dataset_type]
      if (nrow(sel)) dt_sub <- sel
    }
    idx <- match(study_case_ids, dt_sub$case_id)
    study_batches <- dt_sub$batch_id[idx]
    if (anyNA(study_batches) && !is.null(case_batches)) {
      missing_idx <- which(is.na(study_batches))
      if (length(missing_idx)) {
        fallback <- match(study_case_ids[missing_idx], case_batches$case_id)
        study_batches[missing_idx] <- case_batches$batch_id[fallback]
      }
    }
    study_cases_matched <- sum(!is.na(study_batches))
  } else {
    say("Dataset %s has no batch annotation file for %s", dataset, study_key)
  }
  ref_batches <- if (has_reference) reference_batch_labels(colnames(ref_mat)) else character()
  replicate_labels <- c(study_batches, ref_batches)
  fallback_labels <- sample_names
  if (has_reference) fallback_labels[is_reference] <- paste0("reference_", dataset)
  if (all(is.na(replicate_labels))) {
    replicate_labels <- fallback_labels
  } else {
    replicate_labels[is.na(replicate_labels)] <- fallback_labels[is.na(replicate_labels)]
    dup_tmp <- table(replicate_labels)
    if (!any(dup_tmp >= 2)) {
      say("Dataset %s replicate labels from batches have no duplicates; falling back to reference grouping", dataset)
      replicate_labels <- fallback_labels
    }
  }
  say("Dataset %s (%s) samples=%d (batches assigned=%d) refs=%d unique_replicates=%d",
      dataset, dtype, ncol(mat), study_cases_matched, ref_cols, length(unique(replicate_labels)))

  corrected <- mat
  if (ncol(full_mat) < 2) {
    say("Dataset %s has <2 samples; copying input", dataset)
  } else {
    full_mat[is.na(full_mat)] <- 0
    var_source <- if (has_reference && sum(is_reference) >= 2) {
      say("Dataset %s using reference-derived controls", dataset)
      ref_idx <- which(is_reference)
      apply(full_mat[, ref_idx, drop = FALSE], 1, var, na.rm = TRUE)
    } else {
      if (!has_reference) {
        say("Dataset %s lacks reference columns; using fallback controls", dataset)
      } else if (sum(is_reference) < 2) {
        say("Dataset %s has only one reference column; using fallback controls", dataset)
      }
      apply(full_mat, 1, var, na.rm = TRUE)
    }
    control_order <- order(var_source, decreasing = FALSE)
    control_idx <- control_order[seq_len(min(length(control_order), 500))]
    control_idx <- control_idx[control_idx > 0 & !is.na(control_idx)]

    dup_sizes <- table(replicate_labels)
    max_rep_size <- if (length(dup_sizes)) max(dup_sizes) else 1

    if (!length(control_idx)) {
      say("Dataset %s lacks control features; copying input", dataset)
    } else if (max_rep_size < 2) {
      say("Dataset %s has no replicate groups (unique=%d); skipping RUVIII adjustment", dataset, length(dup_sizes))
    } else {
      k_use <- min(k, length(control_idx))
      if (k_use < 1) {
        say("Dataset %s insufficient controls for k=%d; copying input", dataset, k)
      } else {
        Y <- t(full_mat)
        ctl_vec <- control_idx
        M <- ruv::replicate.matrix(replicate_labels)
        fit <- ruv::RUVIII(Y = Y, M = M, ctl = ctl_vec, k = k_use, return.info = TRUE)
        adj <- if (is.list(fit)) fit$newY else fit
        corrected_full <- t(adj)
        rownames(corrected_full) <- rownames(full_mat)
        corrected <- corrected_full[, !is_reference, drop = FALSE]
        say("Dataset %s (%s) corrected with k=%d via RUVIII (%s replicates)",
            dataset, dtype, k_use,
            if (has_reference && sum(is_reference) >= 2) "reference" else "sample duplicates")
        say("Dataset %s (%s) batch correction completed successfully", dataset, dtype)
      }
    }
  }

  out_dir <- file.path(output_root, dtype)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, basename(file))
  out_dt <- data.table(feature = rownames(corrected))
  out_dt <- cbind(out_dt, as.data.table(corrected))
  fwrite(out_dt, out_file)

  if (dtype == "iso_log") {
    eps <- 1e-6
    frac_mat <- plogis(corrected)
    frac_mat <- pmin(pmax(frac_mat, eps), 1 - eps)
    frac_dt <- data.table(feature = rownames(frac_mat))
    frac_dt <- cbind(frac_dt, as.data.table(frac_mat))
    frac_dir <- file.path(output_root, "iso_frac")
    dir.create(frac_dir, recursive = TRUE, showWarnings = FALSE)
    frac_file <- file.path(frac_dir, iso_frac_name(basename(file)))
    fwrite(frac_dt, frac_file)
    say("Derived iso_frac for %s", dataset)
  }
}

for (dtype in dtypes) {
  dtype_root <- file.path(input_root, dtype)
  if (!dir.exists(dtype_root)) {
    say("Missing directory for %s: %s", dtype, dtype_root)
    next
  }
  files <- list.files(dtype_root, pattern = paste0("_", dtype, "\\.csv$"), recursive = TRUE, full.names = TRUE)
  if (!length(files)) {
    say("No %s CSVs under %s", dtype, dtype_root)
    next
  }
  files <- files[!grepl("_reference_", basename(files), ignore.case = TRUE)]
  if (!is.na(target_dataset)) {
    files <- Filter(function(f) dataset_id(f, dtype) == target_dataset, files)
  }
  if (!length(files)) {
    say("No matching datasets for %s", dtype)
    next
  }
  for (expr_file in sort(files)) {
    process_file(expr_file, dtype)
  }
}

#!/usr/bin/env Rscript

# Perform ComBat-based batch correction per study/platform/sample_type using the
# combined per-study CSV layout:
#   <input_root>/<dtype>/<dataset_id>_<dtype>.csv
# The corrected matrices are written to:
#   <output_root>/<dtype>/<dataset_id>_<dtype>.csv
# Optional iso_frac matrices are regenerated from iso_log.

suppressPackageStartupMessages({
  library(data.table)
})
if (!requireNamespace("sva", quietly = TRUE)) {
  stop("Package 'sva' is required. Install via install.packages('sva').")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 5) {
  stop("Usage: 03_combat_batch_correct.R <input_root> <output_root> [comma_separated_dtypes] [target_dataset] [batch_annotation_dir]")
}

input_root <- args[[1]]
output_root <- args[[2]]
dtype_arg <- if (length(args) >= 3) args[[3]] else "gene,iso_log"
target_dataset <- if (length(args) >= 4) args[[4]] else NA_character_
batch_dir <- if (length(args) == 5) args[[5]] else "02_proteomics/data/batch_annotation"

dtypes <- unique(trimws(strsplit(dtype_arg, ",")[[1]]))
dtypes <- dtypes[nzchar(dtypes)]
if (!length(dtypes)) stop("No dtypes specified.")

# Mapping for iso_frac derivation
map_file <- Sys.getenv("PROT_ENST_ENSG_MAP", "02_proteomics/data/raw/ENST-ENSG_mapping.csv")
load_mapping <- function(path) {
  if (!file.exists(path)) stop("[combat] ENST-ENSG mapping not found: ", path)
  dt <- fread(path, header = FALSE, col.names = c("transcript_id", "gene_id"))
  dt[, transcript_id := trimws(transcript_id)]
  dt[, gene_id := trimws(gene_id)]
  dt <- dt[nzchar(transcript_id) & nzchar(gene_id)]
  dt
}
mapping_dt <- load_mapping(map_file)

iso_frac_from_iso_log <- function(mat, mapping_dt) {
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) return(NULL)
  tx <- rownames(mat)
  gene_map <- mapping_dt$gene_id[match(tx, mapping_dt$transcript_id)]
  keep <- !is.na(gene_map)
  if (!any(keep)) return(NULL)
  mat <- mat[keep, , drop = FALSE]
  gene_map <- gene_map[keep]
  res <- mat
  genes <- unique(gene_map)
  for (g in genes) {
    idx <- which(gene_map == g)
    sub <- mat[idx, , drop = FALSE]
    denom <- colSums(sub, na.rm = TRUE)
    frac <- sweep(sub, 2, denom, "/")
    frac[, denom == 0] <- NA_real_
    res[idx, ] <- frac
  }
  res
}

if (!dir.exists(input_root)) stop("Input root not found: ", input_root)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

say <- function(fmt, ...) cat(sprintf(paste0("[combat] ", fmt, "\n"), ...))

canonical_dataset_name <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(x)
  if (grepl("_reference$", x, ignore.case = TRUE)) return(x)
  sub("_(normal|tumor)$", "", x, ignore.case = TRUE)
}

is_reference_dataset <- function(x) {
  !is.null(x) && isTRUE(grepl("_reference$", x, ignore.case = TRUE))
}

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
  if (is_reference_dataset(dataset)) return("reference")
  "tumor"
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
  base <- sub(paste0("_", dtype, "\\.csv$"), "", basename(path), ignore.case = TRUE)
  canonical_dataset_name(base)
}

reference_candidates <- function(file, dtype) {
  base_raw <- sub(paste0("_", dtype, "\\.csv$"), "", basename(file), ignore.case = TRUE)
  base <- canonical_dataset_name(base_raw)
  dirn <- dirname(file)
  cand <- c(
    file.path(dirn, paste0(base, "_reference_", dtype, ".csv")),
    file.path(dirn, paste0(base_raw, "_reference_", dtype, ".csv"))
  )
  unique(cand)
}

target_match_values <- if (!is.na(target_dataset) && nzchar(target_dataset)) {
  unique(na.omit(c(target_dataset, canonical_dataset_name(target_dataset))))
} else {
  character()
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

  say("Starting ComBat batch correction prep for %s (%s)", dataset, dtype)
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

  # batches from annotation
  study_key <- sub("_(normal|tumor|reference)$", "", dataset, ignore.case = TRUE)
  case_batches <- load_case_batches(study_key)
  study_case_ids <- toupper(colnames(mat))
  study_batches <- rep(NA_character_, length(study_case_ids))
  if (!is.null(case_batches)) {
    dt_sub <- case_batches
    if (!is.na(dataset_type) && "sample_type" %in% names(case_batches)) {
      sel <- dt_sub[sample_type == dataset_type]
      if (nrow(sel)) dt_sub <- sel
    }
    idx <- match(study_case_ids, dt_sub$case_id)
    study_batches <- dt_sub$batch_id[idx]
  } else {
    say("Dataset %s has no batch annotation file for %s", dataset, study_key)
  }
  ref_batches <- if (has_reference) reference_batch_labels(colnames(ref_mat)) else character()
  batch_vec <- c(study_batches, ref_batches)
  # fallback batch labels for missing ones (each its own batch)
  if (all(is.na(batch_vec))) {
    batch_vec <- paste0("batch_", seq_along(sample_names))
  } else {
    missing <- is.na(batch_vec) | !nzchar(batch_vec)
    if (any(missing)) batch_vec[missing] <- paste0("batch_", seq_along(sample_names))[missing]
  }

  unique_batches <- unique(batch_vec)
  if (length(unique_batches) < 2) {
    say("Dataset %s (%s) has <2 batches; copying input", dataset, dtype)
    corrected <- mat
  } else {
    say("Dataset %s (%s) batches: %s", dataset, dtype, paste(unique_batches, collapse = ","))
    # simple per-feature median imputation for NAs
    full_mat_imp <- apply(full_mat, 1, function(row) {
      if (all(!is.finite(row))) return(rep(0, length(row)))
      m <- safe_median(row)
      row[!is.finite(row)] <- m
      row
    })
    full_mat_imp <- t(full_mat_imp)
    combat_res <- sva::ComBat(dat = full_mat_imp, 
                              batch = batch_vec, 
                              par.prior = TRUE, 
                              prior.plots = FALSE)
    corrected_full <- combat_res[, !is_reference, drop = FALSE]
    rownames(corrected_full) <- rownames(full_mat_imp)
    if (all(!is.finite(corrected_full))) {
      say("Dataset %s (%s) ComBat returned all non-finite values → keeping uncorrected input",
          dataset, dtype)
      corrected <- mat
    } else {
      corrected <- corrected_full
      say("Dataset %s (%s) batch correction completed via ComBat", dataset, dtype)
    }
  }

  out_dir <- file.path(output_root, dtype)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_name <- sprintf("%s_%s.csv", dataset, dtype)
  out_file <- file.path(out_dir, out_name)
  out_dt <- data.table(feature = rownames(corrected))
  out_dt <- cbind(out_dt, as.data.table(corrected))
  fwrite(out_dt, out_file)

  if (dtype == "iso_log") {
    eps <- 1e-12
    lin <- pmax(2^corrected - 1, 0)  # invert log2(+1)
    frac_mat <- iso_frac_from_iso_log(lin, mapping_dt)
    if (is.null(frac_mat)) {
      warning(sprintf("[combat] Could not derive iso_frac for %s (no transcript->gene mapping)", dataset))
    } else {
      frac_mat <- pmin(pmax(frac_mat, eps), 1 - eps)
      frac_dt <- data.table(feature = rownames(frac_mat))
      frac_dt <- cbind(frac_dt, as.data.table(frac_mat))
      frac_dir <- file.path(output_root, "iso_frac")
      dir.create(frac_dir, recursive = TRUE, showWarnings = FALSE)
      frac_file <- file.path(frac_dir, iso_frac_name(out_name))
      fwrite(frac_dt, frac_file)
      say("Derived iso_frac for %s", dataset)
    }
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
  if (length(target_match_values)) {
    files <- Filter(function(f) dataset_id(f, dtype) %in% target_match_values, files)
  }
  if (!length(files)) {
    say("No matching datasets for %s", dtype)
    next
  }
  for (expr_file in sort(files)) {
    process_file(expr_file, dtype)
  }
}

#!/usr/bin/env Rscript

# ============================================================
# Proteomics univariate CoxPH with batch stratification
# ------------------------------------------------------------
# Mirrors the transcriptomics 03a script but:
#   * works on proteomics expression matrices
#   * optionally stratifies by batch (if batch annotation exists)
#   * collapses duplicate case_ids by averaging replicates
# Usage:
#   Rscript 01_univariate_coxph.R \
#     <expr_matrix.csv> \
#     <clinical_manifest.csv> \
#     <covariates.yaml> \
#     <batch_annotation_dir> \
#     <out_results.csv> \
#     <out_summary.txt>
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(parallel)
  library(yaml)
})

say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }
`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================================================
# Arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  die("Usage: 01_univariate_coxph.R <expr.csv> <clinical.csv> <covariates.yaml> <batch_dir> <out_results.csv> <out_summary.txt>")
}

expr_in   <- args[[1]]
clin_in   <- args[[2]]
cov_yaml  <- args[[3]]
batch_dir <- args[[4]]
out_res   <- args[[5]]
out_sum   <- args[[6]]

dir.create(dirname(out_res), recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helpers
# ============================================================
infer_study_name <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  stem <- sub("_(gene|iso_log|iso_frac)$", "", stem, ignore.case = TRUE)
  stem <- sub("_(normal|tumor)$", "", stem, ignore.case = TRUE)
  stem
}

read_expr_csv <- function(file) {
  dt <- fread(file)
  if (ncol(dt) < 2) return(NULL)
  if (names(dt)[1] == "") setnames(dt, 1, "feature")
  feats <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- feats
  storage.mode(mat) <- "double"
  mat
}

combine_matrices <- function(base, add) {
  if (is.null(base)) return(add)
  if (is.null(add)) return(base)
  feats <- union(rownames(base), rownames(add))
  base_ext <- matrix(NA_real_, length(feats), ncol(base))
  rownames(base_ext) <- feats
  colnames(base_ext) <- colnames(base)
  base_ext[match(rownames(base), feats), ] <- base
  add_ext <- matrix(NA_real_, length(feats), ncol(add))
  rownames(add_ext) <- feats
  colnames(add_ext) <- colnames(add)
  add_ext[match(rownames(add), feats), ] <- add
  cbind(base_ext, add_ext)
}

load_expr_from_file <- function(path) {
  mat <- read_expr_csv(path)
  if (is.null(mat)) die("Failed to read expression matrix: %s", path)
  list(
    expr = mat,
    study = infer_study_name(path),
    platform = NA_character_,
    data_type = sub("^.*_(gene|iso_log|iso_frac)\\.csv$", "\\1", basename(path)),
    cohort = infer_study_name(path)
  )
}

load_expr_from_dir <- function(dir_path) {
  norm_dir <- normalizePath(dir_path, winslash = "/", mustWork = TRUE)
  dtype <- basename(dirname(norm_dir))
  if (!dtype %in% c("gene","iso_log","iso_frac")) {
    die("Could not infer data type from path: %s", dir_path)
  }
  files <- list.files(norm_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (!length(files)) die("No CSV files found under directory: %s", dir_path)
  target_files <- files[basename(files) == paste0(dtype, ".csv")]
  if (!length(target_files)) die("No %s CSV files found under %s", dtype, dir_path)

  rel_paths <- sub(paste0("^", norm_dir, "/?"), "", normalizePath(target_files, winslash = "/", mustWork = TRUE))
  rel_parts <- strsplit(rel_paths, "/")
  platform_vals <- unique(vapply(rel_parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1)))
  platform_vals <- platform_vals[!is.na(platform_vals)]
  platform <- if (length(platform_vals) == 1) platform_vals else NA_character_
  sample_type_vals <- unique(vapply(rel_parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1)))
  sample_type_vals <- sample_type_vals[!is.na(sample_type_vals)]
  sample_type <- if (length(sample_type_vals) == 1) sample_type_vals else NA_character_
  study <- basename(norm_dir)

  expr_mat <- NULL
  for (f in sort(target_files)) {
    mat <- read_expr_csv(f)
    if (!is.null(mat)) expr_mat <- combine_matrices(expr_mat, mat)
  }
  if (is.null(expr_mat)) die("Failed to assemble expression matrix from %s", dir_path)
  cohort_parts <- c(study, platform, sample_type)
  cohort_label <- paste(cohort_parts[!is.na(cohort_parts)], collapse = "/")

  list(
    expr = expr_mat,
    study = study,
    platform = platform,
    sample_type = sample_type,
    data_type = dtype,
    cohort = cohort_label
  )
}

find_batch_file <- function(study, batch_dir) {
  if (!dir.exists(batch_dir)) return(NA_character_)
  files <- list.files(batch_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) return(NA_character_)
  bases <- tools::file_path_sans_ext(basename(files))
  idx <- which(tolower(bases) == tolower(study))
  if (!length(idx)) return(NA_character_)
  files[idx[1]]
}

pick_batch_column <- function(dt) {
  cols <- intersect(c("Folder_name","Batch","BATCH","batch_id",
                      "folder_name","folder","run","Run",
                      "specimen_run","batchname"),
                    names(dt))
  if (length(cols)) return(cols[1])
  NA_character_
}

collapse_duplicates <- function(mat) {
  ids <- colnames(mat)
  dups <- ids[duplicated(ids)]
  if (!length(dups)) return(mat)
  say("[INFO] Averaging duplicate case_ids: %s", paste(unique(dups), collapse = ", "))
  uniq <- unique(ids)
  agg <- vapply(uniq, function(id) {
    cols <- mat[, ids == id, drop = FALSE]
    if (ncol(cols) == 1) cols[, 1] else rowMeans(cols, na.rm = TRUE)
  }, numeric(nrow(mat)))
  rownames(agg) <- rownames(mat)
  agg
}

# ============================================================
# Load data
# ============================================================
expr_meta <- if (dir.exists(expr_in)) load_expr_from_dir(expr_in) else load_expr_from_file(expr_in)
expr <- collapse_duplicates(expr_meta$expr)
mani <- fread(clin_in)
cfg  <- yaml::read_yaml(cov_yaml)

study <- as.character(expr_meta$study)
platform <- expr_meta$platform
data_type <- expr_meta$data_type
cohort_label <- as.character(expr_meta$cohort %||% study)
cancer_map <- cfg$cancers
cancer <- if (!is.null(cancer_map) && study %in% names(cancer_map)) as.character(cancer_map[[study]]) else study
run_id <- sprintf("[prot_cox][%s][%s]", cohort_label, data_type)

say("=== Proteomics CoxPH ===")
say("%s Expression source: %s", run_id, expr_in)
say("%s Clinical file: %s", run_id, clin_in)

# ============================================================
# Clean clinical data
# ============================================================
mani <- mani[complete.cases(mani[, .(case_id, OS_time, OS_event)]), ]
mani <- mani[!duplicated(case_id)]

# Stage ordering if provided
if (!is.null(cfg$stage_levels) && "stage" %in% names(mani)) {
  mani[, stage := factor(stage, levels = cfg$stage_levels, ordered = TRUE)]
}

# Align samples
common <- intersect(colnames(expr), mani$case_id)
if (length(common) < 20) {
  die("%s Too few overlapping samples (n=%d)", run_id, length(common))
}
expr <- expr[, common, drop = FALSE]
mani <- mani[match(common, mani$case_id)]

# ============================================================
# Attach batch annotations
# ============================================================
batch_file <- find_batch_file(study, batch_dir)
if (is.na(batch_file) && !is.null(platform) && !is.na(platform)) {
  batch_file <- find_batch_file(paste(study, platform, sep = "_"), batch_dir)
}
if (!is.na(batch_file)) {
  say("%s Using batch annotation: %s", run_id, batch_file)
  batch_dt <- fread(batch_file)
  colnames(batch_dt) <- tolower(colnames(batch_dt))
  batch_col <- pick_batch_column(batch_dt)
  if (!is.na(batch_col) && "case_id" %in% names(batch_dt)) {
    mani <- merge(
      mani,
      batch_dt[, .(case_id, batch_id = get(batch_col))],
      by = "case_id",
      all.x = TRUE,
      suffixes = c("", "_batch")
    )
  } else {
    say("%s [WARN] batch column not found in %s", run_id, batch_file)
  }
} else {
  say("%s [INFO] No batch annotation found for %s", run_id, study)
}

has_batch <- "batch_id" %in% names(mani) &&
  sum(!is.na(mani$batch_id)) >= 2 &&
  length(unique(na.omit(mani$batch_id))) > 1
if (has_batch) {
  mani[, batch_id := factor(batch_id)]
  say("%s Stratifying by batch (%d levels)", run_id, length(levels(mani$batch_id)))
}

y <- with(mani, Surv(OS_time, OS_event))

# ============================================================
# Covariates
# ============================================================
covariates <- cfg$baseline_covariates
if (!is.null(cfg$conditional_covariates)) {
  for (cond in cfg$conditional_covariates) {
    nm <- cond$name
    thr <- cond$coverage_threshold %||% 0
    if (nm %in% names(mani) && mean(!is.na(mani[[nm]])) >= thr) {
      covariates <- c(covariates, nm)
    }
  }
}
covariates <- unique(na.omit(covariates))
Xcov <- mani[, covariates, with = FALSE]
Xcov <- Xcov[, lapply(.SD, function(x) if (length(unique(na.omit(x))) > 1) x else NULL)]
say("%s Covariates used: %s", run_id,
    if (ncol(Xcov)) paste(names(Xcov), collapse = ", ") else "<none>")

# ============================================================
# Build formula
# ============================================================
rhs_terms <- c("expr")
if (ncol(Xcov)) rhs_terms <- c(rhs_terms, names(Xcov))
if (has_batch) rhs_terms <- c(rhs_terms, "strata(batch_id)")
fml <- as.formula(paste("y ~", paste(rhs_terms, collapse = " + ")))
say("%s Cox formula: %s", run_id, deparse(fml))

# ============================================================
# Run univariate CoxPH
# ============================================================
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
idx <- seq_len(nrow(expr))

res_list <- mclapply(idx, function(i) {
  tryCatch({
    vals <- as.numeric(expr[i, ])
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) return(data.frame())
    vals_z <- as.numeric(scale(vals))
    df <- data.frame(expr = vals_z, Xcov)
    if (has_batch) df$batch_id <- mani$batch_id
    df$y <- y

    fit <- suppressWarnings(coxph(fml, data = df, ties = "efron"))
    if (inherits(fit, "try-error") || isFALSE(fit$converged)) return(data.frame())
    s <- summary(fit)
    if (!"expr" %in% rownames(s$coef)) return(data.frame())

    z_expr        <- s$coef["expr", "z"]
    wald_expr     <- z_expr^2
    delta_LL_expr <- wald_expr / 2

    data.frame(
      feature       = rownames(expr)[i],
      beta          = s$coef["expr","coef"],
      HR            = exp(s$coef["expr","coef"]),
      z             = z_expr,
      p             = s$coef["expr","Pr(>|z|)"],
      se            = s$coef["expr","se(coef)"],
      delta_LL_expr = delta_LL_expr
    )
  }, error = function(e) data.frame())
}, mc.cores = ncores)

res_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, res_list)
if (!length(res_list)) die("%s No valid CoxPH models fitted.", run_id)

res <- rbindlist(res_list, fill = TRUE)
res <- res[is.finite(p)]
res[, FDR := p.adjust(p, "BH")]
res[, study := cohort_label]
res[, data_type := data_type]

out_full <- sub("\\.csv$", "_full.csv", out_res)
fwrite(res, out_full)
fwrite(res[FDR < 0.05 & is.finite(HR) & HR > 0], out_res)

sink(out_sum)
cat("=== Proteomics Cox Summary ===\n")
cat("Study:", cancer, "\n")
cat("Datatype:", data_type, "\n")
cat("Samples:", length(common), "\n")
cat("Features tested:", nrow(res), "\n")
cat("Significant (FDR<0.05):", sum(res$FDR < 0.05, na.rm = TRUE), "\n")
cat("Batch strata:", if (has_batch) length(levels(mani$batch_id)) else "none", "\n")
sink()

say("%s Done. Results: %s | %s", run_id, out_full, out_res)

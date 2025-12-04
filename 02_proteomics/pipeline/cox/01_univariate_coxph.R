#!/usr/bin/env Rscript

# ============================================================
# Proteomics univariate CoxPH
# ------------------------------------------------------------
# Mirrors the transcriptomics 03a script but:
#   * works on proteomics expression matrices
#   * collapses duplicate case_ids by averaging replicates
# Usage:
#   Rscript 01_univariate_coxph.R \
#     <expr_matrix.csv> \
#     <clinical_manifest.csv> \
#     <covariates.yaml> \
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
if (length(args) != 5) {
  die("Usage: 01_univariate_coxph.R <expr.csv> <clinical.csv> <covariates.yaml> <out_results.csv> <out_summary.txt>")
}

expr_in   <- args[[1]]
clin_in   <- args[[2]]
cov_yaml  <- args[[3]]
out_res   <- args[[4]]
out_sum   <- args[[5]]

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

infer_sample_type <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  if (grepl("_reference$", stem, ignore.case = TRUE)) return("reference")
  "tumor"
}

parse_thresholds <- function(env_var, default_vals) {
  raw <- Sys.getenv(env_var, "")
  if (nzchar(raw)) {
    vals <- suppressWarnings(as.numeric(strsplit(raw, ",")[[1]]))
    vals <- vals[is.finite(vals) & vals > 0]
    vals <- unique(vals)
    if (length(vals)) return(sort(vals))
  }
  sort(unique(default_vals))
}

format_threshold <- function(x) {
  sub("0+$", "", sub("(\\.\\d*?)0+$", "\\1", sprintf("%.3f", x)))
}

FDR_THRESHOLDS <- parse_thresholds("COX_UNIV_FDR_THRESHOLDS", c(0.05, 0.10))
P_THRESHOLDS <- parse_thresholds("COX_UNIV_P_THRESHOLDS", c(0.05, 0.10))

write_sig_summary <- function(count_fun, thresholds, label) {
  if (!length(thresholds)) return()
  for (thr in thresholds) {
    count <- count_fun(thr)
    cat(sprintf("Significant (%s<%s): %d\n", label, format_threshold(thr), count))
  }
}

build_cohort_label <- function(study, platform, sample_type) {
  parts <- c(study)
  if (!is.null(platform) && !is.na(platform) && nzchar(platform)) {
    parts <- c(parts, platform)
  }
  if (!is.null(sample_type) && !is.na(sample_type) && nzchar(sample_type) &&
      sample_type != "tumor") {
    parts <- c(parts, sample_type)
  }
  paste(parts, collapse = "/")
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
  sample_type <- infer_sample_type(path)
  cohort_label <- build_cohort_label(infer_study_name(path), NA_character_, sample_type)
  list(
    expr = mat,
    study = infer_study_name(path),
    platform = NA_character_,
    sample_type = sample_type,
    data_type = sub("^.*_(gene|iso_log|iso_frac)\\.csv$", "\\1", basename(path)),
    cohort = cohort_label
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
  sample_type <- if (length(sample_type_vals) == 1) {
    val <- tolower(sample_type_vals)
    if (val == "reference") "reference" else "tumor"
  } else {
    NA_character_
  }
  study <- basename(norm_dir)

  expr_mat <- NULL
  for (f in sort(target_files)) {
    mat <- read_expr_csv(f)
    if (!is.null(mat)) expr_mat <- combine_matrices(expr_mat, mat)
  }
  if (is.null(expr_mat)) die("Failed to assemble expression matrix from %s", dir_path)
  cohort_label <- build_cohort_label(study, platform, sample_type)

  list(
    expr = expr_mat,
    study = study,
    platform = platform,
    sample_type = sample_type,
    data_type = dtype,
    cohort = cohort_label
  )
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

# Quantile-normalize a numeric vector to standard normal scores (per feature)
quantile_normalize_feature <- function(x) {
  res <- rep(NA_real_, length(x))
  idx <- which(is.finite(x))
  if (!length(idx)) return(res)
  vals <- x[idx]
  ranks <- rank(vals, ties.method = "average")
  n <- length(ranks)
  # map ranks to (0,1) then to z-scores
  probs <- (ranks - 0.5) / n
  eps <- 1e-6
  probs <- pmin(pmax(probs, eps), 1 - eps)
  res[idx] <- qnorm(probs)
  res
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
sample_type <- if (!is.null(expr_meta$sample_type) && !is.na(expr_meta$sample_type)) expr_meta$sample_type else "tumor"
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
if (!length(common)) {
  die("%s No overlapping samples between expression and clinical data", run_id)
}
expr <- expr[, common, drop = FALSE]
mani <- mani[match(common, mani$case_id)]

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
    finite_vals <- vals[is.finite(vals)]
    if (length(finite_vals) < 2 || sd(finite_vals, na.rm = TRUE) == 0) return(data.frame())
    vals_qn <- quantile_normalize_feature(vals)
    if (!any(is.finite(vals_qn))) return(data.frame())
    df <- data.frame(expr = vals_qn, Xcov)
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
out_full <- sub("\\.csv$", "_full.csv", out_res)

if (!length(res_list)) {
  say("%s No valid CoxPH models fitted → writing empty outputs", run_id)
  empty <- data.table(
    feature = character(),
    beta = numeric(),
    HR = numeric(),
    z = numeric(),
    p = numeric(),
    se = numeric(),
    delta_LL_expr = numeric(),
    FDR = numeric(),
    study = character(),
    data_type = character()
  )
  fwrite(empty, out_full)
  fwrite(empty, out_res)
  sink(out_sum)
  cat("=== Proteomics Cox Summary ===\n")
  cat("Study:", cancer, "\n")
  cat("Datatype:", data_type, "\n")
  cat("Samples:", length(common), "\n")
  cat("Features tested:", 0, "\n")
  write_sig_summary(function(thr) 0, FDR_THRESHOLDS, "FDR")
  write_sig_summary(function(thr) 0, P_THRESHOLDS, "p")
  sink()
  quit(status = 0)
}

res <- rbindlist(res_list, fill = TRUE)
res <- res[is.finite(p)]
res[, FDR := p.adjust(p, "BH")]
res[, study := cohort_label]
res[, data_type := data_type]
res[, sample_type := sample_type]
fwrite(res, out_full)
top_features <- res[FDR < 0.05 & is.finite(HR) & HR > 0 & !is.na(HR)]
if (!nrow(top_features)) {
  fallback <- res[is.finite(HR) & !is.na(HR)]
  fallback <- fallback[order(FDR, na.last = TRUE)]
  fallback <- fallback[is.finite(FDR)]
  top_features <- fallback[seq_len(min(100L, nrow(fallback)))]
}
fwrite(top_features, out_res)

sink(out_sum)
cat("=== Proteomics Cox Summary ===\n")
cat("Study:", cancer, "\n")
cat("Datatype:", data_type, "\n")
cat("Samples:", length(common), "\n")
cat("Features tested:", nrow(res), "\n")
write_sig_summary(function(thr) sum(res$FDR < thr, na.rm = TRUE), FDR_THRESHOLDS, "FDR")
write_sig_summary(function(thr) sum(res$p < thr, na.rm = TRUE), P_THRESHOLDS, "p")
sink()

say("%s Done. Results: %s | %s", run_id, out_full, out_res)

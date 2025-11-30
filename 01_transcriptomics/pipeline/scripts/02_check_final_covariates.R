#!/usr/bin/env Rscript
# ============================================================
# check_manifest_covariates.R (optimized + progress bar)
# ------------------------------------------------------------
# Prints covariate completeness (age/stage/event) across cancers
# + generates a features-kept summary per data type.
# Uses fast line counting and shows a progress bar during counting.
# ============================================================

suppressPackageStartupMessages(library(data.table))
say <- function(...) message(sprintf(...))

# ANSI color helpers
col_red    <- function(x) paste0("\033[31m", x, "\033[0m")
col_yellow <- function(x) paste0("\033[33m", x, "\033[0m")
col_green  <- function(x) paste0("\033[32m", x, "\033[0m")

# ============================================================
# Paths
# ============================================================
base_dir <- "01_transcriptomics/out/02_norm"
data_dir <- "01_transcriptomics/data"
if (!dir.exists(base_dir)) stop(sprintf("Manifest directory not found: %s", base_dir))

say("[INFO] Scanning manifests in: %s", base_dir)
mani_files <- list.files(
  path = base_dir,
  pattern = "^TCGA_.*\\.sample_manifest\\.csv$",
  full.names = TRUE
)
if (length(mani_files) == 0) {
  cat("[WARN] No manifest files found.\n")
  quit(status = 0)
}

# ============================================================
# Helper: fast line count (no file load)
# ============================================================
fast_linecount <- function(file) {
  if (!file.exists(file)) return(NA_integer_)
  lines <- count.fields(file, sep = ",", blank.lines.skip = TRUE)
  n <- length(lines) - 1L
  if (n < 0) n <- 0L
  n
}

cache_count <- new.env(parent = emptyenv())
get_feature_count <- function(path) {
  if (is.null(cache_count[[path]])) {
    cache_count[[path]] <- fast_linecount(path)
  }
  cache_count[[path]]
}

# ============================================================
# (1) Manifest completeness summary
# ============================================================
mani_summary <- rbindlist(lapply(mani_files, function(f) {
  cancer <- sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", basename(f))
  mani <- fread(f, showProgress = FALSE, nThread = 1L)
  mani[mani == "" | mani == " " | mani == "NA"] <- NA

  n_total  <- nrow(mani)
  n_age    <- if ("age" %in% names(mani)) sum(!is.na(mani$age) & trimws(mani$age) != "") else 0
  n_events <- if ("OS_event" %in% names(mani)) sum(mani$OS_event == 1, na.rm = TRUE) else 0

  n_stage <- n_I <- n_II <- n_III <- n_IV <- 0
  p_I <- p_II <- p_III <- p_IV <- 0
  if ("stage" %in% names(mani)) {
    mani$stage <- trimws(as.character(mani$stage))
    mani$stage <- gsub("^Stage ", "", mani$stage)
    mani$stage <- factor(mani$stage, levels = c("I", "II", "III", "IV"), ordered = TRUE)
    n_stage <- sum(!is.na(mani$stage))
    st_tab <- table(mani$stage)
    n_I <- as.integer(st_tab["I"]); if (is.na(n_I)) n_I <- 0
    n_II <- as.integer(st_tab["II"]); if (is.na(n_II)) n_II <- 0
    n_III <- as.integer(st_tab["III"]); if (is.na(n_III)) n_III <- 0
    n_IV <- as.integer(st_tab["IV"]); if (is.na(n_IV)) n_IV <- 0
    p_I <- round(100 * n_I / n_total, 1)
    p_II <- round(100 * n_II / n_total, 1)
    p_III <- round(100 * n_III / n_total, 1)
    p_IV <- round(100 * n_IV / n_total, 1)
  }

  data.table(
    cancer,
    n_total,
    n_age,
    n_stage,
    n_events,
    n_I, n_II, n_III, n_IV,
    p_I, p_II, p_III, p_IV
  )
}), fill = TRUE)

mani_summary <- mani_summary[, .(
  n_total  = max(n_total, na.rm = TRUE),
  n_age    = max(n_age, na.rm = TRUE),
  n_stage  = max(n_stage, na.rm = TRUE),
  n_events = max(n_events, na.rm = TRUE),
  n_I   = max(n_I, na.rm = TRUE),
  n_II  = max(n_II, na.rm = TRUE),
  n_III = max(n_III, na.rm = TRUE),
  n_IV  = max(n_IV, na.rm = TRUE),
  p_I   = max(p_I, na.rm = TRUE),
  p_II  = max(p_II, na.rm = TRUE),
  p_III = max(p_III, na.rm = TRUE),
  p_IV  = max(p_IV, na.rm = TRUE)
), by = cancer][order(cancer)]

# ============================================================
# (2) Feature summary table with progress bar
# ============================================================
say("[INFO] Counting features (kept / total) using fast linecount...")

pb <- txtProgressBar(min = 0, max = length(mani_files), style = 3)
feature_summary <- rbindlist(lapply(seq_along(mani_files), function(i) {
  f <- mani_files[i]
  setTxtProgressBar(pb, i)

  cancer <- sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", basename(f))
  datatype <- sub("^TCGA_[A-Z0-9]+_(.*?)\\.sample_manifest\\.csv$", "\\1", basename(f))

  expr_norm <- sub("\\.sample_manifest\\.csv$", ".normalized.csv", f)
  expr_in_guess <- list.files(
    data_dir, pattern = paste0("RNA_", cancer, "_", datatype, "\\.csv$"),
    recursive = TRUE, full.names = TRUE
  )[1]

  features_kept  <- get_feature_count(expr_norm)
  features_total <- get_feature_count(expr_in_guess)
  frac_kept <- if (!is.na(features_kept) && !is.na(features_total))
    round(100 * features_kept / features_total, 1) else NA_real_

  data.table(cancer, datatype, features_kept, features_total, frac_kept)
}), fill = TRUE)
close(pb)

# ============================================================
# --- Print summaries
# ============================================================
cat("\n📊 Covariate completeness summary (per cancer type)\n")
cat(sprintf("Found %d manifest files.\n\n", nrow(mani_summary)))
cat(sprintf("%-6s | %7s | %12s | %12s | %14s\n",
            "Cancer", "Samples", "Age valid", "Stage valid", "Events"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(mani_summary))) {
  row <- mani_summary[i]
  n_total <- row$n_total
  age_pct <- if (n_total > 0) round(100 * row$n_age / n_total, 1) else 0
  stage_pct <- if (n_total > 0) round(100 * row$n_stage / n_total, 1) else 0
  event_pct <- if (n_total > 0) round(100 * row$n_events / n_total, 1) else 0

  age_txt   <- sprintf("%4d (%5.1f%%)", row$n_age, age_pct)
  stage_txt <- sprintf("%4d (%5.1f%%)", row$n_stage, stage_pct)
  event_txt <- sprintf("%4d (%5.1f%%)", row$n_events, event_pct)

  if (age_pct < 70) age_txt <- col_red(age_txt)
  else if (age_pct < 90) age_txt <- col_yellow(age_txt)
  else age_txt <- col_green(age_txt)

  if (stage_pct < 70) stage_txt <- col_red(stage_txt)
  else if (stage_pct < 90) stage_txt <- col_yellow(stage_txt)
  else stage_txt <- col_green(stage_txt)

  if (event_pct < 20) event_txt <- col_red(event_txt)
  else if (event_pct < 40) event_txt <- col_yellow(event_txt)
  else event_txt <- col_green(event_txt)

  cat(sprintf("%-6s | %7d | %12s | %12s | %14s\n",
              row$cancer, n_total, age_txt, stage_txt, event_txt))
}
cat(strrep("-", 70), "\n")

cat("\n📈 Feature retention summary (per cancer × data type)\n")
cat(sprintf("%-6s | %-9s | %10s | %10s | %8s\n",
            "Cancer", "Datatype", "Kept", "Total", "Kept %"))
cat(strrep("-", 55), "\n")
for (i in seq_len(nrow(feature_summary))) {
  row <- feature_summary[i]
  pct_col <- if (is.na(row$frac_kept)) "" else
    if (row$frac_kept < 10) col_red(sprintf("%5.1f", row$frac_kept))
    else if (row$frac_kept < 30) col_yellow(sprintf("%5.1f", row$frac_kept))
    else col_green(sprintf("%5.1f", row$frac_kept))
  cat(sprintf("%-6s | %-9s | %10d | %10d | %8s\n",
              row$cancer, row$datatype,
              row$features_kept, row$features_total, pct_col))
}
cat(strrep("-", 55), "\n")

# ============================================================
# Save both tables
# ============================================================
out_cov <- file.path(base_dir, "manifest_covariate_counts.csv")
out_feat <- file.path(base_dir, "manifest_feature_counts.csv")
fwrite(mani_summary, out_cov)
fwrite(feature_summary, out_feat)
cat(sprintf("[SAVED] %s\n", out_cov))
cat(sprintf("[SAVED] %s\n\n", out_feat))

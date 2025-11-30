#!/usr/bin/env Rscript

# 04_collect_univariate_results.R
# Collate all univariate Cox results (per cancer × datatype), coerce numerics safely,
# compute FDR (within-panel & global), and emit master + summary tables.

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

say <- function(...) message(sprintf(...))

IN_DIR  <- "01_transcriptomics/out/03_univariate_coxph"
OUT_DIR <- "01_transcriptomics/out/04_univariate_collect"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# discover result files
files <- list.files(IN_DIR, pattern = "\\.cox_results\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No cox_results files found in: ", IN_DIR)
say("Found %d result files.", length(files))

# parse cancer & dtype from file name: TCGA_<CANCER>_<DTYPE>.cox_results.csv
parse_meta <- function(fp) {
  fn <- basename(fp)
  m <- str_match(fn, "^TCGA_([A-Z0-9]+)_([a-z_]+)\\.cox_results\\.csv$")
  if (is.na(m[1,1])) return(list(cancer=NA_character_, dtype=NA_character_))
  list(cancer = m[1,2], dtype = m[1,3])
}

# read and tag each file
all_dt <- rbindlist(lapply(files, function(f) {
  meta <- parse_meta(f)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt) || !nrow(dt)) return(NULL)

  # Ensure expected columns exist (fill with NA if missing)
  need <- c("feature","n","beta","HR","z","p","cindex","se")
  miss <- setdiff(need, names(dt))
  for (m in miss) dt[[m]] <- NA_real_

  dt[, cancer := meta$cancer][, dtype := meta$dtype][]
}), fill = TRUE)

if (!nrow(all_dt)) stop("After reading, no rows present. Check input files.")

# keep rows with a feature id
all_dt <- all_dt[!is.na(feature)]
setcolorder(all_dt, c("cancer","dtype","feature","n","beta","HR","z","p","cindex","se"))

# --- numeric coercion & sanitization ---
num_cols <- c("n","beta","HR","z","p","cindex","se")
for (cc in num_cols) {
  if (!cc %in% names(all_dt)) all_dt[[cc]] <- NA_real_
  suppressWarnings({
    all_dt[[cc]] <- as.numeric(gsub("[ ,]", "", as.character(all_dt[[cc]])))
  })
}

# Guard against invalid HR for log()
bad_before <- sum(!is.finite(all_dt$HR) | all_dt$HR <= 0, na.rm = TRUE)
all_dt[!is.finite(HR) | HR <= 0, HR := NA_real_]
say("Rows with non-finite or non-positive HR (set to NA): %d", bad_before)

# FDR (within each cancer×dtype) and global
all_dt[, q_within := p.adjust(p, method = "BH"), by = .(cancer, dtype)]
all_dt[, q_global := p.adjust(p, method = "BH")]

# Derived helpers
all_dt[, direction := fifelse(HR > 1, "risk↑", fifelse(HR < 1, "protective↑", "neutral"))]
all_dt[, logHR := ifelse(is.finite(HR) & HR > 0, log(HR), NA_real_)]

# write master table
master_path <- file.path(OUT_DIR, "univariate_master_all.csv")
fwrite(all_dt, master_path)
say("Master table: %s  (%d rows)", master_path, nrow(all_dt))

# index by cancer × dtype
idx <- all_dt[, .(
  n_features = .N,
  n_sig_p05 = sum(p < 0.05, na.rm = TRUE),
  n_sig_q10 = sum(q_within < 0.10, na.rm = TRUE),
  n_sig_q05 = sum(q_within < 0.05, na.rm = TRUE),
  median_cindex = median(cindex, na.rm = TRUE)
), by = .(cancer, dtype)][order(cancer, dtype)]

idx_path <- file.path(OUT_DIR, "index_by_cancer_dtype.csv")
fwrite(idx, idx_path)
say("Index: %s", idx_path)

# filtered “significant” table (tune threshold here)
sig_path <- file.path(OUT_DIR, "significant_q05_within.csv")
fwrite(all_dt[q_within < 0.05], sig_path)
say("Significant (q_within<0.05): %s", sig_path)

# Top-N by raw p within each panel (handy for browsing)
topN <- 50
top_path <- file.path(OUT_DIR, sprintf("top%d_by_p_each_panel.csv", topN))
top_by_panel <- all_dt[order(p)][, head(.SD, topN), by = .(cancer, dtype)]
fwrite(top_by_panel, top_path)
say("Top-N by p (per panel): %s", top_path)

# Optional per-cancer rollups
per_cancer_dir <- file.path(OUT_DIR, "per_cancer_tables")
dir.create(per_cancer_dir, showWarnings = FALSE)
for (cc in sort(unique(all_dt$cancer))) {
  fwrite(all_dt[cancer == cc][order(dtype, p)],
         file.path(per_cancer_dir, sprintf("TCGA_%s_univariate_all.csv", cc)))
  fwrite(all_dt[cancer == cc & q_within < 0.05][order(dtype, p)],
         file.path(per_cancer_dir, sprintf("TCGA_%s_univariate_sig_q05.csv", cc)))
}

say("[DONE]")

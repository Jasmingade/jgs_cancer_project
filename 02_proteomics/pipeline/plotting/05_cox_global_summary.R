#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: 05_cox_global_summary.R <cox_root> <out_dir> [fdr_thresh]")
}

cox_root <- args[[1]]
out_dir  <- args[[2]]
fdr_thr  <- if (length(args) == 3) as.numeric(args[[3]]) else 0.05
if (!is.finite(fdr_thr) || fdr_thr <= 0) fdr_thr <- 0.05

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

say <- function(...) message(sprintf(...))

tx2gene_path <- Sys.getenv("PROT_TX2GENE", "02_proteomics/data/raw/tx2gene.csv")
if (!file.exists(tx2gene_path)) stop("Transcript-to-gene map not found: ", tx2gene_path)
tx2gene <- fread(tx2gene_path)
sel_col <- function(dt, candidates, fallback_idx) {
  for (nm in candidates) if (nm %in% names(dt)) return(nm)
  names(dt)[fallback_idx]
}
iso_col  <- sel_col(tx2gene, c("isoform_id","transcript","tx_id","transcript_id"), 1)
gene_col <- sel_col(tx2gene, c("gene_id","gene","gene_name"), min(2L, ncol(tx2gene)))
tx2gene[, transcript_clean := sub("\\..*$", "", get(iso_col))]
tx2gene[, gene_clean       := sub("\\..*$", "", get(gene_col))]
tx2gene_map <- unique(tx2gene[, .(transcript_clean, gene_clean, gene_name)])
tx2gene_lookup <- setNames(tx2gene_map$gene_clean, tx2gene_map$transcript_clean)
tx2gene_name_lookup <- setNames(tx2gene_map$gene_name, tx2gene_map$transcript_clean)
gene_name_lookup <- tx2gene_map[!duplicated(gene_clean)]
gene_name_lookup <- setNames(gene_name_lookup$gene_name, gene_name_lookup$gene_clean)

map_feature_to_gene <- function(dtype, feature) {
  dtype <- as.character(dtype)
  feature <- as.character(feature)
  res <- rep(NA_character_, length(dtype))
  gene_idx <- !is.na(dtype) & dtype == "gene"
  if (any(gene_idx)) {
    res[gene_idx] <- sub("\\..*$", "", feature[gene_idx])
  }
  iso_idx <- !is.na(dtype) & dtype %in% c("iso_log","iso_frac")
  if (any(iso_idx)) {
    bare <- sub("\\..*$", "", feature[iso_idx])
    res[iso_idx] <- tx2gene_lookup[bare]
  }
  res
}

map_feature_to_gene_name <- function(dtype, feature, gene_ids = NULL) {
  dtype <- as.character(dtype)
  feature <- as.character(feature)
  res <- rep(NA_character_, length(dtype))
  gene_idx <- !is.na(dtype) & dtype == "gene"
  if (any(gene_idx)) {
    ids <- if (is.null(gene_ids)) sub("\\..*$", "", feature[gene_idx]) else gene_ids[gene_idx]
    res[gene_idx] <- gene_name_lookup[ids]
  }
  iso_idx <- !is.na(dtype) & dtype %in% c("iso_log","iso_frac")
  if (any(iso_idx)) {
    bare <- sub("\\..*$", "", feature[iso_idx])
    res[iso_idx] <- tx2gene_name_lookup[bare]
  }
  res
}

files <- list.files(
  cox_root,
  pattern = "cox_results(_full)?\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)
if (!length(files)) stop("No cox_results*.csv found in ", cox_root)

file_dt <- data.table(path = files)
file_dt[, base := sub("_full\\.csv$", ".csv", path)]
file_dt[, is_full := grepl("_full\\.csv$", path)]
setorder(file_dt, base, -is_full)
file_dt <- file_dt[, .SD[1], by = base]
files <- file_dt$path
n_full <- sum(file_dt$is_full)
if (!length(files)) stop("No usable cox_results files after filtering.")
say("[summary] Using %d Cox files (%d prefer *_full.csv when available)",
    length(files), n_full)

parse_meta <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  parts <- strsplit(stem, "\\.")[[1]]
  dataset <- parts[1]
  dtype <- if (length(parts) >= 2) parts[2] else "unknown"
  study <- sub("_(TMT|ITRAQ|LABELFREE).*", "", dataset)
  list(dataset = dataset, dtype = dtype, study_id = study)
}

read_cox_file <- function(path) {
  meta <- parse_meta(path)
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  if (!all(c("feature","beta","HR","p","FDR") %in% names(dt))) return(NULL)

  dt[, feature_id := as.character(feature)]
  dt[, feature := NULL]
  dt[, data_type  := meta$dtype]
  dt[, study      := meta$study_id]
  dt[, cancer     := meta$study_id]

  dt[, HR   := as.numeric(HR)]
  dt[, beta := as.numeric(beta)]
  dt[, pval := as.numeric(p)]
  dt[, qval := as.numeric(FDR)]
  dt <- dt[is.finite(HR) & HR > 0]

  dt[, logHR := log2(HR)]
  dt[, direction := fifelse(logHR > 0, "risk",
                            fifelse(logHR < 0, "protective", "neutral"))]

  dt[, gene_id := map_feature_to_gene(data_type, feature_id)]
  dt[, gene_id := as.character(gene_id)]
  dt[, gene_name := map_feature_to_gene_name(data_type, feature_id, gene_id)]
  dt[, gene_name := fifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)]
  dt[, gene    := gene_name]
  dt[, isoform := fifelse(data_type %in% c("iso_log","iso_frac"), feature_id, NA_character_)]

  if (!"n"      %in% names(dt)) dt[, n      := NA_real_]
  if (!"nevent" %in% names(dt)) dt[, nevent := NA_real_]
  dt[, `:=`(n      = as.numeric(n),
            nevent = as.numeric(nevent))]

  dt[, .(study, cancer, data_type, feature_id, gene, gene_name, isoform,
         gene_id, logHR, beta, HR, pval, qval, direction, n, nevent)]
}

say("[summary] Reading Cox outputs from %s", cox_root)
dt_list <- lapply(files, read_cox_file)
dt_list <- Filter(function(x) !is.null(x) && nrow(x) > 0, dt_list)
if (!length(dt_list)) stop("No usable Cox rows after reading inputs.")

long_dt <- rbindlist(dt_list, fill = TRUE)

long_dt[, HR     := as.numeric(HR)]
long_dt[, logHR  := as.numeric(logHR)]
long_dt[, beta   := as.numeric(beta)]
long_dt[, pval   := as.numeric(pval)]
long_dt[, qval   := as.numeric(qval)]
long_dt[, n      := as.numeric(n)]
long_dt[, nevent := as.numeric(nevent)]

if (file.exists("02_proteomics/config/cancers.yaml")) {
  cancer_cfg <- tryCatch(yaml::read_yaml("02_proteomics/config/cancers.yaml"),
                         error = function(e) NULL)
  if (!is.null(cancer_cfg) && !is.null(cancer_cfg$cancers)) {
    cancer_map <- cancer_cfg$cancers
    long_dt[, cancer := ifelse(study %in% names(cancer_map),
                               paste0(cancer_map[study], ":", study),
                               study)]
  }
}

long_path <- file.path(out_dir, "cox_results_long.csv")
fwrite(long_dt, long_path)
say("[summary] Wrote long table to %s (%d rows)", long_path, nrow(long_dt))

sig_dt <- long_dt[qval < fdr_thr & !is.na(gene_id)]
sig_path <- file.path(out_dir, "cox_results_significant.csv")
fwrite(sig_dt, sig_path)
say("[summary] Wrote significant subset (%d rows) to %s", nrow(sig_dt), sig_path)

if (!nrow(sig_dt)) {
  stop("No significant rows at FDR <", fdr_thr)
}

by_gene <- sig_dt[, .(types = list(sort(unique(data_type)))),
                  by = .(cancer, study, gene_id)]
by_gene[, pattern := vapply(types, function(x) paste(x, collapse = "+"), character(1))]

pattern_levels <- c(
  "gene",
  "iso_log",
  "iso_frac",
  "gene+iso_log",
  "gene+iso_frac",
  "iso_frac+iso_log",
  "gene+iso_log+iso_frac"
)
by_gene[, pattern := factor(pattern, levels = pattern_levels)]

pattern_counts <- by_gene[, .N, by = .(cancer, pattern)]
pattern_counts <- pattern_counts[!is.na(pattern)]
all_cancers <- sort(unique(by_gene$cancer))
pattern_counts <- pattern_counts[
  CJ(cancer = all_cancers,
     pattern = factor(pattern_levels, levels = pattern_levels),
     unique = TRUE),
  on = .(cancer, pattern)]
pattern_counts[is.na(N), N := 0]
pattern_path <- file.path(out_dir, "cox_overlap_pattern_counts.csv")
fwrite(pattern_counts, pattern_path)
say("[summary] Wrote pattern counts to %s", pattern_path)

iso_summary <- sig_dt[data_type %in% c("iso_log","iso_frac"),
                      .(
                        n_iso_sig = .N,
                        n_isoforms = uniqueN(feature_id),
                        n_datatypes = uniqueN(data_type),
                        has_risk = any(direction == "risk"),
                        has_protective = any(direction == "protective")
                      ),
                      by = .(cancer, study, gene_id, gene)]
iso_summary[, heterogeneity := has_risk & has_protective]
iso_path <- file.path(out_dir, "isoform_summary.csv")
fwrite(iso_summary, iso_path)
say("[summary] Wrote isoform summary to %s", iso_path)

say("[summary] Completed Cox global summary at %s", out_dir)
say("[summary] Run 06_cox_global_plots.R against %s to recreate the figures.", out_dir)

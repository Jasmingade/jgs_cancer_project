#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(yaml)
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop("Usage: 04_feature_overlap.R <cox_root> <out_dir> [fdr_threshold]")
}

cox_root <- args[[1]]
out_dir  <- args[[2]]
fdr_thr  <- if (length(args) == 3) as.numeric(args[[3]]) else 0.05
if (!is.finite(fdr_thr) || fdr_thr <= 0) fdr_thr <- 0.05

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

say <- function(...) message(sprintf(...))

tx2gene_path <- Sys.getenv("PROT_TX2GENE", "02_proteomics/data/raw/tx2gene.csv")
if (!file.exists(tx2gene_path)) {
  stop("Transcript-to-gene mapping not found: ", tx2gene_path)
}
tx2gene <- fread(tx2gene_path)
select_col <- function(dt, candidates, fallback_idx) {
  for (nm in candidates) if (nm %in% names(dt)) return(nm)
  names(dt)[fallback_idx]
}
iso_col  <- select_col(tx2gene, c("isoform_id", "transcript", "tx_id", "transcript_id"), 1)
gene_col <- select_col(tx2gene, c("gene_id", "gene", "gene_name"), min(2L, ncol(tx2gene)))
tx2gene[, transcript_clean := sub("\\..*$", "", get(iso_col))]
tx2gene[, gene_clean       := sub("\\..*$", "", get(gene_col))]
tx2gene_lookup <- setNames(tx2gene$gene_clean, tx2gene$transcript_clean)

parse_meta <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  parts <- strsplit(stem, "\\.")[[1]]
  dataset <- parts[1]
  dtype <- if (length(parts) >= 2) parts[2] else "unknown"
  study <- sub("_(TMT|ITRAQ|LABELFREE).*", "", dataset)
  list(dataset = dataset, dtype = dtype, study_id = study)
}

files <- list.files(cox_root, pattern = "cox_results\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(files)) {
  stop("No cox_results.csv found under: ", cox_root)
}

map_feature_to_gene <- function(dtype, feature) {
  if (dtype == "gene") {
    sub("\\..*$", "", feature)
  } else if (dtype %in% c("iso_log", "iso_frac")) {
    bare <- sub("\\..*$", "", feature)
    tx2gene_lookup[bare]
  } else {
    NA_character_
  }
}

dt_list <- vector("list", length(files))
for (i in seq_along(files)) {
  f <- files[[i]]
  meta <- parse_meta(f)
  dt <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(dt)) next
  if (!all(c("feature", "FDR") %in% names(dt))) next
  dt[, FDR := as.numeric(FDR)]
  dt <- dt[is.finite(FDR) & FDR < fdr_thr]
  if (!nrow(dt)) next
  dt[, gene_id := map_feature_to_gene(meta$dtype, feature)]
  dt <- dt[!is.na(gene_id)]
  if (!nrow(dt)) next
  dt[, c("dataset", "data_type", "study_id") := list(meta$dataset, meta$dtype, meta$study_id)]
  dt_list[[i]] <- dt[, .(study = study_id, data_type, gene_id)]
}

dt <- rbindlist(dt_list, fill = TRUE)
if (!nrow(dt)) {
  stop("No significant rows found below FDR threshold ", fdr_thr)
}

if (file.exists("02_proteomics/config/cancers.yaml")) {
  cancer_cfg <- tryCatch(yaml::read_yaml("02_proteomics/config/cancers.yaml"),
                         error = function(e) NULL)
  if (!is.null(cancer_cfg) && !is.null(cancer_cfg$cancers)) {
    cancer_map <- cancer_cfg$cancers
    dt[, cancer := ifelse(study %in% names(cancer_map),
                          paste0(cancer_map[study], ":", study),
                          study)]
  } else {
    dt[, cancer := study]
  }
} else {
  dt[, cancer := study]
}

studies <- unique(dt$study)
summary_rows <- list()
detail_rows <- list()

add_detail <- function(study_id, cancer_label, category, genes, dtype_label) {
  if (!length(genes)) return()
  detail_rows[[length(detail_rows) + 1L]] <<- data.table(
    study = study_id,
    cancer = cancer_label,
    category = category,
    gene_id = genes,
    data_types = dtype_label
  )
}

for (sid in studies) {
  sub <- dt[study == sid]
  cancer_label <- sub$cancer[1]
  S_gene    <- unique(sub[data_type == "gene",    gene_id])
  S_iso_log <- unique(sub[data_type == "iso_log", gene_id])
  S_iso_frac<- unique(sub[data_type == "iso_frac",gene_id])
  S_gene    <- S_gene[!is.na(S_gene)]
  S_iso_log <- S_iso_log[!is.na(S_iso_log)]
  S_iso_frac<- S_iso_frac[!is.na(S_iso_frac)]

  union_iso <- union(S_iso_log, S_iso_frac)
  union_gene <- union(S_gene, union_iso)

  only_gene     <- setdiff(S_gene, union_iso)
  only_iso_log  <- setdiff(S_iso_log, union(S_gene, S_iso_frac))
  only_iso_frac <- setdiff(S_iso_frac, union(S_gene, S_iso_log))
  gene_iso_log  <- setdiff(intersect(S_gene, S_iso_log), S_iso_frac)
  gene_iso_frac <- setdiff(intersect(S_gene, S_iso_frac), S_iso_log)
  iso_log_frac  <- setdiff(intersect(S_iso_log, S_iso_frac), S_gene)
  if (length(S_gene) && length(S_iso_log) && length(S_iso_frac)) {
    all_three <- Reduce(intersect, list(S_gene, S_iso_log, S_iso_frac))
  } else {
    all_three <- character()
  }

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    study = sid,
    cancer = cancer_label,
    n_gene = length(S_gene),
    n_iso_log = length(S_iso_log),
    n_iso_frac = length(S_iso_frac),
    n_gene_only = length(only_gene),
    n_iso_log_only = length(only_iso_log),
    n_iso_frac_only = length(only_iso_frac),
    n_gene_iso_log = length(gene_iso_log),
    n_gene_iso_frac = length(gene_iso_frac),
    n_iso_log_iso_frac = length(iso_log_frac),
    n_all_three = length(all_three),
    n_total_unique = length(union_gene)
  )

  add_detail(sid, cancer_label, "gene_only", only_gene, "gene")
  add_detail(sid, cancer_label, "iso_log_only", only_iso_log, "iso_log")
  add_detail(sid, cancer_label, "iso_frac_only", only_iso_frac, "iso_frac")
  add_detail(sid, cancer_label, "gene_iso_log", gene_iso_log, "gene+iso_log")
  add_detail(sid, cancer_label, "gene_iso_frac", gene_iso_frac, "gene+iso_frac")
  add_detail(sid, cancer_label, "iso_log_iso_frac", iso_log_frac, "iso_log+iso_frac")
  add_detail(sid, cancer_label, "all_three", all_three, "gene+iso_log+iso_frac")
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
detail_dt <- if (length(detail_rows)) rbindlist(detail_rows, fill = TRUE) else data.table(
  study = character(), cancer = character(), category = character(),
  gene_id = character(), data_types = character()
)

summary_path <- file.path(out_dir, "cox_feature_overlap_summary.csv")
detail_path  <- file.path(out_dir, "cox_feature_overlap_genes.csv")
fwrite(summary_dt, summary_path)
fwrite(detail_dt, detail_path)

say("[overlap] Summary written to %s", summary_path)
say("[overlap] Gene lists written to %s", detail_path)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 4) {
  stop("Usage: 03_forest_by_cancer.R <cox_root> <out_dir> [fdr_thresh] [top_n]")
}

cox_root <- args[[1]]
out_dir  <- args[[2]]
fdr_thr  <- if (length(args) >= 3) as.numeric(args[[3]]) else 0.05
top_n    <- if (length(args) == 4) as.integer(args[[4]]) else 25L

if (!is.finite(fdr_thr) || fdr_thr <= 0) fdr_thr <- 0.05
if (!is.finite(top_n)) top_n <- 25L

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tx2gene_path <- Sys.getenv("PROT_TX2GENE", "02_proteomics/data/raw/tx2gene.csv")
if (!file.exists(tx2gene_path)) {
  stop("Transcript-to-gene mapping not found: ", tx2gene_path)
}
tx2gene <- fread(tx2gene_path)
select_col <- function(dt, candidates, fallback_idx) {
  for (nm in candidates) if (nm %in% names(dt)) return(nm)
  names(dt)[fallback_idx]
}
iso_col  <- select_col(tx2gene, c("isoform_id","transcript","tx_id","transcript_id"), 1)
gene_col <- select_col(tx2gene, c("gene_id","gene","gene_name"), min(2L, ncol(tx2gene)))
tx2gene[, transcript_clean := sub("\\..*$", "", get(iso_col))]
tx2gene[, gene_clean       := sub("\\..*$", "", get(gene_col))]
tx2gene_lookup <- setNames(tx2gene$gene_clean, tx2gene$transcript_clean)

say <- function(...) message(sprintf(...))

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

files <- list.files(cox_root, pattern = "cox_results\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(files)) {
  stop("No cox_results.csv found under: ", cox_root)
}

parse_meta <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  parts <- strsplit(stem, "\\.")[[1]]
  dataset <- parts[1]
  dtype <- if (length(parts) >= 2) parts[2] else "unknown"
  study <- sub("_(TMT|ITRAQ|LABELFREE).*", "", dataset)
  list(dataset = dataset, dtype = dtype, study_id = study)
}

dt_list <- vector("list", length(files))
for (i in seq_along(files)) {
  f <- files[[i]]
  meta <- parse_meta(f)
  dt <- fread(f)
  if (!all(c("beta", "se", "HR", "FDR", "feature") %in% names(dt))) next
  dt[, HR := as.numeric(HR)]
  dt <- dt[is.finite(HR) & HR > 0]
  if (!nrow(dt)) next
  dt[, beta := as.numeric(beta)]
  dt[, se   := as.numeric(se)]
  dt[, FDR  := as.numeric(FDR)]
  dt <- dt[is.finite(beta) & is.finite(se) & is.finite(FDR)]
  if (!nrow(dt)) next
  dt <- dt[FDR < fdr_thr]
  if (!nrow(dt)) next
  dt[, gene_id := map_feature_to_gene(meta$dtype, feature)]
  dt <- dt[!is.na(gene_id)]
  if (!nrow(dt)) next
  dt[, c("dataset","data_type","study_id") := list(meta$dataset, meta$dtype, meta$study_id)]
  dt_list[[i]] <- dt
}

dt <- rbindlist(dt_list, fill = TRUE)
if (!nrow(dt)) {
  stop("No rows passed FDR filter (threshold=", fdr_thr, ").")
}

dt[, log2HR := log2(HR)]
dt[, ci_low_log2  := (beta - 1.96 * se) / log(2)]
dt[, ci_high_log2 := (beta + 1.96 * se) / log(2)]
dt[, study := study_id]
dt[, cancer := study_id]
if (file.exists("02_proteomics/config/cancers.yaml")) {
  cancer_cfg <- tryCatch(yaml::read_yaml("02_proteomics/config/cancers.yaml"), error = function(e) NULL)
  if (!is.null(cancer_cfg) && !is.null(cancer_cfg$cancers)) {
    cancer_map <- cancer_cfg$cancers
    dt[, cancer := ifelse(study %in% names(cancer_map),
                          paste0(cancer_map[study], ":", study),
                          study)]
  }
}
dt[, label := gene_id]

type_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00"
)

dt[, order_key := rank(FDR, ties.method = "first"), by = .(cancer, data_type)]
dt <- dt[order(FDR, -abs(log2HR))]

plot_forest_panel <- function(subdt, cancer_id) {
  if (!nrow(subdt)) return(NULL)
  subdt[, type_count := uniqueN(data_type), by = gene_id]
  subdt <- subdt[type_count >= 2]
  if (!nrow(subdt)) return(NULL)
  subdt <- subdt[order(-type_count, FDR, -abs(log2HR))]
  if (top_n > 0L) {
    keep_genes <- unique(subdt$gene_id)
    keep_genes <- keep_genes[seq_len(min(length(keep_genes), top_n))]
    subdt <- subdt[gene_id %in% keep_genes]
  }
  subdt[, label_f := factor(gene_id, levels = rev(unique(gene_id)))]

  ggplot(subdt,
         aes(x = log2HR, y = label_f, color = data_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = ci_low_log2, xmax = ci_high_log2),
                   height = 0.25, linewidth = 0.4,
                   position = position_dodge(width = 0.6)) +
    geom_point(size = 1.8, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = type_colors, name = "Data Type") +
    labs(
      title = sprintf("CoxPH log2(HR) – %s", cancer_id),
      subtitle = sprintf("FDR<%.2g, genes shared by ≥2 data types", fdr_thr),
      x = "log2(HR)",
      y = NULL,
      color = "Data Type"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank()
    )
}

by_cancer <- unique(dt$cancer)
say("[forest] Generating %d forest plots (≥2 datatype overlaps)", length(by_cancer))

for (cid in by_cancer) {
  sub <- dt[cancer == cid]
  p <- plot_forest_panel(sub, cid)
  if (is.null(p)) next
  fname <- sprintf("forest_%s.png", gsub("[^A-Za-z0-9]+", "_", cid))
  fpath <- file.path(out_dir, fname)
  counts <- table(sub$data_type)
  per_panel <- if (top_n > 0L) pmin(top_n, as.numeric(counts)) else as.numeric(counts)
  height <- max(5, 2 + 0.25 * sum(per_panel))
  ggsave(fpath, p, width = 10, height = 10, dpi = 300)
  say("[forest] Wrote %s", fpath)
}

say("[forest] Done.")
tx2gene_path <- Sys.getenv("PROT_TX2GENE", "02_proteomics/data/raw/tx2gene.csv")
if (!file.exists(tx2gene_path)) {
  stop("Transcript-to-gene map not found: ", tx2gene_path)
}
tx2gene <- fread(tx2gene_path)
select_col <- function(dt, candidates, fallback_idx) {
  for (nm in candidates) if (nm %in% names(dt)) return(nm)
  names(dt)[fallback_idx]
}
iso_col  <- select_col(tx2gene, c("isoform_id","transcript","tx_id","transcript_id"), 1)
gene_col <- select_col(tx2gene, c("gene_id","gene","gene_name"), min(2L, ncol(tx2gene)))
tx2gene[, transcript_clean := sub("\\..*$", "", get(iso_col))]
tx2gene[, gene_clean := sub("\\..*$", "", get(gene_col))]
tx2gene_lookup <- setNames(tx2gene$gene_clean, tx2gene$transcript_clean)

map_feature_to_gene <- function(dtype, feature) {
  if (dtype == "gene") {
    sub("\\..*$", "", feature)
  } else if (dtype %in% c("iso_log","iso_frac")) {
    bare <- sub("\\..*$", "", feature)
    tx2gene_lookup[bare]
  } else {
    NA_character_
  }
}

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(biomaRt)
  library(data.table)
})

strip_version <- function(x) sub("\\.\\d+$", "", x)

# ------------------------------------------------------------
# 1) Collect all isoform IDs from 02_norm iso_log / iso_frac
# ------------------------------------------------------------
norm_dir <- "01_transcriptomics/out/02_norm"

files_iso <- list.files(
  norm_dir,
  pattern = "TCGA_.*_(iso_log|iso_frac)\\.normalized\\.csv$",
  full.names = TRUE
)

if (length(files_iso) == 0) {
  stop(sprintf(
    "No iso_log/iso_frac normalized files found in %s.\nCheck path/pattern.",
    norm_dir
  ))
}

cat("[INFO] Found", length(files_iso),
    "isoform expression files in", norm_dir, "\n")

get_features <- function(f) {
  dt <- fread(f, nThread = 1)
  feat_col <- if ("feature" %in% names(dt)) "feature" else names(dt)[1]
  dt[[feat_col]]
}

iso_ids_raw <- unique(unlist(lapply(files_iso, get_features)))
iso_ids_raw <- iso_ids_raw[!is.na(iso_ids_raw) & nzchar(iso_ids_raw)]

iso_ids_clean <- unique(strip_version(iso_ids_raw))

cat("[INFO] Unique isoform IDs (raw):  ", length(iso_ids_raw), "\n")
cat("[INFO] Unique isoform IDs (clean):", length(iso_ids_clean), "\n")

if (length(iso_ids_clean) == 0) {
  stop("No isoform IDs found in iso_log/iso_frac normalized files.")
}

# ------------------------------------------------------------
# 2) Query Ensembl BioMart for transcript → gene mapping
# ------------------------------------------------------------
cat("[INFO] Connecting to Ensembl BioMart...\n")

mart <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl"
  # add mirror / version if needed
)

cat("[INFO] Querying BioMart for", length(iso_ids_clean), "transcript IDs...\n")

tx2gene <- getBM(
  attributes = c(
    "ensembl_transcript_id",
    "ensembl_gene_id",
    "hgnc_symbol"              # optional, nice to have
  ),
  filters   = "ensembl_transcript_id",
  values    = iso_ids_clean,
  mart      = mart
)

if (nrow(tx2gene) == 0) {
  stop("BioMart query returned 0 rows — check ENST IDs / Ensembl version.")
}

setDT(tx2gene)
setnames(
  tx2gene,
  c("ensembl_transcript_id", "ensembl_gene_id", "hgnc_symbol"),
  c("tx_id", "gene_id", "gene_name")
)

tx2gene[, tx_id   := strip_version(tx_id)]
tx2gene[, gene_id := strip_version(gene_id)]

cat("[INFO] Mapping rows:", nrow(tx2gene), "\n")

# ------------------------------------------------------------
# 3) Save mapping
# ------------------------------------------------------------
out_dir  <- "01_transcriptomics/data/raw"
out_file <- file.path(out_dir, "tx2gene.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(tx2gene, out_file)

cat("[DONE] Wrote tx2gene.csv with", nrow(tx2gene),
    "rows →", out_file, "\n")

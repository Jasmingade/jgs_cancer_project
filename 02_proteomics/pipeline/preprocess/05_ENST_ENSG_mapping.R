#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(biomaRt)
  library(data.table)
})

strip_version <- function(x) sub("\\.\\d+$", "", x)

# ------------------------------------------------------------
# 1) Collect all isoform IDs from out/isoform_data iso_log / iso_frac
# ------------------------------------------------------------
iso_dir <- "02_proteomics/out/preprocessed/filtered"

files_iso <- c(
  list.files(file.path(iso_dir, "iso_log"), pattern = "\\.csv$", full.names = TRUE, recursive = TRUE),
  list.files(file.path(iso_dir, "iso_frac"), pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
)

if (length(files_iso) == 0) {
  stop(sprintf(
    "No iso_log/iso_frac files found in %s.\nCheck path/pattern.",
    iso_dir
  ))
}

cat("[INFO] Found", length(files_iso),
    "isoform expression files in", iso_dir, "\n")

get_features <- function(f) {
  dt <- fread(f, nThread = 1)
  feat_col <- if ("feature" %in% names(dt)) "feature" else names(dt)[1]
  dt[[feat_col]]
}

iso_ids_raw <- unique(unlist(lapply(files_iso, get_features)))
iso_ids_raw <- iso_ids_raw[!is.na(iso_ids_raw) & nzchar(iso_ids_raw)]

iso_lookup <- data.table(
  iso_id_full = iso_ids_raw,
  tx_id_clean = strip_version(iso_ids_raw)
)
iso_lookup <- unique(iso_lookup, by = c("iso_id_full", "tx_id_clean"))

iso_ids_clean <- unique(iso_lookup$tx_id_clean)

cat("[INFO] Unique isoform IDs (raw):  ", length(iso_ids_raw), "\n")
cat("[INFO] Unique isoform IDs (clean):", length(iso_ids_clean), "\n")

if (length(iso_ids_clean) == 0) {
  stop("No isoform IDs found in iso_log/iso_frac files.")
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
    "ensembl_transcript_id_version",
    "ensembl_gene_id",
    "hgnc_symbol"
  ),
  filters = "ensembl_transcript_id",
  values  = iso_ids_clean,
  mart    = mart
)

if (nrow(tx2gene) == 0) {
  stop("BioMart query returned 0 rows — check ENST IDs / Ensembl version.")
}

setDT(tx2gene)
setnames(
  tx2gene,
  c("ensembl_transcript_id", "ensembl_transcript_id_version",
    "ensembl_gene_id", "hgnc_symbol"),
  c("tx_id", "tx_id_versioned", "gene_id", "gene_name")
)

tx2gene[, tx_id   := strip_version(tx_id)]
tx2gene[, gene_id := strip_version(gene_id)]

tx2gene <- merge(
  tx2gene,
  iso_lookup,
  by.x = "tx_id",
  by.y = "tx_id_clean",
  all.x = TRUE
)

setnames(tx2gene, "iso_id_full", "isoform_id")
tx2gene[is.na(isoform_id), isoform_id := tx_id_versioned]

cat("[INFO] Mapping rows:", nrow(tx2gene), "\n")

# ------------------------------------------------------------
# 3) Save mapping
# ------------------------------------------------------------
out_dir  <- "02_proteomics/data/raw"
out_file <- file.path(out_dir, "tx2gene.csv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setcolorder(tx2gene, c("isoform_id", "tx_id", "tx_id_versioned", "gene_id", "gene_name"))
fwrite(tx2gene, out_file)

fwrite(data.table(iso_id = unique(iso_lookup$iso_id_full)),
       "02_proteomics/data/raw/iso_ids_used_for_mapping.csv")

cat("[DONE] Wrote tx2gene.csv with", nrow(tx2gene),
    "rows →", out_file, "\n")

#!/usr/bin/env Rscript
# ============================================================
# 03e_prepare_gene_lists.R  (updated)
# ------------------------------------------------------------
# Goal:
#   For each (cancer, mutation_group), find genes that are
#   significant in BOTH:
#     - 03a (expression CoxPH: gene, iso_log, iso_frac)
#     - 03b (mutation-only CoxPH)
#
#   Expression features can be ENSG (gene) or ENST (isoforms).
#   ENST IDs are mapped to ENSG using tx2gene.csv.
#
#   Writes one gene_list.txt per (cancer, mut_group) containing
#   version-stripped ENSG IDs, to be used by 03e_gene_expr_mut.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

say           <- function(...) message(sprintf(...))
strip_version <- function(x) sub("\\.\\d+$", "", x)

root      <- "01_transcriptomics/out"
expr_root <- file.path(root, "03a_univariate_coxph")          # gene + iso_log + iso_frac
mut_root  <- file.path(root, "03b_mutation_univariate_coxph")
out_dir   <- file.path(root, "03e_gene_lists")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# tx2gene mapping (for iso_log / iso_frac → ENSG)
# ------------------------------------------------------------
tx2gene_file <- "01_transcriptomics/data/raw/tx2gene.csv"
tx2gene_map  <- NULL

if (file.exists(tx2gene_file)) {
  tx2gene <- fread(tx2gene_file)
  names(tx2gene) <- tolower(names(tx2gene))

  tx_col <- intersect(
    c("tx_id", "transcript_id", "ensembl_transcript_id"),
    names(tx2gene)
  )[1]
  gene_col <- intersect(
    c("gene_id", "ensembl_gene_id", "gene"),
    names(tx2gene)
  )[1]

  if (is.na(tx_col) || is.na(gene_col)) {
    warning("[WARN] tx2gene file found but no suitable tx/gene columns — isoforms cannot be mapped and will be ignored.")
  } else {
    tx2gene_map <- unique(tx2gene[, .(
      tx_core   = strip_version(get(tx_col)),
      gene_core = strip_version(get(gene_col))
    )])
    setkey(tx2gene_map, tx_core)
    say("[INFO] Loaded tx2gene mapping: %d transcript→gene rows", nrow(tx2gene_map))
  }
} else {
  warning(sprintf("[WARN] tx2gene mapping not found at %s — isoform results will be ignored.", tx2gene_file))
}

# ============================================================
# 1) Load 03a expression (gene + iso_log + iso_frac) results
# ============================================================
say("[INFO] Loading 03a expression results from: %s", expr_root)

expr_files <- list.files(
  expr_root,
  pattern    = "\\.cox_results\\.csv$",
  full.names = TRUE,
  recursive  = TRUE
)

if (length(expr_files) == 0) {
  stop("No 03a expression results found (gene / iso_log / iso_frac).")
}

load_expr_file <- function(f) {
  dt <- fread(f)
  needed <- c("feature", "FDR", "HR", "cancer", "data_type")
  if (!all(needed %in% names(dt))) {
    warning(sprintf("[WARN] Skipping 03a file (missing cols): %s", f))
    return(NULL)
  }

  # keep significant, sensible HR
  dt <- dt[FDR < 0.05 & is.finite(HR) & HR > 0]
  if (nrow(dt) == 0) return(NULL)

  # base core ID from feature
  dt[, feature_core := strip_version(feature)]
  dt[, gene_core := NA_character_]

  # gene-level: assume feature_core is already ENSG
  dt[data_type == "gene", gene_core := feature_core]

  # isoform-level: map ENST → ENSG if mapping available
  if (any(dt$data_type %in% c("iso_log", "iso_frac"))) {
    if (is.null(tx2gene_map)) {
      warning("[WARN] Isoform results detected but tx2gene mapping not loaded; dropping iso_log / iso_frac rows.")
    } else {
      iso_idx <- dt$data_type %in% c("iso_log", "iso_frac")
      dt[iso_idx, gene_core := tx2gene_map[.(feature_core), gene_core]]
    }
  }

  # drop rows without a mapped gene_core
  dt <- dt[!is.na(gene_core)]
  if (nrow(dt) == 0) return(NULL)

  dt
}

expr_list <- lapply(expr_files, load_expr_file)
expr_list <- Filter(Negate(is.null), expr_list)

if (length(expr_list) == 0) {
  stop("No significant expression features (with mapped gene_core) found in 03a.")
}

expr_sig_all <- rbindlist(expr_list, fill = TRUE)

say("[INFO] 03a: %d significant expression rows after filtering (across gene / iso_log / iso_frac).",
    nrow(expr_sig_all))
say("[INFO] 03a counts by data_type:")
print(table(expr_sig_all$data_type))

# ============================================================
# 2) Load 03b mutation significant results
# ============================================================
say("[INFO] Loading 03b mutation results from: %s", mut_root)

mut_files <- list.files(
  mut_root,
  pattern    = "\\.cox_results\\.csv$",
  full.names = TRUE,
  recursive  = TRUE
)

if (length(mut_files) == 0) {
  stop("No 03b mutation results found.")
}

load_mut_file <- function(f) {
  dt <- fread(f)
  needed <- c("feature", "FDR", "HR", "cancer", "mut_group")
  if (!all(needed %in% names(dt))) {
    warning(sprintf("[WARN] Skipping 03b file (missing cols): %s", f))
    return(NULL)
  }

  dt <- dt[FDR < 0.05 & is.finite(HR) & HR > 0]
  if (nrow(dt) == 0) return(NULL)

  # mutation features are already gene-level ENSG
  dt[, gene_core := strip_version(feature)]
  dt
}

mut_list <- lapply(mut_files, load_mut_file)
mut_list <- Filter(Negate(is.null), mut_list)

if (length(mut_list) == 0) {
  stop("No significant mutation features found in 03b.")
}

mut_sig_all <- rbindlist(mut_list, fill = TRUE)

say("[INFO] 03b: %d significant mutation rows after filtering.",
    nrow(mut_sig_all))

# ============================================================
# 3) For each (cancer, mut_group), intersect genes
# ============================================================
setkey(expr_sig_all, cancer, gene_core)
setkey(mut_sig_all,  cancer, gene_core)

combos <- unique(mut_sig_all[, .(cancer, mut_group)])
say("[INFO] Found %d (cancer, mut_group) combinations in 03b.",
    nrow(combos))

n_written <- 0L

for (i in seq_len(nrow(combos))) {
  ca <- combos$cancer[i]
  mg <- combos$mut_group[i]

  # mutation genes for this cancer & mut_group
  mut_sub <- mut_sig_all[cancer == ca & mut_group == mg]
  if (nrow(mut_sub) == 0) next

  # expression genes (from gene / iso_log / iso_frac) for this cancer
  expr_sub <- expr_sig_all[cancer == ca]
  if (nrow(expr_sub) == 0) {
    say("[INFO] No significant expression genes/isoforms for cancer %s → skip %s",
        ca, mg)
    next
  }

  common_core <- intersect(expr_sub$gene_core, mut_sub$gene_core)
  common_core <- sort(unique(common_core))

  if (length(common_core) == 0) {
    say("[INFO] No overlapping significant genes for %s / %s", ca, mg)
    next
  }

  out_file <- file.path(
    out_dir,
    sprintf("TCGA_%s_gene_%s_gene_list.txt", ca, mg)
  )

  fwrite(
    data.table(gene_id = common_core),
    out_file,
    col.names = FALSE
  )

  say("[INFO] Wrote %d genes → %s", length(common_core), out_file)
  n_written <- n_written + 1L
}

if (n_written == 0) {
  say("[WARN] No gene_list files were written (no overlapping sig genes).")
} else {
  say("[DONE] Wrote %d gene_list files into %s", n_written, out_dir)
}

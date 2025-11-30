#!/usr/bin/env Rscript
# ============================================================
# 05a_build_gene_modality_table.R
# ------------------------------------------------------------
# Build a long-format gene × modality table from Cox outputs:
#
# Uses:
#   - 03a_univariate_coxph (M1: expression-only)
#       data_type ∈ {gene, iso_log, iso_frac}
#   - 03b_mutation_univariate_coxph (M2: mutation-only)
#       mut_group ∈ {truncating_LOF, missense_or_inframe, rna, splice, silent, ...}
#   - optionally 03e_gene_expr_mut (M5: gene-matched expr+mut)
#       if directory 01_transcriptomics/out/03e_gene_expr_mut exists
#
# For each (cancer, gene, modality) we summarise:
#   - n_features (# features collapsed into this gene, e.g. isoforms)
#   - n_sig      (# with FDR < 0.05)
#   - HR         (median HR)
#   - logHR      (median log2(HR))
#   - FDR        (min FDR)
#   - delta_LL   (sum of delta log-likelihood contributions, if available)
#
# Columns in final table:
#   cancer, gene_id, model, data_type, term, mut_group, modality,
#   n_features, n_sig, HR, logHR, FDR, delta_LL, source
#
# Output:
#   01_transcriptomics/out/05_gene_modality/gene_modality_long.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

say <- function(...) message(sprintf(...))
strip_version <- function(x) sub("\\.\\d+$", "", x)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
out_root <- "01_transcriptomics/out"
expr_dir <- file.path(out_root, "03a_univariate_coxph")
mut_dir  <- file.path(out_root, "03b_mutation_univariate_coxph")
m5_dir   <- file.path(out_root, "03e_gene_expr_mut")  # optional

gm_out_dir <- file.path(out_root, "05_gene_modality")
dir.create(gm_out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Optional tx2gene mapping (for isoforms → ENSG)
# ------------------------------------------------------------
tx2gene_file <- "01_transcriptomics/data/raw/tx2gene.csv"
tx2gene <- NULL

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
    warning(sprintf(
      "[WARN] tx2gene file %s present but no suitable tx/gene cols; isoforms will NOT be collapsed.",
      tx2gene_file
    ))
    tx2gene <- NULL
  } else {
    tx2gene <- unique(tx2gene[, .(
      tx_id   = strip_version(get(tx_col)),
      gene_id = strip_version(get(gene_col))
    )])
    say("[INFO] Loaded tx2gene mapping: %d transcript→gene rows", nrow(tx2gene))
  }
} else {
  warning(sprintf(
    "[WARN] tx2gene mapping not found at %s — isoform modalities will NOT be collapsed.",
    tx2gene_file
  ))
}

# Small helper to compute median safely
safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x, na.rm = TRUE)
}

# Master list of modality tables
gm_list <- list()

# ============================================================
# 1) Model 1 – Expression-only (03a_univariate_coxph)
# ============================================================
if (dir.exists(expr_dir)) {
  # Prefer *_full outputs if present
  files_03a_full <- list.files(
    expr_dir,
    pattern = "cox_results_full\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(files_03a_full) > 0) {
    files_03a <- files_03a_full
    say("[M1] Using FULL expression results (_full.csv): %d files", length(files_03a))
  } else {
    files_03a <- list.files(
      expr_dir,
      pattern = "cox_results\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    )
    say("[M1] No *_full files; using cox_results.csv: %d files", length(files_03a))
  }

  load_03a <- function(f) {
    dt <- fread(f)
    needed <- c("feature", "beta", "HR", "p")
    if (!all(needed %in% names(dt))) {
      warning(sprintf("[M1] Skipping %s – missing required columns", f))
      return(NULL)
    }

    dt[, HR := as.numeric(HR)]
    dt[, p  := as.numeric(p)]

    if (!"FDR" %in% names(dt)) {
      dt[, FDR := p.adjust(p, "BH")]
    } else {
      dt[, FDR := as.numeric(FDR)]
    }

    # delta_LL_expr if available or derivable
    if (!"delta_LL_expr" %in% names(dt)) {
      if ("z" %in% names(dt)) {
        dt[, delta_LL_expr := (as.numeric(z)^2) / 2]
      } else {
        dt[, delta_LL_expr := NA_real_]
      }
    }

    fname <- basename(f)

    if (!"cancer" %in% names(dt)) {
      dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", fname)]
    }
    if (!"data_type" %in% names(dt)) {
      dt[, data_type :=
            fifelse(grepl("_gene", fname),     "gene",
            fifelse(grepl("_iso_log", fname),  "iso_log",
            fifelse(grepl("_iso_frac", fname), "iso_frac", "unknown")))]
    }

    dt[
      is.finite(HR) & HR > 0 &
      is.finite(p)  & p >= 0 & p <= 1
    ]
  }

  expr_list <- lapply(files_03a, load_03a)
  expr_list <- Filter(Negate(is.null), expr_list)

  if (length(expr_list) > 0) {
    expr_all <- rbindlist(expr_list, fill = TRUE)
    say("[M1] Loaded %d expression rows", nrow(expr_all))
  } else {
    expr_all <- data.table()
    say("[M1] No valid expression tables loaded.")
  }

  # -------------------------
  # M1 – gene-level expression
  # -------------------------
  expr_gene <- expr_all[data_type == "gene"]
  if (nrow(expr_gene) > 0) {
    expr_gene[, gene_id := strip_version(feature)]
    expr_gene[, logHR := log2(HR)]

    gm_M1_gene <- expr_gene[
      ,
      .(
        n_features = .N,
        n_sig      = sum(FDR < 0.05, na.rm = TRUE),
        HR         = safe_median(HR),
        logHR      = safe_median(logHR),
        FDR        = suppressWarnings(min(FDR, na.rm = TRUE)),
        delta_LL   = sum(delta_LL_expr, na.rm = TRUE)
      ),
      by = .(cancer, gene_id)
    ]

    gm_M1_gene[, model     := "M1"]
    gm_M1_gene[, term      := "expr"]
    gm_M1_gene[, data_type := "gene"]
    gm_M1_gene[, mut_group := NA_character_]
    gm_M1_gene[, source    := "03a"]

    gm_list[["M1_gene"]] <- gm_M1_gene
    say("[M1] Added gene-level entries: %d rows", nrow(gm_M1_gene))
  } else {
    say("[M1] No gene-level expression rows.")
  }

  # -------------------------
  # M1 – iso_log expression (collapsed to ENSG)
  # -------------------------
  expr_iso_log <- expr_all[data_type == "iso_log"]
  if (nrow(expr_iso_log) > 0 && !is.null(tx2gene)) {
    expr_iso_log[, tx_id := strip_version(feature)]
    expr_iso_log <- merge(expr_iso_log, tx2gene, by = "tx_id", all.x = TRUE)
    expr_iso_log <- expr_iso_log[!is.na(gene_id)]

    expr_iso_log[, gene_id := strip_version(gene_id)]
    expr_iso_log[, logHR := log2(HR)]

    gm_M1_iso_log <- expr_iso_log[
      ,
      .(
        n_features = .N,
        n_sig      = sum(FDR < 0.05, na.rm = TRUE),
        HR         = safe_median(HR),
        logHR      = safe_median(logHR),
        FDR        = suppressWarnings(min(FDR, na.rm = TRUE)),
        delta_LL   = sum(delta_LL_expr, na.rm = TRUE)
      ),
      by = .(cancer, gene_id)
    ]

    gm_M1_iso_log[, model     := "M1"]
    gm_M1_iso_log[, term      := "expr"]
    gm_M1_iso_log[, data_type := "iso_log"]
    gm_M1_iso_log[, mut_group := NA_character_]
    gm_M1_iso_log[, source    := "03a"]

    gm_list[["M1_iso_log"]] <- gm_M1_iso_log
    say("[M1] Added iso_log entries (collapsed to ENSG): %d rows", nrow(gm_M1_iso_log))
  } else if (nrow(expr_iso_log) > 0) {
    say("[M1] iso_log results exist but tx2gene is missing – skipping iso_log collapse.")
  }

  # -------------------------
  # M1 – iso_frac expression (collapsed to ENSG)
  # -------------------------
  expr_iso_frac <- expr_all[data_type == "iso_frac"]
  if (nrow(expr_iso_frac) > 0 && !is.null(tx2gene)) {
    expr_iso_frac[, tx_id := strip_version(feature)]
    expr_iso_frac <- merge(expr_iso_frac, tx2gene, by = "tx_id", all.x = TRUE)
    expr_iso_frac <- expr_iso_frac[!is.na(gene_id)]

    expr_iso_frac[, gene_id := strip_version(gene_id)]
    expr_iso_frac[, logHR := log2(HR)]

    gm_M1_iso_frac <- expr_iso_frac[
      ,
      .(
        n_features = .N,
        n_sig      = sum(FDR < 0.05, na.rm = TRUE),
        HR         = safe_median(HR),
        logHR      = safe_median(logHR),
        FDR        = suppressWarnings(min(FDR, na.rm = TRUE)),
        delta_LL   = sum(delta_LL_expr, na.rm = TRUE)
      ),
      by = .(cancer, gene_id)
    ]

    gm_M1_iso_frac[, model     := "M1"]
    gm_M1_iso_frac[, term      := "expr"]
    gm_M1_iso_frac[, data_type := "iso_frac"]
    gm_M1_iso_frac[, mut_group := NA_character_]
    gm_M1_iso_frac[, source    := "03a"]

    gm_list[["M1_iso_frac"]] <- gm_M1_iso_frac
    say("[M1] Added iso_frac entries (collapsed to ENSG): %d rows", nrow(gm_M1_iso_frac))
  } else if (nrow(expr_iso_frac) > 0) {
    say("[M1] iso_frac results exist but tx2gene is missing – skipping iso_frac collapse.")
  }

} else {
  say(sprintf("[M1] Expression directory not found: %s", expr_dir))
}

# ============================================================
# 2) Model 2 – Mutation-only (03b_mutation_univariate_coxph)
# ============================================================
if (dir.exists(mut_dir)) {
  files_03b_full <- list.files(
    mut_dir,
    pattern = "cox_results_full\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(files_03b_full) > 0) {
    files_03b <- files_03b_full
    say("[M2] Using FULL mutation results (_full.csv): %d files", length(files_03b))
  } else {
    files_03b <- list.files(
      mut_dir,
      pattern = "cox_results\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    )
    say("[M2] No *_full files; using cox_results.csv: %d files", length(files_03b))
  }

  load_03b <- function(f) {
    dt <- fread(f)
    needed <- c("feature", "beta", "HR", "p")
    if (!all(needed %in% names(dt))) {
      warning(sprintf("[M2] Skipping %s – missing required columns", f))
      return(NULL)
    }

    dt[, HR := as.numeric(HR)]
    dt[, p  := as.numeric(p)]

    if (!"FDR" %in% names(dt)) {
      dt[, FDR := p.adjust(p, "BH")]
    } else {
      dt[, FDR := as.numeric(FDR)]
    }

    # delta_LL_mut if available or derivable
    if (!"delta_LL_mut" %in% names(dt)) {
      if ("z" %in% names(dt)) {
        dt[, delta_LL_mut := (as.numeric(z)^2) / 2]
      } else {
        dt[, delta_LL_mut := NA_real_]
      }
    }

    fname <- basename(f)

    if (!"cancer" %in% names(dt)) {
      dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", fname)]
    }
    if (!"mut_group" %in% names(dt)) {
      dt[, mut_group :=
            sub("^TCGA_[A-Z0-9]+_mutation_(.*)\\.cox_results.*$", "\\1", fname)]
    }

    dt <- dt[
      is.finite(HR) & HR > 0 &
      is.finite(p)  & p >= 0 & p <= 1
    ]

    dt[, gene_id := strip_version(feature)]
    dt[, logHR := log2(HR)]

    dt
  }

  mut_list <- lapply(files_03b, load_03b)
  mut_list <- Filter(Negate(is.null), mut_list)

  if (length(mut_list) > 0) {
    mut_all <- rbindlist(mut_list, fill = TRUE)
    say("[M2] Loaded %d mutation rows", nrow(mut_all))

    gm_M2 <- mut_all[
      ,
      .(
        n_features = .N,
        n_sig      = sum(FDR < 0.05, na.rm = TRUE),
        HR         = safe_median(HR),
        logHR      = safe_median(logHR),
        FDR        = suppressWarnings(min(FDR, na.rm = TRUE)),
        delta_LL   = sum(delta_LL_mut, na.rm = TRUE)
      ),
      by = .(cancer, gene_id, mut_group)
    ]

    gm_M2[, model     := "M2"]
    gm_M2[, term      := "mut"]
    gm_M2[, data_type := "mutation"]
    gm_M2[, source    := "03b"]

    gm_list[["M2_mut"]] <- gm_M2
    say("[M2] Added mutation-only entries: %d rows", nrow(gm_M2))
  } else {
    say("[M2] No valid mutation tables loaded.")
  }

} else {
  say(sprintf("[M2] Mutation directory not found: %s", mut_dir))
}

# ============================================================
# 3) Model 5 – Gene-matched expr + mut (03e_gene_expr_mut, optional)
# ============================================================
if (dir.exists(m5_dir)) {
  files_03e_full <- list.files(
    m5_dir,
    pattern = "cox_results_full\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(files_03e_full) > 0) {
    files_03e <- files_03e_full
    say("[M5] Using FULL gene-matched expr+mut results (_full.csv): %d files", length(files_03e))
  } else {
    files_03e <- list.files(
      m5_dir,
      pattern = "cox_results\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    )
    say("[M5] No *_full files; using cox_results.csv: %d files", length(files_03e))
  }

  load_03e <- function(f) {
    dt <- fread(f)

    needed <- c("feature", "HR_expr", "p_expr", "HR_mut", "p_mut")
    if (!all(needed %in% names(dt))) {
      warning(sprintf("[M5] Skipping %s – missing required columns", f))
      return(NULL)
    }

    dt[, HR_expr := as.numeric(HR_expr)]
    dt[, HR_mut  := as.numeric(HR_mut)]
    dt[, p_expr  := as.numeric(p_expr)]
    dt[, p_mut   := as.numeric(p_mut)]

    # FDR if missing
    if (!"FDR_expr" %in% names(dt)) {
      dt[, FDR_expr := p.adjust(p_expr, "BH")]
    } else {
      dt[, FDR_expr := as.numeric(FDR_expr)]
    }
    if (!"FDR_mut" %in% names(dt)) {
      dt[, FDR_mut := p.adjust(p_mut, "BH")]
    } else {
      dt[, FDR_mut := as.numeric(FDR_mut)]
    }

    # delta_LL_expr / delta_LL_mut if missing
    if (!"delta_LL_expr" %in% names(dt)) {
      if ("z_expr" %in% names(dt)) {
        dt[, delta_LL_expr := (as.numeric(z_expr)^2) / 2]
      } else {
        dt[, delta_LL_expr := NA_real_]
      }
    }
    if (!"delta_LL_mut" %in% names(dt)) {
      if ("z_mut" %in% names(dt)) {
        dt[, delta_LL_mut := (as.numeric(z_mut)^2) / 2]
      } else {
        dt[, delta_LL_mut := NA_real_]
      }
    }

    # cancer / mut_group if not stored in file
    fname <- basename(f)
    if (!"cancer" %in% names(dt)) {
      dt[, cancer := sub("^TCGA_([A-Z0-9]+).*", "\\1", fname)]
    }
    if (!"mut_group" %in% names(dt)) {
      # allow pattern ...ensembl_<group>.cox_results.csv if used
      dt[, mut_group := sub("^.*ensembl_(.*)\\.cox_results.*$", "\\1", fname)]
    }

    dt[, gene_id := strip_version(feature)]

    dt
  }

  m5_list <- lapply(files_03e, load_03e)
  m5_list <- Filter(Negate(is.null), m5_list)

  if (length(m5_list) > 0) {
    m5_all <- rbindlist(m5_list, fill = TRUE)
    say("[M5] Loaded %d gene-matched expr+mut rows", nrow(m5_all))

    # ----- expression term in gene-matched model -----
    m5_expr <- m5_all[
      ,
      .(
        n_features = .N,
        n_sig      = sum(FDR_expr < 0.05, na.rm = TRUE),
        HR         = safe_median(HR_expr),
        logHR      = safe_median(log2(HR_expr)),
        FDR        = suppressWarnings(min(FDR_expr, na.rm = TRUE)),
        delta_LL   = sum(delta_LL_expr, na.rm = TRUE)
      ),
      by = .(cancer, gene_id, mut_group)
    ]

    m5_expr[, model     := "M5"]
    m5_expr[, term      := "expr"]
    m5_expr[, data_type := "gene"]
    m5_expr[, source    := "03e"]

    gm_list[["M5_expr"]] <- m5_expr
    say("[M5] Added gene-matched expr entries: %d rows", nrow(m5_expr))

    # ----- mutation term in gene-matched model -----
    m5_mut <- m5_all[
      ,
      .(
        n_features = .N,
        n_sig      = sum(FDR_mut < 0.05, na.rm = TRUE),
        HR         = safe_median(HR_mut),
        logHR      = safe_median(log2(HR_mut)),
        FDR        = suppressWarnings(min(FDR_mut, na.rm = TRUE)),
        delta_LL   = sum(delta_LL_mut, na.rm = TRUE)
      ),
      by = .(cancer, gene_id, mut_group)
    ]

    m5_mut[, model     := "M5"]
    m5_mut[, term      := "mut"]
    m5_mut[, data_type := "mutation"]
    m5_mut[, source    := "03e"]

    gm_list[["M5_mut"]] <- m5_mut
    say("[M5] Added gene-matched mut entries: %d rows", nrow(m5_mut))

  } else {
    say("[M5] No valid 03e tables loaded.")
  }

} else {
  say(sprintf("[M5] Optional directory not found, skipping 03e: %s", m5_dir))
}

# ============================================================
# Combine all modalities and add 'modality' label
# ============================================================
if (length(gm_list) == 0) {
  stop("No gene × modality entries built – nothing to write.")
}

gene_modality <- rbindlist(gm_list, use.names = TRUE, fill = TRUE)

# Ensure core columns exist
if (!"mut_group" %in% names(gene_modality))
  gene_modality[, mut_group := NA_character_]

# Build human-readable modality label
gene_modality[, modality := ifelse(
  is.na(mut_group) | mut_group == "",
  sprintf("%s:%s:%s", model, data_type, term),
  sprintf("%s:%s:%s:%s", model, data_type, term, mut_group)
)]

# Column order
setcolorder(
  gene_modality,
  c("cancer", "gene_id", "model", "data_type", "term",
    "mut_group", "modality", "n_features", "n_sig",
    "HR", "logHR", "FDR", "delta_LL", "source")
)

# ============================================================
# Write output
# ============================================================
out_file <- file.path(gm_out_dir, "gene_modality_long.csv")
fwrite(gene_modality, out_file)
say("[DONE] Wrote gene × modality table → %s", out_file)

# Done.

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(future.apply)
  library(qvalue)
})

# ============================================================
# Utility: Safe loader for Cox results
# ============================================================
load_cox_file <- function(file, cancer, dtype) {
  if (!file.exists(file)) {
    message("[WARN] Missing file: ", file)
    return(NULL)
  }
  dt <- fread(file)
  dt[, cancer := cancer]
  if (!missing(dtype)) dt[, dtype := dtype]
  return(dt)
}

# ============================================================
# Add FDR (q-value)
# ============================================================
add_fdr <- function(dt, p_col = "pval") {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  dt[, FDR := qvalue::qvalue(get(p_col))$qvalues]
  return(dt)
}

# ============================================================
# Format HR values (log10 scale)
# ============================================================
format_hr <- function(hr) log10(hr)

# ============================================================
# facet labels (pretty formatting)
# ============================================================
pretty_cancer <- function(x) {
  gsub("TCGA_", "", x)
}

# ============================================================
# Multiomics overlap input generator
# Accepts: list of sets, returns binary matrix
# ============================================================
make_overlap_matrix <- function(list_of_sets, all_features) {
  out <- data.table(feature_id = all_features)
  for (name in names(list_of_sets)) {
    out[, (name) := as.integer(feature_id %in% list_of_sets[[name]])]
  }
  return(out)
}

# ============================================================
# Utility: Save plot
# ============================================================
save_plot <- function(p, file, width=12, height=7, dpi=300) {
  ggsave(file, p, width = width, height = height, dpi = dpi, bg = "white")
  message("[SAVED] ", file)
}

# ============================================================
# Utility: Combine many PNGs into a PDF (optional)
# ============================================================
save_pdf_from_pngs <- function(png_files, pdf_file) {
  if (length(png_files) == 0) return(NULL)
  system(
    paste("convert", paste(png_files, collapse=" "), pdf_file)
  )
  message("[PDF CREATED] ", pdf_file)
}

# ============================================================
# Utility: progress message
# ============================================================
say <- function(...) cat("[INFO]", sprintf(...), "\n")

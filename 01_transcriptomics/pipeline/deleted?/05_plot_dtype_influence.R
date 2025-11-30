#!/usr/bin/env Rscript

# 05_plot_dtype_influence.R
# Compare per-cancer influence of datatypes (gene, iso_log, iso_frac) using univariate Cox results.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(stringr)
})

IN_MASTER <- "01_transcriptomics/out/04_univariate_collect/univariate_master_all.csv"
OUT_DIR   <- "01_transcriptomics/out/05_dtype_influence"
SIG_Q     <- 0.10      # use features with q_within < SIG_Q
TOP_N     <- 20        # take top-N (by p) within each cancer×dtype to compute median |logHR|
METRIC    <- "median_abs_logHR_topN"  # choices: "median_abs_logHR_topN", "max_abs_logHR_sig", "n_sig"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dt <- fread(IN_MASTER)
stopifnot(all(c("cancer","dtype","feature","p","q_within","HR","logHR","cindex") %in% names(dt)))

# Keep only significant rows (panel-level FDR)
sig <- dt[q_within < SIG_Q & is.finite(logHR)]

# Guard: if some panels have no significant hits, keep the best p anyway for comparability
fill_best <- dt[!is.finite(logHR), logHR := NA_real_][
  , .SD[order(p)][1], by = .(cancer, dtype)]
fill_best <- fill_best[!is.na(logHR)]
sig_complete <- rbind(sig, fill_best[!.(sig$cancer, sig$dtype), on=.(cancer,dtype)], fill = TRUE)

# Compute metrics per cancer×dtype
# a) median |logHR| over top-N most significant hits (robust)
topN <- sig_complete[order(p)][, head(.SD, TOP_N), by = .(cancer, dtype)]
agg_topN <- topN[, .(
  median_abs_logHR_topN = median(abs(logHR), na.rm = TRUE),
  n_topN = .N
), by = .(cancer, dtype)]

# b) max |logHR| among significant features
agg_max <- sig_complete[, .(
  max_abs_logHR_sig = max(abs(logHR), na.rm = TRUE)
), by = .(cancer, dtype)]

# c) number of significant features
agg_nsig <- dt[, .(n_sig = sum(q_within < SIG_Q, na.rm = TRUE)), by = .(cancer, dtype)]

# Combine
sum_dt <- Reduce(function(x,y) merge(x,y,by=c("cancer","dtype"), all=TRUE),
                 list(agg_topN, agg_max, agg_nsig))

# Choose y metric
ycol <- switch(METRIC,
  "median_abs_logHR_topN" = "median_abs_logHR_topN",
  "max_abs_logHR_sig"     = "max_abs_logHR_sig",
  "n_sig"                 = "n_sig",
  "median_abs_logHR_topN"
)

# Order cancers by overall max influence for nicer plotting
ord <- sum_dt[, .(score = max(get(ycol), na.rm = TRUE)), by = cancer][order(-score)]$cancer
sum_dt[, cancer := factor(cancer, levels = ord)]
sum_dt[, dtype  := factor(dtype, levels = c("gene","iso_log","iso_frac"))]

# Plot
p <- ggplot(sum_dt, aes(x = cancer, y = .data[[ycol]], fill = dtype)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = ifelse(is.finite(n_sig), paste0("n=", n_sig), "")),
            position = position_dodge(width = 0.75), vjust = -0.2, size = 3) +
  labs(
    title = sprintf("Datatype influence by cancer (metric: %s, q<%.2f)", ycol, SIG_Q),
    x = "Cancer",
    y = if (ycol == "n_sig") "Significant features (count)" else "|log(HR)|",
    fill = "Datatype"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

# Save
png_path <- file.path(OUT_DIR, sprintf("datatype_influence_%s.png", ycol))
pdf_path <- file.path(OUT_DIR, sprintf("datatype_influence_%s.pdf", ycol))
ggsave(png_path, p, width = 11, height = 5, dpi = 180)
ggsave(pdf_path, p, width = 11, height = 5)

# Also write the summary table used for plotting
fwrite(sum_dt[order(cancer, dtype)], file.path(OUT_DIR, sprintf("datatype_influence_%s.csv", ycol)))

message("Wrote: ", png_path, " and ", pdf_path)

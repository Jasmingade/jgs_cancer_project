#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------
# CLI + setup
# ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: 05_isoform_lollipops.R <cox_long_csv> <cox_sig_csv> <out_dir>")
}

cox_long_path <- args[[1]]
cox_sig_path  <- args[[2]]
summary_dir   <- args[[3]]

dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

say <- function(...) message(sprintf(...))

say("[init] long = %s", cox_long_path)
say("[init] sig  = %s", cox_sig_path)
say("[init] out  = %s", summary_dir)

# Colours (consistent with earlier scripts)
type_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00"
)

# ------------------------------------------------------------------
# Load long + significant Cox tables
# ------------------------------------------------------------------
long_dt <- fread(cox_long_path)
sig_dt  <- fread(cox_sig_path)

# basic sanity
if (!all(c("cancer","data_type","gene_id","gene","feature_id","logHR","direction") %in% names(long_dt))) {
  stop("long_dt is missing required columns (expected at least: cancer, data_type, gene_id, gene, feature_id, logHR, direction)")
}
if (!all(c("cancer","data_type","gene_id","gene","feature_id","logHR","direction") %in% names(sig_dt))) {
  stop("sig_dt is missing required columns (expected at least: cancer, data_type, gene_id, gene, feature_id, logHR, direction)")
}

# ------------------------------------------------------------------
# Isoform-focused summaries
# ------------------------------------------------------------------
iso_summary <- sig_dt[data_type %in% c("iso_log","iso_frac"),
                      .(
                        n_iso_sig      = .N,
                        n_isoforms     = uniqueN(feature_id),
                        n_datatypes    = uniqueN(data_type),
                        has_risk       = any(direction == "risk"),
                        has_protective = any(direction == "protective")
                      ),
                      by = .(cancer, study, gene_id, gene)]
iso_summary[, heterogeneity := has_risk & has_protective]

# Base per-gene isoform summary (reused in C2/C3)
iso_gene_summary_base <- sig_dt[data_type %in% c("iso_log", "iso_frac"),
                                .(
                                  n_isoforms = uniqueN(feature_id),
                                  n_risk     = sum(direction == "risk"),
                                  n_prot     = sum(direction == "protective")
                                ),
                                by = .(cancer, gene_id, gene)]
iso_gene_summary_base[, total := n_risk + n_prot]

# Gene-level Cox results (reused)
if (!"pval" %in% names(long_dt) && "p" %in% names(long_dt)) {
  long_dt[, pval := as.numeric(p)]
}
gene_info_all <- long_dt[
  data_type == "gene",
  .SD[which.min(pval)],
  by = .(cancer, gene_id, gene)
]

# ------------------------------------------------------------------
# DIAGNOSTICS: how many genes per cancer for each plot flavour
# ------------------------------------------------------------------

# C1: heterogeneous genes with >= min_isoforms_mixed
min_isoforms_mixed <- 3

C1_pre <- iso_summary[
  n_isoforms >= min_isoforms_mixed & heterogeneity == TRUE,
  .(n_genes_hetero = uniqueN(gene_id)),
  by = cancer
][order(-n_genes_hetero)]

C1_post <- merge(
  iso_summary[
    n_isoforms >= min_isoforms_mixed & heterogeneity == TRUE,
    .(cancer, gene_id, gene)
  ],
  gene_info_all[, .(cancer, gene_id, gene)],
  by = c("cancer", "gene_id", "gene"),
  all = FALSE
)[, .(n_genes_plottable = uniqueN(gene_id)), by = cancer][order(-n_genes_plottable)]

fwrite(C1_pre,  file.path(summary_dir, "C1_hetero_genes_per_cancer_premerge.csv"))
fwrite(C1_post, file.path(summary_dir, "C1_hetero_genes_per_cancer_postmerge.csv"))
say("[diag] Wrote C1 diagnostics")

# C2: both directions, stricter
min_isoforms <- 5

iso_gene_summary_C2 <- iso_gene_summary_base[
  total > 0 & n_isoforms >= min_isoforms & n_risk > 0 & n_prot > 0
]

C2_pre <- iso_gene_summary_C2[
  , .(n_genes_mixeddir = uniqueN(gene_id)),
  by = cancer
][order(-n_genes_mixeddir)]

C2_post <- merge(
  iso_gene_summary_C2[, .(cancer, gene_id, gene)],
  gene_info_all[, .(cancer, gene_id, gene)],
  by = c("cancer", "gene_id", "gene"),
  all = FALSE
)[, .(n_genes_plottable = uniqueN(gene_id)), by = cancer][order(-n_genes_plottable)]

fwrite(C2_pre,  file.path(summary_dir, "C2_mixed_genes_per_cancer_premerge.csv"))
fwrite(C2_post, file.path(summary_dir, "C2_mixed_genes_per_cancer_postmerge.csv"))
say("[diag] Wrote C2 diagnostics")

# C3: softer criteria, any direction, classified by prop_risk
iso_gene_summary_C3 <- iso_gene_summary_base[
  total > 0 & n_isoforms >= min_isoforms
]
iso_gene_summary_C3[, prop_risk := n_risk / total]
iso_gene_summary_C3[, category := fifelse(
  prop_risk >= 2/3, "Risk-dominated",
  fifelse(prop_risk <= 1/3, "Protective-dominated", "Balanced mixed")
)]

C3_pre <- iso_gene_summary_C3[
  , .(n_genes_eligible = uniqueN(gene_id)),
  by = cancer
][order(-n_genes_eligible)]

C3_pre_cat <- iso_gene_summary_C3[
  , .(n_genes = uniqueN(gene_id)),
  by = .(cancer, category)
][order(cancer, category)]

fwrite(C3_pre,     file.path(summary_dir, "C3_genes_per_cancer_premerge.csv"))
fwrite(C3_pre_cat, file.path(summary_dir, "C3_genes_per_cancer_category_premerge.csv"))
say("[diag] Wrote C3 diagnostics")

say("[done] Isoform lollipop pipeline finished.")

















# ------------------------------------------------------------------
# Removed plots from 06_cox_global_plots.R (insert back later)
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# Overlap composition plot (Pattern counts per cancer)
# ------------------------------------------------------------------
pattern_levels <- c(
  "gene",
  "iso_log",
  "iso_frac",
  "gene+iso_log",
  "gene+iso_frac",
  "iso_frac+iso_log",
  "gene+iso_log+iso_frac"
)

if (file.exists(pattern_path)) {
  pattern_counts <- fread(pattern_path)
  pattern_counts[, pattern := factor(pattern, levels = pattern_levels)]
  pattern_counts[, cancer := factor(cancer, levels = sort(unique(sig_dt$cancer)))]
} else {
  by_gene <- sig_dt[, .(types = list(sort(unique(data_type)))),
                    by = .(cancer, study, gene_id)]
  by_gene[, pattern := vapply(types, function(x) paste(x, collapse = "+"), character(1))]
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
}

pattern_palette <- c(
  "gene" = "#0072B2",
  "iso_log" = "#009E73",
  "iso_frac" = "#E69F00",
  "gene+iso_log" = "#54B2BE",
  "gene+iso_frac" = "#7DB4D6",
  "iso_frac+iso_log" = "#D6B642",
  "gene+iso_log+iso_frac" = "#9B89C6"
)

plot_overlap <- ggplot(pattern_counts,
                       aes(x = cancer, y = N, fill = pattern)) +
  geom_col(position = "fill", colour = "white") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = pattern_palette, drop = FALSE, na.translate = FALSE) +
  labs(
    x = "Cancer",
    y = "Proportion of significant genes",
    fill = "Significant in",
    title = "Overlap of FDR-significant genes across data types"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(summary_dir, "plot_overlap_composition.png"),
       plot_overlap, width = 15, height = 9, dpi = 300)
say("[plots] Saved overlap composition plot")

# ------------------------------------------------------------------
# Jaccard similarity heatmap
# ------------------------------------------------------------------
jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (!length(a) && !length(b)) return(NA_real_)
  if (!length(union(a, b))) return(NA_real_)
  if (!length(a) || !length(b)) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

jac <- sig_dt[, {
  g  <- unique(gene_id[data_type == "gene"])
  il <- unique(gene_id[data_type == "iso_log"])
  ifr<- unique(gene_id[data_type == "iso_frac"])
  data.table(
    pair = c("gene-iso_log","gene-iso_frac","iso_log-iso_frac"),
    J    = c(jaccard(g, il), jaccard(g, ifr), jaccard(il, ifr))
  )
}, by = cancer]

plot_jaccard <- ggplot(jac, aes(x = pair, y = cancer, fill = J)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Jaccard\nsimilarity", na.value = "grey90") +
  labs(x = "", y = "Cancer", title = "Similarity of significant gene sets") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(summary_dir, "plot_jaccard_heatmap.png"),
       plot_jaccard, width = 7, height = max(4, 0.4 * length(unique(jac$cancer))), dpi = 300)
say("[plots] Saved Jaccard heatmap")











# ------------------------------------------------------------------
# Significant count plots (with/without direction)
# ------------------------------------------------------------------
count_dt <- sig_dt[, .N, by = .(cancer, data_type, direction)]
count_dt[, cancer := factor(cancer, levels = sort(unique(cancer)))]
direction_levels <- c("risk", "protective", "neutral")
count_dt[, direction := factor(direction, levels = direction_levels)]

plot_counts <- ggplot(
  count_dt,
  aes(x = cancer, y = N, fill = data_type, colour = direction)
) +
  geom_col(position = position_dodge(width = 0.85), width = 0.6, linewidth = 0.5) +
  scale_fill_manual(values = type_colors, drop = FALSE, name = "Data type") +
  scale_colour_manual(values = direction_cols, drop = FALSE, name = "Direction") +
  scale_y_log10(labels = comma_format(accuracy = 1)) +
  labs(
    x = "Cancer",
    y = "# of significant features",
    title = sprintf("Significant features per cancer and data type (%s)", fdr_label)
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(summary_dir, "plot_significant_counts.png"),
       plot_counts, width = 18, height = 8, dpi = 300)
say("[plots] Saved directional counts plot")

count_simple <- sig_dt[, .N, by = .(cancer, data_type)]
count_simple <- count_simple[N > 0]
count_simple[, cancer := factor(cancer, levels = sort(unique(cancer)))]

plot_counts_simple <- ggplot(count_simple,
                             aes(x = cancer, y = N, fill = data_type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = type_colors, drop = FALSE, name = "Data type") +
  scale_y_log10(labels = comma_format(accuracy = 1)) +
  labs(
    x = "Cancer",
    y = "# of significant features",
    title = sprintf("Significant features per cancer (%s)", fdr_label)
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(summary_dir, "plot_significant_counts_simple.png"),
       plot_counts_simple, width = 18, height = 8, dpi = 300)
say("[plots] Saved stacked counts plot")

# ------------------------------------------------------------------
# Diverging proportion plot (risk vs protective)
# ------------------------------------------------------------------
sum_dt <- sig_dt[, .(
  risk       = sum(direction == "risk"),
  protective = sum(direction == "protective"),
  neutral    = sum(direction == "neutral")
), by = .(cancer, data_type)]

pyr_dt <- melt(
  sum_dt,
  id.vars = c("cancer", "data_type"),
  variable.name = "direction",
  value.name = "n"
)

pyr_dt[, n_signed := fifelse(direction == "protective", -n,
                      fifelse(direction == "risk", n, 0))]

cancer_order <- sig_dt[, .N, by = cancer][order(-N), cancer]
pyr_dt[, cancer := factor(cancer, levels = cancer_order)]
pyr_dt[, direction := factor(direction, levels = direction_levels)]

pyr_prop <- copy(pyr_dt)
pyr_prop[, total_abs := sum(abs(n_signed)), by = .(cancer, data_type)]
pyr_prop[, prop_signed := ifelse(total_abs > 0, n_signed / total_abs, 0)]

pyr_prop_main <- pyr_prop[direction %in% c("risk", "protective")]

plot_diverging <- ggplot(pyr_prop_main,
                         aes(x = prop_signed, y = cancer, fill = direction)) +
  geom_col(width = 0.6, colour = "white") +
  facet_wrap(~ data_type, ncol = 1) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     name = "Proportion of significant features") +
  scale_fill_manual(values = direction_cols[c("risk", "protective")],
                    drop = FALSE, name = "Direction") +
  labs(y = "Cancer",
       title = "Proportion of risk vs protective features") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(summary_dir, "plot_diverging_props.png"),
       plot_diverging, width = 10, height = 8, dpi = 300)
say("[plots] Saved diverging proportion plot")

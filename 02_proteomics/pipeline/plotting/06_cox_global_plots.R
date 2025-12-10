#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || length(args) > 2) {
  stop("Usage: 06_cox_global_plots.R <summary_dir> [fdr_thresh]")
}

summary_dir <- args[[1]]
fdr_thr <- if (length(args) == 2) as.numeric(args[[2]]) else 0.05
if (!is.finite(fdr_thr) || fdr_thr <= 0) fdr_thr <- 0.05

say <- function(...) message(sprintf(...))

sig_path  <- file.path(summary_dir, "cox_results_significant.csv")
long_path <- file.path(summary_dir, "cox_results_long.csv")
pattern_path <- file.path(summary_dir, "cox_overlap_pattern_counts.csv")

for (p in list(sig_path, long_path)) {
  if (!file.exists(p)) stop("Required summary file missing: ", p)
}

sig_dt  <- fread(sig_path)
long_dt <- fread(long_path)
if (!nrow(sig_dt)) stop("No significant rows available for plotting.")

sig_dt[, cancer := factor(cancer, levels = sort(unique(cancer)))]
long_dt[, cancer := factor(cancer, levels = sort(unique(cancer)))]

fdr_label <- sprintf("FDR<%.3f", fdr_thr)

type_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00"
)

direction_cols <- c(
  risk        = "#c90028",
  protective  = "#E69F00",
  neutral     = "#6C757D"
)

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
# Isoform-focused summaries and lollipop plots
# ------------------------------------------------------------------
iso_summary <- sig_dt[data_type %in% c("iso_log","iso_frac"),
                      .(
                        n_iso_sig    = .N,
                        n_isoforms   = uniqueN(feature_id),
                        n_datatypes  = uniqueN(data_type),
                        has_risk     = any(direction == "risk"),
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
gene_info_all <- long_dt[
  data_type == "gene",
  .SD[which.min(pval)],           # best p-value per cancer x gene
  by = .(cancer, gene_id, gene)
]

# ------------------------------------------------------------------
# Helper: lollipop plot with sticks, isoform medians, gene line
# ------------------------------------------------------------------
make_iso_lolli <- function(iso_plot_dt, gene_info, cap_val, title, outfile, ncol = 3) {
  if (!nrow(iso_plot_dt) || !nrow(gene_info)) {
    say("[plots] Skipping %s: no data", title)
    return(invisible(NULL))
  }

  iso_dt  <- data.table::copy(iso_plot_dt)
  iso_dt[, logHR_cap := pmax(pmin(logHR, cap_val), -cap_val)]

  gene_dt <- data.table::copy(gene_info)
  gene_dt[, logHR_cap := pmax(pmin(logHR, cap_val), -cap_val)]

  # Isoform medians per facet x data_type
  iso_line_dt <- iso_dt[
    , .(iso_median = median(logHR_cap, na.rm = TRUE)),
    by = .(facet_label, data_type)
  ]

  if (!"gene_dir" %in% names(gene_dt)) {
    gene_dt[, gene_dir := ifelse(logHR > 0, "risk", "protective")]
  }

  p <- ggplot(
    iso_dt,
    aes(x = reorder(feature_id, logHR_cap),
        y = logHR_cap,
        colour = data_type)
  ) +
    # baseline HR = 1
    geom_hline(yintercept = 0, linetype = 2, colour = "grey80") +
    # isoform medians
    geom_hline(
      data = iso_line_dt,
      aes(yintercept = iso_median, colour = data_type),
      linetype = "dotted", linewidth = 0.7, show.legend = FALSE
    ) +
    # gene-level line
    geom_hline(
      data = gene_dt,
      aes(yintercept = logHR_cap, linetype = gene_dir),
      colour = type_colors["gene"], linewidth = 1.0
    ) +
    # lollipop sticks
    geom_segment(
      aes(x = reorder(feature_id, logHR_cap),
          xend = reorder(feature_id, logHR_cap),
          y = 0, yend = logHR_cap,
          colour = data_type),
      linewidth = 0.4,
      alpha = 0.6
    ) +
    # isoform points
    geom_point(size = 1.8) +
    coord_flip() +
    scale_colour_manual(
      values = type_colors[c("iso_log","iso_frac")],
      name   = "Data type"
    ) +
    scale_linetype_manual(
      values = c(risk = "solid", protective = "dashed"),
      name   = "Gene-level direction"
    ) +
    labs(
      x = "Isoform (feature)",
      y = sprintf("log2(HR) (capped at ±%d)", cap_val),
      title = title
    ) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = ncol) +
    theme_bw()

  ggsave(outfile, p, width = 18, height = 12, dpi = 300)
  say("[plots] Saved %s", basename(outfile))
}

# ------------------------------------------------------------
# Plot C1: genes with heterogeneous isoform directions
# ------------------------------------------------------------
min_isoforms_mixed <- 10

candidates_C1 <- iso_summary[n_isoforms >= min_isoforms_mixed & heterogeneity == TRUE]

if (nrow(candidates_C1)) {
  # Keep only candidates that have a gene-level result
  candidates_C1 <- merge(
    candidates_C1[, .(cancer, gene_id, gene, n_isoforms)],
    gene_info_all,
    by = c("cancer", "gene_id", "gene"),
    all = FALSE
  )

  if (nrow(candidates_C1)) {
    # Keep up to 6 genes with most significant isoforms
    setorder(candidates_C1, -n_isoforms)
    candidates_C1 <- candidates_C1[seq_len(min(nrow(candidates_C1), 6))]

    # Isoform-level data (only chosen genes, only iso types)
    iso_plot_C1 <- sig_dt[
      cancer %in% candidates_C1$cancer &
        gene_id %in% candidates_C1$gene_id &
        data_type %in% c("iso_log","iso_frac")
    ]
    iso_plot_C1[, gene_label := fifelse(!is.na(gene) & gene != "", gene, gene_id)]
    iso_plot_C1[, facet_label := paste(cancer, gene_label, sep = "\n")]

    gene_info_C1 <- candidates_C1[
      , .(cancer, gene_id, gene_label = gene, logHR)
    ]
    gene_info_C1[, facet_label := paste(cancer, gene_label, sep = "\n")]
    gene_info_C1[, gene_dir := ifelse(logHR > 0, "risk", "protective")]

    make_iso_lolli(
      iso_plot_dt = iso_plot_C1,
      gene_info   = gene_info_C1,
      cap_val     = 20,
      title       = "Isoform-level survival effects (mixed directions)",
      outfile     = file.path(summary_dir, "plot_isoform_lollipop_mixed.png")
    )
  } else {
    say("[plots] No mixed-direction genes overlapping with gene-level results.")
  }
} else {
  say("[plots] No genes with >= %d isoforms and mixed directions", min_isoforms_mixed)
}


# ------------------------------------------------------------
# Plot C1: genes with heterogeneous isoform directions (per cancer)
# ------------------------------------------------------------

min_isoforms_mixed       <- 5   # softer: require at least 5 significant isoforms
max_genes_per_cancer     <- 10   # show at most 20 genes per cancer
min_iso_plotted_per_gene <- 2   # drop genes with < 2 isoforms after all filters

candidates_C1 <- iso_summary[
  n_isoforms >= min_isoforms_mixed &
    heterogeneity == TRUE
]

if (nrow(candidates_C1)) {
  # Keep only candidates that have a gene-level result
  candidates_C1 <- merge(
    candidates_C1[, .(cancer, gene_id, gene, n_isoforms)],
    gene_info_all,
    by = c("cancer", "gene_id", "gene"),
    all = FALSE
  )

  if (nrow(candidates_C1)) {
    # loop over cancers and make one plot per cancer
    cancers_C1 <- sort(unique(candidates_C1$cancer))

    for (cc in cancers_C1) {
      cand_cc <- candidates_C1[cancer == cc]

      if (!nrow(cand_cc)) next

      # within this cancer, keep the genes with most isoforms
      data.table::setorder(cand_cc, -n_isoforms)
      cand_cc <- cand_cc[seq_len(min(nrow(cand_cc), max_genes_per_cancer))]

      # Isoform-level data (only this cancer, chosen genes, iso types)
      iso_plot_cc <- sig_dt[
        cancer == cc &
          gene_id %in% cand_cc$gene_id &
          data_type %in% c("iso_log","iso_frac")
      ]
      if (!nrow(iso_plot_cc)) next

      iso_plot_cc[, gene_label := fifelse(!is.na(gene) & gene != "", gene, gene_id)]
      # per-cancer facets: just gene name (shorter)
      iso_plot_cc[, facet_label := gene_label]

      # optionally drop genes that end up with very few isoforms plotted
      keep_genes_cc <- iso_plot_cc[
        , .(n_iso_plotted = uniqueN(feature_id)),
        by = .(gene_id, gene_label)
      ][n_iso_plotted >= min_iso_plotted_per_gene]

      if (!nrow(keep_genes_cc)) next

      iso_plot_cc <- merge(
        iso_plot_cc,
        keep_genes_cc[, .(gene_id, gene_label)],
        by = c("gene_id", "gene_label")
      )

      # Gene-level direction and HR for those genes in this cancer
      gene_info_cc <- cand_cc[
        gene_id %in% keep_genes_cc$gene_id,
        .(cancer, gene_id, gene_label = gene, logHR)
      ]
      if (!nrow(gene_info_cc)) next

      gene_info_cc[, facet_label := gene_label]
      gene_info_cc[, gene_dir := ifelse(logHR > 0, "risk", "protective")]

      # safe filename from cancer id (no weird chars)
      cc_safe <- gsub("[^A-Za-z0-9]+", "_", cc)

      make_iso_lolli(
        iso_plot_dt = iso_plot_cc,
        gene_info   = gene_info_cc,
        cap_val     = 15,
        title       = sprintf("Isoform-level survival effects (mixed directions) – %s", cc),
        outfile     = file.path(summary_dir,
                                sprintf("plot_isoform_lollipop_mixed_%s.png", cc_safe))
      )
    }
  } else {
    say("[plots] No mixed-direction genes overlapping with gene-level results.")
  }
} else {
  say("[plots] No genes with >= %d isoforms and mixed directions", min_isoforms_mixed)
}


# ------------------------------------------------------------
# Plot C2: representative genes per cancer (balanced / risk / protective)
#   (only cancers with both directions and ≥ min_isoforms)
# ------------------------------------------------------------
min_isoforms <- 10

iso_gene_summary_C2 <- data.table::copy(iso_gene_summary_base)
iso_gene_summary_C2 <- iso_gene_summary_C2[
  total > 0 & n_isoforms >= min_isoforms & n_risk > 0 & n_prot > 0
]

if (nrow(iso_gene_summary_C2)) {
  iso_gene_summary_C2[, prop_risk := n_risk / total]

  # one balanced per cancer
  balanced <- iso_gene_summary_C2[
    , .SD[which.min(abs(prop_risk - 0.5))],
    by = cancer
  ]
  balanced[, category := "Balanced mixed"]

  # one risk-dominated per cancer
  risk_dom <- iso_gene_summary_C2[prop_risk > 0.5,
                                  .SD[which.max(prop_risk)],
                                  by = cancer]
  risk_dom[, category := "Risk-dominated"]

  # one protective-dominated per cancer
  prot_dom <- iso_gene_summary_C2[prop_risk < 0.5,
                                  .SD[which.min(prop_risk)],
                                  by = cancer]
  prot_dom[, category := "Protective-dominated"]

  candidates_C2 <- rbind(balanced, risk_dom, prot_dom, fill = TRUE)
  candidates_C2 <- unique(candidates_C2, by = c("cancer", "gene_id", "category"))

  # Gene-level info
  candidates_C2 <- merge(
    candidates_C2,
    gene_info_all,
    by = c("cancer", "gene_id", "gene"),
    all = FALSE
  )

  if (nrow(candidates_C2)) {
    iso_plot_C2 <- merge(
      sig_dt[data_type %in% c("iso_log", "iso_frac")],
      candidates_C2[, .(cancer, gene_id, gene, category)],
      by = c("cancer", "gene_id", "gene"),
      all = FALSE
    )
    iso_plot_C2[, gene_label := fifelse(!is.na(gene) & gene != "", gene, gene_id)]
    iso_plot_C2[, facet_label := paste(cancer, gene_label, category, sep = "\n")]

    # Optional: drop panels with very few isoforms (e.g. < 3)
    keep_genes_C2 <- iso_plot_C2[
      , .(n_iso_plotted = uniqueN(feature_id)),
      by = .(cancer, gene_id, category)
    ][n_iso_plotted >= 3]

    iso_plot_C2 <- merge(
      iso_plot_C2,
      keep_genes_C2[, .(cancer, gene_id, category)],
      by = c("cancer", "gene_id", "category")
    )

    gene_info_C2 <- candidates_C2[
      , .(cancer, gene_id, gene_label = gene, category, logHR)
    ]
    gene_info_C2 <- merge(
      gene_info_C2,
      keep_genes_C2[, .(cancer, gene_id, category)],
      by = c("cancer", "gene_id", "category"),
      all = FALSE
    )
    gene_info_C2[, facet_label := paste(cancer, gene_label, category, sep = "\n")]
    gene_info_C2[, gene_dir := ifelse(logHR > 0, "risk", "protective")]

    make_iso_lolli(
      iso_plot_dt = iso_plot_C2,
      gene_info   = gene_info_C2,
      cap_val     = 20,
      title       = "Isoform-level survival effects (representative genes)",
      outfile     = file.path(summary_dir, "plot_isoform_lollipop_extended.png")
    )
  } else {
    say("[plots] No candidates with matching gene-level signals for extended lollipops.")
  }
} else {
  say("[plots] No cancers with >= %d isoforms showing both directions", min_isoforms)
}

# ------------------------------------------------------------
# Plot C3: representative genes per cancer (for all cancers)
#   - classify genes by prop_risk using thresholds
#   - choose up to one gene per category per cancer
# ------------------------------------------------------------
iso_gene_summary_C3 <- data.table::copy(iso_gene_summary_base)
iso_gene_summary_C3 <- iso_gene_summary_C3[
  total > 0 & n_isoforms >= min_isoforms
]

if (nrow(iso_gene_summary_C3)) {
  iso_gene_summary_C3[, prop_risk := n_risk / total]

  # classify by proportion of risk isoforms
  iso_gene_summary_C3[, category := fifelse(
    prop_risk >= 2/3, "Risk-dominated",
    fifelse(prop_risk <= 1/3, "Protective-dominated", "Balanced mixed")
  )]

  # pick at most one gene per (cancer, category)
  candidates_C3 <- iso_gene_summary_C3[
    , .SD[which.max(n_isoforms)],
    by = .(cancer, category)
  ]
  candidates_C3 <- unique(candidates_C3, by = c("cancer", "gene_id", "category"))

  # gene-level info
  candidates_C3 <- merge(
    candidates_C3,
    gene_info_all,
    by = c("cancer", "gene_id", "gene"),
    all = FALSE
  )

  # Exploration mode: separate plots per cancer
  if (nrow(candidates_C3)) {
    # list of cancers to loop over
    cancers_C3 <- sort(unique(candidates_C3$cancer))

    for (cc in cancers_C3) {
      cand_cc <- candidates_C3[cancer == cc]

      iso_plot_cc <- merge(
        sig_dt[data_type %in% c("iso_log", "iso_frac")],
        cand_cc[, .(cancer, gene_id, gene, category)],
        by = c("cancer", "gene_id", "gene"),
        all = FALSE
      )
      iso_plot_cc[, gene_label := fifelse(!is.na(gene) & gene != "", gene, gene_id)]
      iso_plot_cc[, facet_label := paste(gene_label, category, sep = "\n")]

      gene_info_cc <- cand_cc[
        , .(cancer, gene_id, gene_label = gene, category, logHR)
      ]
      gene_info_cc[, facet_label := paste(gene_label, category, sep = "\n")]
      gene_info_cc[, gene_dir := ifelse(logHR > 0, "risk", "protective")]

      make_iso_lolli(
        iso_plot_dt = iso_plot_cc,
        gene_info   = gene_info_cc,
        cap_val     = 20,
        title       = sprintf("Isoform-level survival effects (%s)", cc),
        outfile     = file.path(summary_dir,
                                sprintf("plot_isoform_lollipop_%s.png", cc))
      )
    }
  }


  if (nrow(candidates_C3)) {
    iso_plot_C3 <- merge(
      sig_dt[data_type %in% c("iso_log", "iso_frac")],
      candidates_C3[, .(cancer, gene_id, gene, category)],
      by = c("cancer", "gene_id", "gene"),
      all = FALSE
    )
    iso_plot_C3[, gene_label := fifelse(!is.na(gene) & gene != "", gene, gene_id)]
    iso_plot_C3[, facet_label := paste(cancer, gene_label, category, sep = "\n")]

    gene_info_C3 <- candidates_C3[
      , .(cancer, gene_id, gene_label = gene, category, logHR)
    ]
    gene_info_C3[, facet_label := paste(cancer, gene_label, category, sep = "\n")]
    gene_info_C3[, gene_dir := ifelse(logHR > 0, "risk", "protective")]

    make_iso_lolli(
      iso_plot_dt = iso_plot_C3,
      gene_info   = gene_info_C3,
      cap_val     = 20,
      title       = "Isoform-level survival effects (representative genes across cancers)",
      outfile     = file.path(summary_dir, "plot_isoform_lollipop_all.png")
    )
  } else {
    say("[plots] No candidates with matching gene-level signals for extended lollipops.")
  }
} else {
  say("[plots] No genes with >= %d isoforms for representative lollipops", min_isoforms)
}




# ------------------------------------------------------------
# Plot C4: driver-gene isoform lollipops per cancer
#   - focus on known cancer genes per tumour type
#   - one plot per cancer (only genes with sig isoforms)
# ------------------------------------------------------------

# 1) Define known driver genes per tumour type (HGNC symbols)
driver_genes <- list(
  COAD = c("APC","TP53","KRAS","PIK3CA","SMAD4","BRAF","FBXW7",
           "CTNNB1","NRAS","TCF7L2","SOX9","ARID1A"),
  HGSC = c("TP53","BRCA1","BRCA2","RAD51C","RAD51D","PALB2",
           "ATM","BRIP1","CCNE1","RB1","NF1","PTEN","NOTCH3"),
  OV   = c("TP53","BRCA1","BRCA2","PTEN","PIK3CA","ARID1A",
           "CTNNB1","CCNE1","NF1","RB1","PPP2R1A","KRAS"),
  BRCA = c("ESR1","ERBB2","PGR","PIK3CA","PTEN","AKT1","TP53",
           "GATA3","MAP3K1","CDH1","NF1","BRCA1","BRCA2","PALB2","CHEK2"),
  UCEC = c("PTEN","PIK3CA","PIK3R1","ARID1A","CTNNB1","TP53",
           "POLE","FGFR2","MLH1","MSH2","MSH6","PMS2"),
  CCRCC = c("VHL","PBRM1","BAP1","SETD2","KDM5C","MTOR","TSC1","TSC2"),
  LUAD = c("EGFR","ALK","ROS1","MET","ERBB2","RET","NTRK1","NTRK2","NTRK3",
           "KRAS","BRAF","PIK3CA","TP53","STK11","KEAP1","NF1","RB1"),
  GBM  = c("EGFR","PDGFRA","PTEN","PIK3CA","TP53","IDH1","ATRX","NF1","BRAF","TERT"),
  HNSCC = c("TP53","CDKN2A","NOTCH1","FAT1","CASP8","FBXW7","PIK3CA","HRAS","EGFR"),
  LSCC = c("TP53","CDKN2A","PTEN","RB1","NFE2L2","KEAP1","PIK3CA",
           "SOX2","FGFR1","PDGFRA"),
  PDAC = c("KRAS","CDKN2A","TP53","SMAD4","BRCA1","BRCA2","ATM","PALB2",
           "ARID1A","TGFBR2")
)

# 2) Helper to get tumour type from 'cancer' label like "COAD:PDC000116"
get_tumour_type <- function(cancer_label) sub(":.*", "", cancer_label)

# 3) Parameters for selection
min_isoforms_driver <- 2   # min # significant isoforms per driver gene to show
max_genes_driver    <- 6   # max driver genes per cancer to include in one plot

# 4) Loop over cancers and plot driver genes with significant isoforms
cancers_all <- sort(unique(sig_dt$cancer))

for (cc in cancers_all) {
  tumour_type <- get_tumour_type(cc)

  if (!tumour_type %in% names(driver_genes)) {
    next  # no driver list defined for this type
  }

  drivers_here <- driver_genes[[tumour_type]]

  # isoform-level: significant iso_log / iso_frac for driver genes in this cancer
  iso_plot_cc <- sig_dt[
    cancer == cc &
      data_type %in% c("iso_log","iso_frac") &
      !is.na(gene) & gene %in% drivers_here
  ]

  if (!nrow(iso_plot_cc)) {
    say("[C4] No significant driver isoforms for %s (%s)", cc, tumour_type)
    next
  }

  # require at least some isoforms per gene
  gene_counts <- iso_plot_cc[
    , .(n_iso = uniqueN(feature_id)),
    by = .(gene_id, gene)
  ][n_iso >= min_isoforms_driver]

  if (!nrow(gene_counts)) {
    say("[C4] Only single-isoform driver hits for %s (%s); skipping", cc, tumour_type)
    next
  }

  # keep only those genes
  iso_plot_cc <- merge(
    iso_plot_cc,
    gene_counts[, .(gene_id, gene)],
    by = c("gene_id","gene"),
    all = FALSE
  )

  # limit to at most max_genes_driver, prioritising genes with most isoforms
  data.table::setorder(gene_counts, -n_iso)
  gene_counts <- gene_counts[seq_len(min(nrow(gene_counts), max_genes_driver))]

  iso_plot_cc <- iso_plot_cc[gene_id %in% gene_counts$gene_id]

  # build facet labels: gene symbol (and maybe gene_id if helpful)
  iso_plot_cc[, gene_label := fifelse(!is.na(gene) & gene != "", gene, gene_id)]
  iso_plot_cc[, facet_label := gene_label]

  # gene-level Cox info for these drivers in this cancer
  gene_info_cc <- long_dt[
    cancer == cc &
      data_type == "gene" &
      !is.na(gene) & gene %in% iso_plot_cc$gene_label,
    .SD[which.min(pval)],  # best p per gene
    by = .(gene_id, gene)
  ]

  if (!nrow(gene_info_cc)) {
    say("[C4] No gene-level Cox entries for driver genes in %s; skipping", cc)
    next
  }

  # align naming with iso_plot_cc
  gene_info_cc[, gene_label := fifelse(!is.na(gene) & gene != "", gene, gene_id)]
  gene_info_cc <- gene_info_cc[gene_label %in% unique(iso_plot_cc$gene_label)]

  if (!nrow(gene_info_cc)) {
    say("[C4] After alignment, no overlapping driver genes in %s; skipping", cc)
    next
  }

  gene_info_cc[, facet_label := gene_label]
  gene_info_cc[, gene_dir := ifelse(logHR > 0, "risk", "protective")]

  # clean filename
  cc_safe <- gsub("[^A-Za-z0-9]+", "_", cc)

  make_iso_lolli(
    iso_plot_dt = iso_plot_cc,
    gene_info   = gene_info_cc,
    cap_val     = 10,
    title       = sprintf("Driver-gene isoform survival effects – %s (%s)", cc, tumour_type),
    outfile     = file.path(summary_dir,
                            sprintf("plot_isoform_lollipop_drivers_%s.png", cc_safe))
  )
}


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


say("[plots] Finished Cox summary plotting for %s", summary_dir)

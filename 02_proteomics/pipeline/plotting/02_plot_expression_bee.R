#!/usr/bin/env Rscript

# ------------------------------------------------------------
# 02_plot_expression.R
#
# Plot significant model outputs across methods (cox, penalized, glmnet).
# Inputs:
#   <cox_root> <out_dir>
# Expects:
#   - Univariate Cox:      *cox_results.csv (with HR, already FDR-filtered)
#   - Penalized (coxmos):  *.penalized.csv (beta)
#   - glmnet Cox:          *.glmnet.csv (beta)
# Writes per-method outputs under <out_dir>/<method>/:
#   - expression_distribution_stats.csv
#   - expression_beeswarm_gene_coloured.png (cox only)
#   - expression_boxplot.png (penalized / glmnet)
#   - expression_median_iqr.png
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(yaml)
  library(ggbeeswarm)
  library(scales)
  library(parallel)
})

say <- function(...) message(sprintf(...))

top_n_per_dataset <- as.integer(Sys.getenv("PLOT_TOP_N", "0"))
cox_fdr_threshold <- suppressWarnings(as.numeric(Sys.getenv("PLOT_FDR_THRESHOLD", "0.05")))
if (!is.finite(cox_fdr_threshold) || cox_fdr_threshold <= 0) cox_fdr_threshold <- NA_real_
load_cores <- as.integer(Sys.getenv("PLOT_LOAD_CORES", "4"))
if (!is.finite(load_cores) || load_cores < 1) load_cores <- 1

# -------------------- args & dirs --------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 02_plot_expression.R <cox_root> <out_dir>")
}
cox_root <- args[[1]]
out_dir  <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------- meta parsing --------------------
parse_meta <- function(stem) {
  parts <- strsplit(stem, "\\.")[[1]]
  dataset <- parts[1]
  dtype <- if (length(parts) >= 2) parts[2] else "unknown"
  study <- sub("_(TMT|ITRAQ|LABELFREE).*", "", dataset)
  list(dataset = dataset, dtype = dtype, study_id = study)
}

# -------------------- colours --------------------
type_colors <- c(
  gene     = "#0072B2",
  iso_log  = "#009E73",
  iso_frac = "#E69F00",
  unknown  = "grey50"
)

# -------------------- ENST -> ENSG mapping --------------------
# tx2gene.csv is expected to have columns containing transcript (ENST*) and
# gene identifiers (ENSG*). Choose by name if available, falling back to
# positional columns.
tx2gene <- fread("02_proteomics/data/raw/tx2gene.csv")

select_col <- function(dt, candidates, fallback_idx) {
  for (nm in candidates) {
    if (nm %in% names(dt)) return(nm)
  }
  names(dt)[fallback_idx]
}

iso_col  <- select_col(tx2gene, c("isoform_id", "transcript", "tx_id", "transcript_id"), 1)
gene_col <- select_col(tx2gene, c("gene_id", "gene", "gene_name"), min(2L, ncol(tx2gene)))

tx2gene[, transcript_clean := sub("\\..*$", "", get(iso_col))]
tx2gene[, gene_clean       := sub("\\..*$", "", get(gene_col))]
tx2gene_map <- tx2gene[, .(transcript_clean, gene_id = gene_clean)]
tx2gene_lookup <- setNames(tx2gene_map$gene_id, tx2gene_map$transcript_clean)

# -------------------- loader --------------------
load_files_for_method <- function(files, method) {
  if (!length(files)) return(data.table())
  say("[info] %s: loading %d files using %d cores", method, length(files), load_cores)
  res_list <- mclapply(files, function(f) {
    dt <- fread(f)
    meta <- parse_meta(tools::file_path_sans_ext(basename(f)))

    if (method == "cox") {
      if (!"HR" %in% names(dt)) return(NULL)
      dt[, HR := as.numeric(HR)]
      dt <- dt[is.finite(HR) & HR > 0]
      if (!nrow(dt)) return(NULL)
      dt[, value := log2(HR)]
      dt[, value_label := "log2(HR)"]
    } else if (method %in% c("penalized", "glmnet")) {
      if (!"beta" %in% names(dt)) return(NULL)
      dt[, beta := as.numeric(beta)]
    dt <- dt[is.finite(beta)]
    if (!nrow(dt)) return(NULL)
    dt[, value := beta]
    dt[, value_label := "beta"]
  } else {
    return(NULL)
  }

    dt[, c("dataset","data_type","study_id","method") :=
         list(meta$dataset, meta$dtype, meta$study_id, method)]
    dt
  }, mc.cores = load_cores)
  dt <- rbindlist(res_list, fill = TRUE)
  say("[info] %s: loaded %d rows total before filtering", method, nrow(dt))
  if (!nrow(dt)) return(dt)

  if (method == "cox" && !is.na(cox_fdr_threshold) && "FDR" %in% names(dt)) {
    before <- nrow(dt)
    dt <- dt[is.finite(FDR) & FDR < cox_fdr_threshold]
    say("[info] %s: kept %d/%d rows with FDR < %.3f", method, nrow(dt), before, cox_fdr_threshold)
  }
  if (!nrow(dt)) return(dt)

  if (top_n_per_dataset <= 0 || method == "cox") return(dt)

  dt[, abs_value := abs(value)]
  if ("FDR" %in% names(dt)) {
    setorder(dt, dataset, data_type, FDR, -abs_value)
  } else {
    setorder(dt, dataset, data_type, -abs_value)
  }
  dt <- dt[, head(.SD, top_n_per_dataset), by = .(dataset, data_type)]
  say("[info] %s: retained %d rows after top-%d filtering", method, nrow(dt), top_n_per_dataset)
  dt[, abs_value := NULL]
  dt
}

# -------------------- plotting --------------------
plot_method <- function(dt, method_name, cancer_map) {
  if (!nrow(dt)) {
    say("[info] %s: no usable rows", method_name)
    return()
  }
  dt[is.na(data_type),  data_type  := "unknown"]
  dt <- dt[is.finite(value)]
  if (!nrow(dt)) {
    say("[warn] %s: all values non-finite", method_name)
    return()
  }

  # ---------- map study_id -> cancer label ----------
  dt[, cancer := {
    if (!is.null(cancer_map) && length(cancer_map)) {
      mapped <- vapply(study_id, function(s) {
        if (!is.null(cancer_map[[s]])) paste0(as.character(cancer_map[[s]]), ":", s) else s
      }, character(1))
      mapped
    } else {
      study_id
    }
  }]

  # ---------- COX: add gene_id (ENSG) for gene + isoform data ----------
  if (method_name == "cox") {
    if (!"feature" %in% names(dt)) {
      stop("Expected a 'feature' column in Cox results to derive gene_id.")
    }

    # Split into gene-level and isoform-level
    dt_gene <- copy(dt[data_type == "gene"])
    dt_iso  <- copy(dt[data_type %in% c("iso_log", "iso_frac")])

    # Gene-level: assume feature is ENSG.*; strip version
    if (nrow(dt_gene)) {
      dt_gene[, gene_id := sub("\\..*$", "", feature)]
    }

    # Isoform-level: feature = ENST.*; map via tx2gene_map (ENST -> ENSG)
    if (nrow(dt_iso)) {
      dt_iso[, transcript_clean := sub("\\..*$", "", feature)]
      dt_iso[, gene_id := tx2gene_lookup[transcript_clean]]

      # Inspect unmapped isoforms BEFORE dropping
      unmapped <- dt_iso[is.na(gene_id)]

      say("[debug] Cox isoform rows: %d, unmapped: %d",
          nrow(dt_iso), nrow(unmapped))

      if (nrow(unmapped)) {
        cat("  Example unmapped features:\n")
        print(head(unmapped$feature, 20))

        cat("\n  ID types among unmapped (first 10):\n")
        print(table(substr(unmapped$feature, 1, 4)))
      }

      # Drop rows where mapping failed
      dropped <- nrow(unmapped)
      if (dropped > 0L) {
        say("[warn] %s: dropping %d isoform rows without gene_id mapping", method_name, dropped)
        dt_iso <- dt_iso[!is.na(gene_id)]
      }
      say("[info] %s: isoform mapping complete (%d mapped)", method_name, nrow(dt_iso))
    }

    dt_other <- dt[!data_type %in% c("gene", "iso_log", "iso_frac")]
    dt <- rbindlist(list(dt_gene, dt_iso, dt_other), fill = TRUE, use.names = TRUE)

  }


  # ---------- summary stats (per dataset, data_type, cancer) ----------
  dist_stats <- dt[
    ,
    .(
      n_sig   = .N,
      med_val = median(value, na.rm = TRUE),
      iqr_val = IQR(value, na.rm = TRUE),
      min_val = min(value, na.rm = TRUE),
      max_val = max(value, na.rm = TRUE),
      value_label = first(value_label)
    ),
    by = .(dataset, study_id, data_type, cancer)
  ]

  method_out <- file.path(out_dir, method_name)
  dir.create(method_out, recursive = TRUE, showWarnings = FALSE)
  fwrite(dist_stats, file.path(method_out, "expression_distribution_stats.csv"))
  say("[write] %s/expression_distribution_stats.csv", method_name)

  # ---------- trimming ----------
  trim_lo <- quantile(dt$value, 0.01, na.rm = TRUE)
  trim_hi <- quantile(dt$value, 0.99, na.rm = TRUE)
  dt[, value_trim := pmax(pmin(value, trim_hi), trim_lo)]

  # ordering
  dataset_order <- unique(dist_stats[order(med_val)][, dataset])
  dt[, dataset := factor(dataset, levels = dataset_order)]
  dist_stats[, dataset := factor(dataset, levels = dataset_order)]

  dt[, cancer := factor(cancer)]
  cancer_order <- dt[, .(med = median(value_trim, na.rm = TRUE)), by = cancer][order(med)]$cancer
  dt[, cancer := factor(cancer, levels = cancer_order)]

  y_lab <- unique(dt$value_label)
  if (length(y_lab) != 1) y_lab <- method_name

  y_floor <- min(dt$value_trim, na.rm = TRUE) - 0.5
  y_cap   <- max(dt$value_trim, na.rm = TRUE) + 0.5

  # ============================================================
  #   COX: Beeswarm coloured by data_type
  # ============================================================
  if (method_name == "cox") {
    if (!"gene_id" %in% names(dt)) {
      stop("gene_id not present for Cox; mapping likely failed.")
    }

    p_bee_gene <- ggplot(
      dt,
      aes(
        x     = cancer,
        y     = value_trim,
        color = data_type
      )
    ) +
      ggbeeswarm::geom_quasirandom(
        aes(group = interaction(cancer, data_type)),
        dodge.width = 0.6,
        varwidth   = TRUE,
        alpha      = 0.7,
        size       = 0.7
      ) +
      geom_hline(
        yintercept = 0,
        linetype   = "dashed",
        color      = "#c90028",
        linewidth  = 0.6
      ) +
      scale_color_manual(values = type_colors, name = "Data Type") +
      labs(
        title    = sprintf("FDR-significant CoxPH per cancer (%s)", toupper(method_name)),
        subtitle = sprintf("Beeswarm by data type; trimmed to 1–99%% (%s)", y_lab),
        x        = "Cancer Type",
        y        = y_lab
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.major.y = element_line(color = "grey85"),
        legend.position    = "right"
      ) +
      coord_cartesian(ylim = c(y_floor, y_cap))

    ggsave(file.path(method_out, "expression_beeswarm_gene_coloured.png"),
           p_bee_gene, width = 12, height = 8, dpi = 300)
    say("[write] %s/expression_beeswarm_gene_coloured.png", method_name)

    dt_box <- if (!is.na(cox_fdr_threshold) && "FDR" %in% names(dt)) {
      dt[is.finite(FDR) & FDR < cox_fdr_threshold]
    } else {
      dt
    }
    fdr_label <- if (!is.na(cox_fdr_threshold)) {
      sprintf("FDR<%s", format(cox_fdr_threshold, trim = TRUE, scientific = FALSE))
    } else {
      "all rows"
    }
    if (!nrow(dt_box)) {
      say("[warn] %s: no %s rows for boxplot; skipping boxplot output",
          method_name, fdr_label)
    } else {
      p_box_cox <- ggplot(
        dt_box,
        aes(
          x = cancer,
          y = value_trim,
          color = data_type,
          group = interaction(cancer, data_type)
        )
      ) +
        geom_violin(
          aes(fill = data_type),
          position = position_dodge(width = 0.75),
          alpha = 0.25,
          color = NA,
          trim = TRUE
        ) +
        geom_boxplot(
          fill = "white",
          outlier.shape = NA,
          width = 0.4,
          position = position_dodge(width = 0.75)
        ) +
        geom_hline(
          yintercept = 0,
          linetype   = "dashed",
          color      = "#c90028",
          linewidth  = 0.6
        ) +
        geom_point(
          position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
          size = 0.8,
          alpha = 0.6,
          show.legend = FALSE
        ) +
        scale_color_manual(values = type_colors, name = "Data Type") +
        scale_fill_manual(values = type_colors, guide = "none") +
        labs(
          title    = sprintf("CoxPH distribution per cancer (%s)", fdr_label),
          subtitle = sprintf("Violin + boxplot; trimmed to 1–99%% (%s)", y_lab),
          x        = "Cancer Type",
          y        = y_lab
        ) +
        theme_bw(base_size = 12) +
        theme(
          axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
          panel.grid.major.y = element_line(color = "grey85"),
          legend.position    = "right"
        ) +
        coord_cartesian(ylim = c(y_floor, y_cap))

      ggsave(file.path(method_out, "expression_boxplot.png"),
             p_box_cox, width = 12, height = 8, dpi = 300)
      say("[write] %s/expression_boxplot.png", method_name)
    }

  } else {
    # ==========================================================
    #   penalized / glmnet: keep original grouped boxplot
    # ==========================================================
    p_box <- ggplot(
      dt,
      aes(
        x = cancer,
        y = value_trim,
        color = data_type,
        group = interaction(cancer, data_type)
      )
    ) +
      geom_boxplot(
        fill = "white",
        outlier.shape = NA,
        width = 0.65,
        position = position_dodge(width = 0.75)
      ) +
      geom_hline(
        yintercept = 0, linetype = "dashed",
        color = "#c90028", linewidth = 0.6
      ) +
      #facet_grid(~ sample_type, scales = "free_y") +
      scale_color_manual(values = type_colors, name = "Data Type") +
      labs(
        title = sprintf("Significant results per cancer (%s)", toupper(method_name)),
        subtitle = sprintf("Trimmed to 1–99%% (%s)", y_lab),
        x = "Cancer Type",
        y = y_lab
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        panel.grid.major.y = element_line(color = "grey85"),
        legend.position = "right"
      ) +
      coord_cartesian(ylim = c(y_floor, y_cap))

    ggsave(file.path(method_out, "expression_boxplot.png"),
           p_box, width = 12, height = 8, dpi = 300)
    say("[write] %s/expression_boxplot.png", method_name)
  }

  # ---------- scatter plot: median vs IQR (all methods) ----------
  p_scatter <- ggplot(dist_stats,
                      aes(x = med_val, y = iqr_val,
                          color = data_type)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = type_colors) +
    labs(
      title = sprintf("Dataset median vs IQR (%s)", toupper(method_name)),
      x     = sprintf("Median %s", y_lab),
      y     = sprintf("IQR %s", y_lab)
    ) +
    theme_bw()

  ggsave(file.path(method_out, "expression_median_iqr.png"),
         p_scatter, width = 7, height = 5, dpi = 300)
  say("[write] %s/expression_median_iqr.png", method_name)
}

# -------------------- methods & cancer map --------------------
methods <- list(
  list(name = "cox",       pattern = "cox_results\\.csv$",       exclude = "_full\\.csv$", loader = function(files) load_files_for_method(files, "cox")),
  list(name = "penalized", pattern = "\\.penalized\\.csv$",      exclude = NULL,          loader = function(files) load_files_for_method(files, "penalized")),
  list(name = "glmnet",    pattern = "\\.glmnet\\.csv$",         exclude = NULL,          loader = function(files) load_files_for_method(files, "glmnet"))
)

cancer_cfg <- tryCatch(yaml::read_yaml("02_proteomics/config/cancers.yaml"),
                       error = function(e) NULL)
cancer_map <- if (!is.null(cancer_cfg) && "cancers" %in% names(cancer_cfg)) cancer_cfg$cancers else list()

# -------------------- main loop --------------------
for (m in methods) {
  files <- list.files(cox_root, pattern = m$pattern, full.names = TRUE, recursive = TRUE)
  if (!length(files)) {
    say("[info] %s: no files matched", m$name)
    next
  }
  if (!is.null(m$exclude)) files <- files[!grepl(m$exclude, files)]
  if (!length(files)) {
    say("[info] %s: files excluded by pattern", m$name)
    next
  }
  dt <- m$loader(files)
  if (!is.null(dt) && nrow(dt)) {
    plot_method(dt, m$name, cancer_map)
  } else {
    say("[info] %s: no usable rows after loading", m$name)
  }
}

message("[done] Wrote plots and stats to ", out_dir)

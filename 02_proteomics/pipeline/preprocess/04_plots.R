#!/usr/bin/env Rscript

# ------------------------------------------------------------
# 04_plots.R
#
# For each dataset (study_platform_type):
#   1) Boxplot grid (3 x 3)
#      rows  : gene, iso_log, iso_frac
#      cols  : RAW, NORMALIZED, BATCH CORRECTED
#   2) Metrics grid
#      rows  : Batch variance (R^2), Silhouette
#      cols  : data types with available metrics
#
# RAW      : cyclic-loess normalised, not batch corrected
# NORMAL   : log2(TPM+1) + cyclic-loess, not batch corrected
# BATCH    : normalised + ComBat batch corrected
#
# iso_frac is always derived from iso_log at each stage.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(cluster)
})

# ---------------- command line args ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 || length(args) > 6) {
  stop("Usage: 04_plots.R <raw_root> <norm_root> <batch_root> <batch_annotation_dir> <out_dir> [dataset_id]")
}
raw_root   <- args[[1]]
norm_root  <- args[[2]]
batch_root <- args[[3]]
batch_dir  <- args[[4]]
out_dir    <- args[[5]]
target_dataset <- if (length(args) == 6 && nzchar(args[[6]])) args[[6]] else NA_character_

# ---------------- utility functions ----------------
say <- function(fmt, ...) cat(sprintf(paste0("[plot] ", fmt, "\n"), ...))
max_plot_features <- as.integer(Sys.getenv("PREPROCESS_PLOT_MAX_FEATURES", "4000"))

# ---------------- sample type inference ----------------
infer_sample_type <- function(dataset) {
  if (grepl("_reference$", dataset, ignore.case = TRUE)) return("reference")
  "tumor"
}

# ---------------- batch annotation loading ----------------
canonicalize_batch <- function(vals) {
  if (is.null(vals)) return(rep(NA_character_, length(vals)))
  vapply(vals, function(v) {
    if (is.null(v) || is.na(v) || !nzchar(v)) return(NA_character_)
    trimws(v)
  }, character(1), USE.NAMES = FALSE)
}

# ---------------- data loading and processing ----------------
load_case_batches <- function(study_key) {
  if (!nzchar(batch_dir)) return(NULL)
  file <- file.path(batch_dir, paste0(study_key, ".csv"))
  if (!file.exists(file)) {
    cand <- list.files(
      batch_dir,
      pattern = paste0("^", study_key, "\\.csv$"),
      ignore.case = TRUE,
      full.names = TRUE
    )
    if (!length(cand)) return(NULL)
    file <- cand[1]
    say("Using batch annotation file %s for study key %s",
        basename(file), study_key)
  }
  dt <- tryCatch(fread(file), error = function(e) NULL)
  if (is.null(dt) || !all(c("case_id", "folder_name") %in% names(dt))) return(NULL)
  dt[, case_id := toupper(trimws(case_id))]
  if ("sample_type" %in% names(dt)) {
    dt[, sample_type := trimws(tolower(sample_type))]
  } else {
    dt[, sample_type := NA_character_]
  }
  dt[, folder_name := trimws(folder_name)]
  dt <- dt[nzchar(case_id) & nzchar(folder_name)]
  if (!nrow(dt)) return(NULL)
  dt <- dt[order(case_id, sample_type, folder_name)]
  dt <- dt[!duplicated(data.table(case_id, sample_type))]
  dt
}

# ---------------- load expression matrix ----------------
load_matrix <- function(path) {
  if (!file.exists(path)) return(NULL)
  dt <- fread(path)
  if (ncol(dt) < 2) return(NULL)
  feats <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- feats
  storage.mode(mat) <- "double"
  keep <- !grepl("POOLED|QC", colnames(mat), ignore.case = TRUE)
  mat[, keep, drop = FALSE]
}

# ---------------- trim batch labels ----------------
trim_batch <- function(b, width = 7) {
  if (is.null(b)) return(NA_character_)
  b_chr <- as.character(b)
  b_chr[!nzchar(b_chr)] <- NA_character_
  substr(b_chr, 1, width)
}

# ---------------- clip matrix values ----------------
clip_matrix <- function(mat, probs = c(0.1, 0.9)) {
  if (is.null(mat)) return(mat)
  qs <- quantile(as.numeric(mat), probs = probs, na.rm = TRUE, names = FALSE)
  if (all(is.finite(qs))) mat <- pmin(pmax(mat, qs[1]), qs[2])
  mat
}

# ---------------- batch palette ----------------
get_batch_palette <- function(batch_levels) {
  n <- length(batch_levels)
  scales::hue_pal(h = c(0, 360), c = 80, l = 60)(n)
}

# ---------------- gather values for boxplots ----------------
gather_values <- function(mat, stage_label, order_idx, ordered_batches,
                          dataset_name, entries_list) {
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) return(entries_list)
  mat_use <- mat
  if (nrow(mat_use) > max_plot_features) {
    idx <- sample.int(nrow(mat_use), max_plot_features)
    mat_use <- mat_use[idx, , drop = FALSE]
  }
  mat_use <- mat_use[, order_idx, drop = FALSE]
  values        <- as.vector(mat_use)
  sample_order  <- rep(seq_along(order_idx), each = nrow(mat_use))
  batches       <- rep(ordered_batches,      each = nrow(mat_use))
  entries_list[[length(entries_list) + 1L]] <- data.table(
    dataset      = dataset_name,
    stage        = stage_label,
    value        = values,
    sample_order = sample_order,
    batch        = batches
  )
  entries_list
}

# ---------------- PCA and metrics ----------------
do_pca_with_metrics <- function(mat, batch_vec, stage_label,
                                max_features = 4000, n_pcs = 5) {
  valid_batches <- batch_vec[!is.na(batch_vec)]
  if (length(unique(valid_batches)) < 2 || length(valid_batches) < 2) return(NULL)
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) return(NULL)

  mat_use <- mat
  if (nrow(mat_use) > max_features) {
    idx <- sample.int(nrow(mat_use), max_features)
    mat_use <- mat_use[idx, , drop = FALSE]
  }
  row_ok <- apply(mat_use, 1, function(x) any(is.finite(x)))
  mat_use <- mat_use[row_ok, , drop = FALSE]
  if (!nrow(mat_use) || ncol(mat_use) < 2) return(NULL)
  mat_use[!is.finite(mat_use)] <- 0

  pcs <- tryCatch(prcomp(t(mat_use), scale. = TRUE, center = TRUE),
                  error = function(e) NULL)
  if (is.null(pcs) || ncol(pcs$x) < 2) return(NULL)

  scores_dt <- data.table(
    PC1   = pcs$x[, 1],
    PC2   = pcs$x[, 2],
    batch = batch_vec,
    stage = stage_label,
    sample_id = rownames(pcs$x)
  )

  batch_fac <- factor(batch_vec)
  max_pc <- min(n_pcs, ncol(pcs$x))
  r2_vals <- numeric(max_pc)
  for (k in seq_len(max_pc)) {
    r2_vals[k] <- summary(lm(pcs$x[, k] ~ batch_fac))$r.squared
  }
  r2_dt <- data.table(stage = stage_label, PC = seq_len(max_pc), R2 = r2_vals)

  sil <- tryCatch(
    silhouette(as.integer(batch_fac),
               dist(pcs$x[, seq_len(max_pc), drop = FALSE])),
    error = function(e) NULL
  )
  if (is.null(sil)) return(NULL)
  sil_widths <- tryCatch(as.numeric(sil[, "sil_width"]), error = function(e) NA_real_)
  mean_sil <- if (length(sil_widths)) mean(sil_widths, na.rm = TRUE) else NA_real_
  sil_dt <- data.table(
    stage    = stage_label,
    mean_sil = mean_sil,
    n_samples = nrow(pcs$x)
  )
  list(scores = scores_dt, r2 = r2_dt, sil = sil_dt)
}

# ---------------- Boxplots ----------------
make_boxplot <- function(box_dt, stage_label, stage_title, dtype_label,
                         pal, show_row_label = FALSE,
                         y_limits = NULL, y_lab = "Intensity") {

  ggplot(box_dt[stage == stage_label],
         aes(x = sample_order, y = value_plot_clipped,
             fill = batch, color = batch)) +
    geom_boxplot(
      colour        = "#8B8682",
      linewidth     = 0.2,
      fatten        = 1.6,
      outlier.shape = NA,
      outlier.alpha = 0.15,
      width         = 0.55,
      alpha         = 0.90
    ) +
    stat_summary(
      fun.data = function(y) {
        m <- median(y, na.rm = TRUE)
        data.frame(y = m, ymin = m, ymax = m)
      },
      geom  = "crossbar",
      width = 0.5,
      color = "black",
      linewidth = 0.35
    ) +
    labs(
      title = if (show_row_label) dtype_label else NULL,
      subtitle = NULL,
      x = "Samples",
      y = y_lab
    ) +
    scale_fill_manual(values = pal, drop = FALSE) +
    coord_cartesian(ylim = y_limits) +
    theme_bw(base_size = 11) +
    theme(
      panel.border       = element_rect(colour = "#8B8682", linewidth = 0.3),
      axis.line          = element_line(colour = "#8B8682", linewidth = 0.2),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(linewidth = 0.15, colour = "#8B8682"),
      plot.title         = element_text(face = "italic", hjust = 0, size = 10),
      plot.subtitle      = element_text(face = "italic", hjust = 0, size = 10),
      axis.title.x       = element_text(margin = margin(t = 6)),
      axis.title.y       = element_text(margin = margin(r = 6)),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      legend.position    = "none"
    )
}

# ---------------- adorn axes ----------------
adorn_axes <- function(p, show_x = FALSE, show_y = FALSE) {
  if (inherits(p, "ggplot")) {
    p <- p +
      theme(
        axis.title.x = if (show_x) element_text(margin = margin(t = 6)) else element_blank(),
        axis.title.y = if (show_y) element_text(margin = margin(r = 6)) else element_blank()
      )
  }
  p
}

# ---------------- PCA plot ----------------
make_pca_plot <- function(dt, stage_label, panel_label, panel_title, pal) {
  ggplot(dt[stage == stage_label], aes(x = PC1, y = PC2, color = batch)) +
    geom_point(alpha = 0.95, size = 2) +
    labs(
      title  = sprintf("(%s) %s", panel_label, panel_title),
      x = "PC1", y = "PC2", color = "Batch"
    ) +
    scale_colour_manual(values = pal, drop = FALSE) +
    stat_ellipse(level = 0.68, linewidth = 0.3, alpha = 0.6) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0),
      legend.position  = "none",
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.3, "cm")
    )
}

# ---------------- Dataset list / output dirs ----------------
gene_dir <- file.path(raw_root, "gene")
if (!dir.exists(gene_dir)) stop("Missing raw gene directory: ", gene_dir)
datasets <- sub("_gene\\.csv$", "", list.files(gene_dir,
                                               pattern = "_gene\\.csv$",
                                               full.names = FALSE))
datasets <- datasets[!grepl("_reference$", datasets, ignore.case = TRUE)]
if (!is.na(target_dataset)) {
  datasets <- datasets[startsWith(datasets, target_dataset)]
}
if (!length(datasets)) stop("No datasets found for plotting.")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
box_dir <- file.path(out_dir, "boxplots")
met_dir <- file.path(out_dir, "metrics")
dir.create(box_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(met_dir, recursive = TRUE, showWarnings = FALSE)

plot_dtypes <- c("gene", "iso_log", "iso_frac") # for boxplots
metric_dtypes <- c("gene", "iso_log")           # for metrics

# ---------------- iso_frac mapping ----------------
map_file <- Sys.getenv("PROT_ENST_ENSG_MAP",
                       "02_proteomics/data/raw/ENST-ENSG_mapping.csv")

load_mapping <- function(path) {
  if (!file.exists(path)) stop("[plot] ENST-ENSG mapping not found: ", path)
  dt <- fread(path, header = FALSE,
              col.names = c("transcript_id", "gene_id"))
  dt[, transcript_id := trimws(transcript_id)]
  dt[, gene_id      := trimws(gene_id)]
  dt <- dt[nzchar(transcript_id) & nzchar(gene_id)]
  dt
}
mapping_dt <- load_mapping(map_file)

iso_frac_from_iso_log <- function(mat, mapping_dt) {
  if (is.null(mat) || !nrow(mat) || !ncol(mat)) return(NULL)
  tx <- rownames(mat)
  gene_map <- mapping_dt$gene_id[match(tx, mapping_dt$transcript_id)]
  keep <- !is.na(gene_map)
  if (!any(keep)) return(NULL)
  mat <- mat[keep, , drop = FALSE]
  gene_map <- gene_map[keep]
  res <- mat
  genes <- unique(gene_map)
  for (g in genes) {
    idx  <- which(gene_map == g)
    sub  <- mat[idx, , drop = FALSE]
    denom <- colSums(sub, na.rm = TRUE)
    frac  <- sweep(sub, 2, denom, "/")
    frac[, denom == 0] <- NA_real_
    res[idx, ] <- frac
  }
  res
}

# ---------------- data loading per dtype ----------------
load_mats_for_dtype <- function(dataset, dtype) {
  if (dtype == "iso_frac") {
    iso_fname <- paste0(dataset, "_iso_log.csv")
    iso_raw   <- load_matrix(file.path(raw_root,   "iso_log", iso_fname))
    iso_norm  <- load_matrix(file.path(norm_root,  "iso_log", iso_fname))
    iso_batch <- load_matrix(file.path(batch_root, "iso_log", iso_fname))

    raw_mat   <- if (!is.null(iso_raw))
      iso_frac_from_iso_log(iso_raw, mapping_dt) else NULL
    norm_mat  <- if (!is.null(iso_norm))
      iso_frac_from_iso_log(pmax(2^iso_norm - 1, 0), mapping_dt) else NULL
    batch_mat <- if (!is.null(iso_batch))
      iso_frac_from_iso_log(pmax(2^iso_batch - 1, 0), mapping_dt) else NULL
  } else {
    fname <- paste0(dataset, "_", dtype, ".csv")
    raw_mat   <- load_matrix(file.path(raw_root,   dtype, fname))
    norm_mat  <- load_matrix(file.path(norm_root,  dtype, fname))
    batch_mat <- load_matrix(file.path(batch_root, dtype, fname))
  }
  list(raw = raw_mat, norm = norm_mat, batch = batch_mat)
}

# ---------------- main per-dataset loop ----------------
for (dataset in sort(datasets)) {
  say("Processing dataset %s", dataset)

  study_key    <- sub("_(normal|tumor|reference)$", "", dataset, ignore.case = TRUE)
  dataset_type <- infer_sample_type(dataset)
  case_batches <- load_case_batches(study_key)

  panels_by_dtype <- list()
  box_panels_by_dtype <- list()

  for (dtype in plot_dtypes) {

    mats <- load_mats_for_dtype(dataset, dtype)
    available_mats <- Filter(Negate(is.null), mats)
    if (!length(available_mats)) {
      say("Skipping %s (%s) – missing matrices", dataset, dtype)
      next
    }

    common_cols <- Reduce(intersect, lapply(available_mats, colnames))
    if (!length(common_cols)) {
      say("No overlapping samples for %s (%s)", dataset, dtype)
      next
    }

    raw_mat   <- if (!is.null(mats$raw))   mats$raw[,   common_cols, drop = FALSE] else NULL
    norm_mat  <- if (!is.null(mats$norm))  mats$norm[,  common_cols, drop = FALSE] else NULL
    batch_mat <- if (!is.null(mats$batch)) mats$batch[, common_cols, drop = FALSE] else NULL

    sample_ids <- common_cols
    sample_batches <- rep("unknown", length(sample_ids))

    if (!is.null(case_batches)) {
      dt_sub <- case_batches
      if (!is.na(dataset_type) && "sample_type" %in% names(case_batches)) {
        sel <- dt_sub[sample_type == dataset_type]
        if (nrow(sel)) dt_sub <- sel
      }
      idx <- match(toupper(sample_ids), toupper(dt_sub$case_id))
      batch_vals <- dt_sub$folder_name[idx]
      sample_batches[!is.na(batch_vals)] <- batch_vals[!is.na(batch_vals)]
    }

    order_idx <- order(sample_batches, sample_ids)
    ordered_batches <- trim_batch(sample_batches[order_idx])

    # ---------------- boxplots ----------------
    entries <- list()
    entries <- gather_values(raw_mat,   "raw",             order_idx, ordered_batches, dataset, entries)
    entries <- gather_values(norm_mat,  "normalized",      order_idx, ordered_batches, dataset, entries)
    entries <- gather_values(batch_mat, "batch_corrected", order_idx, ordered_batches, dataset, entries)

    box_dt <- rbindlist(entries, fill = TRUE)
    box_dt[, stage := factor(stage,
                             levels = c("raw", "normalized", "batch_corrected"))]

    batch_levels <- sort(unique(box_dt$batch))
    box_dt[, batch := factor(batch, levels = batch_levels)]

    box_dt[, value_plot := value]

    box_dt[, value_plot_clipped := value_plot]
    box_dt[, sample_order := factor(sample_order)]

    q_global <- quantile(box_dt$value_plot,
                         probs = c(0.05, 0.95),
                         na.rm = TRUE,
                         names = FALSE)
    if (all(is.finite(q_global))) {
      box_dt[!is.na(value_plot_clipped),
             value_plot_clipped := pmin(pmax(value_plot_clipped,
                                             q_global[1]),
                                        q_global[2])]
    }

    y_limits <- range(box_dt$value_plot_clipped, na.rm = TRUE)
    if (!all(is.finite(y_limits))) y_limits <- NULL
    y_lab <- if (dtype == "iso_frac") "Fraction" else "Intensity"

    batch_palette <- get_batch_palette(batch_levels)

    has_raw  <- any(box_dt$stage == "raw")
    has_norm <- any(box_dt$stage == "normalized")
    has_batch<- any(box_dt$stage == "batch_corrected")

    p_raw   <- if (has_raw)
      make_boxplot(box_dt, "raw", "", dtype, batch_palette,
                   show_row_label = TRUE,
                   y_limits = y_limits, y_lab = y_lab) else nullGrob()

    p_norm  <- if (has_norm)
      make_boxplot(box_dt, "normalized", "", dtype, batch_palette,
                   show_row_label = FALSE,
                   y_limits = y_limits, y_lab = y_lab) else nullGrob()

    p_batch <- if (has_batch)
      make_boxplot(box_dt, "batch_corrected", "", dtype, batch_palette,
                   show_row_label = FALSE,
                   y_limits = y_limits, y_lab = y_lab) else nullGrob()

    # ---------------- Batch metrics via PCA (only for gene / iso_log) ----------------
    p_r2  <- nullGrob()
    p_sil <- nullGrob()

    if (dtype %in% c("gene", "iso_log")) {   # <--- only these get PCA metrics
      res_norm  <- if (has_norm)
        do_pca_with_metrics(norm_mat[, order_idx, drop = FALSE],
                            ordered_batches, "normalized", max_plot_features)
      else NULL

      res_batch <- if (has_batch)
        do_pca_with_metrics(batch_mat[, order_idx, drop = FALSE],
                            ordered_batches, "batch_corrected", max_plot_features)
      else NULL

      if (!is.null(res_norm) && !is.null(res_batch)) {
        r2_dt  <- rbind(res_norm$r2,  res_batch$r2)
        sil_dt <- rbind(res_norm$sil, res_batch$sil)

        r2_dt[, stage := factor(stage,
                                levels = c("normalized", "batch_corrected"), 
                                labels = c("Norm", "Norm+Batch"))]
        p_r2 <- ggplot(r2_dt, aes(x = factor(PC), y = R2, fill = stage)) +
          geom_col(position = position_dodge(width = 0.7), width = 0.6) +
          geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.2) +
          scale_fill_brewer(palette = "Set2") +
          labs(title = "Batch variance in PCs (R²)",
               x = "Principal component",
               y = expression(R^2*" (batch)"),
               fill = "Stage") +
          theme_bw(base_size = 10) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0))

        sil_dt[, stage := factor(stage,
                                 levels = c("normalized", "batch_corrected"),
                                 labels = c("Norm", "Norm+Batch"))]
        p_sil <- ggplot(sil_dt, aes(x = stage, y = mean_sil, fill = stage)) +
          geom_col(width = 0.5) +
          geom_text(aes(label = sprintf("%.2f", mean_sil)), vjust = -0.2, size = 3) +
          scale_fill_brewer(palette = "Set2", guide = "none") +
          ylim(0, 1) +
          labs(title = "Batch clustering (silhouette)",
               x = NULL, y = "Silhouette width") +
          theme_bw(base_size = 10) +
          theme(panel.grid.minor = element_blank(),
                plot.title = element_text(face = "bold", hjust = 0))
      }
    }

    box_panels_by_dtype[[dtype]] <- list(p_raw, p_norm, p_batch)
    panels_by_dtype[[dtype]]     <- list(p_r2, p_sil)
  }

  if (!length(box_panels_by_dtype)) {
    say("No panels generated for %s", dataset)
    next
  }

  # ---------------- assemble boxplot grid ----------------
  dtype_order   <- intersect(plot_dtypes, names(box_panels_by_dtype))
  stage_headers <- c("RAW", "NORMALIZED", "BATCH CORRECTED")
  grobs_box <- list()

  # header row: four stage headers
  for (s in stage_headers) {
    grobs_box[[length(grobs_box) + 1]] <- textGrob(
      s,
      gp = gpar(fontface = "bold"),
      y  = unit(0.3, "npc")
    )
  }

  # body rows: panels only (row labels are titles of RAW panels)
  for (d in dtype_order) {
    panels <- box_panels_by_dtype[[d]]
    for (s_idx in seq_along(stage_headers)) {
      panel <- if (length(panels) >= s_idx && !is.null(panels[[s_idx]])) {
        adorn_axes(panels[[s_idx]], show_x = FALSE, show_y = FALSE)
      } else {
        nullGrob()
      }
      grobs_box[[length(grobs_box) + 1]] <- panel
    }
  }

  box_layout_mat <- matrix(
    seq_along(grobs_box),
    ncol = length(stage_headers),
    byrow = TRUE
  )

  box_heights <- c(unit(0.3, "lines"),
                   rep(unit(1, "null"), length(dtype_order)))

  box_body <- grid.arrange(
    grobs = grobs_box,
    layout_matrix = box_layout_mat,
    heights = box_heights
  )

  box_final_plot <- arrangeGrob(
    box_body,
    left = textGrob(
      "Intensity",
      rot = 90,
      gp  = gpar(fontface = "bold"),
      x   = unit(1.5, "lines")
    ),
    bottom = textGrob(
      "Samples",
      gp  = gpar(fontface = "bold"),
      y   = unit(1.2, "lines")
    ),
    padding = unit(c(1, 0.5, 1, 0.5), "lines")
  )

  # ----------------- Metrics layout -----------------
  metric_titles <- c("Batch variance (R²)", "Silhouette")
  dtype_metrics <- intersect(metric_dtypes, names(panels_by_dtype))

  if (length(dtype_metrics)) {

    grobs_met <- list()

    # Header row: one header per dtype (gene, iso_log)
    for (d in dtype_metrics) {
      grobs_met[[length(grobs_met) + 1]] <- textGrob(
        d,
        gp = gpar(fontface = "bold"),
        y  = unit(0.3, "npc")
      )
    }

    # Body rows: one row per metric, columns = dtypes
    for (m_idx in seq_along(metric_titles)) {
      for (d in dtype_metrics) {
        panels <- panels_by_dtype[[d]]
        panel  <- if (!is.null(panels) &&
                      length(panels) >= m_idx &&
                      !is.null(panels[[m_idx]])) {
          panels[[m_idx]]
        } else {
          nullGrob()
        }
        grobs_met[[length(grobs_met) + 1]] <- panel
      }
    }

    # Layout: ncol = number of dtypes, first row is header
    met_layout <- matrix(
      seq_along(grobs_met),
      ncol = length(dtype_metrics),
      byrow = TRUE
    )

    met_heights <- c(
      unit(0.5, "lines"),                           # header row (thin)
      rep(unit(1, "null"), length(metric_titles))   # one for each metric row
    )

    metric_body <- grid.arrange(
      grobs         = grobs_met,
      layout_matrix = met_layout,
      heights       = met_heights
    )

    final_plot <- arrangeGrob(
      metric_body,
      left = textGrob(
        "",
        rot = 90,
        gp  = gpar(fontface = "bold"),
        x   = unit(1.5, "lines")
      ),
      bottom = textGrob(
        "",
        gp  = gpar(fontface = "bold"),
        y   = unit(1.2, "lines")
      ),
      padding = unit(c(1, 0.5, 1, 0.5), "lines")
    )

  } else {
    final_plot <- nullGrob()
  }


  # ---------------- save plots ----------------
  out_png <- file.path(box_dir,
                       sprintf("%s_preprocess_boxplots.png", dataset))
  ggsave(out_png, box_final_plot,
         width = 20, height = 8, units = "in", dpi = 300)

  out_png <- file.path(met_dir,
                       sprintf("%s_metrics.png", dataset))
  ggsave(out_png, final_plot,
         width = 10, height = 6, units = "in", dpi = 300)

  say("Saved plot for %s", dataset)
}
say("All done.")

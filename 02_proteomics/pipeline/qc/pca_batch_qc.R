#!/usr/bin/env Rscript

# ------------------------------------------------------------
# PCA & clustering QC plots before/after batch correction
# ------------------------------------------------------------
# This script inspects the matrices under 02_proteomics/out/preprocessed
# and generates PCA scatter plots plus hierarchical clustering trees for
# each data type (gene / iso_log / iso_frac) across the following stages:
#   - Normalized (before batch correction)
#   - Batch-corrected (after correction)
#
# Output directory can be overridden via PROT_PCA_QC_OUT.
# Batch annotation directory can be set via PROT_BATCH_ANNOTATION_DIR.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggdendro)
  library(gridExtra)
  library(yaml)
})

say <- function(...) message(sprintf(...))
`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a)) a else b
normalize_sample_id <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- trimws(x)
  toupper(x)
}
study_from_filename <- function(name) {
  sub("_(gene|iso_log|iso_frac)$", "", name, ignore.case = TRUE)
}

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

parse_data_types <- function(raw) {
  if (!nchar(raw)) return(character())
  toks <- trimws(unlist(strsplit(raw, "[,;]", fixed = FALSE)))
  toks <- toks[nzchar(toks)]
  unique(toks)
}

DATA_TYPES <- {
  env_val <- Sys.getenv("PROT_PCA_QC_DATA_TYPES", "")
  parsed <- parse_data_types(env_val)
  if (length(parsed)) parsed else c("gene", "iso_log")
}
STAGE_DIRS <- list(
  RAW   = "02_proteomics/out/preprocessed/normalization",
  BATCH = "02_proteomics/out/preprocessed/batch_corrected"
)
STAGE_LABEL <- c(
  RAW   = "Pre-batch",
  BATCH = "Batch-corrected"
)

OUT_ROOT <- Sys.getenv(
  "PROT_PCA_QC_OUT",
  "02_proteomics/out/preprocessed/plots/pca_batch_qc"
)
dir.create(OUT_ROOT, recursive = TRUE, showWarnings = FALSE)

BATCH_DIR <- Sys.getenv(
  "PROT_BATCH_ANNOTATION_DIR",
  "02_proteomics/data/batch_annotation"
)

SHOW_LEGEND <- {
  val <- tolower(Sys.getenv("PROT_PCA_QC_SHOW_LEGEND", "true"))
  val %in% c("1", "true", "yes", "y", "on")
}

# ------------------------------------------------------------
# Metadata helpers
# ------------------------------------------------------------

load_cancer_map <- function() {
  cfg_path <- "02_proteomics/config/cancers.yaml"
  if (!file.exists(cfg_path)) return(list())
  cfg <- tryCatch(yaml::read_yaml(cfg_path), error = function(e) NULL)
  if (is.null(cfg) || is.null(cfg$cancers)) return(list())
  cfg$cancers
}

build_metadata <- function(batch_dir) {
  if (!dir.exists(batch_dir)) {
    stop("Batch annotation directory not found: ", batch_dir)
  }
  files <- list.files(batch_dir, pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) {
    stop("No batch annotation CSV files found under: ", batch_dir)
  }
  cancer_map <- load_cancer_map()
  rows <- lapply(files, function(f) {
    dt <- tryCatch(fread(f), error = function(e) NULL)
    if (is.null(dt) || !nrow(dt)) return(NULL)
    col_names <- names(dt)
    lower_names <- tolower(col_names)
    case_idx <- which(lower_names == "case_id")
    if (!length(case_idx)) {
      case_idx <- which(lower_names == "sample_id")
    }
    if (!length(case_idx)) return(NULL)
    case_col <- col_names[case_idx[1]]
    folder_idx <- which(lower_names == "folder_name")
    if (!length(folder_idx)) return(NULL)
    folder_col <- col_names[folder_idx[1]]
    dt[, sample_id := normalize_sample_id(get(case_col))]
    dt[, batch := trimws(as.character(get(folder_col)))]
    dt <- dt[nzchar(sample_id) & nzchar(batch)]
    if (!nrow(dt)) return(NULL)
    dt <- dt[, .SD[1], by = sample_id]
    study_name <- tools::file_path_sans_ext(basename(f))
    study_key <- toupper(sub("_reference$", "", study_name))
    cancer <- if (study_key %in% names(cancer_map)) cancer_map[[study_key]] else study_key
    dt[, `:=`(study = study_key,
              cancer = as.character(cancer))]
    dt[, .(sample_id, cancer, study, batch)]
  })
  rows <- Filter(function(x) !is.null(x) && nrow(x), rows)
  if (!length(rows)) stop("Unable to assemble metadata from batch annotations.")
  meta <- rbindlist(rows, fill = TRUE)
  meta <- unique(meta, by = c("study","sample_id"))
  required <- c("sample_id", "cancer", "study", "batch")
  missing <- setdiff(required, names(meta))
  if (length(missing)) {
    stop("Metadata missing required columns: ", paste(missing, collapse = ", "))
  }
  meta
}

META <- build_metadata(BATCH_DIR)
say("[META] Loaded %d samples spanning %d studies", nrow(META), uniqueN(META$study))
batch_levels <- sort(unique(na.omit(META$batch)))
distinct_palette <- function(n) {
  if (n <= 0) return(character())
  if (requireNamespace("colorspace", quietly = TRUE)) {
    # use evenly spaced hues in HCL space for arbitrarily large n
    hues <- seq(0, 360, length.out = n + 1)[- (n + 1)]
    return(colorspace::hex(colorspace::polarLUV(L = 65, C = 100, H = hues)))
  }
  hues <- seq(0, 360, length.out = n + 1)[- (n + 1)]
  grDevices::hcl(h = hues, c = 100, l = 65)
}
batch_palette <- if (length(batch_levels)) {
  setNames(distinct_palette(length(batch_levels)), batch_levels)
} else {
  character()
}

batch_colors_for <- function(levels) {
  lev <- unique(na.omit(levels))
  if (!length(lev)) return(character())
  cols <- batch_palette[lev]
  missing <- is.na(cols)
  if (any(missing)) {
    cols[missing] <- grDevices::rainbow(sum(missing))
    names(cols)[missing] <- lev[missing]
  }
  cols
}

# ------------------------------------------------------------
# Matrix helpers
# ------------------------------------------------------------

read_matrix <- function(path) {
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || ncol(dt) < 2) return(NULL)
  mat <- as.matrix(dt[, -1, with = FALSE])
  colnames(mat) <- normalize_sample_id(colnames(mat))
  rownames(mat) <- dt[[1]]
  mat
}

run_pca <- function(mat, meta, study_id, dataset_name = NULL) {
  meta_study <- meta[study == study_id]
  if (!nrow(meta_study)) {
    say("[SKIP] %s (%s): no metadata rows", dataset_name %||% "<unknown>", study_id)
    return(NULL)
  }
  samples <- colnames(mat)
  samples <- samples[samples %in% meta_study$sample_id]
  if (length(samples) < 3) {
    say("[SKIP] %s (%s): matched %d samples < 3", dataset_name %||% "<unknown>", study_id, length(samples))
    return(NULL)
  }
  mat_sub <- mat[, samples, drop = FALSE]
  mat_sub[!is.finite(mat_sub)] <- NA_real_
  col_sd <- apply(mat_sub, 2, sd, na.rm = TRUE)
  keep_cols <- which(is.finite(col_sd) & col_sd > 0)
  if (length(keep_cols) < 3) {
    say("[SKIP] %s (%s): <3 informative samples after variance filter", dataset_name %||% "<unknown>", study_id)
    return(NULL)
  }
  mat_sub <- mat_sub[, keep_cols, drop = FALSE]
  samples <- samples[keep_cols]
  row_sd <- apply(mat_sub, 1, sd, na.rm = TRUE)
  keep_rows <- which(is.finite(row_sd) & row_sd > 0)
  if (length(keep_rows) < 3) {
    say("[SKIP] %s (%s): <3 informative features after variance filter", dataset_name %||% "<unknown>", study_id)
    return(NULL)
  }
  mat_sub <- mat_sub[keep_rows, , drop = FALSE]
  if (anyNA(mat_sub)) {
    row_ok <- apply(mat_sub, 1, function(x) all(is.finite(x)))
    mat_sub <- mat_sub[row_ok, , drop = FALSE]
    if (nrow(mat_sub) < 3) {
      say("[SKIP] %s (%s): <3 complete features after NA removal", dataset_name %||% "<unknown>", study_id)
      return(NULL)
    }
    col_ok <- apply(mat_sub, 2, function(x) all(is.finite(x)))
    mat_sub <- mat_sub[, col_ok, drop = FALSE]
    samples <- samples[col_ok]
    if (ncol(mat_sub) < 3) {
      say("[SKIP] %s (%s): <3 complete samples after NA removal", dataset_name %||% "<unknown>", study_id)
      return(NULL)
    }
  }
  mat_t <- t(mat_sub)
  meta_sub <- meta_study[match(samples, meta_study$sample_id)]
  pca <- prcomp(mat_t, center = TRUE, scale. = TRUE)
  data.table(
    sample_id = samples,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    cancer = meta_sub$cancer,
    study = study_id,
    batch = meta_sub$batch
  )
}

build_pca_plot <- function(dt_study, dtype, stage, stage_lab, study_id) {
  dt_plot <- copy(dt_study)
  dt_plot[, batch := as.character(batch)]
  dt_plot[!nzchar(batch) | is.na(batch), batch := "Unknown"]
  cols <- batch_colors_for(dt_plot$batch)
  if (!length(cols)) {
    cols <- setNames("#1f78b4", unique(dt_plot$batch))
  }
  p <- ggplot(dt_plot, aes(x = PC1, y = PC2, colour = batch)) +
    geom_point(size = 1.5, alpha = 1) +
    labs(
      title = sprintf("PCA – %s (%s) – %s", dtype, stage_lab, study_id),
      subtitle = sprintf("Batch colouring (%d samples)", nrow(dt_study)),
      x = "PC1",
      y = "PC2",
      colour = "Batch"
    ) +
    theme_bw()
  unique_cancers <- unique(na.omit(dt_plot$cancer))
  if (length(unique_cancers) > 1) {
    p <- p + aes(shape = cancer) + labs(shape = "Cancer")
  }
  p <- p + scale_colour_manual(values = cols, drop = FALSE)
  if (!SHOW_LEGEND) {
    p <- p + theme(legend.position = "none")
  }
  p
}

plot_pca <- function(dt, dtype, stage, out_dir) {
  if (is.null(dt) || !nrow(dt)) return(invisible(NULL))
  stage_lab <- STAGE_LABEL[[stage]] %||% stage
  studies <- unique(na.omit(dt$study))
  for (study_id in studies) {
    dt_study <- dt[study == study_id]
    if (nrow(dt_study) < 3) next
    p <- build_pca_plot(dt_study, dtype, stage, stage_lab, study_id)
    out_file <- file.path(out_dir, sprintf("pca_%s_%s_%s.png", dtype, tolower(stage), study_id))
    ggsave(out_file, p, width = 6.5, height = 5, dpi = 300)
    say("[QC] %s", out_file)
  }
}

plot_pca_pairs <- function(stage_results, dtype, out_dir) {
  stage_names <- intersect(names(STAGE_DIRS), names(stage_results))
  if (length(stage_names) < 2) return(invisible(NULL))
  combined <- rbindlist(lapply(stage_names, function(stage) {
    dt <- stage_results[[stage]]
    if (is.null(dt) || !nrow(dt)) return(NULL)
    copy_dt <- copy(dt)
    stage_name <- stage
    stage_label <- STAGE_LABEL[[stage_name]]
    if (is.null(stage_label) || is.na(stage_label)) stage_label <- stage_name
    copy_dt[, stage := stage_name]
    copy_dt[, stage_label := stage_label]
    copy_dt
  }), fill = TRUE)
  if (!nrow(combined)) return(invisible(NULL))
  stage_label_levels <- STAGE_LABEL[stage_names]
  stage_label_levels <- ifelse(is.na(stage_label_levels), stage_names, stage_label_levels)
  studies <- unique(combined$study)
  for (study_id in studies) {
    dt_study <- combined[study == study_id]
    if (!nrow(dt_study) || uniqueN(dt_study$stage) < 2) next
    dt_plot <- copy(dt_study)
    dt_plot[, batch := as.character(batch)]
    dt_plot[!nzchar(batch) | is.na(batch), batch := "Unknown"]
    dt_plot[, stage_label := factor(stage_label, levels = stage_label_levels)]
    cols <- batch_colors_for(dt_plot$batch)
    if (!length(cols)) cols <- setNames("#1f78b4", unique(dt_plot$batch))
    p <- ggplot(dt_plot, aes(x = PC1, y = PC2, colour = batch)) +
      geom_point(size = 2.2, alpha = 0.9) +
      labs(
        title = sprintf("PCA pair – %s – %s", dtype, study_id),
        x = "PC1",
        y = "PC2",
        colour = "Batch"
      ) +
      facet_wrap(~stage_label, nrow = 1) +
      theme_bw()
    unique_cancers <- unique(na.omit(dt_plot$cancer))
    if (length(unique_cancers) > 1) {
      p <- p + aes(shape = cancer) + labs(shape = "Cancer")
    }
    p <- p + scale_colour_manual(values = cols, drop = FALSE)
    if (!SHOW_LEGEND) {
      p <- p + theme(legend.position = "none")
    }
    out_file <- file.path(out_dir, sprintf("pca_pair_%s_%s.png", dtype, study_id))
    ggsave(out_file, p, width = 11, height = 5, dpi = 300)
    say("[QC] %s", out_file)
  }
}

prepare_cluster_components <- function(dt) {
  dt_study <- unique(dt, by = "sample_id")
  if (!nrow(dt_study) || nrow(dt_study) < 3) return(NULL)
  dt_study[, batch := as.character(batch)]
  dt_study[!nzchar(batch) | is.na(batch), batch := "Unknown"]
  mat <- as.matrix(dt_study[, .(PC1, PC2)])
  rownames(mat) <- dt_study$sample_id
  if (nrow(mat) < 3) return(NULL)
  hc <- hclust(dist(mat))
  batch_levels <- sort(unique(dt_study$batch))
  batch_cols <- if (length(batch_levels)) batch_colors_for(batch_levels) else NULL
  if (!length(batch_cols)) {
    batch_cols <- setNames("#444444", "Unknown")
  }
  label_cex <- min(0.85, 45 / nrow(dt_study))
  label_size <- max(2, label_cex * 4)
  dend_data <- ggdendro::dendro_data(as.dendrogram(hc), type = "rectangle")
  segments_dt <- as.data.table(dend_data$segments)
  labels_dt <- as.data.table(dend_data$labels)
  labels_dt[, batch := dt_study$batch[match(label, dt_study$sample_id)]]
  labels_dt[is.na(batch) | !nzchar(batch), batch := "Unknown"]
  labels_dt[, colour := batch_cols[batch]]
  labels_dt[is.na(colour), colour := "#444444"]
  max_height <- max(segments_dt$y, na.rm = TRUE)
  labels_dt[, y := y - max_height * 0.03]
  list(
    segments = segments_dt,
    labels = labels_dt,
    bars = labels_dt[, .(x, batch)],
    batch_cols = batch_cols,
    label_size = label_size,
    y_max = max_height
  )
}

build_cluster_plot <- function(comp, title_text) {
  colour_values <- comp$batch_cols
  y_lower <- min(comp$labels$y, na.rm = TRUE)
  y_upper <- comp$y_max * 1.02
  span <- y_upper - y_lower
  bar_top <- y_lower - span * 0.03
  bar_bottom <- y_lower - span * 0.08
  labels_y <- bar_bottom - span * 0.03
  label_data <- copy(comp$labels)
  label_data[, y := labels_y]
  bar_data <- copy(comp$bars)
  bar_data[, `:=`(y = bar_bottom, yend = bar_top)]
  y_min <- labels_y - span * 0.05
  p <- ggplot() +
    geom_segment(
      data = comp$segments,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour = "#555555",
      linewidth = 0.35
    ) +
    geom_segment(
      data = bar_data,
      aes(x = x, xend = x, y = y, yend = yend, colour = batch),
      lineend = "butt",
      linewidth = 3.5
    ) +
    geom_text(
      data = label_data,
      aes(x = x, y = y, label = label, colour = batch),
      angle = 90,
      hjust = 1,
      size = comp$label_size
    ) +
    labs(title = title_text, x = "", y = "Height", colour = "Batch") +
    scale_colour_manual(values = colour_values, drop = FALSE) +
    coord_cartesian(ylim = c(y_min, y_upper), expand = FALSE) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  if (!SHOW_LEGEND) {
    p <- p + theme(legend.position = "none")
  }
  p
}

plot_cluster <- function(dt, dtype, stage, out_dir) {
  if (is.null(dt) || !nrow(dt)) return(invisible(NULL))
  stage_lab <- STAGE_LABEL[[stage]] %||% stage
  studies <- unique(na.omit(dt$study))
  for (study_id in studies) {
    dt_study <- unique(dt[study == study_id], by = "sample_id")
    components <- prepare_cluster_components(dt_study)
    if (is.null(components)) next
    out_file <- file.path(out_dir, sprintf("cluster_%s_%s_%s.png", dtype, tolower(stage), study_id))
    plot_obj <- build_cluster_plot(
      components,
      sprintf("Clustering – %s (%s) – %s", dtype, stage_lab, study_id)
    )
    ggsave(out_file, plot_obj, width = 7, height = 5, dpi = 300)
    say("[QC] %s", out_file)
  }
}

plot_cluster_pairs <- function(stage_results, dtype, out_dir) {
  stage_names <- intersect(names(STAGE_DIRS), names(stage_results))
  if (length(stage_names) < 2) return(invisible(NULL))
  studies <- unique(unlist(lapply(stage_results, function(dt) unique(dt$study))))
  studies <- studies[!is.na(studies)]
  if (!length(studies)) return(invisible(NULL))
  for (study_id in studies) {
    components_list <- lapply(stage_names, function(stage) {
      dt <- stage_results[[stage]]
      if (is.null(dt) || !nrow(dt)) return(NULL)
      prepare_cluster_components(dt[study == study_id])
    })
    valid_idx <- which(vapply(components_list, function(x) !is.null(x), logical(1)))
    if (length(valid_idx) < 2) next
    kept_stages <- stage_names[valid_idx]
    kept_components <- components_list[valid_idx]
    out_file <- file.path(out_dir, sprintf("cluster_pair_%s_%s.png", dtype, study_id))
    plot_list <- lapply(seq_along(kept_components), function(i) {
      stage <- kept_stages[i]
      stage_lab <- STAGE_LABEL[[stage]] %||% stage
      build_cluster_plot(
        kept_components[[i]],
        sprintf("%s – %s", stage_lab, study_id)
      )
    })
    grDevices::png(out_file, width = 1400, height = 700, res = 140)
    gridExtra::grid.arrange(grobs = plot_list, nrow = 1)
    grDevices::dev.off()
    say("[QC] %s", out_file)
  }
}

plot_overview <- function(stage_results, dtype, out_dir) {
  required_stages <- c("RAW", "BATCH")
  if (!all(required_stages %in% names(stage_results))) return(invisible(NULL))
  studies <- unique(unlist(lapply(stage_results[required_stages], function(dt) unique(dt$study))))
  studies <- studies[!is.na(studies)]
  if (!length(studies)) return(invisible(NULL))
  for (study_id in studies) {
    skip_study <- FALSE
    pca_grobs <- list()
    cluster_grobs <- list()
    for (stage in required_stages) {
      dt_stage <- stage_results[[stage]]
      if (is.null(dt_stage) || !nrow(dt_stage)) {
        skip_study <- TRUE
        break
      }
      dt_study <- dt_stage[study == study_id]
      if (nrow(dt_study) < 3) {
        skip_study <- TRUE
        break
      }
      stage_lab <- STAGE_LABEL[[stage]] %||% stage
      pca_plot <- build_pca_plot(dt_study, dtype, stage, stage_lab, study_id)
      components <- prepare_cluster_components(dt_study)
      if (is.null(components)) {
        skip_study <- TRUE
        break
      }
      cluster_plot <- build_cluster_plot(
        components,
        sprintf("Clustering – %s (%s) – %s", dtype, stage_lab, study_id)
      )
      pca_grobs[[stage]] <- pca_plot
      cluster_grobs[[stage]] <- cluster_plot
    }
    if (skip_study || length(pca_grobs) < length(required_stages) ||
        length(cluster_grobs) < length(required_stages)) next
    grobs <- list(
      pca_grobs[["RAW"]],
      pca_grobs[["BATCH"]],
      cluster_grobs[["RAW"]],
      cluster_grobs[["BATCH"]]
    )
    layout <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
    out_file <- file.path(out_dir, sprintf("overview_%s_%s.png", dtype, study_id))
    grDevices::png(out_file, width = 1600, height = 1200, res = 140)
    gridExtra::grid.arrange(grobs = grobs, layout_matrix = layout)
    grDevices::dev.off()
    say("[QC] %s", out_file)
  }
}

# ------------------------------------------------------------
# Main loop
# ------------------------------------------------------------

for (dtype in DATA_TYPES) {
  say("=== Data type: %s ===", dtype)
  out_dir <- file.path(OUT_ROOT, dtype)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  stage_results <- list()
  for (stage in names(STAGE_DIRS)) {
    stage_dir <- file.path(STAGE_DIRS[[stage]], dtype)
    if (!dir.exists(stage_dir)) {
      say("  [%s] Missing directory: %s", stage, stage_dir)
      next
    }
    files <- list.files(stage_dir, pattern = "\\.csv$", full.names = TRUE)
    if (!length(files)) {
      say("  [%s] No CSV files under %s", stage, stage_dir)
      next
    }
    files <- files[!grepl("_reference", basename(files), ignore.case = TRUE)]
    if (!length(files)) {
      say("  [%s] No non-reference CSVs under %s", stage, stage_dir)
      next
    }
    say("  [%s] Found %d non-reference file(s) in %s", stage, length(files), stage_dir)
    pca_list <- lapply(files, function(f) {
      mat <- read_matrix(f)
      dataset_name <- tools::file_path_sans_ext(basename(f))
      if (is.null(mat)) {
        say("[SKIP] %s: could not read matrix", dataset_name)
        return(NULL)
      }
      study_id <- study_from_filename(dataset_name)
      res <- run_pca(mat, META, study_id, dataset_name)
      if (is.null(res)) return(NULL)
      res[, dataset := dataset_name]
      res
    })
    valid_entries <- Filter(function(x) !is.null(x) && nrow(x), pca_list)
    if (!length(valid_entries)) {
      say("  [%s] No usable PCA samples.", stage)
      next
    }
    pca_dt <- rbindlist(valid_entries, fill = TRUE)
    plot_pca(pca_dt, dtype, stage, out_dir)
    plot_cluster(pca_dt, dtype, stage, out_dir)
    stage_results[[stage]] <- pca_dt
  }
  plot_pca_pairs(stage_results, dtype, out_dir)
  plot_cluster_pairs(stage_results, dtype, out_dir)
  plot_overview(stage_results, dtype, out_dir)
}

say("Done. Plots written to %s", OUT_ROOT)

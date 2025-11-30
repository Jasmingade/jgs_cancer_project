#!/usr/bin/env Rscript

# Filter proteomics matrices to keep features with > threshold proportion of
# samples having positive (>0) values.
# Usage:
#   Rscript 02_filter_positive_expression.R \
#     <input_root> <output_root> <min_prop> [comma_separated_dtypes]
# Example:
#   Rscript 02_filter_positive_expression.R \
#     02_proteomics/out/preprocessed \
#     02_proteomics/out/preprocessed_filtered \
#     0.2 \
#     gene,iso_log

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 6) {
  stop("Usage: 02_filter_positive_expression.R <input_root> <output_root> <min_prop> [comma_sep_dtypes] [clinical_dir] [target_dataset]")
}

input_root <- args[[1]]
output_root <- args[[2]]
min_prop <- as.numeric(args[[3]])
if (is.na(min_prop) || min_prop < 0 || min_prop > 1) {
  stop("min_prop must be between 0 and 1 (exclusive).")
}

dtype_arg <- if (length(args) >= 4) args[[4]] else "gene,iso_log"
dtypes <- unique(trimws(strsplit(dtype_arg, ",")[[1]]))
dtypes <- dtypes[nzchar(dtypes)]
if (!length(dtypes)) stop("No datatypes specified.")

clinical_dir <- if (length(args) >= 5) args[[5]] else NA_character_
target_dataset <- if (length(args) == 6) args[[6]] else NA_character_

if (!dir.exists(input_root)) stop("Input root not found: ", input_root)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

say <- function(fmt, ...) cat(sprintf(paste0("[filter] ", fmt, "\n"), ...))

summary_rows <- list()

iso_frac_name <- function(base_name) {
  if (grepl("_iso_log(\\.csv)?$", base_name, ignore.case = TRUE)) {
    sub("_iso_log(\\.csv)?$", "_iso_frac.csv", base_name, ignore.case = TRUE)
  } else if (tolower(base_name) == "iso_log.csv") {
    "iso_frac.csv"
  } else {
    sub("\\.csv$", "_iso_frac.csv", base_name, ignore.case = TRUE)
  }
}

dataset_id_from_file <- function(path, dtype) {
  sub(paste0("_", dtype, "\\.csv$"), "", basename(path), ignore.case = TRUE)
}

for (dtype in dtypes) {
  in_dir <- file.path(input_root, dtype)
  out_dir <- file.path(output_root, dtype)

  if (!dir.exists(in_dir)) {
    say("Skipping %s (missing dir %s)", dtype, in_dir)
    next
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  if (!is.na(target_dataset)) {
    files <- Filter(function(f) dataset_id_from_file(f, dtype) == target_dataset, files)
  }
  if (!length(files)) {
    say("No CSVs found for %s in %s", dtype, in_dir)
    next
  }

  norm_in_dir <- normalizePath(in_dir, winslash = "/", mustWork = TRUE)

  for (expr_file in sort(files)) {
    expr_dt <- fread(expr_file)
    if (ncol(expr_dt) < 2) {
      say("Skipping %s (<2 columns)", expr_file)
      next
    }

    expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
    storage.mode(expr_mat) <- "double"
    total_features <- nrow(expr_mat)
    feature_col <- names(expr_dt)[1]

    if (!total_features) {
      say("%s contains 0 features; copying as-is", expr_file)
      rel_path <- sub(paste0("^", norm_in_dir, "/?"), "", normalizePath(expr_file, winslash = "/", mustWork = TRUE))
      out_file <- file.path(out_dir, rel_path)
      dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
      fwrite(expr_dt, out_file)
      next
    }

    prop_pos <- rowMeans(expr_mat > 0, na.rm = TRUE)
    prop_pos[is.na(prop_pos)] <- 0
    keep_idx <- prop_pos > min_prop
    kept <- sum(keep_idx)

    rel_path <- sub(paste0("^", norm_in_dir, "/?"), "", normalizePath(expr_file, winslash = "/", mustWork = TRUE))
    dataset_name <- gsub("\\.csv$", "", rel_path)
    out_file <- file.path(out_dir, rel_path)
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

    filtered_dt <- expr_dt[keep_idx]
    fwrite(filtered_dt, out_file)

    say("Filtered %s %s: kept %d/%d features (threshold %.2f)",
        dtype, rel_path, kept, total_features, min_prop)

    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      dtype = dtype,
      dataset = dataset_name,
      total_features = total_features,
      kept_features = kept,
      min_prop = min_prop
    )

    if (dtype == "iso_log" && kept > 0) {
      eps <- 1e-6
      frac_mat <- as.matrix(filtered_dt[, -1, with = FALSE])
      storage.mode(frac_mat) <- "double"
      frac_mat <- plogis(frac_mat)
      frac_mat <- pmin(pmax(frac_mat, eps), 1 - eps)
      frac_dt <- cbind(
        filtered_dt[, .SD, .SDcols = feature_col],
        as.data.table(frac_mat)
      )
      setnames(frac_dt, c(feature_col, colnames(frac_mat)))
      frac_dirname <- dirname(rel_path)
      if (frac_dirname == ".") frac_dirname <- ""
      frac_rel <- file.path(frac_dirname, iso_frac_name(basename(rel_path)))
      frac_out <- file.path(output_root, "iso_frac", frac_rel)
      dir.create(dirname(frac_out), recursive = TRUE, showWarnings = FALSE)
      fwrite(frac_dt, frac_out)
      say("Derived iso_frac from filtered %s → %s", rel_path, frac_out)
    }
  }
}

if (length(summary_rows)) {
  summary_dt <- rbindlist(summary_rows)
  summary_file <- file.path(output_root, "filter_summary.csv")
  if (file.exists(summary_file)) {
    existing <- fread(summary_file)
    existing <- existing[!(dtype %in% summary_dt$dtype & dataset %in% summary_dt$dataset)]
    summary_dt <- rbindlist(list(existing, summary_dt), use.names = TRUE, fill = TRUE)
  }
  fwrite(summary_dt, summary_file)
  say("Wrote summary to %s", file.path(output_root, "filter_summary.csv"))

  library(ggplot2)
  plot_dir <- file.path(output_root, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  summary_dt[, study := sub("^(PDC[0-9]+).*", "\\1", dataset)]
  agg <- summary_dt[, .(
    before = sum(total_features),
    after = sum(kept_features)
  ), by = .(dtype, study)]

  agg_long <- melt(
    agg,
    id.vars = c("dtype", "study"),
    variable.name = "stage",
    value.name = "features"
  )
  agg_long[, study := factor(study, levels = rev(unique(study)))]
  retain_plot <- ggplot(agg_long, aes(x = study, y = features, fill = stage)) +
    geom_col(position = "dodge") +
    facet_wrap(~dtype, scales = "free_y") +
    coord_flip() +
    scale_fill_manual(values = c(before = "#bdbdbd", after = "#1b9e77"), labels = c("Before filter", "After filter")) +
    labs(
      title = "Features before vs. after positive-fraction filtering (per study)",
      x = "PDC study",
      y = "Feature count",
      fill = ""
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(file.path(plot_dir, "feature_retention_barplot.pdf"), retain_plot, width = 10, height = 6, units = "in")
  ggsave(file.path(plot_dir, "feature_retention_barplot.png"), retain_plot, width = 10, height = 6, units = "in", dpi = 300)
  say("Saved retention plot to %s/feature_retention_barplot.(pdf|png)", plot_dir)
} else {
  say("No datasets were filtered.")
}

if (!is.na(clinical_dir) && dir.exists(clinical_dir)) {
  clinical_files <- list.files(clinical_dir, pattern = "_clinical\\.csv$", full.names = TRUE)
  if (length(clinical_files)) {
    clin_sum <- rbindlist(lapply(clinical_files, function(cf) {
      dt <- fread(cf)
      event_col <- intersect(c("OS_event", "OS.event", "OS"), names(dt))[1]
      if (is.na(event_col)) {
        say("Warning: %s lacks an OS_event column; counting all as missing", basename(cf))
        event_vals <- rep(NA_real_, nrow(dt))
      } else {
        event_vals <- suppressWarnings(as.numeric(dt[[event_col]]))
      }
      data.table(
        study = tools::file_path_sans_ext(basename(cf)),
        dead = sum(event_vals == 1, na.rm = TRUE),
        alive = sum(event_vals == 0, na.rm = TRUE),
        missing = sum(is.na(event_vals)),
        total = nrow(dt)
      )
    }), fill = TRUE)
    fwrite(clin_sum, file.path(output_root, "clinical_event_summary.csv"))
    say("Clinical event summary written to %s", file.path(output_root, "clinical_event_summary.csv"))
  } else {
    say("No clinical CSVs found in %s; skipping event summary", clinical_dir)
  }
} else {
  say("Clinical dir not provided or missing; skipping event summary.")
}

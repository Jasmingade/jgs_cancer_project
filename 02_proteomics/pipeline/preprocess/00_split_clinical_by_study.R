#!/usr/bin/env Rscript

# Split the merged proteomics clinical metadata into per-study CSVs by case_id.
# Usage:
#   Rscript 00_split_clinical_by_study.R \
#     02_proteomics/data/raw/clinical_all_merged.csv \
#     02_proteomics/data/gene \
#     02_proteomics/data/clinical \
#     02_proteomics/config/cancers.yaml

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(yaml))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop("Usage: 00_split_clinical_by_study.R <clinical_all.csv> <expr_dir> <out_dir> [cancer_config.yaml]")
}

clinical_path <- args[[1]]
expr_dir <- args[[2]]
out_dir <- args[[3]]
config_path <- if (length(args) >= 4) args[[4]] else "02_proteomics/config/cancers.yaml"

if (!file.exists(clinical_path)) stop("Clinical CSV not found: ", clinical_path)
if (!dir.exists(expr_dir)) stop("Expression directory not found: ", expr_dir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

say <- function(fmt, ...) cat(sprintf(paste0("[split] ", fmt, "\n"), ...))

load_cancer_map <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    warning("[split] Cancer config not found: ", path)
    return(character())
  }
  cfg <- yaml::read_yaml(path)
  cancers <- cfg$cancers
  if (is.null(cancers)) {
    warning("[split] No 'cancers' entry in config: ", path)
    return(character())
  }
  if (is.list(cancers) && is.null(names(cancers))) {
    cancers <- unlist(cancers, recursive = TRUE, use.names = TRUE)
  } else if (!is.atomic(cancers)) {
    cancers <- unlist(cancers, recursive = TRUE, use.names = TRUE)
  }
  if (is.null(names(cancers))) {
    warning("[split] Cancer config lacks named entries; ignoring.")
    return(character())
  }
  cancers <- as.character(cancers)
  names(cancers) <- toupper(names(cancers))
  cancers
}

clin <- fread(clinical_path)
if (!"case_id" %in% names(clin)) stop("clinical_all.csv must contain a case_id column.")
clin[, case_id := trimws(as.character(case_id))]

expr_files <- list.files(expr_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(expr_files)) stop("Found no .csv files in ", expr_dir)

legacy_subset_files <- list.files(
  out_dir,
  pattern = "_(normal|tumor)_(clinical\\.csv|missing_case_ids\\.txt)$",
  full.names = TRUE
)

studies <- list()
subset_stats <- list()

root_of <- function(x) sub("_(normal|tumor)$", "", x, perl = TRUE)

study_cancer_map <- load_cancer_map(config_path)
if (length(study_cancer_map)) {
  say("Loaded %d cancer mappings from %s", length(study_cancer_map), config_path)
}

for (expr_file in sort(expr_files)) {
  dataset_name <- sub("_gene\\.csv$", "", basename(expr_file))
  if (dataset_name == basename(expr_file)) {
    say("Skipping %s (expected *_gene.csv)", basename(expr_file))
    next
  }

  study <- root_of(dataset_name)
  hdr <- fread(expr_file, nrows = 1)
  case_ids <- names(hdr)
  if (length(case_ids) <= 1) {
    say("No case columns detected in %s", basename(expr_file))
    next
  }

  case_ids <- unique(trimws(case_ids[-1]))
  case_ids <- case_ids[nzchar(case_ids)]
  if (!length(case_ids)) {
    say("No valid case IDs found in %s", basename(expr_file))
    next
  }

  sub_clin <- clin[case_id %in% case_ids]
  missing_cases <- setdiff(case_ids, sub_clin$case_id)
  say("Processed subset %s (study %s): expr=%d matched=%d missing=%d",
      dataset_name, study, length(case_ids), nrow(sub_clin), length(missing_cases))

  subset_stats[[length(subset_stats) + 1L]] <- data.table(
    study = study,
    subset = dataset_name,
    expr_cases = length(case_ids),
    matched_cases = nrow(sub_clin),
    missing_cases = length(missing_cases)
  )

  if (is.null(studies[[study]])) {
    studies[[study]] <- list(
      data = data.table(),
      order = character(),
      missing = character()
    )
  }

  entry <- studies[[study]]
  entry$order <- unique(c(entry$order, case_ids))
  entry$missing <- unique(c(entry$missing, missing_cases))

  if (nrow(sub_clin)) {
    if (nrow(entry$data)) {
      entry$data <- rbindlist(list(entry$data, sub_clin), fill = TRUE)
      entry$data <- entry$data[!duplicated(case_id)]
    } else {
      entry$data <- copy(sub_clin)
    }
  }

  studies[[study]] <- entry
}

subset_stats_dt <- if (length(subset_stats)) rbindlist(subset_stats) else data.table()

summary_rows <- list()

for (study_name in sort(names(studies))) {
  entry <- studies[[study_name]]
  total_expr_cases <- length(entry$order)
  matched_cases <- nrow(entry$data)
  missing_ids <- unique(entry$missing)
  missing_ids <- missing_ids[nzchar(missing_ids)]
  missing_ids <- sort(missing_ids)
  missing_case_count <- length(missing_ids)

  subset_detail <- ""
  if (nrow(subset_stats_dt)) {
    detail_vec <- subset_stats_dt[study == study_name,
      sprintf("%s(expr=%d,matched=%d,missing=%d)", subset, expr_cases, matched_cases, missing_cases)
    ]
    if (length(detail_vec)) subset_detail <- paste(detail_vec, collapse = "; ")
  }

  study_key <- sub("^(PDC[0-9]+).*", "\\1", study_name)
  study_cancer <- unname(study_cancer_map[study_key])
  if (is.na(study_cancer)) study_cancer <- NA_character_

  out_file <- file.path(out_dir, sprintf("%s_clinical.csv", study_name))
  if (matched_cases > 0) {
    order_cases <- entry$order[entry$order %in% entry$data$case_id]
    entry$data[, case_id := factor(case_id, levels = order_cases)]
    setorder(entry$data, case_id)
    entry$data[, case_id := as.character(case_id)]
    fwrite(entry$data, out_file)
    say("Wrote %s (study %s; matched %d of %d; missing %d)",
        out_file, study_name, matched_cases, total_expr_cases, missing_case_count)
  } else {
    say("No matched clinical rows for study %s (expr cases: %d) → not writing file",
        study_name, total_expr_cases)
    out_file <- NA_character_
  }

  missing_file <- file.path(out_dir, sprintf("%s_missing_case_ids.txt", study_name))
  if (missing_case_count > 0) {
    writeLines(missing_ids, missing_file)
    say("Logged %d missing case IDs for %s to %s", missing_case_count, study_name, missing_file)
  } else if (file.exists(missing_file)) {
    file.remove(missing_file)
  }

  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    study = study_name,
    study_cancer_type = study_cancer,
    subsets = subset_detail,
    total_expr_cases = total_expr_cases,
    matched_case_count = matched_cases,
    missing_case_count = missing_case_count,
    missing_case_ids = paste(missing_ids, collapse = ";"),
    output_file = out_file
  )
}

if (length(summary_rows)) {
  summary_dt <- rbindlist(summary_rows, fill = TRUE)
  summary_file <- file.path(out_dir, "clinical_split_summary.csv")
  fwrite(summary_dt, summary_file)
  say("Summary written to %s", summary_file)
}

if (length(legacy_subset_files)) {
  removed <- legacy_subset_files[file.exists(legacy_subset_files)]
  if (length(removed)) {
    ok <- file.remove(removed)
    if (any(ok)) {
      say("Removed %d legacy subset files (normal/tumor-specific)", sum(ok))
    }
  }
}

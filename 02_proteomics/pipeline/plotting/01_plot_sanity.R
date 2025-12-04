#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(yaml)
})

say <- function(...) message(sprintf(...))

# -------------------------------------------------------------------
# Shared theme & helpers (reuse transcriptomics plotting infra)
# -------------------------------------------------------------------
source("01_transcriptomics/pipeline/plotting/plot_theme.R")
source("01_transcriptomics/pipeline/plotting/helpers_plotting.R")
source("02_proteomics/pipeline/plotting/helpers_proteomics.R", local = TRUE, chdir = TRUE)

# -------------------------------------------------------------------
# Toggles (env vars)
# -------------------------------------------------------------------
to_bool <- function(x) tolower(x) %in% c("1","true","yes")

RUN_PLOT_M0_AGE_FOREST   <- to_bool(Sys.getenv("RUN_PLOT_M0_AGE_FOREST", "true"))
RUN_PLOT_M0_KM_CURVES    <- to_bool(Sys.getenv("RUN_PLOT_M0_KM_CURVES",  "false"))

say("[PROT M0] sanity plotting: age_forest=%s, KM=%s",
    RUN_PLOT_M0_AGE_FOREST, RUN_PLOT_M0_KM_CURVES)

# -------------------------------------------------------------------
# Config & paths
# -------------------------------------------------------------------
# Proteomics cancer mapping
prot_cancer_config <- Sys.getenv("PROT_CANCER_CONFIG", "02_proteomics/config/cancers.yaml")
prot_cfg <- yaml::read_yaml(prot_cancer_config)
cancer_map <- prot_cfg$cancers
cancers  <- if (!is.null(names(cancer_map))) names(cancer_map) else as.character(cancer_map)
cancers  <- as.character(cancers)

# ---- covariates config ----
cov_cfg  <- yaml::read_yaml("02_proteomics/config/covariates.yaml")

out_dir <- "02_proteomics/out/plots/sanity/model0_sanity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


canon_sample_type <- function(x) {
  x <- tolower(x)
  x[is.na(x)] <- ""
  out <- ifelse(grepl("normal", x), "normal", "tumor")
  out
}

# -------------------------------------------------------------------
# 1) Fit Cox(age) for tumor/normal per proteomics cancer/study
# -------------------------------------------------------------------
age_summaries   <- list()
clin_all        <- list()
types           <- c("tumor", "normal")

for (canc in cancers) {

  # -----------------------------------------------------------------
  # Clinical / manifest file for PROTEOMICS
  # -----------------------------------------------------------------
  # Adjust ONLY this line to your actual manifest location:
  clin_glob <- sprintf("02_proteomics/data/clinical/%s_*_clinical.csv", canc)
  clin_matches <- Sys.glob(clin_glob)
  if (!length(clin_matches)) {
    say("[PROT M0] Clinical/manifest file missing for %s → %s", canc, clin_glob)
    next
  }
  clin_file <- clin_matches[[1]]

  if (!file.exists(clin_file)) {
    say("[PROT M0] Clinical/manifest file missing for %s → %s", canc, clin_file)
    next
  }

  say("[PROT M0] Reading proteomics clinical manifest for %s", canc)
  clin <- fread(clin_file)

  # Coerce key columns to numeric and report NAs introduced
  coerce_numeric <- function(dt, col) {
    if (col %in% names(dt)) {
      n_before <- sum(is.na(dt[[col]]))
      suppressWarnings(dt[, (col) := as.numeric(get(col))])
      n_after  <- sum(is.na(dt[[col]]))
      if (n_after > n_before) {
        say("[PROT M0] %s: %d values became NA after numeric coercion", col, n_after - n_before)
      }
    }
  }
  coerce_numeric(clin, "age")
  coerce_numeric(clin, "OS_time")
  coerce_numeric(clin, "OS_event")

  # sample_type handling (prefer batch_annotation; otherwise default tumor)
  annot_lookup <- list()
  annot_matches <- Sys.glob(file.path("02_proteomics/data/batch_annotation",
                                      sprintf("%s_*.csv", canc)))
  if (length(annot_matches) > 0) {
    annot_file <- annot_matches[[1]]
    annot_dt <- try(fread(annot_file), silent = TRUE)
    if (!inherits(annot_dt, "try-error") && "sample_type" %in% names(annot_dt) &&
        "case_id" %in% names(annot_dt)) {
      annot_dt[, sample_type := tolower(sample_type)]
      annot_lookup <- split(annot_dt$sample_type, annot_dt$case_id)
    }
  }
  if (!"sample_type" %in% names(clin)) {
    clin[, sample_type := NA_character_]
  } else {
    clin[, sample_type := tolower(sample_type)]
  }

  needed_cols <- c("case_id", "OS_time", "OS_event", "age")
  if (!all(needed_cols %in% names(clin))) {
    say("[PROT M0] Skipping %s – missing required cols (%s)", canc,
        paste(setdiff(needed_cols, names(clin)), collapse = ", "))
    next
  }

  # Expand rows by available sample types (clinical/sample_type or batch annotation)
  expanded <- list()
  for (i in seq_len(nrow(clin))) {
    row <- clin[i]

    types_case <- unique(na.omit(c(row$sample_type, annot_lookup[[row$case_id]])))
    if (length(types_case) == 0) {
      types_case <- "tumor"
    }

    # collapse everything to "tumor" / "normal"
    types_case <- canon_sample_type(types_case)

    for (tp in types_case) {
      tmp <- copy(row)
      tmp[, sample_type := tp]
      expanded[[length(expanded) + 1]] <- tmp
    }
  }
  clin_expanded <- rbindlist(expanded, use.names = TRUE, fill = TRUE)

  for (s_type in types) {
    clin_use <- clin_expanded[
      !is.na(OS_time) & !is.na(OS_event) & !is.na(age) & sample_type == s_type
    ]

    if (nrow(clin_use) == 0) {
      next
    }

    n_total  <- nrow(clin_use)
    n_event  <- sum(clin_use$OS_event == 1, na.rm = TRUE)

    if (n_total < 20 || n_event < 5) {
      say("[PROT M0] Skipping %s (%s) – too few samples/events (n=%d, events=%d)",
          canc, s_type, n_total, n_event)
      next
    }

    # ---------------- Cox(age) per cancer & type ----------------
    fit_age <- try(coxph(Surv(OS_time, OS_event) ~ age, data = clin_use), TRUE)

    if (inherits(fit_age, "try-error")) {
      say("[PROT M0] Cox(age) failed for %s (%s)", canc, s_type)
    } else {
      s_age <- summary(fit_age)

      # Make sure the age row actually exists
      if (!"age" %in% rownames(s_age$coef)) {
        say("[PROT M0] Age term missing or non-finite for %s (%s) – skipping age summary", canc, s_type)
      } else {
        coef_age <- s_age$coef["age", "coef"]
        ci_age   <- try(suppressWarnings(confint(fit_age)["age", ]), silent = TRUE)

        if (inherits(ci_age, "try-error") || any(!is.finite(ci_age))) {
          say("[PROT M0] confint(age) failed for %s (%s) – skipping age summary", canc, s_type)
        } else {
          age_summaries[[length(age_summaries) + 1]] <- data.table(
            cancer       = canc,
            cancer_label = pretty_cancer(canc),
            sample_type  = s_type,
            n_samples    = n_total,
            n_events     = n_event,
            HR_age       = exp(coef_age),
            HR_age_lo    = exp(ci_age[1]),
            HR_age_hi    = exp(ci_age[2]),
            p_age        = s_age$coef["age", "Pr(>|z|)"]
          )
        }
      }
    }

    clin_use[, cancer := canc]
    clin_use[, sample_type := s_type]
    clin_all[[length(clin_all) + 1]] <- clin_use
  }
}
if (length(age_summaries) == 0) stop("[PROT M0] No valid cancers for age plots.")


age_dt <- rbindlist(age_summaries)
fwrite(age_dt, file.path(out_dir, "prot_m0_age_cox_summary.csv"))
say("[PROT M0] Wrote age-only Cox summary → %s",
    file.path(out_dir, "prot_m0_age_cox_summary.csv"))

# -------------------------------------------------------------------
# 2) FOREST PLOT – age sanity per proteomics cancer
# -------------------------------------------------------------------
if (RUN_PLOT_M0_AGE_FOREST) {

  make_age_plot <- function(dt, suffix) {
    if (nrow(dt) == 0) return(NULL)

    # order by HR (ensure unique levels)
    dt[, cancer_label := factor(
      cancer_label,
      levels = unique(dt[order(HR_age), cancer_label])
    )]

    # significance flag
    dt[, signif := p_age < 0.05]
    dt[, sig_group := factor(
      ifelse(signif, "Significant", "Not significant"),
      levels = c("Not significant", "Significant")
    )]

    # add sample count into the axis label
    dt[, cancer_label_n := sprintf("%s (n=%d)", cancer_label, n_samples)]
    dt[, cancer_label_n := factor(
      cancer_label_n,
      levels = unique(dt[order(HR_age), cancer_label_n])
    )]

    mean_hr_age <- exp(mean(log(dt$HR_age), na.rm = TRUE))

    p_age <- ggplot(dt, aes(x = HR_age, y = cancer_label_n)) +
      geom_point(
        data = dt[sig_group == "Not significant"],
        aes(x = HR_age, y = cancer_label_n, color = sig_group),
        shape = 1, size = 3
      ) +
      geom_point(
        data = dt[sig_group == "Significant"],
        aes(x = HR_age, y = cancer_label_n, color = sig_group),
        shape = 16, size = 3
      ) +
      geom_errorbarh(aes(xmin = HR_age_lo, xmax = HR_age_hi),
                     height = 0.25, alpha = 0.8) +
      geom_vline(xintercept = 1, linetype = "dashed",
                 color = "black", alpha = 0.7) +
      geom_vline(xintercept = mean_hr_age, color = "#c90028", linewidth = 0.7) +
      annotate("label",
               x = mean_hr_age, y = Inf,
               label = sprintf("Mean HR = %.2f", mean_hr_age),
               hjust = -0.1, vjust = 1.3,
               size = 3, fill = "white", color = "#c90028") +
      scale_x_log10() +
      scale_color_manual(
        values = c(
          "Significant"     = "#009E73",
          "Not significant" = "#000000"
        ),
        name = "Significance"
      ) +
      labs(
        title    = sprintf("Proteomics – Age effect (%s)", suffix),
        subtitle = "Per cancer; black dashed = HR=1; red = mean HR",
        x        = "Hazard Ratio (log10 scale)",
        y        = "Cancer Type (with sample size)"
      ) +
      theme_big(13) +
      theme(legend.position = "right", axis.text.y = element_text(size = 10))

    out_age_png <- file.path(out_dir, sprintf("prot_m0_age_forest_%s.png", suffix))
    save_plot(p_age, out_age_png, width = 12, height = 10)
    say("[PROT M0] Saved age sanity forest plot (%s) → %s", suffix, out_age_png)
  }

  for (s_type in types) {
    dt_type <- age_dt[sample_type == s_type]
    if (nrow(dt_type) > 0) {
      make_age_plot(dt_type, s_type)
    } else {
      say("[PROT M0] No age data for %s", s_type)
    }
  }
}
say("[PROT M0] sanity plotting done.")

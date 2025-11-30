#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(yaml)
})

# -------------------------------------------------------------------
# Shared theme & helpers
# -------------------------------------------------------------------
source("01_transcriptomics/pipeline/plotting/plot_theme.R")
source("01_transcriptomics/pipeline/plotting/helpers_plotting.R")

# -------------------------------------------------------------------
# Toggles (env vars)
# -------------------------------------------------------------------
to_bool <- function(x) tolower(x) %in% c("1","true","yes")

RUN_PLOT_M0_AGE_FOREST   <- to_bool(Sys.getenv("RUN_PLOT_M0_AGE_FOREST", "true"))
RUN_PLOT_M0_KM_CURVES    <- to_bool(Sys.getenv("RUN_PLOT_M0_KM_CURVES",  "false"))

say("M0 sanity plotting: age_forest=%s, KM=%s",
    RUN_PLOT_M0_AGE_FOREST, RUN_PLOT_M0_KM_CURVES)

# -------------------------------------------------------------------
# Config & paths
# -------------------------------------------------------------------
cfg      <- yaml::read_yaml("01_transcriptomics/config/cancers.yaml")
cancers  <- paste0("TCGA_", cfg$cancers)

cov_cfg  <- yaml::read_yaml("01_transcriptomics/config/covariates.yaml")

# --- get stage coverage threshold from covariates.yaml ---
stage_cov_threshold <- 0.8   # fallback
if (!is.null(cov_cfg$conditional_covariates)) {
  idx <- vapply(
    cov_cfg$conditional_covariates,
    function(cc) identical(cc$name, "stage"),
    logical(1)
  )
  if (any(idx)) {
    stage_block <- cov_cfg$conditional_covariates[[which(idx)[1]]]
    if (!is.null(stage_block$coverage_threshold)) {
      stage_cov_threshold <- as.numeric(stage_block$coverage_threshold)
    }
  }
}

# --- stage levels (order) from yaml or default ---
if (!is.null(cov_cfg$stage_levels)) {
  stage_levels <- as.character(cov_cfg$stage_levels)
} else {
  stage_levels <- c("I", "II", "III", "IV")
}

say("Stage coverage threshold for inclusion: %.2f", stage_cov_threshold)
say("Stage levels (ordered): %s", paste(stage_levels, collapse = ", "))

out_dir <- "01_transcriptomics/out/05_plots/model0_sanity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# 1) Fit Cox(age) and Cox(stage) per cancer
# -------------------------------------------------------------------
age_summaries   <- list()
stage_summaries <- list()
clin_all        <- list()

for (canc in cancers) {

  # Using the sample_manifest as clinical source
  clin_file <- sprintf("01_transcriptomics/out/02_norm/%s_gene.sample_manifest.csv", canc)

  if (!file.exists(clin_file)) {
    say("Clinical/manifest file missing for %s → %s", canc, clin_file)
    next
  }

  say("Reading clinical manifest for %s", canc)
  clin <- fread(clin_file)

  needed_cols <- c("case_id", "OS_time", "OS_event", "age")
  if (!all(needed_cols %in% names(clin))) {
    say("Skipping %s – missing required cols (%s)", canc,
        paste(setdiff(needed_cols, names(clin)), collapse = ", "))
    next
  }

  # Drop NAs in survival + age
  clin_use <- clin[!is.na(OS_time) & !is.na(OS_event) & !is.na(age)]
  n_total  <- nrow(clin_use)
  n_event  <- sum(clin_use$OS_event == 1, na.rm = TRUE)

  if (n_total < 20 || n_event < 5) {
    say("Skipping %s – too few samples/events (n=%d, events=%d)", canc, n_total, n_event)
    next
  }

  # ---------------- Cox(age) per cancer ----------------
  fit_age <- try(coxph(Surv(OS_time, OS_event) ~ age, data = clin_use), TRUE)
  if (inherits(fit_age, "try-error")) {
    say("Cox(age) failed for %s", canc)
  } else {
    coef_age <- coef(fit_age)["age"]
    ci_age   <- suppressWarnings(confint(fit_age)["age", ])
    s_age    <- summary(fit_age)

    age_summaries[[length(age_summaries) + 1]] <- data.table(
      cancer       = canc,
      cancer_label = pretty_cancer(canc),
      n_samples    = n_total,
      n_events     = n_event,
      HR_age       = exp(coef_age),
      HR_age_lo    = exp(ci_age[1]),
      HR_age_hi    = exp(ci_age[2]),
      p_age        = s_age$coef["age", "Pr(>|z|)"]
    )
  }

  clin_use[, cancer := canc]
  clin_all[[length(clin_all) + 1]] <- clin_use

  # ---------------- Cox(stage) per cancer ----------------
  if (!"stage" %in% names(clin_use)) {
    say("  No 'stage' column in manifest for %s → skipping stage sanity", canc)
    next
  }

  stage_cov <- mean(!is.na(clin_use$stage))
  if (is.na(stage_cov) || stage_cov < stage_cov_threshold) {
    say("  Stage coverage %.2f < %.2f for %s → skipping stage sanity",
        stage_cov, stage_cov_threshold, canc)
    next
  }

  say("  Stage coverage for %s: %.2f ≥ %.2f → fitting Cox(stage)",
      canc, stage_cov, stage_cov_threshold)

  # keep only samples with stage
  clin_stage <- clin_use[!is.na(stage)]

  # map to ordered levels I–IV (anything not in stage_levels becomes NA and dropped)
  clin_stage[, stage_clean := factor(stage, levels = stage_levels, ordered = TRUE)]
  clin_stage <- clin_stage[!is.na(stage_clean)]

  if (nrow(clin_stage) < 20) {
    say("  After stage cleaning: too few samples for %s (n=%d)", canc, nrow(clin_stage))
    next
  }

  # require at least 2 distinct stage levels represented
  if (length(unique(clin_stage$stage_clean)) < 2) {
    say("  Only one stage level present for %s → skipping stage sanity", canc)
    next
  }

  n_stage_event <- sum(clin_stage$OS_event == 1, na.rm = TRUE)
  if (n_stage_event < 5) {
    say("  Too few events for stage model in %s (events=%d)", canc, n_stage_event)
    next
  }

  # numeric ordinal coding (I=1, II=2, ...)
  clin_stage[, stage_ord := as.numeric(stage_clean)]

  fit_stage <- try(coxph(Surv(OS_time, OS_event) ~ stage_ord, data = clin_stage), TRUE)
  if (inherits(fit_stage, "try-error")) {
    say("  Cox(stage_ord) failed for %s", canc)
  } else {
    coef_stage <- coef(fit_stage)["stage_ord"]
    ci_stage   <- suppressWarnings(confint(fit_stage)["stage_ord", ])
    s_stage    <- summary(fit_stage)

    stage_summaries[[length(stage_summaries) + 1]] <- data.table(
      cancer        = canc,
      cancer_label  = pretty_cancer(canc),
      n_samples     = nrow(clin_stage),
      n_events      = n_stage_event,
      HR_stage      = exp(coef_stage),
      HR_stage_lo   = exp(ci_stage[1]),
      HR_stage_hi   = exp(ci_stage[2]),
      p_stage       = s_stage$coef["stage_ord", "Pr(>|z|)"],
      stage_cov     = stage_cov
    )
  }
}

if (length(age_summaries) == 0) stop("No valid cancers for M0 age plots.")

age_dt <- rbindlist(age_summaries)
fwrite(age_dt, file.path(out_dir, "m0_age_only_cox_summary.csv"))
say("Wrote age-only Cox summary → %s",
    file.path(out_dir, "m0_age_only_cox_summary.csv"))

if (length(stage_summaries) > 0) {
  stage_dt <- rbindlist(stage_summaries)
  fwrite(stage_dt, file.path(out_dir, "m0_stage_only_cox_summary.csv"))
  say("Wrote stage-only Cox summary → %s",
      file.path(out_dir, "m0_stage_only_cox_summary.csv"))
} else {
  stage_dt <- NULL
  say("No valid cancers for stage sanity plot (stage_dt is NULL).")
}

# -------------------------------------------------------------------
# 2) FOREST PLOT – age sanity per cancer
# -------------------------------------------------------------------
if (RUN_PLOT_M0_AGE_FOREST) {

  # ---------------- AGE FOREST ----------------
  dt <- copy(age_dt)

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

  # mean HR (overall)
  mean_hr_age <- exp(mean(log(dt$HR_age), na.rm = TRUE))

  p_age <- ggplot(dt, aes(x = HR_age, y = cancer_label_n)) +

    # NON-SIGNIFICANT: hollow black circles
    geom_point(
      data = dt[sig_group == "Not significant"],
      aes(x = HR_age, y = cancer_label_n, color = sig_group),
      shape = 1,
      size  = 3
    ) +
    # SIGNIFICANT: filled green circles
    geom_point(
      data = dt[sig_group == "Significant"],
      aes(x = HR_age, y = cancer_label_n, color = sig_group),
      shape = 16,
      size  = 3
    ) +
    geom_errorbarh(aes(xmin = HR_age_lo, xmax = HR_age_hi),
                   height = 0.25, alpha = 0.8) +

    # HR=1 reference (black)
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "black", alpha = 0.7) +

    # Mean HR annotation (red)
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
      title    = "Clinical Covariate: Age Effect on Survival",
      subtitle = "Per cancer; black dashed = HR=1; red = mean HR",
      x        = "Hazard Ratio (log10 scale)",
      y        = "Cancer Type (with sample size)"
    ) +

    theme_big(13) +
    theme(
      legend.position = "right",
      axis.text.y     = element_text(size = 10)
    )

  out_age_png <- file.path(out_dir, "m0_age_forest.png")
  save_plot(p_age, out_age_png, width = 12, height = 10)
  say("Saved age sanity forest plot → %s", out_age_png)

  # ---------------- STAGE FOREST ----------------
  if (!is.null(stage_dt) && nrow(stage_dt) > 0) {

    ds <- copy(stage_dt)

    # order by HR_stage (ensure unique levels)
    ds[, cancer_label := factor(
      cancer_label,
      levels = unique(ds[order(HR_stage), cancer_label])
    )]

    # significance flag
    ds[, signif := p_stage < 0.05]
    ds[, sig_group := factor(
      ifelse(signif, "Significant", "Not significant"),
      levels = c("Not significant", "Significant")
    )]

    # add sample count into the axis label
    ds[, cancer_label_n := sprintf("%s (n=%d)", cancer_label, n_samples)]
    ds[, cancer_label_n := factor(
      cancer_label_n,
      levels = unique(ds[order(HR_stage), cancer_label_n])
    )]

    # mean HR across cancers (stage effect)
    mean_hr_stage <- exp(mean(log(ds$HR_stage), na.rm = TRUE))

    p_stage <- ggplot(ds, aes(x = HR_stage, y = cancer_label_n)) +

      geom_point(
        data = ds[sig_group == "Not significant"],
        aes(x = HR_stage, y = cancer_label_n, color = sig_group),
        shape = 1,
        size  = 3
      ) +
      geom_point(
        data = ds[sig_group == "Significant"],
        aes(x = HR_stage, y = cancer_label_n, color = sig_group),
        shape = 16,
        size  = 3
      ) +
      geom_errorbarh(aes(xmin = HR_stage_lo, xmax = HR_stage_hi),
                     height = 0.25, alpha = 0.8) +

      geom_vline(xintercept = 1, linetype = "dashed",
                 color = "black", alpha = 0.7) +
      geom_vline(xintercept = mean_hr_stage, color = "#c90028", linewidth = 0.7) +
      annotate("label",
               x = mean_hr_stage, y = Inf,
               label = sprintf("Mean HR = %.2f", mean_hr_stage),
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
        title    = "Clinical Covariate: Stage Effect on Survival",
        subtitle = sprintf(
          "Per cancer; only cancers with ≥ %.0f%% valid stage (I–IV)",
          stage_cov_threshold * 100
        ),
        x        = "Hazard Ratio (per one-stage increase, log10 scale)",
        y        = "Cancer Type (with sample size)"
      ) +

      theme_big(13) +
      theme(
        legend.position = "right",
        axis.text.y     = element_text(size = 10)
      )

    out_stage_png <- file.path(out_dir, "m0_stage_forest.png")
    save_plot(p_stage, out_stage_png, width = 12, height = 10)
    say("Saved stage sanity forest plot → %s", out_stage_png)

  } else {
    say("Stage sanity forest plot skipped: no valid cancers with ≥ %.2f coverage.",
        stage_cov_threshold)
  }
}

say("M0 sanity plotting done.")

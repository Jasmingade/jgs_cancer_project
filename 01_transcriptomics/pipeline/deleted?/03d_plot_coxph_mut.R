#!/usr/bin/env Rscript
# ============================================================
# 03_plot_coxph.R
# ------------------------------------------------------------
# Generates grouped HR boxplots for baseline, combined, and interaction
# and a comparison across analyses + summary statistics.
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

say <- function(...) message(sprintf(...))
winsorize <- function(x, p=0.01) {
  if (is.null(p) || p <= 0) return(x)
  lo <- quantile(x, p, na.rm=TRUE); hi <- quantile(x, 1-p, na.rm=TRUE)
  x[x < lo] <- lo; x[x > hi] <- hi; x
}

# ------------------ Parse Args ------------------
args <- commandArgs(trailingOnly = TRUE)
kv <- list()
for (i in seq(1, length(args), 2)) kv[[sub("^--", "", args[i])]] <- args[i+1]
mode <- kv$mode
if (is.null(mode)) stop("Missing --mode argument")

# ------------------ Output Directories ------------------
out_dir <- kv$out_dir %||% kv$in_dir %||% "01_transcriptomics/out/03_plots_mut"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------ Color Palette ------------------
pal <- c(Expression = "#0072B2", Mutation = "#D55E00")

# ============================================================
# MODE-SPECIFIC PROCESSING
# ============================================================
if (mode %in% c("baseline", "combined", "interaction")) {

  say("[MODE] %s", toupper(mode))
  index_file <- kv$index
  stopifnot(file.exists(index_file))

  idx <- fread(index_file)
  stopifnot(all(c("path","cancer","data_type") %in% names(idx)))

  all_dt <- rbindlist(lapply(seq_len(nrow(idx)), function(i) {
    fp <- idx$path[i]
    if (!file.exists(fp)) return(NULL)
    dt <- try(fread(fp, select=c("HR","p","FDR","beta")), silent=TRUE)
    if (inherits(dt, "try-error")) return(NULL)
    if (!"HR" %in% names(dt) && "beta" %in% names(dt)) dt[, HR := exp(beta)]
    dt[, cancer := idx$cancer[i]]
    dt[, data_type := idx$data_type[i]]
    dt
  }), fill = TRUE)

  if (nrow(all_dt) == 0) stop("No valid results for mode ", mode)

  # Clean & annotate
  all_dt <- all_dt[is.finite(HR) & HR > 0]
  all_dt[, HR_w := winsorize(HR, 0.01)]
  all_dt[, group := ifelse(grepl("^mutation", data_type, ignore.case=TRUE), "Mutation", "Expression")]

  say("[INFO] Loaded %d rows", nrow(all_dt))

  # Plot grouped boxplot
  p <- ggplot(all_dt, aes(x=cancer, y=log2(HR_w), fill=group)) +
    geom_hline(yintercept=0, linetype="dashed", color="gray40") +
    geom_boxplot(outlier.shape=NA, alpha=0.85) +
    scale_fill_manual(values=pal) +
    labs(title=sprintf("Hazard Ratios by Cancer — %s", toupper(mode)),
         subtitle="log2(HR), winsorized at 1%",
         y="log2(HR)", x="Cancer type", fill="Data Group") +
    theme_bw(base_size=12) +
    theme(axis.text.x=element_text(angle=60, hjust=1))

  out_png <- file.path(out_dir, sprintf("cox_boxplot_%s.png", mode))
  ggsave(out_png, p, width=14, height=8, dpi=300)
  say("[OK] Saved: %s", out_png)

  fwrite(all_dt, file.path(out_dir, sprintf("cox_data_%s.csv", mode)))

} else if (mode == "compare") {

  say("[MODE] COMPARISON ACROSS MODES")
  modes <- c("baseline","combined","interaction")
  all_res <- list()

  for (m in modes) {
    fp <- file.path(out_dir, sprintf("cox_data_%s.csv", m))
    if (file.exists(fp)) {
      dt <- fread(fp)
      dt[, mode := m]
      all_res[[m]] <- dt
    }
  }

  dt_all <- rbindlist(all_res, fill=TRUE)
  if (nrow(dt_all) == 0) stop("No mode data found for comparison")

  # Summary stats
  summary_dt <- dt_all[, .(
    median_HR = median(HR, na.rm=TRUE),
    mean_HR   = mean(HR, na.rm=TRUE),
    n = .N
  ), by=.(mode, data_type)]

  fwrite(summary_dt, file.path(out_dir, "cox_summary_statistics.csv"))
  say("[OK] Summary stats written")

  # Comparison plot
  comp_plot <- ggplot(summary_dt, aes(x=data_type, y=log2(median_HR), fill=mode)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot(position=position_dodge(width=0.8), alpha=0.8) +
    labs(title="Comparison of Median HRs Across Analyses",
         subtitle="log2(Median HR)",
         y="log2(Median HR)", x="Data type", fill="Analysis Mode") +
    theme_bw(base_size=13) +
    theme(axis.text.x=element_text(angle=45, hjust=1))
  
  ggsave(file.path(out_dir, "cox_comparison_across_modes.png"), comp_plot, width=12, height=6, dpi=300)
  say("[OK] Comparison plot saved")
}

say("[DONE] Completed mode: %s", mode)

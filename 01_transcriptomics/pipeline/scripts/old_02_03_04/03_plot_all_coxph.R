#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# ----------------------
# Args
# ----------------------
args <- commandArgs(trailingOnly = TRUE)
kv <- list(); i <- 1
while (i <= length(args)) {
  a <- args[[i]]
  if (grepl("^--", a)) {
    key <- sub("^--", "", a); val <- TRUE
    if (i + 1 <= length(args) && !grepl("^--", args[[i+1]])) { val <- args[[i+1]]; i <- i + 1 }
    kv[[key]] <- val
  } else stop(sprintf("Unexpected arg: %s", a))
  i <- i + 1
}
get_opt <- function(k, def=NULL) if (!is.null(kv[[k]])) kv[[k]] else def

index_path <- get_opt("index")
out_path   <- get_opt("out",   "cox_hr_boxplots.png")
title_txt  <- get_opt("title", "Hazard Ratio distributions by cancer")
flt_expr   <- get_opt("filter", "")
winsor_p   <- as.numeric(get_opt("winsor", 0))
clip_str   <- get_opt("clip", "")
use_logy   <- isTRUE(as.logical(get_opt("logy", TRUE)))
order_key  <- get_opt("order", "medianHR")     # "medianHR" | "alphabet"
width_in   <- as.numeric(get_opt("width", 16))
height_in  <- as.numeric(get_opt("height", 9))
dpi_out    <- as.integer(get_opt("dpi", 300))
jitter_n   <- as.integer(get_opt("jitter_n", 0))   # 0 = no dots (clean)
jitter_w   <- as.numeric(get_opt("jitter_width", 0.25))

stopifnot(!is.null(index_path), file.exists(index_path))

# ----------------------
# Helpers
# ----------------------
winsorize <- function(x, p=0) {
  if (is.null(p) || is.na(p) || p <= 0) return(x)
  lo <- quantile(x, p, na.rm=TRUE); hi <- quantile(x, 1-p, na.rm=TRUE)
  x[x < lo] <- lo; x[x > hi] <- hi; x
}
parse_clip <- function(s) {
  if (!nzchar(s)) return(NULL)
  pr <- strsplit(s, ",")[[1]]
  if (length(pr)!=2) stop("--clip needs 'lo,hi' e.g. 0.01,0.99")
  as.numeric(pr)
}

# ----------------------
# Read index + Cox files
# ----------------------
idx <- fread(index_path)
stopifnot(all(c("path","cancer","data_type") %in% names(idx)))

read_one <- function(fp) {
  dt <- try(fread(fp), silent = TRUE); if (inherits(dt, "try-error")) return(NULL)

  # keep only expected cols if present
  keep <- intersect(c("feature","HR","p","FDR","beta","HR_lo","HR_hi","logHR"), names(dt))
  dt <- dt[, ..keep]

  # coerce HR and beta to numeric robustly (strip commas etc.)
  if ("HR" %in% names(dt)) {
    dt[, HR := suppressWarnings(as.numeric(gsub(",", "", as.character(HR))))]
  }
  if ("beta" %in% names(dt)) {
    dt[, beta := suppressWarnings(as.numeric(gsub(",", "", as.character(beta))))]
  }

  # If HR missing but beta present -> derive HR
  if (!"HR" %in% names(dt) && "beta" %in% names(dt)) {
    dt[, HR := exp(beta)]
  }

  # Build/repair logHR: prefer beta (already log(HR)); else log(HR) if HR>0
  if (!"logHR" %in% names(dt)) {
    if ("beta" %in% names(dt) && any(is.finite(dt$beta))) {
      dt[, logHR := beta]
    } else if ("HR" %in% names(dt)) {
      # guard: only log positive finite HR
      dt[!(is.finite(HR) & HR > 0), HR := NA_real_]
      dt[, logHR := suppressWarnings(log(HR))]
    } else {
      return(NULL)
    }
  }

  # Quick sanity check: HR vs exp(beta), only when both are usable
  if (all(c("HR","beta") %in% names(dt))) {
    ok <- is.finite(dt$HR) & dt$HR > 0 & is.finite(dt$beta)
    if (any(ok)) {
      cval <- suppressWarnings(cor(log(dt$HR[ok]), dt$beta[ok], use = "pairwise.complete.obs"))
      if (!is.na(cval) && cval < 0.999) {
        warning(sprintf("HR and beta mismatch in %s (corr=%.3f)", basename(fp), cval))
      }
    }
  }

  dt
}


all_dt <- rbindlist(lapply(seq_len(nrow(idx)), function(i) {
  fp <- idx$path[i]
  dt <- read_one(fp); if (is.null(dt)) return(NULL)
  if (nzchar(flt_expr)) {
    dt <- try(dt[eval(parse(text=flt_expr))], silent=TRUE)
    if (inherits(dt, "try-error")) stop(sprintf("Invalid --filter: %s", flt_expr))
  }
  if (!is.numeric(dt$HR)) dt[, HR := as.numeric(HR)]
  dt[, `:=`(cancer = idx$cancer[i], data_type = idx$data_type[i])]
  dt
}), fill=TRUE)

if (nrow(all_dt) == 0) stop("No results loaded after filtering.")

# enforce numeric & positivity for HR globally, then drop bad rows (very rare)
all_dt[, HR := suppressWarnings(as.numeric(HR))]
bad_n <- all_dt[!(is.finite(HR) & HR > 0), .N]
if (bad_n > 0) {
  warning(sprintf("Dropping %d rows with non-positive or non-finite HR", bad_n))
  all_dt <- all_dt[is.finite(HR) & HR > 0]
}


# ----------------------
# Summary (counts per cancer/type)
# ----------------------
summary_counts <- all_dt[, .(
  total_features = .N,
  sig_p05        = sum(p < 0.05, na.rm = TRUE),
  sig_FDR05      = sum(FDR < 0.05, na.rm = TRUE)
), by = .(cancer, data_type)]
print(summary_counts)

# ----------------------
# Prettify + transforms (order, winsor, clip)
# ----------------------
all_dt[, data_type := factor(data_type,
    levels=c("gene","iso_frac","iso_log",
             "mutation_baseline","mutation_combined","mutation_interaction"))]
# Define color palette
pal <- c(
  gene                 = "#0072B2",  # Blue
  iso_frac             = "#E69F00",  # Orange
  iso_log              = "#009E73",  # Green
  mutation_baseline    = "#D55E00",  # Vermilion
  mutation_combined    = "#CC79A7",  # Pink
  mutation_interaction = "#56B4E9"   # Light blue
)

# Per-group winsor on HR
all_dt[, HR_w := winsorize(HR, winsor_p), by=.(cancer, data_type)]

# Global clip for comparability
clip_q <- parse_clip(clip_str)
if (!is.null(clip_q)) {
  lo <- quantile(all_dt$HR_w, clip_q[1], na.rm=TRUE)
  hi <- quantile(all_dt$HR_w, clip_q[2], na.rm=TRUE)
  all_dt[HR_w < lo, HR_w := lo]
  all_dt[HR_w > hi, HR_w := hi]
}

# Order cancers
if (order_key == "medianHR") {
  med <- all_dt[, .(medianHR = median(HR_w, na.rm=TRUE)), by=cancer]
  ord <- med[order(medianHR)]$cancer
} else {
  ord <- sort(unique(all_dt$cancer))
}
all_dt[, cancer := factor(cancer, levels=ord)]

# Optional jitter sample per (cancer, data_type)
jitter_dt <- NULL
if (!is.na(jitter_n) && jitter_n > 0) {
  set.seed(1L)
  jitter_dt <- all_dt[, {
    n <- .N; if (n <= 0) .SD[0] else .SD[sample.int(n, min(n, jitter_n))]
  }, by=.(cancer, data_type)]
}

# ----------------------
# Plot: white boxes, colored borders, red outliers, optional jitter
# ----------------------
dodge_w <- 0.75
box_w   <- 0.45

base <- ggplot(all_dt, aes(x=cancer, y=HR_w)) +
  geom_hline(yintercept=1, color="red", linetype="dashed", linewidth=0.6) +
  labs(
    title = title_txt, x = "Cancer type",
    y = if (use_logy) "Hazard Ratio (log10 scale)" else "Hazard Ratio",
    color = "Data type"
  ) +
  theme_minimal(base_size=12) +
  theme(
    axis.text.x = element_text(angle=60, hjust=1, vjust=1),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p <- base + geom_boxplot( fill = "white",       # <- white fill for all boxes
    aes(color = data_type),                     # map border color to data type
    width = 0.65,                               # box width
    position = position_dodge(width = 0.75),    # group width
    linewidth = 0.6,                            # thicker border
    outlier.shape = 16, outlier.size = 0.6, outlier.alpha = 0.7, outlier.colour = "black" # colored outliers
    ) + scale_color_manual(values=pal)

if (!is.null(jitter_dt)) {
  p <- p +
    geom_jitter(
      data = jitter_dt,
      aes(x = cancer, y = HR_w, color = data_type),
      inherit.aes = FALSE,
      alpha = 0.22, size = 0.55,
      position = position_jitterdodge(jitter.width = jitter_w, dodge.width = dodge_w)
    )
}

if (use_logy) p <- p + scale_y_log10()

# Save main figure
dir.create(dirname(out_path), showWarnings=FALSE, recursive=TRUE)
ggsave(out_path, p, width=width_in, height=height_in, dpi=dpi_out, limitsize=FALSE)
message(sprintf("[SUMMARY BOXPLOT] Saved: %s", out_path))

# ----------------------
# Summary stats (per cancer/type): medians, IQR, counts
# ----------------------
summary_stats <- all_dt[, .(
  medianHR      = median(HR_w, na.rm = TRUE),
  median_logHR  = median(logHR, na.rm = TRUE),   # uses β if available, else log(HR)
  IQR           = IQR(HR_w, na.rm = TRUE),
  lowerQ        = quantile(HR_w, 0.25, na.rm = TRUE),
  upperQ        = quantile(HR_w, 0.75, na.rm = TRUE),
  total_features = .N,
  sig_p05        = sum(p < 0.05, na.rm = TRUE),
  sig_FDR05      = sum(FDR < 0.05, na.rm = TRUE)
), by = .(cancer, data_type)]

summary_out <- sub("\\.[^.]+$", "_summary.tsv", out_path)
fwrite(summary_stats, summary_out, sep="\t")
message(sprintf("[SUMMARY] Saved: %s", summary_out))













RUN_ALL <- FALSE   # set to TRUE to generate all plots below

if (RUN_ALL) {
# ----------------------
# Heatmap: median log(HR) (symmetric around 0)
# ----------------------
lim <- max(abs(range(summary_stats$median_logHR, na.rm = TRUE)))
heat <- ggplot(summary_stats, aes(x = cancer, y = data_type, fill = median_logHR)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2b83ba", mid = "white", high = "#d7191c",
    midpoint = 0, limits = c(-lim, lim), name = "median log(HR)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Median log(Hazard Ratio) per Cancer and Data Type", x = "Cancer Type", y = "Data Type")

heat_out <- sub("\\.[^.]+$", "_logHR_heatmap.png", out_path)
ggsave(heat_out, heat, width = 12, height = 4, dpi = 300)
message(sprintf("[SUMMARY HEATMAP] Saved: %s", heat_out))

# ----------------------
# Barplot: significant features (FDR < 0.05)
# ----------------------
bar <- ggplot(summary_stats, aes(x = cancer, y = sig_FDR05, fill = data_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = pal, name = "Data type") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "top"
  ) +
  labs(title = "Significant Features (FDR < 0.05)", x = "Cancer Type", y = "Count")

bar_out <- sub("\\.[^.]+$", "_sigFDR_barplot.png", out_path)
ggsave(bar_out, bar, width = 12, height = 4, dpi = 300)
message(sprintf("[SUMMARY BARPLOT (FDR)] Saved: %s", bar_out))


# ----------------------
# Barplot: significant features (p < 0.05)
# ----------------------
bar <- ggplot(summary_stats, aes(x = cancer, y = sig_p05, fill = data_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = pal, name = "Data type") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "top"
  ) +
  labs(title = "Significant Features (p < 0.05)", x = "Cancer Type", y = "Count")

bar_out <- sub("\\.[^.]+$", "_sigp05_barplot.png", out_path)
ggsave(bar_out, bar, width = 12, height = 4, dpi = 300)
message(sprintf("[SUMMARY BARPLOT (p-value)] Saved: %s", bar_out))



# ----------------------
# P-value histograms
# ----------------------
# Keep finite p-values in (0,1]
p_dt <- all_dt[is.finite(p) & p > 0 & p <= 1]

# A) One histogram per data type (aggregated over cancers)
p_hist1 <- ggplot(p_dt, aes(x = p, fill = data_type)) +
  geom_histogram(bins = 40, position = "identity", alpha = 0.5) +
  geom_hline(yintercept = nrow(p_dt)/40, linetype = "dotted", linewidth = 0.4) +  # rough uniform ref
  scale_fill_manual(values = pal, name = "Data Type") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", panel.grid.minor = element_blank()) +
  labs(title = "P-value histograms by data type (all cancers combined)",
       x = "p-value", y = "Count", fill = "Data type")
p_hist1_out <- sub("\\.[^.]+$", "_pHist_datatype.png", out_path)
ggsave(p_hist1_out, p_hist1, width = 10, height = 5, dpi = 300)
message(sprintf("[P-VALUES] Saved: %s", p_hist1_out))

# B) Faceted by cancer and data type (compact): one PDF, many panels
p_hist2 <- ggplot(p_dt, aes(x = p)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  facet_grid(rows = vars(data_type), cols = vars(cancer), scales = "free_y") +
  theme_minimal(base_size = 10) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "P-value histograms per cancer and data type",
       x = "p-value", y = "Count")
p_hist2_out <- sub("\\.[^.]+$", "_pHist_faceted.pdf", out_path)
ggsave(p_hist2_out, p_hist2, width = 22, height = 6, dpi = 300, limitsize = FALSE)
message(sprintf("[P-VALUES] Saved: %s", p_hist2_out))
# ----------------------
} # end RUN_ALL


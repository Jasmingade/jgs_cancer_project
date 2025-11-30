#!/usr/bin/env Rscript
# ============================================================
# 03a_plot_PH_facets.R
# ------------------------------------------------------------
# Combines all *_PH_summary.csv outputs into one dataset
# and generates a single faceted PDF per data type (gene, iso_log, iso_frac)
# showing proportional hazards violation rates.
# ============================================================

#!/usr/bin/env Rscript
# ============================================================
# 03a_investigate_PH_assumption.R
# ------------------------------------------------------------
# Investigates proportional hazards (PH) assumption across
# features (gene/iso_log/iso_frac) within each cancer type.
# Fits unstratified Cox models and tests PH assumption using cox.zph().
#
# Outputs:
#   - *_per_feature.csv : p-values for PH tests per feature
#   - *_PH_summary.csv  : summary of PH violations
#   - *_PH_violations.pdf/png : visual summary plots
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(yaml)
  library(parallel)
  library(ggplot2)
  library(reshape2)
  library(scales)
})

say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }

# ============================================================
# Arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5)
  die("Usage: 03a_investigate_PH_assumption.R <expr_norm.csv> <manifest.csv> <covariates.yaml> <out_summary.csv> <out_plot.pdf>")

expr_in <- args[[1]]
mani_in <- args[[2]]
cov_yaml <- args[[3]]
out_csv  <- args[[4]]
out_pdf  <- args[[5]]

# Ensure correct PDF extension
if (!grepl("\\.pdf$", out_pdf)) {
  out_pdf <- sub("\\.txt$", ".pdf", out_pdf)
}

say("=== Investigating PH assumption ===")
say("Expr: %s", expr_in)
say("Manifest: %s", mani_in)
say("Covariates: %s", cov_yaml)

# ============================================================
# Load input data
# ============================================================
expr <- fread(expr_in)
mani <- fread(mani_in)
cfg  <- yaml::read_yaml(cov_yaml)

feat_col <- if ("feature" %in% names(expr)) "feature" else names(expr)[1]
features <- expr[[feat_col]]
expr[[feat_col]] <- NULL
expr <- as.matrix(expr)
rownames(expr) <- features

cancer <- sub("^.*TCGA_([A-Z0-9]+)_.*$", "\\1", basename(expr_in))
data_type <- sub("^.*TCGA_[A-Z0-9]+_(.*?)\\.normalized\\.csv$", "\\1", basename(expr_in))
say("[INFO] Cancer=%s | DataType=%s", cancer, data_type)

# ============================================================
# Iso_frac logit transform
# ============================================================
if (identical(data_type, "iso_frac")) {
  eps <- 1e-6
  expr <- pmin(pmax(expr, eps), 1 - eps)
  expr <- log(expr / (1 - expr))
}

# ============================================================
# Align expression and survival data
# ============================================================
mani <- mani[complete.cases(mani[, .(OS_time, OS_event)]), ]
common <- intersect(colnames(expr), mani$case_id)
if (length(common) < 20) die("Too few overlapping samples.")
expr <- expr[, common, drop = FALSE]
mani <- mani[match(common, mani$case_id)]
y <- with(mani, Surv(OS_time, OS_event))

# ============================================================
# Covariate setup
# ============================================================
covariates <- cfg$baseline_covariates
if (!is.null(cfg$conditional_covariates)) {
  for (cond in cfg$conditional_covariates) {
    nm <- cond$name
    thr <- cond$coverage_threshold
    if (nm %in% names(mani) && mean(!is.na(mani[[nm]])) >= thr)
      covariates <- c(covariates, nm)
  }
}
covariates <- unique(covariates)
Xcov <- data.frame(row.names = mani$case_id)
for (c in covariates) {
  if (c %in% names(mani)) Xcov[[c]] <- mani[[c]]
}
Xcov <- Xcov[, sapply(Xcov, function(x) length(unique(na.omit(x))) > 1), drop = FALSE]
say("[INFO] Using covariates: %s", paste(names(Xcov), collapse=", "))

# ============================================================
# Sample subset of features for PH test
# ============================================================
n_features <- nrow(expr)
n_sample <- min(300, n_features)
set.seed(1)
sample_idx <- sample(seq_len(n_features), n_sample)
expr <- expr[sample_idx, , drop = FALSE]
features <- features[sample_idx]
say("[INFO] Testing PH assumption on %d randomly sampled features.", n_sample)

# ============================================================
# Run Cox + PH assumption tests
# ============================================================
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "6"))
say("[INFO] Using %d cores for CoxPH fits.", ncores)

ph_results <- mclapply(seq_len(nrow(expr)), function(i) {
  tryCatch({
    vals <- as.numeric(expr[i, ])
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0) return(NULL)
    vals_z <- scale(vals)
    df <- cbind.data.frame(expr = vals_z, Xcov, y = y)
    rhs <- paste(c("expr", names(Xcov)), collapse = " + ")
    fml <- as.formula(paste("y ~", rhs))
    
    fit <- suppressWarnings(try(coxph(fml, data = df, ties = "efron"), silent = TRUE))
    if (inherits(fit, "try-error") || isFALSE(fit$converged)) return(NULL)
    
    ph_test <- suppressWarnings(try(cox.zph(fit), silent = TRUE))
    if (inherits(ph_test, "try-error")) return(NULL)
    tab <- as.data.frame(ph_test$table)
    
    covars_in_model <- intersect(names(Xcov), rownames(tab))
    out <- data.frame(feature = features[i], GLOBAL_p = tab["GLOBAL", "p"])
    for (cv in covars_in_model) {
      out[[paste0(cv, "_p")]] <- tab[cv, "p"]
    }
    out
  }, error = function(e) NULL)
}, mc.cores = ncores)

ph_results <- Filter(Negate(is.null), ph_results)
if (length(ph_results) == 0) die("No valid PH tests completed.")
ph_dt <- rbindlist(ph_results, fill = TRUE)

# ============================================================
# Summarize PH violations
# ============================================================
cols <- setdiff(names(ph_dt), "feature")
ph_summary <- data.table(cancer = cancer, data_type = data_type)
for (col in cols) {
  ph_summary[[paste0(col, "_viol")]] <- mean(ph_dt[[col]] < 0.05, na.rm = TRUE)
}

# Save results
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
fwrite(ph_dt, gsub(".csv$", "_per_feature.csv", out_csv))
fwrite(ph_summary, gsub(".csv$", "_PH_summary.csv", out_csv))
say("[DONE] Saved PH violation summary to %s", out_csv)

# ============================================================
# Prepare data for plotting
# ============================================================
ph_melt <- melt(ph_summary, id.vars = c("cancer", "data_type"))
setDT(ph_melt)
ph_melt[, variable := gsub("_p_viol|_viol|_p", "", variable)]

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(scales)
})

say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }

# ============================================================
# Arguments
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  die("Usage: 03b_plot_PH_facets.R <input_dir> <output_pdf>")

in_dir <- args[[1]]
out_pdf <- args[[2]]
if (!grepl("\\.pdf$", out_pdf)) out_pdf <- paste0(out_pdf, ".pdf")

say("=== Combining PH summaries into facets ===")
say("Input directory: %s", in_dir)
say("Output PDF: %s", out_pdf)

# ============================================================
# Load all *_PH_summary.csv files
# ============================================================
files <- list.files(in_dir, pattern = "_PH_summary\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) die("No *_PH_summary.csv files found in input directory.")

say("[INFO] Found %d summary files", length(files))

all_dt <- rbindlist(lapply(files, function(f) {
  dt <- fread(f)
  dt$cancer <- if ("cancer" %in% names(dt)) dt$cancer else sub("^.*TCGA_([A-Z0-9]+)_.*$", "\\1", f)
  dt$data_type <- if ("data_type" %in% names(dt)) dt$data_type else
    if (grepl("iso_frac", f)) "iso_frac" else if (grepl("iso_log", f)) "iso_log" else "gene"
  dt
}), fill = TRUE)

if (nrow(all_dt) == 0) die("No data loaded from summaries.")

# ============================================================
# Melt and clean variable names
# ============================================================
ph_melt <- melt(all_dt, id.vars = c("cancer", "data_type"))
setDT(ph_melt)
ph_melt[, variable := gsub("_p_viol|_viol|_p", "", variable)]
ph_melt[, cancer := factor(cancer, levels = sort(unique(cancer)))]

# ============================================================
# Create faceted plot
# ============================================================
p <- ggplot(ph_melt, aes(x = variable, y = value, fill = variable)) +
  geom_col(color = "black", alpha = 0.8) +
  facet_wrap(~ cancer, ncol = 6, scales = "free_y") +
  geom_text(aes(label = sprintf("%.1f%%", 100 * value)), vjust = -0.4, size = 2.5) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 9)
  ) +
  labs(
    title = "Proportional Hazards (PH) Assumption Violations Across Cancers",
    subtitle = "Proportion of features violating PH assumption (p < 0.05)",
    x = "Variable",
    y = "Proportion of features",
    fill = "Covariate"
  )

# ============================================================
# Output faceted PDF (and optional PNG)
# ============================================================
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(out_pdf, p, width = 12, height = 9)
png_out <- sub("\\.pdf$", ".png", out_pdf)
ggsave(png_out, p, width = 12, height = 9, dpi = 300)

say("[DONE] Saved faceted PH violation plot:")
say("   %s", out_pdf)
say("   %s", png_out)
say("=== PH facet plot generation complete ✓ ===")

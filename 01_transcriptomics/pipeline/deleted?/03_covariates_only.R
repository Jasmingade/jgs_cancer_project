#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table); library(survival); library(yaml); library(ggplot2) })

`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
say <- function(...) message(sprintf(...))
die <- function(...) { message(sprintf(...)); quit(status = 1) }

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: 03_covariates_only.R <manifest.csv> <covariates_yaml> <out_tsv> <out_plot_png>")
}
mani_in  <- args[[1]]
cov_yaml <- args[[2]]
out_tsv  <- args[[3]]
out_png  <- args[[4]]

# --- load
mani <- fread(mani_in)
cfg  <- yaml::read_yaml(cov_yaml)

if (!all(c("OS_time","OS_event") %in% names(mani))) stop("Manifest missing OS_time/OS_event.")
covariates <- unique(c(cfg$baseline_covariates %||% c("age","sex","stage"),
                       vapply(cfg$conditional_covariates %||% list(), `[[`, "", "name")))
covariates <- covariates[covariates %in% names(mani)]
if (!length(covariates)) stop("No covariates found in manifest.")

# clean & align
mani <- mani[complete.cases(mani[, .(OS_time, OS_event)]), ]
mani <- mani[complete.cases(mani[, ..covariates])]
y <- with(mani, Surv(OS_time, OS_event))

# fit covariates-only model
form <- as.formula(paste("y ~", paste(covariates, collapse = " + ")))
fit  <- coxph(form, data = mani, x = TRUE)
s    <- summary(fit)

# tidy
co <- as.data.table(coef(summary(fit)), keep.rownames = "covariate")
setnames(co, c("coef","exp(coef)","se(coef)","z","Pr(>|z|)"),
              c("beta","HR","se","z","p"))
# CIs
ci <- confint(fit)
co[, `:=`(HR_lo = exp(ci[,1]), HR_hi = exp(ci[,2]), FDR = p.adjust(p, "BH"), logHR = beta)]

# c-index
cidx <- tryCatch({ survival::concordance(fit)$concordance }, error = function(e) NA_real_)
co[, cindex := cidx]

# write
fwrite(co, out_tsv)

# forest plot
pal <- "#2c7fb8"
p <- ggplot(co, aes(x = reorder(covariate, HR), y = HR)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_pointrange(aes(ymin = HR_lo, ymax = HR_hi), color = pal) +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(title = "Covariate-only Cox model", x = "Covariate", y = "Hazard Ratio (HR)")
ggsave(out_png, p, width = 7, height = 5, dpi = 300)

message(sprintf("[COVARIATES] Saved: %s and %s", out_tsv, out_png))

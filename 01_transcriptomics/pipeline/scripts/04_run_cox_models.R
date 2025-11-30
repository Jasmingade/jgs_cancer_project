#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Usage: 04_run_cox_models.R <TCGA_CANCER> <DTYPE> <TX2GENE> <OUT_BASE>")
}

cancer     <- args[1]   # "TCGA_ACC"
dtype      <- args[2]   # "gene", "iso_log", "iso_frac"
tx2gene_fp <- args[3]   # path to tx2gene mapping
out_base   <- args[4]

# ============================================================
# Read toggles from environment variables
# ============================================================
to_bool <- function(x) tolower(x) %in% c("1","true","yes")

RUN_MODEL1_EXPR_ONLY      <- to_bool(Sys.getenv("RUN_MODEL1_EXPR_ONLY", "false"))
RUN_MODEL2_MUT_ONLY       <- to_bool(Sys.getenv("RUN_MODEL2_MUT_ONLY","false"))
RUN_MODEL3_EXPR_PLUS_MUT  <- to_bool(Sys.getenv("RUN_MODEL3_EXPR_PLUS_MUT","false"))
RUN_MODEL4_ISO_MUT_INTERACT <- to_bool(Sys.getenv("RUN_MODEL4_ISO_MUT_INTERACT","false"))

cat("[INFO] Model toggles:\n")
cat("  M1 expr-only: ", RUN_MODEL1_EXPR_ONLY, "\n")
cat("  M2 mut-only:  ", RUN_MODEL2_MUT_ONLY, "\n")
cat("  M3 expr+mut:  ", RUN_MODEL3_EXPR_PLUS_MUT, "\n")
cat("  M4 iso*mut:   ", RUN_MODEL4_ISO_MUT_INTERACT, "\n")

# ============================================================
# Paths
# ============================================================
clin_file <- sprintf("01_transcriptomics/data/clinical/%s_clinical.csv", cancer)

expr_file <- sprintf("01_transcriptomics/out/02_norm/%s_%s.normalized.csv",
                     cancer, dtype)

mut_file_any <- sprintf(
  "01_transcriptomics/out/03_mutation/%s/gene/%s_gene_ensembl_coding_any.csv",
  cancer, cancer
)

out_dir <- file.path(out_base, cancer)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# LOAD CLINICAL
# ============================================================
clin <- fread(clin_file)
clin[, sex := factor(sex)]
clin[, stage := factor(stage)]
if ("subtype" %in% names(clin)) clin[, subtype := factor(subtype)]
setkey(clin, case_id)

covars <- c("age","sex","stage")
if ("subtype" %in% names(clin)) covars <- c(covars, "subtype")

base_form <- function(extra_terms) {
  as.formula(paste(
    "Surv(OS_time, OS_event) ~",
    paste(c(extra_terms, covars), collapse=" + ")
  ))
}

# ============================================================
# LOAD EXPRESSION
# ============================================================
expr <- fread(expr_file)
setnames(expr, 1, "feature_id")
cases_expr <- setdiff(names(expr), "feature_id")

expr_mat <- as.matrix(expr[, ..cases_expr])
rownames(expr_mat) <- expr$feature_id

# ============================================================
# LOAD MUTATION (coding_any only)
# ============================================================
mut_any <- NULL
if (file.exists(mut_file_any)) {
  mut_any <- fread(mut_file_any)
  setnames(mut_any, 1, "feature_id")
  mut_cases <- setdiff(names(mut_any), "feature_id")
  mut_any_mat <- as.matrix(mut_any[, ..mut_cases])
  rownames(mut_any_mat) <- mut_any$feature_id
} else {
  mut_any_mat <- NULL
}

# ============================================================
# ALIGN CASE IDS
# ============================================================
cases <- Reduce(intersect, list(
  clin$case_id,
  cases_expr,
  if (!is.null(mut_any_mat)) mut_cases else cases_expr
))

cases <- sort(cases)
if (length(cases) < 20) stop("Too few overlapping cases")

# Build base clinical frame
make_base <- function(cases) data.table(
  case_id  = cases,
  OS_time  = clin[match(cases, case_id), OS_time],
  OS_event = clin[match(cases, case_id), OS_event],
  age      = clin[match(cases, case_id), age],
  sex      = clin[match(cases, case_id), sex],
  stage    = clin[match(cases, case_id), stage],
  subtype  = if ("subtype" %in% names(clin)) clin[match(cases, case_id), subtype] else NA
)

base <- make_base(cases)

# subset matrices to aligned cases
expr_mat <- expr_mat[, cases, drop=FALSE]
if (!is.null(mut_any_mat)) mut_any_mat <- mut_any_mat[, cases, drop=FALSE]

# ============================================================
# MODEL 1 — Expression only
# ============================================================
if (RUN_MODEL1_EXPR_ONLY) {
  out_file <- sprintf("%s/model1_expr_only_%s.csv", out_dir, dtype)
  res <- list()

  form <- base_form("expr")

  for (i in seq_len(nrow(expr_mat))) {
    feat <- rownames(expr_mat)[i]
    x <- expr_mat[i, ]

    df <- copy(base); df[, expr := x]

    fit <- try(coxph(form, df), TRUE)
    if (inherits(fit, "try-error")) next

    ci <- suppressWarnings(confint(fit)["expr", ])
    res[[length(res)+1]] <- data.table(
      cancer=cancer, dtype=dtype, feature_id=feat,
      HR=exp(coef(fit)["expr"]), HR_lo=exp(ci[1]), HR_hi=exp(ci[2]),
      pval=summary(fit)$coef["expr","Pr(>|z|)"]
    )
  }
  fwrite(rbindlist(res), out_file)
}

# ============================================================
# MODEL 2 — Mutation only
# ============================================================
if (RUN_MODEL2_MUT_ONLY && !is.null(mut_any_mat)) {
  out_file <- sprintf("%s/model2_mut_only.csv", out_dir)
  res <- list()

  form <- base_form("mut")

  for (feat in rownames(mut_any_mat)) {
    z <- mut_any_mat[feat, ]
    if (sum(z) < 5) next

    df <- copy(base); df[, mut := z]

    fit <- try(coxph(form, df), TRUE)
    if (!inherits(fit, "try-error")) {
      ci <- suppressWarnings(confint(fit)["mut", ])
      res[[length(res)+1]] <- data.table(
        cancer=cancer, feature_id=feat,
        HR=exp(coef(fit)["mut"]), HR_lo=exp(ci[1]), HR_hi=exp(ci[2]),
        pval=summary(fit)$coef["mut","Pr(>|z|)"]
      )
    }
  }

  fwrite(rbindlist(res), out_file)
}

# ============================================================
# MODEL 3 — Expression + Mutation
# ============================================================
if (RUN_MODEL3_EXPR_PLUS_MUT && !is.null(mut_any_mat)) {
  out_file <- sprintf("%s/model3_expr_plus_mut.csv", out_dir)

  common <- intersect(rownames(expr_mat), rownames(mut_any_mat))
  res <- list()

  form <- base_form(c("expr","mut"))

  for (feat in common) {
    x <- expr_mat[feat, ]
    z <- mut_any_mat[feat, ]
    if (sum(z) < 5) next

    df <- copy(base); df[,expr:=x]; df[,mut:=z]

    fit <- try(coxph(form, df), TRUE)
    if (!inherits(fit,"try-error")) {
      s <- summary(fit)$coefficients
      res[[length(res)+1]] <- data.table(
        cancer=cancer, feature_id=feat,
        HR_expr=exp(s["expr","coef"]), p_expr=s["expr","Pr(>|z|)"],
        HR_mut=exp(s["mut","coef"]),   p_mut=s["mut","Pr(>|z|)"]
      )
    }
  }

  fwrite(rbindlist(res), out_file)
}

# ============================================================
# MODEL 4 — Isoform × Gene Mutation Interaction
# ============================================================
if (RUN_MODEL4_ISO_MUT_INTERACT && !is.null(mut_any_mat) && file.exists(tx2gene_fp)) {

  out_file <- sprintf("%s/model4_iso_mut_interaction.csv", out_dir)

  # Load mapping
  tx <- fread(tx2gene_fp)
  tx <- tx[feature_id %in% rownames(expr_mat)]
  tx <- tx[gene_id %in% rownames(mut_any_mat)]

  form <- base_form(c("expr","mut","expr:mut"))
  res <- list()

  for (i in seq_len(nrow(tx))) {
    iso  <- tx$feature_id[i]
    gene <- tx$gene_id[i]

    x <- expr_mat[iso, ]
    z <- mut_any_mat[gene, ]

    if (sum(z) < 5) next

    df <- copy(base); df[,expr:=x]; df[,mut:=z]

    fit <- try(coxph(form, df), TRUE)
    if (!inherits(fit,"try-error")) {
      res[[length(res)+1]] <- data.table(
        cancer=cancer, isoform_id=iso, gene_id=gene,
        HR_int=exp(coef(fit)["expr:mut"]),
        p_int=summary(fit)$coef["expr:mut","Pr(>|z|)"]
      )
    }
  }

  fwrite(rbindlist(res), out_file)
}

cat("[DONE] All selected models finished for", cancer, dtype, "\n")

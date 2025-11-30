#!/usr/bin/env Rscript
# ============================================================
# 03e_gene_expr_mut.R  (updated)
# ------------------------------------------------------------
# Gene- or isoform-based expression + mutation Cox model:
#
#   Surv(OS_time, OS_event) ~ expr + mut + covariates
#
# For each cancer + mutation group:
#   - expr comes from an expression matrix:
#       * rows can be ENSG (genes) OR ENST (isoforms)
#       * ENST are mapped to ENSG with tx2gene.csv
#   - mut comes from a gene-level mutation matrix (ENSG)
#
# Matching is done on core ENSG IDs (version stripped).
#
# Optional 7th arg = gene list file (one ENSG per line, WITHOUT
# version) to restrict to specific genes. If omitted → all
# overlapping genes (by core ENSG) are tested.
#
# Env:
#   TX2GENE_FILE can override the default path to tx2gene.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(parallel)
  library(yaml)
})

say           <- function(...) message(sprintf(...))
die           <- function(...) { message(sprintf(...)); quit(status = 1) }
strip_version <- function(x) sub("\\.\\d+$", "", x)

# ------------------------------------------------------------
# tx2gene mapping for ENST → ENSG (optional)
# ------------------------------------------------------------
tx2gene_file <- Sys.getenv("TX2GENE_FILE",
                           "01_transcriptomics/data/raw/tx2gene.csv")
tx2gene_map  <- NULL

if (file.exists(tx2gene_file)) {
  tx2gene <- fread(tx2gene_file)
  names(tx2gene) <- tolower(names(tx2gene))

  tx_col <- intersect(
    c("tx_id", "transcript_id", "ensembl_transcript_id"),
    names(tx2gene)
  )[1]
  gene_col <- intersect(
    c("gene_id", "ensembl_gene_id", "gene"),
    names(tx2gene)
  )[1]

  if (is.na(tx_col) || is.na(gene_col)) {
    warning("[WARN] tx2gene file found but no suitable tx/gene columns — ENST will NOT be mapped.")
  } else {
    tx2gene_map <- unique(tx2gene[, .(
      tx_core   = strip_version(get(tx_col)),
      gene_core = strip_version(get(gene_col))
    )])
    setkey(tx2gene_map, tx_core)
    say("[INFO] Loaded tx2gene mapping in 03e: %d transcript→gene rows",
        nrow(tx2gene_map))
  }
} else {
  say(sprintf("[INFO] No tx2gene mapping at %s; assuming expression features are already gene-level ENSG.",
              tx2gene_file))
}

# ============================================================
# Args
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6 || length(args) > 7) {
  die(paste(
    "Usage:",
    "03e_gene_expr_mut.R <expr_matrix.csv> <mut_gene.csv> <manifest.csv>",
    "<covariates.yaml> <out_results.csv> <out_summary.txt>",
    "[optional_gene_list.txt]"
  ))
}

expr_in   <- args[[1]]  # rows: ENSG or ENST
mut_in    <- args[[2]]  # rows: ENSG
mani_in   <- args[[3]]
cov_yaml  <- args[[4]]
out_res   <- args[[5]]
out_sum   <- args[[6]]
gene_list <- if (length(args) == 7) args[[7]] else NA_character_

dir.create(dirname(out_res), recursive = TRUE, showWarnings = FALSE)

say("=== 03e_gene_expr_mut (expr + mut, gene-mapped) ===")
say("[INFO] expr_in   = %s", expr_in)
say("[INFO] mut_in    = %s", mut_in)
say("[INFO] mani_in   = %s", mani_in)
say("[INFO] cov_yaml  = %s", cov_yaml)
if (!is.na(gene_list)) say("[INFO] gene_list = %s", gene_list)

# Parse cancer + mutation group from filename (like 03b)
cancer    <- sub("^.*TCGA_([A-Z0-9]+).*", "\\1", basename(mut_in))
mut_group <- sub("^.*ensembl_(.*)\\.csv$", "\\1", basename(mut_in))
if (mut_group == basename(mut_in)) mut_group <- "unknown"

say("[INFO] Cancer = %s | MutationGroup = %s", cancer, mut_group)

# ============================================================
# Load expression (rows: ENSG or ENST)
# ============================================================
expr_dt <- fread(expr_in)
feat_col <- if ("feature" %in% names(expr_dt)) "feature" else names(expr_dt)[1]
features_expr <- expr_dt[[feat_col]]
expr_dt[[feat_col]] <- NULL

expr_mat <- as.matrix(expr_dt)
rownames(expr_mat) <- features_expr
case_expr <- colnames(expr_mat)

say("[INFO] Expression matrix: %d features × %d samples",
    nrow(expr_mat), ncol(expr_mat))

# ============================================================
# Load mutation (gene-level ENSG)
# ============================================================
mut_dt <- fread(mut_in)
if (!"feature_id" %in% names(mut_dt)) {
  die("Mutation file must have feature_id as first column.")
}

features_mut <- mut_dt$feature_id
case_mut     <- setdiff(names(mut_dt), "feature_id")

mut_mat <- as.matrix(mut_dt[, ..case_mut])
rownames(mut_mat) <- features_mut

# convert to 0/1
mut_mat[is.na(mut_mat)] <- 0
mut_mat <- (mut_mat > 0) * 1

say("[INFO] Mutation matrix: %d genes × %d samples",
    nrow(mut_mat), ncol(mut_mat))

# ============================================================
# Load manifest
# ============================================================
mani <- fread(mani_in)
if (!all(c("case_id","OS_time","OS_event") %in% names(mani))) {
  die("Manifest must contain: case_id, OS_time, OS_event")
}

# ============================================================
# Align case IDs
# ============================================================
common_samples <- Reduce(
  intersect,
  list(case_expr, case_mut, mani$case_id)
)

if (length(common_samples) < 30) {
  die(sprintf("Too few overlapping samples between expr/mut/manifest: %d",
              length(common_samples)))
}

common_samples <- sort(common_samples)

expr_mat <- expr_mat[, common_samples, drop = FALSE]
mut_mat  <- mut_mat[,  common_samples, drop = FALSE]
mani     <- mani[match(common_samples, mani$case_id)]

y <- with(mani, Surv(OS_time, OS_event))

say("[INFO] Aligned samples: %d", length(common_samples))

# ============================================================
# Covariates (same logic as 03c)
# ============================================================
cfg <- yaml::read_yaml(cov_yaml)

covariates <- cfg$baseline_covariates
if (!is.null(cfg$conditional_covariates)) {
  for (cond in cfg$conditional_covariates) {
    nm  <- cond$name
    thr <- cond$coverage_threshold
    if (nm %in% names(mani) && mean(!is.na(mani[[nm]])) >= thr) {
      covariates <- c(covariates, nm)
    }
  }
}

covariates <- unique(covariates)
covariates <- covariates[covariates %in% names(mani)]

Xcov <- mani[, covariates, with = FALSE]

# drop covariates with no variability
Xcov <- Xcov[, lapply(.SD, function(x) {
  if (length(unique(na.omit(x))) > 1) x else NULL
})]

cov_used <- names(Xcov)
say("[INFO] Covariates used after filtering: %s",
    if (length(cov_used)) paste(cov_used, collapse = ", ") else "<none>")

# Build formula: y ~ expr + mut + cov...
if (length(cov_used) > 0) {
  cov_terms  <- paste(cov_used, collapse = " + ")
  fml_string <- paste("y ~ expr + mut +", cov_terms)
} else {
  fml_string <- "y ~ expr + mut"
}

fml_string <- gsub("\\+\\s*$", "", fml_string)
fml <- as.formula(fml_string)
say("[INFO] Cox formula: %s", fml_string)

# ============================================================
# Candidate genes (map expr features to core ENSG)
# ============================================================
genes_expr_raw <- rownames(expr_mat)   # ENSG or ENST (+version)
genes_mut_raw  <- rownames(mut_mat)    # ENSG (+version)

# start with version-stripped IDs
genes_expr_core0 <- strip_version(genes_expr_raw)
genes_mut_core   <- strip_version(genes_mut_raw)

# if tx2gene is available, map ENST → ENSG; otherwise assume ENSG
if (!is.null(tx2gene_map)) {
  map_dt <- data.table(
    tx_core  = genes_expr_core0,
    expr_row = genes_expr_raw
  )
  # join: will fill gene_core where tx_core is in tx2gene_map
  map_dt <- tx2gene_map[map_dt, on = "tx_core"]
  # for rows without mapping, keep tx_core (assume these are ENSG)
  map_dt[is.na(gene_core), gene_core := tx_core]

  genes_expr_core <- map_dt$gene_core
} else {
  genes_expr_core <- genes_expr_core0
}

# Map from core ENSG → one representative rowname in each matrix
map_expr <- tapply(genes_expr_raw, genes_expr_core, `[`, 1)
map_mut  <- tapply(genes_mut_raw,  genes_mut_core,  `[`, 1)

cores_both <- intersect(names(map_expr), names(map_mut))

if (!is.na(gene_list)) {
  gl_core   <- fread(gene_list, header = FALSE)$V1
  genes_core <- intersect(cores_both, gl_core)
  say("[INFO] Restricting to %d core genes from gene_list", length(genes_core))
} else {
  genes_core <- cores_both
}

if (length(genes_core) == 0) {
  die("No overlapping genes between expr, mut and optional gene_list.")
}

genes_expr_use <- unname(map_expr[genes_core])
genes_mut_use  <- unname(map_mut[genes_core])

valid_idx <- which(!is.na(genes_expr_use) & !is.na(genes_mut_use))
genes_core     <- genes_core[valid_idx]
genes_expr_use <- genes_expr_use[valid_idx]
genes_mut_use  <- genes_mut_use[valid_idx]

if (length(genes_core) == 0) {
  die("No valid gene mappings after resolving core IDs.")
}

say("[INFO] Number of genes to test (core ENSG): %d", length(genes_core))

# ============================================================
# Run gene-matched Cox models
# ============================================================
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "6"))
if (is.na(ncores) || ncores < 1) ncores <- 2
say("[INFO] Using %d cores", ncores)

res_list <- mclapply(seq_along(genes_core), function(i) {
  tryCatch({
    g_core    <- genes_core[i]
    g_expr_id <- genes_expr_use[i]
    g_mut_id  <- genes_mut_use[i]

    expr_vec <- as.numeric(expr_mat[g_expr_id, ])
    mut_vec  <- as.numeric(mut_mat[g_mut_id, ])

    # basic QC: expression must vary
    if (all(is.na(expr_vec)) || sd(expr_vec, na.rm = TRUE) == 0)
      return(data.frame())

    # mutation must be 0/1 and vary
    mut_vec[is.na(mut_vec)] <- 0
    mut_vec <- ifelse(mut_vec > 0, 1L, 0L)

    if (length(unique(mut_vec)) < 2)
      return(data.frame())

    # events in each group
    events_mut <- sum(mani$OS_event[mut_vec == 1] == 1, na.rm = TRUE)
    events_wt  <- sum(mani$OS_event[mut_vec == 0] == 1, na.rm = TRUE)

    # safety: require at least 3 events in each group
    if (events_mut < 3 || events_wt < 3)
      return(data.frame())

    df <- data.frame(
      expr = as.numeric(scale(expr_vec)),
      mut  = mut_vec,
      Xcov,
      y
    )

    fit <- suppressWarnings(try(coxph(fml, data = df, ties = "efron"), silent = TRUE))
    if (inherits(fit, "try-error") || isFALSE(fit$converged))
      return(data.frame())

    s  <- summary(fit)
    co <- s$coef

    get_row <- function(term) {
      if (term %in% rownames(co)) co[term, ] else rep(NA_real_, ncol(co))
    }

    expr_row <- get_row("expr")
    mut_row  <- get_row("mut")

    if (all(is.na(expr_row)) || all(is.na(mut_row)))
      return(data.frame())

    z_expr        <- expr_row["z"]
    wald_expr     <- z_expr^2
    delta_LL_expr <- wald_expr / 2

    z_mut        <- mut_row["z"]
    wald_mut     <- z_mut^2
    delta_LL_mut <- wald_mut / 2

    data.frame(
      feature        = g_core,      # core ENSG (no version)
      feature_expr   = g_expr_id,   # rowname in expr_mat (ENSG or ENST)
      feature_mut    = g_mut_id,    # rowname in mut_mat (ENSG)

      beta_expr      = expr_row["coef"],
      HR_expr        = exp(expr_row["coef"]),
      p_expr         = expr_row["Pr(>|z|)"],

      beta_mut       = mut_row["coef"],
      HR_mut         = exp(mut_row["coef"]),
      p_mut          = mut_row["Pr(>|z|)"],

      delta_LL_expr  = delta_LL_expr,
      delta_LL_mut   = delta_LL_mut
    )
  }, error = function(e) data.frame())
}, mc.cores = ncores)

res_list <- Filter(function(x) is.data.frame(x) && nrow(x) > 0, res_list)

run_id <- sprintf("[03e][%s][%s]", cancer, mut_group)

if (length(res_list) == 0) {
  say("%s [WARN] No valid expr+mut models fitted. Writing empty output.",
      run_id)

  res_empty <- data.table(
    feature        = character(),
    feature_expr   = character(),
    feature_mut    = character(),
    beta_expr      = numeric(),
    HR_expr        = numeric(),
    p_expr         = numeric(),
    beta_mut       = numeric(),
    HR_mut         = numeric(),
    p_mut          = numeric(),
    delta_LL_expr  = numeric(),
    delta_LL_mut   = numeric(),
    FDR_expr       = numeric(),
    FDR_mut        = numeric(),
    cancer         = character(),
    mut_group      = character()
  )

  out_res_full <- sub("\\.cox_results\\.csv$", ".cox_results_full.csv", out_res)
  fwrite(res_empty, out_res_full)
  fwrite(res_empty, out_res)

  sink(out_sum)
  cat("=== Expr+Mut Cox Summary (03e) ===\n")
  cat("Cancer:", cancer, "\n")
  cat("Mutation group:", mut_group, "\n")
  cat("Genes tested (core expr ∩ mut):", length(genes_core), "\n")
  cat("Valid models:", 0, "\n")
  cat("Significant expr term (FDR<0.05):", 0, "\n")
  cat("Significant mut term  (FDR<0.05):", 0, "\n")
  sink()

  quit(status = 0)
}

# ============================================================
# Combine & annotate
# ============================================================
res <- rbindlist(res_list, fill = TRUE)

res[, cancer    := cancer]
res[, mut_group := mut_group]

# FDR per term across genes
res[, FDR_expr := p.adjust(p_expr, "BH")]
res[, FDR_mut  := p.adjust(p_mut,  "BH")]

# Save full results
out_res_full <- sub("\\.cox_results\\.csv$", ".cox_results_full.csv", out_res)
fwrite(res, out_res_full)

# Significant subset: either term FDR<0.05
res_sig <- res[FDR_expr < 0.05 | FDR_mut < 0.05]
fwrite(res_sig, out_res)

# Summary file
sink(out_sum)
cat("=== Expr+Mut Cox Summary (03e) ===\n")
cat("Cancer:", cancer, "\n")
cat("Mutation group:", mut_group, "\n")
cat("Genes tested (core expr ∩ mut):", length(genes_core), "\n")
cat("Valid models:", nrow(res), "\n")
cat("Significant expr term (FDR<0.05):", sum(res$FDR_expr < 0.05, na.rm = TRUE), "\n")
cat("Significant mut term  (FDR<0.05):", sum(res$FDR_mut  < 0.05, na.rm = TRUE), "\n")
sink()

say("%s [DONE] Results: %s (full) | %s (sig)",
    run_id, out_res_full, out_res)

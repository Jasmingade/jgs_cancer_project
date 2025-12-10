#!/usr/bin/env Rscript

# -------------------------------------------------------------------
# Example:
'
Rscript 01_transcriptomics/pipeline/scripts/03_split_mutations.R \
    01_transcriptomics/data/mutation/mc3.v0.2.8.PUBLIC.maf.gz \
    01_transcriptomics/data/raw/TCGA_clinical.csv \
    ensembl
'
# -------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# Arguments
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: 03_split_mutations.R <mutation_file> <clinical_file> [gene_id_type]")
}

mutation_file <- args[1]
clinical_file <- args[2]
gene_id_type  <- ifelse(length(args) >= 3, args[3], "ensembl")

gene_id_type <- match.arg(gene_id_type, c("ensembl", "hugo"))

out_dir <- "01_transcriptomics/out/03_mutation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n========== CLEAN GENE-LEVEL MUTATION SPLIT ==========\n")
cat("[INFO] Mutation file: ", mutation_file, "\n")
cat("[INFO] Clinical file: ", clinical_file, "\n")
cat("[INFO] Output dir:    ", out_dir, "\n")
cat("[INFO] Gene ID type:  ", gene_id_type, "\n\n")

# =============================================================================
# 1) Load clinical → define canonical case list per cancer
# =============================================================================
clin <- fread(clinical_file)

clin[, case_id := substr(bcr_patient_barcode, 1, 12)]
clin[, cancer  := paste0("TCGA_", cancer_type)]

clin_map <- clin[, .(case_id, cancer)]
valid_cases <- sort(unique(clin_map$case_id))

setkey(clin_map, case_id)

cat("[INFO] Clinical patients:", length(valid_cases), "\n")

# =============================================================================
# 2) Define mutation groups
# =============================================================================
class_to_group <- list(
  missense_or_inframe = c(
  "Missense_Mutation",
  "In_Frame_Del",
  "In_Frame_Ins"
  ),
  
  truncating_LOF = c(
    "Frame_Shift_Del", 
    "Frame_Shift_Ins",
    "Nonsense_Mutation", 
    "Nonstop_Mutation"
  ),
  splice = c("Splice_Site"),
  #silent = c("Silent"),
  rna = c("RNA")
)

interesting_classes <- unlist(class_to_group, use.names = FALSE)
group_levels <- names(class_to_group)

# =============================================================================
# 3) Load mutation data
# =============================================================================
mut_cols <- c(
  "Hugo_Symbol", "Gene", "Variant_Classification",
  "Tumor_Sample_Barcode"
)

mut <- fread(mutation_file, select = mut_cols)
mut[, case_id := substr(Tumor_Sample_Barcode, 1, 12)]

# Keep only clinical cases
mut <- mut[case_id %in% valid_cases]

# Keep only 10 interesting classes
mut <- mut[Variant_Classification %in% interesting_classes]

# Attach cancer type
mut <- merge(mut, clin_map, by = "case_id", all.x = TRUE)
mut <- mut[!is.na(cancer)]

# Assign mutation group
mut[, mut_group := NA_character_]
for (grp in group_levels)
  mut[Variant_Classification %in% class_to_group[[grp]], mut_group := grp]

mut <- mut[!is.na(mut_group)]
cat("[INFO] Mutation rows retained:", nrow(mut), "\n")

# Decide whether to use Ensembl or Hugo
gene_col <- if (gene_id_type == "ensembl") "Gene" else "Hugo_Symbol"

# =============================================================================
# Helper: Create full mutation matrix (gene × patient)
# =============================================================================
create_full_matrix <- function(dt, feature_col, case_ids, out_file) {

  if (nrow(dt) == 0) {
    warning("[WARN] No data for ", out_file)
    return(NULL)
  }

  # True mutation pairs
  dtu <- unique(dt[, .(feature_id = get(feature_col), case_id)])
  dtu[, key := paste(feature_id, case_id)]

  # Full (gene × case) grid
  full_grid <- CJ(
    feature_id = unique(dtu$feature_id),
    case_id    = case_ids,
    sorted = TRUE
  )
  full_grid[, key := paste(feature_id, case_id)]

  # Assign 1 if the pair appears in true mutation events
  full_grid[, has_mut := as.integer(key %in% dtu$key)]
  full_grid[, key := NULL]

  # Cast wide
  wide <- dcast(full_grid, feature_id ~ case_id, value.var = "has_mut", fill = 0L)

  fwrite(wide, out_file)
  cat("[INFO] Wrote:", out_file, " (",
      nrow(wide), " genes × ", length(case_ids), " samples)\n")

  return(wide)
}

# =============================================================================
# 4) Split by cancer
# =============================================================================
cancers <- sort(unique(mut$cancer))

for (canc in cancers) {

  cat("\n============== ", canc, " ==============\n")

  mut_c <- mut[cancer == canc]
  cancer_cases <- sort(unique(clin_map[cancer == canc, case_id]))

  gene_dir <- file.path(out_dir, canc, "gene")
  dir.create(gene_dir, recursive = TRUE)

  gene_group_mats <- list()

  # ------------------------------------------
  # 4.1 Generate mutation matrices for each group
  # ------------------------------------------
  for (grp in group_levels) {

    dt_grp <- mut_c[
      mut_group == grp &
      !is.na(get(gene_col)) & get(gene_col) != "",
      .(feature_id = get(gene_col), case_id)
    ]

    out_file <- file.path(
      gene_dir,
      sprintf("%s_gene_%s_%s.csv", canc, gene_id_type, grp)
    )

    gene_group_mats[[grp]] <- create_full_matrix(dt_grp, "feature_id", cancer_cases, out_file)
  }

  # ------------------------------------------
  # 4.2 Build coding_any (LOF OR missense)
  # ------------------------------------------
  if (!is.null(gene_group_mats$truncating_or_splice_LOF) &&
      !is.null(gene_group_mats$missense_or_inframe)) {

    cat("[INFO] Creating coding_any...\n")

    lof <- melt(gene_group_mats$truncating_or_splice_LOF,
                id.vars="feature_id", variable.name="case_id", value.name="lof")
    mis <- melt(gene_group_mats$missense_or_inframe,
                id.vars="feature_id", variable.name="case_id", value.name="mis")

    any_long <- merge(lof, mis, by=c("feature_id", "case_id"), all=TRUE)
    any_long[is.na(lof), lof := 0]
    any_long[is.na(mis), mis := 0]

    any_long[, any := as.integer(lof == 1 | mis == 1)]

    any_wide <- dcast(any_long, feature_id ~ case_id, value.var="any", fill = 0L)

    out_file_any <- file.path(
      gene_dir,
      sprintf("%s_gene_%s_coding_any.csv", canc, gene_id_type)
    )
    fwrite(any_wide, out_file_any)

    cat("[INFO] Wrote coding_any:", out_file_any, "\n")
  }
}

cat("\n[INFO] Mutation splitting complete.\n")

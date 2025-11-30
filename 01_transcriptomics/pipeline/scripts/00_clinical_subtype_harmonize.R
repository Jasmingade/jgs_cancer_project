harmonize_clinical_subtype <- function(subtype_tsv, clinical_csv, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  library(data.table)

  first12 <- function(x) substr(x, 1, 12)

  # Subtype
  SUB <- fread(subtype_tsv, sep = "\t")
  sc <- intersect(c("sampleID","sample","Sample"), names(SUB))[1]
  if (is.na(sc)) stop("No sample column found in subtype TSV.")
  SUB[, case_id := first12(get(sc))]
  subtype_col <- intersect(c("Subtype_Selected","Subtype","mRNA.subtype","Subtype_mRNA"), names(SUB))[1]
  if (!is.na(subtype_col)) {
    SUB <- SUB[!is.na(get(subtype_col))]
    SUB <- SUB[, .(Subtype = names(sort(table(get(subtype_col)), decreasing=TRUE))[1]), by = case_id]
  } else {
    SUB <- unique(SUB[, .(case_id, Subtype = NA_character_)])
  }

  # Clinical
  CLIN <- fread(clinical_csv)
  if (!"bcr_patient_barcode" %in% names(CLIN)) stop("Expected bcr_patient_barcode in clinical.")
  CLIN[, case_id := first12(bcr_patient_barcode)]
  if ("OS.time" %in% names(CLIN)) CLIN[, OS_time := as.numeric(OS.time)]
  if ("OS" %in% names(CLIN)) CLIN[, OS_event := as.integer(OS)]
  CLIN[, age := as.numeric(age_at_initial_pathologic_diagnosis)]
  CLIN[, sex := tolower(as.character(gender))]

  # --- Robust stage harmonization ---
  simp <- function(x) {
    x <- tolower(trimws(x))
    # remove placeholder or invalid entries
    bad <- c("not available", "not applicable", "unknown", 
             "discrepancy", "stage x", "i/ii nos", "na", "is")
    x[grepl(paste(bad, collapse="|"), x)] <- NA

    # extract I, II, III, IV (with or without "stage" prefix)
    stage <- sub(".*?(stage\\s*)?(i{1,3}|iv)\\b.*", "\\2", x, perl=TRUE)
    stage[!grepl("^(i{1,3}|iv)$", stage)] <- NA
    stage <- toupper(stage)
    stage
  }

  if ("ajcc_pathologic_tumor_stage" %in% names(CLIN)) {
    CLIN[, stage := simp(ajcc_pathologic_tumor_stage)]
  } else if ("clinical_stage" %in% names(CLIN)) {
    CLIN[, stage := simp(clinical_stage)]
  }

  # Merge subtype
  M <- merge(CLIN, SUB, by = "case_id", all.x = TRUE)

  # Fall back cancer_type from Subtype prefix if missing
  if (!"cancer_type" %in% names(M) || all(is.na(M$cancer_type))) {
    M[, cancer_type := tstrsplit(Subtype, ".", fixed=TRUE, keep=1)]
  }

  # Keep essentials
  keep <- intersect(c("case_id","cancer_type","OS_time","OS_event","age","sex","stage","Subtype"), names(M))
  M <- M[!is.na(OS_time) & !is.na(OS_event), ..keep]

  # Write per cancer
  for (ct in sort(na.omit(unique(M$cancer_type)))) {
    fwrite(M[cancer_type == ct], file.path(out_dir, sprintf("TCGA_%s_clinical.csv", ct)))
  }
}

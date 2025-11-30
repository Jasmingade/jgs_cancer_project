#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexUpset)
  library(viridis)
  library(patchwork)
})

say <- function(...) message(sprintf(...))

# ============================================================
# PATHS
# ============================================================
expr_root <- "01_transcriptomics/out/03a_univariate_coxph"
mut_root  <- "01_transcriptomics/out/03b_mutation_univariate_coxph"
outdir    <- "01_transcriptomics/out/05_plots/model3_overlap"
tx2gene_path <- "01_transcriptomics/data/raw/tx2gene.csv"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
say("[INFO] Saving outputs into: %s", outdir)

# ============================================================
# Helper: strip ENSEMBL version
# ============================================================
strip_ver <- function(x) sub("\\.\\d+$", "", x)

# ============================================================
# Load tx2gene (ENST → ENSG) for iso_* mapping
# ============================================================
say("[INFO] Loading tx2gene mapping from: %s", tx2gene_path)
tx2gene <- fread(tx2gene_path)
if (!all(c("tx_id", "gene_id") %in% names(tx2gene))) {
  stop("tx2gene file must contain columns: tx_id, gene_id")
}

tx2gene[, tx_id   := strip_ver(tx_id)]
tx2gene[, gene_id := strip_ver(gene_id)]
setkey(tx2gene, tx_id)

say("[INFO] tx2gene loaded: %d rows, %d unique tx_id, %d unique gene_id",
    nrow(tx2gene), uniqueN(tx2gene$tx_id), uniqueN(tx2gene$gene_id))

# ============================================================
# Loader: Expression CoxPH (gene / iso_log / iso_frac)
#   - gene      : feature = ENSG gene_id
#   - iso_*     : feature (ENST) → tx2gene → ENSG gene_id
#   - mutation  : handled in separate loader
# ============================================================
load_expr_file <- function(f){
    say("[DEBUG]   [EXPR] Reading %s", f)
    dt <- fread(f)

    needed <- c("feature", "HR", "FDR")
    if (!all(needed %in% names(dt))) {
        say("[WARN]   [EXPR] Skipping %s — missing one of: %s",
            basename(f), paste(needed, collapse = ", "))
        return(NULL)
    }

    fname <- basename(f)

    # Show a few examples of strip_ver
    ex_raw <- head(dt$feature[!is.na(dt$feature)], 5L)
    if (length(ex_raw) > 0) {
      ex_stripped <- strip_ver(ex_raw)
      say("[DEBUG]   [EXPR] Example feature stripping (first up to 5 rows):")
      for (k in seq_along(ex_raw)) {
        say("[DEBUG]     %s -> %s", ex_raw[k], ex_stripped[k])
      }
    }

    # Strip version from feature ID
    dt[, feature := strip_ver(feature)]
    dt[, HR  := as.numeric(HR)]
    dt[, FDR := as.numeric(FDR)]

    # metadata
    dt[, cancer :=
          sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", fname)]

    dt[, data_type :=
          fifelse(grepl("_gene", fname), "gene",
          fifelse(grepl("_iso_log", fname), "iso_log",
          fifelse(grepl("_iso_frac", fname), "iso_frac", NA_character_)))]

    # ---------- map iso_* (ENST) → gene_id (ENSG) using tx2gene ----------
    if (grepl("_iso_log", fname) || grepl("_iso_frac", fname)) {
        say("[DEBUG]   [EXPR] %s detected as isoform data → mapping tx_id to gene_id", fname)

        dt[, tx_id := feature]  # current feature is ENST (stripped)
        # join to tx2gene
        dt <- merge(
          dt,
          tx2gene[, .(tx_id, gene_id)],
          by = "tx_id",
          all.x = TRUE
        )

        n_before  <- nrow(dt)
        n_missing <- sum(is.na(dt$gene_id))
        if (n_missing > 0) {
          say("[WARN]   [EXPR] %s: %d / %d iso rows had no tx2gene mapping and were dropped",
              fname, n_missing, n_before)
          dt <- dt[!is.na(gene_id)]
        }

        # for overlap, we want gene-level IDs as 'feature'
        dt[, feature := gene_id]

        # (optional) drop tx_id if you don't need it further
        # dt[, tx_id := NULL]

        say("[DEBUG]   [EXPR] %s: %d rows remain after tx2gene mapping",
            fname, nrow(dt))
    }

    # keep only significant rows
    out <- dt[FDR < 0.05 & is.finite(HR) & HR > 0]
    say("[DEBUG]   [EXPR] %s → %d significant rows (FDR < 0.05)",
        fname, nrow(out))
    out
}

# ============================================================
# Loader: Mutation-only CoxPH (already at gene level)
# ============================================================
load_mut_file <- function(f){
    say("[DEBUG]   [MUT] Reading %s", f)
    dt <- fread(f)

    needed <- c("feature", "HR", "FDR")
    if (!all(needed %in% names(dt))) {
        say("[WARN]   [MUT] Skipping %s — missing one of: %s",
            basename(f), paste(needed, collapse = ", "))
        return(NULL)
    }

    # Show a few examples of strip_ver
    ex_raw <- head(dt$feature[!is.na(dt$feature)], 5L)
    if (length(ex_raw) > 0) {
      ex_stripped <- strip_ver(ex_raw)
      say("[DEBUG]   [MUT] Example feature stripping (first up to 5 rows):")
      for (k in seq_along(ex_raw)) {
        say("[DEBUG]     %s -> %s", ex_raw[k], ex_stripped[k])
      }
    }

    dt[, feature := strip_ver(feature)]
    dt[, HR  := as.numeric(HR)]
    dt[, FDR := as.numeric(FDR)]

    fname <- basename(f)

    dt[, cancer := sub("^TCGA_([A-Z0-9]+)_.*$", "\\1", fname)]

    mut_group <- sub("^TCGA_[A-Z0-9]+_mutation_(.*)\\.cox_results\\.csv$", "\\1", fname)
    dt[, data_type := paste0("mut_", mut_group)]

    out <- dt[FDR < 0.05 & is.finite(HR) & HR > 0]
    say("[DEBUG]   [MUT] %s → %d significant rows (FDR < 0.05)",
        fname, nrow(out))
    out
}

# ============================================================
# LOAD ALL DATA
# ============================================================
say("[INFO] Loading expression Cox results from: %s", expr_root)
expr_files <- list.files(expr_root, pattern="cox_results\\.csv$", recursive=TRUE, full.names=TRUE)
say("[INFO] Found %d expression result files", length(expr_files))

if (length(expr_files) == 0) {
    say("[WARN] No expression Cox result files found — expr_root is empty?")
    expr_res <- data.table()
} else {
    expr_res <- rbindlist(
        Filter(Negate(is.null), lapply(expr_files, load_expr_file)),
        fill = TRUE
    )
}
say("[INFO] Loaded %d significant expression rows (combined)", nrow(expr_res))

say("[INFO] Loading mutation Cox results from: %s", mut_root)
mut_files <- list.files(mut_root, pattern="cox_results\\.csv$", recursive=TRUE, full.names=TRUE)
say("[INFO] Found %d mutation result files", length(mut_files))

if (length(mut_files) == 0) {
    say("[WARN] No mutation Cox result files found — mut_root is empty?")
    mut_res <- data.table()
} else {
    mut_res <- rbindlist(
        Filter(Negate(is.null), lapply(mut_files, load_mut_file)),
        fill = TRUE
    )
}
say("[INFO] Loaded %d significant mutation rows (combined)", nrow(mut_res))

# Combine
res_all <- rbindlist(list(expr_res, mut_res), fill=TRUE)
say("[INFO] Total combined rows (expr + mut): %d", nrow(res_all))

all_cancers <- sort(unique(res_all$cancer))
say("[INFO] Found %d cancers with ≥1 significant feature", length(all_cancers))
say("[INFO] Cancers: %s", paste(all_cancers, collapse = ", "))

# Global collector for heatmap
global_list <- list()

# Collector for per-cancer UpSet plots
upset_plots <- list()

# ============================================================
# LOOP PER CANCER — UpSet + overlap matrices
# ============================================================
say("[INFO] Starting per-cancer overlap + UpSet loop")
for (C in all_cancers) {

    say("[INFO] ---- Cancer: %s ----", C)
    sub <- res_all[cancer == C]

    say("[INFO]   %s: %d total significant rows", C, nrow(sub))

    if (nrow(sub) == 0) {
        say("[WARN]   %s: No significant features — skipping", C)
        next
    }

    # Create per-datatype lists
    sets <- split(sub$feature, sub$data_type)
    say("[INFO]   %s: data types present: %s",
        C, paste(names(sets), collapse = ", "))

    # Save for global use
    global_list[[C]] <- sets

    # Build logical membership matrix
    say("[INFO]   %s: Building logical membership matrix", C)
    all_features <- sort(unique(unlist(sets)))
    M <- data.table(feature = all_features)

    for (nm in names(sets)){
        M[, (nm) := feature %in% sets[[nm]]]
    }

    # save matrix
    mat_path <- file.path(outdir, sprintf("%s_overlap_matrix.csv", C))
    fwrite(M, mat_path)
    say("[INFO]   %s: Saved overlap membership matrix → %s", C, mat_path)

    # ============ Per-cancer UpSet plot ===============
    upset_cols <- setdiff(colnames(M), "feature")

    # Skip if fewer than 2 data types
    if (length(upset_cols) < 2) {
        say(sprintf("[SKIP]   %s: only one data_type → no UpSet plot", C))

        # Save a small text note for reproducibility
        writeLines(
          paste("Cancer", C, "has only one data type:", upset_cols),
          file.path(outdir, sprintf("%s_upset_SKIPPED.txt", C))
        )
    } else {

        say("[INFO]   %s: Building UpSet plot", C)
        p_up <- upset(
          M,
          upset_cols,
          base_annotations = list('Intersection size' = intersection_size()),
          wrap = TRUE
        ) +
          ggtitle(C)  # title per panel so you know which cancer is which

        # Instead of saving here, store in list
        upset_plots[[C]] <- p_up
    }

    # ============ Per-cancer feature barplot ==========
    sum_dt <- data.table(
        data_type = names(sets),
        n = sapply(sets, length)
    )

    p_bar <- ggplot(sum_dt, aes(x=data_type, y=n+1, fill=data_type)) +
      geom_col() +
      scale_y_log10() +
      theme_bw(12) +
      theme(axis.text.x=element_text(angle=45,hjust=1)) +
      labs(title=paste("Significant Feature Counts —", C),
           y="Count (log10 scale)", x="Data Type")

    ggsave(
      file.path(outdir, sprintf("%s_feature_counts.png", C)),
      p_bar, width=10, height=6, dpi=300
    )
}

say("[INFO] Finished per-cancer loop. Combining UpSet plots into one figure...")

# ============================================================
# COMBINE ALL PER-CANCER UPSETS INTO ONE GRID
# ============================================================
if (length(upset_plots) > 0) {
  # choose layout
  n_col <- 8
  n_row <- ceiling(length(upset_plots) / n_col)

  combined_up <- wrap_plots(upset_plots, ncol = n_col)

  up_grid_path <- file.path(outdir, "ALL_cancers_upset_grid.png")
  ggsave(
    up_grid_path,
    combined_up,
    width  = 5 * n_col,        # tweak as needed
    height = 4 * n_row,        # tweak as needed
    dpi = 300
  )
  say("[DONE] Saved combined UpSet grid → %s", up_grid_path)
} else {
  say("[WARN] No per-cancer UpSet plots were created → nothing to combine.")
}

# ============================================================
# GLOBAL barplot (grouped by cancer)
# ============================================================
say("[INFO] Building global grouped barplot of significant feature counts")

global_counts <- rbindlist(lapply(names(global_list), function(C){
   data.table(
     cancer = C,
     data_type = names(global_list[[C]]),
     n = sapply(global_list[[C]], length)
   )
}), fill = TRUE)

say("[INFO] global_counts: %d rows, %d cancers × datatypes",
    nrow(global_counts), length(unique(global_counts$cancer)))

p_bar_global <- ggplot(global_counts,
       aes(x=cancer, y=n+1, fill=data_type)) +
   geom_col(position="dodge") +
   scale_y_log10() +
   theme_bw(12) +
   theme(axis.text.x=element_text(angle=60, hjust=1)) +
   labs(title="Significant Feature Counts Across All Cancers",
        y="Count (log10 scale)")

global_bar_path <- file.path(outdir,"GLOBAL_feature_counts.png")
ggsave(global_bar_path,
       p_bar_global, width=16, height=9, dpi=300)
say("[DONE] Saved global feature count barplot → %s", global_bar_path)

# ============================================================
# INTERSECTION HEATMAP (cancer × datatype)
# ============================================================
say("[INFO] Creating intersection heatmap (cancer × data_type)")

n_gc <- nrow(global_counts)
say("[INFO] global_counts has %d rows; computing %d pairwise overlaps",
    n_gc, n_gc * n_gc)

pairs <- CJ(i=1:n_gc, j=1:n_gc)
pairs[, overlap := mapply(function(a,b){
     length(intersect(
       global_list[[global_counts$cancer[a]]][[global_counts$data_type[a]]],
       global_list[[global_counts$cancer[b]]][[global_counts$data_type[b]]]
     ))
}, i, j)]

say("[INFO] Pairwise overlap matrix computed")

heat_df <- data.table(
    group_x = paste(global_counts$cancer[pairs$i],
                    global_counts$data_type[pairs$i], sep="__"),
    group_y = paste(global_counts$cancer[pairs$j],
                    global_counts$data_type[pairs$j], sep="__"),
    overlap = pairs$overlap
)

say("[INFO] heat_df rows: %d", nrow(heat_df))

p_heat <- ggplot(heat_df, aes(group_x, group_y, fill=overlap)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal(base_size=7) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5)) +
  labs(title="Intersection Heatmap of Significant Features",
       x="", y="")

heat_path <- file.path(outdir, "GLOBAL_overlap_heatmap.png")
ggsave(heat_path,
       p_heat, width=14, height=14, dpi=300)
say("[DONE] Saved global intersection heatmap → %s", heat_path)

say("[DONE] All overlap plots generated successfully → %s", outdir)

# ============================================================
# GLOBAL UpSet across all cancers (cancer × data_type)
# with size-based filtering for performance/readability
# ============================================================
say("[INFO] Building global UpSet plot across all cancer×data_type sets")

# ---- knobs for how big the global UpSet is allowed to be ----
MAX_GLOBAL_UPSET_SETS <- 30L   # max number of cancer×modality sets in global UpSet
MIN_FEATURES_PER_SET  <- 20L   # min features required for a set to be included

# merge all sets into one long list: key = "CANCER__datatype"
global_union <- list()
for (C in names(global_list)) {
    for (nm in names(global_list[[C]])) {
        key <- paste(C, nm, sep = "__")
        global_union[[key]] <- unique(global_list[[C]][[nm]])
    }
}

say("[INFO] global_union contains %d cancer×data_type sets", length(global_union))

if (length(global_union) < 2) {
    say("[SKIP] GLOBAL UpSet: fewer than two sets with significant features.")
} else {
    # --- compute set sizes ---
    set_sizes <- vapply(global_union, length, integer(1))
    set_sizes_dt <- data.table(
        set   = names(set_sizes),
        size  = as.integer(set_sizes)
    )[order(-size)]

    say("[INFO] Set size summary (first 10):")
    print(head(set_sizes_dt, 10))

    # --- filter by MIN_FEATURES_PER_SET ---
    keep_sets <- set_sizes_dt[size >= MIN_FEATURES_PER_SET, set]
    say("[INFO] Sets with ≥ %d features: %d / %d",
        MIN_FEATURES_PER_SET, length(keep_sets), length(global_union))

    if (length(keep_sets) == 0) {
        say("[SKIP] GLOBAL UpSet: no sets have ≥ %d features", MIN_FEATURES_PER_SET)
    } else {
        # --- limit to top MAX_GLOBAL_UPSET_SETS by size ---
        if (length(keep_sets) > MAX_GLOBAL_UPSET_SETS) {
            keep_sets <- keep_sets[seq_len(MAX_GLOBAL_UPSET_SETS)]
            say("[INFO] Limiting to top %d sets by size", MAX_GLOBAL_UPSET_SETS)
        }

        say("[INFO] Keeping %d sets in global UpSet:", length(keep_sets))
        print(set_sizes_dt[set %in% keep_sets])

        kept_lists <- global_union[keep_sets]

        # --- build feature universe for kept sets ---
        all_feats <- unique(unlist(kept_lists))
        say("[INFO] GLOBAL UpSet: %d unique features in kept sets (before pruning)",
            length(all_feats))

        # --- drop features that appear in only one set (no intersections) ---
        membership_mat <- sapply(kept_lists, function(v) all_feats %in% v)
        row_membership_count <- rowSums(membership_mat)

        feats_keep_idx <- row_membership_count >= 2
        n_dropped <- sum(!feats_keep_idx)
        all_feats2 <- all_feats[feats_keep_idx]

        say("[INFO] Dropping %d features present in only one set", n_dropped)
        say("[INFO] %d features remain for intersection matrix", length(all_feats2))

        if (length(all_feats2) == 0) {
            say("[SKIP] GLOBAL UpSet: no features belong to ≥2 of the kept sets.")
        } else {
            # rebuild membership matrix on pruned feature set
            Mglob <- data.table(feature = all_feats2)
            for (nm in names(kept_lists)) {
                Mglob[, (nm) := feature %in% kept_lists[[nm]]]
            }

            cols <- setdiff(names(Mglob), "feature")
            say("[INFO] GLOBAL UpSet: final sets used: %s",
                paste(cols, collapse = ", "))

            # build UpSet
            p_up_glob <- upset(
              Mglob,
              cols,
              base_annotations = list('Intersection size' = intersection_size()),
              wrap = TRUE
            )

            global_up_path <- file.path(outdir,"GLOBAL_upset_cancer_by_datatype.png")
            ggsave(global_up_path,
                   p_up_glob, width = 14, height = 10, dpi = 300)
            say("[DONE] Saved global UpSet (cancer×data_type) → %s", global_up_path)
        }
    }
}

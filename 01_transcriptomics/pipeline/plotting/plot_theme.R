#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
})

# ============================================================
# Clean, thesis-ready theme
# ============================================================
theme_big <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.text.x = element_text(size = base_size * 0.8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = base_size * 0.8),
      plot.title = element_text(size = base_size * 1.1, face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = base_size * 0.9, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# ============================================================
# Color palettes (consistent across entire thesis)
# ============================================================
palette_expr_types <- c(
  gene    = "#0072B2",  # blue
  iso_log = "#009E73",  # green
  iso_frac = "#E69F00"  # orange
)

palette_mut_types <- c(
  truncating_or_splice_LOF = "#56B4E9",   # red
  missense_or_inframe      = "#D55E00",   # blue
  rna_other                = "#F0E442",   # green
  coding_any               = "#CC79A7"    # purple
)

# ============================================================
# HR reference line
# ============================================================
hr_reference_line <- function() {
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", alpha = 0.6)
}

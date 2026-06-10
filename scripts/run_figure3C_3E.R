# run_figure3C_3E_GO_KEGG_enrichment.R
#
# Generate GO (Figure 3C) and KEGG (Figure 3E) enrichment chord diagrams for NUMB/NUMBL isoforms.
#
# Inputs:
# - results/processed/Figure_2/differential_fused_NUMB_isoforms_annotated.csv
#
# Outputs:
# - results/figures/Figure_3/Figure_3C_GO_chord_top50.tiff
# - results/figures/Figure_3/Figure_3E_KEGG_chord_top50.tiff
# - results/processed/Figure_3/NUMB_isoforms_GO_enrichment.csv
# - results/processed/Figure_3/GO_raw/GO_*_raw.csv
# - results/tables/NUMB_isoforms_KEGG_enrichment.csv

# ============================================================
# 0) Input / Output paths
# ============================================================

diff_file <- file.path("results", "processed", "Figure_2",
                       "differential_fused_NUMB_isoforms_annotated.csv")

fig3_dir   <- file.path("results", "figures",   "Figure_3")
proc_dir   <- file.path("results", "processed", "Figure_3")
tables_dir <- file.path("results", "tables")

dir.create(fig3_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

out_fig_3C <- file.path(fig3_dir,   "Figure_3C_GO_chord_top50.tiff")
out_fig_3E <- file.path(fig3_dir,   "Figure_3E_KEGG_chord_top50.tiff")

out_table_GO   <- file.path(proc_dir,   "NUMB_isoforms_GO_enrichment.csv")
out_table_KEGG <- file.path(tables_dir, "NUMB_isoforms_KEGG_enrichment.csv")

stopifnot(file.exists(diff_file))

isoforms <- c("p65", "p66", "p71", "p72", "numbl")

# ============================================================
# Libraries
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(circlize)
  library(RColorBrewer)
})

# ------------------------------------------------------------
# Load differential expression dataset
# ------------------------------------------------------------

diff_df <- read.csv(diff_file, stringsAsFactors = FALSE)
diff_df <- diff_df[!is.na(diff_df$hgnc_symbol) & diff_df$hgnc_symbol != "", ]

# ------------------------------------------------------------
# Build filtered gene subsets per isoform and direction
# Threshold: |log2FC| > 0.8 and p_value < 0.001
# ------------------------------------------------------------

filtered_sets <- list()

for (iso in isoforms) {

  log_col  <- paste0("log2FC_", iso)
  pval_col <- paste0("p_value_", iso)

  sub_df <- diff_df[, c("hgnc_symbol", log_col, pval_col)]
  colnames(sub_df) <- c("hgnc_symbol", "log2FC", "p_value")

  up   <- sub_df[sub_df$log2FC >  0.8 & sub_df$p_value < 0.001, ]
  down <- sub_df[sub_df$log2FC < -0.8 & sub_df$p_value < 0.001, ]

  filtered_sets[[paste0(iso, "_up")]]   <- up
  filtered_sets[[paste0(iso, "_down")]] <- down
}

# ------------------------------------------------------------
# Helper: SYMBOL -> ENTREZ ID conversion
# ------------------------------------------------------------

genes_to_entrez <- function(symbols) {
  unique(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )$ENTREZID
  )
}

# ------------------------------------------------------------
# Fixed isoform order for chord diagrams (shared by GO and KEGG)
# ------------------------------------------------------------

isoform_order <- c(
  "p72_up", "p71_up", "p66_up", "p65_up", "numbl_up",
  "numbl_down", "p65_down", "p66_down", "p71_down", "p72_down"
)

# ------------------------------------------------------------
# Colorblind-safe palette (Okabe-Ito) for isoform sectors
# Up-regulated isoforms: orange / teal
# Down-regulated isoforms: light blue / dark blue
# ------------------------------------------------------------

isoform_colors <- c(
  "p72_up"     = "#E69F00",  # orange
  "p71_up"     = "#E69F00",  # orange
  "p66_up"     = "#009E73",  # teal
  "p65_up"     = "#009E73",  # teal
  "numbl_up"   = "#009E73",  # teal
  "numbl_down" = "#56B4E9",  # sky blue
  "p65_down"   = "#56B4E9",  # sky blue
  "p66_down"   = "#56B4E9",  # sky blue
  "p71_down"   = "#0072B2",  # dark blue
  "p72_down"   = "#0072B2"   # dark blue
)

# ============================================================
# Figure 3C — GO enrichment chord diagram
# ============================================================

# ------------------------------------------------------------
# STEP 1: Run GO enrichment and save raw CSVs
# (skipped if files already exist — delete to regenerate)
# ------------------------------------------------------------

go_raw_dir <- file.path(proc_dir, "GO_raw")
dir.create(go_raw_dir, recursive = TRUE, showWarnings = FALSE)

for (name in names(filtered_sets)) {

  archivo_raw <- file.path(go_raw_dir, paste0("GO_", name, "_raw.csv"))
  if (file.exists(archivo_raw)) next

  genes  <- filtered_sets[[name]]$hgnc_symbol
  entrez <- genes_to_entrez(genes)

  if (length(entrez) > 0) {

    go <- enrichGO(
      entrez,
      OrgDb    = org.Hs.eg.db,
      ont      = "BP",
      readable = TRUE
    )

    if (!is.null(go) && nrow(go@result) > 0) {
      df           <- go@result
      df$Isoform   <- name
      df$Direction <- ifelse(grepl("_up$", name), "up", "down")
      df           <- df[order(df$pvalue), ]
      write.csv(df, archivo_raw, row.names = FALSE)
    }
  }
}

# ------------------------------------------------------------
# STEP 2: Load full raw CSVs (unfiltered — auditable)
# ------------------------------------------------------------

go_results <- list()

for (name in names(filtered_sets)) {

  archivo_raw <- file.path(go_raw_dir, paste0("GO_", name, "_raw.csv"))

  if (file.exists(archivo_raw)) {
    df <- read.csv(archivo_raw, stringsAsFactors = FALSE)
    if (nrow(df) > 0) {
      go_results[[name]] <- df
    }
  }
}

# ------------------------------------------------------------
# STEP 2.5: Limit to top 100 most significant terms per condition
#
# enrichGO() already applies pvalueCutoff = 0.05 internally, so all
# terms in the raw CSVs are nominally significant. When a condition
# returns more than 100 terms, only the 100 with the lowest p-value
# are kept. This matches the M&M description: "limiting the 100 most
# significant terms to the total when the total number exceeded this
# threshold."
# ------------------------------------------------------------

max_terms_per_group <- 100

go_final_raw <- do.call(rbind, go_results)

go_final_top100 <- do.call(rbind, lapply(go_results, function(df) {
  df_sorted <- df[order(df$pvalue), ]
  head(df_sorted, max_terms_per_group)
}))

go_final <- go_final_top100[
  , c("ID", "Description", "Isoform", "Direction", "pvalue", "p.adjust", "FoldEnrichment")
]

cat("GO terms loaded (all from enrichGO):         ", nrow(go_final_raw),    "\n")
cat("GO terms after top-100 per condition:        ", nrow(go_final),        "\n")
cat("Unique terms after top-100 filter:           ", length(unique(go_final$ID)), "\n")

# ------------------------------------------------------------
# Diagnostic: terms per isoform/direction combination
# ------------------------------------------------------------

cat("\nGO terms per isoform/direction (after top-100 filter):\n")
print(table(go_final$Isoform))

go_counts_diag <- table(go_final$ID)
go_shared_diag <- names(go_counts_diag[go_counts_diag >= 2])
cat("\nShared GO terms (>= 2 isoform/direction combinations) per isoform:\n")
print(table(go_final[go_final$ID %in% go_shared_diag, ]$Isoform))

# ------------------------------------------------------------
# STEP 3: Build the chord diagram
#
# Selection criteria:
#   1. Top 100 most significant terms per isoform/direction (Step 2.5)
#   2. Term must appear in >= 2 isoform/direction combinations
#   3. Top 50 terms ranked by co-enrichment frequency
#   4. Terms classified as:
#        "up"    — enriched only in up-regulated isoform sets
#        "down"  — enriched only in down-regulated isoform sets
#        "mixed" — enriched in both directions (bidirectional)
# ------------------------------------------------------------

write.csv(go_final, out_table_GO, row.names = FALSE)

go_counts    <- table(go_final$ID)
go_shared    <- names(go_counts[go_counts >= 2])
go_shared_df <- go_final[go_final$ID %in% go_shared, ]

top_n      <- 50
freq       <- sort(table(go_shared_df$ID), decreasing = TRUE)
go_top_ids <- names(freq)[1:top_n]
go_top     <- go_shared_df[go_shared_df$ID %in% go_top_ids, ]

classify_go <- function(df) {
  terms <- unique(df$ID)
  group <- character(length(terms))
  for (i in seq_along(terms)) {
    sub   <- df[df$ID == terms[i], ]
    ups   <- sum(sub$Direction == "up")
    downs <- sum(sub$Direction == "down")
    if      (ups > 0 & downs == 0) group[i] <- "up"
    else if (downs > 0 & ups == 0) group[i] <- "down"
    else                           group[i] <- "mixed"
  }
  data.frame(ID = terms, group, stringsAsFactors = FALSE)
}

groups_go <- classify_go(go_top)

# Order: up terms first, then bidirectional (mixed), then down
order_go_ids <- c(
  groups_go$ID[groups_go$group == "up"],
  groups_go$ID[groups_go$group == "mixed"],
  groups_go$ID[groups_go$group == "down"]
)

go_top$ID <- factor(go_top$ID, levels = order_go_ids)
go_top    <- go_top[order(go_top$ID), ]

relations_go <- go_top[, c("Isoform", "ID")]

# Colorblind-safe palette for GO term sectors
# Using an expanded Okabe-Ito-based gradient
go_colors <- setNames(
  colorRampPalette(c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999"
  ))(length(order_go_ids)),
  order_go_ids
)

all_colors_go <- c(isoform_colors, go_colors)

# ------------------------------------------------------------
# Render and save Figure 3C
# ------------------------------------------------------------

tiff(
  out_fig_3C,
  width = 3000, height = 3000, res = 600, compression = "lzw"
)

circos.clear()
circos.par(gap.degree = 2, start.degree = 90)

chordDiagram(
  relations_go,
  order             = c(isoform_order, order_go_ids),
  grid.col          = all_colors_go,
  transparency      = 0.4,
  annotationTrack   = "grid",
  preAllocateTracks = list(track.height = 0.06)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sn <- get.cell.meta.data("sector.index")
    xl <- get.cell.meta.data("xlim")
    circos.text(mean(xl), 0, sn, facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.5)
  },
  bg.border = NA
)

title(paste0("GO terms shared across NUMB isoform groups (", top_n, " most frequent)"))
dev.off()

# ============================================================
# Figure 3E — KEGG enrichment chord diagram
# ============================================================

# ------------------------------------------------------------
# KEGG-specific gene subsets
# Threshold: |log2FC| > 0.8 and p_value < 0.01 (as per M&M)
# ------------------------------------------------------------

filtered_sets_kegg <- list()

for (iso in isoforms) {
  
  log_col  <- paste0("log2FC_", iso)
  pval_col <- paste0("p_value_", iso)
  
  sub_df <- diff_df[, c("hgnc_symbol", log_col, pval_col)]
  colnames(sub_df) <- c("hgnc_symbol", "log2FC", "p_value")
  
  up   <- sub_df[sub_df$log2FC >  0.8 & sub_df$p_value < 0.01, ]
  down <- sub_df[sub_df$log2FC < -0.8 & sub_df$p_value < 0.01, ]
  
  filtered_sets_kegg[[paste0(iso, "_up")]]   <- up
  filtered_sets_kegg[[paste0(iso, "_down")]] <- down
}

# ------------------------------------------------------------
# STEP 1: Run KEGG enrichment
# pvalueCutoff = 0.05 (default), top 100 per condition
# ------------------------------------------------------------

kegg_results <- list()

for (name in names(filtered_sets_kegg)) {
  
  genes  <- filtered_sets_kegg[[name]]$hgnc_symbol
  entrez <- genes_to_entrez(genes)
  
  if (length(entrez) > 0) {
    
    kegg <- enrichKEGG(
      gene         = entrez,
      organism     = "hsa",
      pvalueCutoff = 0.05
    )
    
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      df           <- kegg@result
      df$Isoform   <- name
      df$Direction <- ifelse(grepl("_up$", name), "up", "down")
      df           <- df[order(df$pvalue), ]
      if (nrow(df) > 100) df <- df[1:100, ]
      df <- df[, c("ID", "Description", "Isoform", "Direction", "pvalue", "FoldEnrichment")]
      kegg_results[[name]] <- df
    }
  }
}

# ------------------------------------------------------------
# STEP 2: Merge, filter shared pathways, select top 50
# ------------------------------------------------------------

kegg_final <- do.call(rbind, kegg_results)

write.csv(kegg_final, out_table_KEGG, row.names = FALSE)

kegg_counts <- table(kegg_final$ID)
shared_ids  <- names(kegg_counts[kegg_counts >= 2])
kegg_shared <- kegg_final[kegg_final$ID %in% shared_ids, ]

top_n        <- 50
freq         <- sort(table(kegg_shared$ID), decreasing = TRUE)
kegg_top_ids <- names(freq)[1:top_n]
kegg_top     <- kegg_shared[kegg_shared$ID %in% kegg_top_ids, ]

classify_kegg <- function(df) {
  terms <- unique(df$ID)
  group <- character(length(terms))
  for (i in seq_along(terms)) {
    sub   <- df[df$ID == terms[i], ]
    ups   <- sum(sub$Direction == "up")
    downs <- sum(sub$Direction == "down")
    if      (ups > 0 & downs == 0) group[i] <- "up"
    else if (downs > 0 & ups == 0) group[i] <- "down"
    else                           group[i] <- "mixed"
  }
  data.frame(ID = terms, group, stringsAsFactors = FALSE)
}

groups_kegg <- classify_kegg(kegg_top)

order_kegg_ids <- c(
  groups_kegg$ID[groups_kegg$group == "up"],
  groups_kegg$ID[groups_kegg$group == "mixed"],
  groups_kegg$ID[groups_kegg$group == "down"]
)

kegg_top$ID <- factor(kegg_top$ID, levels = order_kegg_ids)
kegg_top    <- kegg_top[order(kegg_top$ID), ]

relations <- kegg_top[, c("Isoform", "ID")]

kegg_colors <- setNames(
  colorRampPalette(c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999"
  ))(length(order_kegg_ids)),
  order_kegg_ids
)

all_colors_kegg <- c(isoform_colors, kegg_colors)

# ------------------------------------------------------------
# Render and save Figure 3E
# ------------------------------------------------------------

tiff(
  out_fig_3E,
  width = 3000, height = 3000, res = 600, compression = "lzw"
)

circos.clear()
circos.par(gap.degree = 2, start.degree = 90)

chordDiagram(
  relations,
  order             = c(isoform_order, order_kegg_ids),
  grid.col          = all_colors_kegg,
  transparency      = 0.4,
  annotationTrack   = "grid",
  preAllocateTracks = list(track.height = 0.06)
)

circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sn <- get.cell.meta.data("sector.index")
    xl <- get.cell.meta.data("xlim")
    circos.text(mean(xl), 0, sn, facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.5)
  },
  bg.border = NA
)

title(paste0("KEGG pathways shared across NUMB isoform groups (", top_n, " most frequent)"))
dev.off()
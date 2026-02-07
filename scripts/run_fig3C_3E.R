# run_figure3C_3E_GO_KEGG_enrichment.R
#
# Generate GO (Figure 3C) and KEGG (Figure 3E) enrichment chord diagrams for NUMB/NUMBL isoforms.
#
# Inputs:
# - results/processed/differential_fused_NUMB_isoforms_annotated.csv
#
# Outputs:
# - results/figures/Figure_3/Figure_3C_GO_chord_top40.tiff
# - results/figures/Figure_3/Figure_3E_KEGG_chord_top50.tiff
# - results/processed/Figure_3/NUMB_isoforms_GO_enrichment.csv
# - results/tables/NUMB_isoforms_KEGG_enrichment.csv

# ============================================================
# 0) Input / Output
# ============================================================

diff_file <- file.path("results", "processed", "Figure_2", "differential_fused_NUMB_isoforms_annotated.csv")

fig3_dir   <- file.path("results", "figures", "Figure_3")
proc_dir   <- file.path("results", "processed", "Figure_3")
tables_dir <- file.path("results", "tables")

dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

out_fig_3C <- file.path(fig3_dir, "Figure_3C_GO_chord_top40.tiff")
out_fig_3E <- file.path(fig3_dir, "Figure_3E_KEGG_chord_top50.tiff")

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
# Load differential dataset (once)
# ------------------------------------------------------------

diff_df <- read.csv(diff_file, stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Build filtered gene subsets (shared by GO + KEGG)
# ------------------------------------------------------------

filtered_sets <- list()

for (iso in isoforms) {
  
  log_col  <- paste0("log2FC_", iso)
  pval_col <- paste0("p_value_", iso)
  
  sub_df <- diff_df[, c("hgnc_symbol", log_col, pval_col)]
  colnames(sub_df) <- c("hgnc_symbol", "log2FC", "p_value")
  
  up   <- sub_df[sub_df$log2FC >  0.8 & sub_df$p_value < 0.01, ]
  down <- sub_df[sub_df$log2FC < -0.8 & sub_df$p_value < 0.01, ]
  
  filtered_sets[[paste0(iso, "_up")]]   <- up
  filtered_sets[[paste0(iso, "_down")]] <- down
}

# ------------------------------------------------------------
# SYMBOL → ENTREZ
# ------------------------------------------------------------

genes_to_entrez <- function(symbols) {
  unique(
    bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )$ENTREZID
  )
}

# ============================================================
# Figure 3C — GO chord diagram
# ============================================================

go_results <- list()

for (name in names(filtered_sets)) {
  
  genes <- filtered_sets[[name]]$hgnc_symbol
  entrez <- genes_to_entrez(genes)
  
  if (length(entrez) > 0) {
    
    go <- enrichGO(
      entrez,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      readable = TRUE
    )
    
    if (!is.null(go) && nrow(go@result) > 0) {
      
      df <- go@result
      df$Isoform <- name
      df$Direction <- ifelse(grepl("_up$", name), "up", "down")
      df <- df[order(df$pvalue), ]
      if (nrow(df) > 100) df <- df[1:100, ]
      df <- df[, c("ID", "Description", "Isoform", "Direction", "pvalue", "FoldEnrichment")]
      
      go_results[[name]] <- df
    }
  }
}

go_final <- do.call(rbind, go_results)

write.csv(go_final, out_table_GO, row.names = FALSE)

# Shared ≥2
go_counts <- table(go_final$ID)
go_shared <- names(go_counts[go_counts >= 2])
go_shared_df <- go_final[go_final$ID %in% go_shared, ]

top_n <- 40
freq <- sort(table(go_shared_df$ID), decreasing = TRUE)
go_top_ids <- names(freq)[1:top_n]
go_top <- go_shared_df[go_shared_df$ID %in% go_top_ids, ]

classify_go <- function(df) {
  
  terms <- unique(df$ID)
  group <- character(length(terms))
  
  for (i in seq_along(terms)) {
    
    sub <- df[df$ID == terms[i], ]
    ups <- sum(sub$Direction == "up")
    downs <- sum(sub$Direction == "down")
    
    if (ups > 0 & downs == 0) group[i] <- "up"
    else if (downs > 0 & ups == 0) group[i] <- "down"
    else group[i] <- "mixed"
  }
  
  data.frame(ID = terms, group)
}

groups_go <- classify_go(go_top)

order_go_ids <- c(
  groups_go$ID[groups_go$group == "up"],
  groups_go$ID[groups_go$group == "mixed"],
  groups_go$ID[groups_go$group == "down"]
)

go_top$ID <- factor(go_top$ID, levels = order_go_ids)
go_top <- go_top[order(go_top$ID), ]

relations_go <- go_top[, c("Isoform", "ID")]

isoform_colors_go <- c(
  "p72_up"="#FFA500","p71_up"="#FFA500",
  "p66_up"="#FF0000","numbl_up"="#FF0000","p65_up"="#FF0000",
  "p66_down"="#00BFFF","p65_down"="#00BFFF","numbl_down"="#00BFFF",
  "p71_down"="#0000FF","p72_down"="#0000FF"
)

go_colors <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(length(order_go_ids)),
  order_go_ids
)

all_colors_go <- c(isoform_colors_go, go_colors)

tiff(
  out_fig_3C,
  width = 3000,
  height = 3000,
  res = 600,
  compression = "lzw"
)

circos.clear()
circos.par(gap.degree = 2, start.degree = 90)

chordDiagram(
  relations_go,
  grid.col = all_colors_go,
  transparency = 0.4,
  annotationTrack = "grid",
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
# Figure 3E — KEGG chord diagram
# ============================================================

kegg_results <- list()

for (name in names(filtered_sets)) {
  
  genes <- filtered_sets[[name]]$hgnc_symbol
  entrez <- genes_to_entrez(genes)
  
  if (length(entrez) > 0) {
    
    kegg <- enrichKEGG(gene = entrez, organism = "hsa", pvalueCutoff = 0.05)
    
    if (!is.null(kegg) && nrow(kegg@result) > 0) {
      
      df <- kegg@result
      df$Isoform <- name
      df$Direction <- ifelse(grepl("_up$", name), "up", "down")
      df <- df[order(df$pvalue), ]
      if (nrow(df) > 100) df <- df[1:100, ]
      df <- df[, c("ID", "Description", "Isoform", "Direction", "pvalue", "FoldEnrichment")]
      
      kegg_results[[name]] <- df
    }
  }
}

kegg_final <- do.call(rbind, kegg_results)

write.csv(kegg_final, out_table_KEGG, row.names = FALSE)

kegg_counts <- table(kegg_final$ID)
shared_ids <- names(kegg_counts[kegg_counts >= 2])
kegg_shared <- kegg_final[kegg_final$ID %in% shared_ids, ]

top_n <- 50
freq <- sort(table(kegg_shared$ID), decreasing = TRUE)
kegg_top_ids <- names(freq)[1:top_n]
kegg_top <- kegg_shared[kegg_shared$ID %in% kegg_top_ids, ]

classify_kegg <- function(df) {
  
  terms <- unique(df$ID)
  group <- character(length(terms))
  
  for (i in seq_along(terms)) {
    
    sub <- df[df$ID == terms[i], ]
    ups <- sum(sub$Direction == "up")
    downs <- sum(sub$Direction == "down")
    
    if (ups > 0 & downs == 0) group[i] <- "up"
    else if (downs > 0 & ups == 0) group[i] <- "down"
    else group[i] <- "mixed"
  }
  
  data.frame(ID = terms, group)
}

groups_kegg <- classify_kegg(kegg_top)

order_kegg_ids <- c(
  groups_kegg$ID[groups_kegg$group == "up"],
  groups_kegg$ID[groups_kegg$group == "mixed"],
  groups_kegg$ID[groups_kegg$group == "down"]
)

kegg_top$ID <- factor(kegg_top$ID, levels = order_kegg_ids)
kegg_top <- kegg_top[order(kegg_top$ID), ]

relations <- kegg_top[, c("Isoform", "ID")]

kegg_colors <- setNames(
  colorRampPalette(brewer.pal(9, "Set1"))(length(order_kegg_ids)),
  order_kegg_ids
)

all_colors <- c(isoform_colors_go, kegg_colors)

tiff(
  out_fig_3E,
  width = 3000,
  height = 3000,
  res = 600,
  compression = "lzw"
)

circos.clear()
circos.par(gap.degree = 2, start.degree = 90)

chordDiagram(
  relations,
  grid.col = all_colors,
  transparency = 0.4,
  annotationTrack = "grid",
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

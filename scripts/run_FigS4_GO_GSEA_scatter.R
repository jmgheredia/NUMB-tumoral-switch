# run_FigureS4_GO_GSEA_scatter.R
#
# Supplementary Figure S4 — GSEA (GO Biological Process) NES scatter plots
# comparing NUMB/NUMBL isoforms.
#
# This script:
# 1) Runs GSEA (GO BP) independently for each isoform using ranked log2FC
#    (most extreme transcript per gene, preserving sign).
# 2) Combines results across isoforms.
# 3) Keeps GO terms shared by ≥2 isoforms.
# 4) Generates pairwise NES scatter plots for all isoform combinations.
#
# Inputs:
# - results/processed/differential_fused_NUMB_isoforms_annotated.csv
#
# Outputs:
# - results/figures/Supplementary/S4/NES_scatter_<ISO1>_vs_<ISO2>.tiff
# - results/tables/Supplementary_Table_3_GO_GSEA_shared_terms.csv
# - results/tables/GSEA_GO_summary_all_isoforms.csv
# - results/tables/GSEA_GO_<ISOFORM>_full.csv

# ============================================================
# 0) Input / Output
# ============================================================

diff_file <- file.path(
  "results", "processed", "Figure_2",   "differential_fused_NUMB_isoforms_annotated.csv")

figS4_dir  <- file.path("results", "figures", "Supplementary", "S4")
tables_dir <- file.path("results", "tables")

dir.create(figS4_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(diff_file))

# ============================================================
# Libraries
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(tidyr)
})

# ============================================================
# Load differential dataset
# ============================================================

diff_df <- read.csv(diff_file, stringsAsFactors = FALSE)

isoforms <- c("p72", "p71", "p66", "p65", "numbl")

fc_cols <- paste0("log2FC_", isoforms)
diff_df[fc_cols] <- lapply(diff_df[fc_cols], function(x)
  as.numeric(gsub(",", ".", x))
)

# ============================================================
# Run GSEA (GO BP) per isoform
# ============================================================

gsea_results <- list()

for (iso in isoforms) {
  
  fc_col <- paste0("log2FC_", iso)
  
  df_iso <- diff_df[, c("hgnc_symbol", fc_col)]
  colnames(df_iso)[2] <- "log2FC"
  
  df_iso <- df_iso[!is.na(df_iso$log2FC) & df_iso$hgnc_symbol != "", ]
  
  df_iso_max <- df_iso %>%
    group_by(hgnc_symbol) %>%
    summarise(extreme_log2FC = log2FC[which.max(abs(log2FC))]) %>%
    ungroup()
  
  gene_list <- df_iso_max$extreme_log2FC
  names(gene_list) <- df_iso_max$hgnc_symbol
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea <- gseGO(
    geneList = gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.1,
    verbose = FALSE
  )
  
  gsea_results[[iso]] <- gsea
  
  if (!is.null(gsea) && nrow(gsea@result) > 0) {
    write.csv(
      gsea@result,
      file.path(tables_dir, paste0("GSEA_GO_", iso, "_full.csv")),
      row.names = FALSE
    )
  }
}

# ============================================================
# Combine GSEA results across isoforms
# ============================================================

gsea_summary <- do.call(rbind, lapply(names(gsea_results), function(iso) {
  
  g <- gsea_results[[iso]]
  
  if (!is.null(g) && nrow(g@result) > 0) {
    
    out <- g@result[, c("ID", "Description", "NES", "pvalue", "p.adjust", "qvalue")]
    out$Isoform <- iso
    out
    
  } else {
    NULL
  }
}))

write.csv(
  gsea_summary,
  file.path(tables_dir, "GSEA_GO_summary_all_isoforms.csv"),
  row.names = FALSE
)

gsea_shared <- gsea_summary %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  ungroup()

write.csv(
  gsea_shared,
  file.path(tables_dir, "Supplementary_Table_3_GO_GSEA_shared_terms.csv"),
  row.names = FALSE
)

# ============================================================
# Pairwise NES scatter plots (Supplementary Figure S4)
# ============================================================

gsea_wide <- gsea_shared %>%
  dplyr::select(ID, NES, Isoform) %>%
  pivot_wider(names_from = Isoform, values_from = NES)

pairs <- combn(isoforms, 2, simplify = FALSE)

for (p in pairs) {
  
  iso1 <- p[1]
  iso2 <- p[2]
  
  plot_df <- gsea_wide
  
  g <- ggplot(plot_df, aes_string(x = iso1, y = iso2)) +
    geom_point(na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    theme_minimal() +
    labs(
      title = paste("NES comparison:", iso1, "vs", iso2),
      x = paste0("NES (", iso1, ")"),
      y = paste0("NES (", iso2, ")")
    )
  
  ggsave(
    filename = file.path(figS4_dir,
                         paste0("NES_scatter_", iso1, "_vs_", iso2, ".tiff")),
    plot = g,
    width = 6,
    height = 5,
    units = "in",
    dpi = 300,
    device = "tiff"
  )
}

message("Supplementary Figure S4 (GO GSEA NES scatter) completed successfully.")

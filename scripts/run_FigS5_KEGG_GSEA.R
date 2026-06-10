# run_FigS4_KEGG_GSEA.R
#
# Supplementary Figure S5 — KEGG GSEA analyses across NUMB/NUMBL isoforms.
#
# This script:
# 1) Runs GSEA (KEGG) independently for each isoform using ranked log2FC
#    (most extreme transcript per gene, preserving sign).
# 2) Saves full KEGG GSEA tables per isoform.
# 3) Generates selected KEGG GSEA enrichment plots.
# 4) Computes average NES for two isoform groups (p72/p71 vs p66/p65/NUMBL)
#    and produces the group comparison plot.
#
# Inputs:
# - results/processed/differential_fused_NUMB_isoforms_annotated.csv
#
# Outputs:
# Figures (Supplementary Figure S5):
# - results/figures/Supplementary/S5/<KEGGID>_<ISOFORM>.tiff
# - results/figures/Supplementary/S5/NES_group_comparison.tiff
#
# Tables:
# - results/tables/GSEA_KEGG_<ISOFORM>_full.csv
# - results/tables/GSEA_KEGG_summary_all_isoforms.csv

# ============================================================
# 0) Input / Output
# ============================================================

diff_file <- file.path( "results", "processed", "Figure_2", "differential_fused_NUMB_isoforms_annotated.csv")

figS5_dir  <- file.path("results", "figures", "Supplementary", "S5")
tables_dir <- file.path("results", "tables")

dir.create(figS5_dir, recursive = TRUE, showWarnings = FALSE)
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
# Run KEGG GSEA per isoform
# ============================================================

gsea_results <- list()

for (iso in isoforms) {
  
  fc_col <- paste0("log2FC_", iso)
  
  df_iso <- diff_df[, c("hgnc_symbol", fc_col)]
  colnames(df_iso)[2] <- "log2FC"
  df_iso <- df_iso[!is.na(df_iso$log2FC) & df_iso$hgnc_symbol != "", ]
  
  entrez <- bitr(
    df_iso$hgnc_symbol,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  df_iso <- merge(df_iso, entrez, by.x = "hgnc_symbol", by.y = "SYMBOL")
  
  df_iso_max <- df_iso %>%
    group_by(ENTREZID) %>%
    summarise(extreme_log2FC = log2FC[which.max(abs(log2FC))]) %>%
    ungroup()
  
  gene_list <- df_iso_max$extreme_log2FC
  names(gene_list) <- df_iso_max$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea <- gseKEGG(
    geneList = gene_list,
    organism = "hsa",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.1,
    verbose = FALSE
  )
  
  gsea_results[[iso]] <- gsea
  
  if (!is.null(gsea) && nrow(gsea@result) > 0) {
    
    gsea@result$Description <- gsub(
      " - Homo sapiens \\(human\\)",
      "",
      gsea@result$Description
    )
    
    write.csv(
      gsea@result,
      file.path(tables_dir, paste0("GSEA_KEGG_", iso, "_full.csv")),
      row.names = FALSE
    )
  }
}

# ============================================================
# Combine KEGG GSEA results across isoforms
# ============================================================

gsea_summary <- do.call(rbind, lapply(names(gsea_results), function(iso) {
  
  g <- gsea_results[[iso]]
  
  if (!is.null(g) && nrow(g@result) > 0) {
    
    out <- g@result[, c(
      "ID", "Description",
      "setSize", "enrichmentScore", "NES",
      "pvalue", "p.adjust", "qvalue"
    )]
    out$Isoform <- iso
    out
    
  } else {
    NULL
  }
}))

write.csv(
  gsea_summary,
  file.path(tables_dir, "GSEA_KEGG_summary_all_isoforms.csv"),
  row.names = FALSE
)

# ============================================================
# Selected KEGG GSEA plots
# ============================================================

kegg_terms <- list(
  "hsa04081" = "Hormone signaling",
  "hsa04060" = "Cytokine-cytokine receptor interaction",
  "hsa04510" = "Focal adhesion",
  "hsa04152" = "AMPK signaling pathway",
  "hsa00020" = "Citrate cycle (TCA cycle)",
  "hsa00620" = "Pyruvate metabolism"
)

for (kegg_id in names(kegg_terms)) {
  for (iso in isoforms) {
    
    g <- gsea_results[[iso]]
    
    if (!is.null(g) && kegg_id %in% g@result$ID) {
      
      p <- gseaplot2(
        g,
        geneSetID = kegg_id,
        title = paste(kegg_id, "-", kegg_terms[[kegg_id]], "-", iso)
      )
      
      ggsave(
        filename = file.path(figS5_dir, paste0(kegg_id, "_", iso, ".tiff")),
        plot = p,
        width = 8,
        height = 6,
        units = "in",
        dpi = 300,
        device = "tiff"
      )
    }
  }
}

# ============================================================
# NES group comparison (p72/p71 vs p66/p65/NUMBL)
# ============================================================

group_isoforms <- list(
  Group1 = c("p72", "p71"),
  Group2 = c("p66", "p65", "numbl")
)

gsea_group <- gsea_summary %>%
  filter(Isoform %in% unlist(group_isoforms)) %>%
  group_by(ID) %>%
  filter(any(Isoform %in% group_isoforms$Group1) &
           any(Isoform %in% group_isoforms$Group2)) %>%
  ungroup()

gsea_means <- gsea_group %>%
  mutate(Group = case_when(
    Isoform %in% group_isoforms$Group1 ~ "Group1",
    Isoform %in% group_isoforms$Group2 ~ "Group2"
  )) %>%
  group_by(ID, Description, Group) %>%
  summarise(mean_NES = mean(NES), .groups = "drop")

gsea_means_wide <- gsea_means %>%
  pivot_wider(names_from = Group, values_from = mean_NES) %>%
  drop_na()

df_points <- gsea_means_wide %>%
  pivot_longer(cols = c(Group1, Group2),
               names_to = "Group",
               values_to = "NES")

df_points$Group <- recode(
  df_points$Group,
  "Group1" = "p72/p71",
  "Group2" = "p66/p65/NUMBL"
)

df_lines <- gsea_means_wide %>%
  mutate(Term = ID) %>%
  dplyr::select(Term, NES_Group1 = Group1, NES_Group2 = Group2) %>%
  pivot_longer(cols = starts_with("NES_"),
               names_to = "Group",
               values_to = "NES") %>%
  mutate(Group = recode(Group,
                        "NES_Group1" = "p72/p71",
                        "NES_Group2" = "p66/p65/NUMBL"))

p_group <- ggplot() +
  geom_point(data = df_points, aes(x = Group, y = NES), size = 2) +
  geom_line(data = df_lines,
            aes(x = Group, y = NES, group = Term),
            color = "grey40", alpha = 0.6) +
  scale_x_discrete(limits = c("p72/p71", "p66/p65/NUMBL")) +
  labs(
    title = "Average NES comparison between isoform groups",
    x = "Isoform group",
    y = "Mean NES"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(figS5_dir, "NES_group_comparison.tiff"),
  plot = p_group,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  device = "tiff"
)

message("Supplementary Figure S5 (KEGG GSEA) completed successfully.")

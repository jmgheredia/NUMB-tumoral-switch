# ============================================================
# run_Fig5C_CCLE_NUMB_GO_enrichment.R
#
# Performs GO Biological Process enrichment on genes differentially expressed
# between High vs Low NUMB score CCLE BRCA cell lines.
#
# Input:
#   - results/processed/Figure5/gene_expression_scoreNUMB_high.csv
#   - results/processed/Figure5/gene_expression_scoreNUMB_low.csv
#
# Output:
#   Figures:
#     - results/Figure_5/GO_enrichment_high_NUMB_score.tiff
#     - results/Figure_5/GO_enrichment_low_NUMB_score.tiff
# ============================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# ============================================================
# 0) Input / Output
# ============================================================

high_csv <- file.path("results", "processed", "Figure_5", "gene_expression_scoreNUMB_high.csv")
low_csv  <- file.path("results", "processed", "Figure_5", "gene_expression_scoreNUMB_low.csv")

fig_dir  <- file.path("results", "figures", "Figure_5")
proc_dir <- file.path("results", "processed", "Figure_5")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(high_csv))
  stop("Missing gene_expression_scoreNUMB_high.csv in results/processed/Figure5/")

if (!file.exists(low_csv))
  stop("Missing gene_expression_scoreNUMB_low.csv in results/processed/Figure5/")

# ============================================================
# 2) Load expression matrices
# ============================================================

exp_high <- read.csv(high_csv, header = TRUE, stringsAsFactors = FALSE)
exp_low  <- read.csv(low_csv, header = TRUE, stringsAsFactors = FALSE)

expr_high <- as.matrix(exp_high[, 4:ncol(exp_high)])
expr_low  <- as.matrix(exp_low[, 4:ncol(exp_low)])

# ============================================================
# 3) Compute log2 fold change (High - Low)
# ============================================================

mean_high <- rowMeans(expr_high, na.rm = TRUE)
mean_low  <- rowMeans(expr_low, na.rm = TRUE)

log2FC <- mean_high - mean_low

genes_df <- data.frame(
  gene_name = exp_high$gene_name,
  gene_id = exp_high$gene_id,
  transcript_id = exp_high$transcript_id,
  log2FC = log2FC
)

# ============================================================
# 4) Select up/downregulated genes
# ============================================================

genes_up <- unique(genes_df$gene_name[genes_df$log2FC > 1])
genes_down <- unique(genes_df$gene_name[genes_df$log2FC < -1])

genes_up_entrez <- bitr(
  genes_up,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

genes_down_entrez <- bitr(
  genes_down,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# ============================================================
# 5) GO enrichment (Biological Process)
# ============================================================

ego_up <- enrichGO(
  gene = genes_up_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

ego_down <- enrichGO(
  gene = genes_down_entrez$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# ============================================================
# 6) Dotplots
# ============================================================

p_up <- dotplot(ego_up, showCategory = 15, font.size = 10) +
  ggtitle("GO enrichment – High NUMB score") +
  theme(plot.margin = margin(5, 5, 5, 20, unit = "mm"))

p_down <- dotplot(ego_down, showCategory = 15, font.size = 10) +
  ggtitle("GO enrichment – Low NUMB score") +
  theme(plot.margin = margin(5, 5, 5, 20, unit = "mm"))

# ============================================================
# 7) Save figures
# ============================================================

ggsave(
  filename = file.path(fig_dir, "GO_enrichment_high_NUMB_score.tiff"),
  plot = p_up,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = file.path(fig_dir, "GO_enrichment_low_NUMB_score.tiff"),
  plot = p_down,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

message("Figure 5C GO enrichment completed.")

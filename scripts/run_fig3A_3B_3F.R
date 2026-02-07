# run_figure3A_3B_3F_panels.R
#
# Generate selected panels of Figure 3:
# - Figure 3A: UpSet plot of NUMB/NUMBL isoform intersections
# - Figure 3B: Heatmaps for signaling pathways (Notch, WNT, Hippo, Hedgehog)
# - Figure 3F: Metabolic pathway heatmap (Glycolysis / OXPHOS related)
#
# Inputs:
# - results/processed/differential_fused_NUMB_isoforms_annotated.csv
# - data/metadata/genes_by_pathway.csv
# - data/metadata/genes_by_mitochondrial_pathways.csv
#
# Outputs:
# - results/figures/Figure_3/Figure_3A_UpSet_isoforms_ordered.tiff
# - results/figures/Figure_3/Figure_3B_Heatmap_<PATHWAY>_log2FC_isoforms.tiff
# - results/figures/Figure_3/Figure_3F_Metabolic_heatmap_log2FC.tiff
# - results/tables/Supplementary_Table_2_pathways.csv
# - results/tables/Supplementary_Table_4_Figure3F_heatmap_log2FC.csv

# ============================================================
# 0) Input / Output
# ============================================================

diff_file      <- file.path("results", "processed", "Figure_2", "differential_fused_NUMB_isoforms_annotated.csv")
pathway_file   <- file.path("data", "metadata", "genes_by_pathway.csv")
metabolic_file <- file.path("data", "metadata", "genes_by_mitochondrial_pathways.csv")

output_fig_dir   <- file.path("results", "figures", "Figure_3")
output_table_dir <- file.path("results", "tables")

dir.create(output_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_table_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(diff_file))
stopifnot(file.exists(pathway_file))
stopifnot(file.exists(metabolic_file))

# ============================================================
# Libraries
# ============================================================

suppressPackageStartupMessages({
  library(UpSetR)
  library(dplyr)
  library(tibble)
  library(pheatmap)
})

# ====================================
# Load main differential dataset
# ====================================

diff_df <- read.csv(diff_file, stringsAsFactors = FALSE)

selected_cols <- c("log2FC_p72", "log2FC_p71", "log2FC_p66", "log2FC_p65", "log2FC_numbl")

diff_df[selected_cols] <- lapply(
  diff_df[selected_cols],
  function(x) as.numeric(gsub(",", ".", x))
)

# ===========================
# Figure 3A — UpSet plot
# ===========================

get_gene_list <- function(col, direction = "up", threshold = 1) {
  if (direction == "up") rownames(diff_df[diff_df[[col]] >= threshold, ])
  else rownames(diff_df[diff_df[[col]] <= -threshold, ])
}

gene_sets <- list(
  p65_up     = get_gene_list("log2FC_p65", "up"),
  p65_down   = get_gene_list("log2FC_p65", "down"),
  p66_up     = get_gene_list("log2FC_p66", "up"),
  p66_down   = get_gene_list("log2FC_p66", "down"),
  p71_up     = get_gene_list("log2FC_p71", "up"),
  p71_down   = get_gene_list("log2FC_p71", "down"),
  p72_up     = get_gene_list("log2FC_p72", "up"),
  p72_down   = get_gene_list("log2FC_p72", "down"),
  NUMBL_up   = get_gene_list("log2FC_numbl", "up"),
  NUMBL_down = get_gene_list("log2FC_numbl", "down")
)

all_genes <- unique(unlist(gene_sets))
upset_input <- data.frame(row.names = all_genes)

for (n in names(gene_sets)) {
  upset_input[[n]] <- ifelse(all_genes %in% gene_sets[[n]], 1, 0)
}

upset_input <- upset_input[rowSums(upset_input) > 1, ]

comb_counts <- table(apply(upset_input, 1, function(x) paste0(which(x == 1), collapse = "_")))
valid_combos <- names(comb_counts[comb_counts >= 20])

rows_to_keep <- apply(
  upset_input,
  1,
  function(x) paste0(which(x == 1), collapse = "_") %in% valid_combos
)

upset_final <- upset_input[rows_to_keep, ]

set_order <- c("p72_up","p71_up","p65_down","NUMBL_down","p66_down","p66_up","p72_down","p71_down","p65_up","NUMBL_up")

set_colors <- c(
  p72_up="#E41A1C", p71_up="#D62728", p65_down="#2171B5", NUMBL_down="#6BAED6",
  p66_down="#CCCCCC", p66_up="#969696", p72_down="#FCAE91", p71_down="#FDAE6B",
  p65_up="#9ECAE1", NUMBL_up="#C6DBEF"
)

out_fig_3A <- file.path(output_fig_dir, "Figure_3A_UpSet_isoforms_ordered.tiff")

tiff(
  out_fig_3A,
  width = 3000,
  height = 2000,
  res = 300,
  compression = "lzw"
)

upset(
  upset_final,
  sets = set_order,
  sets.bar.color = set_colors[set_order],
  keep.order = TRUE,
  order.by = "freq",
  main.bar.color = "#333333",
  mainbar.y.label = "Gene Intersections",
  sets.x.label = "Genes per Isoform (up/down)",
  text.scale = c(1.2,1.2,1,1,1,1)
)

dev.off()

message("Figure 3A saved to: ", out_fig_3A)

# ============================================================
# Figure 3B — Signaling pathway heatmaps + Supplementary Table 2
# ============================================================

pathways <- read.csv(pathway_file, stringsAsFactors = FALSE)
if ("gene_symbol" %in% colnames(pathways)) {
  colnames(pathways)[colnames(pathways) == "gene_symbol"] <- "hgnc_symbol"
}

merged <- merge(diff_df, pathways, by = "hgnc_symbol")

supp_table_3B <- merged %>%
  select(pathway, hgnc_symbol, transcript_id, all_of(selected_cols)) %>%
  arrange(pathway, desc(log2FC_p72))

out_table_3B <- file.path(output_table_dir, "Supplementary_Table_2_pathways.csv")
write.csv(supp_table_3B, out_table_3B, row.names = FALSE)

for (pw in unique(supp_table_3B$pathway)) {
  
  mat <- supp_table_3B %>%
    filter(pathway == pw) %>%
    column_to_rownames("transcript_id") %>%
    select(all_of(selected_cols)) %>%
    as.matrix()
  
  out_fig <- file.path(output_fig_dir, paste0("Figure_3B_Heatmap_", pw, "_log2FC_isoforms.tiff"))
  
  tiff(
    out_fig,
    width = 2000,
    height = 1500,
    res = 300,
    compression = "lzw"
  )
  
  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = seq(-3, 3, length.out = 101),
    show_rownames = FALSE,
    main = paste0(pw, " - log2FC"),
    border_color = NA
  )
  
  dev.off()
}

# ============================================================
# Figure 3F — Metabolic heatmap + Supplementary Table 4
# ============================================================

metabolic <- read.csv(metabolic_file, stringsAsFactors = FALSE)
if ("gene_symbol" %in% colnames(metabolic)) {
  colnames(metabolic)[colnames(metabolic) == "gene_symbol"] <- "hgnc_symbol"
}

matched <- merge(diff_df, metabolic, by = "hgnc_symbol")

selected_pathways <- c(
  "Glycolysis",
  "Pyruvate",
  "TCA_cycle",
  "Oxidative_phosphorylation",
  "Carbon_metabolism"
)

matched <- matched[matched$pathway %in% selected_pathways, ]
matched$pathway <- factor(matched$pathway, levels = selected_pathways)
matched <- matched[order(matched$pathway), ]

mat_met <- as.matrix(matched[, selected_cols])
rownames(mat_met) <- make.unique(rep("", nrow(mat_met)))

row_ann <- data.frame(Pathway = matched$pathway)
rownames(row_ann) <- rownames(mat_met)

pal <- c(
  Glycolysis = "#E41A1C",
  Pyruvate = "#377EB8",
  TCA_cycle = "#4DAF4A",
  Oxidative_phosphorylation = "#984EA3",
  Carbon_metabolism = "#FF7F00"
)

out_fig_3F <- file.path(output_fig_dir, "Figure_3F_Metabolic_heatmap_log2FC.tiff")

tiff(
  out_fig_3F,
  width = 2800,
  height = 3200,
  res = 400
)

pheatmap(
  mat_met,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1.5, 1.5, length.out = 100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  annotation_row = row_ann,
  annotation_colors = list(Pathway = pal[unique(row_ann$Pathway)]),
  fontsize_col = 10,
  main = "Heatmap of log2FC – Metabolic pathways",
  border_color = NA
)

dev.off()

supp_table_3F <- data.frame(
  pathway = matched$pathway,
  hgnc_symbol = matched$hgnc_symbol,
  transcript_id = matched$transcript_id,
  mat_met
)

out_table_3F <- file.path(output_table_dir, "Supplementary_Table_4_Figure3F_heatmap_log2FC.csv")
write.csv(supp_table_3F, out_table_3F, row.names = FALSE)

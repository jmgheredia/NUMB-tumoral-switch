# ============================================================
# run_FigureS6B_isoform_log2FC_PDMR.R
#
# Supplementary Figure S6B — Isoform-level log2FC distributions
# for selected signaling pathways (Hedgehog, Hippo, Notch, WNT)
# across PDMR passages.
#
# This script:
# 1) Loads isoform-level log2 expression data.
# 2) Computes log2 fold changes (Resection vs P0, P0 vs P4).
# 3) Generates jitter plots of isoform-level log2FC values
#    for pathway-specific gene sets.
#
# Inputs:
# - data/metadata/metadata_PDMR.csv
# - data/metadata/genes_by_pathway.csv
# - data/processed/pdmr_isoform_expression_log2.csv
#
# Outputs:
# - results/figures/Supplementary/S6/Fig_S6B_*.tiff
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ============================================================
# 0) Input / Output
# ============================================================

metadata_file <- file.path("data", "metadata", "metadata_PDMR.csv")
genes_by_pathway_file <- file.path("data", "metadata", "genes_by_pathway.csv")
expression_file <- file.path("data", "processed", "pdmr_isoform_expression_log2.csv")

out_dir <- file.path("results", "figures", "Supplementary", "S6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(metadata_file),
  file.exists(genes_by_pathway_file),
  file.exists(expression_file)
)

# ============================================================
# 1) Helper: normalize sample IDs (CRITICAL)
# ============================================================

normalize_id <- function(x) {
  x |>
    as.character() |>
    gsub("^X", "", x = _) |>
    gsub("_", ".", x = _) |>
    gsub("-", ".", x = _) |>
    gsub("\\.+", ".", x = _) |>
    trimws()
}

# ============================================================
# 2) Load data
# ============================================================

metadata <- read.csv(metadata_file, sep = ";", stringsAsFactors = FALSE)
genes_by_pathway <- read.csv(genes_by_pathway_file, stringsAsFactors = FALSE)
expression <- read.csv(expression_file, stringsAsFactors = FALSE, check.names = FALSE)

# ============================================================
# 3) Preprocess expression matrix
# ============================================================

expression$transcript_id <- sub("\\..*", "", expression$transcript_id)
expression <- expression[!duplicated(expression$transcript_id), ]
rownames(expression) <- expression$transcript_id

exp_ids <- expression[, c("transcript_id", "gene_id")]
exp_mat <- expression[, !(colnames(expression) %in% c("transcript_id", "gene_id"))]

colnames(exp_mat) <- normalize_id(colnames(exp_mat))

# ============================================================
# 4) Preprocess metadata
# ============================================================

metadata$sample <- normalize_id(metadata$PatientID.SpecimenID.SampleID)

get_samples <- function(passage) {
  metadata %>%
    filter(Passage_of_this_sample == passage) %>%
    pull(sample)
}

samples_R <- get_samples("Resection")
samples_0 <- get_samples("0")
samples_4 <- get_samples("4")

# ============================================================
# 5) Helper: plotting function
# ============================================================

plot_isoform_log2fc <- function(df, fc_col, gene_order,
                                title, ylim, outfile) {
  
  df_plot <- df %>%
    mutate(
      color_fc = ifelse(.data[[fc_col]] >= 0, "positive", "negative"),
      gene_id = factor(gene_id, levels = gene_order)
    )
  
  tiff(
    filename = file.path(out_dir, outfile),
    width = 2000,
    height = 1500,
    res = 300
  )
  
  print(
    ggplot(df_plot, aes(x = gene_id, y = .data[[fc_col]], color = color_fc)) +
      geom_jitter(width = 0.2, size = 2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      coord_cartesian(ylim = ylim) +
      labs(title = title, x = "Gene", y = "log2FC (isoforms)") +
      scale_color_manual(values = c("positive" = "red", "negative" = "blue")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "none",
        panel.grid = element_blank()
      )
  )
  
  dev.off()
}

# ============================================================
# ======================= HEDGEHOG ==========================
# ============================================================

genes_hedgehog <- genes_by_pathway %>% filter(pathway == "Hedgehog")

exp_hh <- expression %>% filter(gene_id %in% genes_hedgehog$gene_symbol)
exp_ids_hh <- exp_hh[, c("transcript_id", "gene_id")]
exp_mat_hh <- exp_hh[, !(colnames(exp_hh) %in% c("transcript_id", "gene_id"))]

log2fc_hh <- data.frame(
  transcript_id = rownames(exp_mat_hh),
  log2FC_R0 = rowMeans(exp_mat_hh[, samples_0], na.rm = TRUE) -
    rowMeans(exp_mat_hh[, samples_R], na.rm = TRUE),
  log2FC_04 = rowMeans(exp_mat_hh[, samples_4], na.rm = TRUE) -
    rowMeans(exp_mat_hh[, samples_0], na.rm = TRUE)
) %>%
  left_join(exp_ids_hh, by = "transcript_id")

gene_order_hh <- unique(genes_hedgehog$gene_symbol)

plot_isoform_log2fc(
  log2fc_hh, "log2FC_R0", gene_order_hh,
  "Isoform-level expression change (Hedgehog, Resection vs P0)",
  c(-7, 4),
  "Fig_S6B_hedgehog_isoforms_log2FC_R0.tiff"
)

plot_isoform_log2fc(
  log2fc_hh, "log2FC_04", gene_order_hh,
  "Isoform-level expression change (Hedgehog, P0 vs P4)",
  c(-7, 4),
  "Fig_S6B_hedgehog_isoforms_log2FC_04.tiff"
)

# ============================================================
# ========================= HIPPO ===========================
# ============================================================

genes_hippo_40 <- c(
  "MST1","MST2","LATS1","LATS2","SAV1","NF2","FRMD6","WWC1","AMOT","MOB1A",
  "YAP1","WWTR1","TEAD1","TEAD2","TEAD3","TEAD4",
  "BIRC5","BIRC2","AREG","MYC","CCND1","CCND2",
  "DLG1","LLGL1","LLGL2","SCRIB","RASSF1","RASSF6",
  "AJUBA","AMOTL1","AMOTL2","DVL1","DVL2","FBXW11",
  "ID1","ID2","SERPINE1","BBC3"
)

exp_hp <- expression %>% filter(gene_id %in% genes_hippo_40)
exp_ids_hp <- exp_hp[, c("transcript_id", "gene_id")]
exp_mat_hp <- exp_hp[, !(colnames(exp_hp) %in% c("transcript_id", "gene_id"))]

log2fc_hp <- data.frame(
  transcript_id = rownames(exp_mat_hp),
  log2FC_R0 = rowMeans(exp_mat_hp[, samples_0], na.rm = TRUE) -
    rowMeans(exp_mat_hp[, samples_R], na.rm = TRUE),
  log2FC_04 = rowMeans(exp_mat_hp[, samples_4], na.rm = TRUE) -
    rowMeans(exp_mat_hp[, samples_0], na.rm = TRUE)
) %>%
  left_join(exp_ids_hp, by = "transcript_id")

plot_isoform_log2fc(
  log2fc_hp, "log2FC_R0", genes_hippo_40,
  "Isoform-level expression change (Hippo, Resection vs P0)",
  c(-7.5, 2.5),
  "Fig_S6B_hippo_isoforms_log2FC_R0_40genes.tiff"
)

plot_isoform_log2fc(
  log2fc_hp, "log2FC_04", genes_hippo_40,
  "Isoform-level expression change (Hippo, P0 vs P4)",
  c(-7.5, 2.5),
  "Fig_S6B_hippo_isoforms_log2FC_04_40genes.tiff"
)

# ============================================================
# ========================== NOTCH ==========================
# ============================================================

genes_notch <- genes_by_pathway %>% filter(pathway == "Notch")

exp_notch <- expression %>% filter(gene_id %in% genes_notch$gene_symbol)
exp_ids_notch <- exp_notch[, c("transcript_id", "gene_id")]
exp_mat_notch <- exp_notch[, !(colnames(exp_notch) %in% c("transcript_id", "gene_id"))]

log2fc_notch <- data.frame(
  transcript_id = rownames(exp_mat_notch),
  log2FC_R0 = rowMeans(exp_mat_notch[, samples_0], na.rm = TRUE) -
    rowMeans(exp_mat_notch[, samples_R], na.rm = TRUE),
  log2FC_04 = rowMeans(exp_mat_notch[, samples_4], na.rm = TRUE) -
    rowMeans(exp_mat_notch[, samples_0], na.rm = TRUE)
) %>%
  left_join(exp_ids_notch, by = "transcript_id")

gene_order_notch <- unique(genes_notch$gene_symbol)

plot_isoform_log2fc(
  log2fc_notch, "log2FC_R0", gene_order_notch,
  "Isoform-level expression change (Notch, Resection vs P0)",
  c(-7.5, 5),
  "Fig_S6B_notch_isoforms_log2FC_R0.tiff"
)

plot_isoform_log2fc(
  log2fc_notch, "log2FC_04", gene_order_notch,
  "Isoform-level expression change (Notch, P0 vs P4)",
  c(-7.5, 5),
  "Fig_S6B_notch_isoforms_log2FC_04.tiff"
)

# ============================================================
# =========================== WNT ===========================
# ============================================================

genes_wnt <- genes_by_pathway %>% filter(pathway == "WNT")

exp_wnt <- expression %>% filter(gene_id %in% genes_wnt$gene_symbol)
exp_ids_wnt <- exp_wnt[, c("transcript_id", "gene_id")]
exp_mat_wnt <- exp_wnt[, !(colnames(exp_wnt) %in% c("transcript_id", "gene_id"))]

log2fc_wnt <- data.frame(
  transcript_id = rownames(exp_mat_wnt),
  log2FC_R0 = rowMeans(exp_mat_wnt[, samples_0], na.rm = TRUE) -
    rowMeans(exp_mat_wnt[, samples_R], na.rm = TRUE),
  log2FC_04 = rowMeans(exp_mat_wnt[, samples_4], na.rm = TRUE) -
    rowMeans(exp_mat_wnt[, samples_0], na.rm = TRUE)
) %>%
  left_join(exp_ids_wnt, by = "transcript_id")

gene_order_wnt <- unique(genes_wnt$gene_symbol)

plot_isoform_log2fc(
  log2fc_wnt, "log2FC_R0", gene_order_wnt,
  "Isoform-level expression change (WNT, Resection vs P0)",
  c(-9, 2.5),
  "Fig_S6B_wnt_isoforms_log2FC_R0.tiff"
)

plot_isoform_log2fc(
  log2fc_wnt, "log2FC_04", gene_order_wnt,
  "Isoform-level expression change (WNT, P0 vs P4)",
  c(-9, 2.5),
  "Fig_S6B_wnt_isoforms_log2FC_04.tiff"
)

message("Supplementary Figure S6B (isoform-level log2FC, PDMR) completed successfully.")

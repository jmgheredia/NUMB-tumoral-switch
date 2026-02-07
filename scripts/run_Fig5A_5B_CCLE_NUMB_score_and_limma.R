# ============================================================
# run_Fig5A_5B_CCLE_NUMB_score_and_limma.R
#
# Computes NUMB isoform-based score in CCLE BRCA cell lines (PDMR context)
# and performs LIMMA analysis.
# Generates Figures 5A and 5B.
#
# Input:
#   - data/processed/BCCL_gene_names_included.csv
#   - data/metadata/genes_by_pathway.csv
#
# Output:
#   Figures:
#     - results/figures/Figure_5/Figure5A_score_NUMB_histogram.tiff
#     - results/figures/Figure_5/Figure5B_score_NUMB_boxplot_filtered.tiff
#
#   Processed results:
#     - results/processed/Figure_5/genes_significant_by_scoreNUMB_pathways.csv
#     - results/processed/Figure_5/BCCL_grouped_by_NUMB_score.csv
#     - results/processed/Figure_5/gene_expression_scoreNUMB_high.csv
#     - results/processed/Figure_5/gene_expression_scoreNUMB_low.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(limma)
})

# ============================================================
# 0) Input / Output
# ============================================================

expr_csv <- file.path("data", "processed", "BCCL_gene_names_included.csv")
genes_by_pathway_csv <- file.path("data", "metadata", "genes_by_pathway.csv")

fig_dir  <- file.path("results", "figures", "Figure_5")
proc_dir <- file.path("results", "processed", "Figure_5")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(expr_csv))
  stop("Missing BCCL_gene_names_included.csv in data/processed/")

if (!file.exists(genes_by_pathway_csv))
  stop("Missing genes_by_pathway.csv in data/metadata/")

# ============================================================
# 2) Load data
# ============================================================

CCLE_breast <- read.csv(expr_csv, header = TRUE, stringsAsFactors = FALSE)
genes_by_pathway <- read.csv(genes_by_pathway_csv, header = TRUE, stringsAsFactors = FALSE)

# ============================================================
# 3) Filter genes expressed in at least one sample
# ============================================================

expression_cols <- colnames(CCLE_breast)[4:ncol(CCLE_breast)]
CCLE_breast <- CCLE_breast[rowSums(CCLE_breast[, expression_cols] != 0) > 0, ]

# ============================================================
# 4) Log2(TPM + 1)
# ============================================================

CCLE_log2 <- CCLE_breast
CCLE_log2[, expression_cols] <- log2(CCLE_log2[, expression_cols] + 1)

CCLE_log2$transcript_base <- sub("\\..*", "", CCLE_log2$transcript_id)

# ============================================================
# 5) Select functional NUMB / NUMBL isoforms
# ============================================================

selected_transcripts <- c(
  "ENST00000356296",
  "ENST00000252891",
  "ENST00000554546",
  "ENST00000557597",
  "ENST00000555238"
)

CCLE_isoforms <- CCLE_log2[CCLE_log2$transcript_base %in% selected_transcripts, ]

# ============================================================
# 6) Build expression matrix
# ============================================================

expr_cols <- colnames(CCLE_isoforms)[4:(ncol(CCLE_isoforms) - 1)]

iso_mat <- CCLE_isoforms[, c("gene_name", "transcript_base", expr_cols)]
rownames(iso_mat) <- paste0(iso_mat$gene_name, "_", iso_mat$transcript_base)

expr_mat <- t(iso_mat[, expr_cols])
expr_mat <- as.data.frame(apply(expr_mat, 2, as.numeric))
rownames(expr_mat) <- rownames(t(iso_mat[, expr_cols]))

# ============================================================
# 7) Compute NUMB score
# ============================================================

score_NUMB <- rowMeans(expr_mat[, c("NUMB_ENST00000555238", "NUMB_ENST00000557597")], na.rm = TRUE) -
  rowMeans(expr_mat[, c("NUMB_ENST00000554546", "NUMB_ENST00000356296", "NUMBL_ENST00000252891")], na.rm = TRUE)

# ============================================================
# 8) Remove 10% most central samples (5 cell lines)
# ============================================================

med <- median(score_NUMB, na.rm = TRUE)
dist <- abs(score_NUMB - med)
central <- names(sort(dist)[1:5])

score_filt <- score_NUMB[!names(score_NUMB) %in% central]

threshold <- median(score_filt, na.rm = TRUE)
group_score <- ifelse(score_filt > threshold, "high", "low")
names(group_score) <- names(score_filt)

group_df <- data.frame(
  cell_line = names(score_filt),
  score_NUMB = score_filt,
  group = group_score
)

# ============================================================
# 9) Figure 5A
# ============================================================

tiff(
  file.path(fig_dir, "Figure5A_score_NUMB_histogram.tiff"),
  width = 800, height = 600, units = "px", res = 150
)

hist(
  score_filt,
  breaks = 20,
  col = "steelblue",
  main = "Distribution of NUMB score",
  xlab = "score_NUMB",
  ylab = "Frequency"
)

abline(v = threshold, col = "red", lwd = 2)
dev.off()

# ============================================================
# 10) Figure 5B
# ============================================================

tiff(
  file.path(fig_dir, "Figure5B_score_NUMB_boxplot_filtered.tiff"),
  width = 800, height = 600, units = "px", res = 150
)

boxplot(
  score_filt ~ group,
  data = group_df,
  main = "score_NUMB by defined group",
  xlab = "score_NUMB group",
  ylab = "score_NUMB",
  col = c("red", "blue")
)

dev.off()

# ============================================================
# 11) LIMMA
# ============================================================

global_expr <- CCLE_log2[, expression_cols]
rownames(global_expr) <- CCLE_log2$transcript_id
global_expr <- global_expr[, names(group_score)]

groups <- factor(group_score, levels = c("low", "high"))
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(global_expr, design)
contrast <- makeContrasts(high_vs_low = high - low, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast))

limma_res <- topTable(fit2, number = Inf, adjust = "BH")
limma_res$transcript_id <- rownames(limma_res)

limma_res <- merge(
  limma_res,
  CCLE_log2[, c("transcript_id", "gene_name")],
  by = "transcript_id",
  all.x = TRUE
)

genes_interest <- genes_by_pathway$gene_symbol
limma_filt <- limma_res[limma_res$gene_name %in% genes_interest, ]

write.csv(limma_filt,
          file.path(proc_dir, "genes_significant_by_scoreNUMB_pathways.csv"),
          row.names = FALSE)

write.csv(group_df,
          file.path(proc_dir, "BCCL_grouped_by_NUMB_score.csv"),
          row.names = FALSE)

# ============================================================
# 12) Save expression matrices by group
# ============================================================

low_lines <- group_df$cell_line[group_df$group == "low"]
high_lines <- group_df$cell_line[group_df$group == "high"]

annotations <- CCLE_log2[, 1:3]

expr_low <- cbind(annotations, CCLE_log2[, low_lines, drop = FALSE])
expr_high <- cbind(annotations, CCLE_log2[, high_lines, drop = FALSE])

write.csv(expr_high,
          file.path(proc_dir, "gene_expression_scoreNUMB_high.csv"),
          row.names = FALSE)

write.csv(expr_low,
          file.path(proc_dir, "gene_expression_scoreNUMB_low.csv"),
          row.names = FALSE)

message("Figure 5A/5B NUMB score and LIMMA analysis completed.")

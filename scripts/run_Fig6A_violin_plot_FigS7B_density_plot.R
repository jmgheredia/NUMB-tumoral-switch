# ============================================================
# run_Fig6A_S7B_NUMB_score_TCGA.R
#
# Performs:
# - Loading isoform expression data from normal and tumor samples
# - Strict ordering of NUMB isoforms of interest
# - NUMB score calculation per sample
# - Comparative visualization (violin plot, density plot)
#
# Input:
#   - data/processed/normal_samples.csv
#   - data/processed/primary_tumor_samples.csv
#
# Output:
#   Figures:
#     - results/Figures/Figure_6/Fig_6A_violin_NUMB_score_TCGA.tiff
#     - results/Figures/Supplementary/S7/Fig_S7B_density_NUMB_score_TCGA.tiff
#
#   Processed data:
#     - results/processed/Figure_6/expr_mat_total_NUMB_score_TCGA.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# 1) Input files
# ============================================================

normal_csv <- file.path("data", "processed", "normal_samples.csv")
tumor_csv  <- file.path("data", "processed", "primary_tumor_samples.csv")

if (!file.exists(normal_csv)) stop("Missing normal_samples.csv in data/processed/")
if (!file.exists(tumor_csv))  stop("Missing primary_tumor_samples.csv in data/processed/")

# ============================================================
# 2) Output directories
# ============================================================

fig6_dir   <- file.path("results", "Figures", "Figure_6")
figS7_dir  <- file.path("results", "Figures", "Supplementary", "S7")
proc6_dir  <- file.path("results", "processed", "Figure_6")

dir.create(fig6_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figS7_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc6_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3) Load expression data
# ============================================================

normal <- read.csv(normal_csv, stringsAsFactors = FALSE)
tumor  <- read.csv(tumor_csv, stringsAsFactors = FALSE)

# ============================================================
# 4) Clean transcript IDs (remove version)
# ============================================================

normal$isoform_id <- sub("\\..*", "", normal$isoform_id)
tumor$isoform_id  <- sub("\\..*", "", tumor$isoform_id)

# ============================================================
# 5) Log2(TPM + 1) transformation
# ============================================================

expr_cols_normal <- setdiff(colnames(normal), "isoform_id")
expr_cols_tumor  <- setdiff(colnames(tumor), "isoform_id")

normal[, expr_cols_normal] <- log2(normal[, expr_cols_normal] + 1)
tumor[, expr_cols_tumor]   <- log2(tumor[, expr_cols_tumor] + 1)

# ============================================================
# 6) Define NUMB isoforms of interest (strict order)
# ============================================================

selected_isoforms <- c(
  "uc001xoa",  # p66 (suppressive)
  "uc002oon",  # NUMBL (suppressive)
  "uc001xob",  # p65 (suppressive)
  "uc001xnz",  # p71 (oncogenic)
  "uc001xny"   # p72 (oncogenic)
)

tumor_isoforms <- tumor[tumor$isoform_id %in% selected_isoforms, ]
tumor_isoforms <- tumor_isoforms[match(selected_isoforms, tumor_isoforms$isoform_id), ]

normal_isoforms <- normal[normal$isoform_id %in% selected_isoforms, ]
normal_isoforms <- normal_isoforms[match(selected_isoforms, normal_isoforms$isoform_id), ]

# ============================================================
# 7) Transpose matrices and assign isoform names
# ============================================================

rownames(tumor_isoforms)   <- tumor_isoforms$isoform_id
rownames(normal_isoforms) <- normal_isoforms$isoform_id

tumor_isoforms$isoform_id  <- NULL
normal_isoforms$isoform_id <- NULL

expr_mat_tumor  <- as.data.frame(t(tumor_isoforms))
expr_mat_normal <- as.data.frame(t(normal_isoforms))

colnames(expr_mat_tumor)  <- c("p66", "NUMBL", "p65", "p71", "p72")
colnames(expr_mat_normal) <- c("p66", "NUMBL", "p65", "p71", "p72")

expr_mat_tumor$Type  <- "Tumor"
expr_mat_normal$Type <- "Normal"

# ============================================================
# 8) Merge matrices and compute NUMB score
# ============================================================

expr_mat_total <- rbind(expr_mat_tumor, expr_mat_normal)

expr_mat_total$NUMB_score <-
  rowMeans(expr_mat_total[, c("p72", "p71")], na.rm = TRUE) -
  rowMeans(expr_mat_total[, c("p65", "p66", "NUMBL")], na.rm = TRUE)

# ============================================================
# 9) Add sample and patient identifiers
# ============================================================

expr_mat_total$sample_id <- rownames(expr_mat_total)
expr_mat_total$sample_id <- gsub("\\.", "-", expr_mat_total$sample_id)

expr_mat_total$patient_id <- sapply(
  strsplit(expr_mat_total$sample_id, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

# ============================================================
# 10) Figure 6A – Violin plot
# ============================================================

tiff(file.path(fig6_dir, "Fig_6A_violin_NUMB_score_TCGA.tiff"),
     width = 1200, height = 1000, res = 150)

ggplot(expr_mat_total, aes(x = Type, y = NUMB_score, fill = Type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, color = "black") +
  geom_jitter(width = 0.2, alpha = 0.4) +
  theme_minimal() +
  labs(
    title = "Distribution of NUMB score by group",
    y = "NUMB score",
    x = ""
  )

dev.off()

# ============================================================
# 11) Figure S7B – Density plot
# ============================================================

tiff(file.path(figS7_dir, "Fig_S7B_density_NUMB_score_TCGA.tiff"),
     width = 1200, height = 1000, res = 150)

ggplot(expr_mat_total, aes(x = NUMB_score, fill = Type)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "NUMB score density",
    x = "NUMB score",
    y = "Density"
  )

dev.off()

# ============================================================
# 12) Wilcoxon test Tumor vs Normal
# ============================================================

wilcox_test <- wilcox.test(NUMB_score ~ Type, data = expr_mat_total)
print(wilcox_test)

# ============================================================
# 13) Save processed data
# ============================================================

write.csv(
  expr_mat_total,
  file.path(proc6_dir, "expr_mat_total_NUMB_score_TCGA.csv"),
  row.names = FALSE
)

message("Figure 6A and Supplementary Figure S7B NUMB score analysis completed.")

# ============================================================
# run_Fig6B_Fig6C_NUMB_score_TCGA_clinical.R
#
# Generates:
# - Figure 6B: NUMB score distribution by AJCC clinical stage
#              including Non-tumor samples
# - Figure 6C: NUMB score distribution by PAM50 molecular subtype
#              including Non-tumor samples
#
# Performs:
# - Integration of NUMB score with TCGA BRCA clinical annotations
# - Grouping of samples by AJCC stage or PAM50 subtype
# - Comparative visualization using violin plots
# - Wilcoxon tests comparing each tumor group vs Non-tumor
#
# Input:
#   - data/metadata/TCGA.BRCA.sampleMap_BRCA_clinicalMatrix
#   - results/processed/Figure_6/expr_mat_total_NUMB_score_TCGA.csv
#
# Output:
#   Figures:
#     - results/Figures/Figure_6/Fig_6B_violin_NUMB_score_AJCC_NonTumor.tiff
#     - results/Figures/Figure_6/Fig_6C_violin_NUMB_score_PAM50_NonTumor.tiff
#
#   Processed data:
#     - results/processed/Figure_6/score_NUMB_with_AJCC_stage_and_normal.csv
#     - results/processed/Figure_6/score_NUMB_combined_PAM50_and_NonTumor.csv
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ============================================================
# 1) Input files
# ============================================================

clinical_csv <- file.path(
  "data", "metadata", "TCGA.BRCA.sampleMap_BRCA_clinicalMatrix"
)

score_csv <- file.path(
  "results", "processed", "Figure_6",
  "expr_mat_total_NUMB_score_TCGA.csv"
)

if (!file.exists(clinical_csv))
  stop("Missing TCGA clinical matrix in data/metadata/")

if (!file.exists(score_csv))
  stop("Missing expr_mat_total_NUMB_score_TCGA.csv in results/processed/Figure_6/")

# ============================================================
# 2) Output directories
# ============================================================

fig6_dir  <- file.path("results", "Figures", "Figure_6")
proc6_dir <- file.path("results", "processed", "Figure_6")

dir.create(fig6_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(proc6_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3) Load data
# ============================================================

clinical_matrix <- read.delim(
  clinical_csv,
  header = TRUE,
  stringsAsFactors = FALSE
)

score_NUMB <- read.csv(
  score_csv,
  stringsAsFactors = FALSE
)

# ============================================================
# FIGURE 6B — AJCC stage + Non-tumor
# ============================================================

# Extract patient ID for AJCC annotation
clinical_matrix$patient_id_ajcc <- substr(
  clinical_matrix$sampleID, 1, 12
)

clinical_stage <- clinical_matrix[, c(
  "patient_id_ajcc",
  "AJCC_Stage_nature2012"
)]

colnames(clinical_stage)[1] <- "patient_id"

clinical_stage$AJCC_Stage_nature2012 <-
  gsub("^Stage[[:space:]]+", "",
       clinical_stage$AJCC_Stage_nature2012)

# Merge NUMB score with AJCC stage
score_NUMB_AJCC <- merge(
  score_NUMB,
  clinical_stage,
  by = "patient_id",
  all.x = TRUE
)

# Define plotting groups
score_NUMB_AJCC$Group_plot <- ifelse(
  score_NUMB_AJCC$Type == "Normal",
  "Non-tumor",
  as.character(score_NUMB_AJCC$AJCC_Stage_nature2012)
)

score_NUMB_AJCC$Group_plot[
  is.na(score_NUMB_AJCC$Group_plot) &
    score_NUMB_AJCC$Type == "Tumor"
] <- "Unknown"

score_NUMB_clean <- score_NUMB_AJCC[
  !is.na(score_NUMB_AJCC$Group_plot) &
    score_NUMB_AJCC$Group_plot != "",
]

# Save processed AJCC data
write.csv(
  score_NUMB_clean,
  file.path(proc6_dir, "score_NUMB_with_AJCC_stage_and_normal.csv"),
  row.names = FALSE
)

# Order AJCC levels
ajcc_levels <- c(
  "Non-tumor",
  "I","IA","IB",
  "II","IIA","IIB",
  "III","IIIA","IIIB","IIIC",
  "IV",
  "X",
  "Unknown"
)

score_NUMB_clean$Group_plot <-
  factor(score_NUMB_clean$Group_plot, levels = ajcc_levels)

# Plot Figure 6B
tiff(
  file.path(fig6_dir, "Fig_6B_violin_NUMB_score_AJCC_NonTumor.tiff"),
  width = 1800, height = 1200, res = 150
)

ggplot(
  score_NUMB_clean,
  aes(x = Group_plot, y = NUMB_score, fill = Group_plot)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "NUMB score distribution by AJCC stage (including Non-tumor)",
    x = "Clinical group (AJCC stage)",
    y = "NUMB score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

dev.off()

# Wilcoxon tests: AJCC vs Non-tumor
df <- score_NUMB_clean
stages_to_test <- setdiff(levels(df$Group_plot), c("Non-tumor", NA))

wilcox_results <- lapply(stages_to_test, function(stg) {
  
  x <- df$NUMB_score[df$Group_plot == "Non-tumor"]
  y <- df$NUMB_score[df$Group_plot == stg]
  
  if (length(x) > 0 && length(y) > 0) {
    wt <- wilcox.test(y, x)
    
    data.frame(
      Stage = stg,
      N_stage = length(y),
      N_nontumor = length(x),
      W = wt$statistic,
      p_value = wt$p.value,
      stringsAsFactors = FALSE
    )
  }
})

df_wilcox <- do.call(rbind, wilcox_results)
df_wilcox$padj_BH <- p.adjust(df_wilcox$p_value, method = "BH")
print(df_wilcox)

# ============================================================
# FIGURE 6C — PAM50 subtype + Non-tumor
# ============================================================

clinical_matrix$patient_id_pam50 <- sapply(
  strsplit(clinical_matrix$sampleID, "-"),
  function(x) paste(x[1:3], collapse = "-")
)

pam50_data <- clinical_matrix[, c(
  "patient_id_pam50",
  "PAM50Call_RNAseq"
)]

colnames(pam50_data)[1] <- "patient_id"

score_NUMB_pam50 <- merge(
  score_NUMB,
  pam50_data,
  by = "patient_id",
  all.x = TRUE
)

score_NUMB_pam50$PAM50Call_RNAseq <-
  gsub("^Normal$", "Normal-like",
       score_NUMB_pam50$PAM50Call_RNAseq)

score_NUMB_pam50$Group_PAM50_Combined <- ifelse(
  score_NUMB_pam50$Type == "Normal",
  "Non-tumor",
  score_NUMB_pam50$PAM50Call_RNAseq
)

score_NUMB_pam50_clean <- score_NUMB_pam50[
  !is.na(score_NUMB_pam50$Group_PAM50_Combined) &
    score_NUMB_pam50$Group_PAM50_Combined != "",
]

pam50_levels <- c(
  "Non-tumor",
  "Basal",
  "Normal-like",
  "LumA",
  "LumB",
  "Her2"
)

score_NUMB_pam50_clean$Group_PAM50_Combined <-
  factor(
    score_NUMB_pam50_clean$Group_PAM50_Combined,
    levels = pam50_levels
  )

# Save processed PAM50 data
write.csv(
  score_NUMB_pam50_clean,
  file.path(proc6_dir, "score_NUMB_combined_PAM50_and_NonTumor.csv"),
  row.names = FALSE
)

# Plot Figure 6C
tiff(
  file.path(fig6_dir, "Fig_6C_violin_NUMB_score_PAM50_NonTumor.tiff"),
  width = 1600, height = 1200, res = 150
)

ggplot(
  score_NUMB_pam50_clean,
  aes(x = Group_PAM50_Combined, y = NUMB_score,
      fill = Group_PAM50_Combined)
) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "NUMB score distribution by PAM50 subtype (including Non-tumor)",
    x = "Clinical group",
    y = "NUMB score"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = guide_legend(title = "Group"))

dev.off()

message("Figure 6B (AJCC) and Figure 6C (PAM50) NUMB score analyses completed.")

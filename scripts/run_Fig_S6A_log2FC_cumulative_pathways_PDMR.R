# ============================================================
# run_Fig_S6A_log2FC_cumulative_pathways_PDMR.R
#
# Supplementary Figure S6A — Cumulative log2FC across PDMR passages
# using gene-level summaries from isoform-level differential expression
# and KEGG pathway gene sets derived from pathfindR outputs.
#
# This script:
# 1) Performs differential expression analysis (limma) between
#    consecutive PDMR passages and resection samples.
# 2) Aggregates transcript-level results to the gene level (max |logFC| per gene).
# 3) Reads pathfindR enrichment outputs and extracts genes for selected KEGG pathways.
# 4) Computes cumulative mean log2FC from Resection across passages.
# 5) Plots cumulative log2FC trajectories and saves a TIFF.
# 6) Computes pathway reactivation percentage from P0 to P4.
#
# Inputs:
# - data/metadata/metadata_PDMR.csv
# - data/processed/pdmr_isoform_expression_log2.csv
# - results/processed/Supplementary/S6/pathfindR_results_<COMPARISON>.csv (precomputed)
#
# Outputs:
# - results/processed/Supplementary/S6/Fig_s6A_combined_DE_results.csv
# - results/figures/Supplementary/S6/Fig_S6A_Cumulative_log2FC_pathways.tiff
# ============================================================

# REQUIRED LIBRARIES
library(dplyr)
library(tidyr)
library(limma)
library(ggplot2)
library(purrr)

# ============================================================
# 0) Input / Output
# ============================================================

metadata_file <- file.path("data", "metadata", "metadata_PDMR.csv")
expression_file <- file.path("data", "processed", "pdmr_isoform_expression_log2.csv")

out_dir_processed <- file.path("results", "processed", "Supplementary", "S6")
out_dir_figures <- file.path("results", "figures", "Supplementary", "S6")

dir.create(out_dir_processed, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figures, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(metadata_file),
  file.exists(expression_file)
)

# 1) LOAD DATA
metadata <- read.csv(metadata_file, sep = ";")
expression <- read.csv(expression_file, sep = ",", stringsAsFactors = FALSE, check.names = FALSE)

# 2) PREPROCESSING
expression$transcript_id <- sub("\\..*", "", expression$transcript_id)
expression <- expression[!duplicated(expression$transcript_id), ]
rownames(expression) <- expression$transcript_id

exp_ids <- expression[, c("transcript_id", "gene_id")]
exp_mat <- expression[, !(colnames(expression) %in% c("transcript_id", "gene_id"))]

### ---------------- COMPARISON 0 vs Resection ----------------
samples_R <- metadata %>% filter(Passage_of_this_sample == "Resection") %>%
  pull(PatientID.SpecimenID.SampleID)
samples_0 <- metadata %>% filter(Passage_of_this_sample == "0") %>%
  pull(PatientID.SpecimenID.SampleID)

group <- factor(ifelse(colnames(exp_mat) %in% samples_0, "P0",
                       ifelse(colnames(exp_mat) %in% samples_R, "R", NA)))
keep <- !is.na(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

exp_sub <- as.matrix(exp_mat[, keep])
fit <- lmFit(exp_sub, design)
contrast.matrix <- makeContrasts(P0 - R, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results_0_R <- topTable(fit2, number = Inf, sort.by = "none") %>%
  mutate(transcript_id = rownames(.)) %>%
  left_join(exp_ids, by = "transcript_id") %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  summarise(
    logFC_0_vs_R = logFC[which.max(abs(logFC))],
    p_value_0_vs_R = P.Value[which.max(abs(logFC))],
    .groups = "drop"
  )

### ---------------- COMPARISON 1 vs 0 ----------------
samples_1 <- metadata %>% filter(Passage_of_this_sample == "1") %>%
  pull(PatientID.SpecimenID.SampleID)

group <- factor(ifelse(colnames(exp_mat) %in% samples_1, "P1",
                       ifelse(colnames(exp_mat) %in% samples_0, "P0", NA)))
keep <- !is.na(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

exp_sub <- as.matrix(exp_mat[, keep])
fit <- lmFit(exp_sub, design)
contrast.matrix <- makeContrasts(P1 - P0, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results_1_0 <- topTable(fit2, number = Inf, sort.by = "none") %>%
  mutate(transcript_id = rownames(.)) %>%
  left_join(exp_ids, by = "transcript_id") %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  summarise(
    logFC_1_vs_0 = logFC[which.max(abs(logFC))],
    p_value_1_vs_0 = P.Value[which.max(abs(logFC))],
    .groups = "drop"
  )

### ---------------- COMPARISON 2 vs 1 ----------------
samples_2 <- metadata %>% filter(Passage_of_this_sample == "2") %>%
  pull(PatientID.SpecimenID.SampleID)

group <- factor(ifelse(colnames(exp_mat) %in% samples_2, "P2",
                       ifelse(colnames(exp_mat) %in% samples_1, "P1", NA)))
keep <- !is.na(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

exp_sub <- as.matrix(exp_mat[, keep])
fit <- lmFit(exp_sub, design)
contrast.matrix <- makeContrasts(P2 - P1, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results_2_1 <- topTable(fit2, number = Inf, sort.by = "none") %>%
  mutate(transcript_id = rownames(.)) %>%
  left_join(exp_ids, by = "transcript_id") %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  summarise(
    logFC_2_vs_1 = logFC[which.max(abs(logFC))],
    p_value_2_vs_1 = P.Value[which.max(abs(logFC))],
    .groups = "drop"
  )

### ---------------- COMPARISON 3 vs 2 ----------------
samples_3 <- metadata %>% filter(Passage_of_this_sample == "3") %>%
  pull(PatientID.SpecimenID.SampleID)

group <- factor(ifelse(colnames(exp_mat) %in% samples_3, "P3",
                       ifelse(colnames(exp_mat) %in% samples_2, "P2", NA)))
keep <- !is.na(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

exp_sub <- as.matrix(exp_mat[, keep])
fit <- lmFit(exp_sub, design)
contrast.matrix <- makeContrasts(P3 - P2, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results_3_2 <- topTable(fit2, number = Inf, sort.by = "none") %>%
  mutate(transcript_id = rownames(.)) %>%
  left_join(exp_ids, by = "transcript_id") %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  summarise(
    logFC_3_vs_2 = logFC[which.max(abs(logFC))],
    p_value_3_vs_2 = P.Value[which.max(abs(logFC))],
    .groups = "drop"
  )

### ---------------- COMPARISON 4 vs 3 ----------------
samples_4 <- metadata %>% filter(Passage_of_this_sample == "4") %>%
  pull(PatientID.SpecimenID.SampleID)

group <- factor(ifelse(colnames(exp_mat) %in% samples_4, "P4",
                       ifelse(colnames(exp_mat) %in% samples_3, "P3", NA)))
keep <- !is.na(group)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

exp_sub <- as.matrix(exp_mat[, keep])
fit <- lmFit(exp_sub, design)
contrast.matrix <- makeContrasts(P4 - P3, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

results_4_3 <- topTable(fit2, number = Inf, sort.by = "none") %>%
  mutate(transcript_id = rownames(.)) %>%
  left_join(exp_ids, by = "transcript_id") %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id) %>%
  summarise(
    logFC_4_vs_3 = logFC[which.max(abs(logFC))],
    p_value_4_vs_3 = P.Value[which.max(abs(logFC))],
    .groups = "drop"
  )

### ---------------- MERGE ALL RESULTS ----------------
combined_results <- list(
  results_0_R,
  results_1_0,
  results_2_1,
  results_3_2,
  results_4_3
) %>%
  reduce(full_join, by = "gene_id") %>%
  rename(Gene.symbol = gene_id)

write.csv(
  combined_results,
  file.path(out_dir_processed, "Fig_s6A_combined_DE_results.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------
# READ pathfindR OUTPUT FILES (PRECOMPUTED)
# ------------------------------------------------------------------

P0_vs_Resection <- read.csv(file.path(out_dir_processed, "pathfindR_results_P0_vs_resection.csv"), stringsAsFactors = FALSE, check.names = FALSE)
P1_vs_P0        <- read.csv(file.path(out_dir_processed, "pathfindR_results_P1_vs_P0.csv"),        stringsAsFactors = FALSE, check.names = FALSE)
P2_vs_P1        <- read.csv(file.path(out_dir_processed, "pathfindR_results_P2_vs_P1.csv"),        stringsAsFactors = FALSE, check.names = FALSE)
P3_vs_P2        <- read.csv(file.path(out_dir_processed, "pathfindR_results_P3_vs_P2.csv"),        stringsAsFactors = FALSE, check.names = FALSE)
P4_vs_P3        <- read.csv(file.path(out_dir_processed, "pathfindR_results_P4_vs_P3.csv"),        stringsAsFactors = FALSE, check.names = FALSE)

# ------------------ WNT PATHWAY (hsa04310) ------------------
wnt_genes_0_R <- c(
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04310", "Up_regulated"], ",")),
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04310", "Down_regulated"], ","))
)
wnt_genes_1_0 <- c(
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04310", "Up_regulated"], ",")),
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04310", "Down_regulated"], ","))
)
wnt_genes_2_1 <- c(
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04310", "Up_regulated"], ",")),
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04310", "Down_regulated"], ","))
)
wnt_genes_3_2 <- c(
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04310", "Up_regulated"], ",")),
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04310", "Down_regulated"], ","))
)
wnt_genes_4_3 <- c(
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04310", "Up_regulated"], ",")),
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04310", "Down_regulated"], ","))
)
genes_WNT_global <- sort(unique(trimws(c(
  wnt_genes_0_R, wnt_genes_1_0, wnt_genes_2_1, wnt_genes_3_2, wnt_genes_4_3
))))
wnt_logfc_df <- combined_results %>% filter(Gene.symbol %in% genes_WNT_global)
wnt_means <- data.frame(
  Comparison = c("P0_vs_R", "P1_vs_P0", "P2_vs_P1", "P3_vs_P2", "P4_vs_P3"),
  Mean_log2FC = c(
    mean(wnt_logfc_df$logFC_0_vs_R, na.rm = TRUE),
    mean(wnt_logfc_df$logFC_1_vs_0, na.rm = TRUE),
    mean(wnt_logfc_df$logFC_2_vs_1, na.rm = TRUE),
    mean(wnt_logfc_df$logFC_3_vs_2, na.rm = TRUE),
    mean(wnt_logfc_df$logFC_4_vs_3, na.rm = TRUE)
  )
)
wnt_cumulative <- data.frame(
  State = c("Resection", "P0", "P1", "P2", "P3", "P4"),
  Cumulative_log2FC = c(0, cumsum(wnt_means$Mean_log2FC))
)
wnt_cumulative$State <- factor(wnt_cumulative$State, levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------ NOTCH PATHWAY (hsa04330) ------------------
notch_genes_0_R <- c(
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04330", "Up_regulated"], ",")),
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04330", "Down_regulated"], ","))
)
notch_genes_1_0 <- c(
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04330", "Up_regulated"], ",")),
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04330", "Down_regulated"], ","))
)
notch_genes_2_1 <- c(
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04330", "Up_regulated"], ",")),
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04330", "Down_regulated"], ","))
)
notch_genes_3_2 <- c(
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04330", "Up_regulated"], ",")),
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04330", "Down_regulated"], ","))
)
notch_genes_4_3 <- c(
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04330", "Up_regulated"], ",")),
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04330", "Down_regulated"], ","))
)
genes_NOTCH_global <- sort(unique(trimws(c(
  notch_genes_0_R, notch_genes_1_0, notch_genes_2_1, notch_genes_3_2, notch_genes_4_3
))))
notch_logfc_df <- combined_results %>% filter(Gene.symbol %in% genes_NOTCH_global)
notch_means <- data.frame(
  Comparison = c("P0_vs_R", "P1_vs_P0", "P2_vs_P1", "P3_vs_P2", "P4_vs_P3"),
  Mean_log2FC = c(
    mean(notch_logfc_df$logFC_0_vs_R, na.rm = TRUE),
    mean(notch_logfc_df$logFC_1_vs_0, na.rm = TRUE),
    mean(notch_logfc_df$logFC_2_vs_1, na.rm = TRUE),
    mean(notch_logfc_df$logFC_3_vs_2, na.rm = TRUE),
    mean(notch_logfc_df$logFC_4_vs_3, na.rm = TRUE)
  )
)
notch_cumulative <- data.frame(
  State = c("Resection", "P0", "P1", "P2", "P3", "P4"),
  Cumulative_log2FC = c(0, cumsum(notch_means$Mean_log2FC))
)
notch_cumulative$State <- factor(notch_cumulative$State, levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------ HIPPO PATHWAY (hsa04390) ------------------
hippo_genes_0_R <- c(
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04390", "Up_regulated"], ",")),
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04390", "Down_regulated"], ","))
)
hippo_genes_1_0 <- c(
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04390", "Up_regulated"], ",")),
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04390", "Down_regulated"], ","))
)
hippo_genes_2_1 <- c(
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04390", "Up_regulated"], ",")),
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04390", "Down_regulated"], ","))
)
hippo_genes_3_2 <- c(
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04390", "Up_regulated"], ",")),
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04390", "Down_regulated"], ","))
)
hippo_genes_4_3 <- c(
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04390", "Up_regulated"], ",")),
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04390", "Down_regulated"], ","))
)
genes_HIPPO_global <- sort(unique(trimws(c(
  hippo_genes_0_R, hippo_genes_1_0, hippo_genes_2_1, hippo_genes_3_2, hippo_genes_4_3
))))
hippo_logfc_df <- combined_results %>% filter(Gene.symbol %in% genes_HIPPO_global)
hippo_means <- data.frame(
  Comparison = c("P0_vs_R", "P1_vs_P0", "P2_vs_P1", "P3_vs_P2", "P4_vs_P3"),
  Mean_log2FC = c(
    mean(hippo_logfc_df$logFC_0_vs_R, na.rm = TRUE),
    mean(hippo_logfc_df$logFC_1_vs_0, na.rm = TRUE),
    mean(hippo_logfc_df$logFC_2_vs_1, na.rm = TRUE),
    mean(hippo_logfc_df$logFC_3_vs_2, na.rm = TRUE),
    mean(hippo_logfc_df$logFC_4_vs_3, na.rm = TRUE)
  )
)
hippo_cumulative <- data.frame(
  State = c("Resection", "P0", "P1", "P2", "P3", "P4"),
  Cumulative_log2FC = c(0, cumsum(hippo_means$Mean_log2FC))
)
hippo_cumulative$State <- factor(hippo_cumulative$State, levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------ HEDGEHOG PATHWAY (hsa04340) ------------------
hh_genes_0_R <- c(
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04340", "Up_regulated"], ",")),
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04340", "Down_regulated"], ","))
)
hh_genes_1_0 <- c(
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04340", "Up_regulated"], ",")),
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04340", "Down_regulated"], ","))
)
hh_genes_2_1 <- c(
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04340", "Up_regulated"], ",")),
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04340", "Down_regulated"], ","))
)
hh_genes_3_2 <- c(
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04340", "Up_regulated"], ",")),
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04340", "Down_regulated"], ","))
)
hh_genes_4_3 <- c(
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04340", "Up_regulated"], ",")),
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04340", "Down_regulated"], ","))
)
genes_HEDGEHOG_global <- sort(unique(trimws(c(
  hh_genes_0_R, hh_genes_1_0, hh_genes_2_1, hh_genes_3_2, hh_genes_4_3
))))
hh_logfc_df <- combined_results %>% filter(Gene.symbol %in% genes_HEDGEHOG_global)
hh_means <- data.frame(
  Comparison = c("P0_vs_R", "P1_vs_P0", "P2_vs_P1", "P3_vs_P2", "P4_vs_P3"),
  Mean_log2FC = c(
    mean(hh_logfc_df$logFC_0_vs_R, na.rm = TRUE),
    mean(hh_logfc_df$logFC_1_vs_0, na.rm = TRUE),
    mean(hh_logfc_df$logFC_2_vs_1, na.rm = TRUE),
    mean(hh_logfc_df$logFC_3_vs_2, na.rm = TRUE),
    mean(hh_logfc_df$logFC_4_vs_3, na.rm = TRUE)
  )
)
hh_cumulative <- data.frame(
  State = c("Resection", "P0", "P1", "P2", "P3", "P4"),
  Cumulative_log2FC = c(0, cumsum(hh_means$Mean_log2FC))
)
hh_cumulative$State <- factor(hh_cumulative$State, levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------ OXIDATIVE PHOSPHORYLATION (hsa00190) ------------------
oxphos_genes_0_R <- c(
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa00190", "Up_regulated"], ",")),
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa00190", "Down_regulated"], ","))
)
oxphos_genes_1_0 <- c(
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa00190", "Up_regulated"], ",")),
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa00190", "Down_regulated"], ","))
)
oxphos_genes_2_1 <- c(
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa00190", "Up_regulated"], ",")),
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa00190", "Down_regulated"], ","))
)
oxphos_genes_3_2 <- c(
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa00190", "Up_regulated"], ",")),
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa00190", "Down_regulated"], ","))
)
oxphos_genes_4_3 <- c(
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa00190", "Up_regulated"], ",")),
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa00190", "Down_regulated"], ","))
)
genes_OXPHOS_global <- sort(unique(trimws(c(
  oxphos_genes_0_R, oxphos_genes_1_0, oxphos_genes_2_1, oxphos_genes_3_2, oxphos_genes_4_3
))))
oxphos_logfc_df <- combined_results %>% filter(Gene.symbol %in% genes_OXPHOS_global)
oxphos_means <- data.frame(
  Comparison = c("P0_vs_R", "P1_vs_P0", "P2_vs_P1", "P3_vs_P2", "P4_vs_P3"),
  Mean_log2FC = c(
    mean(oxphos_logfc_df$logFC_0_vs_R, na.rm = TRUE),
    mean(oxphos_logfc_df$logFC_1_vs_0, na.rm = TRUE),
    mean(oxphos_logfc_df$logFC_2_vs_1, na.rm = TRUE),
    mean(oxphos_logfc_df$logFC_3_vs_2, na.rm = TRUE),
    mean(oxphos_logfc_df$logFC_4_vs_3, na.rm = TRUE)
  )
)
oxphos_cumulative <- data.frame(
  State = c("Resection", "P0", "P1", "P2", "P3", "P4"),
  Cumulative_log2FC = c(0, cumsum(oxphos_means$Mean_log2FC))
)
oxphos_cumulative$State <- factor(oxphos_cumulative$State, levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------ CELL CYCLE (hsa04110) ------------------
cell_cycle_genes_0_R <- c(
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04110", "Up_regulated"], ",")),
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04110", "Down_regulated"], ","))
)
cell_cycle_genes_1_0 <- c(
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04110", "Up_regulated"], ",")),
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04110", "Down_regulated"], ","))
)
cell_cycle_genes_2_1 <- c(
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04110", "Up_regulated"], ",")),
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04110", "Down_regulated"], ","))
)
cell_cycle_genes_3_2 <- c(
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04110", "Up_regulated"], ",")),
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04110", "Down_regulated"], ","))
)
cell_cycle_genes_4_3 <- c(
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04110", "Up_regulated"], ",")),
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04110", "Down_regulated"], ","))
)
genes_cell_cycle <- sort(unique(trimws(c(
  cell_cycle_genes_0_R,
  cell_cycle_genes_1_0,
  cell_cycle_genes_2_1,
  cell_cycle_genes_3_2,
  cell_cycle_genes_4_3
))))
cell_cycle_logfc_df <- combined_results %>% filter(Gene.symbol %in% genes_cell_cycle)

cell_cycle_means <- data.frame(
  Comparison = c("P0_vs_R", "P1_vs_P0", "P2_vs_P1", "P3_vs_P2", "P4_vs_P3"),
  Mean_log2FC = c(
    mean(cell_cycle_logfc_df$logFC_0_vs_R, na.rm = TRUE),
    mean(cell_cycle_logfc_df$logFC_1_vs_0, na.rm = TRUE),
    mean(cell_cycle_logfc_df$logFC_2_vs_1, na.rm = TRUE),
    mean(cell_cycle_logfc_df$logFC_3_vs_2, na.rm = TRUE),
    mean(cell_cycle_logfc_df$logFC_4_vs_3, na.rm = TRUE)
  )
)
cell_cycle_cumulative <- data.frame(
  State = c("Resection", "P0", "P1", "P2", "P3", "P4"),
  Cumulative_log2FC = c(0, cumsum(cell_cycle_means$Mean_log2FC))
)
cell_cycle_cumulative$State <- factor(cell_cycle_cumulative$State, levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------ Th1 AND Th2 CELL DIFFERENTIATION (hsa04658) ------------------
th_genes_0_R <- c(
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04658", "Up_regulated"], ",")),
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa04658", "Down_regulated"], ","))
)
th_genes_1_0 <- c(
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04658", "Up_regulated"], ",")),
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa04658", "Down_regulated"], ","))
)
th_genes_2_1 <- c(
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04658", "Up_regulated"], ",")),
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa04658", "Down_regulated"], ","))
)
th_genes_3_2 <- c(
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04658", "Up_regulated"], ",")),
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa04658", "Down_regulated"], ","))
)
th_genes_4_3 <- c(
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04658", "Up_regulated"], ",")),
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa04658", "Down_regulated"], ","))
)

genes_th <- sort(unique(trimws(c(
  th_genes_0_R,
  th_genes_1_0,
  th_genes_2_1,
  th_genes_3_2,
  th_genes_4_3
))))

th_logfc_df <- combined_results %>%
  filter(Gene.symbol %in% genes_th)

th_means <- data.frame(
  Comparison = c("P0_vs_R", "P1_vs_P0", "P2_vs_P1", "P3_vs_P2", "P4_vs_P3"),
  Mean_log2FC = c(
    mean(th_logfc_df$logFC_0_vs_R, na.rm = TRUE),
    mean(th_logfc_df$logFC_1_vs_0, na.rm = TRUE),
    mean(th_logfc_df$logFC_2_vs_1, na.rm = TRUE),
    mean(th_logfc_df$logFC_3_vs_2, na.rm = TRUE),
    mean(th_logfc_df$logFC_4_vs_3, na.rm = TRUE)
  )
)

th_cumulative <- data.frame(
  State = c("Resection", "P0", "P1", "P2", "P3", "P4"),
  Cumulative_log2FC = c(0, cumsum(th_means$Mean_log2FC))
)

th_cumulative$State <- factor(th_cumulative$State,
                              levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------ PYRUVATE METABOLISM (hsa00620) ------------------
pyruvate_genes_0_R <- c(
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa00620", "Up_regulated"], ",")),
  unlist(strsplit(P0_vs_Resection[P0_vs_Resection$ID == "hsa00620", "Down_regulated"], ","))
)
pyruvate_genes_1_0 <- c(
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa00620", "Up_regulated"], ",")),
  unlist(strsplit(P1_vs_P0[P1_vs_P0$ID == "hsa00620", "Down_regulated"], ","))
)
pyruvate_genes_2_1 <- c(
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa00620", "Up_regulated"], ",")),
  unlist(strsplit(P2_vs_P1[P2_vs_P1$ID == "hsa00620", "Down_regulated"], ","))
)
pyruvate_genes_3_2 <- c(
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa00620", "Up_regulated"], ",")),
  unlist(strsplit(P3_vs_P2[P3_vs_P2$ID == "hsa00620", "Down_regulated"], ","))
)
pyruvate_genes_4_3 <- c(
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa00620", "Up_regulated"], ",")),
  unlist(strsplit(P4_vs_P3[P4_vs_P3$ID == "hsa00620", "Down_regulated"], ","))
)

genes_pyruvate <- sort(unique(trimws(c(
  pyruvate_genes_0_R,
  pyruvate_genes_1_0,
  pyruvate_genes_2_1,
  pyruvate_genes_3_2,
  pyruvate_genes_4_3
))))

pyruvate_logfc_df <- combined_results %>%
  filter(Gene.symbol %in% genes_pyruvate)

pyruvate_means <- c(
  mean(pyruvate_logfc_df$logFC_0_vs_R, na.rm = TRUE),
  mean(pyruvate_logfc_df$logFC_1_vs_0, na.rm = TRUE),
  mean(pyruvate_logfc_df$logFC_2_vs_1, na.rm = TRUE),
  mean(pyruvate_logfc_df$logFC_3_vs_2, na.rm = TRUE),
  mean(pyruvate_logfc_df$logFC_4_vs_3, na.rm = TRUE)
)

pyruvate_cumulative <- data.frame(
  State = c("Resection", "P0", "P1", "P2", "P3", "P4"),
  Cumulative_log2FC = c(0, cumsum(pyruvate_means))
)

pyruvate_cumulative$State <- factor(pyruvate_cumulative$State,
                                    levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------------------------------------------------------
# MERGE CUMULATIVE DATA FRAMES
# ------------------------------------------------------------------
results_all2 <- bind_rows(
  wnt_cumulative  %>% mutate(Pathway = "WNT"),
  notch_cumulative %>% mutate(Pathway = "Notch"),
  hippo_cumulative %>% mutate(Pathway = "Hippo"),
  hh_cumulative    %>% mutate(Pathway = "Hedgehog"),
  oxphos_cumulative %>% mutate(Pathway = "OXPHOS"),
  cell_cycle_cumulative %>% mutate(Pathway = "Cell cycle"),
  th_cumulative %>% mutate(Pathway = "Th1/Th2"),
  pyruvate_cumulative %>% mutate(Pathway = "Pyruvate")
)

results_all2$State <- factor(results_all2$State,
                             levels = c("Resection", "P0", "P1", "P2", "P3", "P4"))

# ------------------------------------------------------------------
# PATHWAY COLOR PALETTE
# ------------------------------------------------------------------
pathway_colors <- c(
  "WNT" = "#1f77b4",
  "Notch" = "#ff7f0e",
  "Hippo" = "#2ca02c",
  "Hedgehog" = "blue",
  "OXPHOS" = "#9467bd",
  "Cell cycle" = "#7f7f7f",
  "Th1/Th2" = "#ff9896",
  "Pyruvate" = "gold"
)

final_pathway_plot <- ggplot(results_all2, aes(x = State, y = Cumulative_log2FC, color = Pathway, group = Pathway)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = pathway_colors) +
  theme_classic(base_size = 12) +
  labs(
    title = "Cumulative log2FC from Resection",
    x = "Tumor state",
    y = "Cumulative log2FC",
    color = "KEGG pathway"
  )

print(final_pathway_plot)

ggsave(
  file.path(out_dir_figures, "Fig_S6A_Cumulative_log2FC_pathways.tiff"),
  plot = final_pathway_plot,
  device = "tiff",
  width = 10,
  height = 6,
  dpi = 300
)

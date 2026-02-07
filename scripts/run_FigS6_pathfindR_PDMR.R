# ============================================================
# run_FigureS5_pathfindR_PDMR.R
#
# Supplementary Figure S6 — Pathway enrichment analysis (pathfindR)
# based on isoform-level differential expression across PDMR passages.
#
# This script:
# 1) Performs differential expression analysis (limma) between
#    consecutive PDMR passages and resection samples.
# 2) Aggregates transcript-level results at the gene level.
# 3) Runs pathfindR enrichment analysis for each comparison.
#
# Inputs:
# - data/metadata/metadata_PDMR.csv
# - data/processed/pdmr_isoform_expression_log2.csv
#
# Outputs:
# - results/processed/Supplementary/S6/pathfindR_results_<COMPARISON>.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(limma)
  library(pathfindR)
})

# ============================================================
# 0) Input / Output
# ============================================================

metadata_file <- file.path("data", "metadata", "metadata_PDMR.csv")
expression_file <- file.path(
  "data", "processed", "pdmr_isoform_expression_log2.csv"
)

out_dir <- file.path("results", "processed", "Supplementary", "S6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(metadata_file),
  file.exists(expression_file)
)

# ============================================================
# 1) Load data
# ============================================================

metadata <- read.csv(metadata_file, sep = ";", stringsAsFactors = FALSE)
expression <- read.csv(expression_file, stringsAsFactors = FALSE, check.names = FALSE)

# ============================================================
# 2) Helper: normalize sample IDs (CRITICAL)
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
# 3) Preprocessing expression matrix
# ============================================================

expression$transcript_id <- sub("\\..*", "", expression$transcript_id)
expression <- expression[!duplicated(expression$transcript_id), ]
rownames(expression) <- expression$transcript_id

exp_ids <- expression[, c("transcript_id", "gene_id")]
exp_mat <- expression[, !(colnames(expression) %in% c("transcript_id", "gene_id"))]

# Normalize column names (sample IDs)
colnames(exp_mat) <- normalize_id(colnames(exp_mat))

# ============================================================
# 4) Preprocessing metadata
# ============================================================

stopifnot(
  "PatientID.SpecimenID.SampleID" %in% colnames(metadata),
  "Passage_of_this_sample" %in% colnames(metadata)
)

metadata$sample <- normalize_id(metadata$PatientID.SpecimenID.SampleID)

# ============================================================
# 5) Differential expression + pathfindR function
# ============================================================

run_comparison <- function(samples_A,
                           samples_B,
                           label_A,
                           label_B,
                           output_filename) {
  
  # Normalizar IDs (para quitar la X inicial si existe)
  normalize_id <- function(x) {
    x <- as.character(x)
    x <- sub("^X", "", x)
    x
  }
  
  samples_A <- normalize_id(samples_A)
  samples_B <- normalize_id(samples_B)
  colnames(exp_mat) <- normalize_id(colnames(exp_mat))
  
  group <- ifelse(
    colnames(exp_mat) %in% samples_A, label_A,
    ifelse(colnames(exp_mat) %in% samples_B, label_B, NA)
  )
  
  if (length(unique(na.omit(group))) < 2) {
    stop(
      sprintf(
        "Comparison %s vs %s has <2 groups after filtering.\nCheck sample name matching.",
        label_A, label_B
      )
    )
  }
  
  group <- factor(group)
  keep <- !is.na(group)
  
  design <- model.matrix(~ 0 + group[keep])
  colnames(design) <- levels(group[keep])
  
  exp_for_limma <- as.matrix(exp_mat[, keep])
  
  fit <- lmFit(exp_for_limma, design)
  
  contrast.matrix <- makeContrasts(
    contrasts = paste0(label_A, " - ", label_B),
    levels = design
  )
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  results <- topTable(fit2, number = Inf, sort.by = "none") %>%
    mutate(transcript_id = rownames(.)) %>%
    left_join(exp_ids, by = "transcript_id") %>%
    filter(!is.na(gene_id)) %>%
    filter(P.Value < 0.05)
  
  gene_summary <- results %>%
    mutate(direction = ifelse(logFC >= 0, "up", "down")) %>%
    group_by(gene_id, direction) %>%
    summarise(
      mean_logFC = mean(logFC),
      mean_p = mean(P.Value),
      n = n(),
      .groups = "drop"
    ) %>%
    group_by(gene_id) %>%
    filter(n == max(n)) %>%
    filter(abs(mean_logFC) == max(abs(mean_logFC))) %>%
    slice(1)
  
  final_df <- gene_summary %>%
    select(
      Gene.symbol = gene_id,
      logFC = mean_logFC,
      p = mean_p
    ) %>%
    filter(
      !is.na(Gene.symbol),
      Gene.symbol != "",
      grepl("^[A-Za-z0-9]+$", Gene.symbol)
    ) %>%
    mutate(across(c(logFC, p), as.numeric)) %>%
    as.data.frame()
  
  path_res <- run_pathfindR(final_df, n_processes = 2)
  
  write.csv(
    path_res,
    file.path(out_dir, output_filename),
    row.names = FALSE
  )
  
  return(path_res)
}

# ============================================================
# 6) Define sample groups (PDMR passages)
# ============================================================

samples_R <- metadata %>%
  filter(Passage_of_this_sample == "Resection") %>%
  pull(sample)

samples_0 <- metadata %>%
  filter(Passage_of_this_sample == "0") %>%
  pull(sample)

samples_1 <- metadata %>%
  filter(Passage_of_this_sample == "1") %>%
  pull(sample)

samples_2 <- metadata %>%
  filter(Passage_of_this_sample == "2") %>%
  pull(sample)

samples_3 <- metadata %>%
  filter(Passage_of_this_sample == "3") %>%
  pull(sample)

samples_4 <- metadata %>%
  filter(Passage_of_this_sample == "4") %>%
  pull(sample)

# ============================================================
# 7) Run comparisons (Supplementary Figure S6)
# ============================================================

run_comparison(
  samples_0, samples_R,
  "P0", "Resection",
  "pathfindR_results_P0_vs_Resection.csv"
)

run_comparison(
  samples_1, samples_0,
  "P1", "P0",
  "pathfindR_results_P1_vs_P0.csv"
)

run_comparison(
  samples_2, samples_1,
  "P2", "P1",
  "pathfindR_results_P2_vs_P1.csv"
)

run_comparison(
  samples_3, samples_2,
  "P3", "P2",
  "pathfindR_results_P3_vs_P2.csv"
)

run_comparison(
  samples_4, samples_3,
  "P4", "P3",
  "pathfindR_results_P4_vs_P3.csv"
)

message("Supplementary Figure S6 (pathfindR, PDMR) completed successfully.")

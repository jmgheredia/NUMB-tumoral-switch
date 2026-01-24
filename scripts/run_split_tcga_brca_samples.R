# run_split_tcga_brca_samples.R
#
# Purpose:
# Split TCGA-BRCA isoform-level expression data into
# Normal and Primary Tumor samples using curated sample metadata.
#
# Inputs:
# - data/external/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt
# - data/metadata/sample_type.csv
#
# Outputs:
# - data/processed/normal_samples.csv
# - data/processed/primary_tumor_samples.csv

# -----------------------------
# 0) Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(dplyr)
})

# -----------------------------
# 1) File paths
# -----------------------------
expr_file <- file.path(
  "data", "external",
  "BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt"
)

sample_file <- file.path(
  "data", "metadata",
  "sample_type.csv"
)

# -----------------------------
# 2) Sanity checks
# -----------------------------
if (!file.exists(expr_file)) {
  stop("Missing expression file: ", expr_file)
}

if (!file.exists(sample_file)) {
  stop("Missing sample metadata file: ", sample_file)
}

# -----------------------------
# 3) Load expression data
# -----------------------------
expr <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Rename first column (Hybridization REF) to isoform_id 
colnames(expr)[1] <- "isoform_id"

# Remove first row (TCGA metadata row)
expr <- expr[-1, ]

# Convert expression columns to numeric
expr[, -1] <- apply(expr[, -1], 2, as.numeric)

# Remove transcripts with zero expression across all samples
expr <- expr[rowSums(expr[, -1] != 0) > 0, ]

# -----------------------------
# 4) Load sample metadata
# -----------------------------
samples <- read.csv(
  sample_file,
  sep = ";",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Clean formatting
samples[, 1] <- gsub("\\s+", "", samples[, 1])
samples[, 2] <- gsub("\\s+", "", samples[, 2])

# -----------------------------
# 5) Define sample groups
# -----------------------------
normal_samples <- samples[samples[, 2] == "Normal", 1]
primary_tumor_samples <- samples[samples[, 2] == "Primary_tumor", 1]

# -----------------------------
# 6) Split expression matrix
# -----------------------------
normal_expr <- expr[, c(1, which(colnames(expr) %in% normal_samples)), drop = FALSE]
primary_tumor_expr <- expr[, c(1, which(colnames(expr) %in% primary_tumor_samples)), drop = FALSE]

# -----------------------------
# 7) Save outputs
# -----------------------------
out_dir <- file.path("data", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  normal_expr,
  file.path(out_dir, "normal_samples.csv"),
  row.names = FALSE
)

write.csv(
  primary_tumor_expr,
  file.path(out_dir, "primary_tumor_samples.csv"),
  row.names = FALSE
)

message("TCGA-BRCA sample split completed successfully.")

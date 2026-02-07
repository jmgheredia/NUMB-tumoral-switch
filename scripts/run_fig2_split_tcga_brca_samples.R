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
# 0) Input / Output
# -----------------------------
# This repository does not redistribute third-party expression data.
# Please download the expression matrix from Firebrowse and place it at:
# data/external/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt

expr_file <- file.path(
  "data", "external",
  "BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_isoforms_normalized__data.data.txt"
)

sample_file <- file.path(
  "data", "metadata",
  "sample_type.csv"
)

output_dir <- file.path("data", "processed")

# -----------------------------
# 1) Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(dplyr)
})

# -----------------------------
# 2) Load input data
# -----------------------------
if (!file.exists(expr_file)) {
  stop(
    paste0(
      "Missing required input file: ", expr_file, "\n",
      "Please download it from Firebrowse\n",
      "and place it in data/external/."
    )
  )
}

if (!file.exists(sample_file)) {
  stop("Missing sample metadata file: ", sample_file)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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
write.csv(
  normal_expr,
  file.path(output_dir, "normal_samples.csv"),
  row.names = FALSE
)

write.csv(
  primary_tumor_expr,
  file.path(output_dir, "primary_tumor_samples.csv"),
  row.names = FALSE
)

message("TCGA-BRCA sample split completed successfully.")


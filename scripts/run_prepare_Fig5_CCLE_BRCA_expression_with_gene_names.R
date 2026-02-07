# run_prepare_Fig5_CCLE_BRCA_expression_with_gene_names.R
#
# Generates CCLE BRCA expression matrix with Ensembl gene annotations
#
# Input:
#   - data/external/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt
#   - data/metadata/ensembl_gene_annotations.csv
#
# Output:
#   - data/processed/BCCL_gene_names_included.csv
#
# NOTE ON DATA AVAILABILITY:
# This repository does NOT redistribute the third-party CCLE expression file.
# Please download it from firebrowse.org and place it at:
#   data/external/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# 0) Input / Output
# ============================================================

ccle_txt <- file.path(
  "data", "external", "CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt"
)

gene_annotation_csv <- file.path(
  "data", "metadata", "ensembl_gene_annotations.csv"
)

out_dir <- file.path("data", "processed")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(ccle_txt)) {
  stop(
    paste0(
      "Missing required CCLE file: ", ccle_txt, "\n",
      "Please download it from firebrowse.org and place it in data/external/."
    )
  )
}

if (!file.exists(gene_annotation_csv)) {
  stop(
    paste0(
      "Missing gene annotation file: ", gene_annotation_csv, "\n",
      "Please place ensembl_gene_annotations.csv in data/metadata/."
    )
  )
}

# -----------------------------
# 2) Load CCLE expression data
# -----------------------------

ccle_raw <- read.delim(
  ccle_txt,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# -----------------------------
# 3) Subset BRCA cell lines
# -----------------------------

# Keep first two annotation columns (gene_id, transcript/isoform info)
base_df <- ccle_raw[, 1:2]

# Identify BREAST columns
columnas_breast <- grep("_BREAST$", colnames(ccle_raw), value = TRUE)

# Subset expression matrix for BREAST cell lines
ccle_brca_expr <- ccle_raw[, columnas_breast, drop = FALSE]

# Remove "_BREAST" suffix from cell line names
colnames(ccle_brca_expr) <- sub("_BREAST$", "", colnames(ccle_brca_expr))

# Combine annotation columns with BRCA expression matrix
ccle_brca <- cbind(
  base_df,
  ccle_brca_expr
)

# -----------------------------
# 4) Load gene annotations
# -----------------------------

gene_annotations <- read.csv(
  gene_annotation_csv,
  header = TRUE,
  sep = ";",
  stringsAsFactors = FALSE
)

required_cols <- c("gene_id", "gene_name")
if (!all(required_cols %in% colnames(gene_annotations))) {
  stop("ensembl_gene_annotations.csv must contain columns: gene_id, gene_name")
}

# -----------------------------
# 5) Merge expression with gene names
# -----------------------------

ccle_brca_annotated <- merge(
  ccle_brca,
  gene_annotations,
  by = "gene_id",
  all.x = TRUE
)

# Reorder columns: gene_name first, then everything else
ccle_brca_annotated <- ccle_brca_annotated[
  , c("gene_name", setdiff(colnames(ccle_brca_annotated), "gene_name"))
]

# -----------------------------
# 6) Save output
# -----------------------------

write.csv(
  ccle_brca_annotated,
  file = file.path(out_dir, "BCCL_gene_names_included.csv"),
  row.names = FALSE
)

message("CCLE BRCA expression matrix with gene names generated.")

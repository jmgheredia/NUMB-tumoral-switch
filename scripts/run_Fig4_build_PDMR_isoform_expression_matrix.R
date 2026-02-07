# ============================================================
# run_Fig4_build_PDMR_isoform_expression_matrix.R
#
# Construct a combined isoform-level expression matrix from RSEM *.results
# files derived from the PDMR resource, followed by log2(x+1) transformation.
#
# This script performs:
# - Recursive discovery of RSEM *.results files
# - Deterministic sample ordering based on filenames
# - Removal of transcripts with zero counts across all samples
# - Log2(x+1) transformation
#
# IMPORTANT:
# This repository does NOT redistribute original RSEM output files.
# Users must manually download the PDMR RSEM archives and place them under:
#
#   data/external/RSEM_archives/
#
# Inputs:
# - data/external/RSEM_archives/   (directory containing RSEM *.results files)
#
# Outputs:
# - data/processed/pdmr_isoform_expression_log2.csv
#
# ============================================================

# ============================================================
# 0) Input / Output
# ============================================================

rsem_root <- file.path("data", "external", "RSEM_archives")
out_dir   <- file.path("data", "processed")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(rsem_root)) {
  stop(
    paste0(
      "Missing required input directory: ", rsem_root, "\n",
      "Please download the PDMR RSEM archives and place them in data/external/RSEM_archives/."
    )
  )
}

# ============================================================
# Libraries
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

# ============================================================
# Parameters
# ============================================================

metric <- "expected_count"

# ============================================================
# Locate RSEM files (recursive + deterministic)
# ============================================================

files <- list.files(
  rsem_root,
  pattern = "\\.results$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No RSEM *.results files found under: ", rsem_root)
}

# Enforce deterministic ordering
files <- sort(files)

sample_names <- basename(files)

sample_names <- sample_names |>
  gsub("^X", "", x = _) |>
  gsub("[~_-]+", ".", x = _) |>
  gsub("\\.+", ".", x = _)

stopifnot(!anyDuplicated(sample_names))

# ============================================================
# Read RSEM files
# ============================================================

read_one <- function(f, sname) {
  
  df <- read_tsv(f, show_col_types = FALSE)
  
  df %>%
    dplyr::select(transcript_id, gene_id, all_of(metric)) %>%
    dplyr::rename(!!sname := all_of(metric))
}

expr_list <- map2(files, sample_names, read_one)

expr <- reduce(expr_list, full_join, by = c("transcript_id", "gene_id"))

# ============================================================
# Drop rows with all-zero counts
# ============================================================

counts <- expr[, -(1:2)]
keep <- rowSums(counts, na.rm = TRUE) > 0
expr <- expr[keep, ]

# ============================================================
# log2(x + 1)
# ============================================================

logmat <- log2(as.matrix(expr[, -(1:2)]) + 1)

expr_final <- as.data.frame(logmat)
expr_final$transcript_id <- expr$transcript_id
expr_final$gene_id <- expr$gene_id

expr_final <- expr_final[, c(
  "transcript_id",
  "gene_id",
  setdiff(colnames(expr_final), c("transcript_id", "gene_id"))
)]

# ============================================================
# Write output
# ============================================================

out_file <- file.path(out_dir, "pdmr_isoform_expression_log2.csv")

write.csv(
  expr_final,
  out_file,
  row.names = FALSE
)

message("DONE: ", out_file)

# ============================================================
# run_prepare_Fig5E_GDSC_CCLE_BRCA.R
#
# Filters GDSC1 and GDSC2 drug response datasets to CCLE BRCA cell lines
# used in the NUMB-score analysis.
#
# Input:
#   - data/external/GDSC1.csv
#   - data/external/GDSC2.csv
#   - data/processed/BCCL_gene_names_included.csv
#
# Output:
#   - data/processed/GDSC1_breast.csv
#   - data/processed/GDSC2_breast.csv
#
# NOTE ON DATA AVAILABILITY:
# This repository does NOT redistribute GDSC drug response datasets.
# Please download GDSC1 and GDSC2 files from:
#   https://www.cancerrxgene.org/
# and place them in data/external/ as GDSC1.csv and GDSC2.csv.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ============================================================
# 0) Input / Output
# ============================================================

gdsc1_csv <- file.path("data", "external", "GDSC1.csv")
gdsc2_csv <- file.path("data", "external", "GDSC2.csv")
ccle_csv  <- file.path("data", "processed", "BCCL_gene_names_included.csv")

out_dir <- file.path("data", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(gdsc1_csv)) stop("Missing GDSC1.csv in data/external/")
if (!file.exists(gdsc2_csv)) stop("Missing GDSC2.csv in data/external/")
if (!file.exists(ccle_csv))  stop("Missing BCCL_gene_names_included.csv in data/processed/")

# ============================================================
# 2) Load data
# ============================================================

GDSC1 <- read.csv(gdsc1_csv, sep = ";", stringsAsFactors = FALSE)
GDSC2 <- read.csv(gdsc2_csv, sep = ";", stringsAsFactors = FALSE)
CCLE  <- read.csv(ccle_csv, stringsAsFactors = FALSE)

# ============================================================
# 3) Filter GDSC by TCGA BRCA
# ============================================================

GDSC1_breast <- GDSC1 %>%
  filter(TCGA_DESC == "BRCA")

GDSC2_breast <- GDSC2 %>%
  filter(TCGA_DESC == "BRCA")

# ============================================================
# 4) Save outputs
# ============================================================

write.csv(GDSC1_breast,
          file.path(out_dir, "GDSC1_breast.csv"),
          row.names = FALSE)

write.csv(GDSC2_breast,
          file.path(out_dir, "GDSC2_breast.csv"),
          row.names = FALSE)

message("GDSC BRCA subsets generated.")

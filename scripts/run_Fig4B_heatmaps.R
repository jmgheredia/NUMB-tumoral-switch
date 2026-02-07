# ============================================================
# run_Figure4B_heatmaps.R
#
# Generate Figure 4B Pearson correlation heatmaps (Resection / 0 / 4 only)
# directly from log2-normalized isoform expression.
#
# Inputs:
# - data/processed/pdmr_isoform_expression_log2.csv
# - data/metadata/metadata_PDMR.csv
# - data/metadata/genes_by_pathway.csv
#
# Outputs:
# - results/Figures/Figure_4/heatmap_correlation_<PATHWAY>_R0_4.tiff
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

# ============================================================
# 0) Input / Output
# ============================================================

expr_file  <- file.path("data", "processed", "pdmr_isoform_expression_log2.csv")
meta_file  <- file.path("data", "metadata", "metadata_PDMR.csv")
genes_file <- file.path("data", "metadata", "genes_by_pathway.csv")

fig_dir <- file.path("results", "Figures", "Figure_4")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(
  file.exists(expr_file),
  file.exists(meta_file),
  file.exists(genes_file)
)

# ============================================================
# Load data
# ============================================================

entire_dataset    <- read.csv(expr_file, header = TRUE, check.names = FALSE)
metadata          <- read.csv(meta_file, sep = ";")
genes_by_pathway <- read.csv(genes_file, row.names = 1)

# ============================================================
# Basic preprocessing
# ============================================================

entire_dataset$transcript_id <- sub("\\..*", "", entire_dataset$transcript_id)

entire_dataset$gene_id <- as.character(entire_dataset$gene_id)
genes_by_pathway$gene_symbol <- as.character(genes_by_pathway$gene_symbol)

metadata$PatientID.SpecimenID.SampleID <-
  as.character(metadata$PatientID.SpecimenID.SampleID)

metadata$PatientID.SpecimenID.SampleID <-
  trimws(metadata$PatientID.SpecimenID.SampleID)

colnames(entire_dataset) <- trimws(colnames(entire_dataset))

stopifnot("Passage_of_this_sample" %in% colnames(metadata))

# ============================================================
# Merge genes and transcripts
# ============================================================

genes_transcripts_by_pathway <- merge(
  genes_by_pathway,
  entire_dataset,
  by.x = "gene_symbol",
  by.y = "gene_id"
)

genes_transcripts_by_pathway <-
  genes_transcripts_by_pathway[, c("gene_symbol", "transcript_id", "pathway")]

notch_df    <- subset(genes_transcripts_by_pathway, pathway == "Notch")
hippo_df    <- subset(genes_transcripts_by_pathway, pathway == "Hippo")
wnt_df      <- subset(genes_transcripts_by_pathway, pathway == "WNT")
hedgehog_df <- subset(genes_transcripts_by_pathway, pathway == "Hedgehog")

# ============================================================
# Isoforms of interest
# ============================================================

numb_isoforms <- list(
  p65   = "uc001xob",
  p66   = "uc001xoa",
  p71   = "uc001xnz",
  p72   = "uc001xny",
  NUMBL = "uc002oon"
)

# ============================================================
# Passages (paper only)
# ============================================================

passages <- c("Resection", "0", "4")

# ============================================================
# Correlation function
# ============================================================

correlations_by_pathway_by_passage <- function(pathway_df,
                                               entire_dataset,
                                               numb_isoforms,
                                               metadata,
                                               specific_passage) {
  
  valid_samples <- metadata %>%
    filter(Passage_of_this_sample == specific_passage) %>%
    pull(PatientID.SpecimenID.SampleID)
  
  expr_long <- entire_dataset %>%
    pivot_longer(
      -c(transcript_id, gene_id),
      names_to = "SampleID",
      values_to = "expression"
    ) %>%
    filter(SampleID %in% valid_samples)
  
  results <- list()
  
  for (iso_name in names(numb_isoforms)) {
    
    iso_id <- numb_isoforms[[iso_name]]
    
    iso_expr <- expr_long %>%
      filter(transcript_id == iso_id) %>%
      dplyr::select(SampleID, expression) %>%
      dplyr::rename(iso_expr = expression)
    
    for (trans_id in unique(pathway_df$transcript_id)) {
      
      trans_expr <- expr_long %>%
        filter(transcript_id == trans_id) %>%
        dplyr::select(SampleID, expression) %>%
        dplyr::rename(trans_expr = expression)
      
      joined <- inner_join(iso_expr, trans_expr, by = "SampleID")
      
      if (nrow(joined) >= 3) {
        corr_val <- cor(joined$iso_expr, joined$trans_expr, method = "pearson")
      } else {
        corr_val <- NA
      }
      
      results[[length(results) + 1]] <- data.frame(
        isoform = iso_name,
        isoform_id = iso_id,
        transcript_id = trans_id,
        correlation = corr_val,
        passage = specific_passage
      )
    }
  }
  
  do.call(rbind, results)
}

# ============================================================
# Process pathway
# ============================================================

process_pathway <- function(df_pathway) {
  
  res <- do.call(rbind, lapply(passages, function(p) {
    correlations_by_pathway_by_passage(
      df_pathway,
      entire_dataset,
      numb_isoforms,
      metadata,
      p
    )
  }))
  
  res <- merge(
    res,
    genes_transcripts_by_pathway,
    by = "transcript_id",
    all.x = TRUE
  )
  
  res[, c("isoform", "isoform_id", "transcript_id",
          "gene_symbol", "correlation", "pathway", "passage")]
}

# ============================================================
# Compute correlations
# ============================================================

notch_corr    <- process_pathway(notch_df)    %>% filter(pathway == "Notch")
hippo_corr    <- process_pathway(hippo_df)    %>% filter(pathway == "Hippo")
wnt_corr      <- process_pathway(wnt_df)      %>% filter(pathway == "WNT")
hedgehog_corr <- process_pathway(hedgehog_df) %>% filter(pathway == "Hedgehog")

# ============================================================
# Heatmap column order (R / 0 / 4)
# ============================================================

filtered_column_order <- c(
  "p72_Resection","p72_0","p72_4",
  "p71_Resection","p71_0","p71_4",
  "p66_Resection","p66_0","p66_4",
  "p65_Resection","p65_0","p65_4",
  "NUMBL_Resection","NUMBL_0","NUMBL_4"
)

# ============================================================
# Heatmap function
# ============================================================

make_heatmap <- function(df_corr, pathway_name) {
  
  df_corr <- df_corr %>% filter(passage %in% passages)
  
  table <- df_corr %>%
    mutate(iso_passage = paste0(isoform, "_", passage)) %>%
    dplyr::select(transcript_id, iso_passage, correlation) %>%
    pivot_wider(names_from = iso_passage, values_from = correlation) %>%
    distinct(transcript_id, .keep_all = TRUE) %>%
    arrange(transcript_id)
  
  mat <- as.matrix(table[, -1])
  rownames(mat) <- table$transcript_id
  
  mat <- mat[, intersect(filtered_column_order, colnames(mat)), drop = FALSE]
  mat <- na.omit(mat)
  
  if (nrow(mat) < 2) {
    warning(paste("Skipping", pathway_name, "- insufficient data"))
    return(invisible(NULL))
  }
  
  outfile <- file.path(
    fig_dir,
    paste0("heatmap_correlation_", toupper(pathway_name), "_R0_4.tiff")
  )
  
  tiff(outfile, width = 2000, height = 2500, res = 300)
  
  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = paste("Correlations -", pathway_name, "(R / 0 / 4)"),
    fontsize_row = 6,
    fontsize_col = 8
  )
  
  dev.off()
}

# ============================================================
# Generate Figure 4B heatmaps
# ============================================================

make_heatmap(notch_corr,    "Notch")
make_heatmap(hippo_corr,    "Hippo")
make_heatmap(wnt_corr,      "WNT")
make_heatmap(hedgehog_corr, "Hedgehog")

message("Figure 4B heatmaps generated.")

# run_isoform_histograms_TCGA_BRCA.R
#
# Generate expression and paired log2FC histograms for NUMB / NUMBL isoforms in TCGA-BRCA, including:
#   - Full cohort Tumor expression histograms (Supplementary S3)
#   - Paired Normal vs Tumor expression histograms (Figure 2A)
#   - Paired log2FC histograms (Figure 2B)
# Additionally:
#   - Dip test for multimodality (paired log2FC)
#   - High/low split (k-means if multimodal; terciles otherwise)
#   - Differential tables, split tables, fused differentials, correlations, and heatmap (Figure 2D)
#
# Inputs:
# - data/processed/normal_samples.csv
# - data/processed/primary_tumor_samples.csv
#
# Outputs:
# - results/figures/Figure_2/
# - results/figures/Supplementary/S3/
# - results/processed/Figure_2/
# - results/processed/Figure_2/diptest_analysis/{differential,splits}

# -----------------------------
# 0) Input / Output
# -----------------------------
normal_file <- file.path("data", "processed", "normal_samples.csv")
tumor_file  <- file.path("data", "processed", "primary_tumor_samples.csv")

fig2_dir <- file.path("results","figures","Figure_2")

fig_expr <- fig2_dir
fig_fc   <- fig2_dir
fig_2d   <- fig2_dir

supp_dir <- file.path("results","figures","Supplementary","S3")

data_dir   <- file.path("results","processed","Figure_2")
dip_dir    <- file.path(data_dir,"diptest_analysis")
diff_dir   <- file.path(dip_dir,"differential")
splits_dir <- file.path(dip_dir,"splits")

# -----------------------------
# 1) Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(diptest)
})

# -----------------------------
# 2) Sanity checks + dirs
# -----------------------------
if (!file.exists(normal_file)) stop("Missing file: ", normal_file)
if (!file.exists(tumor_file))  stop("Missing file: ", tumor_file)

dir.create(fig_expr, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_fc,   recursive = TRUE, showWarnings = FALSE)
dir.create(fig_2d,   recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(diff_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(splits_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 3) Load data
# -----------------------------
normal_df <- read.csv(normal_file, stringsAsFactors = FALSE, check.names = FALSE)
tumor_df  <- read.csv(tumor_file,  stringsAsFactors = FALSE, check.names = FALSE)

colnames(normal_df)[1] <- "isoform_id"
colnames(tumor_df)[1]  <- "isoform_id"

normal_df$isoform_id <- sub("\\..*", "", normal_df$isoform_id)
tumor_df$isoform_id  <- sub("\\..*", "", tumor_df$isoform_id)

# -----------------------------
# 4) Expression columns
# -----------------------------
expr_cols_normal <- setdiff(colnames(normal_df), "isoform_id")
expr_cols_tumor  <- setdiff(colnames(tumor_df),  "isoform_id")

# -----------------------------
# 5) log2(x + 1)
# -----------------------------
normal_df[, expr_cols_normal] <- log2(normal_df[, expr_cols_normal] + 1)
tumor_df[,  expr_cols_tumor]  <- log2(tumor_df[,  expr_cols_tumor]  + 1)

# -----------------------------
# 6) Isoforms
# -----------------------------
isoforms <- list(
  p65   = "uc001xob",
  p66   = "uc001xoa",
  p71   = "uc001xnz",
  p72   = "uc001xny",
  NUMBL = "uc002oon"
)

# -----------------------------
# 7) Styles
# -----------------------------
style_full  <- list(bar = "lightpink", line = "darkred",  border = "white")
style_npair <- list(bar = "lightblue", line = "darkblue", border = "white")
style_tpair <- list(bar = "lightpink", line = "darkred",  border = "white")
style_fc    <- list(bar = "lightblue", line = "darkblue", border = "white")

hist_breaks    <- 20
hist_breaks_fc <- 30

# ============================================================
# 8) PAIRED SAMPLE IDENTIFICATION
# ============================================================
extract_patient <- function(x) sapply(strsplit(x, "-"), function(z) paste(z[1:3], collapse = "-"))
extract_type    <- function(x) sapply(strsplit(x, "-"), function(z) z[4])

meta_normal <- data.frame(column=expr_cols_normal,patient=extract_patient(expr_cols_normal),type=extract_type(expr_cols_normal))
meta_tumor  <- data.frame(column=expr_cols_tumor, patient=extract_patient(expr_cols_tumor), type=extract_type(expr_cols_tumor))

meta_normal <- subset(meta_normal, type=="11A")
meta_tumor  <- subset(meta_tumor, type %in% c("01A","01B"))

paired_patients <- intersect(meta_normal$patient, meta_tumor$patient)

meta_normal <- meta_normal[meta_normal$patient %in% paired_patients,]
meta_tumor  <- meta_tumor[ meta_tumor$patient  %in% paired_patients,]

meta_normal <- meta_normal[!duplicated(meta_normal$patient),]
meta_tumor  <- meta_tumor[ !duplicated(meta_tumor$patient),]

meta_normal <- meta_normal[order(meta_normal$patient),]
meta_tumor  <- meta_tumor[ order(meta_tumor$patient),]

# ============================================================
# 9) PAIRED MATRICES
# ============================================================
normal_p <- normal_df[,c("isoform_id",meta_normal$column)]
tumor_p  <- tumor_df[, c("isoform_id",meta_tumor$column)]

colnames(normal_p)[-1] <- meta_normal$patient
colnames(tumor_p)[-1]  <- meta_tumor$patient

stopifnot(all(colnames(normal_p)[-1]==colnames(tumor_p)[-1]))

write.csv(normal_p,file.path(data_dir,"normal_paired.csv"),row.names=FALSE)
write.csv(tumor_p, file.path(data_dir,"tumor_paired.csv"), row.names=FALSE)

# ============================================================
# 10) FIGURE 2A
# ============================================================
for (iso in names(isoforms)) {
  
  tx <- isoforms[[iso]]
  row_n <- normal_p[normal_p$isoform_id==tx,]
  row_t <- tumor_p[ tumor_p$isoform_id ==tx,]
  if(nrow(row_n)==0||nrow(row_t)==0) next
  
  v_n <- as.numeric(row_n[1,-1])
  v_t <- as.numeric(row_t[1,-1])
  
  out_n <- file.path(fig_expr,paste0("histogram_NORMAL_paired_",iso,"_",tx,".tiff"))
  out_t <- file.path(fig_expr,paste0("histogram_TUMOR_paired_", iso,"_",tx,".tiff"))
  
  tiff(out_n, width = 1600, height = 1200, res = 200, compression = "lzw")
  hist(
    v_n,
    breaks = hist_breaks,
    prob = TRUE,
    main = paste("Normal (paired) -", iso),
    xlab = "log2(TPM + 1)",
    col = style_npair$bar,
    border = style_npair$border
  )
  lines(density(v_n, na.rm = TRUE), col = style_npair$line, lwd = 2)
  rug(v_n)
  dev.off()
  
  tiff(out_t, width = 1600, height = 1200, res = 200, compression = "lzw")
  hist(
    v_t,
    breaks = hist_breaks,
    prob = TRUE,
    main = paste("Tumor (paired) -", iso),
    xlab = "log2(TPM + 1)",
    col = style_tpair$bar,
    border = style_tpair$border
  )
  lines(density(v_t, na.rm = TRUE), col = style_tpair$line, lwd = 2)
  rug(v_t)
  dev.off()
}

# ============================================================
# 11) FIGURE 2B
# ============================================================
common_tx <- intersect(normal_p$isoform_id,tumor_p$isoform_id)
normal_p2 <- normal_p[normal_p$isoform_id%in%common_tx,]
tumor_p2  <- tumor_p[ tumor_p$isoform_id %in%common_tx,]

normal_p2 <- normal_p2[order(normal_p2$isoform_id),]
tumor_p2  <- tumor_p2[ order(tumor_p2$isoform_id),]

log2FC <- tumor_p2; log2FC[,-1] <- tumor_p2[,-1]-normal_p2[,-1]

hist_data <- data.frame()

for(iso in names(isoforms)){
  tx<-isoforms[[iso]]
  row_tx<-log2FC[log2FC$isoform_id==tx,]; if(nrow(row_tx)==0) next
  values<-as.numeric(row_tx[1,-1])
  hist_data<-rbind(hist_data,data.frame(isoform=iso,transcript_id=tx,patient=colnames(row_tx)[-1],log2FC=values))
  out <- file.path(fig_fc, paste0("histogram_log2FC_paired_", iso, "_", tx, ".tiff"))
  tiff(out, width = 1600, height = 1200, res = 200, compression = "lzw")
  hist(
    values,
    breaks = hist_breaks_fc,
    prob = TRUE,
    main = paste("Paired log2FC -", iso),
    xlab = "log2FC",
    col = style_fc$bar,
    border = style_fc$border
  )
  lines(density(values, na.rm = TRUE), col = style_fc$line, lwd = 2)
  rug(values)
  dev.off()
}

write.csv(hist_data,file.path(data_dir,"log2FC_paired_isoforms_fig2B.csv"),row.names=FALSE)

# ============================================================
# 12) SUPPLEMENTARY S3
# ============================================================
for(iso in names(isoforms)){
  tx<-isoforms[[iso]]
  row_t<-tumor_df[tumor_df$isoform_id==tx,]
  if(nrow(row_t)>0){
    values<-as.numeric(row_t[1,expr_cols_tumor])
    out <- file.path(supp_dir, paste0("histogram_TUMOR_full_", iso, "_", tx, ".tiff"))
    tiff(out, width = 1600, height = 1200, res = 200, compression = "lzw")
    hist(
      values,
      breaks = hist_breaks,
      prob = TRUE,
      main = paste("Tumor full -", iso),
      xlab = "log2(TPM + 1)",
      col = style_full$bar,
      border = style_full$border
    )
    lines(density(values, na.rm = TRUE), col = style_full$line, lwd = 2)
    rug(values)
    dev.off()
  }
}

# ============================================================
# 13) DIPTEST ON PAIRED log2FC + HIGH/LOW SPLIT + DIFFERENTIAL
# ============================================================

message("Starting diptest on paired log2FC...")

for (iso in names(isoforms)) {
  
  tx <- isoforms[[iso]]
  
  row_tx <- log2FC[log2FC$isoform_id == tx, ]
  
  if (nrow(row_tx) == 0) next
  
  log2fc_values <- as.numeric(row_tx[1, -1])
  patients <- colnames(row_tx)[-1]
  names(log2fc_values) <- patients
  
  message("Diptest for ", iso, " (", tx, ")")
  
  p_dip <- dip.test(log2fc_values)$p.value
  
  if (p_dip < 0.05) {
    clustering <- kmeans(log2fc_values, centers = 2)
    high_group <- names(log2fc_values)[clustering$cluster == which.max(clustering$centers)]
    low_group  <- names(log2fc_values)[clustering$cluster == which.min(clustering$centers)]
  } else {
    terciles <- quantile(log2fc_values, probs = c(1/3, 2/3), na.rm = TRUE)
    low_group  <- names(log2fc_values[log2fc_values <= terciles[1]])
    high_group <- names(log2fc_values[log2fc_values >= terciles[2]])
  }
  
  message("Patients high: ", length(high_group), " | low: ", length(low_group))
  
  if (length(high_group) < 2 || length(low_group) < 2) {
    message("Skipping ", iso, " due to insufficient groups")
    next
  }
  
  genes <- log2FC$isoform_id
  gene_log2fc <- numeric(length(genes))
  p_vals <- numeric(length(genes))
  
  for (i in seq_along(genes)) {
    
    expr_high <- as.numeric(log2FC[i, high_group])
    expr_low  <- as.numeric(log2FC[i, low_group])
    
    tt <- try(t.test(expr_high, expr_low), silent = TRUE)
    
    if (inherits(tt, "htest")) {
      gene_log2fc[i] <- mean(expr_high, na.rm = TRUE) - mean(expr_low, na.rm = TRUE)
      p_vals[i] <- tt$p.value
    } else {
      gene_log2fc[i] <- NA
      p_vals[i] <- NA
    }
  }
  
  differential_df <- data.frame(
    transcript_id = genes,
    log2FC = gene_log2fc,
    p_value = p_vals,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    differential_df,
    file.path(diff_dir, paste0("diferencial_pareado_", iso, "_", tx, ".csv")),
    row.names = FALSE
  )
  
  split_df <- data.frame(
    patient = c(high_group, low_group),
    group = c(rep("high", length(high_group)), rep("low", length(low_group))),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    split_df,
    file.path(splits_dir, paste0("split_pareado_", iso, "_", tx, ".csv")),
    row.names = FALSE
  )
}

message("Paired diptest + differential analysis completed.")

# ============================================================
# 14) TRANSCRIPTOMIC PROFILE CORRELATION (PAIRED)
# ============================================================

load_log2fc <- function(filename, new_col) {
  df <- read.csv(filename, header = TRUE, stringsAsFactors = FALSE)
  df <- df[, c("transcript_id", "log2FC")]
  colnames(df)[2] <- new_col
  df
}

isoforms_corr <- c("p72", "p71", "p65", "p66", "numbl")
codes_corr    <- c("uc001xny", "uc001xnz", "uc001xob", "uc001xoa", "uc002oon")

analyze_log2fc <- function(type = "PAIRED") {
  
  files <- file.path(
    diff_dir,
    paste0("diferencial_pareado_", isoforms_corr, "_", codes_corr, ".csv")
  )
  
  col_names <- paste0("log2FC_", isoforms_corr)
  
  dfs <- mapply(load_log2fc, files, col_names, SIMPLIFY = FALSE)
  
  merged <- Reduce(function(x, y) merge(x, y, by = "transcript_id"), dfs)
  
  merged[, 2:6] <- lapply(merged[, 2:6], function(x) as.numeric(gsub(",", ".", x)))
  
  filtered <- merged[
    apply(merged[, 2:6], 1, function(x) any(abs(x) > 0.2)),
  ]
  
  cor_matrix <- cor(filtered[, 2:6], method = "pearson", use = "complete.obs")
  
  write.csv(
    cor_matrix,
    file.path(data_dir, paste0("correlation_", type, ".csv")),
    row.names = TRUE
  )
}

analyze_log2fc()

# ============================================================
# 15) FUSE AND FILTER DIFFERENTIAL ISOFORM PROFILES
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

log2fc_threshold <- 0.8

diff_files <- list(
  p65   = file.path(diff_dir, paste0("diferencial_pareado_p65_", isoforms$p65, ".csv")),
  p66   = file.path(diff_dir, paste0("diferencial_pareado_p66_", isoforms$p66, ".csv")),
  p71   = file.path(diff_dir, paste0("diferencial_pareado_p71_", isoforms$p71, ".csv")),
  p72   = file.path(diff_dir, paste0("diferencial_pareado_p72_", isoforms$p72, ".csv")),
  numbl = file.path(diff_dir, paste0("diferencial_pareado_NUMBL_", isoforms$NUMBL, ".csv"))
)

diff_list <- lapply(diff_files, read.csv, stringsAsFactors = FALSE)

no_expression <- Reduce(`&`, lapply(diff_list, function(df) df$log2FC == 0 & is.na(df$p_value)))
diff_list <- lapply(diff_list, function(df) df[!no_expression, ])

not_significant <- Reduce(`&`, lapply(diff_list, function(df) df$p_value > 0.05 | is.na(df$p_value)))
diff_list <- lapply(diff_list, function(df) df[!not_significant, ])

low_log2fc <- Reduce(`&`, lapply(diff_list, function(df) abs(df$log2FC) < log2fc_threshold))
diff_list <- lapply(diff_list, function(df) df[!low_log2fc, ])

diff_p65   <- diff_list$p65
diff_p66   <- diff_list$p66
diff_p71   <- diff_list$p71
diff_p72   <- diff_list$p72
diff_numbl <- diff_list$numbl

diff_p65   <- diff_p65[, c("transcript_id", "log2FC", "p_value")]
diff_p66   <- diff_p66[, c("transcript_id", "log2FC", "p_value")]
diff_p71   <- diff_p71[, c("transcript_id", "log2FC", "p_value")]
diff_p72   <- diff_p72[, c("transcript_id", "log2FC", "p_value")]
diff_numbl <- diff_numbl[, c("transcript_id", "log2FC", "p_value")]

colnames(diff_p65)[-1]   <- c("log2FC_p65", "p_value_p65")
colnames(diff_p66)[-1]   <- c("log2FC_p66", "p_value_p66")
colnames(diff_p71)[-1]   <- c("log2FC_p71", "p_value_p71")
colnames(diff_p72)[-1]   <- c("log2FC_p72", "p_value_p72")
colnames(diff_numbl)[-1] <- c("log2FC_numbl", "p_value_numbl")

diff_fused <- diff_p65 %>%
  full_join(diff_p66,   by = "transcript_id") %>%
  full_join(diff_p71,   by = "transcript_id") %>%
  full_join(diff_p72,   by = "transcript_id") %>%
  full_join(diff_numbl, by = "transcript_id")

write.csv(
  diff_fused,
  file.path(data_dir, "differential_fused.csv"),
  row.names = FALSE
)

# ============================================================
# 16) ANNOTATE FUSED DIFFERENTIAL TABLE
# ============================================================

fused_diff <- read.csv(file.path(data_dir, "differential_fused.csv"), stringsAsFactors = FALSE)

gene_info <- read.csv(
  file.path("data", "metadata", "UcsC_gene_names.csv"),
  stringsAsFactors = FALSE
)

colnames(gene_info)[colnames(gene_info) == "ucsc"] <- "transcript_id"
gene_info$transcript_id <- sub("\\..*", "", gene_info$transcript_id)
gene_info <- gene_info[!duplicated(gene_info$transcript_id), ]

fused_annotated <- merge(
  fused_diff,
  gene_info[, c("transcript_id", "hgnc_symbol")],
  by = "transcript_id",
  all.x = TRUE
)

cols <- colnames(fused_annotated)
fused_annotated <- fused_annotated[, c("transcript_id", "hgnc_symbol", cols[!cols %in% c("transcript_id", "hgnc_symbol")])]
fused_annotated <- fused_annotated[!is.na(fused_annotated$hgnc_symbol), ]

write.csv(
  fused_annotated,
  file.path(data_dir, "differential_fused_NUMB_isoforms_annotated.csv"),
  row.names = FALSE
)

message("Annotated fused differential table saved.")

# ============================================================
# 17) HEATMAP ASSOCIATED WITH NUMB ISOFORMS (Figure 2D)
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(grid)
})

heat_input <- read.csv(
  file.path(data_dir, "differential_fused_NUMB_isoforms_annotated.csv"),
  row.names = 1
)

cols_interest <- c("log2FC_p72","log2FC_p71","log2FC_p66","log2FC_p65","log2FC_numbl")

heat_input[cols_interest] <- lapply(
  heat_input[cols_interest],
  function(x) as.numeric(gsub(",", ".", x))
)

threshold <- 1.5
filtered <- heat_input[
  apply(abs(heat_input[, cols_interest]) > threshold, 1, any),
]

mat <- t(filtered[, cols_interest])

out_file <- file.path(fig_2d, "Figure_2D_heatmap_log2FC_NUMB_isoforms_transposed.tiff")

tiff(out_file, width = 2000, height = 1500, res = 300, compression = "lzw")

ph <- pheatmap(
  mat,
  color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
  breaks = seq(-1.5, 1.5, length.out = 86),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 10,
  main = "Heatmap of log2FC (genes filtered by NUMB isoforms)",
  border_color = NA,
  silent = TRUE
)

grid.newpage()
grid.draw(ph$gtable)
dev.off()

message("Figure 2D heatmap saved to: ", out_file)
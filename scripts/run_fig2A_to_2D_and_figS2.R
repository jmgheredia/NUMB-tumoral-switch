# run_isoform_histograms_TCGA_BRCA.R
#
# Purpose:
# Generate expression and paired log2FC histograms for NUMB / NUMBL isoforms in TCGA-BRCA, including:
#   - Full cohort Normal and Tumor expression histograms (Supplementary)
#   - Paired Normal vs Tumor expression histograms (Figure 2A)
#   - Paired log2FC histograms (Figure 2B)
# Additionally:
#   - Dip test for multimodality (Normal/Tumor; full cohort)
#   - High/low split (k-means if multimodal; terciles otherwise)
#   - Differential tables and split tables (as in original script)
#
# Inputs:
# - data/processed/normal_samples.csv
# - data/processed/primary_tumor_samples.csv
#
# Outputs:
# - results/figures/supplementary/expression_histograms/
# - results/figures/figure2A/expression/
# - results/figures/figure2B/log2FC/
# - results/processed/normal_paired.csv
# - results/processed/tumor_paired.csv
# - results/processed/log2FC_paired_isoforms.csv
# - results/processed/diptest_analysis/differential/
# - results/processed/diptest_analysis/splits/

# -----------------------------
# 0) Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(diptest)
})

# -----------------------------
# 1) File paths
# -----------------------------
normal_file <- file.path("data", "processed", "normal_samples.csv")
tumor_file  <- file.path("data", "processed", "primary_tumor_samples.csv")

# -----------------------------
# 2) Sanity checks
# -----------------------------
if (!file.exists(normal_file)) stop("Missing file: ", normal_file)
if (!file.exists(tumor_file))  stop("Missing file: ", tumor_file)

# -----------------------------
# 3) Load data
# -----------------------------
normal_df <- read.csv(normal_file, stringsAsFactors = FALSE, check.names = FALSE)
tumor_df  <- read.csv(tumor_file,  stringsAsFactors = FALSE, check.names = FALSE)

colnames(normal_df)[1] <- "isoform_id"
colnames(tumor_df)[1]  <- "isoform_id"

# Remove transcript version suffix
normal_df$isoform_id <- sub("\\..*", "", normal_df$isoform_id)
tumor_df$isoform_id  <- sub("\\..*", "", tumor_df$isoform_id)

# -----------------------------
# 4) Expression columns
# -----------------------------
expr_cols_normal <- setdiff(colnames(normal_df), "isoform_id")
expr_cols_tumor  <- setdiff(colnames(tumor_df),  "isoform_id")

# -----------------------------
# 5) log2(x + 1) transformation
# -----------------------------
normal_df[, expr_cols_normal] <- log2(normal_df[, expr_cols_normal] + 1)
tumor_df[,  expr_cols_tumor]  <- log2(tumor_df[,  expr_cols_tumor]  + 1)

# -----------------------------
# 6) Isoforms of interest
# -----------------------------
isoforms <- list(
  p65   = "uc001xob",
  p66   = "uc001xoa",
  p71   = "uc001xnz",
  p72   = "uc001xny",
  NUMBL = "uc002oon"
)

# -----------------------------
# 7) Output directories
# -----------------------------
supp_dir <- file.path("results", "figures", "supplementary", "expression_histograms")
fig_expr <- file.path("results", "figures", "figure2A", "expression")
fig_fc   <- file.path("results", "figures", "figure2B", "log2FC")
data_dir <- file.path("results", "processed")

dip_dir    <- file.path(data_dir, "diptest_analysis")
diff_dir   <- file.path(dip_dir, "differential")
splits_dir <- file.path(dip_dir, "splits")

# Create directories (explicit hierarchy)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_expr, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_fc,   recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

dir.create(dip_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(diff_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(splits_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 7b) Plot styles
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

meta_normal <- data.frame(
  column  = expr_cols_normal,
  patient = extract_patient(expr_cols_normal),
  type    = extract_type(expr_cols_normal),
  stringsAsFactors = FALSE
)

meta_tumor <- data.frame(
  column  = expr_cols_tumor,
  patient = extract_patient(expr_cols_tumor),
  type    = extract_type(expr_cols_tumor),
  stringsAsFactors = FALSE
)

meta_normal <- subset(meta_normal, type == "11A")
meta_tumor  <- subset(meta_tumor,  type %in% c("01A", "01B"))

paired_patients <- intersect(meta_normal$patient, meta_tumor$patient)

meta_normal <- meta_normal[meta_normal$patient %in% paired_patients, ]
meta_tumor  <- meta_tumor[ meta_tumor$patient  %in% paired_patients, ]

meta_normal <- meta_normal[!duplicated(meta_normal$patient), ]
meta_tumor  <- meta_tumor[ !duplicated(meta_tumor$patient),  ]

meta_normal <- meta_normal[order(meta_normal$patient), ]
meta_tumor  <- meta_tumor[ order(meta_tumor$patient),  ]

# ============================================================
# 9) PAIRED EXPRESSION MATRICES
# ============================================================
normal_p <- normal_df[, c("isoform_id", meta_normal$column)]
tumor_p  <- tumor_df[,  c("isoform_id", meta_tumor$column)]

colnames(normal_p)[-1] <- meta_normal$patient
colnames(tumor_p)[-1]  <- meta_tumor$patient

stopifnot(all(colnames(normal_p)[-1] == colnames(tumor_p)[-1]))

write.csv(normal_p, file.path(data_dir, "normal_paired.csv"), row.names = FALSE)
write.csv(tumor_p,  file.path(data_dir, "tumor_paired.csv"),  row.names = FALSE)

# ============================================================
# 10) PAIRED EXPRESSION HISTOGRAMS (Figure 2A)
# ============================================================
for (iso in names(isoforms)) {
  
  tx <- isoforms[[iso]]
  
  row_n <- normal_p[normal_p$isoform_id == tx, ]
  row_t <- tumor_p[ tumor_p$isoform_id  == tx, ]
  
  if (nrow(row_n) == 0 || nrow(row_t) == 0) next
  
  v_n <- as.numeric(row_n[1, -1])
  v_t <- as.numeric(row_t[1, -1])
  
  out_n <- file.path(fig_expr, paste0("histogram_NORMAL_paired_", iso, "_", tx, ".tiff"))
  out_t <- file.path(fig_expr, paste0("histogram_TUMOR_paired_",  iso, "_", tx, ".tiff"))
  
  tiff(out_n, width = 1600, height = 1200, res = 200, compression = "lzw")
  hist(v_n, breaks = hist_breaks, prob = TRUE,
       main = paste("Normal (paired) -", iso),
       xlab = "log2(TPM + 1)",
       col = style_npair$bar, border = style_npair$border)
  lines(density(v_n, na.rm = TRUE), col = style_npair$line, lwd = 2)
  rug(v_n)
  dev.off()
  
  tiff(out_t, width = 1600, height = 1200, res = 200, compression = "lzw")
  hist(v_t, breaks = hist_breaks, prob = TRUE,
       main = paste("Tumor (paired) -", iso),
       xlab = "log2(TPM + 1)",
       col = style_tpair$bar, border = style_tpair$border)
  lines(density(v_t, na.rm = TRUE), col = style_tpair$line, lwd = 2)
  rug(v_t)
  dev.off()
}

# ============================================================
# 11) PAIRED log2FC HISTOGRAMS (Figure 2B)
# ============================================================
common_tx <- intersect(normal_p$isoform_id, tumor_p$isoform_id)

normal_p2 <- normal_p[normal_p$isoform_id %in% common_tx, ]
tumor_p2  <- tumor_p[ tumor_p$isoform_id  %in% common_tx, ]

normal_p2 <- normal_p2[order(normal_p2$isoform_id), ]
tumor_p2  <- tumor_p2[ order(tumor_p2$isoform_id),  ]

log2FC <- tumor_p2
log2FC[, -1] <- tumor_p2[, -1] - normal_p2[, -1]

hist_data <- data.frame()

for (iso in names(isoforms)) {
  
  tx <- isoforms[[iso]]
  row_tx <- log2FC[log2FC$isoform_id == tx, ]
  if (nrow(row_tx) == 0) next
  
  values <- as.numeric(row_tx[1, -1])
  
  hist_data <- rbind(
    hist_data,
    data.frame(
      isoform = iso,
      transcript_id = tx,
      patient = colnames(row_tx)[-1],
      log2FC = values,
      stringsAsFactors = FALSE
    )
  )
  
  out <- file.path(fig_fc, paste0("histogram_log2FC_paired_", iso, "_", tx, ".tiff"))
  tiff(out, width = 1600, height = 1200, res = 200, compression = "lzw")
  hist(values, breaks = hist_breaks_fc, prob = TRUE,
       main = paste("Paired log2FC -", iso),
       xlab = "log2FC (Tumor - Normal)",
       col = style_fc$bar, border = style_fc$border)
  lines(density(values, na.rm = TRUE), col = style_fc$line, lwd = 2)
  rug(values)
  dev.off()
}

write.csv(hist_data, file.path(data_dir, "log2FC_paired_isoforms_fig2B.csv"), row.names = FALSE)

# ============================================================
# 12) FULL COHORT TUMOR HISTOGRAMS (Supplementary)
# ============================================================
for (iso in names(isoforms)) {
  
  tx <- isoforms[[iso]]
  row_t <- tumor_df[tumor_df$isoform_id == tx, ]
  
  if (nrow(row_t) > 0) {
    
    values <- as.numeric(row_t[1, expr_cols_tumor])
    out <- file.path(supp_dir, paste0("histogram_TUMOR_full_", iso, "_", tx, ".tiff"))
    
    tiff(out, width = 1600, height = 1200, res = 200, compression = "lzw")
    hist(values, breaks = hist_breaks, prob = TRUE,
         main = paste("Tumor (full cohort) -", iso),
         xlab = "log2(TPM + 1)",
         col = style_full$bar, border = style_full$border)
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
  
  fila_tx <- log2FC[log2FC$isoform_id == tx, ]
  
  if (nrow(fila_tx) == 0) next
  
  valores_log2FC <- as.numeric(fila_tx[1, -1])
  pacientes <- colnames(fila_tx)[-1]
  names(valores_log2FC) <- pacientes
  
  message("Diptest for ", iso, " (", tx, ")")
  
  p_dip <- dip.test(valores_log2FC)$p.value
  
  if (p_dip < 0.05) {
    clustering <- kmeans(valores_log2FC, centers = 2)
    grupo_alto <- names(valores_log2FC)[clustering$cluster == which.max(clustering$centers)]
    grupo_bajo <- names(valores_log2FC)[clustering$cluster == which.min(clustering$centers)]
  } else {
    terciles <- quantile(valores_log2FC, probs = c(1/3, 2/3), na.rm = TRUE)
    grupo_bajo <- names(valores_log2FC[valores_log2FC <= terciles[1]])
    grupo_alto <- names(valores_log2FC[valores_log2FC >= terciles[2]])
  }
  
  message("Patients high: ", length(grupo_alto), " | low: ", length(grupo_bajo))
  
  if (length(grupo_alto) < 2 || length(grupo_bajo) < 2) {
    message("Skipping ", iso, " due to insufficient groups")
    next
  }
  
  genes <- log2FC$isoform_id
  log2FC_genes <- numeric(length(genes))
  p_vals <- numeric(length(genes))
  
  for (i in seq_along(genes)) {
    
    expr_alto <- as.numeric(log2FC[i, grupo_alto])
    expr_bajo <- as.numeric(log2FC[i, grupo_bajo])
    
    tt <- try(t.test(expr_alto, expr_bajo), silent = TRUE)
    
    if (inherits(tt, "htest")) {
      log2FC_genes[i] <- mean(expr_alto, na.rm = TRUE) - mean(expr_bajo, na.rm = TRUE)
      p_vals[i] <- tt$p.value
    } else {
      log2FC_genes[i] <- NA
      p_vals[i] <- NA
    }
  }
  
  resultados <- data.frame(
    transcript_id = genes,
    log2FC = log2FC_genes,
    p_value = p_vals,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    resultados,
    file.path(diff_dir, paste0("diferencial_pareado_", iso, "_", tx, ".csv")),
    row.names = FALSE
  )
  
  split_df <- data.frame(
    paciente = c(grupo_alto, grupo_bajo),
    grupo = c(rep("alto", length(grupo_alto)), rep("bajo", length(grupo_bajo))),
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
# 14) TRANSCRIPTOMIC PROFILE CORRELATION ASSOCIATED WITH
#     NUMB ISOFORMS (TCGA-BRCA PAIRED SAMPLES)
#     Single output: correlacion_PAREADO.csv
# ============================================================

# Function to load log2FC per isoform
cargar_log2fc <- function(nombre_archivo, nueva_col) {
  df <- read.csv(nombre_archivo, header = TRUE, stringsAsFactors = FALSE)
  df <- df[, c("transcript_id", "log2FC")]
  colnames(df)[2] <- nueva_col
  return(df)
}

# NUMB / NUMBL isoforms (same order as original script)
isoformas_corr <- c("p72", "p71", "p65", "p66", "numbl")
codigos_corr   <- c("uc001xny", "uc001xnz", "uc001xob", "uc001xoa", "uc002oon")

analyze_log2fc <- function(tipo = "PAREADO") {
  
  archivos <- file.path(
    diff_dir,
    paste0("diferencial_pareado_", isoformas_corr, "_", codigos_corr, ".csv")
  )
  
  nombres_col <- paste0("log2FC_", isoformas_corr)
  
  # Load per-isoform differential results
  listas <- mapply(cargar_log2fc, archivos, nombres_col, SIMPLIFY = FALSE)
  
  # Merge by transcript_id
  comparativa <- Reduce(function(x, y) merge(x, y, by = "transcript_id"), listas)
  
  # Explicit numeric conversion (robust to previous formats)
  comparativa[, 2:6] <- lapply(comparativa[, 2:6], function(x) {
    as.numeric(gsub(",", ".", x))
  })
  
  # Filter genes with signal in at least one isoform
  comparativa_filtrada <- comparativa[
    apply(comparativa[, 2:6], 1, function(x) any(abs(x) > 0.2)),
  ]
  
  # Pearson correlation between transcriptomic profiles
  cor_matrix <- cor(
    comparativa_filtrada[, 2:6],
    method = "pearson",
    use = "complete.obs"
  )
  
  # Single output
  write.csv(
    cor_matrix,
    file.path(data_dir, paste0("correlacion_", tipo, ".csv")),
    row.names = TRUE
  )
}

# Run
analyze_log2fc()

# ============================================================
# 15) FUSE AND FILTER DIFFERENTIAL ISOFORM PROFILES
#     Output: differential_fused.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# log2FC threshold (same as original script)
log2fc_threshold <- 0.8

# Differential files generated in Step 13
diff_files <- list(
  p65   = file.path(diff_dir, paste0("diferencial_pareado_p65_", isoforms$p65, ".csv")),
  p66   = file.path(diff_dir, paste0("diferencial_pareado_p66_", isoforms$p66, ".csv")),
  p71   = file.path(diff_dir, paste0("diferencial_pareado_p71_", isoforms$p71, ".csv")),
  p72   = file.path(diff_dir, paste0("diferencial_pareado_p72_", isoforms$p72, ".csv")),
  numbl = file.path(diff_dir, paste0("diferencial_pareado_NUMBL_", isoforms$NUMBL, ".csv"))
)

# Load differential datasets
diff_list <- lapply(diff_files, read.csv, stringsAsFactors = FALSE)

# Remove rows with no expression across all isoforms
no_expression <- Reduce(`&`, lapply(diff_list, function(df) df$log2FC == 0 & is.na(df$p_value)))
message("Rows with no expression: ", sum(no_expression))
diff_list <- lapply(diff_list, function(df) df[!no_expression, ])

# Remove rows with no statistical significance across all isoforms
not_significant <- Reduce(`&`, lapply(diff_list, function(df) df$p_value > 0.05 | is.na(df$p_value)))
message("Rows with no statistical significance: ", sum(not_significant))
diff_list <- lapply(diff_list, function(df) df[!not_significant, ])

# Remove rows with low log2FC across all isoforms
low_log2fc <- Reduce(`&`, lapply(diff_list, function(df) abs(df$log2FC) < log2fc_threshold))
message("Rows with |log2FC| < ", log2fc_threshold, " across all isoforms: ", sum(low_log2fc))
diff_list <- lapply(diff_list, function(df) df[!low_log2fc, ])

# Assign individual datasets
diff_p65   <- diff_list$p65
diff_p66   <- diff_list$p66
diff_p71   <- diff_list$p71
diff_p72   <- diff_list$p72
diff_numbl <- diff_list$numbl

message("Final rows after filtering: ", nrow(diff_p65))

# Keep required columns only
diff_p65   <- diff_p65[, c("transcript_id", "log2FC", "p_value")]
diff_p66   <- diff_p66[, c("transcript_id", "log2FC", "p_value")]
diff_p71   <- diff_p71[, c("transcript_id", "log2FC", "p_value")]
diff_p72   <- diff_p72[, c("transcript_id", "log2FC", "p_value")]
diff_numbl <- diff_numbl[, c("transcript_id", "log2FC", "p_value")]

# Rename columns to identify isoforms
colnames(diff_p65)[-1]   <- c("log2FC_p65", "p_value_p65")
colnames(diff_p66)[-1]   <- c("log2FC_p66", "p_value_p66")
colnames(diff_p71)[-1]   <- c("log2FC_p71", "p_value_p71")
colnames(diff_p72)[-1]   <- c("log2FC_p72", "p_value_p72")
colnames(diff_numbl)[-1] <- c("log2FC_numbl", "p_value_numbl")

# Sequential full join by transcript_id
diff_fused <- diff_p65 %>%
  full_join(diff_p66,   by = "transcript_id") %>%
  full_join(diff_p71,   by = "transcript_id") %>%
  full_join(diff_p72,   by = "transcript_id") %>%
  full_join(diff_numbl, by = "transcript_id")

# Save fused differential table
write.csv(
  diff_fused,
  file.path(data_dir, "differential_fused.csv"),
  row.names = FALSE
)

message("Differential profiles fused and saved to: ", file.path(data_dir, "differential_fused.csv"))

# ============================================================
# 16) ANNOTATE FUSED NUMB ISOFORM DIFFERENTIALS WITH GENE SYMBOLS
#     Output:
#       - differential_fused_NUMB_isoforms_annotated.csv
# ============================================================

# Load fused differential table
fused_diff <- read.csv(file.path(data_dir, "differential_fused.csv"), stringsAsFactors = FALSE)

# Load gene annotation table
gene_info <- read.csv(
  file.path("data", "metadata", "UcsC_gene_names.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

# Rename ucsc column to transcript_id
colnames(gene_info)[colnames(gene_info) == "ucsc"] <- "transcript_id"

# Remove transcript version suffix (.1, .2, etc.)
gene_info$transcript_id <- sub("\\..*", "", gene_info$transcript_id)

# Remove duplicated transcript_id (keep first)
gene_info <- gene_info[!duplicated(gene_info$transcript_id), ]

# Merge annotation into fused differentials
fused_annotated <- merge(
  fused_diff,
  gene_info[, c("transcript_id", "hgnc_symbol")],
  by = "transcript_id",
  all.x = TRUE
)

# Reorder columns: transcript_id | hgnc_symbol | rest
cols <- colnames(fused_annotated)
new_order <- c(
  "transcript_id",
  "hgnc_symbol",
  cols[!cols %in% c("transcript_id", "hgnc_symbol")]
)

fused_annotated <- fused_annotated[, new_order]

# Remove transcripts without gene symbol
fused_annotated <- fused_annotated[!is.na(fused_annotated$hgnc_symbol), ]

# Save final annotated table
write.csv(
  fused_annotated,
  file.path(data_dir, "differential_fused_NUMB_isoforms_annotated.csv"),
  row.names = FALSE
)

message("Annotated fused NUMB isoform differential table saved to: ",
        file.path(data_dir, "differential_fused_NUMB_isoforms_annotated.csv"))


# ============================================================
# 17) TRANSCRIPTOMIC HEATMAP ASSOCIATED WITH NUMB ISOFORMS
#     (Figure 2D)
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(grid)
})

# Read input file (pipeline output)
datos <- read.csv(
  file.path("results", "processed", "differential_fused_NUMB_isoforms_annotated.csv"),
  row.names = 1,
  sep = ",",
  header = TRUE
)

# Columns of interest (unchanged)
cols_interes <- c(
  "log2FC_p72",
  "log2FC_p71",
  "log2FC_p66",
  "log2FC_p65",
  "log2FC_numbl"
)

# Convert to numeric (replace decimal commas if present)
datos[cols_interes] <- lapply(
  datos[cols_interes],
  function(x) as.numeric(gsub(",", ".", x))
)

# Filtering: genes with |log2FC| > 1.5 in at least one isoform
umbral <- 1.5
datos_filtrados <- datos[
  apply(abs(datos[, cols_interes]) > umbral, 1, any),
]

# Original matrix (genes x isoforms)
mat <- datos_filtrados[, cols_interes]

# Transpose: isoforms x genes
mat_t <- t(mat)

# Output file
out_file <- file.path("results", "figures", "figure2D", "Figure_2D_heatmap_log2FC_NUMB_isoforms_transposed.tiff")
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

# Write heatmap to TIFF using base R device (robust on Windows)
tiff(
  out_file,
  width = 2000,
  height = 1500,
  res = 300,
  compression = "lzw"
)

ph <- pheatmap(
  mat_t,
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

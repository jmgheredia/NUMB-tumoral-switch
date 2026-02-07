# run_Figure4A_isoform_expression.R
#
# Generate Figure 4A: Isoform expression across passages with Welch t-tests
# versus Resection baseline.
#
# Inputs:
# - data/processed/pdmr_isoform_expression_log2.csv
# - data/metadata/metadata_PDMR.csv
#
# Outputs:
# - results/Figures/Figure_4/Figure4A_isoform_expression.tiff
#
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(rstatix)
})

# ============================================================
# 0) Input / Output
# ============================================================

input_expression <- file.path("data", "processed", "pdmr_isoform_expression_log2.csv")
metadata_file    <- file.path("data", "metadata", "metadata_PDMR.csv")

figure_dir <- file.path("results", "Figures", "Figure_4")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_expression))
stopifnot(file.exists(metadata_file))

# ============================================================
# Load data
# ============================================================

expression <- read.csv(input_expression, header = TRUE, check.names = FALSE)
metadata   <- read.csv(metadata_file, sep = ";", check.names = FALSE)

# ============================================================
# Helpers
# ============================================================

normalize_id <- function(x) {
  x |>
    as.character() |>
    gsub("^X", "", x = _) |>
    gsub("_", ".", x = _) |>
    gsub("-", ".", x = _) |>
    gsub("\\.+", ".", x = _)
}

# ============================================================
# Preprocessing
# ============================================================

expression$transcript_id <- sub("\\..*", "", expression$transcript_id)

isoforms_numb <- list(
  p65   = "uc001xob",
  p66   = "uc001xoa",
  p71   = "uc001xnz",
  p72   = "uc001xny",
  NUMBL = "uc002oon"
)

expression_interest <- expression %>%
  filter(transcript_id %in% unlist(isoforms_numb))

expression_interest$isoform <- names(isoforms_numb)[
  match(expression_interest$transcript_id, unlist(isoforms_numb))
]

if ("gene_id" %in% colnames(expression_interest)) {
  expression_interest <- expression_interest %>% select(-gene_id)
}

# ============================================================
# Reshape to long format
# ============================================================

expression_long <- expression_interest %>%
  pivot_longer(
    cols = -c(transcript_id, isoform),
    names_to = "sample",
    values_to = "expression"
  )

expression_long$sample <- normalize_id(expression_long$sample)

# ============================================================
# Join with metadata
# ============================================================

stopifnot("PatientID.SpecimenID.SampleID" %in% colnames(metadata))
metadata$sample <- normalize_id(metadata$PatientID.SpecimenID.SampleID)

idx_check <- match(expression_long$sample, metadata$sample)
if (any(is.na(idx_check))) {
  missing_samples <- unique(expression_long$sample[is.na(idx_check)])
  msg <- paste0(
    "Hay muestras sin match (", length(missing_samples), "). Ejemplos: ",
    paste(head(missing_samples, 10), collapse = ", ")
  )
  stop(msg)
}

full_data <- left_join(expression_long, metadata, by = "sample")

passages_of_interest <- c("Resection", "0", "1", "2", "3", "4")
full_data <- full_data %>%
  filter(Passage_of_this_sample %in% passages_of_interest)

full_data$Passage_of_this_sample <- factor(
  full_data$Passage_of_this_sample,
  levels = passages_of_interest
)

if (!("Resection" %in% levels(full_data$Passage_of_this_sample))) {
  stop("El nivel 'Resection' no existe en Passage.of.this.sample (levels).")
}
if (nrow(full_data %>% filter(Passage_of_this_sample == "Resection")) == 0) {
  stop("No hay muestras con Passage_of_this_sample == 'Resection' tras el filtrado.")
}

# ============================================================
# Normalize each isoform to its Resection mean
# ============================================================

resection_means <- full_data %>%
  filter(Passage_of_this_sample == "Resection") %>%
  group_by(isoform) %>%
  summarise(resection_mean = mean(expression, na.rm = TRUE), .groups = "drop")

bad_iso <- resection_means %>%
  filter(is.na(resection_mean) | resection_mean == 0) %>%
  pull(isoform)

if (length(bad_iso) > 0) {
  stop(paste0(
    "Resection mean es NA o 0 para isoform(s): ",
    paste(bad_iso, collapse = ", "),
    ". No se puede normalizar."
  ))
}

full_data <- full_data %>%
  left_join(resection_means, by = "isoform") %>%
  mutate(expression = expression / resection_mean) %>%
  select(-resection_mean)

# ============================================================
# Summary statistics
# ============================================================

summary_df <- full_data %>%
  group_by(isoform, Passage_of_this_sample) %>%
  summarise(
    mean = mean(expression, na.rm = TRUE),
    sd   = sd(expression, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n),
    .groups = "drop"
  )

summary_df$Passage_of_this_sample <- factor(
  summary_df$Passage_of_this_sample,
  levels = passages_of_interest
)

# ============================================================
# Plot
# ============================================================

p <- ggplot(summary_df, aes(x = isoform, y = mean, fill = Passage_of_this_sample)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  labs(
    title = "Isoform expression across passages",
    x = "Isoform",
    y = "Relative expression (normalized to Resection mean)",
    fill = "Passage"
  ) +
  scale_fill_manual(values = c(
    "Resection" = "#E41A1C",
    "0"         = "#FFD92F",
    "1"         = "#4DAF4A",
    "2"         = "#377EB8",
    "3"         = "#984EA3",
    "4"         = "#999999"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

figure_file <- file.path(figure_dir, "Figure4A_isoform_expression.tiff")

ggsave(
  filename = figure_file,
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  compression = "lzw"
)

message("Figure saved to: ", figure_file)

# ============================================================
# Statistical testing: Welch t-test vs Resection (per isoform)
# ============================================================

summary_stats <- full_data %>%
  group_by(isoform, Passage_of_this_sample) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    sd_expression   = sd(expression, na.rm = TRUE),
    n_samples       = n(),
    sem_expression  = sd_expression / sqrt(n_samples),
    .groups = "drop"
  )

stats_vs_resection <- full_data %>%
  group_by(isoform) %>%
  t_test(
    expression ~ Passage_of_this_sample,
    ref.group = "Resection",
    detailed  = TRUE,
    var.equal = FALSE
  ) %>%
  ungroup() %>%
  mutate(significance = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )) %>%
  dplyr::rename(Passage_of_this_sample = group2) %>%
  dplyr::select(isoform, Passage_of_this_sample, p, significance, statistic, df)

summary_with_stats <- summary_stats %>%
  left_join(stats_vs_resection, by = c("isoform", "Passage_of_this_sample"))

message("Statistical analysis completed (summary_with_stats object created).")

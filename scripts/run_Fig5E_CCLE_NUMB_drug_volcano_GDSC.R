# ============================================================
# run_Fig5E_CCLE_NUMB_drug_volcano_GDSC.R
#
# Generates volcano plots of differential drug sensitivity between High vs Low NUMB-score
# BRCA cell lines using GDSC1 and GDSC2 data (Figure 5E).
#
# Input:
#   - data/processed/GDSC1_breast.csv
#   - data/processed/GDSC2_breast.csv
#   - results/processed/Figure_5/BCCL_grouped_by_NUMB_score.csv
#
# Output:
#   Figures:
#     - results/Figures/Figure_5/Figure5E_volcano_NUMB_score.tiff
#     - results/Figures/Figure_5/Figure5E_volcano_NUMB_score_all_labels.tiff
#     - results/Figures/Figure_5/Figure5E_volcano_NUMB_score_no_labels.tiff
#
#   Processed data:
#     - results/processed/Figure_5/Drugs_result_by_scoreNUMB.csv
#     - results/processed/Figure_5/Drugs_significant_NUMB_volcano.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

# ============================================================
# 0) Input / Output
# ============================================================

gdsc1_csv <- file.path("data", "processed", "GDSC1_breast.csv")
gdsc2_csv <- file.path("data", "processed", "GDSC2_breast.csv")
group_csv <- file.path("results", "processed", "Figure_5", "BCCL_grouped_by_NUMB_score.csv")

fig_dir  <- file.path("results", "Figures", "Figure_5")
proc_dir <- file.path("results", "processed", "Figure_5")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(gdsc1_csv)) stop("Missing GDSC1_breast.csv in data/processed/")
if (!file.exists(gdsc2_csv)) stop("Missing GDSC2_breast.csv in data/processed/")
if (!file.exists(group_csv)) stop("Missing BCCL_grouped_by_NUMB_score.csv in results/processed/Figure_5/")

# ============================================================
# 2) Load data
# ============================================================

GDSC1 <- read.csv(gdsc1_csv, stringsAsFactors = FALSE)
GDSC2 <- read.csv(gdsc2_csv, stringsAsFactors = FALSE)
BCCL_score_NUMB <- read.csv(group_csv, stringsAsFactors = FALSE)

# ============================================================
# 2b) Normalize cell line names (GDSC <-> NUMB)
# ============================================================

normalize_name <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

GDSC1 <- GDSC1 %>%
  mutate(CELL_LINE_NAME_NORM = normalize_name(CELL_LINE_NAME))

GDSC2 <- GDSC2 %>%
  mutate(CELL_LINE_NAME_NORM = normalize_name(CELL_LINE_NAME))

BCCL_score_NUMB <- BCCL_score_NUMB %>%
  mutate(cell_line_NORM = normalize_name(cell_line))

# ============================================================
# 3) Merge GDSC1 + GDSC2 and keep NUMB-extreme lines
# ============================================================

GDSC_merged <- bind_rows(GDSC1, GDSC2)

cell_lines_interest_norm <- BCCL_score_NUMB$cell_line_NORM

GDSC_merged <- GDSC_merged %>%
  filter(CELL_LINE_NAME_NORM %in% cell_lines_interest_norm)

# ============================================================
# 4) Add NUMB-score group
# ============================================================

GDSC_with_group <- GDSC_merged %>%
  left_join(
    BCCL_score_NUMB[, c("cell_line_NORM", "group")],
    by = c("CELL_LINE_NAME_NORM" = "cell_line_NORM")
  ) %>%
  dplyr::select(-CELL_LINE_NAME_NORM)

# ============================================================
# 5) Mean IC50/AUC per drug and group
# ============================================================

means_by_group <- GDSC_with_group %>%
  group_by(DRUG_NAME, group) %>%
  summarise(
    AUC = mean(AUC, na.rm = TRUE),
    LN_IC50 = mean(LN_IC50, na.rm = TRUE),
    .groups = "drop"
  )

means_wide <- means_by_group %>%
  pivot_wider(
    names_from = group,
    values_from = c(AUC, LN_IC50),
    names_sep = "_"
  )

drug_freq <- GDSC_with_group %>%
  group_by(DRUG_NAME) %>%
  summarise(Frequency = n(), .groups = "drop")

Drugs <- drug_freq %>%
  left_join(means_wide, by = "DRUG_NAME")

# ============================================================
# 6) Statistics (t-test per drug)
# ============================================================

stats <- GDSC_with_group %>%
  filter(!is.na(group)) %>%
  group_by(DRUG_NAME) %>%
  summarise(
    pval_LN_IC50 = tryCatch(t.test(LN_IC50 ~ group)$p.value, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(pval_LN_IC50_adj = p.adjust(pval_LN_IC50, method = "fdr"))

Drugs <- Drugs %>%
  left_join(stats, by = "DRUG_NAME")

write.csv(
  Drugs,
  file.path(proc_dir, "Drugs_result_by_scoreNUMB.csv"),
  row.names = FALSE
)

# ============================================================
# 7) Prepare volcano plot data
# ============================================================

plot_data <- Drugs %>%
  dplyr::select(DRUG_NAME, LN_IC50_high, LN_IC50_low, pval_LN_IC50) %>%
  filter(!is.na(pval_LN_IC50)) %>%
  mutate(
    delta_IC50 = LN_IC50_low - LN_IC50_high,
    neg_log10_pval = -log10(pval_LN_IC50),
    is_significant = pval_LN_IC50 < 0.05,
    label = ifelse(is_significant & delta_IC50 > 0, DRUG_NAME, NA)
  )

significant_drugs <- plot_data %>%
  filter(is_significant) %>%
  arrange(desc(neg_log10_pval)) %>%
  dplyr::select(
    DRUG_NAME,
    LN_IC50_high,
    LN_IC50_low,
    delta_IC50,
    neg_log10_pval,
    pval_LN_IC50
  )

write.csv(
  significant_drugs,
  file.path(proc_dir, "Drugs_significant_NUMB_volcano.csv"),
  row.names = FALSE
)

# ============================================================
# 8) Volcano plotting function
# ============================================================

base_volcano <- function(data, title_text, show_labels = TRUE, show_legend = TRUE) {
  p <- ggplot(data, aes(x = delta_IC50, y = neg_log10_pval)) +
    geom_point(aes(color = is_significant)) +
    scale_color_manual(values = c("gray70", "red")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    labs(
      title = title_text,
      x = expression(Delta * " LN(IC50) [Low – High NUMB score]"),
      y = expression(-log[10] * "(p-value)"),
      color = if (show_legend) "Significant (p < 0.05)" else NULL
    ) +
    coord_cartesian(ylim = c(0, 4)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "right" else "none"
    )
  
  if (show_labels) {
    p <- p + ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf)
  }
  
  p
}

plot_all_labels <- plot_data %>% mutate(label = ifelse(is_significant, DRUG_NAME, NA))
plot_no_labels  <- plot_data %>% mutate(label = NA)

# ============================================================
# 9) Save figures (Figure 5E)
# ============================================================

tiff(file.path(fig_dir, "Figure5E_volcano_NUMB_score.tiff"),
     width = 8, height = 6, units = "in", res = 300, compression = "lzw")
print(base_volcano(plot_data, "Volcano plot – Differential drug sensitivity"))
dev.off()

tiff(file.path(fig_dir, "Figure5E_volcano_NUMB_score_all_labels.tiff"),
     width = 10, height = 8, units = "in", res = 300, compression = "lzw")
print(base_volcano(plot_all_labels, "Volcano plot – All drugs labelled"))
dev.off()

tiff(file.path(fig_dir, "Figure5E_volcano_NUMB_score_no_labels.tiff"),
     width = 8, height = 6, units = "in", res = 300, compression = "lzw")
print(base_volcano(plot_no_labels, "Volcano plot – No labels", show_labels = FALSE, show_legend = FALSE))
dev.off()

message("Figure 5E GDSC volcano plots completed.")

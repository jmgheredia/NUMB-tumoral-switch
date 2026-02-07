# run_fig1A_switch_summary.R
# Generates Figure 1A: summary of significant isoform switches per cancer type
# Output: results/figures/figure1A_switch_summary.tiff

# -----------------------------
# 0) Input / Output
# -----------------------------
# This repository does not redistribute the third-party Rdata file.
# Download and place it at: data/external/File_S1_switchAnalyzeRlists.Rdata

input_data <- file.path("data", "external", "File_S1_switchAnalyzeRlists.Rdata")
output_dir <- file.path("results", "figures", "Figure_1")

# -----------------------------
# 1) Libraries
# -----------------------------
library(IsoformSwitchAnalyzeR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)  # for unit()

# -----------------------------
# 2) Load input data
# -----------------------------
if (!file.exists(input_data)) {
  stop(
    paste0(
      "Missing required input file: ", input_data, "\n",
      "Please download it from Figshare (DOI: 10.6084/m9.figshare.4924724)\n",
      "and place it in data/external/."
    )
  )
}

load(input_data)

# -----------------------------
# 3) Compute summary per cancer type
# -----------------------------
cancer_types <- gsub(
  "_cancer",
  "",
  unique(tcgaSwitchAnalyzeRlistSignificant$isoformFeatures$condition_2)
)

results_sig <- data.frame()

for (cancer in cancer_types) {
  cancer_condition <- paste0(cancer, "_cancer")
  control_condition <- paste0(cancer, "_ctrl")
  
  # Subset: only genes with functional consequences
  subset_obj <- subsetSwitchAnalyzeRlist(
    tcgaSwitchAnalyzeRlistSignificant,
    tcgaSwitchAnalyzeRlistSignificant$isoformFeatures$condition_2 == cancer_condition &
      tcgaSwitchAnalyzeRlistSignificant$isoformFeatures$condition_1 == control_condition &
      tcgaSwitchAnalyzeRlistSignificant$isoformFeatures$switchConsequencesGene == TRUE
  )
  
  # Extract isoforms with dIF
  iso_df <- subset_obj$isoformFeatures %>%
    filter(!is.na(dIF))
  
  # Counts
  n_genes <- length(unique(iso_df$gene_id))
  n_switches <- nrow(iso_df) / 2  # each switch = 2 isoforms
  n_isoforms <- length(unique(iso_df$isoform_id))
  
  results_sig <- rbind(
    results_sig,
    data.frame(
      Cancer = cancer,
      Genes = n_genes,
      Switches = n_switches,
      Isoforms = n_isoforms
    )
  )
}

results_sig_long <- pivot_longer(
  results_sig,
  cols = c("Genes", "Switches", "Isoforms"),
  names_to = "Feature",
  values_to = "Count"
)

# -----------------------------
# 4) Plot (final style)
# -----------------------------
p <- ggplot(results_sig_long, aes(x = Cancer, y = Count, fill = Feature)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.5) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_blank(),
    
    # Legend
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.key.size = unit(0.8, "cm"),
    legend.spacing.y = unit(0.4, "cm"),
    legend.text = element_text(size = 11)
  ) +
  labs(y = "Number", x = NULL) +
  scale_fill_manual(
    values = c(
      "Isoforms" = "#0072B2",
      "Switches" = "#E69F00",
      "Genes" = "#666666"
    )
  ) +
  guides(fill = guide_legend(byrow = TRUE))

# -----------------------------
# 5) Save output
# -----------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(output_dir, "figure1A_switch_summary.tiff"),
  plot = p,
  dpi = 600,
  width = 7, height = 5,
  units = "in",
  compression = "lzw",
  device = "tiff"
)

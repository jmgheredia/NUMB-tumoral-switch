# run_fig1B_1C_splicing_events.R
# Generates:
# - Figure 1B: Alternative splicing events across cancer types
# - Figure 1C: BRCA splicing gain vs loss (tumor vs control)
#
# Statistical analyses are computed to support interpretation and figure legends

# -----------------------------
# 0) Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(IsoformSwitchAnalyzeR)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(DescTools)
})

# -----------------------------
# 1) Input data
# -----------------------------
# This repository does not redistribute third-party Rdata files.
# Download and place the file at:
# data/external/File_S1_switchAnalyzeRlists.Rdata
input_rdata <- file.path(
  "data", "external", "File_S1_switchAnalyzeRlists.Rdata"
)

if (!file.exists(input_rdata)) {
  stop(
    paste0(
      "Missing required input file: ", input_rdata, "\n",
      "Please download it from Figshare (DOI: 10.6084/m9.figshare.4924724)\n",
      "and place it in data/external/."
    )
  )
}

load(input_rdata)

# -----------------------------
# 2) Alternative splicing analysis (pan-cancer)
# -----------------------------
tcgaSubset <- analyzeAlternativeSplicing(
  tcgaSwitchAnalyzeRlistSignificant,
  showProgress = FALSE
)

splicing_annotated <- tcgaSubset$AlternativeSplicingAnalysis %>%
  left_join(
    tcgaSubset$isoformFeatures[, c("isoform_id", "condition_2")],
    by = "isoform_id"
  )

splicing_counts <- splicing_annotated %>%
  mutate(cancer_type = sub("_.*", "", condition_2)) %>%
  pivot_longer(
    cols = c(ES, MEE, MES, IR, A5, A3, ATSS, ATTS),
    names_to = "event_type",
    values_to = "event_present"
  ) %>%
  filter(event_present == TRUE) %>%
  group_by(cancer_type, event_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    event_type = recode(
      event_type,
      A5   = "AD",
      A3   = "AA",
      ATSS = "AT",
      ATTS = "AP"
    )
  )

# -----------------------------
# 3) Figure 1B: splicing events across cancer types
# -----------------------------
splicing_colors <- c(
  ES = "#CC79A7",
  IR = "#0072B2",
  MES = "#56B4E9",
  AD = "#009E73",
  AA = "#E69F00",
  AT = "#D55E00",
  AP = "#8E6C8A"
)

p1B <- ggplot(
  splicing_counts,
  aes(x = factor(cancer_type), y = count, fill = event_type)
) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = splicing_colors) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1
    )
  ) +
  labs(
    x = NULL,
    y = "# ASEs",
    fill = "Event type"
  )

out_fig_dir <- file.path("results", "figures")
dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(
    out_fig_dir,
    "figure1B_splicing_events_across_cancers.tiff"
  ),
  plot = p1B,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  device = "tiff"
)

# -----------------------------
# 4) BRCA-specific analysis (tumor vs control)
# -----------------------------
brca_subset <- subsetSwitchAnalyzeRlist(
  tcgaSubset,
  tcgaSubset$isoformFeatures$condition_1 == "BRCA_ctrl" &
    tcgaSubset$isoformFeatures$condition_2 == "BRCA_cancer"
)

brca_subset <- analyzeAlternativeSplicing(
  brca_subset,
  showProgress = FALSE
)

brca_annotated <- brca_subset$AlternativeSplicingAnalysis %>%
  left_join(
    brca_subset$isoformFeatures[, c("isoform_id", "dIF")],
    by = "isoform_id"
  ) %>%
  mutate(
    condition = ifelse(dIF > 0, "BRCA_cancer", "BRCA_ctrl")
  )

brca_event_counts <- brca_annotated %>%
  pivot_longer(
    cols = c(A3, A5, ATSS, ATTS, ES, IR, MES),
    names_to = "event_type",
    values_to = "event_present"
  ) %>%
  filter(event_present == TRUE) %>%
  group_by(event_type, condition) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    direction = ifelse(condition == "BRCA_cancer", "+", "-")
  )

# -----------------------------
# 5) Figure 1C: BRCA splicing gain vs loss
# -----------------------------
p1C <- ggplot(
  brca_event_counts,
  aes(x = direction, y = count, fill = condition)
) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.6),
    width = 0.5
  ) +
  facet_wrap(~event_type, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c(
      "BRCA_ctrl"   = "#ff7f00",
      "BRCA_cancer" = "#984ea3"
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_blank()
  ) +
  labs(y = "Number", fill = NULL)

ggsave(
  filename = file.path(
    out_fig_dir,
    "figure1C_BRCA_splicing_gain_loss.tiff"
  ),
  plot = p1C,
  width = 7,
  height = 2.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  device = "tiff"
)

# -----------------------------
# 6) Statistical analyses
# -----------------------------

# 6.1 Global compositional test (2 x K)
tabC <- xtabs(count ~ condition + event_type, data = brca_event_counts)

chisq_global <- suppressWarnings(chisq.test(tabC))
gtest_global  <- GTest(tabC)

# 6.2 Event-level 2x2 tests (proportions)
totals_by_cond <- brca_event_counts %>%
  group_by(condition) %>%
  summarise(total = sum(count), .groups = "drop")

event_2x2_stats <- brca_event_counts %>%
  left_join(totals_by_cond, by = "condition") %>%
  mutate(non_event = total - count) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(count, non_event),
    values_fill = 0
  ) %>%
  rowwise() %>%
  mutate(
    a = count_BRCA_cancer,
    b = non_event_BRCA_cancer,
    c = count_BRCA_ctrl,
    d = non_event_BRCA_ctrl,
    p_value = if (min(c(a, b, c, d)) < 5) {
      fisher.test(
        matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      )$p.value
    } else {
      chisq.test(
        matrix(c(a, b, c, d), nrow = 2, byrow = TRUE),
        correct = FALSE
      )$p.value
    }
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# 6.3 Binomial test (directional bias)
binomial_stats <- brca_event_counts %>%
  select(event_type, condition, count) %>%
  pivot_wider(
    names_from = condition,
    values_from = count,
    values_fill = 0
  ) %>%
  rowwise() %>%                      # ← CLAVE
  mutate(
    total = BRCA_cancer + BRCA_ctrl,
    p_value = binom.test(
      x = BRCA_cancer,
      n = total,
      p = 0.5,
      alternative = "two.sided"
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

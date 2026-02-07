# ============================================================
# run_Fig6E_FigS7C_NUMB_score_TCGA_survival.R
#
# Performs Kaplan–Meier survival analysis of TCGA-BRCA patients
# based on NUMB score groups:
# - WT analysis using a fixed NUMB-score cutoff
# - Sensitivity analysis excluding top 1% influential samples
#   based on DFBETA from a Cox proportional hazards model
#
# Input:
#   - data/metadata/clinical.tsv
#   - results/processed/Figure_6/expr_mat_total_NUMB_score_TCGA.csv
#   - results/processed/Figure_6/score_NUMB_combined_PAM50_and_NonTumor.csv
#
# Output:
#   Figures:
#     - results/Figures/Supplementary/S7/Figure_S7C_KM_NUMBscore_255_WT.tiff
#     - results/Figures/Figure_6/Figure_6E_KM_NUMBscore_255_DFBETA_top1pct_removed.tiff
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# ============================================================
# 1) Input files
# ============================================================

clinical_tsv <- file.path("data", "metadata", "clinical.tsv")

score_csv <- file.path(
  "results", "processed", "Figure_6",
  "expr_mat_total_NUMB_score_TCGA.csv"
)

pam50_csv <- file.path(
  "results", "processed", "Figure_6",
  "score_NUMB_combined_PAM50_and_NonTumor.csv"
)

if (!file.exists(clinical_tsv))
  stop("Missing clinical.tsv in data/metadata/")

if (!file.exists(score_csv))
  stop("Missing expr_mat_total_NUMB_score_TCGA.csv in results/processed/Figure_6/")

if (!file.exists(pam50_csv))
  stop("Missing score_NUMB_combined_PAM50_and_NonTumor.csv in results/processed/Figure_6/")

# ============================================================
# 2) Output directories
# ============================================================

fig6_dir <- file.path("results", "Figures", "Figure_6")
figS7_dir <- file.path("results", "Figures", "Supplementary", "S7")

dir.create(fig6_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figS7_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3) Load data
# ============================================================

clinical <- read_tsv(clinical_tsv, col_types = cols())
numb_scores <- read.csv(score_csv, stringsAsFactors = FALSE)
pam50_groups <- read.csv(pam50_csv, stringsAsFactors = FALSE)

# ============================================================
# 4) Keep tumor samples and merge PAM50
# ============================================================

numb_scores_tumor <- subset(
  numb_scores,
  grepl("Tumor", Type, ignore.case = TRUE)
)

pam50_tumor <- pam50_groups %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  dplyr::select(patient_id, Group_PAM50_Combined)

numb_scores_tumor <- numb_scores_tumor %>%
  left_join(pam50_tumor, by = "patient_id")

# ============================================================
# 5) Merge with clinical data
# ============================================================

clinical_con_score <- clinical %>%
  inner_join(
    numb_scores_tumor,
    by = c("cases.submitter_id" = "patient_id")
  ) %>%
  filter(!duplicated(cases.submitter_id))

# ============================================================
# 6) Define survival variables
# ============================================================

clinical_final <- clinical_con_score %>%
  mutate(
    days_to_death_num =
      as.numeric(as.character(demographic.days_to_death)),
    days_to_followup_num =
      as.numeric(as.character(diagnoses.days_to_last_follow_up)),
    days_to_followup_num =
      ifelse(
        is.na(days_to_followup_num) &
          demographic.vital_status == "Alive",
        4500,
        days_to_followup_num
      ),
    survival_time =
      coalesce(days_to_death_num, days_to_followup_num),
    event = demographic.vital_status == "Dead"
  ) %>%
  filter(!is.na(survival_time), !is.na(NUMB_score))

# ============================================================
# 7) Define NUMB-score groups
# ============================================================

threshold_score <- 2.55

clinical_final <- clinical_final %>%
  mutate(
    score_group = ifelse(
      NUMB_score < threshold_score,
      "Low score",
      "High score"
    )
  )

clinical_final$score_group <- factor(
  clinical_final$score_group,
  levels = c("Low score", "High score")
)

# ============================================================
# FIGURE S7C — WT Kaplan–Meier
# ============================================================

fit_WT <- survfit(
  Surv(survival_time, event) ~ score_group,
  data = clinical_final
)

logrank_WT <- survdiff(
  Surv(survival_time, event) ~ score_group,
  data = clinical_final
)

p_WT <- 1 - pchisq(
  logrank_WT$chisq,
  df = length(logrank_WT$n) - 1
)

cat("WT log-rank p-value:", p_WT, "\n")

tiff(
  file.path(figS7_dir, "Figure_S7C_KM_NUMBscore_255_WT.tiff"),
  width = 2000, height = 1800, res = 300
)

ggsurvplot(
  fit_WT,
  data = clinical_final,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  break.time.by = 1000,
  palette = c("blue", "red"),
  title = "WT Kaplan–Meier (NUMB-score cutoff = 2.55)",
  xlab = "Days",
  ylab = "Overall survival probability",
  xlim = c(0, 3000),
  ylim = c(0.6, 1),
  legend.title = "Group",
  legend.labs = c("Low score", "High score"),
  ggtheme = theme_minimal()
)

dev.off()

# ============================================================
# FIGURE 6E — DFBETA sensitivity (top 1% removed)
# ============================================================

cox_base <- coxph(
  Surv(survival_time, event) ~ score_group,
  data = clinical_final
)

dfb <- residuals(cox_base, type = "dfbeta")
dfb_vec <- if (is.matrix(dfb)) as.numeric(dfb[, 1]) else as.numeric(dfb)

clinical_infl <- clinical_final %>%
  mutate(dfbeta_abs = abs(dfb_vec))

cut_infl <- quantile(
  clinical_infl$dfbeta_abs,
  0.99,
  na.rm = TRUE
)

clinical_sens <- clinical_infl %>%
  filter(dfbeta_abs <= cut_infl)

fit_sens <- survfit(
  Surv(survival_time, event) ~ score_group,
  data = clinical_sens
)

logrank_sens <- survdiff(
  Surv(survival_time, event) ~ score_group,
  data = clinical_sens
)

p_sens <- 1 - pchisq(
  logrank_sens$chisq,
  df = length(logrank_sens$n) - 1
)

cat("DFBETA sensitivity log-rank p-value:", p_sens, "\n")

tiff(
  file.path(
    fig6_dir,
    "Figure_6E_KM_NUMBscore_255_DFBETA_top1pct_removed.tiff"
  ),
  width = 2000, height = 1800, res = 300
)

ggsurvplot(
  fit_sens,
  data = clinical_sens,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  break.time.by = 1000,
  palette = c("blue", "red"),
  title = "Sensitivity Kaplan–Meier (excluding top 1% influential cases by DFBETA)",
  xlab = "Days",
  ylab = "Overall survival probability",
  xlim = c(0, 3000),
  ylim = c(0.6, 1),
  legend.title = "Group",
  legend.labs = c("Low score", "High score"),
  ggtheme = theme_minimal()
)

dev.off()

message("Figure S7B (WT) and Figure 6E (DFBETA sensitivity) Kaplan–Meier analyses completed.")

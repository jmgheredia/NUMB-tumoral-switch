# ============================================================
# run_FigS7A_TCGA_NUMB_isoform_survival_SCAN.R
#
# TCGA-BRCA overall survival analysis by individual NUMB isoforms
# using RNA-seq expression data.
#
# Performs:
# - Optimal cutpoint scan per isoform (maxstat)
# - Minimal group size constraint (minprop = 15%)
# - Kaplan–Meier survival curves
# - Cox proportional hazards models
# - Benjamini–Hochberg correction across isoforms
# - Sensitivity analysis across different minprop thresholds
#
# This script generates Supplementary Figure S7A and
# associated Supplementary Tables.
#
# Input:
#   - data/metadata/clinical.tsv
#   - results/processed/Figure_6/expr_mat_total_NUMB_score_TCGA.csv
#
# Output:
#   Figures (Figure S7A):
#     - results/Figures/Supplementary/S7/Figure_S7A/KM_<isoform>_SCAN_min15pct.tiff
#
#   Tables:
#     - results/Tables/Tables_S7A/TableS_TCGA_Cox_by_isoform_SCAN_min15pct.csv
#     - results/Tables/Tables_S7A/TableS_TCGA_Cox_by_isoform_SCAN_sensitivity_minprop.csv
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

expr_csv <- file.path(
  "results", "processed", "Figure_6",
  "expr_mat_total_NUMB_score_TCGA.csv"
)

if (!file.exists(clinical_tsv))
  stop("Missing clinical.tsv in data/metadata/")

if (!file.exists(expr_csv))
  stop("Missing expr_mat_total_NUMB_score_TCGA.csv in results/processed/Figure_6/")

# ============================================================
# 2) Output directories
# ============================================================

figS7A_dir <- file.path(
  "results", "Figures", "Supplementary", "S7", "Figure_S7A"
)

tableS7A_dir <- file.path(
  "results", "Tables", "Tables_S7A"
)

dir.create(figS7A_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tableS7A_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper: robust cutpoint extraction across survminer versions
# ============================================================

extract_cutpoint_value <- function(cut_obj, varname) {
  
  cp <- cut_obj$cutpoint
  
  if (is.data.frame(cp)) {
    
    if (all(c("variable", "cutpoint") %in% colnames(cp))) {
      return(cp %>% filter(variable == varname) %>% pull(cutpoint) %>% as.numeric())
    }
    
    if (all(c("var", "cutpoint") %in% colnames(cp))) {
      return(cp %>% filter(var == varname) %>% pull(cutpoint) %>% as.numeric())
    }
    
    if (varname %in% colnames(cp)) {
      return(as.numeric(cp[[varname]][1]))
    }
    
    if ("cutpoint" %in% colnames(cp)) {
      return(as.numeric(cp$cutpoint[1]))
    }
    
    return(as.numeric(cp[1, 1]))
  }
  
  as.numeric(cp)
}

# ============================================================
# 3) Load data
# ============================================================

clinical <- read_tsv(clinical_tsv, col_types = cols())

expr_tbl <- read.csv(expr_csv, stringsAsFactors = FALSE)

# ============================================================
# 4) Keep tumour samples only
# ============================================================

expr_tumor <- expr_tbl %>%
  filter(grepl("Tumor", Type, ignore.case = TRUE)) %>%
  distinct(patient_id, .keep_all = TRUE)

# ============================================================
# 5) Minimal clinical cleaning
# ============================================================

clinical <- clinical %>%
  select(
    cases.submitter_id,
    demographic.vital_status,
    demographic.days_to_death,
    diagnoses.days_to_last_follow_up
  )

# ============================================================
# 6) Merge clinical + expression
# ============================================================

df <- clinical %>%
  inner_join(expr_tumor, by = c("cases.submitter_id" = "patient_id")) %>%
  filter(!duplicated(cases.submitter_id))

# ============================================================
# 7) Build survival variables (OS)
# ============================================================

df <- df %>%
  mutate(
    days_to_death_num = as.numeric(as.character(demographic.days_to_death)),
    days_to_followup_num = as.numeric(as.character(diagnoses.days_to_last_follow_up)),
    days_to_followup_num =
      ifelse(is.na(days_to_followup_num) & demographic.vital_status == "Alive",
             4500, days_to_followup_num),
    survival_time = coalesce(days_to_death_num, days_to_followup_num),
    event = demographic.vital_status == "Dead"
  ) %>%
  filter(!is.na(survival_time))

# ============================================================
# 8) Isoforms and parameters
# ============================================================

isoforms <- c("p66", "p65", "p71", "p72", "NUMBL")
minprop <- 0.15

# ============================================================
# 9) Loop: cutpoint scan + KM + Cox
# ============================================================

results <- list()

for (iso in isoforms) {
  
  if (!iso %in% colnames(df)) {
    warning("Isoform column not found: ", iso)
    next
  }
  
  message("[INFO] Processing isoform: ", iso)
  
  df_iso <- df %>%
    filter(!is.na(.data[[iso]])) %>%
    transmute(
      survival_time = survival_time,
      event = as.integer(event),
      value = .data[[iso]]
    )
  
  df_iso2 <- df_iso %>%
    mutate(!!iso := value) %>%
    select(-value)
  
  cut <- surv_cutpoint(
    df_iso2,
    time = "survival_time",
    event = "event",
    variables = iso,
    minprop = minprop
  )
  
  cut_value <- extract_cutpoint_value(cut, iso)
  cut_value_round <- round(cut_value, 3)
  
  df_cat <- as.data.frame(surv_categorize(cut))
  df_cat$group <- factor(df_cat[[iso]], levels = c("low", "high"))
  
  n_low  <- sum(df_cat$group == "low")
  n_high <- sum(df_cat$group == "high")
  
  fit <- survfit(Surv(survival_time, event) ~ group, data = df_cat)
  
  tiff(
    file.path(figS7A_dir, paste0("KM_", iso, "_SCAN_min15pct.tiff")),
    width = 2000, height = 1800, res = 300
  )
  
  print(
    ggsurvplot(
      fit, data = df_cat,
      pval = TRUE, risk.table = TRUE,
      risk.table.height = 0.25,
      break.time.by = 1000,
      palette = c("blue", "red"),
      title = paste0(
        "TCGA-BRCA OS by ", iso,
        " (cutpoint = ", cut_value_round,
        "; min group = 15%)"
      ),
      xlab = "Days",
      ylab = "Overall survival probability",
      xlim = c(0, 3000),
      legend.title = paste0(iso, " group"),
      legend.labs = c("Low", "High"),
      ggtheme = theme_minimal()
    )
  )
  
  dev.off()
  
  cox <- coxph(Surv(survival_time, event) ~ group, data = df_cat)
  cox_sum <- summary(cox)
  
  results[[iso]] <- data.frame(
    isoform = iso,
    cutpoint = cut_value,
    cutpoint_round = cut_value_round,
    minprop = minprop,
    N = cox$n,
    events = cox_sum$nevent,
    n_low = n_low,
    n_high = n_high,
    HR_high_vs_low = exp(coef(cox)),
    CI_low = exp(confint(cox))[1],
    CI_high = exp(confint(cox))[2],
    p_value = cox_sum$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
}

res_df <- bind_rows(results)

# ============================================================
# 10) Multiple testing correction
# ============================================================

res_df <- res_df %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  arrange(p_adj_BH)

write.csv(
  res_df,
  file.path(tableS7A_dir, "TableS_TCGA_Cox_by_isoform_SCAN_min15pct.csv"),
  row.names = FALSE
)

# ============================================================
# 11) Robustness evaluation across minprop values
#     (including group balance assessment)
# ============================================================

minprops <- c(0.10, 0.15, 0.20)
robust_results <- list()

for (mp in minprops) {
  
  message("[INFO] Robustness scan with minprop = ", mp)
  
  mp_results <- list()
  
  for (iso in isoforms) {
    
    if (!iso %in% colnames(df)) next
    
    df_iso <- df %>%
      filter(!is.na(.data[[iso]])) %>%
      transmute(
        survival_time = survival_time,
        event = as.integer(event),
        value = .data[[iso]]
      )
    
    df_iso2 <- df_iso %>%
      mutate(!!iso := value) %>%
      select(-value)
    
    cut <- surv_cutpoint(
      df_iso2,
      time = "survival_time",
      event = "event",
      variables = iso,
      minprop = mp
    )
    
    cut_value <- extract_cutpoint_value(cut, iso)
    
    df_cat <- as.data.frame(surv_categorize(cut))
    df_cat$group <- factor(df_cat[[iso]], levels = c("low", "high"))
    
    if (length(unique(df_cat$group)) < 2) next
    
    # Group sizes
    n_low  <- sum(df_cat$group == "low", na.rm = TRUE)
    n_high <- sum(df_cat$group == "high", na.rm = TRUE)
    N_tot  <- n_low + n_high
    
    # Balance metric: proportion of largest group
    prop_max_group <- max(n_low, n_high) / N_tot
    
    cox <- coxph(Surv(survival_time, event) ~ group, data = df_cat)
    cox_sum <- summary(cox)
    
    mp_results[[iso]] <- data.frame(
      isoform = iso,
      minprop = mp,
      HR = exp(coef(cox)),
      p_value = cox_sum$coefficients[, "Pr(>|z|)"],
      n_low = n_low,
      n_high = n_high,
      prop_max_group = prop_max_group,
      stringsAsFactors = FALSE
    )
  }
  
  robust_results[[as.character(mp)]] <- bind_rows(mp_results)
}

robust_df <- bind_rows(robust_results)

# ============================================================
# 12) Robustness + balance summary per isoform
# ============================================================

robust_summary <- robust_df %>%
  mutate(
    HR_direction = ifelse(HR > 1, "High_worse", "High_better")
  ) %>%
  group_by(isoform) %>%
  summarise(
    n_tested = n(),
    HR_direction_consistent = length(unique(HR_direction)) == 1,
    n_significant = sum(p_value < 0.05, na.rm = TRUE),
    min_p_value = min(p_value, na.rm = TRUE),
    
    # Balance metrics
    median_max_group_prop = median(prop_max_group, na.rm = TRUE),
    min_max_group_prop = min(prop_max_group, na.rm = TRUE),
    max_max_group_prop = max(prop_max_group, na.rm = TRUE),
    
    robust = HR_direction_consistent & n_significant >= 2,
    .groups = "drop"
  ) %>%
  arrange(desc(robust), min_p_value)

# ============================================================
# 13) Save robustness summary table
# ============================================================

write.csv(
  robust_summary,
  file.path(
    tableS7A_dir,
    "TableS_TCGA_Isoform_Robustness_minprop_summary.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 14) Console report (robustness + balance)
# ============================================================

message("\n[ROBUSTNESS & BALANCE SUMMARY]")

for (i in seq_len(nrow(robust_summary))) {
  
  r <- robust_summary[i, ]
  
  message(
    sprintf(
      paste0(
        "Isoform %s | Robust: %s | Significant in %d/%d | ",
        "HR direction consistent: %s | ",
        "Max group proportion (median [range]): %.2f [%.2f–%.2f]"
      ),
      r$isoform,
      r$robust,
      r$n_significant,
      r$n_tested,
      r$HR_direction_consistent,
      r$median_max_group_prop,
      r$min_max_group_prop,
      r$max_max_group_prop
    )
  )
}

message("\nFigure S7A isoform-level robustness and balance evaluation completed.")

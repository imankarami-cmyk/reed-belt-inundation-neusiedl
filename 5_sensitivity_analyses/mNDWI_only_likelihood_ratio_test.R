# ============================================================
# mNDWI-ONLY vs mNDWI+EVI vs NDWI+EVI — COMPARISON
# ------------------------------------------------------------
# PURPOSE: tests whether EVI adds anything to mNDWI by fitting
# a mNDWI-only logistic regression and comparing it against
# the two-predictor models from Sections 5.4 and 5.5.
#
# Follows EXACTLY the same structure as 5_5_re.R:
#   - Same decile labeling
#   - Same full-dataset fit
#   - Same gauge-validation exclusion of training rows
#   - Same Youden J threshold selection
#   - Same AUC / sensitivity / specificity / precision metrics
#
# ADDITIONAL TEST: likelihood ratio test comparing mNDWI-only
# against mNDWI+EVI to formally test whether EVI contributes
# statistically significant information beyond mNDWI alone.
# ============================================================

library(tidyverse)
library(pROC)

BASE      <- "<PROJECT_ROOT>/result/writting/new/Data/"
GAUGE     <- "<PROJECT_ROOT>/Data/"
IN_CSV    <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
GAUGE_CSV <- paste0(GAUGE, "all_stations_1984_2025.csv")
OUT_DIR   <- paste0(BASE, "mNDWI_only_comparison/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix_raw <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI)) %>%
  mutate(.uid = paste(pixel_id, date, sep = "_"))

gauge <- read_csv(GAUGE_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), level_mua = as.numeric(level_mua),
         station = as.character(station)) %>%
  filter(!is.na(date), !is.na(level_mua)) %>%
  distinct(date, station, .keep_all = TRUE)

# ------------------------------------------------------------
# 1. DECILE LABELING — mNDWI labels used for all three models
#    that involve mNDWI (mNDWI-only and mNDWI+EVI share the
#    same labels; NDWI+EVI uses NDWI labels as before)
# ------------------------------------------------------------

label_decile <- function(df, index_col) {
  df %>%
    group_by(pixel_id) %>%
    mutate(
      p10 = quantile(.data[[index_col]], 0.10, na.rm = TRUE),
      p90 = quantile(.data[[index_col]], 0.90, na.rm = TRUE),
      TCL_label_new = case_when(
        .data[[index_col]] >= p90 ~ 1L,
        .data[[index_col]] <= p10 ~ 0L,
        TRUE ~ NA_integer_
      )
    ) %>% ungroup()
}

ndwi_labeled  <- label_decile(pix_raw, "NDWI")
mndwi_labeled <- label_decile(pix_raw, "mNDWI")

# Three models: original two + mNDWI-only
# mNDWI-only uses the same mNDWI decile labels as mNDWI+EVI
model_defs <- list(
  NDWI_EVI   = list(formula = TCL_label_new ~ NDWI + EVI,
                    labeled = ndwi_labeled,
                    label   = "NDWI + EVI"),
  mNDWI_EVI  = list(formula = TCL_label_new ~ mNDWI + EVI,
                    labeled = mndwi_labeled,
                    label   = "mNDWI + EVI"),
  mNDWI_only = list(formula = TCL_label_new ~ mNDWI,
                    labeled = mndwi_labeled,
                    label   = "mNDWI only")
)

# ------------------------------------------------------------
# 2. FIT ALL THREE MODELS ON FULL LABELED SET
# ------------------------------------------------------------

fit_list   <- list()
train_uids <- list()

for (m in names(model_defs)) {
  dat <- model_defs[[m]]$labeled %>% filter(!is.na(TCL_label_new))
  fit_list[[m]]   <- glm(model_defs[[m]]$formula,
                         data = dat, family = binomial(link = "logit"))
  train_uids[[m]] <- dat$.uid
  cat("\n--- Model:", model_defs[[m]]$label, "---\n")
  print(summary(fit_list[[m]])$coefficients)
}

# ------------------------------------------------------------
# 3. LIKELIHOOD RATIO TEST: mNDWI-only vs mNDWI+EVI
#    Tests formally whether EVI adds significant information
#    beyond mNDWI in the improved model.
#    H0: EVI coefficient = 0 (mNDWI-only is sufficient)
#    H1: EVI adds significant information (mNDWI+EVI is better)
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("LIKELIHOOD RATIO TEST: mNDWI-only vs mNDWI+EVI\n")
cat("H0: EVI adds no significant information beyond mNDWI alone\n")
cat("============================================================\n")

# Both models must be fitted on the same data
# Use the mNDWI-labeled dataset for both
dat_mndwi <- mndwi_labeled %>% filter(!is.na(TCL_label_new))

fit_mndwi_only <- glm(TCL_label_new ~ mNDWI,
                       data = dat_mndwi, family = binomial())
fit_mndwi_evi  <- glm(TCL_label_new ~ mNDWI + EVI,
                       data = dat_mndwi, family = binomial())

lrt <- anova(fit_mndwi_only, fit_mndwi_evi, test = "Chisq")
cat("\nAnova (likelihood ratio test):\n")
print(lrt)

cat("\nInterpretation:\n")
cat("If p-value < 0.05: EVI adds statistically significant information\n")
cat("If p-value >= 0.05: mNDWI-only is sufficient — EVI is redundant\n\n")

lrt_result <- data.frame(
  comparison    = "mNDWI-only vs mNDWI+EVI",
  df_diff       = lrt$Df[2],
  deviance_diff = round(lrt$Deviance[2], 4),
  p_value       = round(lrt$`Pr(>Chi)`[2], 6),
  conclusion    = ifelse(lrt$`Pr(>Chi)`[2] < 0.05,
                         "EVI significant — keep mNDWI+EVI",
                         "EVI not significant — mNDWI-only sufficient")
)
write_csv(lrt_result, paste0(OUT_DIR, "likelihood_ratio_test_mNDWI_vs_mNDWI_EVI.csv"))
print(lrt_result)

# ------------------------------------------------------------
# 4. GAUGE VALIDATION FOR ALL THREE MODELS
#    Identical to 5_5_re.R — excludes training rows
# ------------------------------------------------------------

score_full <- function(m) {
  fit    <- fit_list[[m]]
  full_df <- pix_raw
  full_df$pred_prob <- predict(fit, newdata = full_df, type = "response")
  full_df %>%
    transmute(model = m, pixel_id, date, elevation, station,
              pred_prob, .uid,
              used_in_training = .uid %in% train_uids[[m]])
}

scored_all <- bind_rows(lapply(names(model_defs), score_full))

full_validation_all <- scored_all %>%
  left_join(gauge %>% select(date, station, level_mua),
            by = c("date", "station")) %>%
  filter(!is.na(level_mua), !is.na(elevation)) %>%
  mutate(
    water_present = if_else(level_mua >= elevation, 1L, 0L),
    model_label   = sapply(model, function(m) model_defs[[m]]$label)
  )

# Exclude training rows — same as primary validation in 5_5_re.R
full_validation <- full_validation_all %>% filter(!used_in_training)

# ------------------------------------------------------------
# 5. COMPUTE METRICS FOR ALL THREE MODELS
# ------------------------------------------------------------

results_list <- list()

for (m in names(model_defs)) {
  sub <- full_validation %>% filter(model == m)

  if (nrow(sub) < 50 || n_distinct(sub$water_present) < 2) {
    cat("Skipping", m, "— insufficient validation data\n")
    next
  }

  roc_o   <- roc(sub$water_present, sub$pred_prob, quiet = TRUE)
  auc_val <- as.numeric(roc_o$auc)

  youden  <- coords(roc_o, "best", best.method = "youden",
                    ret = c("threshold", "sensitivity", "specificity"))
  opt_thr <- as.numeric(youden["threshold"])
  sens    <- as.numeric(youden["sensitivity"])
  spec    <- as.numeric(youden["specificity"])

  pred_bin <- as.integer(sub$pred_prob >= opt_thr)
  tp <- sum(pred_bin == 1 & sub$water_present == 1)
  fp <- sum(pred_bin == 1 & sub$water_present == 0)
  tn <- sum(pred_bin == 0 & sub$water_present == 0)
  fn <- sum(pred_bin == 0 & sub$water_present == 1)
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  f1   <- if ((2*tp + fp + fn) > 0) 2*tp / (2*tp + fp + fn) else NA_real_

  results_list[[m]] <- tibble(
    model       = m,
    model_label = model_defs[[m]]$label,
    n_rows      = nrow(sub),
    AUC         = round(auc_val, 4),
    threshold   = round(opt_thr, 4),
    sensitivity = round(sens,    4),
    specificity = round(spec,    4),
    precision   = round(prec,    4),
    f1_score    = round(f1,      4)
  )
}

results_table <- bind_rows(results_list)
write_csv(results_table, paste0(OUT_DIR, "comparison_all_three_models.csv"))

cat("\n============================================================\n")
cat("GAUGE VALIDATION COMPARISON — ALL THREE MODELS\n")
cat("(training rows excluded — same as Section 5.5 primary result)\n")
cat("============================================================\n")
print(results_table)

cat("\nKey things to look for:\n")
cat("1. AUC mNDWI-only vs mNDWI+EVI:\n")
cat("   If similar (within ~0.005): EVI is genuinely redundant\n")
cat("   If mNDWI+EVI noticeably higher: EVI still adds value despite\n")
cat("   the non-significant coefficient in the decile-labeled fit\n\n")
cat("2. Likelihood ratio test p-value:\n")
cat("   If p < 0.05 but AUC difference is tiny: EVI is statistically\n")
cat("   significant but practically negligible\n")
cat("   If p >= 0.05: mNDWI-only model is formally sufficient\n\n")
cat("3. Threshold values:\n")
cat("   If mNDWI-only threshold differs substantially from mNDWI+EVI\n")
cat("   threshold, EVI was contributing to calibration even if AUC similar\n\n")

cat("All outputs saved to:", OUT_DIR, "\n")

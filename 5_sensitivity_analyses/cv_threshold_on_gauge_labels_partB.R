# ============================================================
# SECTION 5.6 — REPEATED CROSS-VALIDATION: GAUGE-THRESHOLD VERSION
# ------------------------------------------------------------
# PURPOSE: comparison-only script to test whether switching
# threshold selection in Part B from DECILE LABELS to GAUGE LABELS
# changes the reported results meaningfully.
#
# DIFFERENCE FROM 5_6_re.R (the main script):
#   Original: threshold selected on partB using TCL_label_new
#             (decile-based spectral labels)
#   This version: threshold selected on partB using water_present
#             (gauge water level >= pixel elevation)
#
# AUC should stay similar (AUC is threshold-independent).
# Sensitivity, specificity, precision, and threshold VALUES
# may change because the threshold is now calibrated to the
# gauge criterion rather than the spectral distribution.
#
# DO NOT replace 5_6_re.R with this script. Use only for
# comparison. If results are similar, the original design
# is confirmed adequate. If substantially different, note
# the difference as a sensitivity analysis finding.
# ============================================================

library(tidyverse)
library(pROC)

BASE      <- "<PROJECT_ROOT>/result/writting/new/Data/"
GAUGE     <- "<PROJECT_ROOT>/Data/"
IN_CSV    <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
GAUGE_CSV <- paste0(GAUGE, "all_stations_1984_2025.csv")
OUT_DIR   <- paste0(BASE, "5.6_gauge_threshold_comparison/")
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
# 1. DECILE LABELING (identical to original)
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

model_defs <- list(
  NDWI_EVI  = list(formula = TCL_label_new ~ NDWI + EVI,
                   labeled = ndwi_labeled,  label = "NDWI + EVI"),
  mNDWI_EVI = list(formula = TCL_label_new ~ mNDWI + EVI,
                   labeled = mndwi_labeled, label = "mNDWI + EVI")
)

# ------------------------------------------------------------
# 2. 10-ROUND CV WITH GAUGE-BASED THRESHOLD SELECTION
#    KEY CHANGE: Part B threshold uses water_present (gauge label)
#    instead of TCL_label_new (decile label)
# ------------------------------------------------------------

N_ROUNDS <- 10
set.seed(99)   # same seed as original for fair comparison

run_cv_gauge_thr <- function(m) {
  spec         <- model_defs[[m]]
  pix_model_base <- spec$labeled %>% filter(!is.na(TCL_label_new))
  full_df      <- pix_raw

  cv_rows <- list()
  for (r in seq_len(N_ROUNDS)) {
    n_all <- nrow(pix_model_base)
    idx   <- sample(n_all)
    cutA  <- floor(0.50 * n_all)
    cutB  <- floor(0.75 * n_all)

    partA <- pix_model_base[idx[1:cutA], ]
    partB <- pix_model_base[idx[(cutA + 1):cutB], ]

    # --- Step 1: fit coefficients on Part A (identical to original) ---
    fit_A <- tryCatch(glm(spec$formula, data = partA, family = binomial),
                      error = function(e) NULL)
    if (is.null(fit_A) || !fit_A$converged) next

    # --- Step 2: select threshold on Part B using GAUGE LABELS ---
    # KEY CHANGE: join gauge water level to partB and compute
    # water_present = level_mua >= elevation instead of using
    # TCL_label_new (the decile label)
    partB_gauge <- partB %>%
      left_join(gauge %>% select(date, station, level_mua),
                by = c("date", "station")) %>%
      filter(!is.na(level_mua), !is.na(elevation)) %>%
      mutate(water_present = if_else(level_mua >= elevation, 1L, 0L))

    # Skip this round if Part B has too few gauge-matched rows
    # or only one class present after gauge joining
    if (nrow(partB_gauge) < 50 ||
        n_distinct(partB_gauge$water_present) < 2) next

    prob_B_gauge <- predict(fit_A, partB_gauge, type = "response")

    roc_B <- tryCatch(
      roc(partB_gauge$water_present, prob_B_gauge, quiet = TRUE),
      error = function(e) NULL)
    if (is.null(roc_B)) next

    coords_B <- coords(roc_B, "all",
                       ret = c("threshold", "sensitivity", "specificity"))
    best_idx <- which.max(
      coords_B$sensitivity + coords_B$specificity - 1)[1]
    thr_r <- coords_B$threshold[best_idx]

    # --- Step 3: evaluate on Part C (identical to original) ---
    used_uids_r <- c(partA$.uid, partB$.uid)

    full_df$pred_prob <- predict(fit_A, newdata = full_df, type = "response")
    full_val_r <- full_df %>%
      filter(!(.uid %in% used_uids_r)) %>%
      left_join(gauge %>% select(date, station, level_mua),
                by = c("date", "station")) %>%
      filter(!is.na(level_mua), !is.na(elevation)) %>%
      mutate(water_present = if_else(level_mua >= elevation, 1L, 0L))

    if (nrow(full_val_r) < 100 ||
        n_distinct(full_val_r$water_present) < 2) next

    roc_C <- tryCatch(
      roc(full_val_r$water_present, full_val_r$pred_prob, quiet = TRUE),
      error = function(e) NULL)
    if (is.null(roc_C)) next

    auc_r  <- as.numeric(roc_C$auc)
    pred_C <- as.integer(full_val_r$pred_prob >= thr_r)
    tp <- sum(pred_C == 1 & full_val_r$water_present == 1)
    fn <- sum(pred_C == 0 & full_val_r$water_present == 1)
    tn <- sum(pred_C == 0 & full_val_r$water_present == 0)
    fp <- sum(pred_C == 1 & full_val_r$water_present == 0)
    sens_r <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    spec_r <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    prec_r <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_

    cv_rows[[r]] <- data.frame(
      model        = m,
      round        = r,
      AUC          = round(auc_r,  4),
      sensitivity  = round(sens_r, 4),
      specificity  = round(spec_r, 4),
      precision    = round(prec_r, 4),
      threshold    = round(thr_r,  4),
      n_eval_rows  = nrow(full_val_r),
      n_partB_gauge = nrow(partB_gauge)
    )
  }
  bind_rows(cv_rows)
}

cv_all <- bind_rows(lapply(names(model_defs), run_cv_gauge_thr)) %>%
  mutate(model_label = sapply(model,
                              function(m) model_defs[[m]]$label))

write_csv(cv_all,
          paste0(OUT_DIR, "cv_gauge_threshold_all_rounds.csv"))
print(cv_all)

# ------------------------------------------------------------
# 3. SUMMARY TABLE — mean +/- SD across 10 rounds
# ------------------------------------------------------------

cv_summary <- cv_all %>%
  group_by(model, model_label) %>%
  summarise(
    AUC_mean  = round(mean(AUC,         na.rm = TRUE), 4),
    AUC_sd    = round(sd(AUC,           na.rm = TRUE), 4),
    sens_mean = round(mean(sensitivity, na.rm = TRUE), 4),
    sens_sd   = round(sd(sensitivity,   na.rm = TRUE), 4),
    spec_mean = round(mean(specificity, na.rm = TRUE), 4),
    spec_sd   = round(sd(specificity,   na.rm = TRUE), 4),
    prec_mean = round(mean(precision,   na.rm = TRUE), 4),
    prec_sd   = round(sd(precision,     na.rm = TRUE), 4),
    thr_mean  = round(mean(threshold,   na.rm = TRUE), 4),
    thr_sd    = round(sd(threshold,     na.rm = TRUE), 4),
    n_rounds  = n(),
    .groups   = "drop"
  )

write_csv(cv_summary,
          paste0(OUT_DIR, "cv_gauge_threshold_summary.csv"))

cat("\n=== GAUGE-THRESHOLD CV SUMMARY ===\n")
print(cv_summary)

# ------------------------------------------------------------
# 4. COMPARISON TABLE: original (decile threshold) vs
#    this version (gauge threshold)
#    Paste your original 5_6_re.R summary output here to compare
# ------------------------------------------------------------

original <- tibble(
  model_label = c("NDWI + EVI", "mNDWI + EVI"),
  AUC_mean_original  = c(0.7396, 0.8078),
  AUC_sd_original    = c(0.0008, 0.0003),
  sens_mean_original = c(0.6803, 0.8132),
  thr_mean_original  = c(0.4408, 0.3516)
)

comparison <- cv_summary %>%
  select(model_label, AUC_mean, AUC_sd, sens_mean, thr_mean) %>%
  left_join(original, by = "model_label") %>%
  mutate(
    AUC_change  = round(AUC_mean  - AUC_mean_original,  4),
    sens_change = round(sens_mean - sens_mean_original,  4),
    thr_change  = round(thr_mean  - thr_mean_original,   4)
  )

cat("\n=== COMPARISON: decile-threshold vs gauge-threshold ===\n")
cat("Positive = gauge-threshold version is higher\n\n")
print(comparison %>% select(model_label, AUC_mean, AUC_mean_original,
                             AUC_change, sens_change, thr_change))

write_csv(comparison,
          paste0(OUT_DIR, "cv_comparison_decile_vs_gauge_threshold.csv"))

cat("\nIf AUC_change is close to zero: the original design is adequate.\n")
cat("If AUC_change is substantial (>0.01): consider noting as a\n")
cat("sensitivity finding rather than changing the main results.\n")
cat("\nAll outputs saved to:", OUT_DIR, "\n")

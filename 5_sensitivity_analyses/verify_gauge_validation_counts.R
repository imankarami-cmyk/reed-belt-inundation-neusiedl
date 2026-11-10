# ============================================================
# SECTION 5.5 VERIFICATION
# ------------------------------------------------------------
# Checks:
# 1. Why NDWI+EVI has 252,486 held-out rows and
#    mNDWI+EVI has 252,355 — are these numbers correct?
# 2. Confirms which data the threshold was selected from
# 3. Recomputes AUC, threshold, sensitivity, specificity,
#    precision, accuracy, F1 to verify Table 9 numbers
# ============================================================

library(tidyverse)
library(pROC)

BASE      <- "<PROJECT_ROOT>/result/writting/new/Data/"
GAUGE     <- "<PROJECT_ROOT>/Data/"
IN_CSV    <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
GAUGE_CSV <- paste0(GAUGE, "all_stations_1984_2025.csv")
OUT_DIR   <- paste0(BASE, "5.5_verification/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix_raw <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI)) %>%
  mutate(.uid = paste(pixel_id, date, sep = "_"))

gauge <- read_csv(GAUGE_CSV, show_col_types = FALSE) %>%
  mutate(date      = as.Date(date),
         level_mua = as.numeric(level_mua),
         station   = as.character(station)) %>%
  filter(!is.na(date), !is.na(level_mua)) %>%
  distinct(date, station, .keep_all = TRUE)

cat("=== TOTAL OBSERVATIONS IN pix_raw ===\n")
cat("Total rows:", nrow(pix_raw), "\n\n")

# ------------------------------------------------------------
# 1. DECILE LABELING — same as Section 5.3
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
        TRUE                      ~ NA_integer_
      )
    ) %>% ungroup()
}

ndwi_labeled  <- label_decile(pix_raw, "NDWI")
mndwi_labeled <- label_decile(pix_raw, "mNDWI")

# ------------------------------------------------------------
# 2. CHECK LABELED COUNTS
# ------------------------------------------------------------

ndwi_train  <- ndwi_labeled  %>% filter(!is.na(TCL_label_new))
mndwi_train <- mndwi_labeled %>% filter(!is.na(TCL_label_new))

cat("=== LABELED OBSERVATION COUNTS ===\n")
cat("NDWI+EVI  labeled (training) rows:", nrow(ndwi_train),  "\n")
cat("mNDWI+EVI labeled (training) rows:", nrow(mndwi_train), "\n\n")

cat("Why are these potentially different?\n")
cat("NDWI labels from NDWI 10th/90th percentiles\n")
cat("mNDWI labels from mNDWI 10th/90th percentiles\n")
cat("Boundary cases (exactly at percentile) may differ slightly\n\n")

# ------------------------------------------------------------
# 3. FIT MODELS ON ALL LABELED DATA
# ------------------------------------------------------------

fit_ndwi  <- glm(TCL_label_new ~ NDWI + EVI,
                 data = ndwi_train,  family = binomial())
fit_mndwi <- glm(TCL_label_new ~ mNDWI + EVI,
                 data = mndwi_train, family = binomial())

cat("=== FITTED COEFFICIENTS ===\n")
cat("NDWI+EVI:\n");  print(coef(fit_ndwi))
cat("mNDWI+EVI:\n"); print(coef(fit_mndwi))

# ------------------------------------------------------------
# 4. APPLY TO FULL DATASET AND JOIN GAUGE
# ------------------------------------------------------------

pix_raw$pred_ndwi  <- predict(fit_ndwi,  newdata = pix_raw, type = "response")
pix_raw$pred_mndwi <- predict(fit_mndwi, newdata = pix_raw, type = "response")

full_gauge <- pix_raw %>%
  left_join(gauge %>% select(date, station, level_mua),
            by = c("date", "station")) %>%
  filter(!is.na(level_mua), !is.na(elevation)) %>%
  mutate(water_present = if_else(level_mua >= elevation, 1L, 0L))

cat("\n=== FULL GAUGE-JOINED DATASET ===\n")
cat("Total gauge-matched rows:", nrow(full_gauge), "\n\n")

# ------------------------------------------------------------
# 5. HELD-OUT SETS: EXCLUDE TRAINING ROWS
# ------------------------------------------------------------

held_ndwi  <- full_gauge %>% filter(!(.uid %in% ndwi_train$.uid))
held_mndwi <- full_gauge %>% filter(!(.uid %in% mndwi_train$.uid))

cat("=== HELD-OUT ROW COUNTS ===\n")
cat("NDWI+EVI  held-out rows:", nrow(held_ndwi),  "\n")
cat("mNDWI+EVI held-out rows:", nrow(held_mndwi), "\n")
cat("\nThe difference between these two numbers explains FV20.1:\n")
cat("Each model excludes its OWN labeled training rows.\n")
cat("NDWI and mNDWI label slightly different observations\n")
cat("(boundary cases differ between indices), so the held-out\n")
cat("sets differ slightly in size.\n\n")

# ------------------------------------------------------------
# 6. THRESHOLD SELECTION AND METRICS
# ------------------------------------------------------------

compute_metrics <- function(held_df, pred_col, model_name) {
  roc_obj <- roc(held_df$water_present, held_df[[pred_col]], quiet = TRUE)
  auc_val <- as.numeric(roc_obj$auc)

  # Youden J threshold from held-out gauge data
  coords_all <- coords(roc_obj, "all",
                       ret = c("threshold", "sensitivity", "specificity"))
  best_idx   <- which.max(
    coords_all$sensitivity + coords_all$specificity - 1)[1]
  thr_val    <- coords_all$threshold[best_idx]

  pred_bin <- as.integer(held_df[[pred_col]] >= thr_val)
  tp <- sum(pred_bin == 1 & held_df$water_present == 1)
  fn <- sum(pred_bin == 0 & held_df$water_present == 1)
  tn <- sum(pred_bin == 0 & held_df$water_present == 0)
  fp <- sum(pred_bin == 1 & held_df$water_present == 0)

  sens  <- tp / (tp + fn)
  spec  <- tn / (tn + fp)
  prec  <- tp / (tp + fp)
  acc   <- (tp + tn) / (tp + tn + fp + fn)
  f1    <- 2 * tp / (2 * tp + fp + fn)

  cat("---", model_name, "---\n")
  cat("Threshold (Youden J on gauge data):", round(thr_val, 4), "\n")
  cat("NOTE: threshold selected from SAME held-out gauge data\n")
  cat("used for evaluation — AUC is unaffected (threshold-independent)\n")
  cat("but sensitivity/specificity reflect this threshold choice.\n")
  cat("AUC:        ", round(auc_val, 4), "\n")
  cat("Accuracy:   ", round(acc,     4), "\n")
  cat("Sensitivity:", round(sens,    4), "\n")
  cat("Specificity:", round(spec,    4), "\n")
  cat("Precision:  ", round(prec,    4), "\n")
  cat("F1:         ", round(f1,      4), "\n\n")

  tibble(model = model_name, threshold = round(thr_val, 4),
         AUC = round(auc_val, 4), accuracy = round(acc, 4),
         sensitivity = round(sens, 4), specificity = round(spec, 4),
         precision = round(prec, 4), f1 = round(f1, 4),
         n_held_out = nrow(held_df))
}

cat("=== SECTION 5.5 METRIC VERIFICATION ===\n\n")
results <- bind_rows(
  compute_metrics(held_ndwi,  "pred_ndwi",  "NDWI+EVI"),
  compute_metrics(held_mndwi, "pred_mndwi", "mNDWI+EVI")
)

cat("=== FINAL TABLE FOR SECTION 5.5 ===\n")
print(results)
write_csv(results, paste0(OUT_DIR, "section_5_5_verified_metrics.csv"))

cat("\nSaved to:", OUT_DIR, "\n")

# ============================================================
# SECTION 5.6 — REPEATED CROSS-VALIDATION: GAUGE-BASED THRESHOLD
# ------------------------------------------------------------
# DESIGN (matches LOSO procedure in 5_7_loso.R):
#
#   Part A (50% of labeled data) → fit model coefficients
#
#   Threshold → apply formula to ALL remaining data
#               (full pix_raw MINUS Part A UIDs)
#               join gauge → water_present (level_mua >= elevation)
#               Youden J on GAUGE labels → threshold for this round
#
#   Evaluation → same remaining dataset as threshold step
#               AUC, sensitivity, specificity, precision
#               all computed against GAUGE labels
#
# This eliminates the separate Part B partition entirely.
# Threshold and evaluation share the same remaining dataset,
# exactly as LOSO uses train_full for threshold and held_full
# for evaluation — both gauge-based, no decile labels after
# the fitting step.
#
# This design directly answers FV2 and FV3 from the supervisor:
# the threshold is now gauge-based (not decile-based), and
# there is no contradiction between the threshold method here
# and the threshold method in Section 4.6.2.
# ============================================================

library(tidyverse)
library(pROC)

BASE      <- "<PROJECT_ROOT>/result/writting/new/Data/"
GAUGE     <- "<PROJECT_ROOT>/Data/"
IN_CSV    <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
GAUGE_CSV <- paste0(GAUGE, "all_stations_1984_2025.csv")
OUT_DIR   <- paste0(BASE, "5.6_gauge_threshold/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix_raw <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI)) %>%
  mutate(.uid = paste(pixel_id, date, sep = "_"))

gauge <- read_csv(GAUGE_CSV, show_col_types = FALSE) %>%
  mutate(date       = as.Date(date),
         level_mua  = as.numeric(level_mua),
         station    = as.character(station)) %>%
  filter(!is.na(date), !is.na(level_mua)) %>%
  distinct(date, station, .keep_all = TRUE)

# ------------------------------------------------------------
# 1. DECILE LABELING — used only for Part A (coefficient fitting)
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

model_defs <- list(
  NDWI_EVI  = list(formula = TCL_label_new ~ NDWI + EVI,
                   labeled = ndwi_labeled,
                   label   = "NDWI + EVI"),
  mNDWI_EVI = list(formula = TCL_label_new ~ mNDWI + EVI,
                   labeled = mndwi_labeled,
                   label   = "mNDWI + EVI")
)

# ------------------------------------------------------------
# 2. 10-ROUND REPEATED CROSS-VALIDATION
#    Part A (50% labeled) → fit coefficients via glm()
#    Remaining data       → gauge join → Youden J threshold
#                        → AUC and all metrics vs gauge labels
# ------------------------------------------------------------

N_ROUNDS <- 10
set.seed(99)

run_cv <- function(m) {
  spec           <- model_defs[[m]]
  pix_model_base <- spec$labeled %>% filter(!is.na(TCL_label_new))

  cv_rows <- list()

  for (r in seq_len(N_ROUNDS)) {

    # ── Part A: fit on 50% of labeled observations ──────────
    n_all <- nrow(pix_model_base)
    idx   <- sample(n_all)
    cutA  <- floor(0.50 * n_all)
    partA <- pix_model_base[idx[1:cutA], ]

    fit_A <- tryCatch(
      glm(spec$formula, data = partA, family = binomial),
      error = function(e) NULL)
    if (is.null(fit_A) || !fit_A$converged) next

    # ── Remaining data: full archive minus Part A UIDs ───────
    # Apply formula → predicted probability for every row
    remaining <- pix_raw %>%
      filter(!(.uid %in% partA$.uid)) %>%
      mutate(pred_prob = predict(fit_A, newdata = ., type = "response")) %>%
      left_join(gauge %>% select(date, station, level_mua),
                by = c("date", "station")) %>%
      filter(!is.na(level_mua), !is.na(elevation)) %>%
      mutate(water_present = if_else(level_mua >= elevation, 1L, 0L))

    if (nrow(remaining) < 100 ||
        n_distinct(remaining$water_present) < 2) next

    # ── Threshold: Youden J on GAUGE labels ─────────────────
    # Identical procedure to LOSO threshold selection
    roc_gauge <- tryCatch(
      roc(remaining$water_present, remaining$pred_prob, quiet = TRUE),
      error = function(e) NULL)
    if (is.null(roc_gauge)) next

    coords_g  <- coords(roc_gauge, "all",
                        ret = c("threshold", "sensitivity", "specificity"))
    best_idx  <- which.max(
      coords_g$sensitivity + coords_g$specificity - 1)[1]
    thr_r     <- coords_g$threshold[best_idx]

    # ── Evaluation: AUC and metrics on remaining vs gauge ────
    auc_r  <- as.numeric(roc_gauge$auc)   # reuse same ROC object
    pred_r <- as.integer(remaining$pred_prob >= thr_r)

    tp <- sum(pred_r == 1 & remaining$water_present == 1)
    fn <- sum(pred_r == 0 & remaining$water_present == 1)
    tn <- sum(pred_r == 0 & remaining$water_present == 0)
    fp <- sum(pred_r == 1 & remaining$water_present == 0)

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
      n_eval_rows  = nrow(remaining)
    )
  }
  bind_rows(cv_rows)
}

cv_all <- bind_rows(lapply(names(model_defs), run_cv)) %>%
  mutate(model_label = sapply(model, function(m) model_defs[[m]]$label))

write_csv(cv_all, paste0(OUT_DIR, "cv_all_rounds.csv"))
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

write_csv(cv_summary, paste0(OUT_DIR, "cv_summary.csv"))
cat("\n=== CV SUMMARY (gauge-based threshold) ===\n")
print(cv_summary)

# ------------------------------------------------------------
# 4. FIGURE — AUC stability across rounds
# ------------------------------------------------------------

model_colors <- c("NDWI + EVI" = "#e07b39", "mNDWI + EVI" = "#01696f")

fig_auc <- ggplot(cv_all, aes(x = round, y = AUC, color = model_label)) +
  geom_line() +
  geom_point(size = 2) +
  scale_color_manual(values = model_colors) +
  coord_cartesian(ylim = c(0.5, 0.9)) +
  labs(title = "10-round cross-validation — AUC stability (gauge-based threshold)",
       x = "CV round", y = "AUC", color = "Model") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(paste0(OUT_DIR, "cv_auc_stability.png"),
       fig_auc, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------
# 5. FIGURE — metric distributions across 10 rounds
# ------------------------------------------------------------

cv_long <- cv_all %>%
  select(model_label, round, AUC, sensitivity, specificity, precision) %>%
  pivot_longer(c(AUC, sensitivity, specificity, precision),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(metric = recode(metric,
                         AUC         = "AUC",
                         sensitivity = "Sensitivity",
                         specificity = "Specificity",
                         precision   = "Precision"))

fig_box <- ggplot(cv_long,
                  aes(x = metric, y = value, fill = model_label)) +
  geom_boxplot(alpha = 0.85) +
  scale_fill_manual(values = model_colors) +
  labs(title = "Cross-validation metric distributions — gauge-based threshold",
       x = NULL, y = "Value", fill = "Model") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(paste0(OUT_DIR, "cv_metric_boxplot.png"),
       fig_box, width = 8, height = 6, dpi = 300)

cat("\nAll outputs saved to:", OUT_DIR, "\n")
cat("Run this script and paste the cv_summary.csv numbers into\n")
cat("Table 5.10 in the thesis, replacing the old decile-threshold results.\n")

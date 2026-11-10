# ============================================================
# SECTION 5.6 — REPEATED CROSS-VALIDATION: TABLES + PLOTS (REVISED)
# ------------------------------------------------------------
# FIX APPLIED (this is the most important one in the pipeline):
#
#  In the original script, `full_df <- pix_raw` was fixed OUTSIDE
#  the round loop, and every round's AUC/sensitivity/specificity
#  was computed against that same complete, unchanging dataset -
#  which fully contains the rows just used to fit the model
#  (partA) and pick the threshold (partB) in that same round.
#  That is why the reported SD across 10 rounds was ~0.0001: most
#  of the "test" set never changed round to round, so of course
#  the metric barely moved. It was measuring resampling noise in
#  a model whose evaluation set barely differs from its training
#  set, not generalization.
#
#  This version excludes, in every round, exactly the rows used
#  in that round's partA (fitting) and partB (threshold selection)
#  from the evaluation set, using a unique pixel+date row key
#  (.uid). Each round now genuinely evaluates on rows the model
#  did not see that round. Expect the SD to be larger and more
#  meaningful than before - that is a feature of the fix, not a
#  problem with it.
# ============================================================

library(tidyverse)
library(pROC)

BASE      <- "<PROJECT_ROOT>/result/writting/new/Data/"
GAUGE     <- "<PROJECT_ROOT>/Data/"
IN_CSV    <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
GAUGE_CSV <- paste0(GAUGE, "all_stations_1984_2025.csv")
OUT_DIR   <- paste0(BASE, "5.6/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix_raw <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI)) %>%
  mutate(.uid = paste(pixel_id, date, sep = "_"))   # NEW: unique row key

gauge <- read_csv(GAUGE_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), level_mua = as.numeric(level_mua),
         station = as.character(station)) %>%
  filter(!is.na(date), !is.na(level_mua)) %>%
  distinct(date, station, .keep_all = TRUE)

# ------------------------------------------------------------
# 1. DECILE LABELING (same as Sections 5.3-5.5)
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
  NDWI_EVI  = list(formula = TCL_label_new ~ NDWI + EVI,  labeled = ndwi_labeled,  label = "NDWI + EVI"),
  mNDWI_EVI = list(formula = TCL_label_new ~ mNDWI + EVI, labeled = mndwi_labeled, label = "mNDWI + EVI")
)

# ------------------------------------------------------------
# 2. 10-ROUND REPEATED CROSS-VALIDATION (FIXED)
#    50% fit / 25% threshold selection / evaluation on the
#    REMAINING rows only (fixed: excludes partA + partB uids)
# ------------------------------------------------------------

N_ROUNDS <- 10
set.seed(99)

run_cv <- function(m) {
  spec <- model_defs[[m]]
  pix_model_base <- spec$labeled %>% filter(!is.na(TCL_label_new))
  full_df <- pix_raw   # fixed pool to draw the held-out evaluation rows from

  cv_rows <- list()
  for (r in seq_len(N_ROUNDS)) {
    n_all <- nrow(pix_model_base)
    idx <- sample(n_all)
    cutA <- floor(0.50 * n_all)
    cutB <- floor(0.75 * n_all)

    partA <- pix_model_base[idx[1:cutA], ]
    partB <- pix_model_base[idx[(cutA + 1):cutB], ]

    fit_A <- tryCatch(glm(spec$formula, data = partA, family = binomial),
                      error = function(e) NULL)
    if (is.null(fit_A) || !fit_A$converged) next

    prob_B <- predict(fit_A, partB, type = "response")
    roc_B <- tryCatch(roc(partB$TCL_label_new, prob_B, quiet = TRUE), error = function(e) NULL)
    if (is.null(roc_B)) next
    coords_B <- coords(roc_B, "all", ret = c("threshold", "sensitivity", "specificity"))
    best_idx <- which.max(coords_B$sensitivity + coords_B$specificity - 1)[1]
    thr_r <- coords_B$threshold[best_idx]

    # ---- FIX: build the evaluation set excluding this round's
    #      fitting (partA) and threshold-selection (partB) rows ----
    used_uids_r <- c(partA$.uid, partB$.uid)

    full_df$pred_prob <- predict(fit_A, newdata = full_df, type = "response")
    full_val_r <- full_df %>%
      filter(!(.uid %in% used_uids_r)) %>%                       # <- the fix
      left_join(gauge %>% select(date, station, level_mua), by = c("date", "station")) %>%
      filter(!is.na(level_mua), !is.na(elevation)) %>%
      mutate(water_present = if_else(level_mua >= elevation, 1L, 0L))

    if (nrow(full_val_r) < 100 || n_distinct(full_val_r$water_present) < 2) next

    roc_C <- tryCatch(roc(full_val_r$water_present, full_val_r$pred_prob, quiet = TRUE),
                      error = function(e) NULL)
    if (is.null(roc_C)) next

    auc_r <- as.numeric(roc_C$auc)
    pred_C <- as.integer(full_val_r$pred_prob >= thr_r)
    tp <- sum(pred_C == 1 & full_val_r$water_present == 1)
    fn <- sum(pred_C == 0 & full_val_r$water_present == 1)
    tn <- sum(pred_C == 0 & full_val_r$water_present == 0)
    fp <- sum(pred_C == 1 & full_val_r$water_present == 0)
    sens_r <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    spec_r <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    prec_r <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_

    cv_rows[[r]] <- data.frame(
      model = m, round = r, AUC = round(auc_r, 4),
      sensitivity = round(sens_r, 4), specificity = round(spec_r, 4),
      precision = round(prec_r, 4),
      threshold = round(thr_r, 4),
      n_eval_rows = nrow(full_val_r)
    )
  }
  bind_rows(cv_rows)
}

cv_all <- bind_rows(lapply(names(model_defs), run_cv)) %>%
  mutate(model_label = sapply(model, function(m) model_defs[[m]]$label))

write_csv(cv_all, paste0(OUT_DIR, "table5_16_cv_all_rounds.csv"))
print(cv_all)

# ------------------------------------------------------------
# 3. TABLE 5.17 — CV SUMMARY (mean +/- SD across 10 rounds)
# ------------------------------------------------------------

cv_summary <- cv_all %>%
  group_by(model, model_label) %>%
  summarise(
    AUC_mean = round(mean(AUC, na.rm = TRUE), 4), AUC_sd = round(sd(AUC, na.rm = TRUE), 4),
    sens_mean = round(mean(sensitivity, na.rm = TRUE), 4), sens_sd = round(sd(sensitivity, na.rm = TRUE), 4),
    spec_mean = round(mean(specificity, na.rm = TRUE), 4), spec_sd = round(sd(specificity, na.rm = TRUE), 4),
    prec_mean = round(mean(precision, na.rm = TRUE), 4), prec_sd = round(sd(precision, na.rm = TRUE), 4),
    thr_mean = round(mean(threshold, na.rm = TRUE), 4), thr_sd = round(sd(threshold, na.rm = TRUE), 4),
    n_rounds = n(), .groups = "drop"
  )

write_csv(cv_summary, paste0(OUT_DIR, "table5_17_cv_summary.csv"))
print(cv_summary)

cat("\nIf AUC_sd here is noticeably larger than the ~0.0001-0.0003 reported\n")
cat("originally, that is expected and correct: each round's evaluation set\n")
cat("now genuinely excludes that round's training/threshold rows, so\n")
cat("round-to-round variation reflects real model instability rather than\n")
cat("re-scoring the same fixed dataset ten times.\n\n")

# ------------------------------------------------------------
# 4. FIGURE 5.12 — CV AUC STABILITY ACROSS ROUNDS
# ------------------------------------------------------------

model_colors <- c("NDWI + EVI" = "#e07b39", "mNDWI + EVI" = "#01696f")

fig5_12 <- ggplot(cv_all, aes(x = round, y = AUC, color = model_label)) +
  geom_line() + geom_point(size = 2) +
  scale_color_manual(values = model_colors) +
  coord_cartesian(ylim = c(0.5, 0.9)) +
  labs(title = "10-round cross-validation AUC stability (held-out rows only)",
       x = "CV round", y = "AUC", color = "Model") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_12_cv_stability.png"), fig5_12, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------
# 5. FIGURE 5.13 — CV METRIC DISTRIBUTIONS (boxplot: AUC, sens, spec, prec)
# ------------------------------------------------------------

cv_long <- cv_all %>%
  select(model_label, round, AUC, sensitivity, specificity, precision) %>%
  pivot_longer(c(AUC, sensitivity, specificity, precision), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, AUC = "AUC", sensitivity = "Sensitivity",
                          specificity = "Specificity", precision = "Precision"))

fig5_13 <- ggplot(cv_long, aes(x = metric, y = value, fill = model_label)) +
  geom_boxplot(alpha = 0.85) +
  scale_fill_manual(values = model_colors) +
  labs(title = "Distribution of cross-validation metrics across 10 rounds (held-out rows only)",
       x = NULL, y = "Value", fill = "Model") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_13_cv_metric_boxplot.png"), fig5_13, width = 8, height = 6, dpi = 300)

cat("Section 5.6 outputs saved to:", OUT_DIR, "\n")
cat("NOTE: for a spatial (not just resampling) generalization check, also\n")
cat("run 5_7_loso.R — that tests whether the model generalizes to a gauge\n")
cat("station it has never seen at all, which this section still does not test.\n")

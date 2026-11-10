# ============================================================
# SECTION 5.7 (NEW) — LEAVE-ONE-STATION-OUT CROSS-VALIDATION
# ------------------------------------------------------------
# WHY THIS IS NEEDED:
#   Section 5.6's repeated CV (even after the fix) only ever
#   resamples rows from the SAME six stations / 600 pixels the
#   model was built on. It tells you the fit is numerically
#   stable, but not whether the model would work at a station
#   it has never seen - which is the real claim implied by
#   "this model detects flooding under reed canopy at Lake
#   Neusiedl" rather than "at these six specific gauge sites."
#
#   This script fits each model on FIVE of the six stations and
#   evaluates it purely on the sixth, rotating through all six.
#   Both the model fit AND the Youden-J threshold are derived
#   only from the five training stations - the held-out station's
#   data is never used for anything except final scoring. This is
#   the standard leave-one-group-out design and directly answers
#   "does this generalize to a new location," which the original
#   pipeline never tested.
#
# Suggested thesis placement: new Section 5.6.4 / 5.7, immediately
# after the repeated cross-validation results, presented as a
# stronger, spatially independent complement to Section 5.6.
# ============================================================

library(tidyverse)
library(pROC)

BASE      <- "<PROJECT_ROOT>/result/writting/new/Data/"
GAUGE     <- "<PROJECT_ROOT>/Data/"
IN_CSV    <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
GAUGE_CSV <- paste0(GAUGE, "all_stations_1984_2025.csv")
OUT_DIR   <- paste0(BASE, "5.7_LOSO/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix_raw <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI))

gauge <- read_csv(GAUGE_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), level_mua = as.numeric(level_mua),
         station = as.character(station)) %>%
  filter(!is.na(date), !is.na(level_mua)) %>%
  distinct(date, station, .keep_all = TRUE)

stations <- sort(unique(pix_raw$station))
cat("Stations found for LOSO rotation:", paste(stations, collapse = ", "), "\n\n")
stopifnot(length(stations) == 6)  # sanity check against the six-station design in Section 4.3

# ------------------------------------------------------------
# 1. DECILE LABELING (same rule as Sections 5.3-5.6)
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
# 2. LEAVE-ONE-STATION-OUT LOOP
# ------------------------------------------------------------

run_loso <- function(m) {
  spec <- model_defs[[m]]
  labeled_base <- spec$labeled %>% filter(!is.na(TCL_label_new))

  fold_rows <- list()

  for (held_station in stations) {

    # --- fit + threshold selection: five training stations only ---
    train_labeled <- labeled_base %>% filter(station != held_station)

    fit <- tryCatch(glm(spec$formula, data = train_labeled, family = binomial()),
                     error = function(e) NULL)
    if (is.null(fit) || !fit$converged) {
      warning(paste("Model failed to converge for", m, "excluding", held_station)); next
    }

    # threshold chosen ONLY from the five training stations' gauge-validated data
    train_full <- pix_raw %>%
      filter(station != held_station) %>%
      mutate(pred_prob = predict(fit, newdata = ., type = "response")) %>%
      left_join(gauge %>% select(date, station, level_mua), by = c("date", "station")) %>%
      filter(!is.na(level_mua), !is.na(elevation)) %>%
      mutate(water_present = if_else(level_mua >= elevation, 1L, 0L))

    roc_train <- tryCatch(roc(train_full$water_present, train_full$pred_prob, quiet = TRUE),
                           error = function(e) NULL)
    if (is.null(roc_train)) next
    youden <- coords(roc_train, "best", best.method = "youden",
                      ret = c("threshold", "sensitivity", "specificity"))
    thr <- as.numeric(youden["threshold"])

    # --- evaluation: held-out station only, never seen during fitting
    #     or threshold selection ---
    held_full <- pix_raw %>%
      filter(station == held_station) %>%
      mutate(pred_prob = predict(fit, newdata = ., type = "response")) %>%
      left_join(gauge %>% select(date, station, level_mua), by = c("date", "station")) %>%
      filter(!is.na(level_mua), !is.na(elevation)) %>%
      mutate(water_present = if_else(level_mua >= elevation, 1L, 0L))

    if (nrow(held_full) < 30 || n_distinct(held_full$water_present) < 2) {
      warning(paste("Insufficient gauge-matched data to evaluate", m, "at", held_station))
      next
    }

    roc_held <- tryCatch(roc(held_full$water_present, held_full$pred_prob, quiet = TRUE),
                          error = function(e) NULL)
    auc_held <- if (is.null(roc_held)) NA_real_ else as.numeric(roc_held$auc)

    pred_bin <- as.integer(held_full$pred_prob >= thr)
    tp <- sum(pred_bin == 1 & held_full$water_present == 1)
    fn <- sum(pred_bin == 0 & held_full$water_present == 1)
    tn <- sum(pred_bin == 0 & held_full$water_present == 0)
    fp <- sum(pred_bin == 1 & held_full$water_present == 0)

    fold_rows[[held_station]] <- tibble(
      model = m, model_label = spec$label, held_out_station = held_station,
      n_train_pixels = n_distinct(train_labeled$pixel_id),
      n_eval_rows = nrow(held_full),
      threshold_from_training = round(thr, 4),
      AUC = round(auc_held, 4),
      sensitivity = round(if ((tp+fn)>0) tp/(tp+fn) else NA_real_, 4),
      specificity = round(if ((tn+fp)>0) tn/(tn+fp) else NA_real_, 4),
      precision   = round(if ((tp+fp)>0) tp/(tp+fp) else NA_real_, 4),
      accuracy    = round((tp+tn)/(tp+fp+tn+fn), 4)
    )
  }
  bind_rows(fold_rows)
}

loso_all <- bind_rows(lapply(names(model_defs), run_loso))

write_csv(loso_all, paste0(OUT_DIR, "table5_18_loso_per_station.csv"))
print(loso_all)

# ------------------------------------------------------------
# 3. TABLE 5.19 — LOSO SUMMARY (mean +/- SD across the 6 station-folds)
# ------------------------------------------------------------

loso_summary <- loso_all %>%
  group_by(model, model_label) %>%
  summarise(
    AUC_mean = round(mean(AUC, na.rm = TRUE), 4), AUC_sd = round(sd(AUC, na.rm = TRUE), 4),
    sens_mean = round(mean(sensitivity, na.rm = TRUE), 4), sens_sd = round(sd(sensitivity, na.rm = TRUE), 4),
    spec_mean = round(mean(specificity, na.rm = TRUE), 4), spec_sd = round(sd(specificity, na.rm = TRUE), 4),
    prec_mean = round(mean(precision, na.rm = TRUE), 4), prec_sd = round(sd(precision, na.rm = TRUE), 4),
    n_folds = n(), .groups = "drop"
  )

write_csv(loso_summary, paste0(OUT_DIR, "table5_19_loso_summary.csv"))
print(loso_summary)

cat("\nCompare AUC_sd here to Table 5.17's AUC_sd. A materially larger spread\n")
cat("across station-folds than across resampling rounds would indicate the\n")
cat("model is more sensitive to WHICH station it is applied to than the\n")
cat("resampling-only CV in Section 5.6 was able to reveal - this is exactly\n")
cat("the kind of result worth discussing explicitly in Chapter 6, whichever\n")
cat("way it comes out.\n\n")

# ------------------------------------------------------------
# 4. FIGURE 5.14 — PER-STATION LOSO AUC (bar chart)
# ------------------------------------------------------------

model_colors <- c("NDWI + EVI" = "#e07b39", "mNDWI + EVI" = "#01696f")

fig5_14 <- ggplot(loso_all, aes(x = held_out_station, y = AUC, fill = model_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(data = loso_summary, aes(yintercept = AUC_mean, color = model_label),
             linetype = "dashed", show.legend = FALSE) +
  scale_fill_manual(values = model_colors) +
  scale_color_manual(values = model_colors) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(title = "Leave-one-station-out AUC by held-out station",
       subtitle = "Dashed lines = mean AUC across all six station-folds",
       x = "Held-out station", y = "AUC", fill = "Model") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_14_loso_auc_by_station.png"), fig5_14, width = 9, height = 6, dpi = 300)

# ------------------------------------------------------------
# 5. FIGURE 5.15 — LOSO vs REPEATED-CV COMPARISON (optional, if
#    table5_17_cv_summary.csv from the fixed 5_6.R is available)
# ------------------------------------------------------------

cv_summary_path <- paste0(BASE, "5.6/table5_17_cv_summary.csv")
if (file.exists(cv_summary_path)) {
  cv_summary <- read_csv(cv_summary_path, show_col_types = FALSE) %>%
    mutate(method = "Repeated CV (resampling, Sec. 5.6)")
  loso_compare <- loso_summary %>%
    mutate(method = "LOSO (spatial, Sec. 5.7)")

  compare_df <- bind_rows(
    cv_summary %>% select(model_label, AUC_mean, AUC_sd, method),
    loso_compare %>% select(model_label, AUC_mean, AUC_sd, method)
  )

  fig5_15 <- ggplot(compare_df, aes(x = model_label, y = AUC_mean, fill = method)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.5) +
    geom_errorbar(aes(ymin = AUC_mean - AUC_sd, ymax = AUC_mean + AUC_sd),
                  position = position_dodge(width = 0.6), width = 0.15) +
    coord_cartesian(ylim = c(0.5, 0.9)) +
    labs(title = "Resampling CV vs. spatial (LOSO) generalization",
         x = NULL, y = "Mean AUC (+/- SD)", fill = NULL) +
    theme_minimal(base_size = 13) + theme(legend.position = "bottom")
  ggsave(paste0(OUT_DIR, "figure5_15_cv_vs_loso.png"), fig5_15, width = 8, height = 6, dpi = 300)
} else {
  cat("Skipping Figure 5.15 (CV-vs-LOSO comparison): run the fixed 5_6.R first\n")
  cat("so table5_17_cv_summary.csv exists at:", cv_summary_path, "\n")
}

cat("\nSection 5.7 (LOSO) outputs saved to:", OUT_DIR, "\n")

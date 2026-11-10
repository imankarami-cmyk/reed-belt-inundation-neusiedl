# ============================================================
# SECTION 5.5 — FULL-DATASET GAUGE VALIDATION: TABLES + PLOTS (REVISED)
# ------------------------------------------------------------
# FIXES APPLIED:
#
#  1. The original script scored the ENTIRE cleaned dataset
#     (pix_raw), including the ~79,736 rows used to fit the
#     model, and called the result "independent gauge validation."
#     Those rows are not independent of the model - it saw them
#     during training. This version tags every row with whether
#     it was part of the decile-labeled training set for its
#     model, and reports the validation metrics BOTH with those
#     rows excluded (primary, table5_14/15) and included (for
#     transparency/comparison, table5_14b/15b) so you can see how
#     much the numbers move.
#
#  2. Added a Precision-Recall curve and AUPRC (via the PRROC
#     package) alongside the ROC curve. Given the class imbalance
#     in the gauge-validated ground truth (~13% flooded), AUPRC is
#     a more informative summary of real-world performance than
#     ROC-AUC alone, and reviewers in this field will expect it.
#
#  Models are still fit on the FULL decile-labeled set (this is
#  correct for the final/deployed model - you want to use all
#  available labeled data for the model you actually ship). The
#  fix is in what counts as "validation," not in how the deployed
#  model itself is fit.
# ============================================================

library(tidyverse)
library(pROC)
library(PRROC)

BASE      <- "<PROJECT_ROOT>/result/writting/new/Data/"
GAUGE     <- "<PROJECT_ROOT>/Data/"
IN_CSV    <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
GAUGE_CSV <- paste0(GAUGE, "all_stations_1984_2025.csv")  # adjust path if stored elsewhere
OUT_DIR   <- paste0(BASE, "5.5/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix_raw <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI)) %>%
  mutate(.uid = paste(pixel_id, date, sep = "_"))  # NEW: unique row key, used for leakage exclusion

gauge <- read_csv(GAUGE_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), level_mua = as.numeric(level_mua),
         station = as.character(station)) %>%
  filter(!is.na(date), !is.na(level_mua)) %>%
  distinct(date, station, .keep_all = TRUE)

# ------------------------------------------------------------
# 1. RE-FIT MODELS ON THE FULL DECILE-LABELED SET
#    (this is the deployed/final model - unchanged from original)
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

fit_list <- list()
train_uids <- list()  # NEW: track exactly which rows were used to fit each model

for (m in names(model_defs)) {
  dat <- model_defs[[m]]$labeled %>% filter(!is.na(TCL_label_new))
  fit_list[[m]] <- glm(model_defs[[m]]$formula, data = dat, family = binomial(link = "logit"))
  train_uids[[m]] <- dat$.uid   # every row used in fitting this model
}

# ------------------------------------------------------------
# 2. SCORE FULL DATASET + JOIN GAUGE GROUND TRUTH
# ------------------------------------------------------------

score_full <- function(m) {
  fit <- fit_list[[m]]
  full_df <- pix_raw
  full_df$pred_prob <- predict(fit, newdata = full_df, type = "response")
  full_df %>%
    transmute(model = m, pixel_id, date, elevation, station, pred_prob, .uid,
              used_in_training = .uid %in% train_uids[[m]])   # NEW flag
}

scored_all <- bind_rows(lapply(names(model_defs), score_full))

full_validation_all <- scored_all %>%
  left_join(gauge %>% select(date, station, level_mua), by = c("date", "station")) %>%
  filter(!is.na(level_mua), !is.na(elevation)) %>%
  mutate(
    water_present = if_else(level_mua >= elevation, 1L, 0L),
    model_label = sapply(model, function(m) model_defs[[m]]$label)
  )

# PRIMARY validation set: training rows excluded (FIX)
full_validation <- full_validation_all %>% filter(!used_in_training)

# Secondary set kept for direct comparison to the original (unfixed) approach
full_validation_incl_train <- full_validation_all

write_csv(full_validation, paste0(OUT_DIR, "full_dataset_scored_gauge_joined_excl_train.csv"))
write_csv(full_validation_incl_train, paste0(OUT_DIR, "full_dataset_scored_gauge_joined_incl_train.csv"))

cat("Validation set sizes:\n")
full_validation_all %>%
  group_by(model_label) %>%
  summarise(n_total = n(), n_used_in_training = sum(used_in_training),
            n_excluded_primary = n() - sum(used_in_training), .groups = "drop") %>%
  print()
cat("\n")

# ------------------------------------------------------------
# 3. HELPER: compute ROC/threshold/confusion/PR metrics for a given
#    validation data frame (reused for both excl/incl comparisons)
# ------------------------------------------------------------

compute_metrics <- function(val_df, tag) {
  roc_list <- list(); thr_list <- list(); conf_list <- list(); pr_list <- list()

  for (m in names(model_defs)) {
    sub <- val_df %>% filter(model == m)
    if (nrow(sub) < 50 || n_distinct(sub$water_present) < 2) next

    roc_o <- roc(sub$water_present, sub$pred_prob, quiet = TRUE)
    roc_list[[m]] <- tibble(model = m, tag = tag,
                             fpr = 1 - roc_o$specificities, tpr = roc_o$sensitivities)

    youden <- coords(roc_o, "best", best.method = "youden",
                     ret = c("threshold", "sensitivity", "specificity"))
    opt_thr <- as.numeric(youden["threshold"])

    # NEW: Precision-Recall curve / AUPRC
    pr_o <- pr.curve(scores.class0 = sub$pred_prob[sub$water_present == 1],
                      scores.class1 = sub$pred_prob[sub$water_present == 0],
                      curve = TRUE)
    pr_list[[m]] <- tibble(model = m, tag = tag, auprc = pr_o$auc.integral)

    thr_list[[m]] <- tibble(
      model = m, tag = tag, optimal_threshold = round(opt_thr, 4),
      sensitivity_at_optimal = round(as.numeric(youden["sensitivity"]), 4),
      specificity_at_optimal = round(as.numeric(youden["specificity"]), 4),
      auc_full_dataset = round(as.numeric(roc_o$auc), 4),
      auprc = round(pr_o$auc.integral, 4),
      n_rows = nrow(sub), n_pixels = n_distinct(sub$pixel_id),
      n_gauge_wet = sum(sub$water_present == 1), n_gauge_dry = sum(sub$water_present == 0)
    )

    for (thr_name in c("0.5", "optimal")) {
      thr_val <- if (thr_name == "0.5") 0.5 else opt_thr
      pred_bin <- if_else(sub$pred_prob >= thr_val, 1L, 0L)
      tp <- sum(pred_bin==1 & sub$water_present==1); fp <- sum(pred_bin==1 & sub$water_present==0)
      tn <- sum(pred_bin==0 & sub$water_present==0); fn <- sum(pred_bin==0 & sub$water_present==1)
      conf_list[[paste(m, thr_name)]] <- tibble(
        model = m, tag = tag, threshold_type = thr_name, threshold_value = round(thr_val, 4),
        TP=tp, FP=fp, TN=tn, FN=fn,
        accuracy = round((tp+tn)/(tp+fp+tn+fn), 4),
        sensitivity = round(tp/(tp+fn), 4),
        specificity = round(tn/(tn+fp), 4),
        precision = round(tp/(tp+fp), 4),
        f1_score = round(2*tp/(2*tp+fp+fn), 4)
      )
    }
  }
  list(roc = bind_rows(roc_list), thr = bind_rows(thr_list),
       conf = bind_rows(conf_list), pr = bind_rows(pr_list))
}

res_excl <- compute_metrics(full_validation, "excl_train (primary)")
res_incl <- compute_metrics(full_validation_incl_train, "incl_train (reference only)")

table5_14_thresholds <- bind_rows(res_excl$thr, res_incl$thr) %>%
  mutate(model_label = sapply(model, function(m) model_defs[[m]]$label))
table5_15_confusion <- bind_rows(res_excl$conf, res_incl$conf) %>%
  mutate(model_label = sapply(model, function(m) model_defs[[m]]$label))

write_csv(table5_14_thresholds, paste0(OUT_DIR, "table5_14_optimal_thresholds.csv"))
write_csv(table5_15_confusion, paste0(OUT_DIR, "table5_15_confusion_matrix.csv"))
print(table5_14_thresholds)
print(table5_15_confusion)

cat("\nReport the 'excl_train (primary)' rows as Table 5.14/5.15 in the thesis.\n")
cat("The 'incl_train (reference only)' rows show how much the original\n")
cat("(unfixed) approach was inflated by resubstitution — useful as a one-line\n")
cat("transparency note in the text, not as the headline result.\n\n")

# ------------------------------------------------------------
# 4. FIGURE 5.9 — ROC OVERLAY (primary, excl_train)
# ------------------------------------------------------------

model_colors <- c("NDWI + EVI" = "#e07b39", "mNDWI + EVI" = "#01696f")

roc_all <- res_excl$roc %>%
  left_join(table5_14_thresholds %>% filter(tag == "excl_train (primary)") %>%
              select(model, auc_full_dataset, model_label), by = "model") %>%
  mutate(label = sprintf("%s (AUC = %.3f)", model_label, auc_full_dataset))

fig5_9 <- ggplot(roc_all, aes(x = fpr, y = tpr, color = label)) +
  geom_line(linewidth = 1) +
  geom_abline(linetype = "dashed", color = "grey60") +
  scale_color_manual(values = setNames(model_colors[unique(roc_all$model_label)], unique(roc_all$label))) +
  labs(title = "Full-dataset ROC curves (held-out rows only): NDWI+EVI vs mNDWI+EVI",
       x = "False positive rate", y = "True positive rate", color = "Model") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_9_roc_overlay.png"), fig5_9, width = 7, height = 6.5, dpi = 300)

# ------------------------------------------------------------
# 5. FIGURE 5.9b — PRECISION-RECALL OVERLAY (NEW)
# ------------------------------------------------------------

pr_curve_points <- function(m, tagname) {
  sub <- full_validation %>% filter(model == m)
  pr_o <- pr.curve(scores.class0 = sub$pred_prob[sub$water_present == 1],
                    scores.class1 = sub$pred_prob[sub$water_present == 0],
                    curve = TRUE)
  as_tibble(pr_o$curve[, 1:2]) %>%
    rename(recall = V1, precision = V2) %>%
    mutate(model = m, tag = tagname)
}

pr_points_all <- bind_rows(lapply(names(model_defs), pr_curve_points, tagname = "excl_train")) %>%
  left_join(res_excl$pr, by = "model") %>%
  mutate(model_label = sapply(model, function(m) model_defs[[m]]$label),
         label = sprintf("%s (AUPRC = %.3f)", model_label, auprc))

fig5_9b <- ggplot(pr_points_all, aes(x = recall, y = precision, color = label)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = mean(full_validation$water_present), linetype = "dashed", color = "grey60") +
  scale_color_manual(values = setNames(model_colors[unique(pr_points_all$model_label)], unique(pr_points_all$label))) +
  labs(title = "Precision-Recall curves (held-out rows only): NDWI+EVI vs mNDWI+EVI",
       subtitle = "Dashed line = baseline precision at the observed flood prevalence",
       x = "Recall (sensitivity)", y = "Precision", color = "Model") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_9b_pr_overlay.png"), fig5_9b, width = 7, height = 6.5, dpi = 300)

# ------------------------------------------------------------
# 6. FIGURE 5.10 — CONFUSION MATRIX HEATMAP (optimal threshold, primary)
# ------------------------------------------------------------

confusion_plot_df <- table5_15_confusion %>%
  filter(threshold_type == "optimal", tag == "excl_train (primary)") %>%
  select(model_label, TP, FP, TN, FN) %>%
  pivot_longer(-model_label, names_to = "cell", values_to = "count") %>%
  mutate(cell = factor(cell, levels = c("TP","FN","FP","TN")))

fig5_10 <- ggplot(confusion_plot_df, aes(x = cell, y = model_label, fill = count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = format(count, big.mark = ",")), color = "white", size = 4) +
  scale_fill_gradient(low = "#cedcd8", high = "#01696f") +
  labs(title = "Confusion matrix at optimal threshold (held-out rows only)",
       x = NULL, y = NULL, fill = "Count") +
  theme_minimal(base_size = 13)
ggsave(paste0(OUT_DIR, "figure5_10_confusion_heatmap.png"), fig5_10, width = 7, height = 4.5, dpi = 300)

# ------------------------------------------------------------
# 7. FIGURE 5.11 — PREDICTED PROBABILITY VS GAUGE STATE (BOXPLOT, primary)
# ------------------------------------------------------------

fig5_11 <- full_validation %>%
  mutate(gauge_state = if_else(water_present == 1, "wet (gauge)", "dry (gauge)")) %>%
  ggplot(aes(x = model_label, y = pred_prob, fill = gauge_state)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) +
  scale_fill_manual(values = c("dry (gauge)" = "#c98a3a", "wet (gauge)" = "#01696f")) +
  labs(title = "Predicted flood probability vs. gauge-observed state (held-out rows only)",
       x = "Model", y = "Predicted probability of flooding", fill = "Gauge state") +
  theme_minimal(base_size = 13)
ggsave(paste0(OUT_DIR, "figure5_11_probability_boxplot.png"), fig5_11, width = 8, height = 6, dpi = 300)

cat("Section 5.5 outputs saved to:", OUT_DIR, "\n")

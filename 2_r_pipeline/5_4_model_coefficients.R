# ============================================================
# SECTION 5.4 — MODEL COEFFICIENTS: TABLES + PLOTS  (REVISED)
# ------------------------------------------------------------
# FIXES APPLIED (see notes below):
#
#  1. The original script fit the GLM on the full decile-labeled
#     set and then evaluated AUC/accuracy/AIC/R^2 by predicting
#     on that SAME set — i.e. it reported in-sample fit and
#     called it "training performance." That is not a train/test
#     comparison at all, and because decile labels are the top/
#     bottom 10% of the very index being modeled, the problem is
#     close to trivially separable by construction — which is
#     why AUC/R^2 came out so high and the coefficients so large.
#
#  2. This version splits the decile-labeled set into a genuine,
#     stratified 70% train / 30% held-out partition BEFORE
#     fitting. Table 5.13 now reports BOTH:
#       - "in-sample (train)"  -> kept for transparency, but
#         explicitly flagged as a fit diagnostic, not a
#         performance estimate
#       - "held-out (same decile set, 30%)" -> a fairer, still
#         same-distribution, training-stage performance number
#
#  3. A quasi-complete-separation diagnostic (% of fitted
#     probabilities within 0.01 of 0 or 1) is added, so you can
#     see directly how much of the in-sample AUC is an artifact
#     of the decile construction rather than genuine signal.
#
#  Coefficients themselves (Table 5.12) are still reported from
#  the TRAIN fit only — that part of the original design was
#  fine, it just needs to be paired with an honest performance
#  estimate rather than a resubstituted one.
# ============================================================

library(tidyverse)
library(pROC)

BASE   <- "<PROJECT_ROOT>/result/writting/new/Data/"
IN_CSV <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
OUT_DIR <- paste0(BASE, "5.4/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(42)  # fixed seed so the split is reproducible

pix_raw <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI))

# ------------------------------------------------------------
# 1. DECILE LABELING (same as Section 5.3)
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
    ) %>%
    ungroup()
}

ndwi_labeled  <- label_decile(pix_raw, "NDWI")
mndwi_labeled <- label_decile(pix_raw, "mNDWI")

model_defs <- list(
  NDWI_EVI  = list(formula = TCL_label_new ~ NDWI + EVI,  labeled = ndwi_labeled,  label = "NDWI + EVI"),
  mNDWI_EVI = list(formula = TCL_label_new ~ mNDWI + EVI, labeled = mndwi_labeled, label = "mNDWI + EVI")
)

# ------------------------------------------------------------
# 2. STRATIFIED 70/30 TRAIN / HELD-OUT SPLIT (FIX)
# ------------------------------------------------------------

split_stratified <- function(dat, p = 0.7) {
  dat <- dat %>% mutate(.row_id = row_number())
  train_idx <- dat %>%
    group_by(TCL_label_new) %>%
    slice_sample(prop = p) %>%
    pull(.row_id)
  list(
    train   = dat %>% filter(.row_id %in% train_idx)  %>% select(-.row_id),
    heldout = dat %>% filter(!.row_id %in% train_idx) %>% select(-.row_id)
  )
}

fit_list <- list(); coef_list <- list(); perf_list <- list()

for (m in names(model_defs)) {
  spec <- model_defs[[m]]
  dat  <- spec$labeled %>% filter(!is.na(TCL_label_new))

  parts <- split_stratified(dat, p = 0.7)
  train <- parts$train
  held  <- parts$heldout

  fit <- glm(spec$formula, data = train, family = binomial(link = "logit"))
  fit_list[[m]] <- fit

  # --- in-sample (train) diagnostics ---
  prob_train   <- predict(fit, train, type = "response")
  roc_train    <- roc(train$TCL_label_new, prob_train, quiet = TRUE)
  pred05_train <- if_else(prob_train >= 0.5, 1L, 0L)
  pct_extreme_train <- mean(prob_train < 0.01 | prob_train > 0.99) * 100

  # --- held-out (same-distribution decile set, NOT used in fitting) ---
  prob_held   <- predict(fit, held, type = "response")
  roc_held    <- roc(held$TCL_label_new, prob_held, quiet = TRUE)
  pred05_held <- if_else(prob_held >= 0.5, 1L, 0L)

  ci <- confint(fit)
  coef_list[[m]] <- summary(fit)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term") %>%
    rename(estimate = Estimate, std_error = `Std. Error`,
           z_value = `z value`, p_value = `Pr(>|z|)`) %>%
    mutate(
      odds_ratio = exp(estimate),
      ci_lower = exp(ci[, 1]),
      ci_upper = exp(ci[, 2]),
      model = m,
      model_label = spec$label,
      significance = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      )
    )

  perf_list[[paste0(m, "_train")]] <- tibble(
    model = m, model_label = spec$label, partition = "in-sample (train, 70%) - diagnostic only",
    n_rows = nrow(train), n_pixels = n_distinct(train$pixel_id),
    auc = as.numeric(roc_train$auc),
    accuracy_0_5 = mean(train$TCL_label_new == pred05_train),
    aic = AIC(fit),
    mcfadden_r2 = 1 - fit$deviance / fit$null.deviance,
    pct_extreme_prob = round(pct_extreme_train, 2)
  )

  perf_list[[paste0(m, "_held")]] <- tibble(
    model = m, model_label = spec$label, partition = "held-out (same decile set, 30%) - primary",
    n_rows = nrow(held), n_pixels = n_distinct(held$pixel_id),
    auc = as.numeric(roc_held$auc),
    accuracy_0_5 = mean(held$TCL_label_new == pred05_held),
    aic = NA_real_,
    mcfadden_r2 = NA_real_,
    pct_extreme_prob = NA_real_
  )
}

coef_table <- bind_rows(coef_list) %>%
  mutate(across(where(is.numeric), ~round(., 4)))
perf_summary <- bind_rows(perf_list) %>%
  mutate(across(where(is.numeric) & !c(n_rows, n_pixels), ~round(., 4)))

write_csv(coef_table, paste0(OUT_DIR, "table5_12_model_coefficients.csv"))
write_csv(perf_summary, paste0(OUT_DIR, "table5_13_training_vs_heldout_summary.csv"))

print(coef_table)
print(perf_summary)

cat("\nIMPORTANT — how to read table5_13:\n")
cat("The 'in-sample (train)' rows show fit to the decile-construction problem\n")
cat("and will look inflated by design (decile labels are the extremes of the\n")
cat("very index being modeled) - report these only as a diagnostic, alongside\n")
cat("pct_extreme_prob, which flags possible quasi-complete separation.\n")
cat("The 'held-out' rows are the fairer training-stage comparison to quote as\n")
cat("Table 5.13 in the thesis. Section 5.5/5.6 (gauge validation) remain the\n")
cat("real-world test.\n\n")

# ------------------------------------------------------------
# 3. FIGURE 5.7 / 5.8 — Coefficient / odds-ratio plots (from train fit)
# ------------------------------------------------------------

coef_plot_df <- coef_table %>% filter(term != "(Intercept)")
model_colors <- c("NDWI + EVI" = "#e07b39", "mNDWI + EVI" = "#01696f")

fig5_7 <- ggplot(coef_plot_df, aes(x = term, y = estimate, color = model_label)) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_errorbar(aes(ymin = log(ci_lower), ymax = log(ci_upper)),
                position = position_dodge(width = 0.4), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = model_colors) +
  labs(title = "Logistic regression coefficients with 95% confidence intervals",
       x = NULL, y = "Coefficient estimate (log-odds)", color = "Model") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_7_coefficient_plot.png"), fig5_7, width = 7, height = 6, dpi = 300)

fig5_8 <- ggplot(coef_plot_df, aes(x = term, y = odds_ratio, color = model_label)) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                position = position_dodge(width = 0.4), width = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_y_log10() +
  scale_color_manual(values = model_colors) +
  labs(title = "Odds ratios with 95% confidence intervals",
       x = NULL, y = "Odds ratio (log scale)", color = "Model") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_8_odds_ratio_plot.png"), fig5_8, width = 7, height = 6, dpi = 300)

cat("Section 5.4 outputs saved to:", OUT_DIR, "\n")

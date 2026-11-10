# ============================================================
# SECTION 5.3 — DECILE LABELING: TABLES + PLOTS
# ------------------------------------------------------------
# NOTE: No change from the original. The per-pixel decile rule
# itself is a deliberate design choice, not a bug — the issues
# in this pipeline arise later, in how the resulting labels are
# used for "training performance" (5.4) and validation (5.5/5.6),
# which are fixed in those scripts.
# ============================================================

library(tidyverse)

BASE   <- "<PROJECT_ROOT>/result/writting/5.2/"
IN_CSV <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
OUT_DIR <- paste0(BASE, "5.3/")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix_raw <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI))

# ------------------------------------------------------------
# 1. PER-PIXEL DECILE LABELING (10th / 90th percentile)
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

# ------------------------------------------------------------
# 2. TABLE 5.9 — LABEL COUNT SUMMARY (flooded / dry / ambiguous)
# ------------------------------------------------------------

summarize_labels <- function(df, index_name) {
  n_flood <- sum(df$TCL_label_new == 1, na.rm = TRUE)
  n_dry   <- sum(df$TCL_label_new == 0, na.rm = TRUE)
  n_ambig <- sum(is.na(df$TCL_label_new))
  n_total <- nrow(df)
  tibble(
    index = index_name,
    n_total = n_total,
    n_flooded = n_flood,
    n_dry = n_dry,
    n_ambiguous = n_ambig,
    pct_flooded = round(100 * n_flood / (n_flood + n_dry), 1),
    pct_dry = round(100 * n_dry / (n_flood + n_dry), 1),
    pct_ambiguous = round(100 * n_ambig / n_total, 1)
  )
}

table5_9_label_counts <- bind_rows(
  summarize_labels(ndwi_labeled, "NDWI"),
  summarize_labels(mndwi_labeled, "mNDWI")
)
write_csv(table5_9_label_counts, paste0(OUT_DIR, "table5_9_label_counts.csv"))
print(table5_9_label_counts)

# ------------------------------------------------------------
# 3. TABLE 5.10 — POOLED (GLOBAL) DECILE THRESHOLDS FOR REFERENCE
# ------------------------------------------------------------

table5_10_pooled_thresholds <- tibble(
  index = c("NDWI", "mNDWI"),
  p10_pooled = c(quantile(pix_raw$NDWI, 0.10, na.rm = TRUE),
                 quantile(pix_raw$mNDWI, 0.10, na.rm = TRUE)),
  p90_pooled = c(quantile(pix_raw$NDWI, 0.90, na.rm = TRUE),
                 quantile(pix_raw$mNDWI, 0.90, na.rm = TRUE))
) %>% mutate(across(where(is.numeric), ~round(., 3)))
write_csv(table5_10_pooled_thresholds, paste0(OUT_DIR, "table5_10_pooled_thresholds.csv"))
print(table5_10_pooled_thresholds)

# ------------------------------------------------------------
# 4. TABLE 5.11 — PER-PIXEL THRESHOLD VARIABILITY (spread across pixels)
# ------------------------------------------------------------

pixel_thresholds <- pix_raw %>%
  group_by(pixel_id) %>%
  summarise(
    ndwi_p10 = quantile(NDWI, 0.10, na.rm = TRUE),
    ndwi_p90 = quantile(NDWI, 0.90, na.rm = TRUE),
    mndwi_p10 = quantile(mNDWI, 0.10, na.rm = TRUE),
    mndwi_p90 = quantile(mNDWI, 0.90, na.rm = TRUE),
    .groups = "drop"
  )

table5_11_threshold_variability <- pixel_thresholds %>%
  summarise(
    ndwi_p10_mean = mean(ndwi_p10), ndwi_p10_sd = sd(ndwi_p10),
    ndwi_p90_mean = mean(ndwi_p90), ndwi_p90_sd = sd(ndwi_p90),
    mndwi_p10_mean = mean(mndwi_p10), mndwi_p10_sd = sd(mndwi_p10),
    mndwi_p90_mean = mean(mndwi_p90), mndwi_p90_sd = sd(mndwi_p90)
  ) %>% mutate(across(everything(), ~round(., 3)))
write_csv(table5_11_threshold_variability, paste0(OUT_DIR, "table5_11_threshold_variability.csv"))
print(table5_11_threshold_variability)

# ------------------------------------------------------------
# 5. FIGURE 5.5 — LABEL DISTRIBUTION BAR CHART (NDWI vs mNDWI)
# ------------------------------------------------------------

label_plot_df <- table5_9_label_counts %>%
  select(index, n_flooded, n_dry, n_ambiguous) %>%
  pivot_longer(-index, names_to = "label", values_to = "count") %>%
  mutate(label = recode(label, n_flooded = "Flooded", n_dry = "Dry", n_ambiguous = "Ambiguous"),
         label = factor(label, levels = c("Flooded", "Dry", "Ambiguous")))

fig5_5 <- ggplot(label_plot_df, aes(x = index, y = count, fill = label)) +
  geom_col(position = "stack", width = 0.6) +
  geom_text(aes(label = format(count, big.mark = ",")),
            position = position_stack(vjust = 0.5), size = 3.2, color = "white") +
  scale_fill_manual(values = c("Flooded" = "#01696f", "Dry" = "#e07b39", "Ambiguous" = "#b0b0b0")) +
  labs(title = "Decile-based wet/dry label distribution",
       x = NULL, y = "Number of observations", fill = "Label") +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_5_label_distribution.png"), fig5_5, width = 7, height = 6, dpi = 300)

# ------------------------------------------------------------
# 6. FIGURE 5.6 — PER-PIXEL THRESHOLD SPREAD (density of p10/p90 across pixels)
# ------------------------------------------------------------

thr_long <- pixel_thresholds %>%
  select(ndwi_p10, ndwi_p90, mndwi_p10, mndwi_p90) %>%
  pivot_longer(everything(), names_to = "threshold", values_to = "value") %>%
  mutate(
    index = if_else(str_detect(threshold, "ndwi_"), "NDWI", "mNDWI"),
    percentile = if_else(str_detect(threshold, "p10"), "10th percentile", "90th percentile")
  )

fig5_6 <- ggplot(thr_long, aes(x = value, fill = percentile)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~index, scales = "free_x") +
  scale_fill_manual(values = c("10th percentile" = "#e07b39", "90th percentile" = "#01696f")) +
  labs(title = "Distribution of per-pixel decile thresholds",
       x = "Threshold value", y = "Density", fill = NULL) +
  theme_minimal(base_size = 13) + theme(legend.position = "bottom")
ggsave(paste0(OUT_DIR, "figure5_6_threshold_variability.png"), fig5_6, width = 9, height = 5, dpi = 300)

cat("Section 5.3 outputs saved to:", OUT_DIR, "\n")

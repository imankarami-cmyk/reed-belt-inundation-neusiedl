library(tidyverse)
library(lubridate)

BASE     <- "<PROJECT_ROOT>/Data/"
PIX_PATH <- paste0(BASE, "FLATS_newpixels_1984_2025.csv")
OUT_DIR  <- "<PROJECT_ROOT>/result/writting/5.2/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# LOAD + INITIAL FILTER
# ============================================================

pix_raw <- read_csv(PIX_PATH, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))) %>%
  filter(is.finite(NDWI), is.finite(mNDWI), is.finite(EVI))

cat("══════════════════════════════════════════════════════\n")
cat("STEP 1: Raw data loaded\n")
cat("══════════════════════════════════════════════════════\n")
cat("Total observations loaded:", nrow(pix_raw), "\n")
cat("Unique pixels:", n_distinct(pix_raw$pixel_id), "\n")
cat("Sensors present:", paste(unique(pix_raw$sensor), collapse = ", "), "\n\n")

# ============================================================
# STEP 2 — PHYSICAL-RANGE FILTER (remove index computation artifacts)
# ============================================================
# NDWI and mNDWI are mathematically bounded to [-1, 1] by their ratio
# definitions. Values outside this range indicate a near-zero denominator
# (uncorrected reflectance artifacts, e.g. residual cloud/shadow
# contamination not caught by QA_PIXEL masking). EVI is retained within
# a slightly wider operational bound following Huete et al. (2002),
# since dense canopy can occasionally push EVI marginally above 1.

n_before <- nrow(pix_raw)

pix_flagged <- pix_raw %>%
  mutate(
    invalid_NDWI  = NDWI  < -1   | NDWI  > 1,
    invalid_mNDWI = mNDWI < -1   | mNDWI > 1,
    invalid_EVI   = EVI   < -1.5 | EVI   > 1.5,
    invalid_any   = invalid_NDWI | invalid_mNDWI | invalid_EVI
  )

exclusion_summary <- pix_flagged %>%
  summarise(
    n_total           = n(),
    n_invalid_NDWI     = sum(invalid_NDWI),
    n_invalid_mNDWI    = sum(invalid_mNDWI),
    n_invalid_EVI      = sum(invalid_EVI),
    n_invalid_any      = sum(invalid_any),
    pct_excluded       = round(100 * sum(invalid_any) / n(), 2)
  )
write_csv(exclusion_summary, paste0(OUT_DIR, "table5_3b_outlier_exclusion_summary.csv"))
print(exclusion_summary)

exclusion_by_sensor <- pix_flagged %>%
  group_by(sensor) %>%
  summarise(n_total = n(), n_excluded = sum(invalid_any),
            pct_excluded = round(100 * n_excluded / n_total, 2), .groups = "drop")
write_csv(exclusion_by_sensor, paste0(OUT_DIR, "table5_3c_outlier_exclusion_by_sensor.csv"))
print(exclusion_by_sensor)

exclusion_by_station <- pix_flagged %>%
  group_by(station) %>%
  summarise(n_total = n(), n_excluded = sum(invalid_any),
            pct_excluded = round(100 * n_excluded / n_total, 2), .groups = "drop")
write_csv(exclusion_by_station, paste0(OUT_DIR, "table5_3d_outlier_exclusion_by_station.csv"))
print(exclusion_by_station)

pix <- pix_flagged %>% filter(!invalid_any) %>%
  select(-invalid_NDWI, -invalid_mNDWI, -invalid_EVI, -invalid_any)

cat(sprintf("\nRetained %d of %d observations (%.2f%% excluded) after physical-range filtering.\n\n",
            nrow(pix), n_before, 100 * (n_before - nrow(pix)) / n_before))

# ============================================================
# TABLE 5.4 — Overall predictor summary statistics (pooled, cleaned)
# ============================================================

overall_predictor_stats <- pix %>%
  summarise(
    n_obs = n(),
    across(c(NDWI, mNDWI, EVI),
           list(mean = ~mean(.x, na.rm = TRUE),
                sd   = ~sd(.x, na.rm = TRUE),
                min  = ~min(.x, na.rm = TRUE),
                p05  = ~quantile(.x, 0.05, na.rm = TRUE),
                p25  = ~quantile(.x, 0.25, na.rm = TRUE),
                median = ~median(.x, na.rm = TRUE),
                p75  = ~quantile(.x, 0.75, na.rm = TRUE),
                p95  = ~quantile(.x, 0.95, na.rm = TRUE),
                max  = ~max(.x, na.rm = TRUE)),
           .names = "{.col}_{.fn}")
  )
write_csv(overall_predictor_stats, paste0(OUT_DIR, "table5_4_overall_predictor_stats_clean.csv"))
print(overall_predictor_stats)

# ============================================================
# TABLE 5.5 — Per-station predictor summary (cleaned)
# ============================================================

station_predictor_stats <- pix %>%
  group_by(station) %>%
  summarise(
    n_obs = n(),
    NDWI_mean = round(mean(NDWI), 3), NDWI_sd = round(sd(NDWI), 3),
    mNDWI_mean = round(mean(mNDWI), 3), mNDWI_sd = round(sd(mNDWI), 3),
    EVI_mean = round(mean(EVI), 3), EVI_sd = round(sd(EVI), 3),
    .groups = "drop"
  )
write_csv(station_predictor_stats, paste0(OUT_DIR, "table5_5_station_predictor_stats_clean.csv"))
print(station_predictor_stats)

# ============================================================
# TABLE 5.6 — Per-sensor predictor summary (cleaned, cross-sensor check)
# ============================================================

sensor_predictor_stats <- pix %>%
  group_by(sensor) %>%
  summarise(
    n_obs = n(),
    date_min = min(date), date_max = max(date),
    NDWI_mean = round(mean(NDWI), 3), NDWI_sd = round(sd(NDWI), 3),
    mNDWI_mean = round(mean(mNDWI), 3), mNDWI_sd = round(sd(mNDWI), 3),
    EVI_mean = round(mean(EVI), 3), EVI_sd = round(sd(EVI), 3),
    .groups = "drop"
  )
write_csv(sensor_predictor_stats, paste0(OUT_DIR, "table5_6_sensor_predictor_stats_clean.csv"))
print(sensor_predictor_stats)

# ============================================================
# TABLE 5.7 — Predictor correlation matrix (cleaned)
# ============================================================

cor_matrix <- pix %>%
  select(NDWI, mNDWI, EVI, NDVI) %>%
  cor(use = "complete.obs") %>%
  round(3)
write.csv(cor_matrix, paste0(OUT_DIR, "table5_7_predictor_correlation_matrix_clean.csv"))
print(cor_matrix)

# ============================================================
# FIGURE 5.3 — Density distributions of NDWI, mNDWI, EVI (pooled, cleaned)
# ============================================================

density_long <- pix %>%
  select(NDWI, mNDWI, EVI) %>%
  pivot_longer(everything(), names_to = "index", values_to = "value")

density_plot <- ggplot(density_long, aes(x = value, fill = index)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~index, scales = "free", ncol = 3) +
  scale_fill_manual(values = c(NDWI = "#c98a3a", mNDWI = "#01696f", EVI = "#4472c4")) +
  labs(title = "Pooled distribution of NDWI, mNDWI, and EVI (1984-2025)",
       x = "Index value", y = "Density") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(paste0(OUT_DIR, "figure5_3_index_density_pooled_clean.png"), density_plot,
       width = 10, height = 4.5, dpi = 300)

# ============================================================
# FIGURE 5.4 — NDWI and mNDWI with decile markers (cleaned)
# ============================================================

decile_lines <- pix %>%
  summarise(
    NDWI_p10 = quantile(NDWI, 0.10), NDWI_p90 = quantile(NDWI, 0.90),
    mNDWI_p10 = quantile(mNDWI, 0.10), mNDWI_p90 = quantile(mNDWI, 0.90)
  )
write_csv(decile_lines, paste0(OUT_DIR, "table5_8_pooled_decile_thresholds_clean.csv"))

ndwi_plot <- ggplot(pix, aes(x = NDWI)) +
  geom_histogram(bins = 100, fill = "#c98a3a", alpha = 0.8) +
  geom_vline(xintercept = c(decile_lines$NDWI_p10, decile_lines$NDWI_p90),
             linetype = "dashed", color = "black") +
  labs(title = "NDWI distribution with 10th/90th percentile decile thresholds",
       x = "NDWI", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(paste0(OUT_DIR, "figure5_4a_ndwi_decile_thresholds_clean.png"), ndwi_plot,
       width = 7, height = 5, dpi = 300)

mndwi_plot <- ggplot(pix, aes(x = mNDWI)) +
  geom_histogram(bins = 100, fill = "#01696f", alpha = 0.8) +
  geom_vline(xintercept = c(decile_lines$mNDWI_p10, decile_lines$mNDWI_p90),
             linetype = "dashed", color = "black") +
  labs(title = "mNDWI distribution with 10th/90th percentile decile thresholds",
       x = "mNDWI", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(paste0(OUT_DIR, "figure5_4b_mndwi_decile_thresholds_clean.png"), mndwi_plot,
       width = 7, height = 5, dpi = 300)

# ============================================================
# FIGURE 5.5 — NDWI vs mNDWI scatter (cleaned)
# ============================================================

scatter_plot <- pix %>%
  slice_sample(n = min(50000, nrow(pix))) %>%
  ggplot(aes(x = NDWI, y = mNDWI)) +
  geom_point(alpha = 0.15, size = 0.6, color = "#01696f") +
  geom_smooth(method = "lm", color = "#c98a3a", se = FALSE) +
  labs(title = sprintf("NDWI vs. mNDWI relationship (r = %.3f)", cor_matrix["NDWI","mNDWI"]),
       x = "NDWI", y = "mNDWI") +
  theme_minimal(base_size = 12)
ggsave(paste0(OUT_DIR, "figure5_5_ndwi_vs_mndwi_scatter_clean.png"), scatter_plot,
       width = 7, height = 6, dpi = 300)

# ============================================================
# FIGURE 5.6 — EVI vs each water index (cleaned)
# ============================================================

evi_ndwi_plot <- pix %>%
  slice_sample(n = min(50000, nrow(pix))) %>%
  ggplot(aes(x = EVI, y = NDWI)) +
  geom_point(alpha = 0.15, size = 0.6, color = "#c98a3a") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(title = sprintf("EVI vs. NDWI (r = %.3f)", cor_matrix["EVI","NDWI"]),
       x = "EVI", y = "NDWI") +
  theme_minimal(base_size = 12)

evi_mndwi_plot <- pix %>%
  slice_sample(n = min(50000, nrow(pix))) %>%
  ggplot(aes(x = EVI, y = mNDWI)) +
  geom_point(alpha = 0.15, size = 0.6, color = "#01696f") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(title = sprintf("EVI vs. mNDWI (r = %.3f)", cor_matrix["EVI","mNDWI"]),
       x = "EVI", y = "mNDWI") +
  theme_minimal(base_size = 12)

combined_evi <- gridExtra::grid.arrange(evi_ndwi_plot, evi_mndwi_plot, ncol = 2)
ggsave(paste0(OUT_DIR, "figure5_6_evi_vs_water_indices_clean.png"), combined_evi,
       width = 11, height = 5, dpi = 300)

# ============================================================
# FIGURE 5.7 — Boxplots by sensor (cleaned, cross-sensor consistency)
# ============================================================

sensor_box <- pix %>%
  select(sensor, NDWI, mNDWI, EVI) %>%
  pivot_longer(-sensor, names_to = "index", values_to = "value") %>%
  ggplot(aes(x = sensor, y = value, fill = sensor)) +
  geom_boxplot(alpha = 0.8, outlier.alpha = 0.05) +
  facet_wrap(~index, scales = "free_y", ncol = 3) +
  labs(title = "Predictor values by Landsat sensor generation",
       x = "Sensor", y = "Value") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")
ggsave(paste0(OUT_DIR, "figure5_7_predictor_by_sensor_clean.png"), sensor_box,
       width = 11, height = 4.5, dpi = 300)

# ============================================================
# SAVE CLEANED DATASET FOR DOWNSTREAM USE
# ============================================================

write_csv(pix, paste0(OUT_DIR, "FLATS_newpixels_1984_2025_CLEANED.csv"))

cat("\nAll Section 5.2 tables and figures saved to:\n", OUT_DIR, "\n")
cat("Cleaned dataset (for use in Section 4.5 labeling and 5.3 onward) saved as:\n")
cat(paste0(OUT_DIR, "FLATS_newpixels_1984_2025_CLEANED.csv\n"))
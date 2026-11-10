# ============================================================
# CROSS-SENSOR CONTINUITY CHECK
# Lake Neusiedl | 1984-2025
# ------------------------------------------------------------
# PURPOSE: demonstrate that the four Landsat generations merged in
# Section 4.3.3 form one continuous, unbroken record - not four
# separately-behaving datasets stitched together. Two panels:
#   (A) the full 1984-2025 mean EVI record, colored by sensor, to
#       show no visible jump at any sensor transition
#   (B) three zoomed overlap-period panels (L5/L7, L7/L8, L8/L9),
#       each showing both sensors during their simultaneous
#       coverage window, to show they track each other closely
#       when observing the same conditions at the same time.
#
# This directly supports the claim already made in Section 5.2
# ("predictor behaviour remained broadly stable across sensor
# generations"), with a continuity check rather than only a
# pooled boxplot comparison.
# ============================================================

library(tidyverse)
library(patchwork)

BASE    <- "<PROJECT_ROOT>/result/writting/new/5.2/"
IN_CSV  <- paste0(BASE, "FLATS_newpixels_1984_2025_CLEANED.csv")
OUT_DIR <- BASE
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d")))

sensor_colors <- c(L5 = "#e07b39", L7 = "#01696f", L8 = "#8856a7", L9 = "#2c7fb8")

# ------------------------------------------------------------
# 1. AGGREGATE: mean EVI across all 600 pixels, per scene date
#    (one point per Landsat overpass date, not per pixel-date row)
# ------------------------------------------------------------

daily_mean <- pix %>%
  group_by(date, sensor) %>%
  summarise(mean_EVI = mean(EVI, na.rm = TRUE), n_pixels = n(), .groups = "drop")

# ------------------------------------------------------------
# 2. PANEL A - FULL RECORD, ALL FOUR SENSORS
# ------------------------------------------------------------

fig_full <- ggplot(daily_mean, aes(x = date, y = mean_EVI, color = sensor)) +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_smooth(se = FALSE, span = 0.08, linewidth = 0.8) +
  scale_color_manual(values = sensor_colors) +
  labs(title = "Full 1984-2025 mean EVI record, coloured by sensor",
       x = NULL, y = "Mean EVI (all sampled pixels)", color = "Sensor") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

# ------------------------------------------------------------
# 3. PANEL B - THREE OVERLAP-PERIOD ZOOM PANELS
#    Each shows both sensors during their actual simultaneous
#    coverage window, matching the ranges used in the corrected
#    GEE extraction script.
# ------------------------------------------------------------

make_overlap_panel <- function(data, start_date, end_date, sensors_in_window, panel_title) {
  sub <- data %>%
    filter(date >= as.Date(start_date), date <= as.Date(end_date),
           sensor %in% sensors_in_window)

  ggplot(sub, aes(x = date, y = mean_EVI, color = sensor)) +
    geom_line(stat = "smooth", method = "loess", span = 0.15, se = FALSE, linewidth = 0.9) +
    geom_point(size = 0.6, alpha = 0.35) +
    scale_color_manual(values = sensor_colors) +
    labs(title = panel_title, x = NULL, y = "Mean EVI") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

panel_l5_l7 <- make_overlap_panel(daily_mean, "1999-01-01", "2012-12-31", c("L5", "L7"),
                                    "Landsat 5 / 7 overlap (1999-2012)")
panel_l7_l8 <- make_overlap_panel(daily_mean, "2013-01-01", "2022-12-31", c("L7", "L8"),
                                    "Landsat 7 / 8 overlap (2013-2022)")
panel_l8_l9 <- make_overlap_panel(daily_mean, "2021-01-01", "2025-12-31", c("L8", "L9"),
                                    "Landsat 8 / 9 overlap (2021-2025)")

overlap_panels <- panel_l5_l7 + panel_l7_l8 + panel_l8_l9 + plot_layout(ncol = 3)

# ------------------------------------------------------------
# 4. COMBINE AND SAVE
# ------------------------------------------------------------

final_fig <- fig_full / overlap_panels + plot_layout(heights = c(1.3, 1))

ggsave(paste0(OUT_DIR, "figure5_sensor_continuity.png"), final_fig,
       width = 11, height = 8, dpi = 320, bg = "white")

cat("Saved to:", paste0(OUT_DIR, "figure5_sensor_continuity.png"), "\n")

# ------------------------------------------------------------
# 5. QUANTIFY OVERLAP AGREEMENT (optional but strengthens the claim)
#    Mean difference in EVI between the two sensors during each
#    overlap window, matched to the nearest date within 8 days -
#    a numeric complement to the visual check above.
# ------------------------------------------------------------

quantify_overlap <- function(data, s1, s2, start_date, end_date) {
  d1 <- data %>% filter(sensor == s1, date >= as.Date(start_date), date <= as.Date(end_date))
  d2 <- data %>% filter(sensor == s2, date >= as.Date(start_date), date <= as.Date(end_date))

  matched <- d1 %>%
    rowwise() %>%
    mutate(closest_d2_diff = min(abs(as.numeric(d2$date - date)))) %>%
    ungroup() %>%
    filter(closest_d2_diff <= 8)   # within 8 days = effectively the same conditions

  if (nrow(matched) == 0) return(tibble(sensor_pair = paste(s1, s2, sep = "/"), n_matched = 0, mean_abs_diff = NA))

  diffs <- map_dbl(seq_len(nrow(matched)), function(i) {
    ref_date <- matched$date[i]
    match_val <- d2$mean_EVI[which.min(abs(as.numeric(d2$date - ref_date)))]
    matched$mean_EVI[i] - match_val
  })

  tibble(sensor_pair = paste(s1, s2, sep = "/"),
         n_matched = length(diffs),
         mean_diff = round(mean(diffs, na.rm = TRUE), 4),
         mean_abs_diff = round(mean(abs(diffs), na.rm = TRUE), 4))
}

overlap_summary <- bind_rows(
  quantify_overlap(daily_mean, "L5", "L7", "1999-01-01", "2012-12-31"),
  quantify_overlap(daily_mean, "L7", "L8", "2013-01-01", "2022-12-31"),
  quantify_overlap(daily_mean, "L8", "L9", "2021-01-01", "2025-12-31")
)

write_csv(overlap_summary, paste0(OUT_DIR, "table5_sensor_overlap_agreement.csv"))
cat("\nSensor overlap agreement (mean EVI difference within 8-day matched pairs):\n")
print(overlap_summary)

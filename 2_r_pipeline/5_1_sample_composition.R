library(tidyverse)

BASE    <- "<PROJECT_ROOT>/result/1-pixel samplying/"
PIX_PATH <- paste0(BASE, "FLATS_station_buffered_pixels_6stations.csv")
OUT_DIR  <- BASE

pix <- read_csv(PIX_PATH, show_col_types = FALSE)

# ============================================================
# TABLE 5.1 — Per-station pixel counts and buffer radii
# ============================================================

station_summary <- pix %>%
  group_by(station) %>%
  summarise(
    n_pixels          = n(),
    buffer_radius_m    = unique(buffer_m),
    mean_dist_km       = round(mean(distance_to_station_km), 2),
    min_dist_km        = round(min(distance_to_station_km), 2),
    max_dist_km        = round(max(distance_to_station_km), 2),
    min_elevation_m     = round(min(elevation), 3),
    max_elevation_m     = round(max(elevation), 3),
    n_elev_bins_used    = n_distinct(elev_bin),
    .groups = "drop"
  ) %>%
  arrange(station)

write_csv(station_summary, paste0(OUT_DIR, "table5_1_station_summary.csv"))
print(station_summary)

# ============================================================
# TABLE 5.2 — Elevation-bin distribution (lake-wide)
# ============================================================

bin_summary <- pix %>%
  group_by(elev_bin, bin_lower_m, bin_upper_m) %>%
  summarise(n_pixels = n(), .groups = "drop") %>%
  arrange(elev_bin)

write_csv(bin_summary, paste0(OUT_DIR, "table5_2_elevation_bin_distribution.csv"))
print(bin_summary)

# ============================================================
# TABLE 5.3 — Station × elevation-bin cross-tabulation
# ============================================================

station_bin_matrix <- pix %>%
  count(station, elev_bin) %>%
  pivot_wider(names_from = elev_bin, values_from = n, values_fill = 0) %>%
  arrange(station)

write_csv(station_bin_matrix, paste0(OUT_DIR, "table5_3_station_by_bin_matrix.csv"))
print(station_bin_matrix)

# ============================================================
# Overall dataset statistics (for narrative text)
# ============================================================

overall_stats <- tibble(
  total_pixels        = nrow(pix),
  n_stations          = n_distinct(pix$station),
  pixels_per_station   = nrow(pix) / n_distinct(pix$station),
  total_elev_bins_used  = n_distinct(pix$elev_bin),
  min_elevation_lakewide = round(min(pix$elevation), 3),
  max_elevation_lakewide = round(max(pix$elevation), 3),
  mean_elevation_lakewide = round(mean(pix$elevation), 3),
  sd_elevation_lakewide   = round(sd(pix$elevation), 3)
)

write_csv(overall_stats, paste0(OUT_DIR, "table5_0_overall_stats.csv"))
print(overall_stats)

# ============================================================
# FIGURE 5.1 — Elevation-bin histogram, faceted by station
# ============================================================

bin_plot <- pix %>%
  ggplot(aes(x = elev_bin)) +
  geom_bar(fill = "#01696f") +
  facet_wrap(~station, ncol = 3) +
  labs(title = "Sample pixel distribution across elevation bins by station",
       x = "Elevation bin (10 cm intervals, 114.5-117.0 m a.s.l.)",
       y = "Number of pixels") +
  theme_minimal(base_size = 12)
ggsave(paste0(OUT_DIR, "figure5_1_elevation_bin_by_station.png"), bin_plot,
       width = 10, height = 6.5, dpi = 300)

# ============================================================
# FIGURE 5.2 — Distance-to-station boxplot
# ============================================================

dist_plot <- pix %>%
  ggplot(aes(x = station, y = distance_to_station_km, fill = station)) +
  geom_boxplot(alpha = 0.8, show.legend = FALSE) +
  labs(title = "Distribution of sample pixel distances from each gauge station",
       x = "Station", y = "Distance to station (km)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(paste0(OUT_DIR, "figure5_2_distance_boxplot.png"), dist_plot,
       width = 8, height = 6, dpi = 300)

cat("\nAll Section 5.1 tables and figures saved to:\n", OUT_DIR, "\n")
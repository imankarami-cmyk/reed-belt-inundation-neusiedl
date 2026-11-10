# ============================================================
# STATION BUFFER ELEVATION INVESTIGATION
# Lake Neusiedl | DGM_wetland_cliped.tif
# ------------------------------------------------------------
# PURPOSE
# Directly investigates the elevation distribution within each
# station's sampling buffer, using the actual wetland-clipped DEM.
# Reproduces the elevation range / bin-count numbers already in
# Table 5.2, verifies them independently from the raster itself,
# and tests whether enlarging Apetlon's buffer further would
# actually improve elevation-bin coverage (addressing the "why 8
# bins at Apetlon" question directly with data).
# ============================================================

library(terra)
library(tidyverse)

# ------------------------------------------------------------
# 0. CONFIG
# ------------------------------------------------------------

DEM_PATH <- "<PROJECT_ROOT>/Data/DGM_wetland_cliped.tif"
OUT_DIR  <- "<PROJECT_ROOT>/result/writting/new/elevation_investigation/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Target elevation range and bin width used throughout the thesis
BIN_LOW   <- 114.5
BIN_HIGH  <- 117.0
BIN_WIDTH <- 0.10
N_BINS    <- round((BIN_HIGH - BIN_LOW) / BIN_WIDTH)  # should be 25

# Station coordinates (WGS84) and their CURRENT buffer radii
stations_wgs <- tibble::tribble(
  ~station,        ~lon,        ~lat,        ~buffer_m,
  "Breitenbrunn",  16.770256,   47.916203,   2500,
  "Moerbisch",     16.698516,   47.751756,   2000,
  "Neusiedl",      16.835600,   47.929900,   2500,
  "Rust",          16.695700,   47.801700,   2000,
  "Illmitz",       16.752600,   47.770000,   2500,
  "Apetlon",       16.7490326,  47.6814064,  4000
)

# ------------------------------------------------------------
# 1. LOAD DEM AND CHECK CRS
# ------------------------------------------------------------

dem <- rast(DEM_PATH)
cat("DEM loaded. CRS:\n")
print(crs(dem, describe = TRUE))
cat("DEM resolution:", res(dem), "\n")
cat("DEM value range (raw, before any NoData masking):\n")
print(minmax(dem))

# ------------------------------------------------------------
# 2. MASK OUT LIKELY NoData ARTIFACTS
#    Same lesson learned from the GEE zonal check: 0 or extreme
#    negative values are almost certainly NoData, not real
#    elevation, for a lake surface at ~115 m a.s.l.
# ------------------------------------------------------------

dem_masked <- dem
dem_masked[dem_masked <= 0] <- NA

cat("\nDEM value range AFTER masking out <= 0 (likely NoData):\n")
print(minmax(dem_masked))

# ------------------------------------------------------------
# 3. BUILD STATION BUFFERS (project to match DEM CRS)
# ------------------------------------------------------------

stations_pts <- vect(stations_wgs, geom = c("lon", "lat"), crs = "EPSG:4326")
stations_proj <- project(stations_pts, crs(dem))

# Buffer each station individually by its own radius
station_buffers <- vect(lapply(seq_len(nrow(stations_wgs)), function(i) {
  buffer(stations_proj[i], width = stations_wgs$buffer_m[i])
}))
station_buffers$station <- stations_wgs$station
station_buffers$buffer_m <- stations_wgs$buffer_m

# ------------------------------------------------------------
# 4. EXTRACT ELEVATION VALUES WITHIN EACH BUFFER
# ------------------------------------------------------------

extracted <- terra::extract(dem_masked, station_buffers, ID = TRUE)
names(extracted)[2] <- "elevation"

extracted <- extracted %>%
  mutate(station = station_buffers$station[ID]) %>%
  filter(!is.na(elevation))

# ------------------------------------------------------------
# 5. SUMMARY STATISTICS PER STATION
# ------------------------------------------------------------

station_summary <- extracted %>%
  group_by(station) %>%
  summarise(
    n_valid_pixels = n(),
    elev_min = min(elevation),
    elev_max = max(elevation),
    elev_mean = mean(elevation),
    elev_median = median(elevation),
    elev_sd = sd(elevation),
    .groups = "drop"
  ) %>%
  left_join(stations_wgs %>% select(station, buffer_m), by = "station")

cat("\n--- Per-station elevation summary (from actual DEM) ---\n")
print(station_summary)

# ------------------------------------------------------------
# 6. ELEVATION-BIN COVERAGE PER STATION
#    Reproduces the "elevation bins used" column in Table 5.2,
#    computed independently from the raster itself.
# ------------------------------------------------------------

bin_breaks <- seq(BIN_LOW, BIN_HIGH, by = BIN_WIDTH)

extracted <- extracted %>%
  mutate(
    elev_clamped = pmin(pmax(elevation, BIN_LOW), BIN_HIGH - 0.001),
    bin = cut(elev_clamped, breaks = bin_breaks, include.lowest = TRUE, right = FALSE,
              labels = FALSE)
  )

bin_coverage <- extracted %>%
  filter(elevation >= BIN_LOW, elevation < BIN_HIGH) %>%
  group_by(station) %>%
  summarise(n_bins_covered = n_distinct(bin), .groups = "drop")

station_summary <- station_summary %>%
  left_join(bin_coverage, by = "station") %>%
  arrange(desc(n_bins_covered))

cat("\n--- Elevation-bin coverage (out of", N_BINS, "possible bins) ---\n")
print(station_summary %>% select(station, buffer_m, elev_min, elev_max, n_bins_covered))

write_csv(station_summary, paste0(OUT_DIR, "station_elevation_investigation.csv"))

# ------------------------------------------------------------
# 7. APETLON BUFFER-SIZE SENSITIVITY TEST
#    Directly tests whether enlarging Apetlon's buffer beyond its
#    current 4000 m would improve elevation-bin coverage, or
#    whether coverage has already plateaued (i.e. the constraint
#    is genuinely local geography, not an undersized buffer).
# ------------------------------------------------------------

apetlon_pt <- stations_proj[stations_wgs$station == "Apetlon"]
test_radii <- c(2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 7000, 8000)

apetlon_sensitivity <- map_dfr(test_radii, function(r) {
  buf <- buffer(apetlon_pt, width = r)
  vals <- terra::extract(dem_masked, buf)[, 2]
  vals <- vals[!is.na(vals)]
  vals_binned <- vals[vals >= BIN_LOW & vals < BIN_HIGH]
  bins <- cut(pmin(pmax(vals_binned, BIN_LOW), BIN_HIGH - 0.001),
              breaks = bin_breaks, include.lowest = TRUE, right = FALSE, labels = FALSE)
  tibble(
    buffer_radius_m = r,
    n_valid_pixels = length(vals),
    elev_min = if (length(vals) > 0) min(vals) else NA,
    elev_max = if (length(vals) > 0) max(vals) else NA,
    n_bins_covered = n_distinct(bins)
  )
})

cat("\n--- Apetlon buffer-radius sensitivity test ---\n")
cat("If n_bins_covered plateaus as radius increases, the elevation\n")
cat("ceiling near Apetlon is a genuine local constraint, not an\n")
cat("undersized buffer.\n\n")
print(apetlon_sensitivity)

write_csv(apetlon_sensitivity, paste0(OUT_DIR, "apetlon_buffer_sensitivity.csv"))

# ------------------------------------------------------------
# 8. VISUALIZE ELEVATION DISTRIBUTION PER STATION
# ------------------------------------------------------------

density_plot <- ggplot(extracted, aes(x = elevation, fill = station)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = c(BIN_LOW, BIN_HIGH), linetype = "dashed", color = "grey40") +
  facet_wrap(~station, ncol = 3, scales = "free_y") +
  coord_cartesian(xlim = c(114.8, 117.5)) +
  labs(title = "Elevation distribution within each station buffer",
       subtitle = "Dashed lines mark the 114.5-117.0 m target range; x-axis clipped to exclude rare outlier pixels",
       x = "Elevation (m a.s.l.)", y = "Density") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(paste0(OUT_DIR, "station_elevation_distributions.png"), density_plot,
       width = 12, height = 7, dpi = 300, bg = "white")

# Second version: same data, NOT clipped, so the outlier tail (non-wetland
# pixels captured by the buffer) is visible for diagnostic purposes.
density_plot_full <- ggplot(extracted, aes(x = elevation, fill = station)) +
  geom_density(alpha = 0.7) +
  geom_vline(xintercept = c(BIN_LOW, BIN_HIGH), linetype = "dashed", color = "grey40") +
  facet_wrap(~station, ncol = 3, scales = "free") +
  labs(title = "Elevation distribution within each station buffer (full range, unclipped)",
       subtitle = "Shows the full tail, including likely non-wetland pixels at high elevation",
       x = "Elevation (m a.s.l.)", y = "Density") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave(paste0(OUT_DIR, "station_elevation_distributions_full_range.png"), density_plot_full,
       width = 12, height = 7, dpi = 300, bg = "white")

apetlon_plot <- ggplot(apetlon_sensitivity, aes(x = buffer_radius_m, y = n_bins_covered)) +
  geom_line(color = "#01696f", linewidth = 1) +
  geom_point(color = "#01696f", size = 2) +
  geom_hline(yintercept = N_BINS, linetype = "dashed", color = "grey50") +
  labs(title = "Apetlon: elevation-bin coverage vs. buffer radius",
       subtitle = "Dashed line marks the full possible 25 bins",
       x = "Buffer radius (m)", y = "Elevation bins covered") +
  theme_minimal(base_size = 12)

ggsave(paste0(OUT_DIR, "apetlon_buffer_sensitivity.png"), apetlon_plot,
       width = 8, height = 5, dpi = 300, bg = "white")

cat("\nAll outputs saved to:", OUT_DIR, "\n")

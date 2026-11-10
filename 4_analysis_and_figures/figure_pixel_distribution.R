# ============================================================
# FIGURE 5.1 - SAMPLED PIXEL DISTRIBUTION BY ELEVATION BIN, PER STATION
# ------------------------------------------------------------
# This is the figure already referenced in Section 5.1's text:
# "Figure 5.1 shows how sampled pixels are distributed across
# elevation bins for each station." Uses the ACTUAL 600-pixel
# sample carried through the whole pipeline (not the full buffer
# raster from the earlier investigation script) - one row per
# unique pixel_id, since elevation is static per pixel across
# the whole time series.
# ============================================================

library(tidyverse)

IN_CSV  <- "<PROJECT_ROOT>/result/writting/new/5.2/FLATS_newpixels_1984_2025_CLEANED.csv"
OUT_DIR <- "<PROJECT_ROOT>/result/writting/new/elevation_investigation/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BIN_LOW   <- 114.5
BIN_HIGH  <- 117.0
BIN_WIDTH <- 0.10

# ------------------------------------------------------------
# 1. LOAD AND REDUCE TO ONE ROW PER PIXEL
#    Elevation and station are static per pixel_id - collapsing
#    to distinct pixels avoids double-counting the same pixel
#    across its ~500+ observation dates.
# ------------------------------------------------------------

pix_all <- read_csv(IN_CSV, show_col_types = FALSE)

pix_unique <- pix_all %>%
  distinct(pixel_id, station, elevation, elev_bin) %>%
  filter(!is.na(elevation))

cat("Total unique sampled pixels:", nrow(pix_unique), "(expected 600)\n")
cat("Pixels per station:\n")
print(pix_unique %>% count(station))

# ------------------------------------------------------------
# 2. ASSIGN 10 CM ELEVATION BINS
#    Recomputed directly from elevation, rather than trusting the
#    elev_bin column blindly - a useful cross-check that the
#    original bin assignment in the sampling script was correct.
# ------------------------------------------------------------

bin_breaks <- seq(BIN_LOW, BIN_HIGH, by = BIN_WIDTH)
bin_labels <- paste0(round(bin_breaks[-length(bin_breaks)], 1), "-",
                     round(bin_breaks[-1], 1))

pix_unique <- pix_unique %>%
  mutate(
    bin_recomputed = cut(elevation, breaks = bin_breaks, include.lowest = TRUE,
                         right = FALSE, labels = bin_labels)
  )

# Flag any mismatch between the original elev_bin column and the
# freshly recomputed bin, in case of a labeling inconsistency upstream.
mismatch_check <- pix_unique %>%
  filter(!is.na(elev_bin)) %>%
  mutate(bin_matches = as.character(elev_bin) == as.character(bin_recomputed))

n_mismatch <- sum(!mismatch_check$bin_matches, na.rm = TRUE)
cat("\nRows where original elev_bin differs from recomputed bin:", n_mismatch, "\n")
if (n_mismatch > 0) {
  cat("Inspect these rows before trusting the figure - possible upstream labeling issue.\n")
}

# ------------------------------------------------------------
# 3. PIXEL COUNT PER STATION PER BIN
# ------------------------------------------------------------

bin_counts <- pix_unique %>%
  count(station, bin_recomputed, name = "n_pixels") %>%
  complete(station, bin_recomputed, fill = list(n_pixels = 0))

write_csv(bin_counts, paste0(OUT_DIR, "figure5_1_pixel_counts_by_bin_station.csv"))

# ------------------------------------------------------------
# 4. FIGURE 5.1 - FACETED BAR CHART, ONE PANEL PER STATION
# ------------------------------------------------------------

library(RColorBrewer)

station_order <- c("Breitenbrunn", "Neusiedl", "Illmitz", "Moerbisch", "Rust", "Apetlon")

# Professional, colorblind-safe qualitative palette (Dark2, ColorBrewer),
# standard for scientific figures - replaces ggplot's default rainbow fill.
station_palette <- setNames(brewer.pal(6, "Dark2"), station_order)

fig_5_1 <- bin_counts %>%
  mutate(station = factor(station, levels = station_order)) %>%
  ggplot(aes(x = bin_recomputed, y = n_pixels, fill = station)) +
  geom_col(show.legend = FALSE, color = "grey30", linewidth = 0.15) +
  facet_wrap(~station, ncol = 2) +
  scale_fill_manual(values = station_palette) +
  labs(title = "Figure 5.1 - Sampled pixel count by elevation bin, per station",
       x = "Elevation bin (m a.s.l.)", y = "Number of sampled pixels") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = NA),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(OUT_DIR, "Figure5_1_pixel_distribution_by_station.png"), fig_5_1,
       width = 12, height = 9, dpi = 300, bg = "white")

# ------------------------------------------------------------
# 5. LAKE-WIDE VERSION - ALL SIX STATIONS COMBINED
#    Supports the specific sentence in Section 5.1: "Lake-wide,
#    pixel counts are highest between 115.4 and 115.9 m a.s.l.
#    and decline steadily toward the top of the range."
# ------------------------------------------------------------

lakewide_counts <- pix_unique %>%
  count(bin_recomputed, name = "n_pixels") %>%
  complete(bin_recomputed, fill = list(n_pixels = 0))

fig_lakewide <- ggplot(lakewide_counts, aes(x = bin_recomputed, y = n_pixels)) +
  geom_col(fill = "#01696f") +
  labs(title = "Lake-wide sampled pixel count by elevation bin (all six stations combined)",
       x = "Elevation bin (m a.s.l.)", y = "Number of sampled pixels") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))

ggsave(paste0(OUT_DIR, "Figure5_1b_lakewide_pixel_distribution.png"), fig_lakewide,
       width = 10, height = 6, dpi = 300, bg = "white")

cat("\nFigure 5.1 and supporting data saved to:", OUT_DIR, "\n")

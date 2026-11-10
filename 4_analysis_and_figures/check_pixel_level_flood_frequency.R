# ============================================================
# SANITY CHECK: PIXEL-LEVEL FLOOD FREQUENCY BY STATION
# ------------------------------------------------------------
# Recomputes long-term flood frequency per station DIRECTLY from
# the 600-pixel sample, using the same deployed coefficients and
# thresholds as the GEE map export - completely independent of the
# raster/zonal-statistics pipeline. If this matches the GEE-derived
# ranking (Moerbisch/Illmitz high, Neusiedl low), that's strong
# evidence the raster result is real. If it disagrees, the problem
# is in the raster export, not the underlying finding.
# ============================================================

library(tidyverse)

IN_CSV <- "<PROJECT_ROOT>/result/writting/new/5.2/FLATS_newpixels_1984_2025_CLEANED.csv"

pix <- read_csv(IN_CSV, show_col_types = FALSE)

# Deployed model coefficients and thresholds - MUST match GEE_corrected.ipynb exactly.
N_INT <- 3.936; N_NDWI <- 7.957; N_EVI <- -5.566
M_INT <- 1.917; M_MNDWI <- 7.961; M_EVI <- 0.105
NDWI_EVI_THR  <- 0.776
MNDWI_EVI_THR <- 0.561

pix <- pix %>%
  mutate(
    p_ndwi_evi  = 1 - 1 / (1 + exp(N_INT + N_NDWI * NDWI + N_EVI * EVI)),
    p_mndwi_evi = 1 - 1 / (1 + exp(M_INT + M_MNDWI * mNDWI + M_EVI * EVI)),
    flooded_ndwi_evi  = p_ndwi_evi  >= NDWI_EVI_THR,
    flooded_mndwi_evi = p_mndwi_evi >= MNDWI_EVI_THR
  )

station_freq <- pix %>%
  group_by(station) %>%
  summarise(
    n_obs = n(),
    NDWI_EVI_freq_pixel_level  = mean(flooded_ndwi_evi),
    mNDWI_EVI_freq_pixel_level = mean(flooded_mndwi_evi),
    .groups = "drop"
  ) %>%
  arrange(desc(mNDWI_EVI_freq_pixel_level))

cat("Pixel-level flood frequency by station (independent of raster export):\n")
print(station_freq)

cat("\nCompare this ranking directly against table5_11_zonal_flood_frequency_by_station.csv.\n")
cat("If both agree on which stations are highest/lowest, the raster result is confirmed.\n")
cat("If they disagree, the discrepancy points to a problem in the raster/GEE pipeline\n")
cat("(e.g. wetland mask, band order, or elevation clipping), not in the underlying finding.\n")

write_csv(station_freq, "<PROJECT_ROOT>/result/writting/new/5.2/sanity_check_pixel_level_flood_frequency.csv")

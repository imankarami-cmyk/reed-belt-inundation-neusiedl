# ============================================================
# SEICHE ROBUSTNESS CHECK ACROSS THREE WIND STATIONS
# Lake Neusiedl
# ------------------------------------------------------------
# PURPOSE
# Confirms the north-south seiche signal (demonstrated in detail
# from the Neusiedl wind record) is not an artifact of a single
# wind station. Runs the SAME directional test - correlation
# between north-south wind forcing and the Apetlon-minus-Neusiedl
# gauge gap - for all three available wind stations, each over
# its own coverage window, and outputs a single summary table.
#
# This is a robustness supplement to the full Neusiedl analysis
# and Figure 6.1, NOT a replacement. It produces one small table
# for the text, not three separate figures.
#
# WIND DIRECTION ENCODING (dd32, GeoSphere Austria klima-v2-1d,
# DOI: 10.60669/gs6w-jd70): code x 11.25 = degrees.
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 0. CONFIG
# ------------------------------------------------------------

WATER_CSV <- "<PROJECT_ROOT>/Data/all_stations_1984_2025.csv"
WIND_DIR  <- "<PROJECT_ROOT>/Data/"
OUT_DIR   <- "<PROJECT_ROOT>/result/writting/new/seiche/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SOUTH_GAUGE <- "Apetlon"
NORTH_GAUGE <- "Neusiedl"

# The three wind stations and their files. The gauge GAP is always
# Apetlon (south) - Neusiedl (north); only the WIND source changes,
# to test whether each independent wind record reproduces the same
# directional relationship with that fixed north-south gap.
wind_files <- tibble(
  wind_station = c("Neusiedl", "Illmitz", "Podersdorf"),
  file = c("Neusiedl_wind.csv", "Illmitz_wind.csv", "Podersdorf_wind.csv")
)

# ------------------------------------------------------------
# 1. WATER LEVEL GAP (long format: date, level_mua, station, ...)
# ------------------------------------------------------------

water_raw <- read_csv(WATER_CSV, show_col_types = FALSE)

water <- water_raw %>%
  mutate(date = as.Date(date)) %>%
  filter(station %in% c(SOUTH_GAUGE, NORTH_GAUGE)) %>%
  select(date, station, level_mua) %>%
  pivot_wider(names_from = station, values_from = level_mua) %>%
  rename(south = all_of(SOUTH_GAUGE), north = all_of(NORTH_GAUGE)) %>%
  filter(!is.na(south), !is.na(north)) %>%
  mutate(gap = south - north)

# ------------------------------------------------------------
# 2. FUNCTION: run the directional test for one wind station
# ------------------------------------------------------------

run_station <- function(wind_station, file) {
  wind <- read_csv(paste0(WIND_DIR, file), show_col_types = FALSE) %>%
    mutate(
      date = as.Date(substr(time, 1, 10)),
      wind_speed = as.numeric(vv_mittel),
      dd_code = as.numeric(dd32_i),
      wind_deg = dd_code * 11.25,
      ns_component = cos(wind_deg * pi / 180),      # +1 = northerly, -1 = southerly
      ns_forcing = wind_speed * ns_component
    ) %>%
    filter(!is.na(ns_forcing)) %>%
    select(date, ns_forcing)

  df <- inner_join(water, wind, by = "date")
  if (nrow(df) < 30) {
    return(tibble(wind_station = wind_station, n_days = nrow(df),
                  date_from = NA, date_to = NA, r = NA, p_value = NA))
  }

  ct <- cor.test(df$ns_forcing, df$gap, use = "complete.obs")
  tibble(
    wind_station = wind_station,
    n_days = nrow(df),
    date_from = format(min(df$date), "%Y-%m"),
    date_to = format(max(df$date), "%Y-%m"),
    r = round(unname(ct$estimate), 3),
    p_value = ct$p.value
  )
}

# ------------------------------------------------------------
# 3. RUN ALL THREE AND ASSEMBLE SUMMARY TABLE
# ------------------------------------------------------------

results <- pmap_dfr(wind_files, run_station)

results_fmt <- results %>%
  mutate(p_value = ifelse(is.na(p_value), NA, format.pval(p_value, digits = 2, eps = 1e-16)))

cat("\n=== Three-station seiche robustness check ===\n")
cat("Directional correlation (north-south wind forcing vs. Apetlon-Neusiedl gauge gap)\n")
cat("A consistent POSITIVE r across all three stations indicates the seiche\n")
cat("signal is not an artifact of any single wind station.\n\n")
print(results_fmt)

write_csv(results_fmt, paste0(OUT_DIR, "seiche_three_station_robustness.csv"))
cat("\nSaved to:", paste0(OUT_DIR, "seiche_three_station_robustness.csv"), "\n")

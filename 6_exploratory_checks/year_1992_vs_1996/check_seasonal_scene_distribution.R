# ============================================================
# SEASONAL SCENE DISTRIBUTION: 1992 vs 1996
# ------------------------------------------------------------
# Checks whether the 1992 flood frequency peak relative to 1996
# can be explained by scenes being acquired disproportionately
# during wetter (spring/early summer) months in 1992.
#
# Uses the pixel CSV which has one row per pixel per scene date.
# Since every scene date appears once per pixel, we extract
# unique scene dates per year and examine their monthly distribution.
# ============================================================

library(tidyverse)

IN_CSV  <- "<PROJECT_ROOT>/result/writting/new/Data/FLATS_newpixels_1984_2025_CLEANED.csv"
OUT_DIR <- "<PROJECT_ROOT>/result/writting/new/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

pix <- read_csv(IN_CSV, show_col_types = FALSE) %>%
  mutate(date  = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d")),
         year  = year(date),
         month = month(date))

# Unique scene dates (one row per date, not per pixel)
scenes <- pix %>%
  distinct(date, year, month, sensor) %>%
  arrange(date)

cat("=== TOTAL UNIQUE SCENE DATES ===\n")
cat("1992:", sum(scenes$year == 1992), "scenes\n")
cat("1996:", sum(scenes$year == 1996), "scenes\n\n")

# Monthly distribution for 1992 and 1996
monthly <- scenes %>%
  filter(year %in% c(1992, 1996)) %>%
  count(year, month) %>%
  mutate(month_name = month.abb[month])

cat("=== MONTHLY SCENE DISTRIBUTION ===\n")
print(monthly %>%
        arrange(month_name) %>%
        pivot_wider(names_from = year,
                    values_from = n,
                    values_fill = 0))

# Define "wet season" as months when flooding is most likely
# At Lake Neusiedl, spring (March-May) and early summer (June)
# typically have higher water levels
wet_months   <- 3:6    # March to June
other_months <- c(1:2, 7:12)

season_summary <- scenes %>%
  filter(year %in% c(1992, 1996)) %>%
  mutate(season = if_else(month %in% wet_months,
                          "Spring/early summer (Mar-Jun)",
                          "Rest of year")) %>%
  count(year, season) %>%
  group_by(year) %>%
  mutate(pct = round(n / sum(n) * 100, 1))

cat("\n=== SEASONAL BREAKDOWN ===\n")
print(season_summary)

cat("\n=== INTERPRETATION ===\n")
cat("If 1992 has a higher % of scenes in spring/early summer than 1996,\n")
cat("this supports the hypothesis that scene timing explains the\n")
cat("higher flood frequency in 1992 despite lower annual water level.\n\n")

# Also check mean mNDWI per scene date for both years
# This tells us if 1992 scenes were genuinely wetter on average
scene_mndwi <- pix %>%
  filter(year %in% c(1992, 1996)) %>%
  group_by(year, date, month) %>%
  summarise(mean_mNDWI = mean(mNDWI, na.rm = TRUE),
            .groups = "drop")

cat("=== MEAN mNDWI PER SCENE BY MONTH ===\n")
monthly_mndwi <- scene_mndwi %>%
  group_by(year, month) %>%
  summarise(mean_mNDWI = round(mean(mean_mNDWI), 4),
            n_scenes   = n(),
            .groups    = "drop") %>%
  mutate(month_name = month.abb[month])

print(monthly_mndwi %>%
        select(year, month_name, n_scenes, mean_mNDWI) %>%
        arrange(month_name) %>%
        pivot_wider(names_from  = year,
                    values_from = c(n_scenes, mean_mNDWI)))

# Annual mean mNDWI for both years
cat("\n=== ANNUAL MEAN mNDWI ===\n")
annual_mndwi <- scene_mndwi %>%
  group_by(year) %>%
  summarise(mean_mNDWI_annual = round(mean(mean_mNDWI), 4),
            .groups = "drop")
print(annual_mndwi)

# Save results
write_csv(monthly,        paste0(OUT_DIR, "scene_monthly_1992_1996.csv"))
write_csv(season_summary, paste0(OUT_DIR, "scene_seasonal_1992_1996.csv"))
write_csv(monthly_mndwi,  paste0(OUT_DIR, "mNDWI_monthly_1992_1996.csv"))

cat("\nAll outputs saved to:", OUT_DIR, "\n")
cat("\nBring these numbers to the thesis only after checking the results.\n")
cat("If spring/early summer scenes are proportionally more common in 1992\n")
cat("AND mean mNDWI is higher in those months, the hypothesis is confirmed.\n")
cat("If not, remove the seasonal-timing explanation from the thesis text.\n")

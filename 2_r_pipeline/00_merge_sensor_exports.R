# ============================================================
# MERGE PER-SENSOR EXPORTS INTO ONE PIXEL-DATE DATASET
# Lake Neusiedl | 1984-2025
# ------------------------------------------------------------
# FIXES / ADDITIONS relative to the original version:
#   1. Reports which columns get dropped by the common-column
#      intersection (the original silently kept only shared columns
#      with no visibility into what was excluded, if anything).
#   2. Explicitly coerces key columns to their expected type instead
#      of relying on bind_rows() to guess - catches silent NA
#      introduction from type mismatches between the four files.
#   3. Validates that each sensor's dates actually fall inside its
#      expected filter window - catches a file assigned to the
#      wrong variable, which bind_rows() would not otherwise flag.
#   4. Adds per-sensor and per-station summary tables, matching the
#      reporting style already used in 5_2.R.
#   5. Compares the new merged dataset against the previous one
#      (FLATS_newpixels_1984_2025.csv, produced by the old GEE
#      script), so you can see exactly what changed after switching
#      to the corrected cloud thresholds and date ranges.
# ============================================================

library(readr)
library(dplyr)

base_dir <- "<PROJECT_ROOT>/result/writting/new/Data"

file_l5 <- file.path(base_dir, "FLATS_newpixels_L5_1984_2012.csv")
file_l7 <- file.path(base_dir, "FLATS_newpixels_L7_1999_2022.csv")
file_l8 <- file.path(base_dir, "FLATS_newpixels_L8_2013_2025.csv")
file_l9 <- file.path(base_dir, "FLATS_newpixels_L9_2021_2025.csv")

outfile <- file.path(base_dir, "FLATS_newpixels_1984_2025_merged.csv")
dupfile <- file.path(base_dir, "FLATS_newpixels_duplicates_check.csv")

# Path to the PREVIOUS merged dataset, for comparison. Adjust this if it
# lives somewhere else - set to NULL to skip the comparison entirely.
old_file <- "<PROJECT_ROOT>/Data/FLATS_newpixels_1984_2025.csv"

# Expected sensor date windows, matching GEE_pixel_extraction_corrected.txt -
# used below to validate that each file actually contains what its name claims.
expected_ranges <- list(
  L5 = c(as.Date("1984-01-01"), as.Date("2012-12-31")),
  L7 = c(as.Date("1999-01-01"), as.Date("2022-12-31")),
  L8 = c(as.Date("2013-01-01"), as.Date("2025-12-31")),
  L9 = c(as.Date("2021-01-01"), as.Date("2025-12-31"))
)

# ------------------------------------------------------------
# 1. CHECK FILES EXIST
# ------------------------------------------------------------

files <- c(file_l5, file_l7, file_l8, file_l9)
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"))
}

# ------------------------------------------------------------
# 2. LOAD EACH FILE
# ------------------------------------------------------------

l5 <- read_csv(file_l5, show_col_types = FALSE)
l7 <- read_csv(file_l7, show_col_types = FALSE)
l8 <- read_csv(file_l8, show_col_types = FALSE)
l9 <- read_csv(file_l9, show_col_types = FALSE)

cat("Rows loaded per sensor (before any merging or filtering):\n")
cat("  L5:", nrow(l5), "\n  L7:", nrow(l7), "\n  L8:", nrow(l8), "\n  L9:", nrow(l9), "\n\n")

# ------------------------------------------------------------
# 3. COLUMN CONSISTENCY CHECK (FIX #1)
# ------------------------------------------------------------

all_cols     <- Reduce(union, list(names(l5), names(l7), names(l8), names(l9)))
common_cols  <- Reduce(intersect, list(names(l5), names(l7), names(l8), names(l9)))
dropped_cols <- setdiff(all_cols, common_cols)

if (length(dropped_cols) > 0) {
  cat("WARNING - the following columns are NOT present in all four files\n")
  cat("and will be DROPPED from the merged dataset:\n")
  cat(" ", paste(dropped_cols, collapse = ", "), "\n")
  cat("Check whether this is expected before proceeding.\n\n")
} else {
  cat("Column check OK - all four files share the same columns.\n\n")
}

l5 <- l5 %>% select(all_of(common_cols))
l7 <- l7 %>% select(all_of(common_cols))
l8 <- l8 %>% select(all_of(common_cols))
l9 <- l9 %>% select(all_of(common_cols))

# ------------------------------------------------------------
# 4. EXPLICIT TYPE COERCION (FIX #2)
#    Done explicitly, one column at a time, so any coercion failure
#    (e.g. a non-numeric value landing in a numeric column) shows up
#    as a visible NA count below, rather than silently inside bind_rows().
# ------------------------------------------------------------

coerce_types <- function(df) {
  df %>%
    mutate(
      pixel_id = as.character(pixel_id),
      sensor   = as.character(sensor),
      station  = as.character(station),
      date     = as.Date(date),
      year     = as.integer(year),
      month    = as.integer(month),
      elevation = as.numeric(elevation),
      NDWI  = as.numeric(NDWI),
      mNDWI = as.numeric(mNDWI),
      EVI   = as.numeric(EVI)
    )
}

l5 <- coerce_types(l5); l7 <- coerce_types(l7); l8 <- coerce_types(l8); l9 <- coerce_types(l9)

merged <- bind_rows(l5, l7, l8, l9) %>% arrange(pixel_id, date, sensor)

na_report <- merged %>%
  summarise(
    n_rows = n(),
    na_date = sum(is.na(date)),
    na_NDWI = sum(is.na(NDWI)),
    na_mNDWI = sum(is.na(mNDWI)),
    na_EVI = sum(is.na(EVI))
  )
cat("NA check after type coercion (should be 0, or explainable):\n")
print(na_report)
cat("\n")

# ------------------------------------------------------------
# 5. VALIDATE SENSOR DATE RANGES (FIX #3)
#    Confirms each file's dates actually fall inside the range its
#    filename claims - catches a file mix-up that bind_rows() alone
#    would not detect.
# ------------------------------------------------------------

cat("Sensor date-range validation:\n")
for (s in names(expected_ranges)) {
  sub <- merged %>% filter(sensor == s)
  rng <- expected_ranges[[s]]
  out_of_range <- sub %>% filter(date < rng[1] | date > rng[2])
  cat(sprintf("  %s: observed %s to %s | expected %s to %s | out-of-range rows: %d\n",
              s, min(sub$date, na.rm = TRUE), max(sub$date, na.rm = TRUE),
              rng[1], rng[2], nrow(out_of_range)))
}
cat("\n")

# ------------------------------------------------------------
# 6. DEDUPLICATION
# ------------------------------------------------------------

dup_check <- merged %>%
  count(pixel_id, date, sensor, name = "n") %>%
  filter(n > 1)

merged_clean <- merged %>%
  distinct(pixel_id, date, sensor, .keep_all = TRUE) %>%
  arrange(pixel_id, date, sensor)

write_csv(merged_clean, outfile)
write_csv(dup_check, dupfile)

# ------------------------------------------------------------
# 7. PER-SENSOR AND PER-STATION SUMMARY (FIX #4)
# ------------------------------------------------------------

sensor_summary <- merged_clean %>%
  group_by(sensor) %>%
  summarise(n_obs = n(), n_pixels = n_distinct(pixel_id),
            date_min = min(date), date_max = max(date), .groups = "drop")
cat("Per-sensor summary (final merged dataset):\n")
print(sensor_summary)
cat("\n")

station_summary <- merged_clean %>%
  group_by(station) %>%
  summarise(n_obs = n(), n_pixels = n_distinct(pixel_id), .groups = "drop")
cat("Per-station summary (final merged dataset):\n")
print(station_summary)
cat("\n")

# ------------------------------------------------------------
# 8. COMPARISON AGAINST THE PREVIOUS DATASET (FIX #5)
# ------------------------------------------------------------

if (!is.null(old_file) && file.exists(old_file)) {
  old <- read_csv(old_file, show_col_types = FALSE) %>%
    mutate(date = as.Date(date))

  cat("================ COMPARISON: OLD vs NEW MERGED DATASET ================\n")
  cat("Old dataset rows:", nrow(old), " | New dataset rows:", nrow(merged_clean),
      " | Difference:", nrow(merged_clean) - nrow(old), "\n")
  cat("Old unique pixels:", n_distinct(old$pixel_id),
      " | New unique pixels:", n_distinct(merged_clean$pixel_id), "\n")

  if ("sensor" %in% names(old)) {
    old_sensor_summary <- old %>% group_by(sensor) %>% summarise(n_obs = n(), .groups = "drop")
    cat("\nOld per-sensor row counts:\n"); print(old_sensor_summary)
    cat("\nNew per-sensor row counts:\n"); print(sensor_summary %>% select(sensor, n_obs))
  }
  cat("=========================================================================\n\n")
} else {
  cat("No previous dataset found at old_file path - skipping comparison.\n")
  cat("Set old_file to NULL at the top of this script to silence this message.\n\n")
}

# ------------------------------------------------------------
# 9. FINAL REPORT
# ------------------------------------------------------------

cat("Merged file written to:\n", outfile, "\n")
cat("Duplicate check written to:\n", dupfile, "\n")
cat("Rows before deduplication:", nrow(merged), "\n")
cat("Rows after deduplication :", nrow(merged_clean), "\n")
cat("Duplicate rows removed   :", nrow(merged) - nrow(merged_clean), "\n")
cat("Columns:", ncol(merged_clean), "\n")
cat("Unique pixels:", n_distinct(merged_clean$pixel_id), "\n")
cat("Date range:", as.character(min(merged_clean$date, na.rm = TRUE)),
    "to", as.character(max(merged_clean$date, na.rm = TRUE)), "\n")

cat("\nIMPORTANT: 5_2.R expects its input file to be named\n")
cat("'FLATS_newpixels_1984_2025.csv', not 'FLATS_newpixels_1984_2025_merged.csv'.\n")
cat("Either rename this output, or update PIX_PATH in 5_2.R to point here,\n")
cat("before running the cleaning step.\n")

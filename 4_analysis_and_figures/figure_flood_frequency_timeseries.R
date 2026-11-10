# ============================================================
# FIGURE: Annual Flood Frequency (NDWI+EVI vs mNDWI+EVI) with
#          Landsat scene counts and mean annual water level
# Lake Neusiedl, Austria | 1984-2025
# ------------------------------------------------------------
# Reads:
#   FloodMap_annual_timeseries.csv   (from your GEE pipeline)
#   water_mean_level/*-<year>.txt    (one file per year, daily
#                                      "date  level" rows)
# Produces:
#   Figure5_8_annual_flood_timeseries.png  (top: frequency + scenes,
#                                            bottom: water level)
# ============================================================

library(tidyverse)
library(patchwork)

# ------------------------------------------------------------
# 0. PATHS - adjust only if your folders move
# ------------------------------------------------------------

FLOOD_CSV     <- "<PROJECT_ROOT>/result/writting/GEE/FloodMap_annual_timeseries.csv"
WATERLEVEL_DIR <- "<PROJECT_ROOT>/Data/water_mean_level"
OUT_DIR        <- "<PROJECT_ROOT>/result/writting/GEE/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

START_YEAR <- 1984
END_YEAR   <- 2025

# -- Extreme years, corrected against the full ranked gauge record --
# (see rank_water_level_years.R - 1996 = highest mean annual level of all
# 42 years; 2022 and 2023 = the two lowest, 0.18 m below the next-driest
# year, i.e. a genuine two-year drought trough, not a marginal pick)
WETTEST_YEAR <- 1996
DRIEST_YEARS <- c(2022, 2023)

# -- Shaded background periods --
# CRISIS_START/END covers the historically documented 2000-2004 lake level
# crisis (Section 3.1) - this is a separate, independently cited event and
# is NOT the same as the statistically-ranked driest years above, so it is
# left unchanged even though it no longer contains a DRIEST_YEARS marker.
CRISIS_START <- 2000; CRISIS_END <- 2004

# DRYING_START/END widened to run through 2023 so the shaded "lake drying"
# period visually contains both corrected driest-year markers, since 2022
# and 2023 are the culmination of this same post-2015 drying trend
# (Section 3.1; Sailer and Maracek 2019), not a separate event.
DRYING_START <- 2019.5; DRYING_END <- 2023.5

# ------------------------------------------------------------
# 1. LOAD FLOOD FREQUENCY + SCENE COUNT DATA
# ------------------------------------------------------------

flood <- read_csv(FLOOD_CSV, show_col_types = FALSE) %>%
  filter(year >= START_YEAR, year <= END_YEAR)

# ------------------------------------------------------------
# 2. LOAD AND AGGREGATE DAILY WATER LEVEL FILES (ONE PER YEAR)
#    Matches any file ending in "-<year>.txt", regardless of the
#    export-hash prefix, so this doesn't break if the prefix differs
#    slightly between yearly downloads.
# ------------------------------------------------------------

read_year_level <- function(year) {
  pattern <- paste0("-", year, "\\.txt$")
  f <- list.files(WATERLEVEL_DIR, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) {
    warning(paste("No water level file found for", year))
    return(tibble(year = year, mean_level = NA_real_, n_days = 0))
  }
  if (length(f) > 1) {
    warning(paste("Multiple files matched for", year, "- using the first:", f[1]))
    f <- f[1]
  }
  df <- read.table(f, header = FALSE, col.names = c("date", "level"),
                    colClasses = c("character", "numeric"))
  tibble(year = year, mean_level = mean(df$level, na.rm = TRUE), n_days = nrow(df))
}

waterlevel <- map_dfr(START_YEAR:END_YEAR, read_year_level)

cat("Water level files loaded for", sum(!is.na(waterlevel$mean_level)),
    "of", length(START_YEAR:END_YEAR), "years.\n")
if (any(is.na(waterlevel$mean_level))) {
  cat("Missing years:", paste(waterlevel$year[is.na(waterlevel$mean_level)], collapse = ", "), "\n")
}

# ------------------------------------------------------------
# 3. MERGE
# ------------------------------------------------------------

full_data <- flood %>% left_join(waterlevel, by = "year")

# ------------------------------------------------------------
# 4. STYLE CONSTANTS (kept consistent with earlier thesis figures)
# ------------------------------------------------------------

col_baseline <- "#e07b39"   # NDWI + EVI
col_improved <- "#01696f"   # mNDWI + EVI
col_bar      <- "grey80"
col_crisis   <- "#e07b39"
col_drying   <- "#5b8fa8"

# Scale factor to draw scene-count bars on the same axis as frequency (0-1).
# Capped so the tallest bar never exceeds the lowest frequency value in the
# series - otherwise, in years with many scenes but low flood frequency
# (common in your data post-2018), the bars visually bury the lines.
max_scenes <- max(full_data$n_scenes, na.rm = TRUE)
min_freq   <- min(c(full_data$NDWI_EVI_freq, full_data$mNDWI_EVI_freq), na.rm = TRUE)
scene_scale <- (min_freq * 0.85) / max_scenes

# ------------------------------------------------------------
# 5. TOP PANEL - FLOOD FREQUENCY + SCENE COUNT
# ------------------------------------------------------------

# ------------------------------------------------------------
# 4b. COMPUTE CORRELATION FOR IN-PLOT ANNOTATION
#     Same calculation as Section 5.8 / 4.7 - reusing it here so
#     the figure displays the exact value already reported in text,
#     rather than a separately re-derived number.
# ------------------------------------------------------------

r_ndwi  <- cor(full_data$mean_level, full_data$NDWI_EVI_freq,  use = "complete.obs")
r_mndwi <- cor(full_data$mean_level, full_data$mNDWI_EVI_freq, use = "complete.obs")

corr_label <- sprintf("Pearson r (vs. water level):\nNDWI+EVI: %.3f   mNDWI+EVI: %.3f", r_ndwi, r_mndwi)

p_top <- ggplot(full_data, aes(x = year)) +

  # background shading for known lake events
  annotate("rect", xmin = CRISIS_START, xmax = CRISIS_END, ymin = -Inf, ymax = Inf,
           fill = col_crisis, alpha = 0.12) +
  annotate("rect", xmin = DRYING_START, xmax = DRYING_END, ymin = -Inf, ymax = Inf,
           fill = col_drying, alpha = 0.15) +
  annotate("text", x = (CRISIS_START + CRISIS_END) / 2, y = 0.97,
           label = "Lake level\ncrisis", size = 3.4, color = col_crisis, lineheight = 0.9) +
  annotate("text", x = (DRYING_START + DRYING_END) / 2, y = 0.97,
           label = "Lake\ndrying", size = 3.4, color = col_drying, lineheight = 0.9) +

  # scene-count bars (rescaled), drawn first so lines sit on top
  geom_col(aes(y = n_scenes * scene_scale), fill = col_bar, width = 0.65, alpha = 0.75) +

  # flood frequency lines
  geom_line(aes(y = NDWI_EVI_freq, color = "Baseline (NDWI + EVI)"),
             linewidth = 0.9, linetype = "dashed") +
  geom_line(aes(y = mNDWI_EVI_freq, color = "Improved (mNDWI + EVI)"),
             linewidth = 1) +

  # wettest / driest year markers
  geom_vline(xintercept = WETTEST_YEAR, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = DRIEST_YEARS, linetype = "dotted", color = "grey40", linewidth = 0.4) +

  scale_color_manual(values = c("Baseline (NDWI + EVI)" = col_baseline,
                                 "Improved (mNDWI + EVI)" = col_improved),
                      name = NULL) +
  scale_x_continuous(breaks = seq(START_YEAR, END_YEAR, 4)) +
  scale_y_continuous(
    name = "Mean flood frequency",
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    sec.axis = sec_axis(~ . / scene_scale, name = "Number of Landsat scenes",
                         breaks = seq(0, 120, by = 40))
  ) +
  labs(title = "Lake Neusiedl (Austria) - Annual Flood Frequency, 1984-2025") +
  annotate("label", x = START_YEAR + 1.5, y = 0.99, label = corr_label,
           hjust = 0, vjust = 1, size = 3.2, label.size = 0.2,
           fill = "white", color = "grey20") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 13)
  )

# ------------------------------------------------------------
# 6. BOTTOM PANEL - MEAN ANNUAL WATER LEVEL
# ------------------------------------------------------------

p_bottom <- ggplot(full_data, aes(x = year, y = mean_level)) +
  annotate("rect", xmin = CRISIS_START, xmax = CRISIS_END, ymin = -Inf, ymax = Inf,
           fill = col_crisis, alpha = 0.12) +
  annotate("rect", xmin = DRYING_START, xmax = DRYING_END, ymin = -Inf, ymax = Inf,
           fill = col_drying, alpha = 0.15) +
  geom_vline(xintercept = WETTEST_YEAR, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = DRIEST_YEARS, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  geom_line(color = "#33465c", linewidth = 0.9) +
  geom_point(color = "#33465c", size = 1.4) +
  scale_x_continuous(breaks = seq(START_YEAR, END_YEAR, 4)) +
  labs(x = "Year", y = "Mean annual\nwater level (m a.s.l.)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------
# 7. COMBINE AND SAVE
# ------------------------------------------------------------

final_plot <- p_top / p_bottom + plot_layout(heights = c(2.4, 1))

ggsave(paste0(OUT_DIR, "Figure5_8_annual_flood_timeseries.png"), final_plot,
       width = 10, height = 8, dpi = 320, bg = "white")

cat("\nSaved to:", paste0(OUT_DIR, "Figure5_8_annual_flood_timeseries.png"), "\n")

# ------------------------------------------------------------
# 8. CORRELATION: annual flood frequency vs. mean annual water level
# ------------------------------------------------------------
# Diagnostic check BEFORE running cor.test() - confirms exactly how many
# years are actually going into the correlation, and flags if any expected
# years silently dropped out due to missing water level files or missing
# flood frequency rows. cor.test() drops incomplete pairs automatically
# but does not print how many rows it used, so this makes that explicit.

cat("\n--- Correlation input diagnostics ---\n")
cat("Total rows in full_data:", nrow(full_data), "\n")
cat("Rows with non-missing mean_level:", sum(!is.na(full_data$mean_level)), "\n")
cat("Rows with non-missing NDWI_EVI_freq:", sum(!is.na(full_data$NDWI_EVI_freq)), "\n")
cat("Rows with non-missing mNDWI_EVI_freq:", sum(!is.na(full_data$mNDWI_EVI_freq)), "\n")

years_missing_level <- full_data$year[is.na(full_data$mean_level)]
years_missing_freq  <- full_data$year[is.na(full_data$NDWI_EVI_freq) | is.na(full_data$mNDWI_EVI_freq)]

if (length(years_missing_level) > 0) {
  cat("Years dropped due to missing water level:", paste(years_missing_level, collapse = ", "), "\n")
} else {
  cat("No years dropped due to missing water level.\n")
}
if (length(years_missing_freq) > 0) {
  cat("Years dropped due to missing flood frequency:", paste(years_missing_freq, collapse = ", "), "\n")
} else {
  cat("No years dropped due to missing flood frequency.\n")
}

n_used <- sum(!is.na(full_data$mean_level) & !is.na(full_data$NDWI_EVI_freq) & !is.na(full_data$mNDWI_EVI_freq))
cat("Years actually used in both correlation tests below: n =", n_used,
    "(expected up to", length(START_YEAR:END_YEAR), ")\n\n")

cor_ndwi  <- cor.test(full_data$mean_level, full_data$NDWI_EVI_freq,  method = "pearson")
cor_mndwi <- cor.test(full_data$mean_level, full_data$mNDWI_EVI_freq, method = "pearson")

cat("NDWI+EVI  vs water level: r =", round(cor_ndwi$estimate, 3),
    " p =", format.pval(cor_ndwi$p.value, digits = 3),
    " n =", cor_ndwi$parameter + 2, "\n")
cat("mNDWI+EVI vs water level: r =", round(cor_mndwi$estimate, 3),
    " p =", format.pval(cor_mndwi$p.value, digits = 3),
    " n =", cor_mndwi$parameter + 2, "\n")

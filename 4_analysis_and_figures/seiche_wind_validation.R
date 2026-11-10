# ============================================================
# SEICHE VALIDATION: does wind widen the north-south gauge gap?
# Lake Neusiedl | 1984-2025
# ------------------------------------------------------------
# PURPOSE
# Tests, using the study's own data, the seiche assumption that
# Section 4.3.2 / 6.6 currently support only by citation (BMLUK 2020):
# that strong wind along the lake's axis raises water level at one
# shore relative to the other. If confirmed, this converts a cited
# assumption into a demonstrated one.
#
# METHOD
# 1. Compute the daily water-level DIFFERENCE between a southern gauge
#    (Apetlon) and a northern gauge (Neusiedl).
# 2. Join the Neusiedl daily wind record (speed + direction).
# 3. Test whether that difference depends on wind speed and direction
#    in the physically expected way (strong N wind pushes water south,
#    raising Apetlon relative to Neusiedl; strong S wind does the reverse).
#
# WIND DIRECTION ENCODING (dd32, GeoSphere Austria klima-v2-1d,
# DOI: 10.60669/gs6w-jd70): code x 11.25 = degrees. Station reports
# 8-point resolution: 0/32=N, 4=NE, 8=E, 12=SE, 16=S, 20=SW, 24=W, 28=NW.
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 0. CONFIG
# ------------------------------------------------------------

WATER_CSV <- "<PROJECT_ROOT>/Data/all_stations_1984_2025.csv"
WIND_CSV  <- "<PROJECT_ROOT>/Data/Neusiedl_wind.csv"
OUT_DIR   <- "<PROJECT_ROOT>/result/writting/new/seiche/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Which two gauges define the north-south axis. Apetlon is the southern
# shore (lat 47.68), Neusiedl the northern (lat 47.93).
SOUTH_GAUGE <- "Apetlon"
NORTH_GAUGE <- "Neusiedl"

# ------------------------------------------------------------
# 1. LOAD AND RESHAPE WATER LEVEL DATA
#    Adjust column handling to match the actual structure of
#    all_stations_1984_2025.csv - this assumes a 'date' column plus
#    one column per station. If it is in long format instead, see
#    the alternative block below.
# ------------------------------------------------------------

water_raw <- read_csv(WATER_CSV, show_col_types = FALSE)
cat("Water level file columns:\n"); print(names(water_raw))
cat("\nFirst rows:\n"); print(head(water_raw, 3))

# --- WIDE format assumption (one column per station) ---
# If your file has columns like date, Apetlon, Neusiedl, Illmitz, ...
water <- water_raw %>%
  mutate(date = as.Date(.[[1]])) %>%          # assumes first column is the date
  select(date, all_of(c(SOUTH_GAUGE, NORTH_GAUGE))) %>%
  rename(south = all_of(SOUTH_GAUGE), north = all_of(NORTH_GAUGE)) %>%
  filter(!is.na(south), !is.na(north)) %>%
  mutate(gap = south - north)   # positive = southern shore higher

# --- If instead your file is LONG format (date, station, level), use: ---
# water <- water_raw %>%
#   mutate(date = as.Date(date)) %>%
#   filter(station %in% c(SOUTH_GAUGE, NORTH_GAUGE)) %>%
#   select(date, station, level) %>%
#   pivot_wider(names_from = station, values_from = level) %>%
#   rename(south = all_of(SOUTH_GAUGE), north = all_of(NORTH_GAUGE)) %>%
#   filter(!is.na(south), !is.na(north)) %>%
#   mutate(gap = south - north)

cat("\nWater-level gap (south - north) summary:\n")
print(summary(water$gap))

# ------------------------------------------------------------
# 2. LOAD AND DECODE WIND DATA
# ------------------------------------------------------------

wind <- read_csv(WIND_CSV, show_col_types = FALSE) %>%
  mutate(
    date = as.Date(substr(time, 1, 10)),
    wind_speed = as.numeric(vv_mittel),
    dd_code = as.numeric(dd32_i)
  ) %>%
  filter(!is.na(wind_speed), !is.na(dd_code)) %>%
  select(date, wind_speed, dd_code)

# Decode dd32 code into 8-point compass label and degrees
compass_lookup <- c("0"="N","4"="NE","8"="E","12"="SE",
                    "16"="S","20"="SW","24"="W","28"="NW","32"="N")

wind <- wind %>%
  mutate(
    wind_deg = dd_code * 11.25,
    wind_dir = compass_lookup[as.character(dd_code)],
    # North-south component: +1 = pure northerly (pushes water south),
    # -1 = pure southerly. cos(deg) is +1 at N (0 deg), -1 at S (180 deg).
    ns_component = cos(wind_deg * pi / 180),
    # Signed "north-south wind forcing" = speed scaled by how N-S aligned it is
    ns_forcing = wind_speed * ns_component
  )

cat("\nWind direction frequency:\n")
print(wind %>% count(wind_dir, sort = TRUE))

# ------------------------------------------------------------
# 3. JOIN AND ANALYSE
# ------------------------------------------------------------

df <- inner_join(water, wind, by = "date")
cat("\nDays with both water gap and wind data:", nrow(df), "\n")

# --- Test 1: does the gap correlate with north-south wind forcing? ---
# Physical expectation: stronger northerly forcing (ns_forcing > 0)
# should push water toward the south shore, INCREASING south - north gap.
# So we expect a POSITIVE correlation between ns_forcing and gap.
r_forcing <- cor.test(df$ns_forcing, df$gap, use = "complete.obs")
cat("\n--- Test 1: N-S wind forcing vs. water-level gap ---\n")
cat("Pearson r =", round(r_forcing$estimate, 3),
    " p =", format.pval(r_forcing$p.value, digits = 3), "\n")
cat("(Positive r supports the seiche mechanism: northerly wind raises the south shore.)\n")

# --- Test 2: compare the gap on strong-N vs strong-S wind days ---
strong <- quantile(df$wind_speed, 0.75, na.rm = TRUE)  # "windy" = top quartile
df_dir <- df %>%
  filter(wind_speed >= strong) %>%
  mutate(regime = case_when(
    wind_dir %in% c("N", "NE", "NW") ~ "Strong northerly",
    wind_dir %in% c("S", "SE", "SW") ~ "Strong southerly",
    TRUE ~ "Cross-axis (E/W)"
  ))

regime_summary <- df_dir %>%
  group_by(regime) %>%
  summarise(n_days = n(), mean_gap = mean(gap), sd_gap = sd(gap), .groups = "drop")
cat("\n--- Test 2: mean gap by wind regime (windy days only, top quartile speed) ---\n")
print(regime_summary)
cat("Expectation: 'Strong northerly' mean_gap > 'Strong southerly' mean_gap.\n")

# --- Test 3: calm-day baseline for comparison ---
calm <- quantile(df$wind_speed, 0.25, na.rm = TRUE)
calm_gap <- df %>% filter(wind_speed <= calm) %>% summarise(mean_gap = mean(gap), sd_gap = sd(gap))
cat("\n--- Test 3: gap on calm days (bottom quartile speed) ---\n")
print(calm_gap)
cat("Expectation: calm-day gap is near zero with small SD; windy days show larger, direction-dependent gaps.\n")

# ------------------------------------------------------------
# 4. VISUALISE
# ------------------------------------------------------------

p1 <- ggplot(df, aes(x = ns_forcing, y = gap)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "lm", color = "#01696f") +
  labs(title = "Water-level gap (Apetlon - Neusiedl) vs. north-south wind forcing",
       subtitle = "Positive forcing = northerly wind; positive gap = southern shore higher",
       x = "North-south wind forcing (m/s, + = northerly)",
       y = "Water-level gap: Apetlon - Neusiedl (m)") +
  theme_minimal(base_size = 12)

ggsave(paste0(OUT_DIR, "seiche_wind_vs_gap.png"), p1, width = 9, height = 6, dpi = 300, bg = "white")

library(patchwork)

df_dir <- df_dir %>%
  mutate(regime = factor(regime,
    levels = c("Strong southerly", "Cross-axis (E/W)", "Strong northerly")))

regime_cols <- c("Strong southerly" = "#c0392b",
                 "Cross-axis (E/W)" = "#95a5a6",
                 "Strong northerly" = "#2471a3")

# Date range and total-day context, computed from the actual joined data
date_min <- format(min(df$date), "%Y")
date_max <- format(max(df$date), "%Y")
n_total_days <- nrow(df)
n_windy_days <- nrow(df_dir)
date_caption <- sprintf(
  "Neusiedl wind record %s-%s | %s days with paired wind and gauge data; %s windy days (top-quartile speed) shown",
  date_min, date_max, format(n_total_days, big.mark = ","), format(n_windy_days, big.mark = ","))

# --- Panel A: distribution of the gap by wind regime, with the
#     10 cm elevation-bin width shown for scale ---
bin_width <- 0.10

pA <- ggplot(df_dir, aes(x = regime, y = gap, fill = regime)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey60") +
  geom_violin(alpha = 0.35, color = NA, show.legend = FALSE) +
  geom_boxplot(width = 0.25, outlier.alpha = 0.08, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3.5,
               color = "black", show.legend = FALSE) +
  scale_fill_manual(values = regime_cols) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = -bin_width/2, ymax = bin_width/2,
           alpha = 0.12, fill = "black") +
  annotate("text", x = 0.6, y = bin_width/2 + 0.005,
           label = "10 cm elevation bin", hjust = 0, size = 3, color = "grey30") +
  labs(title = "A. Wind direction shifts the north-south water-level gap",
       subtitle = "Windy days only (top-quartile wind speed); diamond = mean",
       x = "Wind regime (Neusiedl station)",
       y = "Water-level gap: Apetlon (S) - Neusiedl (N), m") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# --- Panel B: mean gap per regime as a bar, annotated with the
#     physical interpretation, the swing magnitude, and n days ---
mean_n <- regime_summary$mean_gap[regime_summary$regime == "Strong northerly"]
mean_s <- regime_summary$mean_gap[regime_summary$regime == "Strong southerly"]
swing_cm <- round((mean_n - mean_s) * 100, 1)

pB <- regime_summary %>%
  mutate(regime = factor(regime,
    levels = c("Strong southerly", "Cross-axis (E/W)", "Strong northerly"))) %>%
  ggplot(aes(x = regime, y = mean_gap, fill = regime)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_gap - sd_gap/sqrt(n_days),
                    ymax = mean_gap + sd_gap/sqrt(n_days)),
                width = 0.15) +
  geom_text(aes(label = sprintf("%+.1f cm", mean_gap * 100),
                vjust = ifelse(mean_gap >= 0, -0.8, 1.6)), size = 3.5) +
  # n days shown as a label just inside the base of each bar
  geom_text(aes(y = 0, label = paste0("n = ", format(n_days, big.mark = ","))),
            vjust = ifelse(regime_summary$mean_gap[match(
              c("Strong southerly","Cross-axis (E/W)","Strong northerly"),
              regime_summary$regime)] >= 0, 1.6, -0.8),
            size = 3, color = "grey25") +
  scale_fill_manual(values = regime_cols) +
  labs(title = "B. Mean gap by regime",
       subtitle = paste0("North-to-south swing: ", swing_cm, " cm"),
       x = NULL, y = "Mean water-level gap (m)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1))

seiche_fig <- (pA + pB + plot_layout(widths = c(1.7, 1))) +
  plot_annotation(
    caption = date_caption,
    theme = theme(plot.caption = element_text(size = 9, color = "grey40", hjust = 0))
  )

ggsave(paste0(OUT_DIR, "seiche_gap_by_regime.png"), seiche_fig,
       width = 13, height = 6.3, dpi = 300, bg = "white")

write_csv(regime_summary, paste0(OUT_DIR, "seiche_regime_summary.csv"))

cat("\nDate range in joined data:", date_min, "to", date_max, "\n")
cat("Total paired days:", n_total_days, "| windy days shown:", n_windy_days, "\n")
cat("All outputs saved to:", OUT_DIR, "\n")

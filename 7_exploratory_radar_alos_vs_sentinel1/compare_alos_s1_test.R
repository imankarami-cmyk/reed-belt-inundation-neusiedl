# ============================================================
# ALOS PALSAR-1 (L-band) vs SENTINEL-1 (C-band):
# which detects water under reed canopy better?
# Exploratory comparison (not part of the thesis).
#
# FAIRNESS DESIGN
# - The sensors do not overlap in time (2007-2010 vs 2015-2025), so they
#   cannot be compared on the same observations.
# - Both are scored on the same 600 pixels against the same gauge-elevation
#   label, with the same metric (AUC), WITHOUT fitting anything
#   (single-variable AUC), so neither sensor gets a modelling advantage.
# - Each is also compared with the Landsat mNDWI+EVI model (thesis
#   coefficients, never fitted to gauge labels) in its OWN period, so that
#   "2007-2010 vs 2015-2025 was a different hydrology" is controlled for.
# - Co-pol is compared with co-pol (HH vs VV), cross with cross (HV vs VH).
# - Confidence intervals come from a bootstrap over PIXELS (observations
#   within a pixel are correlated).
# ============================================================

library(tidyverse)

RADAR_CSV   <- "<PROJECT_ROOT>/Data/RADAR_ALOS1_vs_S1_600pixels.csv"
LANDSAT_CSV <- "<PROJECT_ROOT>/result/writting/new/Data/FLATS_newpixels_1984_2025_CLEANED.csv"
GAUGE_CSV   <- "<PROJECT_ROOT>/Data/all_stations_1984_2025.csv"
OUT_DIR     <- "<PROJECT_ROOT>/result/writting/new/ALOS_vs_S1/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# "doy"   : ALOS date band is day of year (about 1-366)
# "epoch" : ALOS date band is days since 1970-01-01 (about 13000-15000)
# "none"  : no usable date band -> ALOS labels use the annual MEAN gauge level
#           (weak; the AUC for ALOS is then only indicative)
DATE_FORMAT   <- "doy"
NEAR_THR_M    <- 0.15      # "hard cases": |gauge level - pixel elevation| <= 15 cm
B_BOOT        <- 300       # bootstrap repeats
S1_MAX_ROWS   <- 150000    # random subsample of S1 rows to keep runtime reasonable
set.seed(1)

# ── 1. LOAD ─────────────────────────────────────────────────
rad <- read_csv(RADAR_CSV, show_col_types = FALSE,
                col_types = cols(date = col_character())) %>%
  rename(station = nearest_station) %>%
  mutate(station = as.character(station)) %>%
  filter(is.finite(co_dB), is.finite(cross_dB))

gauge <- read_csv(GAUGE_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date), level_mua = as.numeric(level_mua),
         station = as.character(station)) %>%
  filter(!is.na(date), !is.na(level_mua)) %>%
  distinct(date, station, .keep_all = TRUE)

cat("In radar file but not in gauge file:",
    setdiff(unique(rad$station), unique(gauge$station)), "\n\n")

# ── 2. OBSERVATION DATES ────────────────────────────────────
rad <- rad %>%
  mutate(obs_date = case_when(
    sensor == "S1" ~ as.Date(date),
    sensor == "ALOS1" & DATE_FORMAT == "doy" & date_band > 0 ~
      as.Date(paste0(year, "-01-01")) + date_band - 1,
    sensor == "ALOS1" & DATE_FORMAT == "epoch" & date_band > 0 ~
      as.Date(date_band, origin = "1970-01-01"),
    TRUE ~ as.Date(NA)))

# ── 3. GAUGE-ELEVATION LABEL ────────────────────────────────
annual_mean <- gauge %>% mutate(year = year(date)) %>%
  group_by(station, year) %>% summarise(level_year = mean(level_mua), .groups = "drop")

rad <- rad %>%
  left_join(gauge %>% select(obs_date = date, station, level_day = level_mua),
            by = c("obs_date", "station")) %>%
  left_join(annual_mean, by = c("station", "year")) %>%
  mutate(level_mua = coalesce(level_day, level_year),   # annual mean only if no date
         label_src = if_else(!is.na(level_day), "daily gauge", "annual-mean gauge"),
         margin = level_mua - elevation,
         water_present = if_else(margin >= 0, 1L, 0L),
         season = case_when(is.na(obs_date) ~ NA_character_,
                            month(obs_date) %in% 4:10 ~ "growing",
                            TRUE ~ "leafoff")) %>%
  filter(!is.na(margin))

# S1: remove geometry effects, same "own history" idea as the thesis decile labels
rad <- rad %>%
  group_by(sensor, pixel_id, rel_orbit, pass) %>%
  mutate(co_anom = co_dB - median(co_dB), cross_anom = cross_dB - median(cross_dB)) %>%
  ungroup()

s1_all <- filter(rad, sensor == "S1")
if (nrow(s1_all) > S1_MAX_ROWS) s1_all <- slice_sample(s1_all, n = S1_MAX_ROWS)
rad <- bind_rows(s1_all, filter(rad, sensor == "ALOS1"))

cat("Radar observations with label:\n")
print(count(rad, sensor, label_src))
cat("\n")

# ── 4. LANDSAT BENCHMARK (thesis model, thesis exclusion rule) ──
ls <- read_csv(LANDSAT_CSV, show_col_types = FALSE) %>%
  mutate(date = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d")),
         station = as.character(station)) %>%
  filter(is.finite(mNDWI), is.finite(EVI)) %>%
  group_by(pixel_id) %>%
  mutate(p10 = quantile(mNDWI, 0.10), p90 = quantile(mNDWI, 0.90),
         is_train = mNDWI >= p90 | mNDWI <= p10) %>%       # rows the model was trained on
  ungroup() %>%
  left_join(gauge %>% select(date, station, level_mua), by = c("date", "station")) %>%
  filter(!is.na(level_mua), !is.na(elevation)) %>%
  mutate(P = plogis(1.917 + 7.961 * mNDWI + 0.105 * EVI),   # deployed coefficients, no refit
         margin = level_mua - elevation,
         water_present = if_else(margin >= 0, 1L, 0L),
         season = if_else(month(date) %in% 4:10, "growing", "leafoff"))

bench <- list(ALOS1 = filter(ls, year(date) %in% 2007:2010),
              S1    = filter(ls, year(date) %in% 2015:2025))

# ── 5. FAST AUC + PIXEL-CLUSTER BOOTSTRAP ───────────────────
fast_auc <- function(score, y) {
  n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(sum(y == 0))  # numeric: n1*n0 overflows as integer above ~2.1e9
  if (n1 < 5 || n0 < 5) return(NA_real_)
  (sum(rank(score)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

boot_auc <- function(df, vars, B = B_BOOT) {
  y <- df$water_present
  signs <- map_dbl(vars, ~ sign(mean(df[[.x]][y == 1]) - mean(df[[.x]][y == 0])))
  dir_auc <- function(k, ii) {
    a <- fast_auc(df[[vars[k]]][ii], y[ii]); if (signs[k] < 0) 1 - a else a
  }
  est <- map_dbl(seq_along(vars), ~ dir_auc(.x, seq_len(nrow(df))))
  by_pix <- split(seq_len(nrow(df)), df$pixel_id)
  pix <- names(by_pix)
  boots <- replicate(B, {
    ii <- unlist(by_pix[sample(pix, length(pix), replace = TRUE)], use.names = FALSE)
    map_dbl(seq_along(vars), ~ dir_auc(.x, ii))
  })
  boots <- matrix(boots, nrow = length(vars), dimnames = list(vars, NULL))
  list(est = setNames(est, vars), boots = boots, n = nrow(df), sign = setNames(signs, vars))
}

# NA-safe 95% interval: bootstrap repeats where a class vanished are skipped.
# If fewer than 30 repeats are usable, the interval is NA ("cannot tell").
ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 30) c(NA_real_, NA_real_) else quantile(x, c(0.025, 0.975), names = FALSE)
}

usable <- function(d) nrow(d) >= 50 && sum(d$water_present == 1) >= 5 && sum(d$water_present == 0) >= 5

subsets <- list(
  all      = function(d) d,
  near_thr = function(d) filter(d, abs(margin) <= NEAR_THR_M),
  growing  = function(d) filter(d, season == "growing"),
  leafoff  = function(d) filter(d, season == "leafoff")
)

results <- list(); store <- list()
for (sen in c("ALOS1", "S1")) {
  vars <- if (sen == "S1") c("co_dB", "cross_dB", "co_anom", "cross_anom") else c("co_dB", "cross_dB")
  for (ss in names(subsets)) {
    r <- subsets[[ss]](filter(rad, sensor == sen))
    o <- subsets[[ss]](bench[[sen]])
    if (!usable(r) || !usable(o)) next
    br <- boot_auc(r, vars)
    bo_heldout <- boot_auc(filter(o, !is_train), "P")
    bo_all     <- boot_auc(o, "P")
    store[[paste(sen, ss)]] <- list(rad = br, ben = bo_heldout)
    results[[paste(sen, ss)]] <- bind_rows(
      tibble(sensor = sen, subset = ss, variable = vars, AUC = br$est,
             lo = apply(br$boots, 1, function(x) ci(x)[1]),
             hi = apply(br$boots, 1, function(x) ci(x)[2]), n = br$n,
             sign = br$sign, valid_boots = rowSums(is.finite(br$boots))),
      tibble(sensor = sen, subset = ss, variable = "Landsat mNDWI+EVI (thesis held-out rows)",
             AUC = bo_heldout$est, lo = ci(bo_heldout$boots)[1], hi = ci(bo_heldout$boots)[2],
             n = bo_heldout$n, sign = NA_real_,
             valid_boots = sum(is.finite(bo_heldout$boots))),
      tibble(sensor = sen, subset = ss, variable = "Landsat mNDWI+EVI (all rows)",
             AUC = bo_all$est, lo = ci(bo_all$boots)[1], hi = ci(bo_all$boots)[2],
             n = bo_all$n, sign = NA_real_,
             valid_boots = sum(is.finite(bo_all$boots))))
  }
}

tab <- bind_rows(results) %>% mutate(across(c(AUC, lo, hi), ~ round(.x, 3)))
cat("=== AUC vs gauge-elevation label (95% pixel-bootstrap CI) ===\n")
cat("sign: +1 = higher backscatter when flooded (double-bounce), -1 = lower (specular)\n\n")
print(tab, n = Inf)
write_csv(tab, paste0(OUT_DIR, "auc_table.csv"))

# ── 6. THE ACTUAL COMPARISON ────────────────────────────────
dq <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 30) rep(NA_real_, 3)
  else round(quantile(x, c(0.5, 0.025, 0.975), names = FALSE), 3)
}
cmp <- list()
for (ss in names(subsets)) {
  a <- store[[paste("ALOS1", ss)]]; s <- store[[paste("S1", ss)]]
  if (is.null(a) || is.null(s)) next
  for (v in c("co_dB", "cross_dB")) {
    d_raw <- s$rad$boots[v, ] - a$rad$boots[v, ]
    d_rel <- (s$rad$boots[v, ] - s$ben$boots[1, ]) - (a$rad$boots[v, ] - a$ben$boots[1, ])
    cmp[[paste(ss, v)]] <- tibble(
      subset = ss, variable = v,
      diff_S1_minus_ALOS = dq(d_raw)[1], lo1 = dq(d_raw)[2], hi1 = dq(d_raw)[3],
      rel_to_Landsat_S1_minus_ALOS = dq(d_rel)[1], lo2 = dq(d_rel)[2], hi2 = dq(d_rel)[3])
  }
}
cmp <- bind_rows(cmp)
cat("\n=== S1 minus ALOS (positive = Sentinel-1 better; CI excluding 0 = clear difference) ===\n")
cat("diff_S1_minus_ALOS            : raw AUC difference\n")
cat("rel_to_Landsat_S1_minus_ALOS  : each sensor first measured against Landsat in its OWN period\n\n")
print(cmp, n = Inf)
write_csv(cmp, paste0(OUT_DIR, "sensor_comparison.csv"))

# ── 7. OPTIONAL: leave-one-area-out logistic models (favours radar) ──
# Radar models are fitted on gauge labels; the Landsat formula never was.
loso <- function(df, formula) {
  map_dfr(unique(df$station), function(st) {
    tr <- filter(df, station != st); te <- filter(df, station == st)
    if (!usable(te)) return(tibble(held_out = st, AUC = NA_real_))
    fit <- tryCatch(glm(formula, data = tr, family = binomial()), error = function(e) NULL)
    if (is.null(fit)) return(tibble(held_out = st, AUC = NA_real_))
    tibble(held_out = st, AUC = fast_auc(predict(fit, te, type = "response"), te$water_present))
  })
}
lo_tab <- bind_rows(
  loso(filter(rad, sensor == "ALOS1"), water_present ~ co_dB + cross_dB) %>% mutate(model = "ALOS1: HH + HV"),
  loso(filter(rad, sensor == "S1"),    water_present ~ co_dB + cross_dB) %>% mutate(model = "S1: VV + VH"),
  loso(filter(rad, sensor == "S1"),    water_present ~ co_anom + cross_anom) %>% mutate(model = "S1: VV + VH (own-median anomaly)")
) %>% group_by(model) %>%
  summarise(AUC_mean = round(mean(AUC, na.rm = TRUE), 3),
            AUC_sd = round(sd(AUC, na.rm = TRUE), 3),
            folds = sum(!is.na(AUC)), .groups = "drop")
cat("\n=== Leave-one-area-out logistic AUC (fitted on gauge labels, so optimistic) ===\n")
print(lo_tab)
write_csv(lo_tab, paste0(OUT_DIR, "loso_table.csv"))

cat("\nOutputs saved to:", OUT_DIR, "\n")

# ============================================================
# RANKING YEARS BY MEAN ANNUAL GAUGE WATER LEVEL
# Lake Neusiedl, Austria | 1984-2025
# ------------------------------------------------------------
# Reads every yearly water_mean_level/*-<year>.txt file, computes
# each year's mean annual water level, and ranks all years from
# lowest to highest. Use this to confirm (or correct) which years
# actually qualify as "driest" and "wettest" before naming them
# in the thesis - don't rely on assumption or memory.
# ============================================================

library(tidyverse)

WATERLEVEL_DIR <- "<PROJECT_ROOT>/Data/water_mean_level"
OUT_DIR        <- "<PROJECT_ROOT>/result/writting/GEE/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

START_YEAR <- 1984
END_YEAR   <- 2025

# ------------------------------------------------------------
# 1. READ + AGGREGATE EVERY YEARLY FILE
# ------------------------------------------------------------

read_year_level <- function(year) {
  pattern <- paste0("-", year, "\\.txt$")
  f <- list.files(WATERLEVEL_DIR, pattern = pattern, full.names = TRUE)
  if (length(f) == 0) {
    warning(paste("No water level file found for", year))
    return(tibble(year = year, mean_level = NA_real_, min_level = NA_real_,
                  max_level = NA_real_, n_days = 0))
  }
  if (length(f) > 1) {
    warning(paste("Multiple files matched for", year, "- using the first:", f[1]))
    f <- f[1]
  }
  df <- read.table(f, header = FALSE, col.names = c("date", "level"),
                    colClasses = c("character", "numeric"))
  tibble(
    year       = year,
    mean_level = mean(df$level, na.rm = TRUE),
    min_level  = min(df$level, na.rm = TRUE),
    max_level  = max(df$level, na.rm = TRUE),
    n_days     = nrow(df)
  )
}

waterlevel <- map_dfr(START_YEAR:END_YEAR, read_year_level)

cat("Loaded", sum(!is.na(waterlevel$mean_level)), "of",
    length(START_YEAR:END_YEAR), "years.\n")
if (any(is.na(waterlevel$mean_level))) {
  cat("WARNING - missing years (not included in ranking):",
      paste(waterlevel$year[is.na(waterlevel$mean_level)], collapse = ", "), "\n")
}

# Flag any year with unusually low day coverage - a "driest year" ranking
# based on an incomplete year would be misleading
waterlevel <- waterlevel %>%
  mutate(days_flag = if_else(n_days < 300, "INCOMPLETE YEAR - check this", ""))

if (any(waterlevel$days_flag != "" & !is.na(waterlevel$mean_level))) {
  cat("\nWARNING - years with fewer than 300 days of gauge data (ranking may be unreliable):\n")
  print(waterlevel %>% filter(days_flag != "") %>% select(year, n_days))
}

# ------------------------------------------------------------
# 2. RANK ALL YEARS, DRIEST TO WETTEST
# ------------------------------------------------------------

ranked <- waterlevel %>%
  filter(!is.na(mean_level)) %>%
  arrange(mean_level) %>%
  mutate(rank_driest_to_wettest = row_number())

write_csv(ranked, paste0(OUT_DIR, "annual_water_level_ranked.csv"))

cat("\n================ 10 DRIEST YEARS (lowest mean water level) ================\n")
print(ranked %>% slice_head(n = 10) %>%
        select(rank_driest_to_wettest, year, mean_level, min_level, max_level, n_days))

cat("\n================ 10 WETTEST YEARS (highest mean water level) ================\n")
print(ranked %>% slice_tail(n = 10) %>% arrange(desc(mean_level)) %>%
        select(year, mean_level, min_level, max_level, n_days))

cat("\nFull ranked table saved to:",
    paste0(OUT_DIR, "annual_water_level_ranked.csv"), "\n")

# ------------------------------------------------------------
# 3. QUICK VISUAL CHECK - ALL YEARS RANKED, LOWEST TO HIGHEST
# ------------------------------------------------------------

fig <- ggplot(ranked, aes(x = reorder(as.character(year), mean_level), y = mean_level)) +
  geom_col(fill = "#33465c") +
  geom_hline(yintercept = mean(ranked$mean_level), linetype = "dashed", color = "grey50") +
  coord_flip() +
  labs(title = "Mean Annual Water Level by Year (ranked, lowest to highest)",
       subtitle = "Dashed line = 1984-2025 overall mean",
       x = NULL, y = "Mean annual water level (m a.s.l.)") +
  theme_minimal(base_size = 9)

ggsave(paste0(OUT_DIR, "annual_water_level_ranked.png"), fig,
       width = 7, height = 12, dpi = 300, bg = "white")

cat("Ranked bar chart saved to:", paste0(OUT_DIR, "annual_water_level_ranked.png"), "\n")

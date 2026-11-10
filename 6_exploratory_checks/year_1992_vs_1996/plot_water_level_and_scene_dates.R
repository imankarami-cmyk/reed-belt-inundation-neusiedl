# ============================================================
# FIGURE: Daily water level 1992 vs 1996
# with actual Landsat scene dates from full wetland GEE export
# ============================================================

library(tidyverse)

DATA_DIR    <- "<PROJECT_ROOT>/Data/water_mean_level/"
SCENE_1992  <- "<PROJECT_ROOT>/Data/scene_dates_1992_full_wetland.csv"
SCENE_1996  <- "<PROJECT_ROOT>/Data/scene_dates_1996_full_wetland.csv"
OUT_DIR     <- "<PROJECT_ROOT>/result/writting/new/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_1992 <- "#2a78d6"
COL_1996 <- "#1d9e75"

# ── 1. LOAD WATER LEVEL DATA ────────────────────────────────
read_wl <- function(year) {
  files <- list.files(DATA_DIR,
                      pattern = paste0(".*", year, "\\.txt"),
                      full.names = TRUE)
  if (length(files) == 0) stop(paste("No file for year", year))
  read_table(files[1], col_names = c("date", "level_m"),
             col_types = "Dd") %>%
    mutate(date = as.Date(date), year = as.integer(year))
}

wl_1992 <- read_wl(1992)
wl_1996 <- read_wl(1996)
wl_both <- bind_rows(wl_1992, wl_1996)

# Annual means
mean_1992 <- round(mean(wl_1992$level_m), 4)
mean_1996 <- round(mean(wl_1996$level_m), 4)

# ── 2. LOAD SCENE DATES AND JOIN WATER LEVEL ────────────────
load_scenes <- function(path, yr) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(date = as.Date(date), year = as.integer(yr)) %>%
    distinct(date, year, month, doy)
}

sc_1992 <- load_scenes(SCENE_1992, 1992) %>%
  left_join(wl_1992 %>% select(date, level_m), by = "date")

sc_1996 <- load_scenes(SCENE_1996, 1996) %>%
  left_join(wl_1996 %>% select(date, level_m), by = "date")

scenes <- bind_rows(sc_1992, sc_1996)

# Mean WL on scene dates
scene_mean_1992 <- round(mean(sc_1992$level_m, na.rm=TRUE), 4)
scene_mean_1996 <- round(mean(sc_1996$level_m, na.rm=TRUE), 4)

cat("=== KEY NUMBERS ===\n")
cat("Annual mean 1992:", mean_1992, "m |",
    "Annual mean 1996:", mean_1996, "m\n")
cat("Mean WL on 1992 scenes:", scene_mean_1992, "m |",
    "Mean WL on 1996 scenes:", scene_mean_1996, "m\n")
cat("Satellite flood frequency — 1992: 64.0% | 1996: 49.7%\n\n")
cat("KEY FINDING: 1996 scenes saw", round(scene_mean_1996-scene_mean_1992,3),
    "m HIGHER water than 1992 — yet lower flood frequency.\n")
cat("This shows canopy obscuration matters more than water level alone.\n")

# ── 3. PLOT ─────────────────────────────────────────────────
wl_plot <- wl_both %>%
  mutate(doy = yday(date), year_fac = factor(year))

scenes_plot <- scenes %>%
  mutate(doy = yday(date), year_fac = factor(year))

y_min <- min(wl_both$level_m) - 0.01
y_max <- max(wl_both$level_m) + 0.02

p <- ggplot() +

  # Spring shading
  annotate("rect", xmin=60, xmax=181,
           ymin=y_min, ymax=y_max,
           fill="lightblue", alpha=0.10) +
  annotate("text", x=120, y=y_max-0.003,
           label="Spring", color="grey50",
           size=3.2, fontface="italic") +

  # Annual mean dotted lines
  geom_hline(yintercept=mean_1992, color=COL_1992,
             linetype="dotted", linewidth=0.7) +
  geom_hline(yintercept=mean_1996, color=COL_1996,
             linetype="dotted", linewidth=0.7) +

  # Mean on scene dates — dashed lines
  geom_hline(yintercept=scene_mean_1992, color=COL_1992,
             linetype="dashed", linewidth=0.9) +
  geom_hline(yintercept=scene_mean_1996, color=COL_1996,
             linetype="dashed", linewidth=0.9) +

  # Daily water level curves
  geom_line(data=wl_plot,
            aes(x=doy, y=level_m, color=year_fac),
            linewidth=1.0) +

  # Scene date circles
  geom_point(data=scenes_plot,
             aes(x=doy, y=level_m, color=year_fac),
             size=4.5, shape=21, fill="white", stroke=2.0) +

  # Labels — annual means (left margin)
  annotate("text", x=3, y=mean_1992+0.005,
           label=paste0("Annual mean: ", mean_1992, " m"),
           color=COL_1992, size=2.8, hjust=0, fontface="italic") +
  annotate("text", x=3, y=mean_1996-0.005,
           label=paste0("Annual mean: ", mean_1996, " m"),
           color=COL_1996, size=2.8, hjust=0, fontface="italic") +

  # Labels — scene means (right margin)
  annotate("label", x=340, y=scene_mean_1992,
           label=paste0("Mean on 1992\nscene dates:\n",
                        scene_mean_1992, " m"),
           color=COL_1992, fill="white", size=2.8,
           label.size=0.2, label.padding=unit(0.25,"lines"),
           hjust=0) +
  annotate("label", x=340, y=scene_mean_1996,
           label=paste0("Mean on 1996\nscene dates:\n",
                        scene_mean_1996, " m"),
           color=COL_1996, fill="white", size=2.8,
           label.size=0.2, label.padding=unit(0.25,"lines"),
           hjust=0) +

  # Key annotation box
  annotate("label",
           x=180, y=y_min+0.025,
           label=paste0(
             "1996 scenes captured water ",
             round(scene_mean_1996-scene_mean_1992, 2),
             " m HIGHER than 1992 scenes\n",
             "Yet satellite flood frequency: 1992 = 64.0%  vs  1996 = 49.7%\n",
             "→ Higher water level does not guarantee stronger spectral flood signal\n",
             "  when water is hidden beneath dense reed canopy"
           ),
           color="grey20", fill="#fff9e6", size=3.0,
           label.size=0.3, label.padding=unit(0.4,"lines"),
           hjust=0.5) +

  scale_color_manual(
    values=c("1992"=COL_1992, "1996"=COL_1996),
    name=NULL,
    labels=c(
      paste0("1992  ·  ", nrow(sc_1992),
             " Landsat scenes  ·  flood frequency: 64.0%"),
      paste0("1996  ·  ", nrow(sc_1996),
             " Landsat scenes  ·  flood frequency: 49.7%")
    )
  ) +

  scale_x_continuous(
    breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
    labels=c("Jan","Feb","Mar","Apr","May","Jun",
             "Jul","Aug","Sep","Oct","Nov","Dec"),
    expand=c(0.01,0),
    limits=c(1, 380)
  ) +
  scale_y_continuous(
    limits=c(y_min, y_max),
    breaks=scales::pretty_breaks(n=8)
  ) +

  labs(
    title="Daily water level 1992 vs 1996 — Lake Neusiedl mean gauge",
    subtitle=paste0(
      "Circles = Landsat acquisition dates across the full reed belt  ·  ",
      "Dashed = mean WL on scene dates  ·  Dotted = annual mean\n",
      "1996 was wetter by every measure — yet the satellite detected",
      " more flooding in 1992. The explanation lies in canopy obscuration."
    ),
    x=NULL,
    y="Lake mean water level (m a.s.l.)"
  ) +

  theme_minimal(base_size=13) +
  theme(
    plot.title       = element_text(face="bold", size=13),
    plot.subtitle    = element_text(color="grey40", size=10,
                                    margin=margin(b=8)),
    legend.position  = "bottom",
    legend.text      = element_text(size=11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color="grey90", linewidth=0.4),
    axis.line.x      = element_line(color="grey70"),
    plot.margin      = margin(10, 80, 10, 10)
  )

OUT_FILE <- paste0(OUT_DIR, "water_level_1992_vs_1996_FINAL.png")
ggsave(OUT_FILE, p, width=14, height=7, dpi=300, bg="white")
cat("\nSaved:", OUT_FILE, "\n")

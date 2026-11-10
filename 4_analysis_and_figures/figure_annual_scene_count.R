# ============================================================
# FIGURE: Annual Landsat scene count 1984–2025
# Shows how many cloud-free scenes were available per year
# to compute the flood frequency maps.
# ============================================================

library(tidyverse)

IN_CSV  <- "<PROJECT_ROOT>/result/writting/FloodMap_annual_timeseries.csv"
OUT_DIR <- "<PROJECT_ROOT>/result/writting/new/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df <- read_csv(IN_CSV, show_col_types = FALSE)

# Confirm the scene count column name
cat("Columns:", names(df), "\n")
cat("Years:", range(df$year), "\n")
cat("Scene count range:", range(df$n_scenes), "\n")

p <- ggplot(df, aes(x = year, y = n_scenes)) +
  geom_col(width = 0.85, fill = "#888780", color = "white", linewidth = 0.2) +

  # Highlight 1996 and 2022 with annotation
  geom_vline(xintercept = 1996, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 2022, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  annotate("text", x = 1996, y = max(df$n_scenes) * 0.97,
           label = "1996\n(32 scenes)", hjust = 1.1,
           size = 3.2, color = "grey30", fontface = "italic") +
  annotate("text", x = 2022, y = max(df$n_scenes) * 0.97,
           label = "2022\n(135 scenes)", hjust = -0.1,
           size = 3.2, color = "grey30", fontface = "italic") +

  scale_x_continuous(breaks = seq(1985, 2025, by = 5),
                     expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, max(df$n_scenes) * 1.08)) +

  labs(
    title    = "Number of cloud-free Landsat scenes per year (1984–2025)",
    subtitle = "1996 and 2022 highlighted — the wettest and driest years compared in Section 5.8",
    x        = "Year",
    y        = "Number of cloud-free scenes"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(color = "grey40", size = 10.5,
                                      margin = margin(b = 8)),
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
    axis.line.x        = element_line(color = "grey70", linewidth = 0.4),
    plot.margin        = margin(10, 15, 10, 10)
  )

OUT_FILE <- paste0(OUT_DIR, "annual_scene_count_1984_2025.png")
ggsave(OUT_FILE, p, width = 12, height = 5.5, dpi = 300, bg = "white")
cat("\nSaved:", OUT_FILE, "\n")

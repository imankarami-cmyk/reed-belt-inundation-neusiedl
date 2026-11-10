# ============================================================
# FIGURE: Per-pixel wet/dry labeling — ALL individual observations
# Pixel: Rust_1_2_15 | Station: Rust | Elevation: 115.617 m
# ============================================================

library(tidyverse)

CSV_PATH  <- "<PROJECT_ROOT>/Data/FLATS_newpixels_1984_2025_CLEANED.csv"
TARGET_ID <- "Rust_1_2_15"
INDEX_COL <- "mNDWI"    # change to "NDWI" for the baseline version
OUT_DIR   <- "<PROJECT_ROOT>/result/writting/new/"
OUT_FILE  <- paste0(OUT_DIR, "slide12_pixel_labeling_", INDEX_COL, "_allobs.png")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

COL_FLOODED <- "#2a78d6"
COL_DRY     <- "#eb6834"
COL_AMBIG   <- "#b4b2a9"

# ------------------------------------------------------------
# 1. LOAD
# ------------------------------------------------------------

pix <- read_csv(CSV_PATH, show_col_types = FALSE) %>%
  filter(pixel_id == TARGET_ID)

if (nrow(pix) == 0) {
  cat("Pixel not found. Available Rust_1_2_* IDs:\n")
  read_csv(CSV_PATH, show_col_types = FALSE) %>%
    filter(str_detect(pixel_id, "Rust_1_2")) %>%
    distinct(pixel_id) %>% print()
  stop("Fix TARGET_ID and re-run.")
}

cat("Found:", nrow(pix), "observations for", TARGET_ID, "\n")

pix <- pix %>%
  mutate(
    date_parsed = as.Date(date, tryFormats = c("%m/%d/%Y", "%Y-%m-%d")),
    index_val   = .data[[INDEX_COL]]
  ) %>%
  filter(!is.na(index_val), !is.na(date_parsed)) %>%
  arrange(date_parsed)

# ------------------------------------------------------------
# 2. THRESHOLDS from this pixel's own raw history (Section 4.5.1)
# ------------------------------------------------------------

p10 <- quantile(pix$index_val, 0.10, na.rm = TRUE)
p90 <- quantile(pix$index_val, 0.90, na.rm = TRUE)

cat(INDEX_COL, "10th pct:", round(p10, 4), " | 90th pct:", round(p90, 4), "\n")

pix <- pix %>%
  mutate(
    label = case_when(
      index_val >= p90 ~ "Flooded",
      index_val <= p10 ~ "Dry",
      TRUE             ~ "Ambiguous"
    ),
    label = factor(label, levels = c("Flooded", "Ambiguous", "Dry"))
  )

cat("Flooded:", sum(pix$label=="Flooded"),
    "| Dry:", sum(pix$label=="Dry"),
    "| Ambiguous:", sum(pix$label=="Ambiguous"),
    "| Total:", nrow(pix), "\n")

# ------------------------------------------------------------
# 3. PLOT
# ------------------------------------------------------------

y_pad <- (max(pix$index_val) - min(pix$index_val)) * 0.15
y_lo  <- min(pix$index_val) - y_pad
y_hi  <- max(pix$index_val) + y_pad

p <- ggplot(pix, aes(x = date_parsed, y = index_val)) +

  # Event bands
  annotate("rect", xmin=as.Date("2000-01-01"), xmax=as.Date("2004-12-31"),
           ymin=y_lo, ymax=y_hi, fill=COL_DRY, alpha=0.07) +
  annotate("text", x=as.Date("2002-06-01"), y=y_hi-(y_hi-y_lo)*0.05,
           label="2000–04 crisis", color=COL_DRY, size=3.2, hjust=0.5, fontface="italic") +

  annotate("rect", xmin=as.Date("2022-01-01"), xmax=as.Date("2025-12-31"),
           ymin=y_lo, ymax=y_hi, fill=COL_DRY, alpha=0.07) +
  annotate("text", x=as.Date("2023-07-01"), y=y_hi-(y_hi-y_lo)*0.05,
           label="2022–25\ndrought", color=COL_DRY, size=3.2, hjust=0.5, fontface="italic") +

  # Ambiguous band
  annotate("rect", xmin=as.Date("1984-01-01"), xmax=as.Date("2026-06-01"),
           ymin=p10, ymax=p90, fill=COL_AMBIG, alpha=0.12) +
  annotate("text", x=as.Date("2006-01-01"), y=(p10+p90)/2,
           label="Middle 80%  —  excluded from training",
           color="grey45", size=3.3, fontface="italic") +

  # Threshold lines
  geom_hline(yintercept=p90, color=COL_FLOODED, linetype="dashed", linewidth=0.75) +
  geom_hline(yintercept=p10, color=COL_DRY,     linetype="dashed", linewidth=0.75) +

  # Threshold labels
  annotate("label", x=as.Date("2026-08-01"), y=p90,
           label=paste0("90th pct: ", round(p90,3), "\n→ labeled FLOODED"),
           color=COL_FLOODED, fill="white", size=3, hjust=0,
           label.size=0.3, label.padding=unit(0.25,"lines")) +
  annotate("label", x=as.Date("2026-08-01"), y=p10,
           label=paste0("10th pct: ", round(p10,3), "\n→ labeled DRY"),
           color=COL_DRY, fill="white", size=3, hjust=0,
           label.size=0.3, label.padding=unit(0.25,"lines")) +

  # All individual points
  geom_point(aes(fill=label, size=label, alpha=label),
             shape=21, color="white", stroke=0.5) +

  scale_fill_manual(values=c("Flooded"=COL_FLOODED,"Dry"=COL_DRY,"Ambiguous"=COL_AMBIG), name=NULL) +
  scale_size_manual(values=c("Flooded"=3.5,"Dry"=3.5,"Ambiguous"=1.8), name=NULL) +
  scale_alpha_manual(values=c("Flooded"=1,"Dry"=1,"Ambiguous"=0.5), name=NULL) +

  scale_x_date(date_breaks="5 years", date_labels="%Y",
               limits=c(as.Date("1984-01-01"), as.Date("2030-01-01")),
               expand=c(0,0)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=7), limits=c(y_lo, y_hi)) +

  labs(
    title    = paste0("Pixel ", TARGET_ID, "  ·  All ", nrow(pix),
                      " observations  ·  ", INDEX_COL, "  (1984–2025)"),
    subtitle = paste0("Elevation: ", round(unique(pix$elevation)[1],3),
                      " m a.s.l.  ·  Labeled: ",
                      sum(pix$label=="Flooded"), " flooded  |  ",
                      sum(pix$label=="Dry"), " dry  |  ",
                      sum(pix$label=="Ambiguous"), " ambiguous (excluded from training)"),
    x = "Date",
    y = paste0(INDEX_COL, " value")
  ) +

  theme_minimal(base_size=13) +
  theme(
    plot.title         = element_text(face="bold", size=13.5, margin=margin(b=4)),
    plot.subtitle      = element_text(color="grey40", size=10.5, margin=margin(b=10)),
    legend.position    = "bottom",
    legend.text        = element_text(size=11),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color="grey90", linewidth=0.4),
    axis.line          = element_line(color="grey70", linewidth=0.4),
    plot.margin        = margin(10, 95, 10, 10)
  ) +
  guides(
    fill  = guide_legend(override.aes=list(size=4.5, alpha=1)),
    size  = "none",
    alpha = "none"
  )

ggsave(OUT_FILE, p, width=14, height=6.5, dpi=300, bg="white")
cat("\nSaved:", OUT_FILE, "\n")
cat("Tip: set INDEX_COL <- 'NDWI' and re-run for the baseline comparison.\n")

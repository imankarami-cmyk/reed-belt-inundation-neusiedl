# Exploratory test: ALOS PALSAR-1 (L-band) vs Sentinel-1 (C-band)

> **Status: exploratory. Not part of the thesis results.** This folder answers one side question and does not feed into any thesis table or figure.

## Question

Can radar backscatter detect water under dense reed canopy, and if so, which sensor does it better: L-band (ALOS PALSAR-1) or C-band (Sentinel-1)? The optical model of the thesis (Landsat mNDWI + EVI) is used as a benchmark.

## Why a head-to-head on the same observations is not possible

The two sensors never overlap in time (ALOS-1 mosaics: 2007–2010; Sentinel-1: 2015–2025). The comparison therefore uses a design that keeps it as fair as the data allow:

- Both sensors are scored on the **same 600 thesis pixels** against the **same gauge-elevation label** (flooded if the gauge level is at or above the pixel elevation on that date).
- The metric is the **AUC of a single radar variable, with no model fitting**, so neither sensor gets a modelling advantage. Co-polarisation is compared with co-polarisation (HH vs VV) and cross with cross (HV vs VH).
- Each sensor is also compared with the **Landsat mNDWI+EVI model in its own period**, using the thesis coefficients (never refitted to gauge labels). This controls for 2007–2010 and 2015–2025 being hydrologically different.
- Confidence intervals are 95% bootstrap intervals **over pixels** (300 repeats), because observations within a pixel are correlated.
- Subsets: all dates, "hard cases" (gauge level within 15 cm of the pixel elevation), growing season (Apr–Oct), leaf-off (Nov–Mar).
- For Sentinel-1 only, an "own-median anomaly" (value minus that pixel's median for the same orbit and pass) isolates the within-pixel change from between-pixel differences.

## Files

| File | Purpose |
|---|---|
| `compare_alos_s1_extract.js` | GEE script. Samples ALOS-1 and Sentinel-1 at the 600 pixels and exports one CSV with common columns (`co_dB` = HH or VV, `cross_dB` = HV or VH) |
| `compare_alos_s1_test.R` | R script. Joins the gauge label, computes AUC with pixel-bootstrap intervals, compares the two sensors and the Landsat benchmark, and runs an optional leave-one-area-out model |
| `phase1_alos_hh_raster_export.ipynb` | The original Phase 1 notebook: exports annual ALOS HH mosaics (2007–2010) as rasters. It does not sample the pixels, and is kept for reference only |
| `results/auc_table.csv` | AUC and 95% intervals per sensor, variable and subset |
| `results/sensor_comparison.csv` | Sentinel-1 minus ALOS differences, raw and relative to Landsat |
| `results/loso_table.csv` | Leave-one-area-out logistic AUC |

## How to run

**Prerequisites.** The thesis pixel table must exist as a GEE asset (`sample_pixels_with_nearest_station`), and you need the thesis files `FLATS_newpixels_1984_2025_CLEANED.csv` and `all_stations_1984_2025.csv` (gauge table). Replace `YOUR_GEE_PROJECT` and `<PROJECT_ROOT>` (see the main README).

1. Open `compare_alos_s1_extract.js` in the GEE Code Editor and run it. Read the console printout of the ALOS mosaic bands and the value of the `date` band at one pixel.
2. Run the export task. If it times out, set `STATION` to one sampling area, or set `S1_REL_ORBIT` to one Sentinel-1 orbit track, and export in several parts.
3. Save the CSV as `<PROJECT_ROOT>/Data/RADAR_ALOS1_vs_S1_600pixels.csv`.
4. In `compare_alos_s1_test.R`, set `DATE_FORMAT` from the console printout in step 1: `"doy"` (values about 1–366), `"epoch"` (about 13000–15000, days since 1970-01-01) or `"none"` (no date band; ALOS then uses the annual mean gauge level, which is a weak label).
5. Run `compare_alos_s1_test.R`. Outputs go to `<PROJECT_ROOT>/result/writting/new/ALOS_vs_S1/`.

Other settings at the top of the R script: `NEAR_THR_M` (hard-case margin), `B_BOOT` (bootstrap repeats), `S1_MAX_ROWS` (random subsample of Sentinel-1 rows, for runtime).

## Results (final run)

AUC against the gauge-elevation label, all dates, 95% pixel-bootstrap interval. All ALOS observations fall in the growing season, so "all dates" and "growing season" are the same for ALOS.

| | AUC (95% interval) |
|---|---|
| ALOS HH | 0.712 (0.672–0.746) |
| ALOS HV | 0.717 (0.674–0.750) |
| Landsat mNDWI+EVI in the ALOS years (2007–2010) | 0.766 (0.732–0.794) |
| Sentinel-1 VV | 0.587 (0.566–0.613) |
| Sentinel-1 VH | 0.610 (0.585–0.639) |
| Sentinel-1 VV, own-median anomaly | 0.513 (0.505–0.520) |
| Sentinel-1 VH, own-median anomaly | 0.505 (0.500–0.512) |
| Landsat mNDWI+EVI in the Sentinel-1 years (2015–2025) | 0.838 (0.810–0.865) |

Sentinel-1 minus ALOS (negative means ALOS is better; an interval that excludes 0 is a clear difference):

| Subset | Co-pol (VV − HH) | Cross-pol (VH − HV) |
|---|---|---|
| All dates | −0.124 (−0.169 to −0.079) | −0.106 (−0.150 to −0.056) |
| Hard cases | −0.044 (−0.109 to 0.026) | −0.074 (−0.143 to −0.007) |

Leave-one-area-out logistic AUC (models fitted on gauge labels, six folds): ALOS HH+HV 0.763 ± 0.142; Sentinel-1 VV+VH 0.668 ± 0.105; Sentinel-1 anomaly 0.569 ± 0.061.

## Interpretation

- **Sentinel-1 C-band does not detect flooding under this canopy.** Raw AUC is 0.59–0.61, against 0.84 for Landsat in the same period. Compared with each pixel's own history it is about 0.5, which is chance.
- **L-band is clearly better than C-band**, consistent with deeper canopy penetration at longer wavelengths. It comes close to the Landsat benchmark but does not beat it. In the hard cases ALOS and Landsat are about equal within the (wide) intervals.
- **The sign is unexpected.** In ALOS, backscatter is *lower* when flooded for both HH and HV. A double-bounce under flooded reed would raise HH. A plausible but unproven explanation is that the raw signal reflects differences between pixels (habitat) rather than flooding events.

## Limitations

- ALOS has about four observations per pixel, so its AUC mostly reflects between-pixel contrast, and the within-pixel signal cannot be tested. Only 249 ALOS observations fall in the hard-case subset.
- The two sensors cover different years, seasons and resolutions (25 m vs about 30 m after smoothing 10 m pixels).
- The within-pixel (anomaly) check was done for radar only, not for mNDWI+EVI, so the comparison is not perfectly even.
- Labels carry noise from LiDAR bias in dense reed and from wind-driven water-level differences of up to about 18 cm across the lake.
- The leave-one-area-out models are fitted on gauge labels, which favours radar, because the Landsat formula was never fitted on gauge labels. Use the single-variable AUC for the verdict.
- The test of adding SAR variables to the mNDWI+EVI model was **not** run.

## Change log

- A first run left the Sentinel-1 "all dates" rows empty because of an integer overflow in the AUC helper (`n1 * n0` exceeds the integer range with 150,000 rows). Counts are now numeric.
- Bootstrap intervals skip repeats in which a class disappeared, and report the number of usable repeats (`valid_boots`).

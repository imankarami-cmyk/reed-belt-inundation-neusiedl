# Long-term remote sensing of reed-belt inundation at Lake Neusiedl

Code for my master's thesis at BOKU Vienna: *Long-Term Remote Sensing of Austrian Reed Belt Inundation Regimes at Lake Neusiedl*.

The thesis adapts the FLATS/FLATS+ approach (a water index paired with EVI in a logistic regression, originally built for tidal salt marshes) to a freshwater, non-tidal reed belt. It tests whether replacing NDWI with mNDWI improves detection of flooding hidden under dense reed canopy, and then maps inundation across 41 years of Landsat imagery (1984–2025).

> **This repository contains code only.** The input data (gauge records, LiDAR terrain model, wind records, Landsat exports) are not redistributed. See [Data](#data).

## Method in brief

1. **Sampling.** 600 pixels are drawn around six gauge stations, stratified by 10 cm elevation bands from a 5 m airborne-laser terrain model (100 pixels per station).
2. **Extraction.** The full Landsat 5/7/8/9 Collection 2 surface-reflectance time series is extracted for each pixel in Google Earth Engine (GEE). NDWI, mNDWI and EVI are computed.
3. **Labels without field ground truth.** For each pixel, observations in the top 10% of its own index history are labelled *flooded* and the bottom 10% *dry*. The middle 80% is discarded.
4. **Models.** Two logistic regressions: baseline (NDWI + EVI) and improved (mNDWI + EVI).
5. **Validation against physical data.** Predictions are scored against a gauge-elevation criterion (flooded if the gauge water level is at or above the pixel's terrain elevation). Thresholds are chosen with Youden's J. Three tests: full-dataset gauge validation, repeated cross-validation, and leave-one-sampling-area-out.
6. **Deployment.** The fixed coefficients and thresholds are applied to every pixel in the wetland and every Landsat scene in GEE, producing annual and long-term flood-frequency maps.

The headline results are in the thesis. They are deliberately not repeated here so they cannot drift out of sync with the document.

## Repository structure

```
.
├── 1_data_extraction_gee/           Pixel sampling and Landsat time-series extraction
├── 2_r_pipeline/                    Merge, clean, label, fit, validate (thesis Chapter 5)
├── 3_deployment_gee/                Apply the models to the whole wetland; zonal statistics
├── 4_analysis_and_figures/          Figures, checks, seiche (wind) analysis
├── 5_sensitivity_analyses/          Alternative designs tested during review
├── 6_exploratory_checks/            Side investigations (1992 vs 1996, shoreline/CORINE)
├── 7_exploratory_radar_alos_vs_sentinel1/   Radar test, NOT part of the thesis (own README)
└── config/station_buffers.csv       Gauge-station coordinates and buffer radii
```

## Pipeline and thesis sections

Run in this order. File names carry the thesis section they produce.

| Step | File | What it does | Thesis |
|---|---|---|---|
| 1 | `1_data_extraction_gee/pixel_sampling.ipynb` | Station buffers, 10 cm elevation bins, stratified pixel draw; exports the sampled pixels | 4.3 |
| 2 | `1_data_extraction_gee/extract_landsat_pixels.js` | Landsat 5/7/8/9 extraction at the sampled pixels (cloud masking, scale factors, NDWI/mNDWI/EVI) | 4.2, 4.4 |
| 3 | `2_r_pipeline/00_merge_sensor_exports.R` | Merges the four per-sensor exports into one pixel-date table. Set the four input file names at the top to match your GEE export names. It writes `..._merged.csv`; rename it to `FLATS_newpixels_1984_2025.csv`, which step 5 reads | 4.3 |
| 4 | `2_r_pipeline/5_1_sample_composition.R` | Sample composition tables and plots | 5.1 |
| 5 | `2_r_pipeline/5_2_predictor_distributions.R` | Physical-range filtering, predictor and per-sensor statistics; writes the `_CLEANED` table | 5.2 |
| 6 | `2_r_pipeline/5_3_decile_labeling.R` | Per-pixel decile labels and thresholds | 5.3 |
| 7 | `2_r_pipeline/5_4_model_coefficients.R` | Logistic regression coefficients | 5.4 |
| 8 | `2_r_pipeline/5_5_gauge_validation.R` | Gauge validation, ROC, Youden-J thresholds | 5.5 |
| 9 | `2_r_pipeline/5_6_repeated_cv.R` | 10-round repeated cross-validation | 5.6 |
| 10 | `2_r_pipeline/5_7_loso.R` | Leave-one-sampling-area-out validation | 5.7 |
| 11 | `3_deployment_gee/deploy_flood_maps.ipynb` | Annual and long-term flood-frequency maps, time series, GeoTIFF export | 4.7, 5.8 |
| 12 | `3_deployment_gee/zonal_stats_by_sampling_area.js` | Mean flood frequency per sampling area | 5.8 |

Supporting scripts in `4_analysis_and_figures/`:

| File | Purpose |
|---|---|
| `rank_water_level_years.R` | Ranks years by mean annual lake water level (picks wettest and driest years) |
| `figure_pixel_distribution.R` | Sampled pixels per elevation bin and sampling area |
| `figure_flood_frequency_timeseries.R` | Annual flood frequency with scene counts and water level |
| `figure_sensor_continuity.R` | Cross-sensor continuity across Landsat generations |
| `figure_annual_scene_count.R` | Number of cloud-free scenes per year |
| `figure_pixel_labeling_example.R` | Decile labelling shown for one pixel |
| `check_elevation_by_sampling_area.R` | Elevation distribution within each station buffer |
| `check_pixel_level_flood_frequency.R` | Independent recomputation of flood frequency from the 600 pixels |
| `seiche_wind_validation.R`, `seiche_three_wind_stations.R` | Wind-driven north–south water-level difference and its robustness across three wind stations |

`5_sensitivity_analyses/` holds the alternative designs examined during review (threshold chosen on decile labels vs gauge labels, an mNDWI-only model with a likelihood-ratio test, and a check of the validation row counts). They are kept for transparency; the final design is in `2_r_pipeline/`.

## Requirements

**R** (4.x): `tidyverse`, `lubridate`, `pROC`, `PRROC`, `patchwork`, `terra`, `RColorBrewer`

```r
install.packages(c("tidyverse", "lubridate", "pROC", "PRROC", "patchwork", "terra", "RColorBrewer"))
```

**Python** (notebooks): see `requirements.txt`. Versions are not pinned because the original environment was not recorded.

**Google Earth Engine:** an Earth Engine account and a Cloud project.

## Configuration before running

The scripts were written on a local Windows machine. Two placeholders replace personal values:

| Placeholder | Replace with |
|---|---|
| `<PROJECT_ROOT>` | Your local project folder (it appears inside file-path strings in the R, JS and notebook files) |
| `YOUR_GEE_PROJECT` | Your Earth Engine Cloud project ID |

Find every occurrence with:

```bash
grep -rn "<PROJECT_ROOT>\|YOUR_GEE_PROJECT" .
```

The scripts assume this layout under `<PROJECT_ROOT>`. A few scripts write to differently named result sub-folders; check the path variables at the top of each file.

```
<PROJECT_ROOT>/
├── Data/                 merged inputs (gauge table, water_mean_level/, wind data, DEM)
└── result/               script outputs (tables, figures, exports)
```

## Data

Not included. All inputs come from third parties:

| Data | Used for | Provider |
|---|---|---|
| Landsat 5/7/8/9 Collection 2 Level-2 surface reflectance | Spectral indices | USGS, accessed through Google Earth Engine |
| Airborne laser scanning terrain model, 5 m | Pixel elevation, stratified sampling | Land Burgenland survey authority |
| Daily water levels, six gauge stations and the lake-wide mean level | Validation criterion, extreme-year selection | Hydrographischer Dienst Burgenland |
| Daily wind records (Neusiedl, Illmitz, Podersdorf) | Seiche analysis | GeoSphere Austria |
| CORINE Land Cover 2018 | Study-area boundary (classes 411 and 324) | Copernicus Land Monitoring Service |

Gauge table: the R scripts read one merged file, `all_stations_1984_2025.csv`, with the columns `date`, `station` and `level_mua` (daily water level, m a.s.l.). It was assembled from the six per-station daily files outside these scripts. Station names must match the `station` column of the pixel table.

Study-area mask: the deployment uses a GEE asset named `DGM_wetland_AT` (terrain model clipped to the CORINE-derived wetland boundary). It was prepared outside the scripts in this repository.

Sampled pixels: the CSV produced by `pixel_sampling.ipynb` (with longitude, latitude, elevation and nearest station) is uploaded to GEE as a table asset named `sample_pixels_with_nearest_station`.

Please obtain the data from the providers under their terms of use.

## Reproducibility notes

- Random seeds: `set.seed(42)` in `5_4_model_coefficients.R` (train/held-out split); `set.seed(99)` in the cross-validation scripts.
- Model coefficients in `3_deployment_gee/deploy_flood_maps.ipynb` (Cell 1) were copied from the output of `5_4_model_coefficients.R`; the thresholds were copied from `5_5_gauge_validation.R`. If either script is re-run, update Cell 1 by hand.
- The scene-level cloud-cover filter is set by `CLOUD_THR` in `extract_landsat_pixels.js`.
- Scripts in `6_exploratory_checks/` and `7_exploratory_radar_alos_vs_sentinel1/` answer side questions and do not feed into the thesis results.

## Citation

> Karami, I. *Long-Term Remote Sensing of Austrian Reed Belt Inundation Regimes at Lake Neusiedl.* Master's thesis, University of Natural Resources and Life Sciences (BOKU), Vienna. [year, DOI or repository link to be added]

## License

Code is released under the MIT License (see `LICENSE`). The data are not covered by this license.

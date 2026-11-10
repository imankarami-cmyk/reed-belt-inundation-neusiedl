// ============================================================
// ZONAL STATISTICS: LONG-TERM MEAN FLOOD FREQUENCY BY STATION
// Lake Neusiedl | 1984-2025
// ------------------------------------------------------------
// PURPOSE
// Extracts the mean long-term flood frequency (both models) within
// each station's sampling buffer, from the long-term mean flood
// frequency raster already produced in Section 4.7 / GEE_corrected.
// This turns the qualitative spatial claim in Chapter 6.4 ("the
// most persistently flooded zones concentrate around Moerbisch and
// Illmitz") into an actual number per station, directly comparable
// to the EVI-based station ranking already in Table 5.5.
//
// IMPORTANT: fill in STATION_COORDS below with the exact same
// station_lon / station_lat values already used throughout your
// pixel sampling pipeline (pixel_sampling_map.ipynb) - do not
// re-derive them, since the buffers must match exactly what was
// used for pixel sampling and model validation.
// ============================================================


// ── 0. CONFIG ────────────────────────────────────────────────

// Path to the long-term mean flood frequency image, two bands:
// one per model. Replace with your actual asset/export path.
var LONGTERM_MEAN_IMG = ee.Image('projects/YOUR_GEE_PROJECT/assets/FloodMap_longterm_mean');

// Print the actual band names before assuming anything - this is what
// caught the mismatch last run (asset used generic b1/b2, not descriptive
// names). Check the console output and confirm which band is which model
// before trusting the BAND_NDWI_EVI / BAND_MNDWI_EVI assignment below.
print('Actual bands in LONGTERM_MEAN_IMG:', LONGTERM_MEAN_IMG.bandNames());

// Band names - confirmed from the asset's actual band list (b1, b2).
// IMPORTANT: verify which of b1/b2 corresponds to which model by checking
// the band order used when this image was originally exported (Section 4.7
// / GEE_corrected.ipynb) - band order there determines which is which here.
// As a default assumption, b1 = NDWI+EVI (baseline, exported first),
// b2 = mNDWI+EVI (improved, exported second) - CONFIRM before trusting output.
var BAND_NDWI_EVI  = 'b1';
var BAND_MNDWI_EVI = 'b2';

var SCALE = 30; // match native Landsat resolution, as used throughout

// Station coordinates and buffer radii - MUST MATCH Section 4.3.1 exactly.
// Fill in the real lon/lat values from your own pixel sampling script.
var STATION_COORDS = [
  {name: 'Breitenbrunn', lon: 16.770256,  lat: 47.916203, radius_m: 2500},
  {name: 'Moerbisch',    lon: 16.698516,  lat: 47.751756, radius_m: 2000},
  {name: 'Neusiedl',     lon: 16.835600,  lat: 47.929900, radius_m: 2500},
  {name: 'Rust',         lon: 16.695700,  lat: 47.801700, radius_m: 2000},
  {name: 'Illmitz',      lon: 16.752600,  lat: 47.770000, radius_m: 2500},
  {name: 'Apetlon',      lon: 16.7490326, lat: 47.6814064, radius_m: 4000}  // enlarged from 3500
];


// ── 1. BUILD STATION BUFFER FEATURES ─────────────────────────

var stationBuffers = ee.FeatureCollection(
  STATION_COORDS.map(function(s) {
    var pt = ee.Geometry.Point([s.lon, s.lat]);
    return ee.Feature(pt.buffer(s.radius_m), {station: s.name, buffer_radius_m: s.radius_m});
  })
);

print('Station buffers:', stationBuffers);
Map.centerObject(stationBuffers, 11);
Map.addLayer(stationBuffers, {color: 'red'}, 'Station buffers');
Map.addLayer(LONGTERM_MEAN_IMG.select(BAND_MNDWI_EVI), 
             {min: 0, max: 1, palette: ['red', 'orange', 'yellow', 'cyan', 'blue']}, 
             'mNDWI+EVI long-term mean');


// ── 2. ZONAL MEAN EXTRACTION ─────────────────────────────────
// reduceRegions computes the mean of both bands within each
// station buffer polygon in a single pass.

var zonalStats = LONGTERM_MEAN_IMG.select([BAND_NDWI_EVI, BAND_MNDWI_EVI])
  .reduceRegions({
    collection: stationBuffers,
    reducer: ee.Reducer.mean(),
    scale: SCALE
  });

print('Zonal statistics by station:', zonalStats);

// ── 2b. SANITY CHECK: MEAN ELEVATION PER BUFFER ──────────────
// Uses the SAME DEM used for pixel stratification in Section 4.3.1.
// If a station's buffer sits at systematically higher elevation than
// others, that alone would explain a low flood frequency there,
// independent of vegetation/EVI condition - this is the first thing
// to rule out before trusting the flood-frequency ranking.

// Correct DEM asset, confirmed by the thesis author - the same wetland-
// clipped ALS DEM used throughout Section 4.3.1 pixel sampling.
var DEM_RAW = ee.Image('projects/YOUR_GEE_PROJECT/assets/DGM_wetland_AT');

// CRITICAL FIX: the raster's NoData value is stored as literal 0, not as
// a proper masked/NoData flag. A lake surface at ~115 m a.s.l. can never
// legitimately read 0, so every 0-value pixel is a NoData artifact being
// silently averaged in as if it were real elevation - this is what was
// dragging every station's mean down toward 0 in previous runs.
// Masking out exact-0 pixels removes these artifacts before any statistics
// are computed.
var DEM = DEM_RAW.updateMask(DEM_RAW.gt(0));

print('Sanity check - min/max AFTER masking out 0-value NoData:',
  DEM.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: stationBuffers.geometry(),
    scale: 5,
    maxPixels: 1e9
  }));

var elevStats = DEM.reduceRegions({
  collection: stationBuffers,
  reducer: ee.Reducer.mean().combine({
    reducer2: ee.Reducer.count(),
    sharedInputs: true
  }),
  scale: 5 // match native ALS DEM resolution
});

print('Mean elevation and valid pixel count per buffer:', elevStats);


// ── 3. FORMAT AND EXPORT ──────────────────────────────────────

var zonalTable = zonalStats.map(function(f) {
  return ee.Feature(null, {
    station: f.get('station'),
    buffer_radius_m: f.get('buffer_radius_m'),
    NDWI_EVI_longterm_mean: f.get(BAND_NDWI_EVI),
    mNDWI_EVI_longterm_mean: f.get(BAND_MNDWI_EVI)
  });
});

print('Formatted zonal table:', zonalTable);

Export.table.toDrive({
  collection: zonalTable,
  description: 'table5_11_zonal_flood_frequency_by_station',
  folder: 'GEE_Exports',
  fileNamePrefix: 'table5_11_zonal_flood_frequency_by_station',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: elevStats,
  description: 'sanity_check_elevation_by_station',
  folder: 'GEE_Exports',
  fileNamePrefix: 'sanity_check_elevation_by_station',
  fileFormat: 'CSV'
});

print(' ');
print('READY: run BOTH export tasks in the Tasks tab.');
print('Compare mean elevation against flood frequency per station -');
print('if Neusiedl sits noticeably higher than Moerbisch, elevation');
print('alone may explain the low Neusiedl flood frequency.');

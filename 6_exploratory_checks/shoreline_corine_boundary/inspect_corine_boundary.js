// NOTE: study-area classes follow thesis Section 4.1 (CORINE 411 + 324).
// The deployed maps use the asset DGM_wetland_AT (DEM clipped to that boundary);
// this script rebuilds the mask directly from CORINE for exploration only.
// ============================================================
// GEE SCRIPT: Inspect CORINE boundary vs lake water edge
// Purpose: check whether the eastern shore red pixels are
// caused by the CORINE mask stopping before the actual
// water edge, or by high elevation at the water's edge
// ============================================================

// ── 1. LOAD CORINE 2018 ─────────────────────────────────────
var corine = ee.Image('COPERNICUS/CORINE/V20/100m/2018');

// Your wetland classes (from thesis Section 4.1)
var wetland = corine.eq(411).or(corine.eq(324));
var openWater = corine.eq(512);

// ── 2. STUDY AREA BOUNDARY ──────────────────────────────────
// Lake Neusiedl approximate bounding box
var lakeRegion = ee.Geometry.Rectangle([16.6, 47.7, 17.0, 47.95]);

Map.setCenter(16.80, 47.82, 12);
Map.setOptions('SATELLITE');

// ── 3. VISUALISE CORINE CLASSES ─────────────────────────────
Map.addLayer(
  wetland.selfMask().clip(lakeRegion),
  {palette: ['1d9e75'], opacity: 0.6},
  '1. CORINE wetland (Classes 411+324)'
);

Map.addLayer(
  openWater.selfMask().clip(lakeRegion),
  {palette: ['2a78d6'], opacity: 0.5},
  '2. CORINE open water (Class 512)'
);

// ── 4. WETLAND BOUNDARY EDGE ────────────────────────────────
// Find the edge of the wetland mask — pixels that are wetland
// but adjacent to open water or non-wetland
var wetlandEdge = wetland
  .focal_max(1, 'square', 'pixels')
  .subtract(wetland)
  .selfMask()
  .clip(lakeRegion);

Map.addLayer(
  wetlandEdge,
  {palette: ['ff0000'], opacity: 0.9},
  '3. Wetland boundary edge (red = edge pixels)'
);

// ── 5. DISTANCE FROM WETLAND EDGE TO OPEN WATER ─────────────
// This shows how far the CORINE wetland boundary sits from
// the classified open water — if there is a gap, that is the
// area excluded from the study
var distToWater = openWater
  .distance(ee.Kernel.euclidean(500, 'meters'))
  .clip(lakeRegion);

Map.addLayer(
  distToWater.updateMask(wetland),
  {min: 0, max: 500, palette: ['ff0000', 'ffff00', '00ff00']},
  '4. Distance from wetland pixels to open water (red=near, green=far)'
);

// ── 6. LOAD ALS ELEVATION ───────────────────────────────────
// Use FABDEM as proxy if Austrian ALS not available in GEE
// FABDEM removes buildings and vegetation — better than SRTM
var fabdem = ee.ImageCollection('projects/sat-io/open-datasets/FABDEM')
  .filterBounds(lakeRegion)
  .mosaic();

// Clip to wetland area only
var elevWetland = fabdem.clip(lakeRegion).updateMask(wetland);

Map.addLayer(
  elevWetland,
  {min: 114.5, max: 117.5,
   palette: ['2a78d6', '74c6e8', 'ffffcc', 'eb6834', 'c0392b']},
  '5. Elevation within CORINE wetland (blue=low, red=high)'
);

// ── 7. EAST vs WEST SHORE PROFILES ──────────────────────────
// Draw two transect lines to sample elevation:
// West shore transect (Mörbisch area) — perpendicular to shore
// East shore transect — perpendicular to shore on opposite side

var westTransect = ee.Geometry.LineString([
  [16.685, 47.795], [16.735, 47.795]  // west to east through Mörbisch reed belt
]);

var eastTransect = ee.Geometry.LineString([
  [16.85, 47.82], [16.90, 47.82]     // through eastern reed belt if it exists
]);

Map.addLayer(
  ee.FeatureCollection([
    ee.Feature(westTransect, {name: 'West transect'}),
    ee.Feature(eastTransect, {name: 'East transect'})
  ]),
  {color: 'yellow'},
  '6. Elevation transects (west and east)'
);

// Sample elevation along west transect
var westProfile = fabdem.sample({
  region: westTransect,
  scale: 30,
  numPixels: 100,
  geometries: true
});

// Sample elevation along east transect
var eastProfile = fabdem.sample({
  region: eastTransect,
  scale: 30,
  numPixels: 100,
  geometries: true
});

print('West transect elevation profile (sample):', westProfile.limit(10));
print('East transect elevation profile (sample):', eastProfile.limit(10));

// ── 8. EXPORT ELEVATION WITHIN WETLAND TO DRIVE ─────────────
Export.image.toDrive({
  image: elevWetland.rename('elevation'),
  description: 'wetland_elevation_FABDEM',
  folder: 'GEE_Exports',
  fileNamePrefix: 'wetland_elevation_FABDEM_30m',
  region: lakeRegion,
  scale: 30,
  crs: 'EPSG:4326',
  maxPixels: 1e9
});

// ── 9. EXPORT CORINE CLASS MAP TO DRIVE ─────────────────────
Export.image.toDrive({
  image: corine.clip(lakeRegion).rename('corine_class'),
  description: 'corine_2018_lake_neusiedl',
  folder: 'GEE_Exports',
  fileNamePrefix: 'corine_2018_lake_neusiedl',
  region: lakeRegion,
  scale: 100,
  crs: 'EPSG:4326',
  maxPixels: 1e9
});

print('✅ Run the export tasks in the Tasks tab');
print('Then use the R script to analyse east vs west elevation profiles');

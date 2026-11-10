// NOTE: study-area classes follow thesis Section 4.1 (CORINE 411 + 324).
// The deployed maps use the asset DGM_wetland_AT (DEM clipped to that boundary);
// this script rebuilds the mask directly from CORINE for exploration only.
// ============================================================
// GEE SCRIPT: Export CORINE boundary context map
// Shows why western shore has blue pixels at water edge
// and eastern shore does not
// ============================================================

// ── 1. LOAD CORINE 2018 ─────────────────────────────────────
var corine = ee.Image('COPERNICUS/CORINE/V20/100m/2018');

var lakeRegion = ee.Geometry.Rectangle([16.60, 47.70, 17.00, 47.98]);

// Define CORINE classes
var wetland   = corine.eq(411).or(corine.eq(324));  // 411 inland marshes + 324 transitional woodland-shrub (thesis Section 4.1)
var openWater = corine.eq(512);                       // water bodies
var other     = wetland.not().and(openWater.not());   // everything else

// ── 2. CREATE RGB CORINE CLASSIFICATION MAP ──────────────────
// Encode three classes as RGB image for export
// Wetland (411/324) = teal #1d9e75  → R:29  G:158 B:117
// Open water (512)  = blue #2a78d6  → R:42  G:120 B:214
// Other             = light grey    → R:220 G:220 B:220

var corineRGB = ee.Image.cat([
  // Red channel
  wetland.multiply(29)
    .add(openWater.multiply(42))
    .add(other.multiply(220)),
  // Green channel
  wetland.multiply(158)
    .add(openWater.multiply(120))
    .add(other.multiply(220)),
  // Blue channel
  wetland.multiply(117)
    .add(openWater.multiply(214))
    .add(other.multiply(220))
]).rename(['R', 'G', 'B'])
  .clip(lakeRegion);

// ── 3. VISUALISE ON MAP ─────────────────────────────────────
Map.setCenter(16.80, 47.82, 12);
Map.setOptions('SATELLITE');

Map.addLayer(
  openWater.selfMask().clip(lakeRegion),
  {palette: ['2a78d6'], opacity: 0.7},
  'Open water (Class 512)'
);

Map.addLayer(
  wetland.selfMask().clip(lakeRegion),
  {palette: ['1d9e75'], opacity: 0.7},
  'Wetland (Classes 411+324)'
);

// Show where wetland directly touches open water
// These are the pixels responsible for blue zones in flood map
var wetlandBuffer = wetland.focal_max(100, 'circle', 'meters');
var wetlandTouchingWater = wetlandBuffer.and(openWater);

Map.addLayer(
  wetlandTouchingWater.selfMask().clip(lakeRegion),
  {palette: ['ff0000'], opacity: 0.9},
  'Where wetland touches open water (red = direct contact zone)'
);

// ── 4. ANNOTATE KEY AREAS ───────────────────────────────────
// Western shore contact zone
var westContact = ee.Geometry.Point([16.695, 47.795]);
// Eastern shore — where gap exists
var eastGap = ee.Geometry.Point([16.88, 47.82]);

var annotations = ee.FeatureCollection([
  ee.Feature(westContact, {label: 'West: wetland directly touches lake'}),
  ee.Feature(eastGap,     {label: 'East: gap between lake and wetland'})
]);

Map.addLayer(
  annotations,
  {color: 'yellow'},
  'Key locations'
);

print('West shore contact zone:', westContact);
print('East shore gap area:', eastGap);

// ── 5. EXPORT CORINE CONTEXT MAP ───────────────────────────
// Export as GeoTIFF for use in thesis/presentation
Export.image.toDrive({
  image: corineRGB.toUint8(),
  description: 'CORINE_boundary_context_map',
  folder: 'GEE_Exports',
  fileNamePrefix: 'CORINE_boundary_wetland_water_contact',
  region: lakeRegion,
  scale: 100,
  crs: 'EPSG:32633',  // UTM zone 33N — same as Austrian ALS data
  maxPixels: 1e9
});

// Also export just the contact zone layer
Export.image.toDrive({
  image: wetlandTouchingWater.toUint8().clip(lakeRegion),
  description: 'Wetland_water_contact_zone',
  folder: 'GEE_Exports',
  fileNamePrefix: 'wetland_openwater_contact_zone',
  region: lakeRegion,
  scale: 100,
  crs: 'EPSG:32633',
  maxPixels: 1e9
});

// Export higher resolution version for presentation figure
Export.image.toDrive({
  image: corineRGB.toUint8(),
  description: 'CORINE_map_presentation_30m',
  folder: 'GEE_Exports',
  fileNamePrefix: 'CORINE_boundary_presentation_30m',
  region: lakeRegion,
  scale: 30,
  crs: 'EPSG:32633',
  maxPixels: 1e9
});

print('✅ Run three export tasks in the Tasks tab');
print('');
print('Map interpretation:');
print('TEAL = CORINE wetland (Classes 411 and 324) — your study area');
print('BLUE = CORINE open water (Class 512) — excluded from study');
print('RED  = Where wetland directly touches open water boundary');
print('');
print('The red zones on the western shore explain why blue pixels');
print('appear there in your flood frequency maps — those wetland');
print('pixels sit at the lake edge at the lowest elevation.');
print('');
print('The absence of red on the eastern shore explains why no');
print('blue pixels appear there — the wetland boundary never');
print('touches the open water class on that side.');

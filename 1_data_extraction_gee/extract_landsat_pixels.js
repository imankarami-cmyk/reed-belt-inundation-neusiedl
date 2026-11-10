// ============================================================
// FLATS+ — CORRECT PIXEL EXTRACTION FOR NEW SAMPLE
// Lake Neusiedl | 1984–2025
//
// PURPOSE
// - Extract harmonized Landsat reflectance bands
// - Keep pixel geometry exactly from uploaded lon/lat
// - Export nearest station for every pixel
// - Compute mNDWI + Pheno in the SAME band logic as thesis R code
//
// IMPORTANT
// 1. Upload sample_pixels_with_nearest_station-2.csv to GEE Assets
// 2. Replace PIXELS_ASSET with your uploaded asset path
// 3. Run all 4 export tasks
// ============================================================


// ── 0. CONFIG ────────────────────────────────────────────────
var PIXELS_ASSET = 'projects/YOUR_GEE_PROJECT/assets/sample_pixels_with_nearest_station'; // change if needed
var START_DATE   = '1984-01-01';
var END_DATE     = '2025-12-31';
var SCALE        = 30;
var TILE_SCALE   = 8;
var CLOUD_THR    = 70;


// ── 1. LOAD PIXEL TABLE AND REBUILD GEOMETRY FROM lon/lat ───
// This is safer than trusting CSV auto-geometry import.
// It guarantees sampleRegions uses your uploaded coordinates.

var ptsTable = ee.FeatureCollection(PIXELS_ASSET);

var pts = ptsTable.map(function(f) {
  var lon = ee.Number.parse(f.get('lon'));
  var lat = ee.Number.parse(f.get('lat'));

  return ee.Feature(
    ee.Geometry.Point([lon, lat]),
    f.toDictionary([
      'pixel_id',
      'elev_bin',
      'bin_lower_m',
      'bin_upper_m',
      'elevation',
      'lon',
      'lat',
      'nearest_station',
      'station_lon',
      'station_lat',
      'distance_to_station_km'
    ])
  );
});

print('Total pixels loaded:', pts.size());
print('Sample pixel:', pts.first());

Map.setCenter(16.75, 47.80, 11);
Map.setOptions('SATELLITE');
Map.addLayer(pts.style({color: 'ffff00', pointSize: 4}), {}, 'New sample pixels');

var stationCounts = pts.reduceColumns({
  reducer: ee.Reducer.frequencyHistogram(),
  selectors: ['nearest_station']
});
print('Pixels per nearest station:', stationCounts);


// ── 2. CLOUD MASKS ───────────────────────────────────────────
// Same core logic as your old scripts, but keep snow masked too.

var maskL57 = function(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0)   // dilated cloud
    .and(qa.bitwiseAnd(1 << 3).eq(0))      // cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))      // cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));     // snow
  return image.updateMask(mask);
};

var maskL89 = function(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0)   // dilated cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))      // cirrus
    .and(qa.bitwiseAnd(1 << 3).eq(0))      // cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))      // cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0));     // snow
  return image.updateMask(mask);
};


// ── 3. HARMONIZE BANDS TO THESIS NAMING ──────────────────────
// Final harmonized band meaning used everywhere below:
// B1 = Blue
// B2 = Green
// B3 = Red
// B4 = NIR
// B5 = SWIR1
// B7 = SWIR2
//
// This matches the old GEE export structure and allows consistent
// downstream index calculation. [file:2913]
//
// IMPORTANT FOR YOUR THESIS R SCRIPT:
// In your training R script, mNDWI and Pheno were defined from
// harmonized B2/B3/B5 logic:
// mNDWI = (B2 - B5) / (B2 + B5)
// Pheno = (B2 - B3) / (B2 + B3)
// even though the comments there describe Green/Blue/NIR naming. [file:2845]
//
// To stay 100% consistent with the dataset format used before,
// we export the harmonized raw bands AND compute the indices in GEE.
// Then in R, use the exported mNDWI and Pheno directly.

var applyScaleOffset = function(image) {
  return image.multiply(2.75e-05).add(-0.2);
};

var prepL57 = function(image) {
  var scaled = applyScaleOffset(
    maskL57(image).select(
      ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
      ['B1',    'B2',    'B3',    'B4',    'B5',    'B7']
    )
  );

  return scaled.copyProperties(image, ['system:time_start']);
};

var prepL89 = function(image) {
  var scaled = applyScaleOffset(
    maskL89(image).select(
      ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
      ['B1',    'B2',    'B3',    'B4',    'B5',    'B7']
    )
  );

  return scaled.copyProperties(image, ['system:time_start']);
};


// ── 4. ADD INDICES EXACTLY FOR YOUR MODEL TEST ───────────────
// This follows the variable structure used in your model workflow:
// mNDWI, Pheno, plus seasonal sin/cos. [file:2845]

var addPredictors = function(image) {
  var B1 = image.select('B1'); // Blue
  var B2 = image.select('B2'); // Green
  var B3 = image.select('B3'); // Red
  var B4 = image.select('B4'); // NIR
  var B5 = image.select('B5'); // SWIR1

  var mNDWI = B2.subtract(B5)
    .divide(B2.add(B5).add(1e-10))
    .rename('mNDWI');

  var Pheno = B2.subtract(B3)
    .divide(B2.add(B3).add(1e-10))
    .rename('Pheno');

  var NDWI = B2.subtract(B4)
    .divide(B2.add(B4).add(1e-10))
    .rename('NDWI');

  var EVI = B4.subtract(B3).multiply(2.5)
    .divide(B4.add(B3.multiply(6)).subtract(B1.multiply(7.5)).add(1).add(1e-10))
    .rename('EVI');

  var date = ee.Date(image.get('system:time_start'));
  var doy = ee.Number(date.getRelative('day', 'year')).add(1);
  var pi2 = 2 * Math.PI / 365.25;

  var sin_doy = ee.Image.constant(doy.multiply(pi2).sin()).rename('sin_doy');
  var cos_doy = ee.Image.constant(doy.multiply(pi2).cos()).rename('cos_doy');

  return image.addBands([mNDWI, Pheno, NDWI, EVI, sin_doy, cos_doy])
              .set('doy', doy)
              .set('year', date.get('year'))
              .set('month', date.get('month'));
};


// ── 5. BUILD COLLECTIONS ─────────────────────────────────────
var aoi = pts.geometry().bounds().buffer(500);

var L5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(aoi)
  .filterDate('1984-01-01', '2012-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_THR))
  .map(prepL57)
  .map(addPredictors)
  .map(function(img) { return img.set('sensor', 'L5'); });

var L7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(aoi)
  .filterDate('1999-01-01', '2013-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_THR))
  .map(prepL57)
  .map(addPredictors)
  .map(function(img) { return img.set('sensor', 'L7'); });

var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(aoi)
  .filterDate('2013-01-01', '2022-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_THR))
  .map(prepL89)
  .map(addPredictors)
  .map(function(img) { return img.set('sensor', 'L8'); });

var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
  .filterBounds(aoi)
  .filterDate('2022-01-01', END_DATE)
  .filter(ee.Filter.lt('CLOUD_COVER', CLOUD_THR))
  .map(prepL89)
  .map(addPredictors)
  .map(function(img) { return img.set('sensor', 'L9'); });

print('L5 scenes:', L5.size());
print('L7 scenes:', L7.size());
print('L8 scenes:', L8.size());
print('L9 scenes:', L9.size());


// ── 6. SAMPLE FUNCTION ───────────────────────────────────────
// filter(notNull(['B1'])) is kept so only valid sampled pixels remain,
// same logic as your older working extraction code. [file:2913]

var sampleSensor = function(collection, sensorLabel) {
  return collection.map(function(image) {
    var date = ee.Date(image.get('system:time_start'));

    return image.sampleRegions({
      collection: pts,
      properties: [
        'pixel_id',
        'elev_bin',
        'bin_lower_m',
        'bin_upper_m',
        'elevation',
        'lon',
        'lat',
        'nearest_station',
        'station_lon',
        'station_lat',
        'distance_to_station_km'
      ],
      scale: SCALE,
      tileScale: TILE_SCALE,
      geometries: false
    })
    .filter(ee.Filter.notNull(['B1']))
    .map(function(f) {
      return f.set({
        'date': date.format('YYYY-MM-dd'),
        'year': date.get('year'),
        'month': date.get('month'),
        'doy': image.get('doy'),
        'sensor': sensorLabel
      });
    });
  }).flatten();
};


// ── 7. RUN SAMPLING ──────────────────────────────────────────
var l5_data = sampleSensor(L5, 'L5');
var l7_data = sampleSensor(L7, 'L7');
var l8_data = sampleSensor(L8, 'L8');
var l9_data = sampleSensor(L9, 'L9');

print('Sampled L5 rows:', l5_data.size());
print('Sampled L7 rows:', l7_data.size());
print('Sampled L8 rows:', l8_data.size());
print('Sampled L9 rows:', l9_data.size());


// ── 8. QUICK VALUE CHECK ─────────────────────────────────────
// Just to confirm indices exist and are numeric.

print('Example sampled row:', l8_data.first());


// ── 9. EXPORT ────────────────────────────────────────────────
var selectors = [
  'pixel_id',
  'nearest_station',
  'station_lon',
  'station_lat',
  'distance_to_station_km',
  'lon',
  'lat',
  'elevation',
  'elev_bin',
  'bin_lower_m',
  'bin_upper_m',
  'date',
  'year',
  'month',
  'doy',
  'sensor',
  'B1',
  'B2',
  'B3',
  'B4',
  'B5',
  'B7',
  'mNDWI',
  'Pheno',
  'NDWI',
  'EVI',
  'sin_doy',
  'cos_doy'
];

Export.table.toDrive({
  collection: l5_data,
  description: 'FLATS_newpixels_L5_1984_2012',
  folder: 'GEE_Exports',
  fileNamePrefix: 'FLATS_newpixels_L5_1984_2012',
  fileFormat: 'CSV',
  selectors: selectors
});

Export.table.toDrive({
  collection: l7_data,
  description: 'FLATS_newpixels_L7_1999_2013',
  folder: 'GEE_Exports',
  fileNamePrefix: 'FLATS_newpixels_L7_1999_2013',
  fileFormat: 'CSV',
  selectors: selectors
});

Export.table.toDrive({
  collection: l8_data,
  description: 'FLATS_newpixels_L8_2013_2022',
  folder: 'GEE_Exports',
  fileNamePrefix: 'FLATS_newpixels_L8_2013_2022',
  fileFormat: 'CSV',
  selectors: selectors
});

Export.table.toDrive({
  collection: l9_data,
  description: 'FLATS_newpixels_L9_2022_2025',
  folder: 'GEE_Exports',
  fileNamePrefix: 'FLATS_newpixels_L9_2022_2025',
  fileFormat: 'CSV',
  selectors: selectors
});

print(' ');
print('✅ READY: run all 4 tasks in the Tasks tab');
print('Output includes mNDWI, Pheno, NDWI, EVI, sin_doy, cos_doy');
print('Output also keeps nearest_station for every pixel');
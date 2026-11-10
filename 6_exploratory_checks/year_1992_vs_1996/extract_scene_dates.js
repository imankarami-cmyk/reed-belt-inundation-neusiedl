// ============================================================
// GEE SCRIPT: Extract all cloud-free scene dates for 1992
// and 1996 across the full CORINE wetland boundary
// Simplified version — extracts dates directly without
// pixel counting inside map() to avoid type errors
// ============================================================

// ── 1. STUDY AREA ────────────────────────────────────────────
var aoi = ee.Geometry.Rectangle([16.65, 47.72, 16.98, 47.97]);

Map.setCenter(16.80, 47.82, 11);
Map.setOptions('SATELLITE');

// ── 2. CLOUD MASKING ─────────────────────────────────────────
var maskL57 = function(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 1).eq(0)
    .and(qa.bitwiseAnd(1 << 3).eq(0))
    .and(qa.bitwiseAnd(1 << 4).eq(0))
    .and(qa.bitwiseAnd(1 << 5).eq(0));
  return image.updateMask(mask);
};

// ── 3. BUILD COLLECTIONS ─────────────────────────────────────
var L5_1992 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(aoi)
  .filterDate('1992-01-01', '1992-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', 70))
  .map(maskL57);

var L5_1996 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(aoi)
  .filterDate('1996-01-01', '1996-12-31')
  .filter(ee.Filter.lt('CLOUD_COVER', 70))
  .map(maskL57);

print('1992 scenes (cloud < 70%):', L5_1992.size());
print('1996 scenes (cloud < 70%):', L5_1996.size());

// ── 4. EXTRACT SCENE DATES AS FEATURE COLLECTION ─────────────
// Simple approach: convert each image's metadata to a Feature
// No pixel counting inside map() — just extract the date info

var getSceneDates = function(collection, yearLabel) {
  return ee.FeatureCollection(
    collection.map(function(image) {
      var date = ee.Date(image.get('system:time_start'));
      return ee.Feature(null, {
        date:        date.format('YYYY-MM-dd'),
        year:        date.get('year'),
        month:       date.get('month'),
        doy:         date.getRelative('day', 'year').add(1),
        cloud_cover: image.get('CLOUD_COVER'),
        scene_id:    image.get('system:index'),
        sensor:      'L5'
      });
    })
  );
};

var sceneDates_1992 = getSceneDates(L5_1992, 1992);
var sceneDates_1996 = getSceneDates(L5_1996, 1996);

print('Scene dates 1992 (first 10):', sceneDates_1992.limit(10));
print('Scene dates 1996 (first 10):', sceneDates_1996.limit(10));

// ── 5. MONTHLY DISTRIBUTION ───────────────────────────────────
// Count scenes per month for each year
var monthCount = function(fc) {
  return ee.List.sequence(1, 12).map(function(m) {
    var count = fc.filter(ee.Filter.eq('month', m)).size();
    return ee.Feature(null, {month: m, n_scenes: count});
  });
};

print('1992 monthly scene distribution:',
      ee.FeatureCollection(monthCount(sceneDates_1992)));
print('1996 monthly scene distribution:',
      ee.FeatureCollection(monthCount(sceneDates_1996)));

// ── 6. EXPORT TO DRIVE ───────────────────────────────────────
Export.table.toDrive({
  collection: sceneDates_1992,
  description: 'scene_dates_1992_full_wetland',
  folder: 'GEE_Exports',
  fileNamePrefix: 'scene_dates_1992_full_wetland',
  fileFormat: 'CSV',
  selectors: ['date', 'year', 'month', 'doy', 'cloud_cover', 'scene_id', 'sensor']
});

Export.table.toDrive({
  collection: sceneDates_1996,
  description: 'scene_dates_1996_full_wetland',
  folder: 'GEE_Exports',
  fileNamePrefix: 'scene_dates_1996_full_wetland',
  fileFormat: 'CSV',
  selectors: ['date', 'year', 'month', 'doy', 'cloud_cover', 'scene_id', 'sensor']
});

print('');
print('✅ Run both export tasks in the Tasks tab');
print('Files will appear in GEE_Exports folder in Google Drive');

// ============================================================
// ALOS PALSAR-1 (L-band) vs SENTINEL-1 (C-band)
// Step 1: extract both at the SAME 600 pixels into ONE table.
// Exploratory comparison (not part of the thesis).
// Common columns: co_dB (HH or VV), cross_dB (HV or VH).
// ============================================================

var PIXELS_ASSET = 'projects/YOUR_GEE_PROJECT/assets/sample_pixels_with_nearest_station';
var STATION      = null;          // e.g. 'Illmitz' to export one area at a time; null = all
var S1_START     = '2015-01-01';
var S1_END       = '2025-12-31';
var S1_REL_ORBIT = null;          // e.g. 117 to keep ONE orbit track (fairer, smaller)
var ALOS_YEARS   = [2007, 2008, 2009, 2010];

// ── 1. THE 600 THESIS PIXELS ─────────────────────────────────
var ptsTable = ee.FeatureCollection(PIXELS_ASSET);
var pts = ptsTable.map(function(f) {
  var lon = ee.Number.parse(f.get('lon'));
  var lat = ee.Number.parse(f.get('lat'));
  return ee.Feature(
    ee.Geometry.Point([lon, lat]),
    f.toDictionary(['pixel_id', 'elev_bin', 'elevation', 'lon', 'lat',
                    'nearest_station', 'distance_to_station_km'])
  );
});
if (STATION !== null) pts = pts.filter(ee.Filter.eq('nearest_station', STATION));
print('Pixels:', pts.size());
var aoi = pts.geometry().bounds().buffer(500);

// ── 2. SENTINEL-1: every acquisition, speckle-reduced ────────
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(aoi)
  .filterDate(S1_START, S1_END)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'));
if (S1_REL_ORBIT !== null) {
  s1 = s1.filter(ee.Filter.eq('relativeOrbitNumber_start', S1_REL_ORBIT));
}
print('Sentinel-1 scenes:', s1.size());

var s1Fc = s1.map(function(img) {
  // average in linear units over about 30 m, back to dB
  var lin = ee.Image(10).pow(img.select(['VV', 'VH']).divide(10))
              .focal_mean({radius: 1, kernelType: 'square', units: 'pixels'});
  var db = lin.log10().multiply(10).rename(['co_dB', 'cross_dB'])
              .addBands(img.select('angle'));
  var date = ee.Date(img.get('system:time_start'));
  return db.reduceRegions({collection: pts, reducer: ee.Reducer.first(),
                           scale: 10, tileScale: 4})
    .filter(ee.Filter.notNull(['co_dB', 'cross_dB']))
    .map(function(f) {
      return f.set({
        sensor:    'S1',
        date:      date.format('YYYY-MM-dd'),
        year:      date.get('year'),
        pass:      img.get('orbitProperties_pass'),
        rel_orbit: img.get('relativeOrbitNumber_start')
      });
    });
}).flatten();

// ── 3. ALOS PALSAR-1 yearly mosaics ──────────────────────────
// Print this once to learn whether a 'date' band exists (see R script).
var alosCol = ee.ImageCollection('JAXA/ALOS/PALSAR/YEARLY/SAR');
var m0 = ee.Image(alosCol.filterDate('2007-01-01', '2008-01-01').first());
print('ALOS bands:', m0.bandNames());
print('ALOS values at first pixel:', m0.reduceRegion({
  reducer: ee.Reducer.first(), geometry: pts.first().geometry(), scale: 25}));

var alosFc = ee.FeatureCollection(ALOS_YEARS.map(function(y) {
  var m  = ee.Image(alosCol.filterDate(y + '-01-01', (y + 1) + '-01-01').first());
  var hh = m.select('HH').pow(2).log10().multiply(10).subtract(83.0).rename('co_dB');
  var hv = m.select('HV').pow(2).log10().multiply(10).subtract(83.0).rename('cross_dB');
  var db = hh.addBands(hv).focal_median({radius: 1, units: 'pixels'});
  var dateBand = ee.Image(ee.Algorithms.If(
    m.bandNames().contains('date'),
    m.select('date').rename('date_band'),
    ee.Image.constant(-9999).rename('date_band')));
  return db.addBands(dateBand).reduceRegions({
      collection: pts, reducer: ee.Reducer.first(), scale: 25})
    .filter(ee.Filter.notNull(['co_dB', 'cross_dB']))
    .map(function(f) { return f.set({sensor: 'ALOS1', year: y}); });
})).flatten();

// ── 4. ONE TABLE, ONE EXPORT ─────────────────────────────────
var all = s1Fc.merge(alosFc);

Export.table.toDrive({
  collection: all,
  description: 'RADAR_ALOS1_vs_S1_600pixels',
  folder: 'GEE_Exports',
  fileNamePrefix: 'RADAR_ALOS1_vs_S1_600pixels',
  fileFormat: 'CSV',
  selectors: ['sensor', 'pixel_id', 'nearest_station', 'elevation', 'elev_bin',
              'year', 'date', 'pass', 'rel_orbit',
              'co_dB', 'cross_dB', 'angle', 'date_band']
});

print('Run the export in the Tasks tab. If it times out, set STATION to one');
print('area, or set S1_REL_ORBIT to one track, and export in several parts.');

/**
 * HAT-Coast: Selected Google Earth Engine functions for coastal-water reconstruction
 *
 * This file provides selected core functions from the HAT-Coast framework,
 * which supports biennial global mainland shoreline mapping from 2000 to 2020.
 * The included functions implement locally adaptive optical and SAR water
 * extraction, mapping-period assignment, hierarchical monthly gap filling,
 * annual Water Inundation Frequency (WIF) calculation, and persistent-water classification.
 *
 * Monthly water observations are reconstructed hierarchically in the order
 * JRC -> GLAD -> locally adaptive optical/SAR water mask -> GLC_FCS30D.
 * Valid JRC observations are preserved, and each supplementary source is
 * applied only to pixels that remain unresolved after the preceding step.
 * Annual WIF is then calculated from the reconstructed monthly observations,
 * and persistent water is identified using the validation-derived threshold
 * of 71%.
 *
 * Input class conventions: JRC, prepared GLAD, and prepared GLC_FCS30D: 0 = unresolved, 1 = land, 2 = water.
 *   Locally adaptive optical and SAR masks: 0 = non-water, 1 = water.
 *   Final persistent-water classification:
 *     0 = land or intermittently inundated area;
 *     1 = persistent water.
 */


// -----------------------------------------------------------------------------
// 1. Locally adaptive optical water extraction (2000-2004)
// -----------------------------------------------------------------------------
/**
 * Extract water from a Landsat composite within one coastal grid cell.
 * Extreme AWEI values outside the 2nd-98th percentile range are excluded only
 * from threshold estimation. The derived threshold is applied to the complete
 * AWEI image before terrain screening.
 */
function extractOpticalWater(startDate, endDate, region, outputName) {
  var composite = buildLandsatComposite(startDate, endDate, region);
  var awei = calculateAWEInsh(composite);

  var percentiles = awei.reduceRegion({
    reducer: ee.Reducer.percentile([2, 98]),
    geometry: region,
    scale: 30,
    maxPixels: 1e13,
    bestEffort: true
  });

  var p2 = ee.Number(percentiles.get('AWEI_p2'));
  var p98 = ee.Number(percentiles.get('AWEI_p98'));
  var thresholdSample = awei.updateMask(awei.gte(p2).and(awei.lte(p98)));

  var histogram = thresholdSample.reduceRegion({
    reducer: ee.Reducer.histogram(255)
      .combine('mean', null, true)
      .combine('variance', null, true),
    geometry: region,
    scale: 30,
    maxPixels: 1e13,
    bestEffort: true
  });

  var threshold = ee.Number(otsuThreshold(histogram.get('AWEI_histogram')));
  var initialWater = awei.gt(threshold);

  var elevation = SRTM.select('elevation');
  var slope = ee.Terrain.slope(elevation);
  var terrainMask = elevation.lt(5).and(slope.lt(5)).clip(region);
  var screenedWater = initialWater.and(terrainMask);

  return {
    water: screenedWater.rename(outputName).toByte(),
    threshold: threshold,
    awei: awei,
    initialWater: initialWater.rename('initial_' + outputName).toByte()
  };
}


// -----------------------------------------------------------------------------
// 2. Locally adaptive SAR water extraction (2006-2020)
// -----------------------------------------------------------------------------

/**
 * Extract water from the PALSAR/PALSAR-2 HH mosaic for a source year within
 * one coastal grid cell. Low HH backscatter is classified as water, followed
 * by the same 5 m and 5 degree terrain screening used in the supplied workflow.
 */
function extractSarWater(sourceYear, region, outputName) {
  var startDate = ee.Date.fromYMD(sourceYear, 1, 1);
  var endDate = startDate.advance(1, 'year');

  var sarImage = ee.Image(
    ee.ImageCollection('JAXA/ALOS/PALSAR/YEARLY/SAR')
      .filterDate(startDate, endDate)
      .first()
  );

  var Db = convertPalsarToDb(sarImage);

  var percentiles = Db.reduceRegion({
      reducer: ee.Reducer.percentile([2, 98]),
      geometry: region,
      scale: 30,
      maxPixels: 1e13,
      bestEffort: true
    });

  var p2 = ee.Number(percentiles.get('HH_dB_p2'));
  var p98 = ee.Number(percentiles.get('HH_dB_p98'));
    
  var dbForHistogram = Db.updateMask(Db.gte(p2).and(Db.lte(p98)));
	
  var histogram = dbForHistogram.reduceRegion({ 
      reducer: ee.Reducer.histogram(255).combine('mean', null, true).combine('variance', null, true), 
      geometry: region, 
      scale: 30,
      maxPixels: 1e13,
      bestEffort: true 
    });

  var histJSON = ee.String.encodeJSON(histogram.get('HH_dB_histogram'));

  var threshold = ee.Number(otsuThreshold(histogram.get('HH_dB_histogram')));
  var initialWater = Db.lt(threshold);

  var elevation = SRTM.select('elevation');
  var slope = ee.Terrain.slope(elevation);
  var terrainMask = elevation.lt(5).and(slope.lt(5)).clip(region);
  var screenedWater = initialWater.and(terrainMask);

  return {
    water: screenedWater.rename(outputName).toByte(),
    threshold: threshold,
    backscatterDb: Db,
    initialWater: initialWater.rename('initial_' + outputName).toByte()
  };
}


/**
 * Generate the annual adaptive water mask for one biennial mapping label.
 * The mapping year must be a client-side even integer from 2000 to 2020.
 */
function getAnnualAdaptiveWater(mappingYear, region) {
  var opticalPeriods = {
  2000: ['2000-01-01', '2002-01-01'],
  2002: ['2002-01-01', '2004-01-01'],
  2004: ['2004-01-01', '2006-01-01']
};

  var regularSarYears = [2008, 2010, 2016, 2018, 2020];
  var water;

  if (opticalPeriods[mappingYear]) {
    var period = opticalPeriods[mappingYear];
    water = extractOpticalWater(
      period[0],
      period[1],
      region,
      String(mappingYear)
    ).water;
  } else if (mappingYear === 2006) {
    // The 2007 PALSAR mosaic represents the 2006-2007 mapping period.
    water = extractSarWater(2007, region, '2006').water;
  } else if (mappingYear === 2012) {
    // No SAR mosaic was available for 2012-2013. The mean > 0.5 criterion
    // retains only pixels classified as water in both 2010 and 2015.
    var water2010 = extractSarWater(2010, region, '2010').water;
    var water2015 = extractSarWater(2015, region, '2015').water;
    water = water2010.add(water2015).divide(2).gt(0.5);
  } else if (mappingYear === 2014) {
    // The 2015 PALSAR-2 mosaic represents the 2014-2015 mapping period.
    water = extractSarWater(2015, region, '2014').water;
  } else if (regularSarYears.indexOf(mappingYear) !== -1) {
    water = extractSarWater(
      mappingYear,
      region,
      String(mappingYear)
    ).water;
  } else {
    throw new Error('mappingYear must be an even integer from 2000 to 2020.');
  }

  return water.rename('annual_water').toByte().clip(region);
}


// -----------------------------------------------------------------------------
// 3. Hierarchical monthly gap filling
// -----------------------------------------------------------------------------
/**
 * Fill one monthly JRC image in the order JRC -> GLAD -> adaptive water mask
 * -> GLC_FCS30D. Only unresolved pixels (value 0) are modified at each step.
 * The GLAD image must have been prepared as 0 = unresolved, 1 = land,
 * and 2 = water before it is supplied to this function.
 */
function fillMonthlyWater(jrcMonth,preparedGladMonth,annualAdaptiveWater,preparedGlcImage,region) {
  var jrcWater = jrcMonth.select('water');
  var gladWater = preparedGladMonth.select('water');
  var glcWater = preparedGlcImage.select('water');

  // Step 1: temporally corresponding GLAD observation.
  var afterGlad = jrcWater.where(jrcWater.eq(0),gladWater.unmask(jrcWater));

  // Step 2: locally adaptive optical/SAR mask. Only water classifications
  // supported by the adaptive mask are inserted at this stage.
  var adaptiveWaterMask = afterGlad.eq(0).and(annualAdaptiveWater.eq(1).unmask(0));
  var afterAdaptiveWater = afterGlad.where(adaptiveWaterMask, 2);

  // Step 3: annual land-cover information assigns the residual gaps.
  var glcLandMask = afterAdaptiveWater.eq(0).and(glcWater.eq(1));
  var glcWaterMask = afterAdaptiveWater.eq(0).and(glcWater.eq(2));

  return afterAdaptiveWater.where(glcLandMask, 1).where(glcWaterMask, 2)
    .rename('water').clip(region)
    .copyProperties(jrcMonth, jrcMonth.propertyNames());
}


// -----------------------------------------------------------------------------
// 4. Annual water inundation frequency and persistent-water classification
// -----------------------------------------------------------------------------

/** Calculate annual WIF from a collection of reconstructed monthly images. */
function calculateAnnualWif(filledMonthlyCollection, region) {
  var monthlyIndicators = filledMonthlyCollection.map(function(image) {
    var validObservation = image.select('water').gt(0).rename('observation');
    var waterObservation = image.select('water').eq(2).rename('water');
    return validObservation.addBands(waterObservation)
      .copyProperties(image, image.propertyNames());
  });

  var totalObservations = monthlyIndicators
    .select('observation')
    .sum()
    .toFloat();
  var totalWater = monthlyIndicators.select('water').sum().toFloat();

  return totalWater
    .divide(totalObservations)
    .multiply(100)
    .updateMask(totalObservations.gt(0))
    .rename('WIF')
    .clip(region);
}


/** Apply the validation-derived 71% WIF threshold. */
function classifyPersistentWater(wifImage, mappingYear, region) {
  return wifImage.gte(71).rename(String(mappingYear)).toByte().clip(region).set('year', mappingYear)
    .set(
      'system:time_start',
      ee.Date.fromYMD(mappingYear, 1, 1).millis()
    );
}
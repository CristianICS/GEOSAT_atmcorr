import os
import math
import json
import subprocess
from datetime import datetime
from xml.dom import minidom
from osgeo_utils import gdal_calc # type: ignore
from osgeo import gdal # type: ignore
gdal.UseExceptions()

import ee # type: ignore
from Py6S import * # type: ignore

def utcdate_to_datetime(UTCtime: str):
    """Convert string with UTC date from dim files into python datetime"""
    date = datetime.strptime(UTCtime, "%Y-%m-%dT%H:%M:%S")
    return date

class DIM:
    """Functions to retrieve metadata values from .dim files."""
    def __init__(self, dim_path):
        self._path = dim_path
        # Open DIM file
        self.file = minidom.parse(dim_path)

    def _get_crs(self):
        """Function to retrieve CRS from the DIM metadata.
        It is not used.

        from pyproj import CRS
        # Fetch the PROJ4 string
        proj4 = self.file.getElementsByTagName('Projection_OGCWKT')

        # Check the data
        try:
            if proj4[0].localName == 'Projection_OGCWKT':
                ogc_wkt_string = proj4[0].firstChild.data
                # Get EPSG code from OGC WKT string
                crs = CRS.from_wkt(ogc_wkt_string)
                return crs.source_crs.to_epsg()
            else:
                raise ValueError(' '.join([
                'The CRS code is not inside <Projection_OGCWKT tag>. Insert',
                  'it manually.'
                ]))
        except:
            raise ValueError(' '.join([
              'The CRS code is not inside <Projection_OGCWKT tag>. Insert it',
              'manually.'
            ]))
        """
  
    def get_bandKeys(self):
        # Init the band list
        bands = []
        spec_info = self.file.getElementsByTagName('Spectral_Band_Info')
        for i in spec_info:
            band_key = i.getElementsByTagName('BAND_INDEX')[0].firstChild.data
            bands.append(int(band_key))

        return bands

    def get_spectral(self, keys: list, band: int):
        """
        Get values from XML Spectral_Band_Info (.dim metadata file)
        Return a dict with XML key: val by band

        :param keys: List with the keys to retrieve.
        :param band: Band index which contains the keys to retrieve.
        """
        params = {} # init dict
        # Get band spectral info
        spec_info = self.file.getElementsByTagName('Spectral_Band_Info')
        # Iterate over spectral info childs
        # Store only the band_key == band
        for c in spec_info:
            band_key = c.getElementsByTagName('BAND_INDEX')[0].firstChild.data
            if int(band_key) == band:
                for key in keys:
                    val = c.getElementsByTagName(key)[0].firstChild.data
                    params[key] = val

        return params

    def get_azimuth(self):
        """Retrieve SATELLITE_AZIMUTH key from .dim file"""
        # Get quality parameter keys
        qp = self.file.getElementsByTagName('Quality_Parameter')
        # Iterate over its items
        for child in qp:
            qp_desc = child.getElementsByTagName('QUALITY_PARAMETER_DESC')
            # Get the key of the first child
            key = qp_desc[0].firstChild.data
            if key == 'SATELLITE_AZIMUTH':
                values = child.getElementsByTagName('QUALITY_PARAMETER_VALUE')
                azimuth = values[0].firstChild.data

                return float(azimuth)
        # When azimuth value has not been found, return None
        return None
    
    def get_firstOrder_keys(self, keys: list):
        """
        Get values from XML by keys

        It's only for unique keys (appears 1 in the XML tree).

        Note: Integers values are returned as strings.
        """
        params = {}
        for key in keys:
            val = self.file.getElementsByTagName(key)[0].firstChild.data
            params[key] = val

        return params
    
    def get_bbox(self):
        """
        Retrieve BBOX coords.

        The first 4 FRAME_LON and FRAME_LAT are the BBOX
        NW, NE, SE, SW coordinates in WGS84 (GEE accepted EPSG)

        Coordinates must be inserted into GEE with this format:
        ['NWLONG','NWLAT','SELONG','SELAT']
        """
        coords_long = self.file.getElementsByTagName('FRAME_LON')
        coords_lat = self.file.getElementsByTagName('FRAME_LAT')
        
        # Retrieve NWLONG and NWLAT
        nwlong = float(coords_long[0].firstChild.data)
        nwlat = float(coords_lat[0].firstChild.data)
        # Retrieve SELONG and SELAT
        selong = float(coords_long[2].firstChild.data)
        selat = float(coords_lat[2].firstChild.data)

        return([nwlong, nwlat, selong, selat])

class Atmospheric():
  """
  atmospheric.py, Sam Murphy (2016-10-26)
  
  Atmospheric water vapour, ozone and AOT from GEE
  
  Usage
  H2O = Atmospheric.water(geom,date)
  O3 = Atmospheric.ozone(geom,date)
  AOT = Atmospheric.aerosol(geom,date)
  
  """
  def round_date(date,xhour):
    """
    rounds a date of to the closest 'x' hours
    """
    y = date.get('year')
    m = date.get('month')
    d = date.get('day')
    H = date.get('hour')
    HH = H.divide(xhour).round().multiply(xhour)
    return date.fromYMD(y,m,d).advance(HH,'hour')
  
  def round_month(date):
    """
    round date to closest month
    """
    # start of THIS month
    m1 = date.fromYMD(date.get('year'),date.get('month'),ee.Number(1))
    
    # start of NEXT month
    m2 = m1.advance(1,'month')
      
    # difference from date
    d1 = ee.Number(date.difference(m1,'day')).abs()
    d2 = ee.Number(date.difference(m2,'day')).abs()
    
    # return closest start of month
    return ee.Date(ee.Algorithms.If(d2.gt(d1),m1,m2))
  
  
  
  def water(geom,date):
    """
    Water vapour column above target at time of image aquisition.
    
    (Kalnay et al., 1996, The NCEP/NCAR 40-Year Reanalysis Project. Bull. 
    Amer. Meteor. Soc., 77, 437-471)
    """
    
    # Point geometry required
    centroid = geom.centroid()
    
    # H2O datetime is in 6 hour intervals
    H2O_date = Atmospheric.round_date(date,6)
    
    # filtered water collection
    water_ic = ee.ImageCollection('NCEP_RE/surface_wv').filterDate(H2O_date, H2O_date.advance(1,'month'))
    
    # water image
    water_img = ee.Image(water_ic.first())
    
    # water_vapour at target
    water = water_img.reduceRegion(reducer=ee.Reducer.mean(), geometry=centroid).get('pr_wtr')
                                        
    # convert to Py6S units (Google = kg/m^2, Py6S = g/cm^2)
    water_Py6S_units = ee.Number(water).divide(10)                                   
    
    return water_Py6S_units
  
  
  
  def ozone(geom,date):
    """
    returns ozone measurement from merged TOMS/OMI dataset
    
    OR
    
    uses our fill value (which is mean value for that latlon and day-of-year)
  
    """
    
    # Point geometry required
    centroid = geom.centroid()
       
    def ozone_measurement(centroid,O3_date):
      
      # filtered ozone collection
      ozone_ic = ee.ImageCollection('TOMS/MERGED').filterDate(O3_date, O3_date.advance(1,'month'))
      
      # ozone image
      ozone_img = ee.Image(ozone_ic.first())
      
      # ozone value IF TOMS/OMI image exists ELSE use fill value
      ozone = ee.Algorithms.If(ozone_img,\
      ozone_img.reduceRegion(reducer=ee.Reducer.mean(), geometry=centroid).get('ozone'),\
      ozone_fill(centroid,O3_date))
      
      return ozone
      
    def ozone_fill(centroid,O3_date):
      """
      Gets our ozone fill value (i.e. mean value for that doy and latlon)
      
      you can see it
      1) compared to LEDAPS: https://code.earthengine.google.com/8e62a5a66e4920e701813e43c0ecb83e
      2) as a video: https://www.youtube.com/watch?v=rgqwvMRVguI&feature=youtu.be
      
      """
      
      # ozone fills (i.e. one band per doy)
      ozone_fills = ee.ImageCollection('users/samsammurphy/public/ozone_fill').toList(366)
      
      # day of year index
      jan01 = ee.Date.fromYMD(O3_date.get('year'),1,1)
      doy_index = date.difference(jan01,'day').toInt()# (NB. index is one less than doy, so no need to +1)
      
      # day of year image
      fill_image = ee.Image(ozone_fills.get(doy_index))
      
      # return scalar fill value
      return fill_image.reduceRegion(reducer=ee.Reducer.mean(), geometry=centroid).get('ozone')
     
    # O3 datetime in 24 hour intervals
    O3_date = Atmospheric.round_date(date,24)
    
    # TOMS temporal gap
    TOMS_gap = ee.DateRange('1994-11-01','1996-08-01')  
    
    # avoid TOMS gap entirely
    ozone = ee.Algorithms.If(TOMS_gap.contains(O3_date),ozone_fill(centroid,O3_date),ozone_measurement(centroid,O3_date))
    
    # fix other data gaps (e.g. spatial, missing images, etc..)
    ozone = ee.Algorithms.If(ozone,ozone,ozone_fill(centroid,O3_date))
    
    #convert to Py6S units 
    ozone_Py6S_units = ee.Number(ozone).divide(1000)# (i.e. Dobson units are milli-atm-cm )                             
    
    return ozone_Py6S_units
 

  def aerosol(geom,date):
    """
    Aerosol Optical Thickness.
    
    try:
      MODIS Aerosol Product (monthly)
    except:
      fill value
    """
    
    def aerosol_fill(date):
      """
      MODIS AOT fill value for this month (i.e. no data gaps)
      """
      return ee.Image('users/samsammurphy/public/AOT_stack')\
               .select([ee.String('AOT_').cat(date.format('M'))])\
               .rename(['AOT_550'])
               
               
    def aerosol_this_month(date):
      """
      MODIS AOT original data product for this month (i.e. some data gaps)
      """
      # image for this month
      img =  ee.Image(\
                      ee.ImageCollection('MODIS/061/MOD08_M3')\
                        .filterDate(Atmospheric.round_month(date))\
                        .first()\
                     )
      
      # fill missing month (?)
      img = ee.Algorithms.If(img,\
                               # all good
                               img\
                               .select(['Aerosol_Optical_Depth_Land_Mean_Mean_550'])\
                               .divide(1000)\
                               .rename(['AOT_550']),\
                              # missing month
                                aerosol_fill(date))

      return img


    def get_AOT(AOT_band,geom):
      """
      AOT scalar value for target
      """  
      return ee.Image(AOT_band).reduceRegion(reducer=ee.Reducer.mean(),\
                                 geometry=geom.centroid())\
                                .get('AOT_550')
                                

    after_modis_start = date.difference(ee.Date('2000-03-01'),'month').gt(0)
    
    AOT_band = ee.Algorithms.If(after_modis_start, aerosol_this_month(date), aerosol_fill(date))
    
    AOT = get_AOT(AOT_band,geom)
    
    AOT = ee.Algorithms.If(AOT,AOT,get_AOT(aerosol_fill(date),geom))
    # i.e. check reduce region worked (else force fill value)
    
    return AOT
  
class Correction:
    """Atmospheric Correction Functions"""

    def get_metadata(band_key: int, nd_min: int, dim: DIM):
        """Return band metadata to perform DOS and COST formulas."""
        # Get ESUN
        e = dim.get_spectral('ESUN', band_key)['ESUN']
        # Retrieve the rest of the keys
        keys = [
           'SUN_AZIMUTH',
           'EARTH_SUN_DISTANCE',
           'SUN_ELEVATION'
        ]
        metadata = dim.get_firstOrder_keys(keys)
        metadata['ESUN'] = e
        # Compute solar zenith angle
        zen_sol_ang = round(90 - float(metadata['SUN_ELEVATION']), 2)
        metadata['zen_sol_ang'] = zen_sol_ang
        # Dark object computation
        lmin = Correction.ARC(band_key, '', dim, nd_min)
        # Reflectance (TOA) from dark object (L1%)
        d = metadata['EARTH_SUN_DISTANCE']
        lone = (0.01 * float(e) * abs(math.cos(zen_sol_ang))) / (float(d)**2 * math.pi)
        # Convert Reflectance to Radiance
        lone_radiance = ((float(d)**2) * math.pi) / (lone * abs(math.cos(zen_sol_ang)))
        # Path radiance (Lhaze)
        lhaze = lmin - lone_radiance

        metadata['lhaze'] = lhaze
        metadata['lonea_radiance'] = lone_radiance
        metadata['lone'] = lone
        return metadata
        

    def ARC(band_key: int, band_letter: str, dim: DIM, dn_min = None):
        """
        Compute the formula to perform Absolute Radiometric Correction (ARC)
        Formula: L = DN * Gain + Bias
        
        Returns a string with the formula to compute it inside gdal_calc.

        The required parameters are inside gesosat image .dim metadata file.

        :param band_key: Same key as the one represented the band in dim files.
        :param band_letter: Letter which will link band in gdal_calc comand.
        :param dim: DIM object.
        :param dn_min: If a DN min value is passed, the formula calcs radiance
        for dn_min value. It's used in DOS and COST functions.
        """
        # Return band metadata
        keys = ['PHYSICAL_GAIN', 'PHYSICAL_BIAS']
        band_metadata = dim.get_spectral(keys, band_key)
        gain = band_metadata['PHYSICAL_GAIN']
        bias = band_metadata['PHYSICAL_BIAS']

        # Return radiance value from ND min
        if dn_min != None:
            radiancia_min = abs(float(dn_min) * float(gain) + float(bias))
            return radiancia_min
        else:
            # Retrun ARC formula
            # Note: DN values corresponds with image pixel
            formula = f"abs( {band_letter} * {gain} + {bias} )"

            # Write log messages
            logs = [
                f'<ARC band={band_key}>',
                f'[formula] {formula}'
            ]

            return (formula, logs)

    def DOS(band_key: int, band_letter: str, dim: DIM, dn_min: int):
        """
        Compute Dark Object Substraction method (DOS)
        
        Formula: REF(BOA) = pi * (L - Lhaze) * d**2 / Elambda * cos(theta)
        Lhaze = Lmin - L1%
        
        The formula unknows are inside DIM metadata file.

        IMPORTANT: L1% is obtained in TOA reflectance. It's mandatory to
        transform it into radiance again. The L1% formula is:

        L1% (TOA) = (0.01 * Radiance * cos(cenSolAng)) / ((d^2) * pi)

        Then get the inverse ecuation isolating the radiance:

        Radiance = ((d^2) * pi) / (L1% * cos(cenSolAng))

        With the second formula we obtain the L1% in Radiance units.

        :param band_key: Same key as the one represented the band in dim files.
        :param band_letter: Letter which represents band in gdal_calc comand.
        :param dim: Object of class DIM.
        :param dn_min: Dark object value (transformed to radiance with ARC)
        """
        # Store required band metadata
        band_metadata = Correction.get_metadata(band_key, dn_min, dim)
        # Zenith solar angle
        zensol = band_metadata['zen_sol_ang']
        # ESUN
        e = band_metadata['ESUN']
        d = band_metadata['EARTH_SUN_DISTANCE']
        # Lhaze
        lhaze = band_metadata['lhaze']
        
        # Write DOS formula
        formula = (f"(pi * (({band_letter} - {lhaze}) * ({d}**2))) /" +
        f"({e} * abs(cos({zensol})))")
        # Save logs
        logs = [
            f'<DOS band={band_key}>',
            f'[formula] f{formula}',
            f"[NDmin units=radiance] {band_metadata['lmin']}",
            f"[L1% units=radiance] {band_metadata['lone_radiance']}"
        ]

        return (formula, logs)
    
    def COST(band_key: int, band_letter: str, dim: DIM, dn_min: str):
        """
        Compute COST
        REF(BOA) = pi * (L - Lhaze) * d**2 / Elambda * cos(theta) * TAUz
        Lhaze = Lmin - L1%

        The TAUz is computed by COST method proposed by Chavez (1996)
        """
        # Store required band metadata
        band_metadata = Correction.get_metadata(band_key, dn_min, dim)
        # Zenith solar angle
        zensol = band_metadata['zen_sol_ang']
        # ESUN
        e = band_metadata['ESUN']
        d = band_metadata['EARTH_SUN_DISTANCE']
        # Lhaze
        lhaze = band_metadata['lhaze']

        # Compute TAUz
        TZ = math.radians(zensol)
        TAUz = 1 - (TZ**2/2) + (TZ**4/(4*3*2)) - (TZ**6/(6*5*4*3*2))

        # Write formula
        formula = (f"(pi * (({band_letter} - {lhaze}) * ({d}**2))) /" +
            f"({e} * abs(cos({zensol})) * {TAUz})")
        # Save logs
        logs = [
            f'<COST band={band_key}>',
            f'[formula] f{formula}',
            f"[TAUz] {TAUz}"
        ]

        return (formula, logs)
    
    def sixS(band_key: int, band_letter: str, dim: DIM, gee_id: str, is_pan=False):
        """
        Compute 6S atmospheric correction model.

        Formula: REF(BOA) = pi * (L - Lp) / t * ( Edir + Edif )

        Lp = Path radiance
        t = transmissivity (absorption transmissivity (TAUv?) * scattering transmissivity (TAUz?))
        Edir = Direct Solar Irradiance (Eo?)
        Edif = Diffuse solar irradiance (Edown?)
        
        Terms follow by ? are the closest approx to 6S formula terms (from Sam
        Murphy repo) in Chavez 1996 radiance to BOA reflectance function.

        :band_key: Same as the key represented the band in metadata.
        :band_letter: Letter which represent band in gdal_calc comand.
        :dim: Object of class DIM.
        :gee_id: Cloud project id from which GEE account is linked.
        :is_pan: When panchromatic image, set different wavelength limits.
        """
        # Open GEE python API
        ee.Initialize(project=gee_id)

        # Retrieve image BBOX
        coords = dim.get_bbox() # EPSG:4326
        # Test footprint
        # var coords = [
        # -0.7785354999436517,
        # 42.3782895196768905,
        # -0.5899609097725235,
        # 42.2822183047948812
        # ]
        # var footprint = ee.Geometry.Rectangle(coords)
        # Map.addLayer(footprint)

        # Create geom poligon with image footprint
        footprint = ee.Geometry.Rectangle(coords)
        # Get center point
        geom = footprint.centroid()

        # 6S object
        # Core class of Py6S. It allows to define the input parameters,
        # to run the radiative transfer code and to access the outputs which
        # are required to convert radiance to surface reflectance.
        s = SixS() # type: ignore

        # Geometric conditions with Geometry class
        # - Month
        # - Day
        # - Sun zenith and azimuth angles
        # - Sensor/View zenith and azimuth angles
        s.geometry = Geometry.User() # type: ignore
        keys = [
            'START_TIME',
            'VIEWING_ANGLE',
            'SUN_ELEVATION',
            'SUN_AZIMUTH'
        ]
        img_metadata = dim.get_firstOrder_keys(keys)

        # Date
        strDate = img_metadata['START_TIME']
        captureDate = utcdate_to_datetime(strDate) # UTC format
        s.geometry.month = captureDate.month # Month (used for Earth-Sun distance)
        s.geometry.day = captureDate.day # Day (used for Earth-Sun distance)
        
        # Sensor zenith angle (NADIR)
        s.geometry.view_z = 0
        # Sensor azimuth angle
        s.geometry.view_a =  dim.get_azimuth()
        # Solar zenith angle
        solar_z = round(90 - float(img_metadata['SUN_ELEVATION']), 2) 
        s.geometry.solar_z = solar_z
        # Solar azimuth angle
        s.geometry.solar_a = float(img_metadata['SUN_AZIMUTH'])

        # Atmospheric profile and Aerosol Profile
        # Import manually the atmospheric constitutients
        # AOT at 550nm
        # Atmospheric constituents
        captureDate_str = f"{captureDate.year}-{captureDate.month}-{captureDate.day}"
        date = ee.Date(captureDate_str)
        h2o = Atmospheric.water(geom,date).getInfo()
        o3 = Atmospheric.ozone(geom,date).getInfo()
        aot = Atmospheric.aerosol(geom,date).getInfo()

        # Atmospheric constituents
        s.atmos_profile = AtmosProfile.UserWaterAndOzone(h2o, o3) # type: ignore
        # s.aero_profile = AeroProfile.Continental
        s.aero_profile = AeroProfile.Continental # type: ignore
        s.aot550 = aot

        # Altitudes. Allows the specification of target and sensor altitudes.
        # Altitude - Shuttle Radar Topography mission
        SRTM = ee.Image('CGIAR/SRTM90_V4') 
        alt = SRTM.reduceRegion(
            reducer = ee.Reducer.mean(),
            geometry = geom.centroid()
        ).get('elevation').getInfo()
        km = alt/1000 # i.e. Py6S uses units of kilometers
        s.altitudes.set_target_custom_altitude(km)
        #s.altitudes.set_sensor_custom_altitude(sensor_altitude)
        # https://py6s.readthedocs.io/en/stable/params.html#altitudes
        # s.altitudes.set_sensor_satellite_level()
        s.altitudes.set_sensor_custom_altitude(99) # 100km is the maximun altitude

        # Wavelenght conditions
        # Wavelenght spectral response is computed with Wavelength()
        # passing the start and end band wavelength (micrometers)
        # - Wavelength function: 
        # https://github.com/robintw/Py6S/blob/master/Py6S/Params/wavelength.py
        # - WV3 Spectral Response:
        # https://earth.esa.int/eogateway/missions/geosat-2#instruments-section
        # Soruce: GEOSAT-2-Imagery-User-Guide.pdf p. 4
        if is_pan:
            lowerWav = 0.560
            upperWav = 0.900
        else:
            wavelengts_geosat = {
                # band-key: [lowerBandEdge, upperBandEdge]
                1: [0.770, 0.892],
                2: [0.640, 0.697],
                3: [0.532, 0.599],
                4: [0.466, 0.525],
            }
            lowerWav = wavelengts_geosat[band_key][0]
            upperWav = wavelengts_geosat[band_key][1]

        s.wavelength = Wavelength(lowerWav, upperWav) # type: ignore

        # run 6S for this waveband
        s.run()

        # extract 6S outputs
        Edir = s.outputs.direct_solar_irradiance # direct solar irradiance
        Edif = s.outputs.diffuse_solar_irradiance # diffuse solar irradiance
        Lp = s.outputs.atmospheric_intrinsic_radiance # path radiance
        # absorption transmissivity
        absorb = s.outputs.trans['global_gas'].upward 
        # scattering transmissivity
        scatter = s.outputs.trans['total_scattering'].upward 
        tau2 = absorb*scatter # total transmissivity

        # Compute surface reflectance formula
        formula = f'(pi*({band_letter} - {Lp}))/({tau2}*({Edir}+{Edif}))'
        
        # Save logs
        logs = [
            f'<6S band={band_key}>',
            f'[BBOX NWLONG,NWLAT,SELONG,SELAT]{coords}',
            "[6S inputs]", f"H2O: {h2o}", f"O3: {o3}", f"AOT: {aot}",
            f"Altitude: {km}", f"Solar zenith angle: {solar_z}",
            "[6S outputs]", f"EDIR: {Edir}", f"EDIF: {Edif}",
            f"Lp: {Lp}", f"Absorb: {absorb}", f"Scatter: {scatter}",
            f"Transmissivity (total): {tau2}", f'[formula] {formula}']

        return (formula, logs)

class Geosat:
    """Handle GEOSAT image correction process."""
    def __init__(self, folder_path):
        # Handle GEOSAT image folder
        self.test_path(folder_path)
        self.folder = folder_path
        # Init the log file
        self.create_logfile()

        # Get images names inside the folder
        self.img_names = self.get_original_images()
        # Get dim paths
        self.dim_paths = [self.get_dim(img) for img in self.img_names]

    def test_path(self, path: str):
        """Check if a path exists."""
        if os.path.exists(path):
            return True
        else:
            raise ValueError(f'{path} path cannot be resolved.')

    def get_original_images(self):
        """Write original GEOSAT image names inside the root folder.
        
        Folder name:
        DE2_PM4_L1C_000000_20210629T102606_20210629T102608_DE2_38083_838B
        
        Images names:
        DE2_[MS4|PAN]_L1C_000000_20210629T102606_20210629T102608_DE2_38083_838B
        """
        mul_img = os.path.basename(self.folder).replace('PM4', 'MS4') + '.tif'
        pan_img = os.path.basename(self.folder).replace('PM4', 'PAN') + '.tif'
        # Add extension

        # Test if images exist
        self.test_path(os.path.join(self.folder, mul_img))
        self.test_path(os.path.join(self.folder, pan_img))

        return [mul_img, pan_img]

    def get_dim(self, img_path: str):
        """Retrieve image dim file with the metadata info."""
        dim_name = img_path.replace('tif', 'dim')
        dim_path = os.path.join(self.folder, dim_name)
        if self.test_path(dim_path):
            return dim_path

    def create_logfile(self):
        """Create the log file path.
        The atmospheric correction's metadata will be written in it.
        """
        self.log_path = os.path.join(self.folder, 'atm_correction.log')
    
    def write_logs(self, lines: str | list):
        """Open the log file and write the log info in it."""
        log = open(self.log_path, 'a')
        # Insert current datetime
        c = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        # Check if there are multiple lines to write
        if type(lines) != type([]):
            lines = [lines]
        for line in lines:
            # Insert log line
            log.write(c + ' - ' + line + '\n')

    def get_minDn(self, img_path):
        """Get image bands min values (GDAL is required)"""
        # Get the stats by band in JSON format
        gdalinfo = f'gdalinfo -json -hist "{img_path}"'
        # Handle JSON response
        infojson = json.loads(subprocess.check_output(gdalinfo, shell=True))

        # Store min digital numbers
        min_dn_values = []
        # Min digital numbers are inside ['bands'] attribute
        for band in infojson['bands']:
            dn_value = band['histogram']['min'] # Float type
            # Handle errors
            if dn_value > 0:
                min_dn_values.append(str(dn_value))
            else:
                min_dn_values.append('1')
        return min_dn_values

    def find_ARC(self, img_path):
        """Find the ARC image"""
        img_arc_path = img_path.replace('.tif', ('_ARC.tif'))
        if os.path.exists(img_arc_path):
            return img_arc_path
        else:
            raise ValueError("ARC image is not exist and it's mandatory to perform the desired correction.")
    
    def is_duplicated(self, path):
        """Check if an image is duplicated."""
        if os.path.exists(path):
            err_txt = f'{os.path.basename(path)} is already exist. Stop the image creation.'
            raise ValueError(err_txt)

    def get_atm_formulas(self, fun: Correction, band_keys: list, **args):
        """
        Create atmcorr formulas by band to insert them inside gdal_calc
        
        If DOS/COST/6S is performed, an image corrected with ARC formula must
        exist.

        :fun: One of the correction functins inside Correction() class.
        :band_keys: Image bands keys
        :dim: DIM object
        :**args: Custom arguments to include inside fun parameter.
        """
        # Assign default letter to band keys
        abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        letters = {}
        for i in range(0, len(band_keys)):
            letters[band_keys[i]] = abc[i]

        # Calculate gdal_calc parameters by band
        formulas = []
        selected_letters = {}

        # Create ATM CORR formulas by band
        for key in band_keys:

            formula, logs = fun(key, letters[key], **args)
            self.write_logs(logs)
            # Add final band formula to the list of formulas
            formulas.append(formula)
            # Add the selected letter
            selected_letters[key] = letters[key]

        return(formulas, selected_letters)
    
    def run_gdal_calc(self, formulas, selected_letters, img_path, out_path, gdaltype = None):
        """Run gdal_calc.Calc function with custom params.

        The function must contain each parameters as bands have the image to
        correct.

        https://github.com/OSGeo/gdal/blob/master/swig/python/gdal-utils/osgeo_utils/gdal_calc.py

        :formulas: List with atmcorr formulas by band.
        :selected_letters: Dict with the band key and its letter which is its
        name inside the formula.
        :img_path: Image to correct path.
        :gdaltype: Data type of the output image.
        """
        # Init the object to store the Calc function arguments
        args = {}
        
        # Include the image band parameters (letter and band number)
        for key, letter in selected_letters.items():
            # Write param to identify the letter with the image.
            # The image is the same in all bands, but it must be specified
            # in each letter.
            args[letter] = img_path
            # Write the band letter (to identify the band with the letter)
            args[f'{letter}_band']= key
        
        # Write creation options
        args['creation_options'] = ['COMPRESS=DEFLATE', 'PREDICTOR=2']

        if type(gdaltype) != type(None):
            args['type'] = gdaltype

        # Construct the final gdal_calc command
        try:
            gdal_calc.Calc(calc=formulas, outfile=out_path, **args)
        except Exception as ex:
            ex_msg = f'The image {os.path.basename(out_path)} is already exist. Stop the image creation.'
            self.write_logs(ex_msg)
            print(ex_msg)

    def atm_corr(self, atm_key: str, **args):
        """Perform the atmospheric correction.
        
        Handle each of the atmcorr types and run gdal_calc.Calc to do the task.

        :atm_key: Correction type.
        :**args: Custom arguments. Current only the gee_id param.
        """
        # Correct each image
        for img, dim in zip(self.img_names, self.dim_paths):

            # Construct the image path
            img_path = os.path.join(self.folder, img)
            # Open DIM object
            dim = DIM(dim)
            # Get band keys
            band_keys = dim.get_bandKeys()

            self.write_logs(f'Init the correction of {img}')

            if atm_key == 'ARC':
                # Check if output image exists
                img_output = img_path.replace('.tif', (f'_{atm_key}.tif'))
                self.is_duplicated(img_output)
                # Create band formulas
                formulas, lettrs = self.get_atm_formulas(
                    Correction.ARC, band_keys, dim=dim)
                # Run the gdal_calc .py script from GDAL
                self.run_gdal_calc(formulas, lettrs, img_path, img_output)

            elif atm_key == 'DOS':
                # Check for ARC image
                img_arc_path = self.find_ARC(img_path)
                # Check if output image exists
                img_output = img_arc_path.replace('.tif', (f'_{atm_key}.tif'))
                self.is_duplicated(img_output)
                # Create band formulas
                formulas, lettrs = self.get_atm_formulas(
                    Correction.DOS, band_keys, dim=dim)
                # Run the gdal_calc .py script from GDAL
                self.run_gdal_calc(formulas, lettrs, img_arc_path, img_output)

            elif atm_key == 'COST':
                # Check for ARC image
                img_arc_path = self.find_ARC(img_path)
                # Check if output image exists
                img_output = img_arc_path.replace('.tif', (f'_{atm_key}.tif'))
                self.is_duplicated(img_output)
                # Create band formulas
                formulas, lettrs = self.get_atm_formulas(
                    Correction.COST, band_keys, dim=dim)
                # Run the gdal_calc .py script from GDAL
                self.run_gdal_calc(formulas, lettrs, img_arc_path, img_output)

            elif atm_key == '6S':
                # Check for ARC image
                img_arc_path = self.find_ARC(img_path)
                # Check if output image exists
                img_output = img_arc_path.replace('.tif', (f'_{atm_key}.tif'))
                self.is_duplicated(img_output)
                # Check if image is a panchromatic one
                is_pan = img.find('PAN') > 0
                # Create band formulas
                formulas, lettrs = self.get_atm_formulas(
                    Correction.sixS, band_keys, dim=dim, is_pan=is_pan, **args)
                # Run the gdal_calc .py script from GDAL
                self.run_gdal_calc(formulas, lettrs, img_arc_path, img_output, 'Float32')

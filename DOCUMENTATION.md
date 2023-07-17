<div align='center'>
<img src="https://github.com/kylejgillett/sounderpy/assets/100786530/2e9477c9-e36a-4163-accb-fe46780058dd" width="200">
</div>

# SOUNDERPY DOCUMENTATION
  > LATEST VERSION: v0.0.2  |  RELEASED: July 9, 2023  |  COPYRIGHT Kyle J Gillett, 2023
<br>

This script is used to access vertical profile data for calculations or plotting of a vertical profile (sounding). 

## SounderPy is currently capable of accessing and processing data from:

- ECMWF CDS ERA5 reanalysis [1940-present] *note: you must set up an account through the CDS to unlock ERA5 data. (see: https://cds.climate.copernicus.eu/api-how-to)
- UNIDATA THREDDS TDS RAP reanalysis [2005-present]
- UNIDATA THREDDS TDS RUC reanalysis [2005-2020]
- The University of Wyoming RAOB archive [1973-present, depending on station]
- Iowa State University's RAOB archive [1945-present, depending on station]

## SounderPy is used for:
- Accessing and loading raw vertical profile data from the sources listed above
- Parsing these raw data into a clean & simple-to-use format for calculations and plotting
- SounderPy offers a built-in quick MetPy plotting function, but sounderPy itself is meant as a source for accessing data -- not for plotting (that may be a future package ;) )
-------
## How to use SounderPy:
1. #### sounderpy is available through pip!
   ```
   pip install sounderpy
   ```
2. Make sure your environment has the required dependencies:
   -  metpy
   -  xarray
   -  urllib
   -  numpy
   -  siphon
   -  cartopy
   -  datetime
   -  cdsapi
3. Finally, its fun to run: `import sounderpy as spy`
------
## Tools included in SounderPy:
- #### `get_docs()`
    - function: print a quick look at the documentation
    - ##### EXAMPLES
          spy.get_docs()
-------- 
- #### `get_model_data(method, latlons, year, month, day, hour, domain='point')`
  - function used to access model profile data from the ERA5, RAP or RUC
  - ##### PARAMETERS
          name      dtype  required  example         explanation 
          method:    str,    yes,    'rap'           model data source: 'rap-ruc', 'era5' 
          latlon:    list,   yes,    [42.3, -97.3]   ordered list of lat, lon points 
          year:      str,    yes,    '2014'          4 digit year (UTC)
          month:     str,    yes,    '06'            2 digit month (UTC)
          day:       str,    yes,    '16'            2 digit day (UTC)
          hour:      str,    yes,    '20'            2 digit hour (UTC)
          domain:    str,    no,     'point'         default='point', single lat/lon 'point' or 'map' for lat/lon box *under dev*
  - ##### RETURNS
          raw_data: for ERA5, an xarray.Dataset of ERA5 reanalysis data or
          raw_data: for RAP/RUC, a netCDF4._netCDF4.Dataset of RAP/RUC reanalysis data
  - ##### EXAMPLES
          raw_data = spy.get_model_data('rap', 'point', [42.3, -97.3], '2014', '06', '16', '18')
          raw_data = spy.get_model_data('era5', 'point', [42.3, -97.3], '2014', '06', '16', '18')
    *note: method 'rap-ruc' will try for RAP data first, if the data is not available, then  it will try for RUC data*
--------
- #### `get_obs_data(station, year, month, day, hour)`
  - function used to access *and* parse RAOB profile data
  - this function will search both UW & ISU for desired station and date
  - ##### PARAMETERS
          name      dtype  required  example         explanation 
          station:   str,    yes,    'OAX'           3 digit RAOB station identifier
          year:      str,    yes,    '2014'          4 digit year (UTC)
          month:     str,    yes,    '06'            2 digit month (UTC)
          day:       str,    yes,    '16'            2 digit day (UTC)
          hour:      str,    yes,    '18'            2 digit hour (UTC)
  - ##### RETURNS
          clean_data: a dict of cleaned profile data ready to use for plotting
  - ##### EXAMPLES
          clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')
    *note: this function will first attempt to get data from the UW archive, if it is not available from UW then it will search the ISU IEM archive*
--------
- #### `parse_data(raw_data)`
  - function used to parse and clean up raw model data *NOTE: only needed for *model data*
  - ##### PARAMETERS
          name       dtype    required   example      explanation 
          raw_data: dataset,    yes,    raw_data      raw_data is produced by get_model_data()
  - ##### RETURNS
          clean_data: a dict of cleaned profile data ready to use for plotting
  - ##### EXAMPLES
          clean_data = spy.parse_data(raw_data)
--------
- #### `metpy_sounding(clean_data)`
  - this function inputs cleaned profile data through a slightly modified MetPy sounding plot script. As mentioned above, sounderPy is a tool for accessing and processing sounding data and as such this is only a very simple built-in plotting tool
  - ##### PARAMETERS
          name         dtype    required   example      explanation 
          clean_data:  dict,      yes,    clean_data    a dict of cleaned profile data acquired from the above functions. 
  - ##### RETURNS
          img: MetPy Sounding Plot
  - ##### EXAMPLES
          spy.metpy_sounding(clean_sounding)

------
## SounderPy Helper Functions 
- Several helper functions are included within SounderPy that a user may find useful when working with sounding data
  - #### `metar_latlon(metar_site)`
    - ##### PARAMETERS
            name        dtype    required   example   explanation 
            metar_site  str      yes        'KMBS'    4 digit METAR site id
    - ##### RETURNS
            latlon:     list     lat lon pair for requested METAR site in a list
    - ##### EXAMPLES
            latlon = spy.metar_latlon('KMBS')
        *note: you can use this lat/lon pair list when calling get_model_data()!*
      
   - #### `raob_latlon(raob_site)`
      - ##### PARAMETERS
              name        dtype    required   example   explanation 
              raob_site  str      yes        'DTX'     3-4 digit RAOB site id
      - ##### RETURNS
              latlon:     list     lat lon pair for requested RAOB site in a list
      - ##### EXAMPLES
              latlon = spy.raob_latlon('OUN')
          *note: you can use this lat/lon pair list when calling get_model_data()!*
        
    - #### `interp_data(variable, heights, step=100)`
      - ##### PARAMETERS
              name        dtype    required   explanation 
              variable:   array    yes        the data to be interpolated. Must be same length as height array
              heights:    array    yes        heights corresponding to the vertical profile used to interpolate
              step:       int      no         resolution of interpolation. Default is 100 (recommended value is 100)
          
      - ##### RETURNS
              varinterp:  array,  array of interpolated data 
      - ##### EXAMPLES
              interp_temp = spy.interp_data(temp, heights, step=100)
  
  
    - #### `find_nearest(array, value)`
      - ##### PARAMETERS
              name         dtype                        required    explanation 
              array        array of ints or floats      yes         the array of data search through to 'find the nearest'
              value        int or float (same as array) yes         the value used to compare against the array of data
          
      - ##### RETURNS
              nearest_idx:  int,  index of the data array that corresponds with the nearest value to the given value.  
      - ##### EXAMPLES
              z_equals_500m = spy.interp_data(z, 500)
  
------
## The future of SounderPy
- Adding more data sources, such as NWP forecasts (HRRR, RAP, NAM, GFS)
- RAP real-time analysis

-------









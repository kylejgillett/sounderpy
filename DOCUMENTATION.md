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
1. #### Download `sounderpy.py` and add it to your file directory
2. `import sounderpy as spy`
------
## Tools included in SounderPy:
- #### `get_docs()`
    - function: print a quick look at the documentation
    - ##### EXAMPLES
          spy.get_docs()
-------- 
- #### `get_model_data(method, domain, latlons, year, month, day, hour)`
  - function used to access model profile data from the ERA5, RAP or RUC
  - ##### PARAMETERS
          name      dtype  required  example         explanation 
          method:    str,    yes,    'rap'           model data source: 'rap', 'era5' 
          domain:    str,    no,     'point'         single 'point' or 'map' *under dev.
          latlon:    list,   yes,    [42.3, -97.3]   ordered list of lat, lon points 
          year:      str,    yes,    '2014'          4 digit year (UTC)
          month:     str,    yes,    '06'            2 digit month (UTC)
          day:       str,    yes,    '16'            2 digit day (UTC)
          hour:      str,    yes,    '20'            2 digit hour (UTC)
  - ##### RETURNS
          raw_data: for ERA5, an xarray.Dataset of ERA5 reanalysis data or
          raw_data: for RAP/RUC, a netCDF4._netCDF4.Dataset of RAP/RUC reanalysis data
  - ##### EXAMPLES
          raw_data = spy.get_obs_data('rap', 'point', [42.3, -97.3], '2014', '06', '16', '18')
          raw_data = spy.get_obs_data('era5', 'point', [42.3, -97.3], '2014', '06', '16', '18')    
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
## The future of SounderPy
- Adding more data sources, such as NWP forecasts (HRRR, RAP, NAM, GFS)
- RAP real-time analysis










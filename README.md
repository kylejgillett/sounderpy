<div align="center">
<img src="https://github.com/kylejgillett/sounderpy/assets/100786530/2e9477c9-e36a-4163-accb-fe46780058dd" width="250">
</div>

# sounderpy

------
[![PyPI Package](https://img.shields.io/pypi/v/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)]([https://pypi.python.org/pypi/ansicolortags/](https://github.com/kylejgillett/sounderpy/blob/main/LICENSE.txt))
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![GitHub commits](https://badgen.net/github/commits/kylejgillett/sounderpy)](https://GitHub.com/kylejgillett/sounderpy/commit/)

[![Maintainer](https://img.shields.io/badge/maintainer-kylejgillett-blue)](https://github.com/kylejgillett)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)



## WELCOME! thank you for visiting SounderPy!
This script is used to access vertical profile data for calculations or plotting of a vertical profile (sounding). 

## SounderPy News
### A new release of SounderPy will be out soon! (Early August)
#### Here is a look at what's coming to SounderPy v1.1.0:
+ A few minor bug fixes
+ Functionality supporting access to the IGRAv2 radiosonde archive
+ Functionality for accessing most-recent RAP analysis data
+ Built in functionality for saving a .png of the 'simple metpy sounding' function
+ New helper functions!
   + A function that retrieves the latitude and longitude of a US buoy/CMAN site for marine profile applications
   + A function that retrieves the latitude and longitude of an IGRAv2 radiosonde site
   + A function that saves parsed sounding data in a .csv file
   + A function that saves parsed sounding data in a CM1 input_sounding file
#### Some of the planned additions to SounderPy:
+ Give users the option to chose between pressure levels and model levels when using the ERA5
+ Add functionality supporting access to the IEM Bufkit archive
+ Increasing SounderPy's plotting capabilities

-----

## Why SounderPy?
+ Sometimes data is tough to find, and often times is even tougher to get it in the format you like. SounderPy gets you this data!
+ The code needed for loading and parsing vertical data (especially from models) can be large and messy. SounderPy keeps it hidden away in a PyPi package -- just import and call sounderPy functions to keep your code clean!

## SounderPy is used for:
- Accessing and loading raw vertical profile data from the sources listed below
- Parsing these raw data into a clean & simple-to-use format for calculations and plotting
- SounderPy offers a built-in quick MetPy plotting function, but sounderPy itself is meant as a source for accessing data -- not for plotting (that may be a future package ;) )
-------

## SounderPy is currently capable of accessing and processing data from:
- ECMWF CDS ERA5 reanalysis [1940-present] *note: you must set up an account through the CDS to unlock ERA5 data. (see: https://cds.climate.copernicus.eu/api-how-to)
- UNIDATA THREDDS TDS RAP reanalysis [2005-present]
- UNIDATA THREDDS TDS RUC reanalysis [2005-2020]
- The University of Wyoming RAOB archive [1973-present, depending on station]
- Iowa State University's RAOB archive [1945-present, depending on station]
-------
## How to use SounderPy:
1. Make sure your environment has the required dependencies:
   -  metpy
   -  xarray
   -  urllib
   -  numpy
   -  siphon
   -  cartopy
   -  datetime
   -  cdsapi
   -  pandas
2. ```
   pip install sounderpy
   ```
   Find it at https://pypi.org/project/sounderpy/1.0.0/
3. ```
   import sounderpy as spy
   ```
4. ```
   year  = '2014'
   month = '06'
   day   = '16'
   hour  = '20'
   latlon = [41.9, -97.01]
   method = 'rap'
   ```
5. ```
   raw_data = spy.get_model_data(method, latlon, year, month, day, hour)
   ```
6. ```
   clean_data = spy.parse_data(raw_data)
   ```
------
  and boom! Now you have a callable dictionary of vertical profile data including... 
1. Temperature
2. Dewpoint
3. Relative Humidity
4. Pressure
5. Height 
6. Height AGL
7. Vertical Velocity
8. U-component Wind 
9. V-component Wind

You can make a quick plot of the data using built-in MetPy plotting functions!, just call...
`spy.metpy_sounding(clean_data)`
<div align="center">
<img src="https://raw.githubusercontent.com/kylejgillett/sounderpy/main/images/example_RAP_0427201122z.png" width="600">
</div>

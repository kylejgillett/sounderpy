# SOUNDERPY DOCUMENTATION 
  > LATEST VERSION: v2.0.0 |  RELEASED: August 6, 2023  |  COPYRIGHT Kyle J Gillett, 2023

<div align="center">
<img src="https://github.com/kylejgillett/sounderpy/assets/100786530/2e9477c9-e36a-4163-accb-fe46780058dd" width="250">

## WELCOME! thank you for visiting SounderPy!
### [VISIT SOUNDERPY DOCUMENTATION HERE](https://github.com/kylejgillett/sounderpy/wiki)
This script is used to access vertical profile data for calculations or plotting of a vertical profile (sounding). 

[![PyPI Package](https://img.shields.io/pypi/v/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)](https://github.com/kylejgillett/sounderpy/blob/main/LICENSE.txt)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![GitHub commits](https://badgen.net/github/commits/kylejgillett/sounderpy)](https://GitHub.com/kylejgillett/sounderpy/commit/)
[![Maintainer](https://img.shields.io/badge/maintainer-kylejgillett-blue)](https://github.com/kylejgillett)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

</div>



-----
## What is SounderPy:

- SounderPy is a Python package used to access vertical profile data for calculations or plotting of a vertical profile (sounding). SounderPy's main use is for getting the data, but some basic plotting tools are included.

## SounderPy is currently capable of accessing and processing data from:

- ECMWF CDS ERA5 reanalysis [1940-present] *note: you must set up an account through the CDS to unlock ERA5 data. (see: https://cds.climate.copernicus.eu/api-how-to)
- UNIDATA THREDDS TDS RAP reanalysis [2005-present]
- UNIDATA THREDDS TDS RUC reanalysis [2005-2020]
- The University of Wyoming RAOB archive [1973-present, depending on station]
- Iowa State University's RAOB archive [1945-present, depending on station]
- The IGRAv2 Observed profile archive [1905-present, depending on station]
- Iowa State University's BUFKIT archive [2011-present, depending on station & model]
- Penn State University's BUFKIT archive [most-recent run, depending on station & model]
- UNIDATA THREDDS TDS RAP-most-recent-analysis [now, most recent analysis only]

## Why SounderPy?
- Sometimes data is tough to find, and often times is even tougher to get it in the format you like. SounderPy gets you this data!
- The code needed for loading and parsing vertical data (especially from models) can be large and messy. SounderPy keeps it hidden away in a PyPi package -- just import and call sounderPy functions to keep your code clean!

-------

## How to use SounderPy:
1. Make sure your environment has the required dependencies:
   - cdsapi>=0.6.1
   - matplotlib>=3.3.0, <=3.7.1
   - metpy>=1.5.1
   - netcdf4>=1.6.4
   - numpy>=1.20.0
   - pandas>=1.2.0
   - siphon>=0.9
   - scipy>= 1.10.1
   - xarray>=0.18.0

2. ```
   pip install sounderpy
   ```
   Find it at https://pypi.org/project/sounderpy/2.0.0/
3. ```
   import sounderpy as spy
   ```
4. ```
   year  = '2011'
   month = '04'
   day   = '27'
   hour  = '22'
   latlon = [33.19, -87.46]
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

# sounderpy

## WELCOME! thank you for visiting SounderPy!
This script is used to access vertical profile data for calculations or plotting of a vertical profile (sounding). 

### Why SounderPy?
+ Sometimes data is tough to find, and often times is even tougher to get in the format you like. SounderPy gets you this data!
+ The code needed for loading and parsing vertical data (especially from models) can be large and messy. SounderPy keeps it hidden away in a separate file -- just import and call sounderPy functions to keep your code clean!

## SounderPy is used for:
- Accessing and loading raw vertical profile data from the sources listed above
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

You can render a quick plot of the data using built-in MetPy plotting functions!, just call...
`spy.metpy_show_sounding(clean_data)`
...or you can save the plot to a file using...
`spy.metpy_save_sounding(clean_data, file_name)`
<div align="center">
<img src="https://raw.githubusercontent.com/kylejgillett/sounderpy/main/images/example_RAP_0427201122z.png" width="600">
</div>
<br>
<br>
<div align="center">
<img src="https://github.com/kylejgillett/sounderpy/assets/100786530/2e9477c9-e36a-4163-accb-fe46780058dd" width="250">
</div>

# sounderpy

## WELCOME! thank you for visiting SounderPy!
This script is used to access vertical profile data for calculations or plotting of a vertical profile (sounding). 

Thank you for visiting SounderPy. At this time [7/8/2023] this script is still under development and testing. This plan is to have the script uploaded to GitHub available for download by August 1st, 2023. 

### Why SounderPy?
+ Sometimes data is tough to find, and often times is even tougher to get in the format you like. SounderPy gets you this data!
+ The code needed for loading and parsing vertical data (especially from models) can be large and messy. SounderPy keeps it hidden away in a separate file -- just import and call sounderPy functions to keep your code clean!

### Using sounderPy is simple!
1. Download this script and add it to your directory (the goal is to have this loaded to GitHub as a true package in the future)
2. ```
   import sounderpy as spy
    ```
3. ```
   year  = '2014'
   month = '06'
   day   = '16'
   hour  = '20'
   latlon = [41.9, -97.01]
   method = 'rap'
   ```
4. ```
   raw_data = spy.get_model_data(method, 'point', latlon, year, month, day, hour)
   ```
5. ```
   clean_data = spy.parse_data(raw_data)
   ```
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
<img src="https://github.com/kylejgillett/sounderpy/assets/100786530/2e9477c9-e36a-4163-accb-fe46780058dd" width="300">
</div>

<div align="center">
<img src="https://github.com/kylejgillett/sounderpy/assets/100786530/2e9477c9-e36a-4163-accb-fe46780058dd" width="300">
</div>

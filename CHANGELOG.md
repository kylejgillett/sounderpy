
# SounderPy Change Log
All notable changes to this project will be documented in this file.

## [2.0.6] - LATEST RELEASE | 2023-10-04 - Release #4
 
SounderPy v2.0.6 will feature a few bug fixes to the plotting functions from v2.0.5. 
 
### Fixed
  - Erroneous dewpoint data from the NCEP-FNL and BUFKIT archived forecasts have been removed. I.e., bad data (considered less than -130C in this case) are set to nans to improve calculations.
  - MetPy CAPE/CIN calculations were changed from MetPy's `surface_based_cape_cin()` to MetPy's `cape_cin()` which seems to both perform better and handle possibly erroneous data better. This was done for SB, ML & MU CAPE/CIN.
  - An improved way to set the temperature axis bounds was created to ensure that the profile plots somewhat in the 'middle' of the skew-t for the best possible readability.

-----


## [2.0.5] - 2023-09-29 - Release #3
 
SounderPy v2.0.5 will feature a number of new tools, improvements to existing functionality, and a few bug fixes from v2.0.4. 
 
### New Features
- Access to Aircraft Communications, Addressing and Reporting System (ACARS) vertical profile data with `acars_data()`
     - Listing recent profiles for a given date and time: `acars_data(year, month, day, hour).list_profiles()`
     - Getting a profile after using `.list_profiles`: `acars_data(year, month, day, hour).get_profile(profile)`
- NCEP FNL 0.25deg Gridded Reanalysis Dataset using the `get_model_data()` function by setting the `method` kwarg to '`ncep`'.
- Ability to output SounderPy data to a SHARPPY input file. 

 
### Changed
- Finding station lat-lon data is now condensed into a single function: `get_latlon(station_type, station_id)` where `station_type` can be `buoy`, `raob`, `igra`, `metar`, or `bufkit`.
     - This change depreciates the following functions: `buoy_latlon()`, `metar_latlon()`, `raob_latlon()`, `igra_latlon()`, `metar_latlon()` and `bufkit_latlon()`.
- Outputting data to a file is now condensed into a single function: `to_file(file_type, filename)`, where `file_type` can be `csv`, `cm1`, or `sharppy`.
     - This change depreciates the following functions: `to_csv()` and `to_cm1`. 
- MetPy sounding and hodograph plots are now more advanced. 


### Fixed
  - Improved Docs
  - Corrected Bufkit data output lat-lon data

-----


 
## [2.0.4] - 2023-08-06 - Release #2
 
SounderPy v2.0.4 will feature a number of new tools, improvements to existing functionality, and minor bug fixes. 
 
### New Features
- IGRAv2 Archive Access via `get_obs_data()` -- just specify an IRGA station ID for the kwarg `station`!
- Most-recent RAP analysis data access via `get_model_data()` -- specify the `method` kwarg as `'rap-now'`
- Most-recent & archive BUFKIT data access via `get_bufkit_data()` -- GFS, NAM, NAMNEST, RAP, HRRR, SREF, & HIRESW data
- Ability to find a lat/lon pair of a US buoy/CMAN site -- `buoy_latlon('site-id')`
- Ability to find a lot/lon pair for a IGRA site -- `igra_latlon('site-id')`
- Ability to save plots to a file -- `metpy_sounding(clean_data, 'save')`
- Ability to save parsed data as a csv -- `to_csv(clean_data)`
- Ability to save parsed data to CM1 input file -- `to_cm1(clean_data)`
- Ability to plot profile data on a MetPy Hodograph! `metpy_hodograph(clean_data, 'show')`
- GitHub Wiki Documentation

 
### Changed
  - `metpy_sounding()` function now offers two rendering options, a user can specify kwarg `method` as `'show'` to display the plot inline, or as `'save'` to save the plot as a .png image. If `'save'` is chosen, the kwarg `filename` can be set to a user-specified file location and name.
    - Example: `metpy_sounding(clean_data, 'save', '/your-file-path/your-file-name')`
  - all `clean_data` dicts now include the key `site-info` which include site information, model/data information and time information.
 
### Fixed
  - Improved Docs
  - Made minor corrections to `get_docs()` function
  - Added `CHANGLOG.md` 
  - fixed `requirements.txt`, thus allowing dependencies to automatically load upon installing SounderPy

-----

## [1.0.0] - 2023-07-15 - Release #1
  
The initial release of SounderPy
 

## All Releases
- [sounderpy-2.0.6](https://pypi.org/project/sounderpy/2.0.6/)
  Minor: bug fixes to plotting functions from v2.0.5
- [sounderpy-2.0.5](https://pypi.org/project/sounderpy/2.0.5/)
  Major: new functions, improved existing functionality, repaired minor bugs
- [sounderpy-2.0.4](https://pypi.org/project/sounderpy/2.0.4/)
  Major: new functions, improved existing functionality, repaired minor bugs
- [sounderpy-1.0.0](https://pypi.org/project/sounderpy/1.0.0/)
  Initial release

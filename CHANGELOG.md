
# SounderPy Change Log
All notable changes to this project will be documented in this file. Apologies for the ridiculous versioning -- blame that on me not using tests and PyPi uploads correctly 🙃


## [3.0.1] PRE-RELEASE | Upcoming Release
Version 3.0.1 of SounderPy is bringing a number of minor and a couple of major bug fixes to release 3.0.0. This release will also feature a few new tools and improvements to existing functions. 

### The Scoop:

### Additions
   1. VAD Hodograph plot function that plots NEXRAD VWP VAD data on a hodograph: `build_vad_hodograph()`.
   2. Ability to target specific RAP-RUC datasets when using `get_model_data()` using the `dataset` kwarg.
   3. Ability to determine a 'box average' size when using `get_model_data()` using the `box_avg_size` kwarg. See the **Changed** section for more info.
   4. *IN THE WORKS*: New ECAPE operations, allowing for the plotting of different ECAPE parcel traces 
 
### Changed
   1. Model reanalysis data now uses an 'area-average' or 'box-sounding' approach to build a vertical profile of reanalysis data. I.e., by default, a user must pass a lat/lon point to `get_model_data()`, SounderPy will extract data for a 0.1x0.1 degree box and average the data together to build an averaged-point profile.

### Bug Fixes 
   1. Corrected EIL-SRH plot on hodographs
   2. Corrected 0-1km Streamwise Vorticity value on sounding and hodograph plots
   3. Improved Freezing-Point calculation using interpolation scheme.
   4. Corrected last-forecast-hour BUFKIT data retrieval issue and created an error message for when an invalid forecast hour is requested.
   5. Corrected RAOB site ID issue -- users can now request data using a WMO ID or ICAO ID.

### Removed 
  None


-------
-------
-------

## [3.0.0] LATEST RELEASE | 2024-01-08 - Release #5

## Summary
Version 3.0.0 of SounderPy brings exciting to features and improvements to existing functionality. This release will provide users with enhanced plotting capabilities, an improved user interface, improved parameter calculations via Sharppy utilities, and cleaner code overall. 

Some major upgrades include... a number of new and improved plotting tools that create one-of-a-kind sounding and hodograph plots, new parameter calculations utilizing [Sharrpy](https://github.com/sharppy/SHARPpy/blob/main/README.md) functionality and a new [documentation site](https://kylejgillett.github.io/sounderpy/).

### The Scoop:

### Additions
   1. New plot creation abilities using the new `build_sounding()` and `build_hodograph()` functions.
   2. Composite sounding plots using the new `build_composite()` function. This will allow users to analyze several profiles at once.
   3. A 'Dark Mode' setting for sounding, hodograph, and composite plots.
   4. A color-blind friendly setting for sounding plots that turns the green dewpoint trace to blue.
   5. Redesigned package structure, including the base sounderpy.py module, a calc.py module, and a plot.py module.
   6. A new stand-alone documentation site that offers clear, easy-to-find documentation for SounderPy's funtionality. [Available, under construction, here.](https://kylejgillett.github.io/sounderpy/)
   7. A print out of general thermodynamic and kinematic parameters of a given profile to the console when data is retrieved by SounderPy.  
 
### Changed
   1. Functions `build_sounding()` and `build_hodograph()` replace `metpy_sounding()` and `metpy_hodograph()`
   2. `get_model_data()` keyword argument `method` changed to `model`
   3.  Some warnings were changed to exceptions and other warnings were made more concise to address specific issues and improve the code's user-interface. 
   4. Computed thermodynamic and kinematic properties of a profile were changed to the widely accepted Sharppy calculation methods. This allows for increased reliability in functioning calculations and for more parameters that can be included in plots. SounderPy now directly uses the Sharrpy package for calculations of profile properties. 

### Bug Fixes 
   1. Minor bugs were addressed throughout the code. Some of these include...
        + The NWS-hosted METAR site list was removed from its online location a couple months ago. SounderPy now has its only METAR site list for finding METAR lat/lon sites. 
        + Some data retrieval functions would only like station ID's in upper-case ('APX'), this bug was fixed so lines like the following will be accepted: `spy.get_obs_data('apx', '2022', '05', '20', '18')`

### Removed 
  1. The `parse_data()` function is depreciated. Its functionality was simply added to the `get_model_data()` function. 

------

## Examples of new plots available with SounderPy functionality:

### New Full Sounding Plot:
![example-sounding_light](https://github.com/kylejgillett/sounderpy/assets/100786530/70d7209f-bc99-45c7-8999-952f26e60dba)

### New Composite Sounding Plot using 'Dark Mode'
![example-composite_dark](https://github.com/kylejgillett/sounderpy/assets/100786530/12d545d1-1518-4cbf-82d8-7062626c42e7)

### New Hodograph Plot in 'Dark Mode'
![example-hodograph_dark](https://github.com/kylejgillett/sounderpy/assets/100786530/3d3e3995-e915-49b5-850d-03c7981eccb6)


-----
-----
-----

## [2.0.6] - 2023-10-04 - Release #4
 
SounderPy v2.0.6 will feature a few bug fixes to the plotting functions from v2.0.5. 
 
### Fixed
  - Erroneous dewpoint data from the NCEP-FNL and BUFKIT archived forecasts have been removed. I.e., bad data (considered less than -130C in this case) are set to nans to improve calculations.
  - MetPy CAPE/CIN calculations were changed from MetPy's `surface_based_cape_cin()` to MetPy's `cape_cin()` which seems to both perform better and handle possibly erroneous data better. This was done for SB, ML & MU CAPE/CIN.
  - An improved way to set the temperature axis bounds was created to ensure that the profile plots somewhat in the 'middle' of the skew-t for the best possible readability.

-----
-----
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
-----
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
-----
-----

## [1.0.0] - 2023-07-15 - Release #1
  
The initial release of SounderPy
 

## All Releases
- [sounderpy-3.0.0](https://pypi.org/project/sounderpy/3.0.0/)
  Major: Major improvements to plot style, calculations, function names, repaired many bugs
- [sounderpy-2.0.6](https://pypi.org/project/sounderpy/2.0.6/)
  Minor: bug fixes to plotting functions from v2.0.5
- [sounderpy-2.0.5](https://pypi.org/project/sounderpy/2.0.5/)
  Major: new functions, improved existing functionality, repaired minor bugs
- [sounderpy-2.0.4](https://pypi.org/project/sounderpy/2.0.4/)
  Major: new functions, improved existing functionality, repaired minor bugs
- [sounderpy-1.0.0](https://pypi.org/project/sounderpy/1.0.0/)
  Initial release


# SounderPy Change Log
All notable changes to this project will be documented in this file.
 
## [Unreleased] - 2023-08-06 - Release #2
 
SounderPy v2.0.0 will feature a number of new tools, improvements to existing functionality, and minor bug fixes. 
 
### New Features
- IGRAv2 Archive Access via `get_obs_data()` -- just specify an IRGA station ID for the kwarg `station`!
- Most-recent RAP analysis data access via `get_model_data()` -- specify the `method` kwarg as `'rap-now'`
- Most-recent & archive BUFKIT data access via `get_bufkit_data()` -- GFS, NAM, NAMNEST, RAP, HRRR, SREF, & HIRESW data
- Ability to find a lat/lon pair of a US buoy/CMAN site -- `buoy_latlon('site-id')`
- Ability to find a lot/lon pair for a IGRA site -- `igra_latlon('site-id')`
- Ability to save plots to a file -- `metpy_sounding(clean_data, 'save')`
- Ability to save parsed data as a csv -- `to_csv(clean_data)`
- Ability to save parsed data to CM1 input file -- `to_cm1(clean_data)`


 
### Changed
  - `metpy_sounding()` function now offers two rendering options, a user can specify kwarg `method` as `'show'` to display the plot inline, or as `'save'` to save the plot as a .png image. If `'save'` is chosen, the kwarg `filename` can be set to a user-specified file location and name.
    - Example: `metpy_sounding(clean_data, 'save', '/your-file-path/your-file-name')`
  - Links to station lists (originally linked back to GitHub) are now packaged within the source-code directory
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
- [sounderpy-1.1.0]()
  Major: new functions, improved existing functionality, repaired minor bugs
- [sounderpy-1.0.0](https://pypi.org/project/sounderpy/1.0.0/)
  Initial release

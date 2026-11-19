"""
    SOUNDERPY DATA

    Purpose of module:

    The modules that get, decode, read, and write profile data:

    acars_data.py            : ACARS aircraft profiles from OU
    bufkit_data.py           : BUFKIT model forecast profiles from PSU and IEM
    ncep_bufr.py             : RRFS, HRRR, and RAP BUFR profiles from NCEP NOMADS
    obs_data.py              : RAOB and IGRA observed profiles
    model_data.py            : ERA5, NCEP-FNL, and RAP analysis profiles from THREDDS and CDS
    rapruc_grib_retriever.py : RAP/RUC GRIB downloads from AWS and NCEI
    rapruc_grib_parser.py    : RAP/RUC GRIB parsing
    radar_utils.py           : radar mosaic and single-site radar data for the map inset
    sfc_obs.py               : the surface observation nearest a sounding, from the IEM ASOS archive
    wrf_utils.py             : WRF output profiles
    cm1_utils.py             : CM1 input_sounding profiles
    file_io.py               : reading and writing SHARPpy, CSV, and CM1 files

    Docs:
        getting data:               https://sounderpy.app/docs/reference/gettingdata.html
        custom data sources:        https://sounderpy.app/docs/reference/customdatasources.html
        the clean_data dictionary:  https://sounderpy.app/docs/reference/clean_data_scheme.html

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""


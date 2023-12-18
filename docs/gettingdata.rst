Getting Data with SounderPy
===========================

SounderPy is currently capable of accessing and processing data from:

- ECMWF CDS ERA5 reanalysis [1940-present] *note: you must set up an account through the CDS to unlock ERA5 data. (see: https://cds.climate.copernicus.eu/api-how-to)*
- UNIDATA THREDDS TDS RAP reanalysis [2005-present]
- UNIDATA THREDDS TDS RUC reanalysis [2005-2020]
- The University of Wyoming RAOB archive [1973-present, depending on station]
- Iowa State University's RAOB archive [1945-present, depending on station]
- The IGRAv2 Observed profile archive [1905-present, depending on station]
- Iowa State University's BUFKIT archive [2011-present, depending on station & model]
- Penn State University's BUFKIT archive [most recent run, depending on station & model]
- UNIDATA THREDDS TDS RAP-most-recent-analysis [now, most recent analysis only]
- OU Aircraft Communications, Addressing and Reporting System (ACARS) [2019-present]
- NCEP FNL 0.25deg Gridded Reanalysis Dataset [2021-present]

|

***************************************************************


Model Reanalysis Data | RAP, ERA5, NCEP 
----------------------------------

- SounderPy hosts a simple function used to access model reanalysis profile data from the ERA5, RAP / RUC, & NCEP FNL datasets. 
- BEWARE: This data is *reanalysis*, therefore not a forecast & not entirely representative of the actual atmosphere. Understanding the caveats of using reanalysis model data is important when utilizing this function. 

We can use the simple ``spy.get_model_data()`` function:

.. py:function:: spy.get_model_data(model, latlon, year, month, day, hour)


   Return a ``dict`` of 'cleaned up' model reanalysis data from a given model, for a given location, date, and time

   :param model: required, the requested model to use (rap-ruc, era5, ncep)
   :type model: str
   :param latlon: required, the latitude & longitude pair for sounding ([44.92, -84.72])
   :type latlon: list
   :param year: required, valid year
   :type year: str
   :param month: required, valid month
   :type month: str
   :param day: required, valid day
   :type day: str
   :param hour: required, valid hour
   :type hour: str
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict


***************************************************************


|
|



Model Forecast Data | BUFKIT 
-----------------------------

- A function used to access BUFKIT model forecast vertical profile data
- This data is *model forecast* data. Users must note that BUFKIT data is model data loaded for *specific designated BUFKIT sites*. 
- To learn more about BUFKIT check out: https://meteor.geol.iastate.edu/~ckarsten/bufkit/bufkit.html

AVAILABLE BUFKIT SITES:

- Iowa State University's IEM BUFKIT warehouse: https://meteor.geol.iastate.edu/~ckarsten/bufkit/data/
- Penn State University's BUFKIT warehouse: http://www.meteo.psu.edu/bufkit/DomainNAMRAP_NAM_12.html


AVAILABLE MODEL DATA:

- Most recent model runs: 
    - GFS, NAM, NAMNEST, RAP, HRRR, SREF & HIRESW
    - via Penn State's BUFKIT Warehouse
- Archive model runs: 
    - GFS, NAM, NAMNEST, RAP, HRRR
    - via Iowa State's BUFKIT Warehouse


.. py:function:: spy.get_bufkit_data(model, station, fcst_hour, run_year=None, run_month=None, run_day=None, run_hour=None)

   Return a ``dict`` of 'cleaned up' model forecast data from a given model, for a given BUFKIT site identifier, forecast hour, & model-run-date

   :param model: required, the requested model to use (such as hrrr, nam, gfs, etc)
   :type model: str

   :param station: required, a 3-4 digit BUFKIT site identifier
   :type station: str

   :param fcst_hour: required, valid forecast hour
   :type fcst_hour: int

   :param year: required, valid year
   :type year: str

   :param month: required, valid month
   :type month: str

   :param day: required, valid day
   :type day: str

   :param hour: required, valid hour
   :type hour: str

   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict


***************************************************************

|
|

Observed Data | RAOB & IGRAv2
-----------------------------

- A function used to access and parse RAOB & IGRAv2 profile data. 
- This function will determine which dataset the user would like to access (RAOB from the University of Wyoming, or IGRAv2 from the IGRAv2 dataset) based on the provided station identifier, then search the appropriate dataset. 


.. py:function:: spy.get_obs_data(station, year, month, day, hour)

   Return a ``dict`` of 'cleaned up' model forecast data from a given model, for a given BUFKIT site identifier, forecast hour, & model-run-date

   :param station: required, a three digit RAOB identifier (such as: 'DTX') or 11 digit IGRAv2 identifier (such as: 'GMM00010393')
   :type station: str

   :param year: required, launch year
   :type year: str

   :param month: required, launch month
   :type month: str

   :param day: required, launch day
   :type day: str

   :param hour: required, launch hour
   :type hour: str

   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict

***************************************************************


|
|

Observed Data | ACARS Aircraft Obs
----------------------------------

- NOTE: this is a python class, not a function like the tools above. 
- This ``Class`` sets up a 'connection' to the ACARS data dataset. 
   - After setting up a 'connection' to the data, you can search for available profiles
   - 
- To learn more about ACARS, check out the 'AIRCRAFT' section of this webpage: NOAA Observing Systems


.. py:function:: spy.get_obs_data(station, year, month, day, hour)

   Return a ``dict`` of 'cleaned up' model forecast data from a given model, for a given BUFKIT site identifier, forecast hour, & model-run-date

   :param station: required, a three digit RAOB identifier (such as: 'DTX') or 11 digit IGRAv2 identifier (such as: 'GMM00010393')
   :type station: str

   :param year: required, launch year
   :type year: str

   :param month: required, launch month
   :type month: str

   :param day: required, launch day
   :type day: str

   :param hour: required, launch hour
   :type hour: str

   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict


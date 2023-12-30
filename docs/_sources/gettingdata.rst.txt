Tools for Getting Data
=======================

Getting data is tough, then getting into a usable format can often be as bad or worse. Thankfully this is SounderPy's main purpose and strong suit!

**SounderPy is currently capable of accessing and processing data from:**

+--------------------------------+-------------------+------------------+------------------+
|         **DATA**               |   **FUNCTION**    |  **TYPE**        |  **TIME RANGE**  |
+================================+===================+==================+==================+
|    ECMWF CDS ERA5 reanalysis*  | get_model_data()  | Reanalysis       |  1940-present    |
+--------------------------------+-------------------+------------------+------------------+
|     UNIDATA THREDDS TDS RAP    | get_model_data()  | Reanalysis       |  2005-present    |
+--------------------------------+-------------------+------------------+------------------+
|     UNIDATA THREDDS TDS RUC    | get_model_data()  | Reanalysis       |  2005-2020       |
+--------------------------------+-------------------+------------------+------------------+
|    UNIDATA THREDDS NCEP-FNL    | get_model_data()  | Reanalysis       |  2005-2020       |
+--------------------------------+-------------------+------------------+------------------+
|    ISU's BUFKIT archive        | get_bufkit_data() | Model Forecast   |  2011-present    |
+--------------------------------+-------------------+------------------+------------------+
|     PSU's BUFKIT feed          | get_bufkit_data() | Model Forecast   | Most recent runs |
+--------------------------------+-------------------+------------------+------------------+
|  UNIDATA THREDDS TDS RAP       | get_model_data()  | Model Analysis   | Most recent run  |
+--------------------------------+-------------------+------------------+------------------+
|   OU ACARS Archive             | acars_data()      | Observations     | 2019-present     |
+--------------------------------+-------------------+------------------+------------------+
|  The Unv. of WY RAOB Archive   | get_obs_data()    | Observations     | 1973-present     |
+--------------------------------+-------------------+------------------+------------------+
|  IGRAv2 Observation archive    | get_obs_data()    | Observations     |  1905-present    |
+--------------------------------+-------------------+------------------+------------------+

                    

***************************************************************

.. _modeldata:

Model Reanalysis Data | RAP, ERA5, NCEP 
----------------------------------------

**SounderPy hosts a simple function used to access model reanalysis profile data from the ERA5, RAP / RUC, & NCEP FNL datasets** 

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

.. note::
   BEWARE: This data is *reanalysis*, therefore not a forecast & not entirely representative of the actual atmosphere. Understanding the caveats of using reanalysis model data is important when utilizing this function. 

.. tip::
   You must set up an account through the CDS to unlock ERA5 data. See: https://cds.climate.copernicus.eu/api-how-to 

.. tip::
   **Is data access taking forever?** Sometimes the NCEP (RAP-RUC, NCEP-FNL) & ECMWF CDS (ERA5) servers are down and not able to be accessed. Sometimes these issues are resolved within hours, other times possibly a few days. 

***************************************************************




.. _bufkitdata:

Model Forecast Data | BUFKIT 
-----------------------------

**A function used to access BUFKIT model forecast vertical profile data**

.. tip:: 
   - This data is *model forecast* data. Users must note that BUFKIT data is model data loaded for *specific designated BUFKIT sites*
   - To learn more about BUFKIT check out: `IEM BUFKIT page <https://meteor.geol.iastate.edu/~ckarsten/bufkit/bufkit.html>`_

AVAILABLE BUFKIT SITES:
^^^^^^^^^^^^^^^^^^^^^^^^
- Iowa State University's `IEM BUFKIT warehouse <https://meteor.geol.iastate.edu/~ckarsten/bufkit/data/>`_
- Penn State University's `PSU BUFKIT warehouse <http://www.meteo.psu.edu/bufkit/DomainNAMRAP_NAM_12.html>`_

AVAILABLE MODEL DATA:
^^^^^^^^^^^^^^^^^^^^^^
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



.. _obsdata:

Observed Data | RAOB & IGRAv2
------------------------------

**A function used to access and parse RAOB & IGRAv2 profile data**
- This function will determine which dataset the user would like to access (RAOB from the University of Wyoming, or IGRAv2 from the IGRAv2 dataset) based on the provided station identifier, then search the appropriate dataset. 


.. py:function:: spy.get_obs_data(station, year, month, day, hour)

   Return a ``dict`` of 'cleaned up' observed profile data

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

.. note::
   Some data in these archives may be missing, incomplete or on occasion mislabled. 

***************************************************************



.. _acarsdata:


Observed Data | ACARS Aircraft Obs
-----------------------------------

- NOTE: this is a Python ``Class``, not a function like the tools above. 
   - This ``Class`` sets up a 'connection' to the ACARS data dataset. 
   - After setting up a 'connection' to the data, you can search for available profiles using the class's function, ``.list_profiles()``
   - Then you may select one of the listed profiles and use it as an argument for the class's function, ``.get_profile()``. See below.

- To learn more about ACARS, check out the 'AIRCRAFT' section of this webpage: `NOAA Observing Systems <https://www.weather.gov/about/observation-equipment>`_

.. class:: acars_data()

   :param year: required, observation year
   :type year: str
   :param month: required, observation month
   :type month: str
   :param day: required, observation day
   :type day: str
   :param hour: required, observation hour
   :type hour: str

   .. py:function:: spy.list_profiles()

      Return a list of strings that represents ACARS profiles for a given date and hour.

   .. py:function:: spy.get_profile(profile)

      Return a ``dict`` of 'cleaned up' ACARS observation profile data. Do so by selecting one of the profile strings listed by ``list_profiles()`` and pasting it as an argument in ``get_profile()``

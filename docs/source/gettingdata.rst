🌐 Tools for Getting Data
==========================

Getting data can be tough, and then getting it into a usable format can often be as bad or worse. Thankfully this is SounderPy's main purpose and strong suit!

**SounderPy is currently capable of accessing and processing data from:**

+--------------------------------+-----------------------+------------------+------------------+
|         **DATA**               |   **FUNCTION**        |  **TYPE**        |  **TIME RANGE**  |
+================================+=======================+==================+==================+
|    ECMWF CDS ERA5 reanalysis*  | get_model_data()      | Reanalysis       |  1940-present    |
+--------------------------------+-----------------------+------------------+------------------+
|     UNIDATA THREDDS TDS RAP    | get_model_data()      | Reanalysis       |  2005-present    |
+--------------------------------+-----------------------+------------------+------------------+
|     UNIDATA THREDDS TDS RUC    | get_model_data()      | Reanalysis       |  2005-2020       |
+--------------------------------+-----------------------+------------------+------------------+
|    UNIDATA THREDDS NCEP-FNL    | get_model_data()      | Reanalysis       |  2005-2020       |
+--------------------------------+-----------------------+------------------+------------------+
|    ISU's BUFKIT archive        | get_bufkit_data()     | Model Forecast   |  2011-present    |
+--------------------------------+-----------------------+------------------+------------------+
|     PSU's BUFKIT feed          | get_bufkit_data()     | Model Forecast   | Most recent runs |
+--------------------------------+-----------------------+------------------+------------------+
|  UNIDATA THREDDS TDS RAP       | get_model_data()      | Model Analysis   | Most recent run  |
+--------------------------------+-----------------------+------------------+------------------+
|   OU ACARS Archive             | acars_data()          | Observations     | 2019-06/2024     |
+--------------------------------+-----------------------+------------------+------------------+
|  The Unv. of WY RAOB Archive   | get_obs_data()        | Observations     | 1973-present     |
+--------------------------------+-----------------------+------------------+------------------+
|  IGRAv2 Observation Archive    | get_obs_data()        | Observations     |  1905-present    |
+--------------------------------+-----------------------+------------------+------------------+

                    

**********************************************************************

.. _modeldata:

Model Reanalysis Data | RAP, ERA5, NCEP 
----------------------------------------

**SounderPy hosts a simple function used to access model reanalysis profile data from the ERA5, RAP / RUC, & NCEP FNL datasets** 

This tool accesses pressure-level and surface model data, parses it using a 'box average approach', and creates a Python dictionary (referred to as a SounderPy clean_data dict in this documentation), of the vertical profile. 

When accessing RAP-RUC data, this tool will search for the given date/time in a list of RAP & RUC datasets available through NCEI -- each dataset does have varying output. Because of this, a new argument, `dataset`, allows you to target a specific dataset instead of searching for the first dataset with the desired date/time.

We can use the simple ``spy.get_model_data()`` function:

.. py:function:: spy.get_model_data(model, latlon, year, month, day, hour, dataset=None, box_avg_size=0.10, hush=False, clean_it=True)


   Return a ``dict`` of 'cleaned up' model reanalysis data from a given model, for a given location, date, and time

   :param model: the requested model to use (rap-ruc, era5, ncep)
   :type model: str, required
   :param latlon: the latitude & longitude pair for sounding ([44.92, -84.72])
   :type latlon: list, required
   :param year: valid year
   :type year: str, required
   :param month: valid month
   :type month: str, required
   :param day: valid day
   :type day: str, required
   :param hour: required, valid hour
   :type hour: str, required
   :param dataset: optional, target a specific dataset instead of searching for the first one with data.
   :type dataset: str, optional
   :param box_avg_size: optional, determine an area-averaged box size in degrees, default is 0.10 degrees.
   :type box_avg_size: int, optional
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting a data.
   :type hush: bool, optional, default is `False`
   :param clean_it: whether to return the raw_data object or a clean_data dict.
   :type clean_it: bool, optional, default is `True`
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict



.. _reanalysiskeys:

Model key names 
^^^^^^^^^^^^^^^
  - ``'era5'``: ECMWF ReAnalysis v5 (ERA5), reanalysis 

  - ``'rap'``, or ``'rap-ruc'``: RAPid refresh (RAP) / Rapid Update Cycle (RUC), reanalysis 

  - ``'ncep'``: NCEP Global Data Assimilation System/Final 0.25 degree, reanalysis 

  - ``'rap-now'``: RAPid refresh, latest analysis


.. _datasetkeys:

Dataset key names 
^^^^^^^^^^^^^^^^^^
  - ``'RAP_25km'``
  - ``'RAP_25km_old'``
      
  - ``'RAP_25km_anl'``
  - ``'RAP_25km_anl_old'``
          
  - ``'RAP_13km'``
  - ``'RAP_13km_old'``
      
  - ``'RAP_13km_anl'``
  - ``'RAP_13km_anl_old'``
      
  - ``'RUC_13km'``
  - ``'RUC_13km_old'``
      
  - ``'RUC_25km'``
  - ``'RUC_25km_old'``


.. _latlonpairs:

Latitude-Longitude pairs
^^^^^^^^^^^^^^^^^^^^^^^^^
  - A `list` of `floats`: ``[44.92, -84.72]``



.. note::
   **BEWARE** This data is *reanalysis*, therefore not a forecast & not entirely representative of the actual atmosphere. Understanding the caveats of using reanalysis model data is important when utilizing this function. 

.. tip::
   **To access ERA5 data** you must set up an account through the CDS to use their API: See: https://cds.climate.copernicus.eu/api-how-to 

.. tip::
   **Is data access taking forever?** Sometimes the NCEP (RAP-RUC, NCEP-FNL) & ECMWF CDS (ERA5) servers are down and not able to be accessed. Sometimes these issues are resolved within hours, other times possibly a few days. 






***************************************************************




.. _bufkitdata:

Model Forecast Data | BUFKIT 
-----------------------------

**A function used to access BUFKIT model forecast vertical profile data**

.. py:function:: spy.get_bufkit_data(model, station, fcst_hour, run_year=None, run_month=None, run_day=None, run_hour=None)

   Return a ``dict`` of 'cleaned up' model forecast data from a :ref:`given model<forecastmodels>`, for a given BUFKIT :ref:`site identifier<forecastsites>`, forecast hour, & model-run-date

   :param model: the model :ref:`'key' <forecastkeys>` name to request data from 
   :type model: str, required
   :param station: a 3-4 digit BUFKIT site identifier
   :type station: str, required
   :param fcst_hour: valid forecast hour
   :type fcst_hour: int, required
   :param year: valid year
   :type year: str, required
   :param month: valid month
   :type month: str, required
   :param day: valid day
   :type day: str, required
   :param hour: valid hour
   :type hour: str, required
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
   :type hush: bool, optional, default is `False`
   :return: :ref:`clean_data<datadescription>`, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict

.. _forecastsites:

Available BUFKIT Sites:
^^^^^^^^^^^^^^^^^^^^^^^^
.. raw:: html

    <embed>
        <iframe src="https://kylejgillett.github.io/sounderpy/bufkit_map" width="100%" height="500"></iframe>
    </embed>

-

.. _forecastmodels:

Available Models:
^^^^^^^^^^^^^^^^^^
- Most recent model runs: 
    - GFS, NAM, NAMNEST, RAP, HRRR, SREF & HIRESW
    - via Penn State's BUFKIT Warehouse
- Archive model runs: 
    - GFS, NAM, NAMNEST, RAP, HRRR
    - via Iowa State's BUFKIT Warehouse


.. _forecastkeys:

Model key names 
^^^^^^^^^^^^^^^
  - ``hrrr``: High Resolution Rapid Refresh, analysis (F00) & forecast
  - ``rap``: RAPid refresh, analysis (F00) & forecast
  - ``nam``: North American Mesoscale model, analysis (F00) & forecast
  - ``namnest``: Nested North American Mesoscale model, analysis (F00) & forecast
  - ``gfs``: Global Forecast System, analysis (F00) & forecast 
  - ``sref``: Short Range Ensemble Forecast, analysis (F00) & forecast
  - ``hiresw``: High RESolution Window forecast system, analysis (F00) & forecast


.. tip:: 
   - This data is *model forecast* data. Users must note that BUFKIT data is model data loaded for *specific designated BUFKIT sites*
   - To learn more about BUFKIT check out: `IEM BUFKIT page <https://meteor.geol.iastate.edu/~ckarsten/bufkit/bufkit.html>`_






***************************************************************



.. _obsdata:

Observed Data | RAOB & IGRAv2
------------------------------

**A function used to access and parse RAOB & IGRAv2 profile data**
- This function will determine which dataset the user would like to access (RAOB from the University of Wyoming, or IGRAv2 from the IGRAv2 dataset) based on the provided station identifier, then search the appropriate dataset. 


.. py:function:: spy.get_obs_data(station, year, month, day, hour)

   Return a ``dict`` of 'cleaned up' observed profile data

   :param station: a three digit RAOB identifier (such as: 'DTX') or 11 digit IGRAv2 identifier (such as: 'GMM00010393')
   :type station: str, required
   :param year: launch year
   :type year: str, required
   :param month: launch month
   :type month: str, required
   :param day: launch day
   :type day: str, required
   :param hour: launch hour
   :type hour: str, required
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
   :type hush: bool, optional, default is `False`
   :return: :ref:`clean_data<datadescription>`, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict

.. note::
   Some data in these archives may be missing, incomplete or on occasion mislabled. 


.. _raobsites:

Available RAOB Sites:
^^^^^^^^^^^^^^^^^^^^^^^^
.. raw:: html

    <embed>
        <iframe src="https://kylejgillett.github.io/sounderpy/raob_map" width="100%" height="500"></iframe>
    </embed>

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

   :param year: observation year
   :type year: str, required
   :param month: observation month
   :type month: str, required
   :param day: observation day
   :type day: str, required
   :param hour: observation hour
   :type hour: str, required

   .. py:function:: .list_profiles()

      Return a :ref:`list of strings<acarslists>` that represents ACARS profiles for a given date and hour.

   .. py:function:: .get_profile(profile)

      Return a ``dict`` of 'cleaned up' ACARS observation profile data. Do so by selecting one of the profile string :ref:`"IDs"<acarslists>` listed by ``list_profiles()`` and pasting it as an argument in ``get_profile()``

      :param profile: profile :ref:`"ID"<acarslists>`
      :type profile: str, required
      :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
      :type hush: bool, optional, default is `False`
      :return: :ref:`clean_data<datadescription>`, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
      :rtype: dict

.. _acarslists:

ACARS Data Retrieval Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Here is a simple example of the ACARS data retrieval functionality:

.. code-block:: python
   :linenos:

   # Start by setting up an 'ACARS connection'
   acars_conn = spy.acars_data('2023', '12', '30', '14')

   # List profiles
   acars_conn.list_profiles()

   '''
   `.list_profiles()` will return a list of all profiles available
   during the date/time entered in `acars_data()`, like this:
   ['ATL_1450',
   'AUS_1410',
   'AUS_1430',
   'AUS_1450',
   'BNA_1420',
   'BWI_1430']
   '''

   # To now get the data for a profile, 
   # copy the 'profile ID' and add it to `.get_profile()`:
   acars_conn.get_profile('AUS_1450')


.. note::
   ACARS data is aircraft observation data, thus these profiles are typically not 'full' profiles (i.e., up to 100 hPa). Often times these profiles extend to only 500 hPa or less. They may also contain various errors such as unreasonably dry dewpoints and unreasonably high wind velocities. 



****************************************


Observed Data | VAD radar data
-------------------------------

.. _vaddata:

.. py:function:: pyart_radar_profile(nexrad_site, scan_dt, from_file=False, data_file='none')

   Return a ``dict`` of 'cleaned up' radar VAD data. This radar data loader and VWP creator function is powered by PyArt 
   (https://arm-doe.github.io/pyart/)

   :param nexrad_site: station ID (``'KDTX'``)
   :type nexrad_site: str, required
   :param scan_dt: the date and time of the requested scan (``datetime(2021, 12, 11, 4, 24)``)
   :type scan_dt: datetime obj, required
   :param from_file: whether or not to search the NEXRAD AWS database or look for a local file, default is False
   :type from_file: bool, optional
   :param data_file: the filename of the local radar file to use
   :type data_file: str, optional




******************************************



.. _datadescription:

What does the data look like?
------------------------------
When using the data-retrevial functions above, they return 'clean_data', which is a Python Dictionary of vertical profile data and profile metadata. 

The profile data this `dict` contains...
   + ``clean_data['p']``: an `array` of pressure data 
   + ``clean_data['z']``: an `array` of height data 
   + ``clean_data['T']``: an `array` of temperature data
   + ``clean_data['Td']``: an `array` of dewpoint data
   + ``clean_data['u']``: an `array` of u-component of wind data 
   + ``clean_data['v']``: an `array` of v-component of wind data

The profile metadata this `dict` contains (via `clean_data['site_info']`)...
   + ``clean_data['site_info']['site-name']`` 
         - a `str` representing the name of a profile site, if available (e.g. 'DTX')
   + ``clean_data['site_info']['site-lctn']`` 
         - a `str` representing additional site location information (e.g. 'MI US')
   + ``clean_data['site_info']['site-latlon']`` 
         - a latitude-longitude pair of `floats` in a `list`
   + ``clean_data['site_info']['site-elv']`` 
         - elevation of the profile 
   + ``clean_data['site_info']['source']`` 
         - a `str` representing the data source name (e.g. 'RAOB OBSERVED PROFILE')
         - other sources are... 'ACARS OBSERVED AIRCRAFT PROFILE', 'BUFKIT FORECAST PROFILE', 'MODEL REANALYSIS PROFILE', 'RAOB OBSERVED PROFILE'
   + ``clean_data['site_info']['model']`` 
         - a `str` representing the model name, if available (e.g., 'no-model' or 'hrrr')
   + ``clean_data['site_info']['fcst-hour']`` 
         - if a model is used, the forecast hour of the model run as a `str` (e.g. 'no-fcst-hour' or 'F01')
   + ``clean_data['site_info']['run-time']`` 
         - if a model is used, the model run time as a `list` of `strs`
   + ``clean_data['site_info']['valid-time']`` 
         - the data's valid time as a `list` of `strs`

Below is an example:

      .. code-block:: python

         {'p': array([944. , 926.4, 925. , 894.5, 863.5, 850. , 848. , 833.4, 804.1,
                 795. , 775.7, 774. , 748. , 721.2, 720. , 700. , 685. , 674. ,
                 670. , 651. , 645.1, 630. , 621.3, 621. , 598.2, 591. , 587. ,
                 583. , 572. , 554. , 509. , 500. , 473.6, 471. , 446. , 442. ,
                 425. , 418. , 402. , 400. , 399. , 395. , 386. , 382. , 370. ,
                 354.7, 354. , 336. , 311.2, 300. , 297. , 279. , 250. , 241. ,
                 239. , 237.6, 232. , 200. , 194. , 190. , 188. , 170. , 168.9,
                 165. , 162. , 161. , 160.6, 155. , 152.8, 150. , 138.2, 135. ,
                 131.3, 131. , 130. , 127. , 125. , 124.8, 122. , 118.7, 118. ,
                 113. , 112. , 111. , 108. , 103. , 102. , 101. , 100. ]) <Unit('hectopascal')>,
          'z': array([  446,   610,   623,   914,  1219,  1356,  1376,  1524,  1829,
                  1926,  2134,  2152,  2438,  2743,  2757,  2990,  3168,  3300,
                  3349,  3584,  3658,  3850,  3962,  3966,  4267,  4364,  4419,
                  4473,  4625,  4877,  5542,  5680,  6096,  6137,  6550,  6617,
                  6911,  7035,  7323,  7360,  7378,  7452,  7620,  7696,  7925,
                  8230,  8243,  8612,  9144,  9400,  9470,  9900, 10640, 10880,
                 10935, 10973, 11128, 12070, 12260, 12389, 12454, 13067, 13106,
                 13246, 13356, 13394, 13411, 13628, 13716, 13830, 14326, 14466,
                 14630, 14645, 14690, 14831, 14927, 14935, 15075, 15240, 15278,
                 15544, 15599, 15655, 15826, 16123, 16184, 16247, 16310]) <Unit('meter')>,
          'T': array([ 26. ,  24.3,  24.2,  21.7,  19.2,  18. ,  17.4,  16.4,  14.3,
                  13.6,  13.4,  13.4,  10.9,   8.3,   8.2,   6.4,   5.2,   5.8,
                   6. ,   4.2,   3.6,   2. ,   2.4,   2.4,   0.2,  -0.5,  -0.7,
                  -0.5,  -0.9,  -3.1,  -9.1, -10.3, -14.1, -14.5, -16.9, -16.7,
                 -19.1, -19.7, -22.3, -22.5, -22.5, -22.5, -23.9, -24.5, -26.7,
                 -29.6, -29.7, -32.7, -36.1, -37.7, -37.9, -41.1, -47.1, -49.3,
                 -49.5, -49.8, -51.1, -59.3, -61.3, -62.3, -62.9, -68.1, -68.3,
                 -68.9, -68.3, -63.9, -63.8, -62.1, -62.8, -63.7, -68.5, -69.9,
                 -70.3, -70.3, -68.1, -66.5, -65.3, -65.3, -65.3, -63.8, -63.5,
                 -62.9, -61.9, -60.5, -59.1, -58.7, -57.5, -55.3, -55.3]) <Unit('degree_Celsius')>,
          'Td': array([ 17. ,  16.3,  16.2,  15.6,  15. ,  14.7,  14.8,  14.2,  13. ,
                  12.6,  10. ,   9.8,   8.7,   7.5,   7.4,   5.3,   4.1,  -1.2,
                  -3. ,  -3.8,  -3.6,  -3. ,  -4.5,  -4.6,  -4.4,  -4.3,  -5.3,
                  -8.5, -12.9, -14.1, -17.1, -17.3, -17.4, -17.4, -20.1, -22.7,
                 -26.1, -29.7, -31.3, -31.5, -31.5, -35.5, -37.6, -38.5, -36.8,
                 -34.5, -34.4, -36.4, -39.8, -41.4, -41.5, -45.7, -50.8, -53. ,
                 -54.3, -54.7, -56.1, -64.3, -66.3, -67.3, -66.9, -72. , -72.2,
                 -72.9, -72.5, -68.5, -68.4, -67.1, -67.8, -68.7, -73.5, -74.9,
                 -75.3, -75.3, -74.1, -74.5, -74.3, -74.4, -76.3, -76.5, -76.5,
                 -78.9, -78.9, -78.5, -79.1, -83.7, -83.5, -83.3, -83.3]) <Unit('degree_Celsius')>,
          'u': array([ 10.7246222 ,  10.60660172,  10.60660172,  17.        ,
                  22.36948102,  26.99707961,  26.99707961,  27.63986722,
                  31.81980515,  34.37362398,  39.83431104,  39.83431104,
                  42.13244437,  45.05336244,  45.05336244,  39.83431104,
                  39.99960775,  40.12982058,  40.22445359,  40.28302882,
                  40.30508653,  40.30508653,  40.30508653,  40.30508653,
                  55.92124435,  56.73165519,  56.73165519,  57.52478501,
                  57.50175672,  58.97894719,  60.00171105,  60.62177826,
                  64.08587988,  64.08587988,  58.51531863,  58.51531863,
                  55.35225748,  53.05840464,  49.9682747 ,  49.9682747 ,
                  49.9682747 ,  48.32997061,  44.23421039,  43.93899135,
                  44.16729559,  50.78742675,  50.78742675,  50.78742675,
                  51.60657879,  51.09549882,  51.09549882,  53.85980316,
                  57.09739058,  55.28477501,  54.37846722,  54.37846722,
                  55.28477501,  61.62892952,  64.34785288,  67.06677624,
                  67.97308403,  77.94246969,  78.84877747,  91.15018422,
                  99.6074178 , 102.42649567, 102.42649567,  80.39200027,
                  71.59831518,  69.53725394,  67.61480784,  52.13005469,
                  33.7059555 ,  34.47199994,  37.03650542,  45.28821067,
                  51.09549882,  51.09549882,  45.033321  ,  37.23909236,
                  37.60864741,  37.74069899,  38.27679749,  37.58770483,
                  37.48920614,  36.5444686 ,  36.63991854,  35.80278823,
                  35.86300913]) <Unit('knot')>,
          'v': array([ 8.99902654, 10.60660172, 10.60660172, 29.44486373, 31.94692973,
                 32.17386661, 32.17386661, 32.93991105, 31.81980515, 32.05392292,
                 33.4249557 , 33.4249557 , 35.35331853, 31.546704  , 31.546704  ,
                 33.4249557 , 34.77112854, 36.13305274, 37.5099098 , 38.90086875,
                 40.30508653, 40.30508653, 40.30508653, 40.30508653, 46.92349551,
                 45.94038855, 45.94038855, 44.9432877 , 43.33068167, 41.29750342,
                 36.05266524, 35.        , 37.        , 37.        , 36.56442923,
                 36.56442923, 35.94617631, 35.78834582, 34.98816262, 34.98816262,
                 34.98816262, 33.84100974, 30.97312756, 29.63722388, 25.5       ,
                 35.56173905, 35.56173905, 35.56173905, 36.13531549, 29.5       ,
                 29.5       , 28.63776533, 26.62495049, 25.77971397, 25.3570957 ,
                 25.3570957 , 25.77971397, 28.7380418 , 30.00589658, 31.27375137,
                 31.69636963, 36.34517051, 36.76778877, 33.1759539 , 36.25413519,
                 37.28019562, 37.28019562, 35.79282459, 33.38684268, 25.30949061,
                 18.11733316, 25.42552651, 28.28265483, 28.92544244, 28.93608934,
                 29.41050789, 29.5       , 29.5       , 26.        , 21.5       ,
                 20.84681367, 16.01997627, 14.69308593, 13.68080573, 10.74985688,
                  5.78807521,  5.14940474,  3.76302468,  3.13760674]) <Unit('knot')>,
          'site_info': {'site-id': 'KAPX',
           'site-name': 'GAYLORD',
           'site-lctn': 'MI US',
           'site-latlon': [44.92, -84.72],
           'site-elv': 446.0,
           'source': 'RAOB OBSERVED PROFILE',
           'model': 'no-model',
           'fcst-hour': 'no-fcst-hour',
           'run-time': ['no-run-time'],
           'valid-time': ['2022', '05', '20', '18']}}
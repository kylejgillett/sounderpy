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
|   OU ACARS Archive             | acars_data()          | Observations     | 2019-present     |
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

   :param model: the requested model to use ("rap-ruc", "era5", "ncep")
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
   :param dataset: target a specific dataset instead of searching for the first one with data ("rap-ruc" only).
   :type dataset: str, optional, default is `None`
   :param box_avg_size: determine an area-averaged box size, in degrees, by which gridded model data will be averaged to find a single vertical porfile.
   :type box_avg_size: int, optional, Default is `0.10`
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
   :type hush: bool, optional, default is `False`
   :param clean_it: whether to return the raw_data object or a clean_data dict.
   :type clean_it: bool, optional, default is `True`
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, omega & model information
   :rtype: dict



.. _reanalysiskeys:

Model key names 
^^^^^^^^^^^^^^^
  - ``'era5'``: ECMWF renalysis v5 (ERA5), reanalysis

  - ``'rap'``, or ``'rap-ruc'``: NCEP Rapid Refresh model (RAP) / Rapid Update Cycle model (RUC), reanalysis

  - ``'ncep'``: NCEP Global Data Assimilation System/Final 0.25 degree (ncep-fnl), reanalysis

  - ``'rap-now'``: NCEP Rapid Refresh model, latest analysis


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
   **To access ERA5 data** you -must- set API access to the ECMWF Climate Data Store (CDS). This includes....

    - creating a CDS API account
    - Setting up a CDS API personal access token
    - Creating a $HOME/.cdsapirc file

   Follow the instructions on the CDSAPI "how to" documentation -- See: https://cds.climate.copernicus.eu/how-to-api

.. tip::
   **Is data access taking forever?** Sometimes the NCEP (RAP-RUC, NCEP-FNL) & ECMWF CDS (ERA5) servers are down and not able to be accessed. Sometimes these issues are resolved within hours, other times possibly a few days. 






***************************************************************




.. _bufkitdata:

Model Forecast Data | BUFKIT 
-----------------------------

**A function used to access BUFKIT model forecast vertical profile data for a given BUFKIT site**

.. py:function:: spy.get_bufkit_data(model, station, fcst_hour, run_year=None, run_month=None, run_day=None, run_hour=None, hush=False, clean_it=True)

   Return a ``dict`` of 'cleaned up' model forecast data from a :ref:`given model<forecastmodels>`, for a given BUFKIT :ref:`site identifier<forecastsites>`, forecast hour, & model-run-date

   :param model: the model :ref:`'key' <forecastkeys>` name to request data from 
   :type model: str, required
   :param station: a 3-4 digit BUFKIT site identifier
   :type station: str, required
   :param fcst_hour: valid forecast hour
   :type fcst_hour: int, required
   :param run_year: valid year
   :type run_year: str, optional, Default=None
   :param run_month: valid month
   :type run_month: str, optional, Default=None
   :param run_day: valid day
   :type run_day: str, optional, Default=None
   :param run_hour: valid hour
   :type run_hour: str, optional, Default=None
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
   :type hush: bool, optional, default is `False`
   :param clean_it: whether to return the raw_data object or a clean_data dict.
   :type clean_it: bool, optional, default is `True`
   :return: :ref:`clean_data<datadescription>`, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, omega, & model information
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
  - ``hrrr``: High Resolution Rapid Refresh, analysis (F00) & forecast; out to forecast hour 48
  - ``rap``: Rapid Refresh Model, analysis (F00) & forecast; out to forecast hour 51
  - ``nam``: North American Mesoscale Model, analysis (F00) & forecast; out to forecast hour 48
  - ``namnest``: Nested North American Mesoscale model, analysis (F00) & forecast; out to forecast hour 60
  - ``gfs``: Global Forecast System, analysis (F00) & forecast; out to forecast hour 180
  - ``sref``: Short Range Ensemble Forecast, analysis (F00) & forecast; out to forecast hour 84
  - ``hiresw``: High Resolution Window Forecast System, analysis (F00) & forecast; out to forecast hour 48


.. tip::
   Running the ``get_bufkit_data()`` function without date kwargs will return the latest available forecast.
   Example:

   .. code-block:: python
       :linenos:

       # RAP model for site KGFK at forecast hour 5
       spy.get_bufkit_data('rap', 'kgfk', 5)



.. tip:: 
   - This data is *model forecast* data. Users must note that BUFKIT data is model data loaded for *specific designated BUFKIT sites*
   - To learn more about BUFKIT check out: `IEM BUFKIT page <https://meteor.geol.iastate.edu/~ckarsten/bufkit/bufkit.html>`_






***************************************************************



.. _obsdata:

Observed Data | RAOB & IGRAv2
------------------------------

**A function used to access and parse RAOB & IGRAv2 profile data**
- This function will determine which dataset the user would like to access (RAOB from the University of Wyoming, or IGRAv2 from the IGRAv2 dataset) based on the provided station identifier, then search the appropriate dataset. 


.. py:function:: spy.get_obs_data(station, year, month, day, hour, hush=False, clean_it=True)

   Return a ``dict`` of 'cleaned up' observed profile data

   :param station: may be a three digit RAOB identifier (such as: 'DTX'), 5 digit WMO identifier (such as: '72317'), or 11 digit IGRAv2 identifier (such as: 'GMM00010393')
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
   :param clean_it: whether to return the raw_data object or a clean_data dict.
   :type clean_it: bool, optional, default is `True`
   :return: :ref:`clean_data<datadescription>`, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & profile information
   :rtype: dict

.. note::
   Some data in these archives may be missing, incomplete or on occasion mislabled. If you can't find a RAOB you know for sure exists, try increasing or decreasing the launch_hour by 1 hour.


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

   .. py:function:: .get_profile(profile, hush=False, clean_it=True)

      Return a ``dict`` of 'cleaned up' ACARS observation profile data. Do so by selecting one of the profile string :ref:`"IDs"<acarslists>` listed by ``list_profiles()`` and pasting it as an argument in ``get_profile()``

      :param profile: profile :ref:`"ID"<acarslists>`
      :type profile: str, required
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
   :type hush: bool, optional, default is `False`
   :param clean_it: whether to return the raw_data object or a clean_data dict.
   :type clean_it: bool, optional, default is `True`
      :return: :ref:`clean_data<datadescription>`, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & profile/flight information
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
   + ``clean_data['omega']``: an `array` of vertical velocity -- model data only

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

**New to v3.0.5, profile metadata also contains pre-built plot titles (via `clean_data['titles']`). This will make creating titles manually for custom data sources easier.**



Below is an example:

      .. code-block:: python

        {'p': <Quantity([984.7 981.8 976.4 967.1 953.4 935.8 914.8 892.  867.2 839.6 809.2 775.8
          739.6 700.4 658.4 613.4 566.1 520.1 478.7 441.6 408.1 377.8 350.5 325.8
          303.6 283.6 265.5 248.7 233.5 219.6 205.7 191.7 177.7 163.8 149.8 135.9
          121.9 108.   94.   81.   70.6  62.2  54.3], 'hectopascal')>,
         'z': <Quantity([  262.29   287.38   334.72   417.39   540.9    701.94   897.92  1116.39
           1360.38  1638.85  1954.13  2311.53  2712.99  3165.6   3673.36  4247.56
           4889.99  5558.47  6202.    6817.56  7410.43  7981.44  8527.44  9050.51
           9548.   10021.78 10474.5  10917.67 11339.77 11745.44 12172.28 12626.51
          13109.32 13623.36 14184.85 14794.86 15466.97 16205.52 17052.44 17961.73
          18802.02 19579.02 20415.56], 'meter')>,
         'T': <Quantity([ 14.84  17.54  20.04  21.34  21.44  20.44  21.44  22.74  21.64  19.94
           17.54  14.84  11.94   8.64   5.24   1.84  -1.86  -6.06 -10.56 -14.66
          -18.36 -22.36 -26.66 -30.66 -34.06 -37.26 -40.06 -43.06 -45.96 -48.76
          -51.56 -54.46 -56.86 -58.26 -58.66 -59.66 -64.26 -65.26 -64.26 -64.66
          -63.86 -63.36 -62.16], 'degree_Celsius')>,
         'Td': <Quantity([   9.3     8.79    8.32    7.88    7.01    6.38   -2.2   -13.24  -14.49
           -24.86  -17.51   -9.57   -7.31   -9.9   -12.11  -12.99  -16.41  -20.62
           -25.37  -30.93  -36.12  -40.09  -42.47  -45.57  -49.39  -52.06  -54.3
           -55.81  -58.64  -62.27  -64.93  -65.45  -68.93  -74.24  -74.84 -101.09
          -101.6  -102.16 -102.8  -103.49 -104.11 -104.68 -105.28], 'degree_Celsius')>,
         'u': <Quantity([-2.33041109 -2.91567705 -2.3323691   0.77626972  3.5007069   5.63835178
           9.13497097 10.6914363  13.41432206 20.79403892 23.12721208 21.57710061
          20.41478107 21.38502534 22.54815854 23.51844683 24.29878342 23.71105505
          20.98924152 22.35724617 28.57257978 28.18931303 25.65927737 28.56994968
          28.57863291 29.73848546 31.87877692 31.48766321 29.74415098 27.98817845
          27.59730449 26.43360693 23.32832754 21.76679295 23.12850533 25.85409643
          33.62524652 19.4383033   8.94253878  9.13760323  6.0283743   6.02452813
           0.38886559], 'knot')>,
         'v': <Quantity([ 3.88459575  9.13602361 12.24993691 13.79818123 13.99892679 13.41310886
          10.68811982  9.33002089  9.33014275  6.99796725  2.91344493  1.36129697
          -1.35952707 -4.85855855 -5.83135029 -5.44263343 -4.85978647 -4.08282603
          -4.27412453 -5.44184192 -4.46923758 -0.19680173 -0.19257433 -4.0813693
          -2.52543878 -1.75116608 -3.49692463 -2.13557144  0.97112424  2.71960789
           2.13803291  0.5813991   1.55355538  5.2458388   6.21890192 20.41100923
          13.40831073  9.91285856 15.74811418  9.3310132   9.72277241  7.38416961
           6.99920592], 'knot')>,
         'omega': <Quantity([ 0.   0.   0.   0.   0.   0.1  0.1  0.1 -0.1 -0.3 -0.4 -0.2 -0.1  0.
           0.   0.1  0.1  0.   0.   0.   0.   0.   0.   0.1  0.1  0.   0.   0.
          -0.1  0.   0.   0.1  0.1  0.   0.   0.   0.   0.   0.   0.   0.   0.
           0. ], 'pascal / second')>,
         'site_info': {'site-id': 'KGFK',
          'site-name': 'GRAND FORKS INTL',
          'site-lctn': 'ND',
          'site-latlon': [47.95, -97.18],
          'site-elv': 257,
          'source': 'BUFKIT FORECAST PROFILE',
          'model': 'RAP',
          'fcst-hour': 'F00',
          'run-time': ['2024', '09', '28', '04'],
          'valid-time': ['2024', '09', '28', '04']},
         'titles': {'top_title': 'BUFKIT MODEL FORECAST PROFILE | 04Z RAP F00',
          'left_title': ' RUN: 09/28/2024 04Z  |  VALID: 09/28/2024 04Z',
          'right_title': 'KGFK - GRAND FORKS INTL, ND | 47.95, -97.18    '}}
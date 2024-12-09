✏️ Custom Data Sources
=======================

A neat and very useful aspect of SounderPy is the "``clean_data`` dictionary" that stores sounding data. This dictionary is fed into a SounderPy plotting function. The plotting functions are designed in such a way that *any* data can be passed to it, so long as the data maintains the SounderPy "``clean_data`` dictionary" format.

***************************************************************

Built-in Custom Data Ingestion Tools
-------------------------------------
As of v3.0.6, SounderPy has some integrated custom data ingestion tools. These include ``make_wrf_profile()`` & ``make_cm1_profile()``.


Ingesting WRF Output Data
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:function:: spy.make_wrf_profile(ds, latlon, model_name='WRF-ARW', run_name='CONTROL RUN')


   Return a ``dict`` of 'cleaned up' WRF output profile data for a given location

   :param ds: a netcdf4 dataset object of your wrfout*.nc file
   :type ds: netCDF4.Dataset(), required
   :param latlon: the latitude & longitude pair for sounding (ex: [44.92, -84.72])
   :type latlon: list, required
   :param model_name: the name of your model for plot titles
   :type model_name: str, recommended, default is "WRF-ARW"
   :param run_name: the name of your run for plot titles
   :type run_name: str, optional, default is "CONTROL"
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict


.. _wrfexample:
.. code-block:: python

   import sounderpy as spy
   import netcdf4

   # open your `wrfout*.nc` file as a netCDF4 dataset
   wrf_ds = netCDF4.Dataset(filename_of_wrf_data)

   # use `make_wrf_profile()` to create a `clean_data` dictionary of sounding data
   clean_data = make_wrf_profile(wrf_ds, [48.978, -100.861])

   # Now, pass the 'clean_data' dictionary into the SounderPy `build_sounding` function
   spy.build_sounding(clean_data)



Ingesting CM1 input_sounding Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:function:: spy.make_cm1_profile(filename, meta_data_dict)

   Return a ``dict`` of CM1 `input_sounding` data

   :param filename: the filename of an `input_sounding` file
   :type filename: str, required
   :param meta_data_dict: a simple dictionary of metadata for the profile (see :ref:`'CM1 Example' <cm1example>`)
   :type meta_data_dict: dict, required
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict


.. _cm1example:
.. code-block:: python

   import sounderpy as spy

   # create this simple dict of information for your profile
   # latitude, longitude, elevation, and top, right & left titles are required
   meta_data_dict = {
    'latlon': [45.100, -100.89],
    'elev': 358,
    'top_title': f"CUSTOM CM1 COMPOSITE SOUNDING",
    'left_title': f"Near-storm inflow RUC sounding for the El-Reno EF3 tornado",
    'right_title': f"May 24, 2011"}

   # use `make_cm1_profile()` to create a `clean_data` dictionary of sounding data
   clean_data = make_cm1_profile(input_sounding_filename, meta_data_dict)

   # Now, pass the 'clean_data' dictionary into the SounderPy `build_sounding` function
   spy.build_sounding(clean_data)

***************************************************************




Ingesting your own data into SounderPy
--------------------------------------

If you have custom data that you want to plot with SounderPy, you can just manually create a "``clean_data`` dictionary" and pass it to a plot!

Your data could be field campaign observations, custom model output, university soundings, and anything else in between. So long as you set it up right, it can be integrated into SounderPy.


Understanding the SounderPy data format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

What's contained in the ``clean_data`` dictionary...

The actual profile data...
   + ``clean_data['p']``: an `array` of pressure data
   + ``clean_data['z']``: an `array` of height data
   + ``clean_data['T']``: an `array` of temperature data
   + ``clean_data['Td']``: an `array` of dewpoint data
   + ``clean_data['u']``: an `array` of u-component of wind data
   + ``clean_data['v']``: an `array` of v-component of wind data
   + ``clean_data['omega']``: an `array` of vertical velocity -- model data only

The profile metadata (via `clean_data['site_info']`)...
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

And finally, ``clean_data`` also contains pre-built plot titles (via `clean_data['titles']`).
   + ``clean_data['titles']['top_title']``
         - the primary, or "main" title that appears above everything else
   + ``clean_data['titles']['left_title']``
         - the left title that typically contains time/date information
   + ``clean_data['titles']['right_title']``
         - the right title that typically contains location information


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



**************************************************



Building a custom clean_data dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After importing your custom sounding data, (perhaps using ``pandas`` or ``numpy``) you can define a ``clean_data`` dictionary with it.

Here is an example. ``raw_data`` represents some obj with your data. Using ``raw_data[2]`` is simply referencing some "column of data", in whatever format ``raw_data`` is. I.e., if the 2nd index in raw_data holds pressure values, ``raw_data[2]`` goes into ``clean_data['p']``.


      .. code-block:: python

        # declare the clean_data dictionary, leave it empty for the moment
        clean_data = {}

        # add profile data | make sure you have p, z, T, Td, u, & v
        # use metpy.units to add units to each array -- make sure they are in the same
        # units as show below!
        clean_data['p']  = np.array(raw_data[2])*units.hPa
        clean_data['z']  = np.array(raw_data[7])*units.m
        clean_data['T']  = np.array(raw_data[3])*units.degC
        clean_data['Td'] = np.array(raw_data[11])*units.degC
        clean_data['u']  = np.array(raw_data[5])*units.kts
        clean_data['v']  = np.array(raw_data[6])*units.kts

        # declare some profile metadata
        clean_data['site_info'] = {
                    'site-id'     : 'UND',                              # could be a station ID, site ID, launch ID, mission ID, etc
                    'site-name'   : 'GRAND FORKS'                       # a location's "name", usually the city or town
                    'site-lctn'   : 'ND',                               # could be another name, or the state
                    'site-latlon' : [47.9213, -97.087]                  # location lat/lon, list of floats
                    'site-elv'    : 257,                                # the profile's elevation in meters (int or float) (easily found on google)
                    'source'      : 'UND AEROSPACE'                     # the 'source' which will be the main title component of the plot
                    'model'       : 'none',                             # model name if a model was involved
                    'fcst-hour'   : f'none',                            # forecast hour if a model was involved
                    'run-time'    : ['none', 'none', 'none', 'none'],   # model run date if a model was involved
                    'valid-time'  : ["2024", "09", "28", "16:15"]}      # the profile's valid date/time.

        # declare the plot titles
        clean_data['titles'] = {
                    'top_title': 'UNIVERSITY OF NORTH DAKOTA AEROSPACE | OBSERVED SOUNDING',
                    'left_title': 'VALID: 09/28/2024 - 16:15Z',
                    'right_title': 'GRAND FORKS, ND [47.9213, -97.087]    '}



Done! Now your custom data source is integrated into a SounderPy ``clean_data`` dictionary.

Finally, you can pass it into the ``build_sounding()`` to visualize your data.

.. code-block:: python

    spy.build_sounding(clean_data *kwargs)

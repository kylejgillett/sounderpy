✏️ Custom Data Sources
=======================

A neat and very useful aspect of SounderPy is the "``clean_data`` dictionary" that stores sounding data. This dictionary is fed into a SounderPy plotting function. The plotting functions are designed in such a way that *any* data can be passed to it, so long as the data maintains the SounderPy "``clean_data`` dictionary" format.

***************************************************************


Built-in Custom Data Ingestion Tools
-------------------------------------
Have model data or your own observations? Pipe them into SounderPy with these custom data tools and methods!

.. note::
   ``from_wrf()`` and ``from_cm1()`` were called ``make_wrf_profile()`` and ``make_cm1_profile()`` before SounderPy 3.3.0. The old names still work, with a ``FutureWarning``, and will be removed in a later release.


Ingesting WRF Output Data
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:function:: spy.from_wrf(ds, latlon, model_name='WRF-ARW', run_name='CONTROL RUN', box_avg_size=0.10, time_idx=0)


   Return a ``dict`` of 'cleaned up' WRF output profile data for a given location

   :param ds: a netcdf4 dataset object of your wrfout*.nc file (an opened dataset, not a file path)
   :type ds: netCDF4.Dataset(), required
   :param latlon: the latitude & longitude pair for sounding (ex: [44.92, -84.72])
   :type latlon: list, required
   :param model_name: the name of your model for plot titles
   :type model_name: str, recommended, default is "WRF-ARW"
   :param run_name: the name of your run for plot titles
   :type run_name: str, optional, default is "CONTROL RUN"
   :param box_avg_size: the width, in degrees, of the box of grid points averaged into the profile. A box smaller than one grid cell uses the nearest grid point.
   :type box_avg_size: float, optional, default is ``0.10``
   :param time_idx: the index of the output time to use, for ``wrfout`` files with more than one time
   :type time_idx: int, optional, default is ``0``
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict


.. _wrfexample:
.. code-block:: python

   import sounderpy as spy
   import netcdf4

   # open your `wrfout*.nc` file as a netCDF4 dataset
   wrf_ds = netCDF4.Dataset(filename_of_wrf_data)

   # use `from_wrf()` to create a `clean_data` dictionary of sounding data
   clean_data = spy.from_wrf(wrf_ds, [48.978, -100.861])

   # Now, pass the 'clean_data' dictionary into the SounderPy `build_sounding` function
   spy.build_sounding(clean_data)

******************************************



Ingesting CM1 input_sounding Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:function:: spy.from_cm1(filename, meta_data_dict)

   Return a ``dict`` of CM1 ``input_sounding`` data

   :param filename: the filename of an ``input_sounding`` file
   :type filename: str, required
   :param meta_data_dict: a simple dictionary of metadata for the profile (see :ref:`'CM1 Example' <cm1example>`)
   :type meta_data_dict: dict, required
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict


.. _cm1example:
.. code-block:: python

   import sounderpy as spy

   # create this simple dict of information for your profile
   # latitude & longitude are required. elevation (m) defaults to 0 m and the
   # titles default to "CM1 INPUT SOUNDING" and blank left & right titles
   meta_data_dict = {
    'latlon': [45.100, -100.89],
    'elev': 358,
    'top_title': f"CUSTOM CM1 COMPOSITE SOUNDING",
    'left_title': f"Near-storm inflow RUC sounding for the El-Reno EF3 tornado",
    'right_title': f"May 24, 2011"}

   # use `from_cm1()` to create a `clean_data` dictionary of sounding data
   clean_data = spy.from_cm1(input_sounding_filename, meta_data_dict)

   # Now, pass the 'clean_data' dictionary into the SounderPy `build_sounding` function
   spy.build_sounding(clean_data)


.. note::
   CM1 input_sounding height data *should* be "AGL", or "above ground level". This means you'll need an elevation value to plot your profile correctly. Find the elevation of where this profile is "supposed to be" and add it to the ``meta_data_dict``. Without one, SounderPy warns and places the surface at 0 m.


***************************************************************



Ingesting SHARPPY Files
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. py:function:: spy.from_sharppy(filepath)


   Return a ``dict`` of 'cleaned up' WRF output profile data for a given location

   :param filepath: filepath to your SHARPPY data file
   :type filepath: str, required
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind 
   :rtype: dict


.. note::
   To plot your SHARPPY file data with SounderPy, you'll have to manually provide additional site meta-data. Be sure to include a site ID, name, location, lat/lon, and elevation. You must also create plot titles for the sounding plot.


.. _sharppyexample:
.. code-block:: python

   import sounderpy as spy

   # set file path
   file_path = "Downloads/sounding_data/09082024_04Z"

   # open your data
   clean_data = spy.from_sharppy(file_path)

   # manually set profile meta-data
   data['site_info'] = {'site-id': 'KGFK',
                        'site-name': 'GRAND FORKS INTL',
                        'site-lctn': 'ND',
                        'site-latlon': [47.95, -97.18],
                        'site-elv': 257,
                        'source': 'BUFKIT FORECAST PROFILE',
                        'model': 'RAP',
                        'fcst-hour': 'F00',
                        'run-time': ['2024', '09', '28', '04'],
                        'valid-time': ['2024', '09', '28', '04']}

   data['titles'] = {'top_title': 'BUFKIT MODEL FORECAST PROFILE | 04Z RAP F00',
                     'left_title': ' RUN: 09/28/2024 04Z  |  VALID: 09/28/2024 04Z',
                     'right_title': 'KGFK - GRAND FORKS INTL, ND | 47.95, -97.18    '}

   # Now, pass the 'clean_data' dictionary into the SounderPy `build_sounding` function
   spy.build_sounding(clean_data)


***************************************************************



Ingesting GRAW, iMet, and Windsond Radiosonde Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SounderPy can read the files saved by three common research radiosonde systems. These readers are based on conversion scripts contributed by Alec Sczepanski (University of North Dakota). Each returns the same ``clean_data`` dict as ``from_sharppy``, keeps only the balloon's ascent (from launch to burst), and fills in whatever site meta-data the file holds.

.. py:function:: spy.from_graw(filepath, date=None)


   Return a ``dict`` of 'cleaned up' profile data from a GRAW radiosonde file

   :param filepath: filepath to the tab-separated table saved by GRAW's software
   :type filepath: str, required
   :param date: the launch date, as ``'YYYY-MM-DD'``, ``'YYYYMMDD'``, or a ``datetime.date``. GRAW files hold the time of day but not the date; without ``date``, it is read from a file name that starts with ``YYYYMMDD`` (e.g. ``20260125_2100.txt``)
   :type date: str or datetime.date, optional
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind
   :rtype: dict

.. py:function:: spy.from_imet(filepath)


   Return a ``dict`` of 'cleaned up' profile data from an iMet radiosonde file (InterMet's SkySonde software)

   :param filepath: filepath to the SkySonde flight data .csv
   :type filepath: str, required
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind
   :rtype: dict

.. py:function:: spy.from_windsond(filepath, elevation=None)


   Return a ``dict`` of 'cleaned up' profile data from a Windsond radiosonde file

   :param filepath: filepath to the Windsond .sounding.csv
   :type filepath: str, required
   :param elevation: the launch site's elevation (m). Windsond heights are above ground level; with ``elevation``, they are put above sea level like other SounderPy data
   :type elevation: int or float, optional
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind
   :rtype: dict


.. note::
   - Rows are kept from launch to burst: iMet rows before launch (ascent rate under 1 m/s) and after the burst are dropped, and so is any level whose pressure doesn't fall below the levels under it.
   - GRAW's first row is the ground station's surface observation. iMet and Windsond dewpoints are calculated from temperature and relative humidity.
   - Missing values (e.g. iMet winds before the GPS locks on) are kept as NaN.
   - As with SHARPpy files, add the site meta-data the file doesn't hold (each reader prints what's still needed) and plot titles before plotting.


.. code-block:: python

   import sounderpy as spy

   # open your data (the date comes from the file name)
   clean_data = spy.from_graw("Downloads/20260125_2100.txt")

   # add the site meta-data the file doesn't have, and plot titles
   clean_data['site_info'].update({'site-id': 'UND',
                                   'site-name': 'UND GRAND FORKS',
                                   'site-lctn': 'ND'})

   clean_data['titles'] = {'top_title': 'GRAW RADIOSONDE OBSERVED VERTICAL PROFILE',
                           'left_title': 'LAUNCH: 01/25/2026 20:26Z',
                           'right_title': 'UND - UND GRAND FORKS, ND | 47.92, -97.09    '}

   spy.build_sounding(clean_data, mode='winter')



---------------------------------------------------------------------



Creating Custom ``clean_data``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given SounderPy's data structure (See :ref:`clean_data schema <clean_data_scheme>`), it's possible to ingest any data so long as its in the correct format. You can just manually create the entire "``clean_data`` dictionary". Your data could be field campaign observations, custom model output, university soundings, and anything else in between. So long as you set it up right, it can be integrated into SounderPy. When constructing a profile manually, the safest approach is to begin with the
six required atmospheric variables and ``site_info``.

At minimum:


.. code-block:: python

   clean_data = {
       "p": ...,
       "z": ...,
       "T": ...,
       "Td": ...,
       "u": ...,
       "v": ...,
       "site_info": {...},
   }

Before passing custom data into SounderPy, verify that:

* all six profile arrays are one-dimensional;
* all six arrays have equal length;
* all six arrays contain physical units;
* pressure is physically valid;
* the profile is ordered upward;
* ``site_info`` is a dictionary.



**Here is a full example**. Here, ``raw_data`` represents a data object with your data. Using ``raw_data[2]`` is simply referencing some "column of data", in whatever format ``raw_data`` is. I.e., if the 2nd index in raw_data holds pressure values, ``raw_data[2]`` goes into ``clean_data['p']``.


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
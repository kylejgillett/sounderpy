🛠️ Helper Tools
================

SounderPy includes several helper functions for profile calculations,
interpolation, file export, station lookup, and manipulation of ``clean_data``
profiles.

See :ref:`clean_data schema <clean_data_scheme>` for the common profile
structure used by these tools.



*************************


Returning a dictionary of profile parameters
--------------------------------------------

Calculate SounderPy's thermodynamic, kinematic, and interpolated profile parameters for a given ``clean_data`` profile.

The calculation returns a dictionary with four dictionaries inside it:

- ``general``: general profile information and derived quantities;
- ``thermo``: thermodynamic parameters;
- ``kinem``: kinematic parameters; and
- ``intrp``: an interpolated version of the profile.

.. class:: sounding_params(clean_data, storm_motion='right_moving', modify_sfc=None, include_all_parcels=False, parcel_depth=None)

   :param clean_data: the dictionary of profile data to calculate profile parameters for (see :doc:`/reference/gettingdata`)
   :type clean_data: dict, required
   :param storm_motion: the storm motion used for calculations. Default is 'right_moving'. Custom storm motions are accepted as a ``list`` of ``floats`` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional
   :param modify_sfc: a ``dict`` in the format ``{'T': 25, 'Td': 21, 'ws': 20, 'wd': 270}`` to modify the surface values of the ``clean_data`` dict.
   :type modify_sfc: None or dict, optional, default is None
   :param include_all_parcels: whether to lift a parcel from every level of the profile, which adds CAPE, CIN, and 3CAPE profiles (``cape_profile``, ``cin_profile``, and ``3cape_profile``) to ``thermo``. This is slower.
   :type include_all_parcels: bool, optional, default is ``False``
   :param parcel_depth: with ``include_all_parcels``, only lift parcels from levels up to this height (m AGL). Levels above it are ``nan``. This is much faster when only the lowest few kilometers are needed.
   :type parcel_depth: int or float, optional, default is ``None``, which lifts a parcel from every level


   .. py:function:: .calc()

      :return: A dictionary with the keys ``general``, ``thermo``, ``kinem``, and ``intrp``, each holding a dictionary of profile data.
      :rtype: dict

Example:

.. code-block:: python

   params = spy.sounding_params(clean_data).calc()

   params['thermo']['mlcape']
   params['kinem']['srh_0_to_1000']

.. versionchanged:: 3.3.0
   ``.calc()`` returns a dictionary. It used to return a tuple, which was unpacked with ``general, thermo, kinem, intrp = spy.sounding_params(clean_data).calc()``.

.. note::
   These calculated parameters are the same used on SounderPy plots. This function simply returns the master dictionary of all of these values.


****************************************************************


Returning a dictionary of winter parameters
-------------------------------------------

Calculate everything ``sounding_params`` does, plus a fifth dictionary, ``winter``, with the values behind the winter sounding (``build_sounding(mode='winter')``):

- precipitation type: Bourgouin's method (``ptype``) and SHARPpy's best guess (``sharppy_ptype``), the precipitation source level, and the phase of falling precipitation from the source to the surface (``phase_trace``);
- melting and refreezing layers and their energies (``layers``), freezing levels, and the warm nose;
- cloud and icing layers;
- snow growth: the dendritic growth zone in cloud (``dgz_z``, ``dgz_p``, ``dgz_depth``), its RH and omega, OPRH, and the snow-to-liquid ratio (``slr``);
- the mixed layer: depth, lapse rate, and mean wind (``ml_depth``, ``ml_lapse_rate``, ``mean_wind_ml``);
- surface and hazard values: wind chill, PBL max gust, the snow squall parameter, surface inversion, and freezing fog and freezing drizzle potential (``fzfg``, ``fzdz``);
- winds: 10 m, ML and DGZ mean winds, the 0-2 km max wind, and the winds at 850, 700, and 500 hPa.

.. class:: winter_params(clean_data, storm_motion='right_moving', modify_sfc=None, include_all_parcels=False, parcel_depth=None)

   :param clean_data: the dictionary of profile data to calculate profile parameters for (see :doc:`/reference/gettingdata`)
   :type clean_data: dict, required
   :param storm_motion: the storm motion used for calculations, as in ``sounding_params``.
   :type storm_motion: str or list of floats, optional
   :param modify_sfc: a ``dict`` in the format ``{'T': -2, 'Td': -3, 'ws': 20, 'wd': 20}`` to modify the surface values of the ``clean_data`` dict. The modified surface is used for both the standard and the winter values.
   :type modify_sfc: None or dict, optional, default is None
   :param include_all_parcels: as in ``sounding_params``.
   :type include_all_parcels: bool, optional, default is ``False``
   :param parcel_depth: as in ``sounding_params``.
   :type parcel_depth: int or float, optional, default is ``None``


   .. py:function:: .calc()

      :return: A dictionary with the keys ``general``, ``thermo``, ``kinem``, ``intrp``, and ``winter``, each holding a dictionary of profile data.
      :rtype: dict

Example:

.. code-block:: python

   params = spy.winter_params(clean_data).calc()

   params['winter']['ptype']        # e.g. 'FZRA'
   params['winter']['dgz_depth']    # m
   params['thermo']['mlcape']       # every sounding_params value is here too

.. versionadded:: 3.3.0

.. note::
   These are the values drawn on the winter sounding. Values that can't be calculated for a profile, such as omega values for data without omega (RAOBs), are masked.


****************************************************************




Saving data to a file
----------------------

.. py:function:: spy.to_file(file_type, clean_data, filename=None, convert_to_AGL=True)


   Create a file of 'cleaned' SounderPy data

   :param file_type: a ``str`` representing the file type you'd like to export data to.
   :type file_type: str, required
   :param clean_data: 'cleaned' SounderPy data ``dict``
   :type clean_data: dict, required
   :param filename: output filename. If omitted, SounderPy uses its defaultoutput filename.
   :type filename: str, optional, default is ``None``
   :param convert_to_AGL: whether or not to convert height values to "above ground level" when saving to a file. Useful for CM1.
   :type convert_to_AGL: bool, optional, default is ``True``
   :return: a file of SounderPy data.

Example:

- Supported ``file_type`` values are:

  - ``"csv"``
  - ``"sharppy"``
  - ``"cm1"``

.. code-block:: python

   spy.to_file("csv", clean_data, filename="sounding.csv")

   spy.to_file("sharppy", clean_data, filename="sounding.snd")

   spy.to_file("cm1", clean_data, filename="input_sounding")

See :ref:`tutorial-exporting-data` for complete export examples.

***************************************************************


Merging two profiles via weighted averaging
--------------------------------------------

Merge two profiles together using weighted averaging of values along a common height axis to return a single profile.

.. py:function:: spy.merge_profiles(profile_a, profile_b, weight_a)

   :param profile_a: a SounderPy ``clean_data`` dictionary containing variable arrays for profile_a.
   :type profile_a: dict, required
   :param profile_b: a SounderPy ``clean_data`` dictionary containing variable arrays for profile_b.
   :type profile_b: dict, required
   :param weight_a: the weight for ``profile_a`` in the averaging scheme, default is 0.5 (50%)
   :type weight_a: float, optional
   :return: a ``clean_data`` dictionary with weighted-averaged values for z, T, Td, u, v, and p **without metadata** (see example below)
   :rtype: dict


Example:

.. code-block:: python

    # load in some sounding data
    profile_a = spy.get_obs_data("OAX", "2014", "06", "16", "18")
    profile_b = spy.get_obs_data("OAX", "2014", "06", "17", "00")

    # apply `merge_profiles()` function
    merged_profile = spy.merge_profiles(profile_a, profile_b, weight_a=0.6)

    # reapply profile metadata and plot titles
    merged_profile['site_info'] = {'site-id': 'KOAX',
                              'site-name': 'OMAHA VALLEY',
                              'site-lctn': 'NE US',
                              'site-latlon': [41.32, -96.37],
                              'site-elv': 350,
                              'source': 'RAOB OBSERVED PROFILE',
                              'model': 'none',
                              'fcst-hour': 'none',
                              'run-time': 'none',
                              'valid-time': 'none',
                              'box_area': 'none'}
    merged_profile['titles'] =  {'top_title':   'COMPOSITE VERTICAL PROFILE',
                              'left_title':  'Profile-a (60%): 18z OAX 06/16/2014 | Profile-b (40%): 00z OAX 06/16/2014',
                              'right_title': 'KOAX, OMAHA VALLEY, NE | 41.32, -96.37    '}

    # build a sounding
    spy.build_sounding(merged_profile, radar=None, special_parcels='simple')

.. image:: /_static/images/merged_profile_example.png
   :width: 500 px
   :align: center

\ 

.. important::
   ``merge_profiles()`` returns the merged profile variables but does not
   automatically create new ``site_info`` or ``titles`` metadata. These should
   be supplied before producing a fully annotated SounderPy plot.

**************************************************




Smooth a profile
-----------------

Apply a simple ``SciPy`` gaussian smoother to a SounderPy ``clean_data`` profile

.. py:function:: spy.smooth_profile(clean_data, sigma)

   :param clean_data: a SounderPy ``clean_data`` dictionary.
   :type clean_data: dict, required
   :param sigma: standard deviation of the Gaussian kernel. Larger values produce stronger smoothing.
   :type sigma: float, optional, default is ``0.5``
   :return: a ``clean_data`` dictionary with smoothed values for z, T, Td, u, v, and p. Metadata remains unchanged.
   :rtype: dict


Here is an (extreme) example:

.. code-block:: python

    # load in some sounding data
    clean_data = spy.get_obs_data("OAX", "2014", "06", "16", "18")

    smoothed_profile = spy.smooth_profile(clean_data, 3)
    spy.build_sounding(smoothed_profile, radar=None, special_parcels='simple')

.. image:: /_static/images/smoothed_profile_example.png
   :width: 500 px
   :align: center

\ 
\ 
.. note::
   Smoothing alters the original observed or modeled profile. Use this tool
   carefully when calculating quantities that depend on sharp vertical
   gradients, such as lapse rates, inversions, shear, or parcel properties.

**************************************************



Finding a station latitude/longitude
--------------------------------------

.. py:function:: spy.get_latlon(station_type, station_id)


   Return the latitude and longitude associated with a supported station
   identifier.

   :param station_type: the station 'type' that corresponds with the given station ID
   :type station_type: str, required
   :param station_id: the station ID for the given station type
   :type station_id: str, required
   :return: lat/lon float pair
   :rtype: list

Example:

.. code-block:: python

   spy.get_latlon("metar", "KMOP")
   spy.get_latlon("bufkit", "APX")
   spy.get_latlon("raob", "OUN")
   spy.get_latlon("buoy", "45210")


.. tip::
   The returned ``[latitude, longitude]`` list can be passed directly to
   :ref:`get_model_data <modeldata>`.


***************************************************************



Interpolating a vertical profile
---------------------------------

.. py:function:: spy.interp_data(variable, heights, step=100)


   Interpolate a 1D array of data (such as a temperature profile) over a given interval (step) based on a corresponding array of height values. 

   :param variable: an array of data to be interpolated. Must be same length as height array.
   :type variable: array-like, required
   :param heights: heights corresponding to the vertical profile used to interpolate. Must be same length as variable array.
   :type heights: array-like, required
   :param step: the resolution of interpolation. Default is 100 (recommended value is 100)
   :type step: int, optional
   :return: interp_var, an array of interpolated data.
   :rtype: numpy.ndarray

Example:

.. code-block:: python

   spy.interp_data(temperature_array, height_array, step=100)  


***************************************************************




Finding a 'nearest' value
--------------------------

.. py:function:: spy.find_nearest(array, value)


   Return the index of the array element closest to a requested value.

   :param array: an array of data to be searched through
   :type array: arr, required
   :param value: Target value to locate within the array.
   :type value: int or float, required
   :return: nearest_idx, index of the data array that corresponds with the nearest value to the given value
   :rtype: int

Example:

.. code-block:: python

   idx_500m = spy.find_nearest(z, 500)

   print(z[idx_500m])

****************************************************************





Printing data to the console
-----------------------------

Print a compact summary of commonly used thermodynamic and kinematic parameters to the console.

.. py:function:: spy.print_variables(clean_data)

   :param clean_data: the dictionary of profile data to calculate profile parameters for (see :doc:`/reference/gettingdata`)
   :type clean_data: dict, required
   :return: ``None``. Calculated parameters, using the right-moving storm motion, are written to standard output.
   :rtype: None

.. code-block:: python

   spy.print_variables(clean_data)

.. code-block:: python

   > THERMODYNAMICS ---------------------------------------------
   --- SBCAPE: 2090.8 | MUCAPE: 2090.8 | MLCAPE: 1878.3 | MUECAPE: 1651.9
   --- MU 0-3: 71.1 | MU 0-6: 533.0 | SB 0-3: 71.1 | SB 0-6: 533.0

   > KINEMATICS -------------------------------------------------
   --- 0-500 SRW: 35.0 knot | 0-500 SWV: 0.019 | 0-500 SHEAR: 21.8 | 0-500 SRH: 186.2
   --- 1-3km SRW: 20.9 knot | 1-3km SWV: 0.005 | 1-3km SHEAR: 14.1 | | 1-3km SRH: 54.0

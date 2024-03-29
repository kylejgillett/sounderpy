📈 Tools for Plotting Data
===========================

.. image:: _static/images/example-sounding_dark3.png
   :width: 300 px
   :align: right


Initially, SounderPy was a tool meant for getting and parsing data with not as much focus on plotting said data. However, recent releases have put more emphasis on plotting capabilites. Version 3.0.0+ features a number of significant & exciting upgrades to SounderPy's plotting abilities. 

SounderPy can create general :ref:`sounding<soundings>` and :ref:`hodograph<hodographs>` plots as well as :ref:`composite sounding<compsoundings>` plots! 

The `full` sounding plots that SounderPy creates are complex figures with unique design geared towards severe convective storm enviroment analysis. 

These sounding plots are 'my baby' and I hope you find them useful! :)

Check out SounderPy's :doc:`examplegallery` & :doc:`examplescripts`

********************************************************************


.. tip::
   **Do your plots have funky scaling?** This is a common issue for smaller screen sizes. To fix this, use the `save=True` kwarg in your plot function.

.. _soundings:

Building Soundings
----------------------------------

We can use the simple ``spy.build_sounding()`` function:

.. py:function:: spy.build_sounding(clean_data, style='full', color_blind=False, dark_mode=False, storm_motion='right_moving', special_parcels=None, save=False, filename='sounderpy_sounding')


   Return a full sounding plot of SounderPy data, ``plt`` 

   :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
   :type clean_data: dict, required
   :param style: may be `simple` or `full`. Default is `full`.
   :type style: str, optional
   :param color_blind: whether or not to change the dewpoint trace line from green to blue for improved readability for color deficient users/readers. Default is ``False``
   :type color_blind: bool, optional
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
   :type dark_mode: bool, optional
   :param storm_motion: the storm motion used for plotting and calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional
   :param special_parcels: a nested list of special parcels from the ``ecape_parcels`` library. The nested list should be a list of two lists (`[[a, b], [c, d]]`) where the first list should include 'highlight parcels' and second list should include 'background parcels'. For more details, see the :ref:`parcels_logic` section. Another option is 'simple', which removes all advanced parcels making the plot quicker.
   :type special_parcels: nested `list` of two `lists`, optional
   :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
   :type save: bool, optional
   :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
   :type filename: str, optional
   :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
   :rtype: plt

Examples
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')

   # build the sounding! 
   spy.build_sounding(clean_data)

.. image:: _static/images/example-sounding_light.png
   :width: 800 px

**************************************************







.. _hodographs:

Building Hodographs
----------------------------------

Very similarly to soundings, we can use the simple ``spy.build_hodograph()`` function:

.. py:function:: spy.build_hodograph(clean_data, dark_mode=False, storm_motion='right_moving', sr_hodo=False, save=False, filename='sounderpy_sounding')


   Return a full hodograph plot of SounderPy data, ``plt`` 

   :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
   :type clean_data: dict, required
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
   :type dark_mode: bool, optional
   :param storm_motion: the storm motion used for plotting and calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional
   :param sr_hodo: transform the hodograph from ground relative to storm relative 
   :type sr_hodo: bool, optional, default is ``False``
   :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
   :type save: bool, optional
   :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
   :type filename: str, optional
   :return: plt, a SounderPy hodograph built with Matplotlib, MetPy, SharpPy, & SounderPy.
   :rtype: plt



Examples
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')

   # build the hodograph! 
   spy.build_hodograph(clean_data)

.. image:: _static/images/example-hodograph_light.png
   :width: 800 px

**************************************************




.. _compsoundings:

Building Composite Soundings
----------------------------------

Sometimes we want to compare two or more profiles against each other. Perhaps at different locations or times, or we may want to compare different models or model run-times. SounderPy allows you to do this!

To do so, a list of :ref:`'clean_data' dicts<datadescription>` is needed. If you want to customize the look of each profile, you can create equal length lists with alphas, linestyles, linewidths, & colors. See below:

.. py:function:: spy.build_composite(data_list, cmap='viridis', colors_to_use='none', shade_between=False, alphas_to_use='none', ls_to_use='none', lw_to_use='none', dark_mode=False, save=False, filename='sounderpy_sounding')


   Return a composite sounding plot of multiple profiles, ``plt`` 

   :param data_list: a list of data dictionaries for each profile to be plotted
   :type data_list: list of dicts, required
   :param shade_between: Lightly shade between the dewpoint & temperature trace. In many cases, this improves readability. Default is ``True``.
   :type shade_between: bool, optional
   :param cmap: a linear colormap, may be any custom or matplotlib cmap. Default is 'viridis'. If ``colors_to_use`` kwarg is provided, ``colors_to_use`` will be used instead.
   :type cmap: `matplotlib.colors.LinearSegmentedColormap` or `str` representing the name of a matplotlib cmap, optional
   :param colors_to_use: A list of custom matplotlib color name stings. List length must match the number of profiles listed in ``data_list``. Default is 'none'.
   :type colors_to_use: list of strings, optional
   :param alphas_to_use: A list of custom alphas (0.0-1.0). List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default alpha is 1.
   :type alphas_to_use: list of floats, optional
   :param ls_to_use: A list of custom matplotlib linestyles. List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default linestyle is '-'.
   :type ls_to_use: list of stings, optional
   :param lw_to_use: A list of custom linewidths. List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default linewidth is 3.
   :type lw_to_use: list of floats, optional
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
   :type dark_mode: bool, optional
   :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
   :type save: bool, optional
   :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
   :type filename: str, optional
   :return: plt, a SounderPy composite sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
   :rtype: plt


Examples
^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   # this example looks at 3 profiles from OAX on Pilger-day.
   clean_data1 = spy.get_obs_data('oax', '2014', '06', '16', '12')
   clean_data2 = spy.get_obs_data('oax', '2014', '06', '16', '18')
   clean_data3 = spy.get_obs_data('oax', '2014', '06', '17', '00')
     
   # add each dict of data to a list
   data_list = [clean_data1, clean_data2, clean_data3]

   # build the composite! 
   spy.build_composite(data_list)


.. image:: _static/images/example-composite_light.png
   :width: 800 px

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!

   data_list = []
   for hour in ['00', '01', '02', '03', '04', '05', '06']:
       cd = spy.get_bufkit_data('hrrr', 'dtx', 0, '2024', '02', '28', hour, hush=True)
       data_list.append(cd)

   # and make it dark-mode for fun! 
   spy.build_composite(data_list, dark_mode=True, lw_to_use=[4 for cd in data_list])

.. image:: _static/images/example-composite_dark.png
   :width: 800 px

**************************************************






.. _vadhodographs:

Building VAD Hodographs
----------------------------------

*Coming soon -- function is still in development*

SounderPy now offers the ability to plot NEXRAD radar VAD data on a hodograph using the ``spy.build_vad_hodograph()`` function:

.. py:function:: spy.build_vad_hodograph(vad_data, dark_mode=False, storm_motion='right_moving', sr_hodo=False, save=False, filename='sounderpy_sounding')

   Return a full hodograph plot of SounderPy VAD data, ``plt`` 

   :param vad_data: the dictionary of VAD data to be plotted
   :type vad_data: dict, required
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
   :type dark_mode: bool, optional
   :param storm_motion: the storm motion used for plotting and calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional
   :param sr_hodo: transform the hodograph from ground relative to storm relative 
   :type sr_hodo: bool, optional, default is ``False``
   :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
   :type save: bool, optional
   :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
   :type filename: str, optional
   :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
   :rtype: plt


Examples
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
   import datetime as dt
     
   # get VAD data using `pyart_radar_profile()`
   vad_data = spy.pyart_radar_profile('kgrr', dt.datetime(2023, 8, 25, 0, 20, 30))

   # build the hodograph! 
   spy.build_vad_hodograph(vad_data)


.. image:: _static/images/example_vad-hodograph.png
   :width: 800 px



**************************************************



.. _storm_motions:

Storm Motion Logic
-------------------

Users can define custom storm motions or choose from a number of 'storm motion keys' to change the storm motion considered by kinematic and thermodynamic parameters during calculations and plotting. All parameters that consider storm motion will be affected by the ``storm_motion`` kwarg. 

Storm Motion Keys 
^^^^^^^^^^^^^^^^^^
   - ``right_moving``: Bunkers Right Moving supercell (default)
   - ``left_moving``: Bunkers Left Moving supercell
   - ``mean_wind``: 0-6km mean wind.

   Example: 

   .. code-block:: python 

      storm_motion='left_moving'

Custom Storm Motions 
^^^^^^^^^^^^^^^^^^^^^
   Custom storm motions must be given in a `list` including direction in degrees and speed in knots. Note: degrees must be in the meteorological convention of 'from', i.e. 'northeast' would be 225 degrees, not 45 degrees.

   Example: 

   .. code-block:: python 

      # 250 degrees at 45 knots
      storm_motion=[250, 45]




**************************************************




.. _parcels_logic:

Parcel Logic
------------

New to v3.0.2+, the 'parcel-update', is a complex scheme for computing and plotting advanced parcels using various adiabatic ascent schemes and entrainment schemes. This toolkit comes from `Amelia Urquhart's <https://github.com/a-urq/ecape-parcel-py>`_ ``ecape-parcels`` Python package, which is based on work by `Peters et. al. 2022 <https://journals.ametsoc.org/view/journals/atsc/79/3/JAS-D-21-0118.1.xml>`_. 

When plotting soundings, users can choose from a number of parcel types to compute and plot, such as...

   - Pseudoadiabatic non-entraining ascent CAPE
   - Pseudoadiabatic entraining ascent CAPE
   - Irreversible Adiabatic non-entraining ascent CAPE
   - Irreversible Adiabatic entraining ascent CAPE


Each of these parcel types can be computed and plotted from a/the...

   - Surface-based parcel
   - Most-Unstable parcel
   - Mixed Layer parcel


How to use this feature 
^^^^^^^^^^^^^^^^^^^^^^^^

When plotting a `full` sounding using the ``build_sounding()`` function, use the kwarg `special_parcels` to choose which parcels you'd like to plot. This kwarg is a nested `list` (``[[a, b], [c, d]]``), where the first `list` contains 'highlight' parcels and the second `list` contains 'background' parcels. I.e., 'highlighted' parcels are darker and on top of 'background' parcels, which appear faded and behind the 'highlight' parcels. 

   - Example:

   .. code-block:: python
      
      special_parcels = [["sb_ia_ecape"], ["sb_ps_ecape", "sb_ps_cape"]]


   By default, SounderPy will plot normal MU/ML/SB-CAPE parcels and an mu_ia_ecape parcel. You can override this by setting ``special_parcels`` to 'simple', which only plots the common MU/ML/SB-CAPE parcels. This is greatly reduce the plot-time!


Parcel Keys 
^^^^^^^^^^^^^^^^^

Note the struture of the 'parcel key': ``sb_ia_ecape``. This is broken into three components: 'parcel-type', 'ascent-scheme', and 'entrainment-scheme'. You can make any parcel you like using this specific nomenclature: ``parcel-type_ascent-scheme_entrainment-scheme``.

  - PARCEL-TYPES
   - ``sb``: surface-based parcel
   - ``mu``: most-unstable parcel
   - ``ml``: mixed-layer parcel


  - ASCENT-SCHEMES
   - ``ps``: Pseudoadiabatic ascent
   - ``ia``: - Irreversible adiabatic ascent


  - ENTRAINMENT-SCHEMES
   - ``cape``: non-entraining convective available potential energy
   - ``ecape``: entraining convective available potential energy


  - Examples:
   - ``'sb_ia_ecape'``: surface-based irreversible adiabatic entraining CAPE
   - ``'mu_ps_cape'``: most-unstable pseudoadiabatic CAPE
   - ``'ml_ia_cape'``: mixed-layer irreversible adiabatic CAPE
   - ``'sb_ps_ecape'``: surface-based pseudoadiabatic entraining CAPE





**************************************************



Printing data to the console
-----------------------------

.. py:function:: spy.print_variables(clean_data, storm_motion='right_moving')

   :param clean_data: the dictionary of profile data to calculate profile parameters for (see :doc:`gettingdata`)
   :type clean_data: dict, required
   :param storm_motion: the storm motion used for calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional
   :return: prints a number of thermodynamic and kinematic variables to the console.
   :rtype: data print out to the console

.. code-block:: python

   > THERMODYNAMICS --------------------------------------------- 
   --- SBCAPE: 2090.8 | MUCAPE: 2090.8 | MLCAPE: 1878.3 | MUECAPE: 1651.9
   --- MU 0-3: 71.1 | MU 0-6: 533.0 | SB 0-3: 71.1 | SB 0-6: 533.0
    
   > KINEMATICS ------------------------------------------------- 
   --- 0-500 SRW: 35.0 knot | 0-500 SWV: 0.019 | 0-500 SHEAR: 21.8 | 0-500 SRH: 186.2
   --- 1-3km SRW: 20.9 knot | 1-3km SWV: 0.005 | 1-3km SHEAR: 14.1 | | 1-3km SRH: 54.0




**************************************************





About These Plots
-----------------

This plot style has been developed in a way that acts to provide as much information to the user as possible with attributes designed specifically for the analysis of severe convective environments, and supercells/tornadoes in particular. You will also find that this particular plot style does not host many of the common and popular severe weather composite indices – that was intentional. Most, if not all, of the data provided on this plot, are considered, for the lack of a better word, ‘true’ observations of the atmosphere though most are still subject to heavy assumptions. 

The data on these plots are considered, by most, to be useful in determining critical characteristics of the atmosphere related to supercellular storm mode and tornadogenesis.

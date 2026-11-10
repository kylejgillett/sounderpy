📈 Tools for Plotting Data
===========================

.. image:: /_static/images/example-sounding_dark3.png
   :width: 300 px
   :align: right


SounderPy can create general :ref:`sounding<soundings>` and :ref:`hodograph<hodographs>` plots as well as :ref:`composite sounding<compsoundings>` plots! 

The ``full`` sounding plot that SounderPy creates is a complex figure with unique design geared towards severe convective storm environment analysis. As of v3.3.0, a long-overdue winter version of the full sounding figure is available.

These sounding plots are 'my baby' and I hope you find them useful! :)

Check out SounderPy's :doc:`/gallery/examplegallery` & :ref:`example notebooks <example-notebooks>`

********************************************************************


.. tip::
   **Do your plots have funky scaling?** This is a common issue for smaller screen sizes. To fix this, use the ``save=True`` kwarg in your plot function.

.. _soundings:

Building Soundings
----------------------------------

While this function may appear complicated at first glance, there is actually only *one* required argument for this function (your sounding data: `clean_data`). The other arguments are simply optional settings for the plot!

.. py:function:: spy.build_sounding(clean_data, color_blind=False, dark_mode=False, storm_motion='right_moving', special_parcels=None, radar='mosaic', radar_time='sounding', map_zoom=2, modify_sfc=None, show_theta=False, hodo_boundary=None, dpi=100, save=False, filename='sounderpy_sounding', max_levels=250, ylim=None, theme=None, sr_hodo=False, mode='severe', show_obs=False)


   Return a full sounding plot of SounderPy data, ``plt`` 

   :param clean_data: the dictionary of data to be plotted (see :doc:`/reference/gettingdata`)
   :type clean_data: dict, required
   :param color_blind: whether or not to change the dewpoint trace line from green to blue for improved readability for color deficient users/readers.
   :type color_blind: bool, optional, Default is ``False``
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding.
   :type dark_mode: bool, optional, Default is ``False``.
   :param storm_motion: the storm motion used for plotting and calculations. Custom storm motions are accepted as a ``list`` of ``floats`` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional, Default is 'right_moving'.
   :param special_parcels: a nested list of special parcels from the ``ecape_parcels`` library. The nested list should be a list of two lists (``[[a, b], [c, d]]``) where the first list should include 'highlight parcels' and second list should include 'background parcels'. For more details, see the :ref:`parcels_logic` section. Another option is 'simple', which removes all advanced parcels making the plot quicker.
   :type special_parcels: nested ``list`` of two ``lists``, optional, Default is None
   :param radar: whether to display mosaic reflectivity, single-site reflectivity, or no radar in the map inset.
   :type radar: str, ``"mosaic"``, ``"single-site"``, or ``None``. Optional. Default is ``mosaic``.
   :param radar_time: radar mosaic data valid time. May be ``"sounding"`` (uses the valid time of the sounding data), ``"now"`` (current time/date), or a ``datetime`` obj **Note: radar mosaic data only goes back 1 month from current date**
   :type radar_time: str or ``datetime`` obj, optional, Default is ``"sounding"``.
   :param map_zoom: a 'zoom' level for the map inset as an ``int``. **Note: Setting** ``map_zoom=0`` **will hide the map**
   :type map_zoom: int, optional, Default is ``2``.
   :param modify_sfc: a ``dict`` in the format ``{'T': 25, 'Td': 21, 'ws': 20, 'wd': 270}`` to modify the surface values of the ``clean_data`` dict. See the :ref:`modify_sfc` section for more details.
   :type modify_sfc: None or dict, optional, default is None
   :param show_theta: whether to replace the piecewise CAPE plot with theta & theta-e. Default is ``False``
   :type show_theta: bool, optional
   :param hodo_boundary: plot a "boundary" (straight line) on the hodograph axis. Do so by providing an angle and a color. The angle represents the angle between zero-degrees (north) and the upper-half of the boundary line, to the right. Multiple boundaries can be plotted.
   :type hodo_boundary: dict of lists, structured as ``{'angle':[], 'color':[]}``, optional, Default is ``None``
   :param dpi: Dots per inch resolution of the figure.
   :type dpi: int, optional, Default is 100.
   :param save: whether to show the plot inline or save to a file.
   :type save: bool, optional, Default is ``False``.
   :param filename: the filename by which a file should be saved to if ``save = True``.
   :type filename: str, optional, Default is ``sounderpy_sounding``.
   :param max_levels: the most vertical levels used to calculate and plot the sounding. Profiles with more levels are thinned to this many evenly spaced levels. Set to ``None`` to use every level.
   :type max_levels: int or None, optional, Default is ``250``.
   :param ylim: the skew-t's pressure range as (bottom, top) in hPa, such as ``(1000, 700)``, to look at one part of the atmosphere. The temperature axis and pressure ticks fit to that range, while the rest of the figure still uses the full profile.
   :type ylim: tuple of two numbers, optional, Default is ``None``, which shows the full profile.
   :param theme: changes to the figure's colors, such as ``{'border': 'k', 'temperature': 'navy', 'parcel': {'mu': 'black'}}``. See :ref:`plot_themes` for every setting.
   :type theme: dict, optional, Default is ``None``.
   :param sr_hodo: draw the hodograph storm relative, shifting the winds by the storm motion from ``storm_motion``. If the profile has no storm motion, the hodograph stays ground relative.
   :type sr_hodo: bool, optional, Default is ``False``.
   :param mode: which sounding to draw. ``"severe"`` is the standard severe-weather sounding. ``"winter"`` draws the winter sounding, built around precipitation type, snow growth, and winter hazards (see :doc:`/tutorials/winter_soundings`). ``"auto"`` draws the winter sounding when the surface wet-bulb temperature is below 35 degrees F, and the severe sounding otherwise. The winter sounding doesn't use ``special_parcels``, ``show_theta``, ``sr_hodo``, or ``hodo_boundary``; if any are set, it warns and ignores them.
   :type mode: str, ``"severe"``, ``"winter"``, or ``"auto"``, optional, Default is ``"severe"``.
   :param show_obs: winter sounding only: draw a station model of the latest surface observation (METAR) near the profile, valid within an hour of the profile's time, from the Iowa Environmental Mesonet. Needs internet access. Values missing from the report read "NA".
   :type show_obs: bool, optional, Default is ``False``.
   :return: plt, a SounderPy sounding figure.
   :rtype: plt

Examples
^^^^^^^^^^

.. container:: spy-compare

   .. container:: spy-compare-col

      **Severe sounding**

      .. image:: /_static/images/example-sounding_light.png

      .. code-block:: python

         import sounderpy as spy

         # get data | Note: any sounderpy data will work!
         clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')

         # build the sounding!
         spy.build_sounding(clean_data)

   .. container:: spy-compare-col

      **Winter sounding**

      .. image:: /_static/images/example-winter-sounding_light.png

      .. code-block:: python

         import sounderpy as spy

         # get data | Note: any sounderpy data will work!
         clean_data = spy.get_bufkit_data('hrrr', 'kgfk', 1,'2024', '01', '24', '12')

         # build the winter sounding, with station plot and radar
         spy.build_sounding(clean_data, mode='winter', show_obs=True, radar='single-site')

See :doc:`/tutorials/winter_soundings` for a guide to the winter sounding.

**************************************************







.. _hodographs:

Building Hodographs
----------------------------------

Very similarly to soundings, we can use the simple ``spy.build_hodograph()`` function:

.. py:function:: spy.build_hodograph(clean_data, dark_mode=False, storm_motion='right_moving', sr_hodo=False, modify_sfc=None, radar='mosaic', radar_time='sounding', map_zoom=2, hodo_boundary=None, dpi=100, save=False, filename='sounderpy_hodograph', max_levels=250, theme=None)


   Return a full hodograph plot of SounderPy data, ``plt`` 

   :param clean_data: the dictionary of data to be plotted (see :doc:`/reference/gettingdata`)
   :type clean_data: dict, required
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding.
   :type dark_mode: bool, optional, Default is ``False``.
   :param storm_motion: the storm motion used for plotting and calculations. Custom storm motions are accepted as a ``list`` of ``floats`` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional, Default is 'right_moving'.
   :param sr_hodo: transform the hodograph from ground relative to storm relative 
   :type sr_hodo: bool, optional, default is ``False``
   :param modify_sfc: a ``dict`` in the format ``{'T': 25, 'Td': 21, 'ws': 20, 'wd': 270}`` to modify the surface values of the ``clean_data`` dict.
   :type modify_sfc: None or dict, optional, default is None
   :param radar: whether to display mosaic reflectivity, single-site reflectivity, or no radar in the map inset.
   :type radar: str, ``"mosaic"``, ``"single-site"``, or ``None``. Optional. Default is ``mosaic``.
   :param radar_time: radar mosaic data valid time. May be ``"sounding"`` (uses the valid time of the sounding data), ``"now"`` (current time/date), or a ``datetime`` obj **Note: radar mosaic data only goes back 1 month from current date**
   :type radar_time: str or ``datetime`` obj, optional, Default is ``"sounding"``.
   :param map_zoom: a 'zoom' level for the map inset as an ``int``. **Note: Setting** ``map_zoom=0`` **will hide the map**
   :type map_zoom: int, optional, Default is ``2``.
   :param hodo_boundary: plot a "boundary" (straight line) on the hodograph axis. Do so by providing an angle and a color. The angle represents the angle between zero-degrees (north) and the upper-half of the boundary line, to the right. Multiple boundaries can be plotted.
   :type hodo_boundary: dict of lists, structured as ``{'angle':[], 'color':[]}``, optional, Default is ``None``
   :param dpi: Dots per inch resolution of the figure.
   :type dpi: int, optional, Default is 100.
   :param save: whether to show the plot inline or save to a file.
   :type save: bool, optional, Default is ``False``.
   :param filename: the filename by which a file should be saved to if ``save = True``.
   :type filename: str, optional, Default is ``sounderpy_hodograph``.
   :param max_levels: the most vertical levels used to calculate and plot the hodograph. Profiles with more levels are thinned to this many evenly spaced levels. Set to ``None`` to use every level.
   :type max_levels: int or None, optional, Default is ``250``.
   :param theme: changes to the figure's colors, such as ``{'border': 'k', 'temperature': 'navy', 'parcel': {'mu': 'black'}}``. See :ref:`plot_themes` for every setting.
   :type theme: dict, optional, Default is ``None``.
   :return: plt, a SounderPy hodograph figure
   :rtype: plt


Example
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')

   # build the hodograph! 
   spy.build_hodograph(clean_data, radar='single')

.. image:: /_static/images/example-hodograph_light.png
   :width: 800 px

**************************************************








.. _compsoundings:

Building Composite Soundings
----------------------------------

Sometimes we want to compare two or more profiles against each other. Perhaps at different locations or times, or we may want to compare different models or model run-times. SounderPy allows you to do this!

To do so, a list of :ref:`'clean_data' dicts<datadescription>` is needed. If you want to customize the look of each profile, you can create equal length lists with alphas, linestyles, linewidths, & colors. See below:

.. py:function:: spy.build_composite(data_list, shade_between=True, cmap='viridis', colors=None, linestyles=None, alphas=None, linewidths=None, dark_mode=False, save=False, filename='sounderpy_composite', show_params=True, storm_motion='right_moving', max_levels=250, dpi=100, ylim=None, theme=None, style='full', sr_hodo=False)


   Return a composite sounding plot of multiple profiles, ``plt`` 

   :param data_list: a list of ``clean_data`` dictionaries for each profile to be plotted. Profiles may come from any mix of sources, such as observed, model, BUFKIT, ACARS, and custom data.
   :type data_list: list of dicts, required
   :param shade_between: Lightly shade between the dewpoint & temperature trace. In many cases, this improves readability.
   :type shade_between: bool, optional, Default is ``True``.
   :param cmap: a linear colormap, may be any custom or matplotlib cmap. The default uses only the part of 'viridis' that stands out from the figure's background. If ``colors`` kwarg is provided, ``colors`` will be used instead.
   :type cmap: ``matplotlib.colors.LinearSegmentedColormap`` or ``str`` representing the name of a matplotlib cmap, optional, Default is 'viridis'.
   :param colors: A list of custom matplotlib color name strings. List length must match the number of profiles listed in ``data_list``.
   :type colors: list of strings, optional, Default is ``None``.
   :param alphas: A list of custom alphas (0.0-1.0). List length must match the number of profiles listed in ``data_list``.
   :type alphas: list of floats, optional, Default is ``None`` (sets alpha to 1)
   :param linestyles: A list of custom matplotlib linestyles. List length must match the number of profiles listed in ``data_list``.
   :type linestyles: list of strings, optional, Default is ``None`` (sets linestyle to '-')
   :param linewidths: A list of custom linewidths. List length must match the number of profiles listed in ``data_list``.
   :type linewidths: list of floats, optional, Default is ``None`` (sets linewidth to 3).
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding.
   :type dark_mode: bool, optional, Default is ``False``.
   :param save: whether to show the plot inline or save to a file.
   :type save: bool, optional, Default is ``False``
   :param filename: the filename by which a file should be saved to if ``save = True``.
   :type filename: str, optional, Default is ``sounderpy_composite``.
   :param show_params: whether to compute each profile's parameters and show them in the table. Setting this to ``False`` is faster.
   :type show_params: bool, optional, Default is ``True``.
   :param storm_motion: the storm motion used to compute every profile's parameters. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional, Default is 'right_moving'.
   :param max_levels: the most vertical levels used to calculate and plot each profile. Profiles with more levels are thinned to this many evenly spaced levels. Set to ``None`` to use every level.
   :type max_levels: int or None, optional, Default is ``250``.
   :param dpi: Dots per inch resolution of the figure.
   :type dpi: int, optional, Default is 100.
   :param ylim: the skew-t's pressure range as (bottom, top) in hPa, such as ``(1000, 700)``, to look at one part of the atmosphere. The temperature axis and pressure ticks fit to that range, while the rest of the figure still uses the full profile.
   :type ylim: tuple of two numbers, optional, Default is ``None``, which shows the full profile.
   :param theme: changes to the figure's colors, such as ``{'border': 'k', 'temperature': 'navy', 'parcel': {'mu': 'black'}}``. See :ref:`plot_themes` for every setting.
   :type theme: dict, optional, Default is ``None``.
   :param style: the panels to show. ``'full'`` shows the skew-T and hodograph, ``'sounding'`` shows only the skew-T, and ``'hodograph'`` shows only the hodograph. Each style's table lists the parameters that go with its panels. ``ylim`` is not used with ``'hodograph'``.
   :type style: str, optional, Default is ``'full'``.
   :param sr_hodo: draw the hodograph storm relative, shifting each profile's winds by its own storm motion from ``storm_motion``. Not used with ``style='sounding'``. If any profile has no storm motion, the hodograph stays ground relative.
   :type sr_hodo: bool, optional, Default is ``False``.
   :return: plt, a SounderPy composite sounding figure
   :rtype: plt

.. note::
   ``colors``, ``linestyles``, ``alphas``, and ``linewidths`` replace the older ``colors_to_use``, ``ls_to_use``, ``alphas_to_use``, and ``lw_to_use``. The old names still work for now, but show a ``FutureWarning`` and will be removed in a future release.


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


.. image:: /_static/images/example-composite_light.png
   :width: 800 px

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!

   data_list = []
   for hour in ['00', '01', '02', '03', '04', '05', '06']:
       cd = spy.get_bufkit_data('hrrr', 'dtx', 0, '2024', '02', '28', hour, hush=True)
       data_list.append(cd)

   # and make it dark-mode for fun! 
   spy.build_composite(data_list, dark_mode=True, linewidths=[4 for cd in data_list])

.. image:: /_static/images/example-composite_dark.png
   :width: 800 px

**************************************************






.. _vadhodographs:

Building VAD Hodographs
----------------------------------

.. caution::
   VAD hodograph support is experimental and may be more sensitive to missing,
   sparse, or unusual radar-derived wind profiles than the standard sounding
   hodograph workflow.

SounderPy now offers the ability to plot NEXRAD radar VAD data on a hodograph using the ``spy.build_vad_hodograph()`` function:

.. py:function:: spy.build_vad_hodograph(vad_data, save=False, dark_mode=False, storm_motion='right_moving', sr_hodo=False, filename='vad_hodograph', radar='mosaic', radar_time='sounding', map_zoom=2, hodo_boundary=None, dpi=100, theme=None)

   Return a full hodograph plot of SounderPy VAD data, ``plt`` 

   :param vad_data: the dictionary of VAD data to be plotted
   :type vad_data: dict, required
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding.
   :type dark_mode: bool, optional,  Default is ``False``.
   :param storm_motion: the storm motion used for plotting and calculations. Custom storm motions are accepted as a ``list`` of ``floats`` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
   :type storm_motion: str or list of floats, optional, Default is 'right_moving'.
   :param sr_hodo: transform the hodograph from ground relative to storm relative 
   :type sr_hodo: bool, optional, default is ``False``
   :param save: whether to show the plot inline or save to a file.
   :type save: bool, optional,  Default is ``False``
   :param filename: the filename by which a file should be saved to if ``save = True``.
   :type filename: str, optional, Default is ``vad_hodograph``.
   :param radar: whether to display mosaic reflectivity, single-site reflectivity, or no radar in the map inset.
   :type radar: str, ``"mosaic"``, ``"single-site"``, or ``None``. Optional. Default is ``mosaic``.
   :param radar_time: the time of radar to show: ``"sounding"`` for the data's valid time, or ``"now"``.
   :type radar_time: str, optional, Default is ``"sounding"``.
   :param map_zoom: a 'zoom' level for the map inset as an ``int``. **Note: Setting** ``map_zoom=0`` **will hide the map**
   :type map_zoom: int, optional, Default is ``2``.
   :param hodo_boundary: plot a "boundary" (straight line) on the hodograph axis, structured as ``{'angle':[], 'color':[]}``.
   :type hodo_boundary: dict of lists, optional, Default is ``None``
   :param dpi: Dots per inch resolution of the figure.
   :type dpi: int, optional, Default is 100.
   :param theme: changes to the figure's colors, such as ``{'border': 'k', 'temperature': 'navy', 'parcel': {'mu': 'black'}}``. See :ref:`plot_themes` for every setting.
   :type theme: dict, optional, Default is ``None``.
   :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
   :rtype: plt


Examples
^^^^^^^^^^

.. image:: /_static/images/example_vad-hodograph.png
   :width: 800 px


.. _plot_themes:

Figure Colors
----------------------------------

Every plotting function accepts a ``theme`` dictionary that changes individual colors in the figure. Only the settings you include change; everything else keeps the light or dark color scheme. ``color_blind=True`` and ``dark_mode=True`` are applied first, so a ``theme`` can adjust them further.

.. code-block:: python

   spy.build_sounding(clean_data, theme={'border': 'k',
                                         'temperature': 'navy',
                                         'parcel': {'mu': 'black'}})

Parcel colors can be changed one at a time, as with ``'parcel': {'mu': 'black'}`` above. An invalid setting name raises an error that lists every valid name.

Theme settings
^^^^^^^^^^^^^^

- **Figure**: ``text``, ``background``, ``border``, ``marker``, ``barb``, ``skew_line``, ``shade_alpha``
- **Skew-T profiles**: ``temperature``, ``dewpoint``, ``virtual``, ``wetbulb``
- **Skew-T background**: ``adiabat``, ``isotherm``, ``shade``
- **Parcels**: ``parcel`` (a dict with ``sb``, ``ml``, and ``mu``), ``special_parcel``, ``secondary_parcel``, ``downdraft``
- **Hodograph**: ``hodo_layers``, ``height_marker``, ``height_marker_edge``, ``mcs``, ``dtm``
- **Skew-T annotations**: ``pbl``, ``eil``, ``dgz``, ``hgz``, ``lapse_max``
- **Shared warm and cool values**: ``warm``, ``cool``, ``ascent``
- **Tables**: ``moisture``, ``moisture_strong``, ``kinematics``, ``kinematics_strong``, ``kinematics_layers``, ``ncape``, ``lapse``, ``composite_index``, ``ehi``
- **Accessory panels**: ``bar_edge``, ``pressure_line``, ``low_fill``, ``omega_zero``, ``omega_reference``, ``zero_line``, ``streamwise``, ``theta``, ``u_wind``, ``v_wind``, ``cin``, ``cape``
- **Watermark**: ``watermark``

``hodo_layers`` is a sequence of ``(bottom_m, top_m, color)`` layers, such as ``((0, 1000, 'purple'), (1000, 3000, 'red'))``, and ``kinematics_layers`` is a sequence of five colors.


.. _watermark:

Watermark
----------------------------------

SounderPy figures by default show the SounderPy logo and credit text. In a new addition, ``spy.set_watermark()`` allows you to change it for every figure made afterwards in the same Python session. This is especially useful for field campaign projects or university observations.

.. py:function:: spy.set_watermark(**settings)

   Change the watermark on SounderPy figures, for this Python session. Pass only the settings to change. Calling this with no arguments restores the default watermark.

   :param show: whether to show the watermark text and logo
   :type show: bool, optional, Default is ``True``.
   :param title: the top line of text, on every figure. By default, each figure uses its own text
   :type title: str, optional
   :param subtitle: the bottom line of text, on every figure. By default, each figure uses its own text
   :type subtitle: str, optional
   :param logo: a path to an image to use as the logo, or ``None`` for no logo. By default, the SounderPy logo is used
   :type logo: str or None, optional
   :param color: the color of the text. By default, the figure's ``watermark`` theme color is used
   :type color: str, optional
   :param alpha: the transparency of the text
   :type alpha: float, optional, Default is ``0.6``.

Example
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy

   # use your own text and logo on every figure
   spy.set_watermark(title='MY RESEARCH GROUP', subtitle='Made with SounderPy', logo='my_logo.png')

   # hide the watermark
   spy.set_watermark(show=False)

   # restore the default watermark
   spy.set_watermark()

.. note::
   Watermark settings are shared by everything running in the same Python process.


**************************************************








.. _modify_sfc:

Surface Modification Logic
---------------------------

Users may override the surface values of the ``clean_data`` dict, when creating a plot or returning sounding parameters, using the ``modify_sfc`` kwarg. This argument is a python dictionary. New to v3.3.0, the surface modification logic now includes three different "methods" of modification: `"mix"`, `"blend"`, and `"replace"`.

With the default ``'method': 'mix'``, a surface temperature warmer than the profile is mixed upward along its dry adiabat until it meets the original profile, the depth a heated surface mixes to (Holzworth 1964). Through that layer, a new surface dewpoint is held at a constant mixing ratio, as in a well-mixed boundary layer. A cooler surface temperature, a dewpoint-only change, and wind changes will change only the surface level (i.e., no mixing). Dewpoint is never allowed to exceed temperature.

With ``'method': 'blend'``, the temperature and dewpoint changes taper linearly from the surface to zero at ``depth`` meters above ground (500 m by default), so a warmer surface may give a superadiabatic layer of that depth while the profile above it is unchanged. For example, a surface 4°C warmer blended over 250 m adds 16 °C/km to that layer's lapse rate. This option is usually unrealistic in the warm season, but may be useful in the cool season (e.g., blowing snow, lake effect, etc)

With ``'method': 'replace'``, only the surface level changes, so a superadiabatic surface layer is kept. This suits cases where the lowest few hundred meters matter most.


The ``modify_sfc`` dict
^^^^^^^^^^^^^^^^^^^^^^^^
- ``T``: Temperature, degrees Celsius
- ``Td``: Dewpoint, degrees Celsius
- ``ws``: Wind speed, knots
- ``wd``: Wind direction, meteorological degrees (north=0)
- ``method``: ``'mix'`` (default), ``'replace'``, or ``'blend'``
- ``depth``: Depth of the blend, meters above ground (``'blend'`` only, default 500 (meters))

Examples:

.. code-block:: python

   spy.build_sounding(clean_data, modify_sfc={'T':21, 'Td':19, 'ws': 30, 'wd':270})

   spy.build_sounding(clean_data, modify_sfc={'T':21, 'Td':19})

   spy.build_sounding(clean_data, modify_sfc={'ws': 30, 'wd':270})

   spy.build_sounding(clean_data, modify_sfc={'T':-8, 'Td':-11, 'method':'replace'})

   spy.build_sounding(clean_data, modify_sfc={'T':-8, 'Td':-11, 'method':'blend', 'depth':250})

.. caution::
   Surface modification should be used carefully. Replacing the surface
   thermodynamic or wind state can materially alter parcel and kinematic
   calculations.

**************************************************








.. _storm_motions:

Storm Motion Logic
-------------------

Users can define custom storm motions or choose from a number of 'storm motion keys' to change the storm motion considered by kinematic and thermodynamic parameters during calculations and plotting. All parameters that consider storm motion will be affected by the ``storm_motion`` kwarg. 

Storm Motion Keys 
^^^^^^^^^^^^^^^^^^
- ``"right_moving"`` — Bunkers right-moving supercell motion (default). Also ``"rm"``, ``"right"``, or ``"right_mover"``
- ``"left_moving"`` — Bunkers left-moving supercell motion. Also ``"lm"``, ``"left"``, or ``"left_mover"``
- ``"mean_wind"`` — 0-6km mean wind. Also ``"mw"`` or ``"mean"``

Keys are not case sensitive. Any other value raises a ``ValueError``.

Example: 

.. code-block:: python 

   storm_motion='left_moving'


Custom Storm Motions 
^^^^^^^^^^^^^^^^^^^^^
Custom storm motions are supplied as:

.. code-block:: python

   [direction_degrees, speed_knots]

Direction follows the meteorological ``from`` convention. For example, 250 degrees at 45kts:

.. code-block:: python

   storm_motion = [250, 45]




**************************************************




.. _parcels_logic:

Parcel Logic
------------

SounderPy supports several parcel-ascent configurations, including pseudoadiabatic and irreversible-adiabatic ascent with entraining and non-entraining formulations. This toolkit comes from `Amelia Urquhart's <https://github.com/a-urq/ecape-parcel-py>`_ ``ecape-parcels`` Python package, which is based on work by `Peters et al. 2022 <https://journals.ametsoc.org/view/journals/atsc/79/3/JAS-D-21-0118.1.xml>`_. 

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

When plotting a ``full`` sounding using the ``build_sounding()`` function, use the kwarg ``special_parcels`` to choose which parcels you'd like to plot. This kwarg is a nested ``list`` (``[[a, b], [c, d]]``), where the first ``list`` contains 'highlight' parcels and the second ``list`` contains 'background' parcels. I.e., 'highlighted' parcels are darker and on top of 'background' parcels, which appear faded and behind the 'highlight' parcels. 

- Example:

.. code-block:: python
   
   special_parcels = [["sb_ia_ecape"], ["sb_ps_ecape", "sb_ps_cape"]]


By default, SounderPy will plot normal MU/ML/SB-CAPE parcels and an mu_ia_ecape parcel. You can override this by setting ``special_parcels`` to 'simple', which only plots the common MU/ML/SB-CAPE parcels. This can greatly reduce the plot-time!


Parcel Keys 
^^^^^^^^^^^^^^^^^

Note the structure of the 'parcel key': ``sb_ia_ecape``. This is broken into three components: 'parcel-type', 'ascent-scheme', and 'entrainment-scheme'. You can make any parcel you like using this specific nomenclature: ``parcel-type_ascent-scheme_entrainment-scheme``.

- PARCEL-TYPES

  - ``sb``: surface-based parcel
  - ``mu``: most-unstable parcel
  - ``ml``: mixed-layer parcel


- ASCENT-SCHEMES

  - ``ps``: Pseudoadiabatic ascent
  - ``ia``: Irreversible adiabatic ascent


- ENTRAINMENT-SCHEMES

  - ``cape``: non-entraining convective available potential energy
  - ``ecape``: entraining convective available potential energy


- Examples:

  - ``'sb_ia_ecape'``: surface-based irreversible adiabatic entraining CAPE
  - ``'mu_ps_cape'``: most-unstable pseudoadiabatic CAPE
  - ``'ml_ia_cape'``: mixed-layer irreversible adiabatic CAPE
  - ``'sb_ps_ecape'``: surface-based pseudoadiabatic entraining CAPE


Setting ``special_parcels="simple"`` plots only the standard SB/ML/MU CAPE parcels and skips the additional MU ECAPE parcel, reducing calculation and plotting time.


**************************************************


.. _radar_logic:

Radar Logic
------------
Knowing the location of a vertical profile relative to nearby meteorological phenomenon is vital in understanding what that vertical profile means for the atmosphere at a given time. To help with this, SounderPy employs the use of mapping & radar data to provide this information to a reader.

Using the Map
^^^^^^^^^^^^^

The map inset is controlled with ``map_zoom``.

- ``map_zoom=0`` disables the map.
- Larger positive values display a broader geographic area.


Using Radar
^^^^^^^^^^^^^^^^^
SounderPy offers two sources of radar data to the user. First is a CONUS reflectivity mosaic, and the other is single-site NEXRAD WSR-88D reflectivity. The CONUS mosaic is default because its fastest and easiest to plot, but data availability is only up to 30 days before the current date. Single-site radar goes back much further but data availability varies by site and date.

To plot radar on a sounding or hodograph, there are two kwargs to understand:

- ``radar``

  - ``"mosaic"``: plots the CONUS reflectivity mosaic
  - ``"single"``: plots single-site WSR-88D reflectivity closest to the latitude/longitude coordinates of the sounding.

- ``radar_time``

  - ``"now"``: plots the most recent radar scan
  - ``"sounding"``: plots the radar scan closest the time of the sounding
  - ``datetime(1999, 5, 4, 0, 27)``: or you may define your own time by providing a datetime obj. This will plot the nearest radar scan to the date/time provided.

**************************************************

About These Plots
-----------------

SounderPy's full sounding layout is designed to summarize thermodynamic and
kinematic characteristics commonly examined in severe-convective-environment
analysis.

The figure intentionally emphasizes physically interpretable profile
quantities and parcel/kinematic diagnostics rather than displaying a large
number of composite severe-weather indices.

The interpretation of any SounderPy figure depends on the source profile.
Observed soundings, aircraft profiles, model forecasts, analyses, and
reanalyses each have different sampling characteristics and limitations.
Users should consider those limitations when interpreting calculated
parameters and environmental diagnostics.

The data on these plots are considered, by most, to be useful in determining critical characteristics of the atmosphere related to supercellular storm mode and tornadogenesis.

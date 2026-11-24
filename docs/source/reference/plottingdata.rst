📈 Tools for Plotting Data
===========================

SounderPy draws three kinds of figures from any profile it reads: its own **operational** figures, **custom**
figures you build on top of them, and clean **publication** figures for journals. The examples use the HRRR's 3-hour
forecast for Aberdeen, SD, valid 00 UTC 29 August 2024, from the :doc:`Mound City case study </case_studies/mound_city_20240828>`.

.. container:: home-card-grid

   .. container:: home-card

      .. image:: /_static/gallery/custom_figures/operational_sounding.png

      **📈 Operational Figures**

      SounderPy's full sounding, hodograph, composite, and VAD hodograph figures, built for forecasting and
      analysis, in one line.

      :ref:`Reference <soundings>` · :doc:`Tutorials </tutorials/plot_appearance>`

   .. container:: home-card

      .. image:: /_static/gallery/custom_figures/return_fig_sounding.png

      **✏️ Custom Figures**

      SounderPy's figures with your own touches: return the figure and its panels to annotate, or draw your own
      data on the map.

      :ref:`Reference <returned_figures>` · :ref:`Tutorial <custom-figures-drawing>`

   .. container:: home-card

      .. image:: /_static/gallery/custom_figures/publication_default.png

      **📝 Publication Figures**

      A clean, journal-sized skew-t, hodograph, map, and parameter table of one or more profiles (experimental).

      :ref:`Reference <publicationfigures>` · :ref:`Tutorial <custom-figures-publication>`

.. code-block:: python

   import sounderpy as spy

   clean_data = spy.get_bufkit_data('hrrr', 'kabr', 3, '2024', '08', '28', '21')

   spy.build_sounding(clean_data)                                # operational
   fig, axes = spy.build_sounding(clean_data, return_fig=True)   # custom: draw on axes, then fig.savefig(...)
   fig, axes = spy.build_publication_sounding(clean_data)        # publication

The ``full`` sounding plot that SounderPy creates is a complex figure with unique design geared towards severe
convective storm environment analysis. As of v3.3.0, a long-overdue winter version of the full sounding figure is
available. These sounding plots are 'my baby' and I hope you find them useful! :)

Check out SounderPy's :doc:`/gallery/examplegallery` & :ref:`example notebooks <example-notebooks>`

The settings these figures share (colors and themes, the watermark, surface modification, storm motions, parcels,
and the map and radar) are on :doc:`/reference/plotsettings`, and every option of the publication figure is on
:doc:`/reference/publicationfigures`.

********************************************************************


.. tip::
   **Do your plots have funky scaling?** This is a common issue for smaller screen sizes. To fix this, use the ``save=True`` kwarg in your plot function.

.. _soundings:

Building Soundings
----------------------------------

While this function may appear complicated at first glance, there is actually only *one* required argument for this function (your sounding data: `clean_data`). The other arguments are simply optional settings for the plot!

.. py:function:: spy.build_sounding(clean_data, color_blind=False, dark_mode=False, storm_motion='right_moving', special_parcels=None, radar='mosaic', radar_time='sounding', map_zoom=2, modify_sfc=None, show_theta=False, hodo_boundary=None, dpi=100, save=False, filename='sounderpy_sounding', max_levels=250, ylim=None, theme=None, sr_hodo=False, mode='severe', show_obs=False, map_layers=None, map_label=None, return_fig=False, readout=False)


   Build a full sounding figure of SounderPy data

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
   :type radar: str, ``"mosaic"``, ``"single"`` (or ``"single-site"``), or ``None``. Optional. Default is ``mosaic``.
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
   :param map_layers: your own data for the map inset: a function, or a list of functions, that each take the map's axes (a cartopy ``PlateCarree`` axes) and draw on it, such as ``lambda ax: ax.contourf(lon, lat, values, transform=ccrs.PlateCarree(), zorder=3)``. They're drawn after the radar and before the profile's location marker; ``zorder=3`` puts filled data above the satellite image and below the radar and the county and state lines. Use ``radar=None`` to show only your data. See :ref:`custom_map_data`.
   :type map_layers: function or list of functions, optional, Default is ``None``.
   :param map_label: a label for the map, such as ``"HRRR MLCAPE (J/kg)"``, shown in a box at the top of the map.
   :type map_label: str, optional, Default is ``None``.
   :param return_fig: return the figure instead of showing it, so you can draw on it or save it yourself. A figure that's saved and not returned is closed.
   :type return_fig: bool, optional, Default is ``False``.
   :param readout: EXPERIMENTAL: with ``save=True``, also write a readout file next to the figure (``kokc.png`` -> ``kokc.readout.json``) for hover readouts and profile editing on a web page. See :doc:`/developers/interactive_figures`.
   :type readout: bool, optional, Default is ``False``.
   :return: ``None``, or with ``return_fig=True``, ``(fig, axes)``: the Matplotlib figure, and a dict of its panels: ``'skewt'`` (a MetPy ``SkewT``), ``'hodograph'`` (a MetPy ``Hodograph``), and ``'map'`` (a cartopy ``GeoAxes``, when the map is drawn).
   :rtype: None or tuple

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
         spy.build_sounding(clean_data, color_blind=True, radar='single', map_zoom=1)

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

.. _drawing_on_figures:

**Drawing on the figure**

``return_fig=True`` hands back the figure and its panels instead of showing it, so you can add your own
annotations and save it in any format. This works the same way for ``build_hodograph()``, ``build_composite()``,
and ``build_vad_hodograph()``. See :ref:`returned_figures` for each function's panels.

.. code-block:: python

   import sounderpy as spy
   import matplotlib.pyplot as plt

   clean_data = spy.get_bufkit_data('hrrr', 'kabr', 3, '2024', '08', '28', '21')

   fig, axes = spy.build_sounding(clean_data, radar='single', map_zoom=1, return_fig=True)

   # the skew-t is a MetPy SkewT: draw on its axes in pressure (hPa) and temperature (°C)
   axes['skewt'].ax.axhspan(850, 700, color='gold', alpha=0.2, zorder=0)
   axes['skewt'].ax.text(-28, 775, 'EML', fontsize=16, weight='bold', va='center')

   fig.savefig('kabr_sounding.png', bbox_inches='tight')
   plt.close(fig)

.. image:: /_static/gallery/custom_figures/return_fig_sounding.png
   :width: 100%

**************************************************







.. _hodographs:

Building Hodographs
----------------------------------

Very similarly to soundings, we can use the simple ``spy.build_hodograph()`` function:

.. py:function:: spy.build_hodograph(clean_data, dark_mode=False, storm_motion='right_moving', sr_hodo=False, modify_sfc=None, radar='mosaic', radar_time='sounding', map_zoom=2, hodo_boundary=None, dpi=100, save=False, filename='sounderpy_hodograph', max_levels=250, theme=None, map_layers=None, map_label=None, return_fig=False, readout=False)


   Build a hodograph figure of SounderPy data

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
   :type radar: str, ``"mosaic"``, ``"single"`` (or ``"single-site"``), or ``None``. Optional. Default is ``mosaic``.
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
   :param map_layers: your own data for the map inset: a function, or a list of functions, that each take the map's axes (a cartopy ``PlateCarree`` axes) and draw on it, such as ``lambda ax: ax.contourf(lon, lat, values, transform=ccrs.PlateCarree(), zorder=3)``. They're drawn after the radar and before the profile's location marker; ``zorder=3`` puts filled data above the satellite image and below the radar and the county and state lines. Use ``radar=None`` to show only your data. See :ref:`custom_map_data`.
   :type map_layers: function or list of functions, optional, Default is ``None``.
   :param map_label: a label for the map, such as ``"HRRR MLCAPE (J/kg)"``, shown in a box at the top of the map.
   :type map_label: str, optional, Default is ``None``.
   :param return_fig: return the figure instead of showing it, so you can draw on it or save it yourself. A figure that's saved and not returned is closed.
   :type return_fig: bool, optional, Default is ``False``.
   :param readout: EXPERIMENTAL: with ``save=True``, also write a readout file next to the figure (``kokc.png`` -> ``kokc.readout.json``) for hover readouts and profile editing on a web page. See :doc:`/developers/interactive_figures`.
   :type readout: bool, optional, Default is ``False``.
   :return: ``None``, or with ``return_fig=True``, ``(fig, axes)``: the Matplotlib figure, and a dict of its panels: ``'hodograph'`` (a MetPy ``Hodograph``) and ``'map'`` (a cartopy ``GeoAxes``, when the map is drawn).
   :rtype: None or tuple


Example
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')

   # build the hodograph! 
   spy.build_hodograph(clean_data, radar='single', map_zoom=1)

.. image:: /_static/images/example-hodograph_light.png
   :width: 800 px

**************************************************








.. _compsoundings:

Building Composite Soundings
----------------------------------

Sometimes we want to compare two or more profiles against each other. Perhaps at different locations or times, or we may want to compare different models or model run-times. SounderPy allows you to do this!

To do so, a list of :ref:`'clean_data' dicts<datadescription>` is needed. If you want to customize the look of each profile, you can create equal length lists with alphas, linestyles, linewidths, & colors. See below:

.. py:function:: spy.build_composite(data_list, shade_between=True, cmap='viridis', colors=None, linestyles=None, alphas=None, linewidths=None, dark_mode=False, save=False, filename='sounderpy_composite', show_params=True, storm_motion='right_moving', max_levels=250, dpi=100, ylim=None, theme=None, style='full', sr_hodo=False, return_fig=False)


   Build a composite sounding figure of multiple profiles

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
   :param return_fig: return the figure instead of showing it, so you can draw on it or save it yourself. A figure that's saved and not returned is closed.
   :type return_fig: bool, optional, Default is ``False``.
   :return: ``None``, or with ``return_fig=True``, ``(fig, axes)``: the Matplotlib figure, and a dict of its panels: ``'skewt'`` (a MetPy ``SkewT``) and ``'hodograph'`` (a MetPy ``Hodograph``), each when drawn.
   :rtype: None or tuple

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

.. py:function:: spy.build_vad_hodograph(vad_data, save=False, dark_mode=False, storm_motion='right_moving', sr_hodo=False, filename='vad_hodograph', radar='mosaic', radar_time='sounding', map_zoom=2, hodo_boundary=None, dpi=100, theme=None, map_layers=None, map_label=None, return_fig=False)

   Build a hodograph figure of SounderPy VAD data

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
   :type radar: str, ``"mosaic"``, ``"single"`` (or ``"single-site"``), or ``None``. Optional. Default is ``mosaic``.
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
   :param map_layers: your own data for the map inset: a function, or a list of functions, that each take the map's axes (a cartopy ``PlateCarree`` axes) and draw on it, such as ``lambda ax: ax.contourf(lon, lat, values, transform=ccrs.PlateCarree(), zorder=3)``. They're drawn after the radar and before the profile's location marker; ``zorder=3`` puts filled data above the satellite image and below the radar and the county and state lines. Use ``radar=None`` to show only your data. See :ref:`custom_map_data`.
   :type map_layers: function or list of functions, optional, Default is ``None``.
   :param map_label: a label for the map, such as ``"HRRR MLCAPE (J/kg)"``, shown in a box at the top of the map.
   :type map_label: str, optional, Default is ``None``.
   :param return_fig: return the figure instead of showing it, so you can draw on it or save it yourself. A figure that's saved and not returned is closed.
   :type return_fig: bool, optional, Default is ``False``.
   :return: ``None``, or with ``return_fig=True``, ``(fig, axes)``: the Matplotlib figure, and a dict of its panels: ``'hodograph'`` (a MetPy ``Hodograph``) and ``'map'`` (a cartopy ``GeoAxes``, when the map is drawn).
   :rtype: None or tuple


VAD data
^^^^^^^^^^

``vad_data`` is a dict of NumPy arrays (without units) from your own VAD source:

- ``'z'``: height, m AGL
- ``'u'``, ``'v'``: wind components, kts
- ``'site_info'``: a dict with ``'source'`` (the top title), ``'site-id'``, ``'site-name'``, ``'site-latlon'``
  (``[lat, lon]``, which also centers the map), and ``'valid-time'`` (``[year, month, day, hour]`` as strings)

Its parameters come from ``spy.vad_params(vad_data, storm_motion='right_moving').calc()``, which returns a dict with
``'kinem'`` (shear, SRH, storm motions) and ``'intrp'`` (the profile at regular heights). Layers above the VAD's top
height are left out, with a warning.

Examples
^^^^^^^^^^

.. code-block:: python

   import numpy as np
   import sounderpy as spy

   vad_data = {'z': np.array([0, 250, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 9000]),
               'u': np.array([-5, 0, 5, 10, 15, 18, 22, 25, 28, 30, 35, 38]),
               'v': np.array([10, 20, 25, 28, 28, 27, 25, 24, 23, 22, 20, 18]),
               'site_info': {'source': 'NEXRAD VAD HODOGRAPH', 'site-id': 'KOAX', 'site-name': 'OMAHA',
                             'site-latlon': [41.32, -96.37], 'valid-time': ['2014', '06', '16', '20']}}

   spy.build_vad_hodograph(vad_data, radar='single')

.. image:: /_static/images/example_vad-hodograph.png
   :width: 800 px


.. _returned_figures:

Returned Figures & Panels
----------------------------------

``build_sounding()``, ``build_hodograph()``, ``build_composite()``, and ``build_vad_hodograph()`` show or save their
figure and return nothing, unless you pass ``return_fig=True``. Then they return ``(fig, axes)`` instead: the
Matplotlib figure, and a dictionary of its panels, so you can draw your own annotations and save it however you like.
``build_publication_sounding()`` always returns ``(fig, axes)``.

=================================  =============================================================================
Function                           Panels in ``axes``
=================================  =============================================================================
``build_sounding()``               ``'skewt'``, ``'hodograph'``, and ``'map'`` (unless ``map_zoom=0``)
``build_hodograph()``              ``'hodograph'`` and ``'map'`` (unless ``map_zoom=0``)
``build_vad_hodograph()``          ``'hodograph'`` and ``'map'`` (unless ``map_zoom=0``)
``build_composite()``              ``'skewt'`` and ``'hodograph'``, each when ``style`` draws it
``build_publication_sounding()``   ``'skewt'``, and ``'hodograph'``, ``'map'``, and ``'parameters'`` when drawn
=================================  =============================================================================

**SounderPy's figures** (``return_fig=True``):

=================  ==========================================================================
Key                Panel
=================  ==========================================================================
``'skewt'``        a MetPy ``SkewT``: draw on ``axes['skewt'].ax`` in hPa and °C
``'hodograph'``    a MetPy ``Hodograph``: draw on ``axes['hodograph'].ax`` in knots
``'map'``          a cartopy ``GeoAxes``: pass ``transform=ccrs.PlateCarree()`` with lon, lat
=================  ==========================================================================

**Publication figures:**

=================  ==========================================================================
Key                Panel
=================  ==========================================================================
``'skewt'``        a MetPy ``SkewT``: draw on ``axes['skewt'].ax`` in hPa and °C
``'hodograph'``    a Matplotlib axes, in the hodograph's units (``wind_units``)
``'map'``          a cartopy ``GeoAxes``: pass ``transform=ccrs.PlateCarree()`` with lon, lat
``'parameters'``   a Matplotlib axes holding the table (its axis lines are hidden)
=================  ==========================================================================


Good to know:

* A returned figure isn't shown, so you can finish it first. Show it with ``plt.show()``, or let a notebook show it
  at the end of the cell.
* ``save=True`` still works with ``return_fig=True``: the figure is saved, then returned.
* A figure that's saved and *not* returned is closed, so making many figures in a loop doesn't use up memory.
* A publication figure saves at exactly its width; don't add ``bbox_inches='tight'``, which trims it.

See :ref:`Drawing on the figure <drawing_on_figures>` for an example, and :doc:`/tutorials/custom_figures` for publication figures.


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

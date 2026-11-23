.. _publication_reference:

=========================
📝 Publication Figures
=========================

.. _publicationfigures:

Building Publication Figures
----------------------------------

``spy.build_publication_sounding()`` draws a clean, journal-style figure of one or more profiles: a skew-t, and
optionally a hodograph, a location map, and a table of parameters, sized for AMS and AGU journals and fitted so the
saved file is exactly the width you choose. It's **experimental** in SounderPy 3.3.0. For a walkthrough with examples,
see :doc:`/tutorials/custom_figures`.

.. py:function:: spy.build_publication_sounding(profiles, panels=('skewt', 'hodograph', 'parameters'), parameters=('mllcl', 'mlcape', 'mlcin', 'srh01', 'srh03', 'stp'), width='double', panel_labels=('(a)', '(b)', '(c)', '(d)'), colors=None, labels=None, storm_motion='right_moving', max_levels=250, ylim=None, parcel='ml', show_adiabats=False, show_barbs=True, hodo_heights=(0, 1000, 3000, 6000, 12000), hodo_max_height=12000, wind_units='m/s', hodo_limits=None, storm_motion_profile=0, sr_vector_heights=(), map_zoom=1.5, map_counties=True, map_radar=None, radar_time='sounding', parameter_profiles=None, style=None, watermark=None, dpi=300, save=False, filename='sounderpy_publication.pdf')

   EXPERIMENTAL: a clean, journal-style figure of one or more profiles: a skew-t, and any of
   a hodograph, a location map, and a table of parameters, laid out for the panels asked for
   (see :ref:`publication_layouts`).

   :param profiles: one SounderPy data dictionary, or a list of them, drawn in ``colors``.
   :type profiles: dict or list of dicts, required
   :param panels: ``'skewt'`` and any of ``'hodograph'``, ``'map'`` (each profile needs ``site-latlon``), and ``'parameters'``, labeled in that order (see :ref:`publication_layouts`).
   :type panels: tuple of str, optional, Default is ``('skewt', 'hodograph', 'parameters')``.
   :param parameters: the parameter table's columns (with the ``'parameters'`` panel): built-in names or ``(label, group, key, multiplier, unit, decimals)`` tuples. About 8 fit at the default width and font.
   :type parameters: list or None, optional, Default is ``('mllcl', 'mlcape', 'mlcin', 'srh01', 'srh03', 'stp')``.
   :param width: ``'single'``, ``'double'``, a journal's width (see :ref:`publication_widths`), or inches.
   :type width: str or float, optional, Default is ``'double'`` (6.5 in).
   :param panel_labels: labels for the panels, in the order skew-t, hodograph, map, and table, or ``None``.
   :type panel_labels: tuple of str, optional, Default is ``('(a)', '(b)', '(c)', '(d)')``.
   :param colors: a color for each profile.
   :type colors: list of str, optional, Default is ``None`` (the style's colors).
   :param labels: a legend name for each profile.
   :type labels: list of str, optional, Default is ``None`` (no legend).
   :param storm_motion: the storm motion for every profile's parameters.
   :type storm_motion: str or list of floats, optional, Default is ``'right_moving'``.
   :param max_levels: the most vertical levels used from each profile.
   :type max_levels: int, optional, Default is ``250``.
   :param ylim: the skew-t's pressure range as ``(bottom, top)`` in hPa.
   :type ylim: tuple, optional, Default is ``None`` (1050 to 100 hPa).
   :param parcel: the parcel trace: ``'sb'``, ``'ml'``, ``'mu'``, or ``None``.
   :type parcel: str or None, optional, Default is ``'ml'``.
   :param show_adiabats: draw dry and moist adiabats.
   :type show_adiabats: bool, optional, Default is ``False``.
   :param show_barbs: draw wind barb columns beside the skew-t.
   :type show_barbs: bool, optional, Default is ``True``.
   :param hodo_heights: heights (m AGL) to mark on the hodograph.
   :type hodo_heights: tuple, optional, Default is ``(0, 1000, 3000, 6000, 12000)``.
   :param hodo_max_height: where the hodograph stops (m AGL).
   :type hodo_max_height: int or float, optional, Default is ``12000``.
   :param wind_units: the hodograph's units: ``'m/s'`` or ``'kt'``.
   :type wind_units: str, optional, Default is ``'m/s'``.
   :param hodo_limits: the hodograph's limits as ``((u_min, u_max), (v_min, v_max))``.
   :type hodo_limits: tuple, optional, Default is ``None`` (fits the data).
   :param storm_motion_profile: which profile's storm motion to draw, or ``None``.
   :type storm_motion_profile: int or None, optional, Default is ``0``.
   :param sr_vector_heights: heights (m AGL) to draw storm-relative wind vectors to.
   :type sr_vector_heights: tuple, optional, Default is ``()`` (none).
   :param map_zoom: half the map's shorter side, in degrees; it grows to fit every profile.
   :type map_zoom: float, optional, Default is ``1.5``.
   :param map_counties: draw county lines on the map.
   :type map_counties: bool, optional, Default is ``True``.
   :param map_radar: ``'single'``, ``'mosaic'``, or ``None``.
   :type map_radar: str or None, optional, Default is ``None``.
   :param radar_time: ``'sounding'`` or ``'now'``.
   :type radar_time: str, optional, Default is ``'sounding'``.
   :param parameter_profiles: which profiles get a row in the table.
   :type parameter_profiles: list of int, optional, Default is ``None`` (all).
   :param style: changes to ``PUBLICATION_STYLE`` for this figure.
   :type style: dict, optional, Default is ``None``.
   :param watermark: ``None``, ``True``, or watermark settings for this figure.
   :type watermark: None, bool, or dict, optional, Default is ``None``.
   :param dpi: the figure's resolution (dots per inch).
   :type dpi: int, optional, Default is ``300``.
   :param save: save the figure to ``filename``.
   :type save: bool, optional, Default is ``False``.
   :param filename: the file to save to; its extension sets the format.
   :type filename: str, optional, Default is ``'sounderpy_publication.pdf'``.
   :return: ``(fig, axes)``: the Matplotlib figure, and a dict of its panels (``'skewt'``, ``'hodograph'``, ``'map'``, ``'parameters'``, each when drawn).
   :rtype: tuple

Example
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy

   raob = spy.get_obs_data('ABR', '2024', '08', '29', '00')                     # observed at Aberdeen, SD
   ashley = spy.get_bufkit_data('hrrr', 'kasy', 3, '2024', '08', '28', '21')    # HRRR forecast for Ashley, ND

   fig, axes = spy.build_publication_sounding([raob, ashley], labels=['00Z ABR RAOB', 'HRRR F03 Ashley, ND'],
                                              panels=('skewt', 'hodograph', 'map', 'parameters'),
                                              parameters=('mllcl', 'mlcape', 'mlcin', 'mucape', 'srh01', 'srh03', 'bwd06', 'stp'), map_radar='single',
                                              map_zoom=1)
   fig.savefig('abr_kasy_20240829.pdf')

.. image:: /_static/gallery/custom_figures/publication_quickstart.png
   :width: 100%


.. _publication_widths:

Publication Widths
^^^^^^^^^^^^^^^^^^^^

``width`` takes a name from ``PUBLICATION_WIDTHS``, or inches:

=====================  ==========  ==================================================
``width``              Inches      Fits
=====================  ==========  ==================================================
``'single'``           3.2         one column in AMS and AGU journals
``'double'``           6.5         two columns (full width) in AMS and AGU journals
``'ams_single'``       3.2         AMS, one column (19 picas)
``'ams_two_thirds'``   4.5         AMS, two-thirds of the page (27 picas)
``'ams_two_column'``   5.5         AMS, two columns (33 picas)
``'ams_full'``         6.5         AMS, full width (39 picas)
``'agu_single'``       3.35        AGU, one column (85 mm)
``'agu_double'``       6.7         AGU, two columns (170 mm)
a number               that many   any width you need
=====================  ==========  ==================================================


Every width uses the same layouts, scaled to fit, and the figure's height follows from them. ``'single'`` (3.2 in)
suits a skew-t alone; two or more panels need a two-column width, and a figure with a panel under 1.4 in across gives
a warning.


.. _publication_layouts:

Publication Layouts
^^^^^^^^^^^^^^^^^^^^^

You choose the panels with ``panels``; SounderPy places and sizes them. Each set of panels has its own predefined
layout, made to fill the figure: the skew-t, hodograph, and map each fill their space, whatever its shape (the
skew-t's isotherms stay at 45°, and the hodograph keeps one unit the same length both ways), with room kept for the
skew-t's wind barbs and the map's radar colorbar. Panels are labeled in the order skew-t, hodograph, map, and table,
whatever order they're listed in.

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - ``panels``
     - Layout
   * - ``('skewt',)``
     - the skew-t alone
   * - ``('skewt', 'hodograph')``, ``('skewt', 'map')``
     - the skew-t beside the other panel, as tall as it is
   * - ``('skewt', 'parameters')``
     - a wide skew-t over the table
   * - ``('skewt', 'hodograph', 'parameters')``, ``('skewt', 'map', 'parameters')``
     - the skew-t beside the other panel, with the table across the bottom
   * - ``('skewt', 'hodograph', 'map')``
     - a tall skew-t beside the hodograph over the map
   * - ``('skewt', 'hodograph', 'map', 'parameters')``
     - a tall skew-t beside the hodograph over the map, with the table across the bottom

The table is a strip with a column for each parameter, and the legend runs in one row across the bottom of the
figure.

.. image:: /_static/gallery/custom_figures/publication_layouts.png
   :width: 100%

*Every layout at the* ``'double'`` *width, drawn at the same scale, with the HRRR's and RAP's 3-hour forecasts for
Aberdeen, SD, valid 00 UTC 29 August 2024.*



.. _publication_parameters:

Publication Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

``parameters`` takes names from ``PUBLICATION_PARAMETERS``:

=============  ===============  ==============================================  ===============
Name           Label            Value                                           Unit
=============  ===============  ==============================================  ===============
``mllcl``      MLLCL            mixed-layer LCL height                          km
``mlcape``     MLCAPE           mixed-layer CAPE                                J kg\ :sup:`-1`
``mlcin``      MLCIN            mixed-layer CIN (negative)                      J kg\ :sup:`-1`
``sbcape``     SBCAPE           surface-based CAPE                              J kg\ :sup:`-1`
``mucape``     MUCAPE           most-unstable CAPE                              J kg\ :sup:`-1`
``muecape``    MU SR-ECAPE      most-unstable storm-relative entraining CAPE    J kg\ :sup:`-1`
``srh01``      0–1 km SRH       0–1 km storm-relative helicity                  m\ :sup:`2` s\ :sup:`-2`
``srh03``      0–3 km SRH       0–3 km storm-relative helicity                  m\ :sup:`2` s\ :sup:`-2`
``bwd06``      0–6 km BWD       0–6 km bulk wind difference                     m s\ :sup:`-1`
``stp``        STP              effective-layer significant tornado parameter   none
``scp``        SCP              effective-layer supercell composite parameter   none
=============  ===============  ==============================================  ===============


Any other value ``spy.sounding_params()`` calculates can be a column, given as a tuple of six items. Names and tuples
can be mixed in one list.

.. code-block:: python

   ('0–3 km SBCAPE', 'thermo', 'sb3cape', 1, r'J kg$^{-1}$', 0)
   #  label          group     key       multiplier  unit     decimals

* **label**: the column's name.
* **group** and **key**: where the value is in ``spy.sounding_params(clean_data).calc()``. The groups are
  ``'general'``, ``'thermo'``, ``'kinem'``, and ``'intrp'``; every value, by group and key, is listed on the
  :doc:`Parameters & References </reference/references>` page.
* **multiplier**: a unit conversion, such as ``0.001`` for m to km, or ``0.514444`` for knots to m/s. Use ``1`` to
  keep the value as is.
* **unit**: the unit shown under the name; Matplotlib math text such as ``r'm$^2$ s$^{-2}$'`` works.
* **decimals**: how many decimal places to show.


.. _publication_style:

Publication Style
^^^^^^^^^^^^^^^^^^^^

Every part of a publication figure's look is a setting in ``PUBLICATION_STYLE``. Change any of them for one figure
with ``style={...}``, or for every figure in a script by changing ``PUBLICATION_STYLE`` itself
(``from sounderpy.plot.publication import PUBLICATION_STYLE``).

.. list-table:: Style settings
   :header-rows: 1
   :widths: 28 24 48

   * - Setting
     - Default
     - What it sets
   * - ``font_family``
     - ``'sans-serif'``
     - the font, such as ``'serif'`` or a font's name
   * - ``fontsize``
     - ``8``
     - text size (pt); AMS and most AGU journals ask for about 8 pt, GRL for 6 pt
   * - ``panel_label_fontsize``
     - ``10``
     - the size of ``(a)``, ``(b)``, ``(c)``
   * - ``colors``
     - red, blue, black, orange, purple, green
     - the profiles' colors, in order
   * - ``linewidth``
     - ``1.1``
     - temperature and dewpoint lines
   * - ``parcel_linewidth`` / ``parcel_linestyle``
     - ``0.9`` / dashed
     - the parcel trace
   * - ``temperature_range``
     - ``(-40, 40)``
     - the skew-t's temperature range (°C); also used when zoomed, if you set it
   * - ``grid_color`` / ``grid_linewidth``
     - ``'#bdbdbd'`` / ``0.6``
     - isobars, isotherms, and adiabats
   * - ``axis_linewidth``
     - ``0.8``
     - axis lines and ticks
   * - ``barb_length`` / ``barb_linewidth`` / ``barb_spacing``
     - ``4.5`` / ``0.5`` / ``0.17``
     - wind barbs; the spacing is between profiles' columns, in inches
   * - ``hodograph_linewidth``
     - ``1.1``
     - hodograph traces
   * - ``height_marker_size`` / ``height_label_fontsize`` / ``height_label_weight``
     - ``7`` / ``5.5`` / ``'bold'``
     - height markers and their labels (for GRL's 6 pt, use ``8`` / ``6``)
   * - ``height_marker_fill``
     - see-through white
     - the inside of each height marker
   * - ``surface_label``
     - ``'S'``
     - the surface marker's label, such as ``'SFC'``
   * - ``hodograph_rings`` / ``ring_color`` / ``ring_linewidth``
     - ``True`` / ``'#d4d4d4'`` / ``0.5``
     - range rings on the hodograph
   * - ``storm_motion_color`` / ``storm_motion_linewidth``
     - ``'black'`` / ``1.2``
     - the storm motion arrow
   * - ``sr_vector_color`` / ``sr_vector_linewidth``
     - ``'gray'`` / ``1.2``
     - storm-relative wind vectors
   * - ``map_line_color`` / ``map_marker_size``
     - ``'#555555'`` / ``30``
     - the map's lines and profile dots
   * - ``radar_colorbar`` / ``radar_range``
     - ``True`` / ``(10, 75)``
     - the radar colorbar, and its dBZ range
   * - ``watermark_fontsize`` / ``watermark_color``
     - ``6`` / ``'#9a9a9a'``
     - the watermark's text

See :ref:`returned_figures` for drawing on the returned panels, and :doc:`/tutorials/custom_figures` for a
walkthrough with examples.

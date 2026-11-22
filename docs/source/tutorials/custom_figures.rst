.. _tutorial-custom-figures:

==================================================
📝 Building Publication & Custom Figures
==================================================

SounderPy gives you two ways to make a figure of your own:

* **Customize a SounderPy figure.** Every ``build_*`` function can hand back its
  Matplotlib figure and panels with ``return_fig=True``, so you can draw on it and
  save it however you like.
* **Build a publication figure.** ``spy.build_publication_sounding()`` draws a clean,
  journal-style figure of one or more profiles: a skew-t, and optionally a hodograph,
  a location map, and a table of parameters. It's sized for AMS and AGU journals,
  every part of its look can be changed, and it returns its panels for annotating.

This tutorial covers:

* drawing on SounderPy's own figures with ``return_fig``;
* making a publication figure, and choosing its panels and width;
* the skew-t, hodograph, parameter table, and map options;
* changing the style, adding a watermark, and saving for a journal;
* drawing your own annotations on the returned panels.

.. image:: /_static/gallery/custom_figures/publication_quickstart.png
   :width: 100%

*The observed Aberdeen, SD sounding and the HRRR's 3-hour forecast for Ashley, ND, both valid 00 UTC 29 August 2024
(see the* :doc:`Mound City case study </case_studies/mound_city_20240828>`\ *); the map shows where each profile is.*

.. note::

   ``build_publication_sounding()`` is **experimental** in SounderPy 3.3.0: its
   options may change in a later release.


***************************************************************

.. _custom-figures-drawing:

Part 1: Drawing on SounderPy's Figures
======================================

``build_sounding()``, ``build_hodograph()``, ``build_composite()``, and
``build_vad_hodograph()`` hand back their figure and its panels with ``return_fig=True``,
so you can draw your own annotations on them and save the figure however you like. The
plotting reference covers it:

* :ref:`Drawing on the figure <drawing_on_figures>`: an example that shades the elevated mixed layer on this
  case's HRRR sounding;
* :ref:`returned_figures`: the panels each function returns, and when a figure is shown,
  saved, or closed.

SounderPy's figures can also be changed without drawing on them: ``theme`` sets their
colors (see :ref:`plot_themes`), and ``map_layers`` draws your own data on the map (see
:ref:`custom_map_data` and the :doc:`custom_map_data` tutorial).


***************************************************************

.. _custom-figures-publication:

Part 2: Publication Figures
===========================

Quick Start
-----------

Pass one profile, or a list of them. Any SounderPy data works: observed soundings,
model forecasts, BUFKIT, ACARS, or your own data.

.. code-block:: python

   import sounderpy as spy

   raob = spy.get_obs_data('ABR', '2024', '08', '29', '00')                     # observed at Aberdeen, SD
   ashley = spy.get_bufkit_data('hrrr', 'kasy', 3, '2024', '08', '28', '21')    # HRRR forecast for Ashley, ND

   fig, axes = spy.build_publication_sounding([raob, ashley], labels=['00Z ABR RAOB', 'HRRR F03 Ashley, ND'],
                                              panels=('skewt', 'hodograph', 'map', 'parameters'),
                                              parameters=('mllcl', 'mlcape', 'mlcin', 'mucape', 'srh01', 'srh03', 'bwd06', 'stp'),
                                              map_radar='single')

   fig.savefig('abr_kasy_20240829.pdf')

That's the figure at the top of this page. Each profile gets its own color, on
every panel, and the map shows where each one is.

The rest of this tutorial compares the HRRR's and RAP's forecasts for Aberdeen. With one
profile and the defaults you get a skew-t, a hodograph, and a short parameter table:

.. code-block:: python

   hrrr = spy.get_bufkit_data('hrrr', 'kabr', 3, '2024', '08', '28', '21')
   rap = spy.get_bufkit_data('rap', 'kabr', 3, '2024', '08', '28', '21')

   fig, axes = spy.build_publication_sounding(hrrr)

.. image:: /_static/gallery/custom_figures/publication_default.png
   :width: 100%


Panels & Layout
---------------

``panels`` names the panels to draw: ``'skewt'`` and any of ``'hodograph'``, ``'map'``, and
``'parameters'`` (the default is the skew-t, hodograph, and parameter table). Each set of
panels has its own layout, and each panel fills its space, whatever its shape:

* the skew-t on the left, with the hodograph or map beside it, or both stacked (the
  hodograph on top);
* the table in a strip across the bottom, with a column for each parameter (the figure
  at the top of this page has all four).

You only choose the panels; SounderPy places and sizes them. Each layout is listed in
:ref:`publication_layouts`.

.. code-block:: python

   fig, axes = spy.build_publication_sounding(hrrr, panels=('skewt', 'hodograph'))

.. image:: /_static/gallery/custom_figures/publication_pair.png
   :width: 100%

Panels are labeled ``(a)``, ``(b)``, ``(c)``, ``(d)`` in the order skew-t, hodograph, map,
and table, whatever order they're listed in. Change the labels with
``panel_labels=('a)', 'b)', 'c)', 'd)')``, or leave them off with ``panel_labels=None``.


Width
-----

``width`` sets the figure's width to fit a journal page. Every width uses the same
layouts, scaled to fit; the figure's height follows from the layout, and everything drawn
is fitted inside the figure, so the saved file is exactly this wide.

``'double'`` (6.5 in, the default) fits two columns in AMS and AGU journals and suits any
layout; ``'single'`` (3.2 in) fits one column and suits a skew-t alone. Every width,
including each journal's own, is listed in :ref:`publication_widths`.

.. tip::

   Panels under 1.4 in across give a warning; use a wider figure, or fewer panels.


Profiles, Colors & Legend
-------------------------

* ``colors``: a color for each profile, such as ``['black', 'darkorange']``. The
  default colors are red, blue, black, orange, purple, and green, in that order.
* ``labels``: a name for each profile, such as ``['18Z RAOB', 'HRRR F03']``. They're
  shown in a legend in one row across the bottom of the figure, under the panels (more
  rows only if one row is wider than the figure).
* ``storm_motion``: the storm motion used for every profile's parameters:
  ``'right_moving'`` (default), ``'left_moving'``, ``'mean_wind'``, or ``[direction, speed]`` in
  degrees and knots (see :ref:`storm_motions`).
* ``max_levels``: the most vertical levels used from each profile (default ``250``).


Skew-T Options
--------------

* ``ylim``: the pressure range as ``(bottom, top)`` in hPa, such as ``(1000, 500)``.
  A zoomed skew-t fits its temperature range to the data.
* ``parcel``: the parcel trace drawn for each profile: ``'ml'`` (default), ``'sb'``,
  ``'mu'``, or ``None``. A parcel is drawn only when it has CAPE.
* ``show_adiabats``: draw dry and moist adiabats (default ``False``).
* ``show_barbs``: draw a column of wind barbs for each profile beside the skew-t
  (default ``True``). Barbs are in knots.

.. code-block:: python

   fig, axes = spy.build_publication_sounding(
       [hrrr, rap], labels=['HRRR F03', 'RAP F03'], panels=('skewt', 'parameters'),
       ylim=(1000, 500), parcel='sb',
       parameters=('sbcape', 'mlcin', ('0–3 km SBCAPE', 'thermo', 'sb3cape', 1, r'J kg$^{-1}$', 0), 'mllcl'))

.. image:: /_static/gallery/custom_figures/publication_zoom_parameters.png
   :width: 100%


Hodograph Options
-----------------

* ``wind_units``: ``'m/s'`` (default) or ``'kt'``.
* ``hodo_heights``: heights (m AGL) to mark with a labeled circle. Default
  ``(0, 1000, 3000, 6000, 12000)``; ``0`` is labeled ``S`` (the surface), the rest in
  km. Markers that would cover each other move aside, with a thin line to their point.
* ``hodo_max_height``: where the hodograph trace stops (m AGL, default ``12000``).
* ``hodo_limits``: the axis limits as ``((u_min, u_max), (v_min, v_max))``. The default
  fits every profile and storm motion.
* ``storm_motion_profile``: which profile's storm motion arrow to draw (default ``0``,
  the first), or ``None`` for none.
* ``sr_vector_heights``: heights (m AGL) to draw **storm-relative wind vectors** to.
  Each is an arrow from the tip of the storm motion arrow to the wind at that height:
  the wind a storm moving with that motion feels there. A vector shorter than its
  height marker (a storm-relative wind near zero) is left out. Default: none.

.. code-block:: python

   fig, axes = spy.build_publication_sounding(
       hrrr, wind_units='kt', hodo_heights=(0, 500, 1000, 3000, 6000, 9000), hodo_max_height=9000,
       sr_vector_heights=(500, 1000, 3000), parameters=('srh01', 'srh03', 'bwd06', 'scp'))

.. image:: /_static/gallery/custom_figures/publication_hodograph.png
   :width: 100%


Parameter Table
---------------

``parameters`` lists the table's columns, in order: each has the parameter's name and
unit on top, then each profile's value in its own color; a missing value shows ``NA``. The
default is ``('mllcl', 'mlcape', 'mlcin', 'srh01', 'srh03', 'stp')``.

Use the **built-in names** listed in :ref:`publication_parameters`.

**Any other value** that ``spy.sounding_params()`` calculates can be a column, given as a
tuple of six items:

.. code-block:: python

   ('0–3 km SBCAPE', 'thermo', 'sb3cape', 1, r'J kg$^{-1}$', 0)
   #  label          group     key       multiplier  unit     decimals

Each item is explained in :ref:`publication_parameters`.

Names and tuples can be mixed in one list. Other options:

* Leave ``'parameters'`` out of ``panels`` to leave out the table.
* ``parameter_profiles``: which profiles get a row, such as ``[0]`` for the first
  only. Default: all of them.

Each column is as wide as its widest entry, and the room left over is shared between
them, so a longer list fills the strip out: about 8 fit at the default width and font. A
list that doesn't fit gives an error.


Map & Radar
-----------

The ``'map'`` panel draws a map of where each profile is, a dot in each profile's color, with
county and state lines and no satellite image. Each profile needs its location in
``site_info['site-latlon']``.

* ``map_zoom``: half the map's shorter side, in degrees (default ``1.5``); it grows to fit every
  profile.
* ``map_counties``: draw county lines (default ``True``).
* ``map_radar``: radar reflectivity on the map: ``'single'`` (the NEXRAD site nearest
  the first profile), ``'mosaic'`` (the CONUS mosaic; only for about the last month),
  or ``None`` (default). A dBZ colorbar goes under the map, with the scan time.
* ``radar_time``: the radar's time: ``'sounding'`` (the first profile's valid time,
  default) or ``'now'``.

.. code-block:: python

   fig, axes = spy.build_publication_sounding([hrrr, rap], labels=['HRRR F03', 'RAP F03'],
                                              panels=('skewt', 'hodograph', 'map'), map_radar='single')

.. image:: /_static/gallery/custom_figures/publication_radar.png
   :width: 100%

The map's layers, bottom to top, are: radar (``zorder`` 4), county lines (5), state
lines (6), and profile dots (10). A field you draw on ``axes['map']`` sits under the
radar and lines at the default ``zorder``; use ``zorder=4.5`` to put it over the radar.


Style
-----

Everything about how the figure looks is in one dictionary of settings,
``PUBLICATION_STYLE``. Change any of them for one figure with ``style``:

.. code-block:: python

   fig, axes = spy.build_publication_sounding(
       [hrrr, rap], labels=['HRRR F03', 'RAP F03'], colors=['black', 'darkorange'],
       show_adiabats=True, parcel='mu', sr_vector_heights=(1000,),
       style={'font_family': 'serif', 'fontsize': 9, 'linewidth': 1.4,
              'hodograph_rings': False, 'sr_vector_color': 'purple'})

.. image:: /_static/gallery/custom_figures/publication_style.png
   :width: 100%

Every setting, with its default, is listed in :ref:`publication_style`.

To change a setting for every figure in a script, change ``PUBLICATION_STYLE`` itself:

.. code-block:: python

   from sounderpy.plot.publication import PUBLICATION_STYLE

   PUBLICATION_STYLE['fontsize'] = 6   # every figure from here on, for GRL


Watermark
---------

Publication figures have no watermark unless you ask for one:

* ``watermark=True``: SounderPy's credit text (or your own, set with
  ``spy.set_watermark()``; see :ref:`watermark`) in small gray text at the bottom right,
  in a row of its own under the panels and legend (a logo goes at the bottom left).
* ``watermark={...}``: your own, for this figure only, such as
  ``{'title': 'UND Atmospheric Sciences', 'subtitle': '', 'logo': 'logo.png'}``. Leave a
  setting out to use the ``spy.set_watermark()`` one. A logo is only drawn when you give
  one.


Saving for a Journal
--------------------

.. code-block:: python

   fig.savefig('figure1.pdf')                 # or .eps, .svg, .png, .tif

* The file is exactly the figure's size: the ``width`` you chose. Don't add
  ``bbox_inches='tight'``, which trims the figure and changes its size.
* PDF, EPS, and SVG files keep text as real text (Type 42 fonts), as AMS and AGU ask,
  even when you save the figure later.
* ``dpi`` (default ``300``) sets the resolution of PNG and TIFF files, and of a figure
  shown in a notebook. Use ``600`` or more for line art in a raster format.
* Or save as part of the call: ``save=True, filename='figure1.pdf'``.


Drawing on the Figure
---------------------

``build_publication_sounding()`` always returns ``(fig, axes)``: the figure, and a
dictionary of its panels (see :ref:`returned_figures`).

.. code-block:: python

   import cartopy.crs as ccrs

   fig, axes = spy.build_publication_sounding([hrrr, rap], labels=['HRRR F03', 'RAP F03'],
                                              panels=('skewt', 'hodograph', 'map'), map_zoom=2)

   # shade a layer on the skew-t
   axes['skewt'].ax.axhspan(850, 700, color='gold', alpha=0.25, zorder=0)
   axes['skewt'].ax.text(-35, 775, 'EML', fontsize=8, va='center')

   # a note on the hodograph, in m/s
   axes['hodograph'].annotate('0–3 km SRH\n342 m$^2$ s$^{-2}$', xy=(18, -6), fontsize=6,
                              ha='center', va='top', style='italic')

   # a place on the map
   axes['map'].plot(-100.07, 45.73, marker='*', color='black', markersize=7,
                    transform=ccrs.PlateCarree(), zorder=11)
   axes['map'].text(-99.9, 45.73, 'Mound City', fontsize=6, va='center',
                    transform=ccrs.PlateCarree(), zorder=11)

   fig.savefig('figure1.pdf')

.. image:: /_static/gallery/custom_figures/publication_annotated.png
   :width: 100%

.. tip::

   Draw your annotations inside the panels. The figure is fitted to its panels when
   it's built, so text placed far outside them can be cut off when it's saved.


***************************************************************

Reference
=========

Every option of ``build_publication_sounding()``, its settings, and the panels each
figure returns are in the plotting reference: :ref:`publicationfigures` and
:ref:`returned_figures`. For a full study built on a publication figure, see the
:doc:`CM1 surface heating case study </case_studies/cm1_surface_heating>`.

.. _plot_settings:

==================
🎨 Plot Settings
==================

The settings SounderPy's figures share: their colors and themes, the watermark, surface modification, storm motions,
parcels, and the map and radar (including your own map data). Each figure's own options are on
:doc:`/reference/plottingdata`.

***************************************************************

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
- **Hodograph**: ``hodo_layers``, ``height_marker``, ``height_marker_edge``, ``mcs``, ``dtm``, ``mean_wind``
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
.. _surface-modification-logic:

Surface Modification
---------------------

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

To get the modified profile itself, such as to save it or compare it with the original, use
``spy.modify_surface(clean_data, modify_sfc)`` (see :ref:`modify_surface`).

.. caution::
   Surface modification should be used carefully. Replacing the surface
   thermodynamic or wind state can materially alter parcel and kinematic
   calculations.

**************************************************








.. _storm_motions:
.. _storm-motion-logic:

Storm Motions
--------------

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
.. _parcel-logic:

Parcels
-------

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


By default, SounderPy will plot normal MU/ML/SB-CAPE parcels and an mu_ps_ecape parcel. You can override this by setting ``special_parcels`` to 'simple', which only plots the common MU/ML/SB-CAPE parcels. This can greatly reduce the plot-time!


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

Map & Radar
-----------
Knowing the location of a vertical profile relative to nearby meteorological phenomenon is vital in understanding what that vertical profile means for the atmosphere at a given time. To help with this, SounderPy employs the use of mapping & radar data to provide this information to a reader.

.. _plot_settings_map:

Using the Map
^^^^^^^^^^^^^

The map inset is controlled with ``map_zoom`` and ``map_background``.

- ``map_zoom=0`` disables the map.
- Larger positive values display a broader geographic area.
- ``map_background='map'`` (the default) draws land and water with city areas, lakes, rivers, major highways, and
  the names of up to five of the largest cities. It's drawn from `Natural Earth <https://www.naturalearthdata.com/>`_
  data, which cartopy downloads the first time it's used.
- ``map_background='satellite'`` draws satellite imagery instead. It's downloaded for each figure, so it's slower.

The map's colors can be changed with ``theme``: ``map_land``, ``map_water``, ``map_urban``, ``map_river``,
``map_road``, and ``map_city``.

.. code-block:: python

   spy.build_sounding(clean_data, map_background='satellite')


Using Radar
^^^^^^^^^^^^^^^^^
SounderPy offers two sources of radar data to the user. First is a CONUS reflectivity mosaic, and the other is single-site NEXRAD WSR-88D reflectivity. The CONUS mosaic is default because it's the fastest and easiest to plot, but data availability is only up to 30 days before the current date. Single-site radar goes back much further but data availability varies by site and date.

To plot radar on a sounding or hodograph, there are two kwargs to understand:

- ``radar``

  - ``"mosaic"``: plots the CONUS reflectivity mosaic
  - ``"single"`` (or ``"single-site"``): plots single-site WSR-88D reflectivity closest to the latitude/longitude coordinates of the sounding.
  - ``None``: no radar.

- ``radar_time``

  - ``"now"``: plots the most recent radar scan
  - ``"sounding"``: plots the radar scan closest the time of the sounding
  - ``datetime(1999, 5, 4, 0, 27)``: or you may define your own time by providing a datetime obj. This will plot the nearest radar scan to the date/time provided.


.. _custom_map_data:

Custom Map Data
^^^^^^^^^^^^^^^^^

The map inset can show your own data, such as surface observations, a model field, or anything else with a latitude
and longitude. ``build_sounding()`` (severe and winter), ``build_hodograph()``, and ``build_vad_hodograph()`` take two
kwargs for it (from Python only, not the command line tool):

- ``map_layers``: a function, or a list of functions, that each take the map's axes (a cartopy ``PlateCarree``
  axes) and draw on it with any Matplotlib or cartopy call: ``contourf``, ``contour``, ``pcolormesh``, ``scatter``,
  ``text``, ``barbs``, and so on.
- ``map_label``: a short label, such as ``'RAP SBCAPE (J/kg)'``, shown in a box at the top of the map.

.. code-block:: python

   import cartopy.crs as ccrs
   import sounderpy as spy

   clean_data = spy.get_obs_data('ABR', '2024', '08', '29', '00')

   # lon, lat, cape, and mslp are your own 2-D arrays
   def plot_cape(ax):
       ax.contourf(lon, lat, cape, levels=range(250, 5001, 250), cmap='magma_r', alpha=0.6,
                   transform=ccrs.PlateCarree(), zorder=3)

   def plot_mslp(ax):
       ax.contour(lon, lat, mslp, levels=range(980, 1041, 2), colors='black', linewidths=1,
                  transform=ccrs.PlateCarree(), zorder=7)

   spy.build_sounding(clean_data, radar=None, map_zoom=3, map_layers=[plot_cape, plot_mslp],
                      map_label='RAP SBCAPE (J/kg) & MSLP (2 hPa)')

.. image:: /_static/gallery/custom_map/rap_cape_mslp.png
   :width: 100%

Pass ``transform=ccrs.PlateCarree()`` so your data's x and y are read as longitude and latitude, and use
``zorder`` to stack it with the map. Your layers are drawn after the radar and before the profile's location marker:

==========================================  ==========
Map part                                    zorder
==========================================  ==========
Land, or the satellite image                1
City areas and lakes                        1.2 - 1.4
Highways and rivers; satellite shading      2 - 2.5
Radar                                       4
County lines                                5
State lines                                 6
City names                                  7
Location marker                             9.9 - 11
==========================================  ==========

``zorder=3`` puts a filled field above the map background and under the radar and borders; ``zorder=7`` or ``8``
puts lines, points, and text on top of the borders. Use ``radar=None`` to show only your data. ``map_layers`` and
``map_label`` need the map, so leave ``map_zoom`` above ``0``.

A publication figure's map (``build_publication_sounding(panels=(..., 'map'))``) is drawn on directly instead: it's returned
as ``axes['map']`` (see :ref:`returned_figures`).

For worked examples with METAR observations from the Iowa Environmental Mesonet and a RAP field from NOAA's archive
on AWS, see the :doc:`Custom Map Data tutorial </tutorials/custom_map_data>`.

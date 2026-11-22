.. _plot-gallery:

================
📊 Plot Gallery
================

SounderPy can create sounding, hodograph, composite, and VAD-based visualizations
from observed, forecast, reanalysis, model-output, and custom vertical profile
data.

Most comparison examples use the same observed OAX sounding from
16 June 2014 at 18 UTC so that visual differences reflect plotting options
rather than differences in the atmospheric profile.

**Jump to:** :ref:`Appearance <gallery-appearance>` |
:ref:`Parcels <gallery-parcels>` |
:ref:`Profile Tools <gallery-profile-tools>` |
:ref:`Surface Modification <gallery-surface>` |
:ref:`Hodographs <gallery-hodographs>` |
:ref:`Maps & Radar <gallery-map-radar>` |
:ref:`Your Own Data <gallery-customization>` |
:ref:`Data Sources <gallery-data-sources>` |
:ref:`Winter Soundings <gallery-winter>` |
:ref:`Composites <gallery-composites>` |
:ref:`Publication Figures <gallery-publication>`


.. _gallery-appearance:

Appearance & Accessibility
==========================

The same profile in light and dark mode, with color-deficiency-friendly colors (the dewpoint
trace turns from green to blue), and with a custom ``theme``.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/themes/sounding_light.png
          :alt: Light-mode SounderPy sounding
          :width: 100%

       **Light Mode**

       Default SounderPy appearance.

       .. code-block:: python

          spy.build_sounding(
              data,
              dark_mode=False,
          )

     - .. image:: /_static/gallery/themes/sounding_dark.png
          :alt: Dark-mode SounderPy sounding
          :width: 100%

       **Dark Mode**

       Enable dark mode with ``dark_mode=True``.

       .. code-block:: python

          spy.build_sounding(
              data,
              dark_mode=True,
          )

   * - .. image:: /_static/gallery/accessibility/sounding_colorblind.png
          :alt: SounderPy sounding using color-deficiency-friendly colors
          :width: 100%

       **Color-Deficiency-Friendly**

       The dewpoint trace turns blue.

       .. code-block:: python

          spy.build_sounding(
              data,
              color_blind=True,
          )

     - .. image:: /_static/gallery/themes/sounding_custom_theme.png
          :alt: SounderPy sounding drawn with a custom parchment theme
          :width: 100%

       **Custom Theme**

       A parchment palette to match these docs, on the Mound City case (HRRR F03 for Ashley, ND). See :ref:`theme settings <plot_themes>` for every name.

       .. code-block:: python

          pastoral = {
              'background': '#f5edda',
              'text': '#241c15',
              'temperature': '#a3241a',
              'dewpoint': '#3c6a1c',
              'parcel': {'mu': '#5a2a73'},
              # ... full palette in
              # make_gallery_figures.py
          }
          spy.build_sounding(
              data,
              theme=pastoral,
          )


.. _gallery-parcels:

Parcels & Thermodynamics
========================

Parcel options, the classic low-level theta and theta-e profiles, and a skew-t zoomed into the
lower atmosphere with ``ylim``.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/parcels/parcels_simple.png
          :alt: SounderPy sounding using simple parcel settings
          :width: 100%

       **Simple Parcels**

       The simplified parcel configuration.

       .. code-block:: python

          spy.build_sounding(
              data,
              special_parcels="simple",
          )

     - .. image:: /_static/gallery/parcels/parcels_custom.png
          :alt: SounderPy sounding using custom parcel settings
          :width: 100%

       **Custom Parcels**

       Highlight some parcel types and draw others in the background.

       .. code-block:: python

          custom_parcels = [
              ["sb_ia_ecape"],
              ["sb_ps_ecape", "sb_ps_cape"],
          ]

          spy.build_sounding(
              data,
              special_parcels=custom_parcels,
          )

   * - .. image:: /_static/gallery/options/show_theta.png
          :alt: SounderPy sounding displaying theta and theta-e
          :width: 100%

       **Theta and Theta-e**

       The classic low-level theta and theta-e profiles, with ``show_theta=True``.

       .. code-block:: python

          spy.build_sounding(
              data,
              show_theta=True,
          )

     - .. image:: /_static/gallery/profile_tools/zoom_ylim.png
          :alt: SounderPy sounding with the skew-t zoomed to 1000-500 hPa
          :width: 100%

       **Zoomed Skew-T**

       Zoom the skew-t into a layer with ``ylim`` (bottom, top) in hPa.

       .. code-block:: python

          spy.build_sounding(
              data,
              ylim=(1000, 500),
          )


.. _gallery-profile-tools:

Profile Tools
=============

Smooth a profile with a Savitzky-Golay or Gaussian filter, merge two profiles by weighted
averaging, or thin a profile to fewer levels with ``max_levels``. See :doc:`/reference/helpertools`.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/profile_tools/smooth_savgol.png
          :alt: SounderPy sounding smoothed with a Savitzky-Golay filter
          :width: 100%

       **Savitzky-Golay Smoothing**

       The default filter; a wider window smooths more.

       .. code-block:: python

          smoothed = spy.smooth_profile(
              data,
              method='savgol',
              window_length=15,
          )
          spy.build_sounding(smoothed)

     - .. image:: /_static/gallery/profile_tools/smooth_gaussian.png
          :alt: SounderPy sounding smoothed with a Gaussian filter
          :width: 100%

       **Gaussian Smoothing**

       A larger ``sigma`` smooths more.

       .. code-block:: python

          smoothed = spy.smooth_profile(
              data,
              sigma=3,
          )
          spy.build_sounding(smoothed)

   * - .. image:: /_static/gallery/profile_tools/merged.png
          :alt: SounderPy sounding of two merged profiles
          :width: 100%

       **Merged Profiles**

       60% of the 18Z OAX sounding and 40% of the 00Z one.

       .. code-block:: python

          merged = spy.merge_profiles(
              oax_18z, oax_00z,
              weight_a=0.6,
          )
          # add site_info and titles, then
          spy.build_sounding(merged)

     - .. image:: /_static/gallery/profile_tools/max_levels.png
          :alt: SounderPy sounding thinned to 30 levels
          :width: 100%

       **Fewer Levels**

       Thin a profile to evenly spaced levels with ``max_levels``.

       .. code-block:: python

          spy.build_sounding(
              data,
              max_levels=30,
          )


.. _gallery-surface:

Surface Modification
====================

The 12Z OAX sounding's surface inversion, with its surface warmed to the 18Z observation
(28 °C, dewpoint 24 °C) by each ``modify_sfc`` method, zoomed to 1000-600 hPa. See
:ref:`surface modification <modify_sfc>`.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/surface_modification/sfc_none.png
          :alt: Unmodified SounderPy sounding with a surface inversion
          :width: 100%

       **Unmodified**

       The morning sounding as observed.

       .. code-block:: python

          spy.build_sounding(
              data,
              ylim=(1000, 600),
          )

     - .. image:: /_static/gallery/surface_modification/sfc_mix.png
          :alt: SounderPy sounding with the surface mixed up its dry adiabat
          :width: 100%

       **Mix (default)**

       The warmer surface mixes up its dry adiabat until it meets the profile.

       .. code-block:: python

          spy.build_sounding(
              data,
              modify_sfc={'T': 28, 'Td': 24},
              ylim=(1000, 600),
          )

   * - .. image:: /_static/gallery/surface_modification/sfc_blend.png
          :alt: SounderPy sounding with the surface change blended over 500 m
          :width: 100%

       **Blend**

       The change tapers to zero at ``depth`` m above ground (500 by default).

       .. code-block:: python

          spy.build_sounding(
              data,
              modify_sfc={'T': 28, 'Td': 24,
                          'method': 'blend'},
              ylim=(1000, 600),
          )

     - .. image:: /_static/gallery/surface_modification/sfc_replace.png
          :alt: SounderPy sounding with only the surface level replaced
          :width: 100%

       **Replace**

       Only the surface level changes, keeping a superadiabatic layer.

       .. code-block:: python

          spy.build_sounding(
              data,
              modify_sfc={'T': 28, 'Td': 24,
                          'method': 'replace'},
              ylim=(1000, 600),
          )


.. _gallery-hodographs:

Hodographs
==========

Ground- or storm-relative hodographs (``sr_hodo``), a custom storm motion, and a boundary line
to show the flow relative to a front or outflow boundary. Each works on hodograph and sounding plots.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/options/hodograph_ground_relative.png
          :alt: Ground-relative SounderPy hodograph
          :width: 100%

       **Ground Relative**

       The default.

       .. code-block:: python

          spy.build_hodograph(
              data,
              sr_hodo=False,
          )

     - .. image:: /_static/gallery/options/hodograph_storm_relative.png
          :alt: Storm-relative SounderPy hodograph
          :width: 100%

       **Storm Relative**

       Winds relative to the storm motion.

       .. code-block:: python

          spy.build_hodograph(
              data,
              sr_hodo=True,
          )

   * - .. image:: /_static/gallery/options/hodograph_custom_storm_motion.png
          :alt: SounderPy hodograph using custom storm motion
          :width: 100%

       **Custom Storm Motion**

       ``storm_motion = [direction, speed]``.

       .. code-block:: python

          spy.build_hodograph(
              data,
              storm_motion=[250.0, 45.0],
          )

     - .. image:: /_static/gallery/options/hodograph_boundary.png
          :alt: SounderPy hodograph with a custom boundary line
          :width: 100%

       **Hodograph Boundary**

       A line at an angle from north, here 45°.

       .. code-block:: python

          spy.build_hodograph(
              data,
              hodo_boundary={
                  "angle": [45],
                  "color": ["tab:purple"],
              },
          )


.. _gallery-map-radar:

Maps & Radar
============

The map inset with or without radar, a wider map extent, and single-site or mosaic radar.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/map_radar/map_only.png
          :alt: SounderPy sounding with map enabled and radar disabled
          :width: 100%

       **Map Only**

       Display the map without radar.

       .. code-block:: python

          spy.build_sounding(
              data,
              radar=None,
              map_zoom=2,
          )

     - .. image:: /_static/gallery/map_radar/map_zoom_4.png
          :alt: SounderPy sounding with a larger map extent
          :width: 100%

       **Map Zoom**

       A larger ``map_zoom`` shows more of the map; ``map_zoom=0`` removes it.

       .. code-block:: python

          spy.build_sounding(
              data,
              radar=None,
              map_zoom=4,
          )

   * - .. image:: /_static/gallery/map_radar/radar_single.png
          :alt: SounderPy sounding with single-site radar
          :width: 100%

       **Single-Site Radar**

       WSR-88D data from the nearest radar site.

       .. code-block:: python

          spy.build_sounding(
              data,
              radar="single",
              radar_time="sounding",
          )

     - .. image:: /_static/gallery/map_radar/radar_mosaic.png
          :alt: SounderPy sounding with CONUS radar mosaic
          :width: 100%

       **Radar Mosaic**

       Recent CONUS mosaic reflectivity.

       .. code-block:: python

          # a sounding from the last ~30 days
          spy.build_sounding(
              recent_data,
              radar="mosaic",
              radar_time="sounding",
          )


.. _gallery-customization:

Your Own Data & Annotations
===========================

Draw your own data on the map inset with ``map_layers``, or take the figure back with
``return_fig=True`` and annotate it. See :doc:`/tutorials/custom_map_data` and :doc:`/tutorials/custom_figures`.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/custom_map/rap_cape_mslp.png
          :alt: SounderPy sounding with RAP CAPE and isobars on the map inset
          :width: 100%

       **RAP CAPE & Pressure**

       A RAP field from NOAA's archive on AWS, filled and contoured.

       .. code-block:: python

          spy.build_sounding(
              data,
              radar=None,
              map_zoom=3,
              map_layers=[plot_cape, plot_mslp],
              map_label='RAP SBCAPE & MSLP',
          )

     - .. image:: /_static/gallery/custom_map/metar_dewpoints.png
          :alt: SounderPy sounding with METAR dewpoints on the map inset
          :width: 100%

       **METAR Dewpoints**

       Surface observations from the Iowa Environmental Mesonet.

       .. code-block:: python

          spy.build_sounding(
              data,
              radar=None,
              map_zoom=3,
              map_layers=plot_metars,
              map_label='METAR DEWPOINTS (F)',
          )

   * - .. image:: /_static/gallery/custom_map/everything.png
          :alt: SounderPy sounding with radar, RAP CAPE, isobars, and METAR dewpoints on the map inset
          :width: 100%

       **Everything**

       Radar, CAPE, isobars, and METARs, all on one map.

       .. code-block:: python

          spy.build_sounding(
              data,
              radar='single',
              map_zoom=3,
              map_layers=[plot_cape, plot_mslp,
                          plot_metars],
          )

     - .. image:: /_static/gallery/custom_figures/return_fig_sounding.png
          :alt: SounderPy sounding with a shaded, labeled EML drawn on afterward
          :width: 100%

       **Annotate the Figure**

       ``return_fig=True`` hands back the figure and its panels.

       .. code-block:: python

          fig, axes = spy.build_sounding(
              data, return_fig=True)
          axes['skewt'].ax.axhspan(850, 700,
              color='gold', alpha=0.2, zorder=0)
          axes['skewt'].ax.text(-28, 775, 'EML',
              fontsize=16, weight='bold')


.. _gallery-data-sources:

Different Data Sources
======================

One of SounderPy's core design features is that vertical profiles from many
different sources are converted into the same ``clean_data`` structure. Once a
profile has been retrieved or converted, the same plotting functions can be
used regardless of its source.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/data_sources/data_raob.png
          :alt: SounderPy sounding using observed RAOB data
          :width: 100%

       **Observed RAOB**

       .. code-block:: python

          data = spy.get_obs_data(
              "OAX",
              "2014", "06", "16", "18",
          )

     - .. image:: /_static/gallery/data_sources/data_rap_ruc.png
          :alt: SounderPy sounding using RAP/RUC reanalysis data
          :width: 100%

       **RAP/RUC Reanalysis**

       .. code-block:: python

          data = spy.get_model_data("rap-ruc",
              [44.58, -100.82],
              "2024", "08", "28", "18",
              box_avg_size=0.25)

   * - .. image:: /_static/gallery/data_sources/data_bufkit.png
          :alt: SounderPy sounding using BUFKIT forecast data
          :width: 100%

       **BUFKIT Forecast**

       .. code-block:: python

          data = spy.get_bufkit_data(
              "gfs",
              "KMOP",
              6,
              "2023", "08", "05", "12",
          )

     - .. image:: /_static/gallery/data_sources/data_acars.png
          :alt: SounderPy sounding using ACARS aircraft data
          :width: 100%

       **ACARS Aircraft Profile**

       .. code-block:: python

          conn = spy.acars_data(
              "2024", "05", "21", "18",
          )

          profiles = conn.list_profiles()
          data = conn.get_profile(profiles[0])

The same approach also applies to ERA5, NCEP, WRF, CM1, and user-created custom
profiles.

See :doc:`Getting Data </reference/gettingdata>` and
:doc:`Custom Data Sources </reference/customdatasources>` for additional information.


.. _gallery-winter:

Winter Soundings
================

``mode='winter'`` draws a sounding built around precipitation type, snow
growth, and winter hazards. These examples add the latest surface observation
(``show_obs=True``) and single-site radar. See
:doc:`/tutorials/winter_soundings` for a tour of the figure.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/winter/winter_bufkit_fzdz.png
          :alt: Winter sounding of a BUFKIT HRRR forecast during freezing drizzle at Grand Forks, ND
          :width: 100%

       **BUFKIT Forecast: Freezing Drizzle**

       HRRR at Grand Forks, ND, 24 Jan 2024 12Z F01.

       .. code-block:: python

          data = spy.get_bufkit_data(
              'hrrr', 'kgfk', 1,
              '2024', '01', '24', '12')
          spy.build_sounding(data, mode='winter',
                             show_obs=True,
                             radar='single-site')

     - .. image:: /_static/gallery/winter/winter_raob_ice.png
          :alt: Winter sounding of the Detroit, MI RAOB during an ice storm
          :width: 100%

       **Observed Sounding: Ice Storm**

       Detroit, MI RAOB, 23 Feb 2023 00Z.

       .. code-block:: python

          data = spy.get_obs_data(
              'DTX', '2023', '02', '23', '00')
          spy.build_sounding(data, mode='winter',
                             show_obs=True,
                             radar='single-site')

   * - .. image:: /_static/gallery/winter/winter_raob_dark.png
          :alt: Dark-mode winter sounding of the Omaha, NE RAOB during a blizzard
          :width: 100%

       **Dark Mode: Blizzard**

       Omaha, NE RAOB, 12 Jan 2024 12Z.

       .. code-block:: python

          data = spy.get_obs_data(
              'OAX', '2024', '01', '12', '12')
          spy.build_sounding(data, mode='winter',
                             dark_mode=True,
                             show_obs=True,
                             radar='single-site')

     - .. image:: /_static/gallery/winter/winter_acars.png
          :alt: Winter sounding of an ACARS aircraft profile at Denver, CO
          :width: 100%

       **ACARS Profile: Arctic Air**

       Denver, CO aircraft profile, 15 Jan 2024 14:20Z.

       .. code-block:: python

          conn = spy.acars_data(
              '2024', '01', '15', '14')
          data = conn.get_profile('DEN_1420')
          spy.build_sounding(data, mode='winter',
                             show_obs=True,
                             radar='single-site')


.. _gallery-composites:

Composite Soundings
===================

Composite plots are useful for comparing temporal evolution, model forecasts,
observations versus reanalysis, or other groups of profiles.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/composites/composite_light.png
          :alt: Light-mode SounderPy composite sounding
          :width: 100%

       **Light Composite**

       A standard multi-profile comparison.

     - .. image:: /_static/gallery/composites/composite_dark.png
          :alt: Dark-mode SounderPy composite sounding
          :width: 100%

       **Dark Composite**

       With ``dark_mode=True``.

   * - .. image:: /_static/gallery/composites/composite_sounding.png
          :alt: SounderPy composite with only a skew-T
          :width: 100%

       **Skew-T Only**

       With ``style="sounding"``.

     - .. image:: /_static/gallery/composites/composite_hodograph.png
          :alt: SounderPy composite with only a hodograph
          :width: 100%

       **Hodograph Only**

       With ``style="hodograph"``.


.. _gallery-publication:

Publication Figures
===================

``spy.build_publication_sounding()`` (experimental) draws clean, journal-sized figures of one or more
profiles, sized for AMS and AGU journals. See :doc:`/tutorials/custom_figures`.

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/custom_figures/publication_quickstart.png
          :alt: Four-panel publication figure with a skew-T, hodograph, radar map, and parameter table
          :width: 100%

       **All Four Panels**

       A RAOB and an HRRR forecast with a skew-T, hodograph, radar map, and parameter table.

       .. code-block:: python

          spy.build_publication_sounding(
              [raob, ashley],
              labels=['00Z ABR RAOB',
                      'HRRR F03 Ashley, ND'],
              panels=('skewt', 'hodograph',
                      'map', 'parameters'),
              map_radar='single',
          )

     - .. image:: /_static/gallery/custom_figures/publication_radar.png
          :alt: Publication figure with a skew-T, hodograph, and radar map
          :width: 100%

       **Map and Radar**

       Single-site radar under the profiles' locations.

       .. code-block:: python

          spy.build_publication_sounding(
              [hrrr, rap],
              labels=['HRRR F03', 'RAP F03'],
              panels=('skewt', 'hodograph',
                      'map'),
              map_radar='single',
          )

   * - .. image:: /_static/gallery/custom_figures/publication_style.png
          :alt: Publication figure restyled with serif fonts and custom colors
          :width: 100%

       **Custom Style**

       Colors, fonts, line widths, and hodograph rings, all set with ``style``.

       .. code-block:: python

          spy.build_publication_sounding(
              [hrrr, rap],
              colors=['black', 'darkorange'],
              show_adiabats=True, parcel='mu',
              style={'font_family': 'serif',
                     'hodograph_rings': False},
          )

     - .. image:: /_static/gallery/custom_figures/publication_annotated.png
          :alt: Publication figure with a shaded EML and a labeled map location
          :width: 100%

       **Annotated**

       The returned panels are ready to draw on: a shaded EML and a place on the map.

       .. code-block:: python

          fig, axes = (
              spy.build_publication_sounding(
                  [hrrr, rap],
                  panels=('skewt', 'hodograph',
                          'map')))
          axes['skewt'].ax.axhspan(
              850, 700, color='gold', alpha=0.25)
          axes['map'].plot(
              -100.07, 45.73, marker='*',
              transform=ccrs.PlateCarree())


Continue Learning
=================

For detailed explanations and additional examples, see:

* :doc:`📈 Plotting Data </reference/plottingdata>`
* :doc:`🎨 Plot Settings </reference/plotsettings>`
* :doc:`📝 Publication Figures </reference/publicationfigures>`
* :doc:`🌎 Getting Data </reference/gettingdata>`
* :doc:`✏️ Custom Data Sources </reference/customdatasources>`
* :doc:`💻 Command Line Tool </reference/cli>`

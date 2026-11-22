.. _tutorial-custom-map-data:

====================
🧭 Custom Map Data
====================

The map inset on SounderPy's soundings and hodographs can show your own data,
such as surface observations, a model field, or anything else with a latitude
and longitude. You get the data, and SounderPy hands you the map to draw it on.

This tutorial covers:

* how ``map_layers`` and ``map_label`` work;
* plotting METAR observations on the map;
* plotting a gridded model field (RAP CAPE and pressure) on the map;
* tips for layering, labels, and other figures.

.. image:: /_static/gallery/custom_map/rap_cape_mslp.png
   :width: 100%

*The Aberdeen, SD RAOB at 00 UTC on 29 August 2024 (see the* :doc:`Mound City case study </case_studies/mound_city_20240828>`\ *), with RAP surface-based CAPE and sea-level pressure on the map.*

.. note::

   Custom map data works from Python only (not the command line tool), and needs
   the map, so leave ``map_zoom`` above ``0``.


***************************************************************

How It Works
============

``map_layers`` takes a function, or a list of functions. SounderPy calls each one
with the map's axes, ``ax``, and your function draws on it with any Matplotlib or
cartopy plotting call: ``contourf``, ``contour``, ``pcolormesh``, ``scatter``,
``text``, ``barbs``, and so on.

``map_label`` puts a short label in a box at the top of the map.

.. code-block:: python

   import cartopy.crs as ccrs

   def my_layer(ax):
       ax.contourf(lon, lat, values, transform=ccrs.PlateCarree(), zorder=3)

   spy.build_sounding(clean_data, radar=None, map_layers=my_layer, map_label='MY DATA')

Three things to know:

* **Pass** ``transform=ccrs.PlateCarree()``. The map is a cartopy map, and this
  tells it your data's x and y are longitude and latitude.
* **Use** ``zorder`` **to stack your data with the map.** Your layers are drawn
  after the radar and before the profile's location marker; each map part's
  ``zorder``, and which values put your data under or over the borders, are listed in
  :ref:`custom_map_data`. The examples below use ``zorder=3`` for filled fields, and
  ``7`` or ``8`` for lines, points, and text.
* **Use** ``radar=None`` to show only your data. Radar can stay on too; its time
  stays in the box at the bottom of the map.


***************************************************************

Example 1: METAR Observations
=============================

Start by fetching a sounding: here's Aberdeen, SD RAOB at 00 UTC on 29 August 2024

.. code-block:: python

   import numpy as np
   import pandas as pd
   import cartopy.crs as ccrs
   import sounderpy as spy

   clean_data = spy.get_obs_data('ABR', '2024', '08', '29', '00')


Next is your custom data:

The `Iowa Environmental Mesonet <https://mesonet.agron.iastate.edu/request/download.phtml>`_
(IEM) serves archived METARs as a plain CSV file, so ``pandas`` can read them
straight from a URL. This asks for the dewpoint (``dwpf``) at every ASOS station in
South Dakota, North Dakota, Nebraska, and Minnesota within 30 minutes of 00 UTC,
then keeps each station's report closest to 00 UTC:

.. code-block:: python

   url = ('https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?'
          'network=SD_ASOS&network=ND_ASOS&network=NE_ASOS&network=MN_ASOS'
          '&data=dwpf&sts=2024-08-28T23:30Z&ets=2024-08-29T00:30Z&tz=Etc/UTC'
          '&format=onlycomma&latlon=yes&missing=empty&report_type=3')

   metars = pd.read_csv(url, parse_dates=['valid']).dropna(subset=['dwpf'])
   metars['offset'] = (metars['valid'] - pd.Timestamp('2024-08-29 00:00')).abs()
   metars = metars.sort_values('offset').drop_duplicates('station')

The layer function draws each station as a dot colored by its dewpoint, with the
value above it. Two steps keep the map readable:

* ``ax.get_extent()`` gives the map's west, east, south, and north edges, so only
  stations inside the map are drawn (with a little room for their labels);
* MetPy's ``reduce_point_density`` drops stations that are too close together,
  so the labels don't overlap.

.. code-block:: python

   from metpy.calc import reduce_point_density

   def plot_metars(ax):
       # keep the stations inside the map, then thin them so labels don't overlap
       west, east, south, north = ax.get_extent(ccrs.PlateCarree())
       obs = metars[metars['lon'].between(west + 0.2, east - 0.2) & metars['lat'].between(south + 0.3, north - 0.3)]
       obs = obs[reduce_point_density(np.column_stack([obs['lon'], obs['lat']]), 0.35)]

       ax.scatter(obs['lon'], obs['lat'], c=obs['dwpf'], cmap='BrBG', vmin=40, vmax=75, s=60,
                  edgecolor='black', transform=ccrs.PlateCarree(), zorder=8)
       for _, ob in obs.iterrows():
           ax.text(ob['lon'], ob['lat'] + 0.12, f"{ob['dwpf']:.0f}", ha='center', va='bottom', fontsize=8,
                   weight='bold', transform=ccrs.PlateCarree(), zorder=8)

   spy.build_sounding(clean_data, radar=None, map_zoom=3, map_layers=plot_metars,
                      map_label='METAR DEWPOINTS (F), 00 UTC')

.. image:: /_static/gallery/custom_map/metar_dewpoints.png
   :width: 100%

.. tip::

   The same request works for other fields: add ``&data=tmpf`` (temperature),
   ``&data=sknt&data=drct`` (wind speed and direction), or ``&data=mslp``, and
   change the ``network=`` list to the states you need (``KS_ASOS``, ``IA_ASOS``,
   and so on).


***************************************************************

Example 2: A Gridded Model Field (RAP)
======================================

NOAA keeps RAP model runs on Amazon Web Services (AWS) as GRIB2 files, from 2021
on. The analysis (forecast hour 0) on the 13-km grid is one 18 MB file per hour;
this downloads the 00 UTC 29 August 2024 file:

.. code-block:: python

   import fsspec
   import xarray as xr

   fs = fsspec.filesystem('s3', anon=True)
   fs.get('s3://noaa-rap-pds/rap.20240829/rap.t00z.awp130pgrbf00.grib2', 'rap.grib2')

``xarray`` opens GRIB2 files with the ``cfgrib`` engine. A GRIB2 file holds many
fields on many kinds of levels, so ``filter_by_keys`` picks out the ones you
want: here, surface CAPE (``'shortName': 'cape'`` on the ``'surface'`` level) and
sea-level pressure (the only field on the ``'meanSea'`` level):

.. code-block:: python

   def open_rap(keys):
       return xr.open_dataset('rap.grib2', engine='cfgrib',
                              backend_kwargs={'filter_by_keys': keys, 'indexpath': ''})

   cape = open_rap({'typeOfLevel': 'surface', 'stepType': 'instant', 'shortName': 'cape'})['cape']
   mslp = open_rap({'typeOfLevel': 'meanSea'})['mslma'] / 100   # Pa to hPa
   lat, lon = cape.latitude, cape.longitude

Then one layer fills the CAPE and a second draws the isobars on top of it:

.. code-block:: python

   def plot_cape(ax):
       ax.contourf(lon, lat, cape, levels=np.arange(250, 5001, 250), cmap='magma_r', alpha=0.6,
                   transform=ccrs.PlateCarree(), zorder=3)

   def plot_mslp(ax):
       ax.contour(lon, lat, mslp, levels=np.arange(980, 1041, 2), colors='black', linewidths=1,
                  transform=ccrs.PlateCarree(), zorder=7)

   spy.build_sounding(clean_data, radar=None, map_zoom=3, map_layers=[plot_cape, plot_mslp],
                      map_label='RAP SBCAPE (J/kg) & MSLP (2 hPa)')

.. image:: /_static/gallery/custom_map/rap_cape_mslp.png
   :width: 100%

.. tip::

   To see what else a GRIB2 file holds, open it with a level and print the
   result, e.g. ``print(open_rap({'typeOfLevel': 'heightAboveGround', 'level': 2}))``
   lists the 2-m temperature, dewpoint, and humidity fields.


***************************************************************

Tips
====

* **Several layers** are drawn in list order, so put fills first and lines,
  points, and text last.
* **Any figure with a map** takes ``map_layers`` and ``map_label``:
  ``build_sounding()`` (including ``mode='winter'``), ``build_hodograph()``, and
  ``build_vad_hodograph()``.
* **Zoom out** with ``map_zoom`` to show more of a field. ``ax.get_extent()``
  always gives the edges of the map as shown.
* **Data in another projection** (for example, a Lambert conformal grid in x/y
  meters) needs that projection's ``transform``, e.g.
  ``transform=ccrs.LambertConformal(...)``, instead of ``ccrs.PlateCarree()``.


***************************************************************

And I guess if you really wanted to, you could add.....everything:
==================================================================

.. code-block:: python

   spy.build_sounding(clean_data, radar='single', map_zoom=3,
                      map_layers=[plot_cape, plot_mslp, plot_metars],
                      map_label='RAP SBCAPE & MSLP + METAR Td + RADAR')

.. image:: /_static/gallery/custom_map/everything.png
   :width: 100%

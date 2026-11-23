.. _tutorial-profile-tools:

===============================
🧰 Profile Tools & Parameters
===============================

Once a profile is in SounderPy, it's more than a figure: you can pull out its
parameters, look up values at any height, follow them through a forecast, and
change the profile itself by modifying its surface, smoothing it, or merging it
with another.

This tutorial covers:

* calculating parameters with ``sounding_params()``;
* the interpolated profile and looking up values by height;
* ``interp_data()`` and ``find_nearest()``;
* parameters through every hour of a forecast;
* modifying the surface with ``modify_sfc`` and ``modify_surface()``;
* smoothing a profile with ``smooth_profile()``;
* merging two profiles with ``merge_profiles()``;
* thinning a profile and zooming the skew-T.

Every function is listed in :doc:`Helper Tools </reference/helpertools>`.


***************************************************************

Retrieve the Example Profile
============================

Most examples use the OAX sounding at 18Z on 16 June 2014, the afternoon of the
Pilger, NE tornadoes:

.. code-block:: python

   import copy

   import sounderpy as spy

   data = spy.get_obs_data(
       "OAX",
       "2014", "06", "16", "18",
       hush=True,
   )


***************************************************************

Calculate Parameters
====================

``sounding_params()`` calculates every value drawn on a SounderPy sounding and
returns them in a dictionary:

.. code-block:: python

   params = spy.sounding_params(data).calc()

   params.keys()
   # dict_keys(['general', 'thermo', 'kinem', 'intrp'])

* ``general``: surface values, and values at each level such as wet-bulb
  temperature, relative humidity, and theta-e;
* ``thermo``: parcels (CAPE, CIN, LCL, LFC, EL), ECAPE, lapse rates, and
  thermodynamic indices;
* ``kinem``: storm motions, bulk shear, storm-relative helicity and wind, and
  composite parameters such as STP and SCP;
* ``intrp``: the profile interpolated to 100 m levels (see below).

Pick out a value by its key:

.. code-block:: python

   thermo, kinem = params["thermo"], params["kinem"]

   print(f"MLCAPE:        {thermo['mlcape']:.0f} J/kg")
   print(f"0-1 km SRH:    {kinem['srh_0_to_1000']:.0f} m2/s2")
   print(f"0-6 km shear:  {kinem['shear_0_to_6000']:.0f} kt")
   print(f"Effective STP: {kinem['eil_stp']:.1f}")

.. code-block:: text

   MLCAPE:        4311 J/kg
   0-1 km SRH:    357 m2/s2
   0-6 km shear:  56 kt
   Effective STP: 9.5

Every key, its units, and its source are listed in
:ref:`Parameters & References <params-severe>`.

``sounding_params()`` takes the same ``storm_motion`` and ``modify_sfc``
settings as the figures, so the numbers match the plot you draw with them:

.. code-block:: python

   params_lm = spy.sounding_params(data, storm_motion="left_moving").calc()

.. tip::
   ``spy.winter_params(data).calc()`` returns all of the above plus a
   ``winter`` dictionary: precipitation type, snow growth, and winter hazards.
   See "Winter values in code" in the :doc:`Winter Soundings </tutorials/winter_soundings>` tutorial.

To print a quick summary of the main values instead:

.. code-block:: python

   spy.print_variables(data)


***************************************************************

The Interpolated Profile
========================

Observed and model levels fall at different heights in every profile.
``params["intrp"]`` holds the profile interpolated to every 100 m above ground,
so the same index means the same height in any profile:

.. code-block:: python

   intrp = params["intrp"]

   intrp["zINTRP"][:5]
   # array([  0., 100., 200., 300., 400.])

It holds temperature (``tINTRP``), wet-bulb temperature (``wbINTRP``), relative
humidity (``rhINTRP``), mixing ratio (``mrINTRP``), theta and theta-e
(``thetaINTRP``, ``thetaeINTRP``), pressure (``pINTRP``), and wind
(``uINTRP``, ``vINTRP``).

``intrp["hgt_lvls"]`` gives the index of common heights: ``"h0"`` is the
surface, ``"h05"`` is 500 m, ``"h1"`` is 1 km, and so on to ``"h14"``:

.. code-block:: python

   lvl = intrp["hgt_lvls"]

   t_1km = intrp["tINTRP"][lvl["h1"]]    # temperature 1 km above ground (°C)
   u_3km = intrp["uINTRP"][lvl["h3"]]    # u-wind 3 km above ground (kt)


***************************************************************

Interpolate and Look Up Values Yourself
=======================================

``interp_data()`` interpolates any array to evenly spaced heights. The heights
must be **above ground**, starting at 0, so subtract the surface height first:

.. code-block:: python

   z_agl = data["z"] - data["z"][0]

   # theta-e every 250 m above ground
   theta_e = params["general"]["theta_e"]
   theta_e_250m = spy.interp_data(theta_e, z_agl, step=250)

``find_nearest()`` returns the index of the value closest to the one you ask
for, which finds the observed level nearest a height or pressure:

.. code-block:: python

   i = spy.find_nearest(z_agl.m, 3000)

   print(z_agl[i], data["p"][i], data["T"][i])
   # 2985.0 meter 682.0 hectopascal 10.4 degree_Celsius


***************************************************************

Parameters Through a Forecast
=============================

``get_bufkit_data()`` with ``all_fcst_hours=True`` downloads every forecast
hour of a model run at once, as a list of profiles. Calculate parameters for
each to see how they change.

This is the 12Z HRRR run for Mobridge, SD on 28 August 2024, the day of the
Mound City tornadoes (see the :doc:`case study </case_studies/mound_city_20240828>`),
from 12Z to 00Z:

.. code-block:: python

   import matplotlib.pyplot as plt

   forecast = spy.get_bufkit_data(
       "hrrr", "KMBG", 0,
       "2024", "08", "28", "12",
       all_fcst_hours=True,
       hush=True,
   )

   hours, mlcape, srh = [], [], []
   for profile in forecast[:13]:                         # F00 to F12
       params = spy.sounding_params(profile).calc()
       hours.append(profile["site_info"]["fcst-hour"])   # 'F00', 'F01', ...
       mlcape.append(params["thermo"]["mlcape"])
       srh.append(params["kinem"]["srh_0_to_1000"])

   fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
   ax1.plot(hours, mlcape, marker="o")
   ax1.set_title("MLCAPE (J/kg)")
   ax2.plot(hours, srh, marker="o")
   ax2.set_title("0-1 km SRH (m²/s²)")
   plt.show()

.. image:: /_static/gallery/profile_tools/params_trend.png
   :alt: MLCAPE and 0-1 km storm-relative helicity at each forecast hour of a HRRR run
   :width: 82%
   :align: center

The instability builds through the afternoon while the low-level helicity
dips, then recovers toward evening. ``fcst_hour`` (the ``0``) is ignored with
``all_fcst_hours=True``. A recent run works the same way without the date:
``spy.get_bufkit_data("hrrr", "KMBG", 0, all_fcst_hours=True)``.


***************************************************************

Modify the Surface
==================

A morning sounding often has a surface inversion that the afternoon will mix
out. ``modify_sfc`` gives a profile new surface values for a figure or a
calculation. This is the 12Z OAX sounding with its surface warmed to the 18Z
observation:

.. code-block:: python

   morning = spy.get_obs_data(
       "OAX",
       "2014", "06", "16", "12",
       hush=True,
   )

   spy.build_sounding(morning, ylim=(1000, 600))

   spy.build_sounding(
       morning,
       modify_sfc={"T": 28, "Td": 24},
       ylim=(1000, 600),
   )

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/surface_modification/sfc_none.png
          :alt: The 12Z OAX sounding with its surface inversion
          :width: 100%

       **As observed**

     - .. image:: /_static/gallery/surface_modification/sfc_mix.png
          :alt: The 12Z OAX sounding with its surface warmed and mixed up its dry adiabat
          :width: 100%

       **Warmed to 28 °C, dewpoint 24 °C**

``T`` and ``Td`` are in °C, ``ws`` in knots, and ``wd`` in degrees; leave out
any you don't want to change. The ``method`` decides how far up the change
reaches:

* ``"mix"`` (default): a warmer surface mixes up its dry adiabat until it meets
  the profile;
* ``"blend"``: the change tapers to zero at ``depth`` m above ground (500 by
  default);
* ``"replace"``: only the surface level changes.

See :ref:`Surface modification <modify_sfc>` for how each works, and the
:ref:`gallery <gallery-surface>` for all three on this sounding.

The same ``modify_sfc`` works in ``sounding_params()``:

.. code-block:: python

   before = spy.sounding_params(morning).calc()
   after = spy.sounding_params(morning, modify_sfc={"T": 28, "Td": 24}).calc()

   print(before["thermo"]["mlcape"], after["thermo"]["mlcape"])
   # about 60 and 3,300 J/kg

``modify_sfc`` changes the profile only for that figure or calculation.
``modify_surface()`` returns the modified profile as a new ``clean_data``
dictionary, leaving the original as it is, so you can save it, compare it, or
use it like any other profile:

.. code-block:: python

   heated = spy.modify_surface(morning, {"T": 28, "Td": 24})

   spy.to_file("csv", heated, filename="oax_12z_heated.csv")

The :doc:`CM1 surface heating case study </case_studies/cm1_surface_heating>`
uses ``modify_surface()`` to make new CM1 input soundings.


***************************************************************

Smooth a Profile
================

``smooth_profile()`` smooths the pressure, height, temperature, dewpoint, and
wind arrays. It changes the profile you give it, so pass a copy to keep the
original:

.. code-block:: python

   # Savitzky-Golay (the default): a wider window smooths more
   smoothed = spy.smooth_profile(copy.deepcopy(data), method="savgol", window_length=15)

   # Gaussian: a larger sigma smooths more
   smoothed = spy.smooth_profile(copy.deepcopy(data), sigma=3)

   spy.build_sounding(smoothed)

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/profile_tools/smooth_savgol.png
          :alt: The OAX sounding smoothed with a Savitzky-Golay filter
          :width: 100%

       **Savitzky-Golay**, ``window_length=15``

     - .. image:: /_static/gallery/profile_tools/smooth_gaussian.png
          :alt: The OAX sounding smoothed with a Gaussian filter
          :width: 100%

       **Gaussian**, ``sigma=3``

The Savitzky-Golay filter keeps peaks and inversions sharper than the Gaussian
filter at the same amount of smoothing.

.. note::
   Smoothing changes the data. Parameters that depend on sharp vertical
   changes, such as lapse rates, inversions, shear, and parcel properties, will
   differ from the original profile's.


***************************************************************

Merge Two Profiles
==================

``merge_profiles()`` averages two profiles at common heights, weighting the
first by ``weight_a``. It can blend two soundings in time, or two models. This
is 60% of the 18Z OAX sounding and 40% of the 00Z one:

.. code-block:: python

   later = spy.get_obs_data(
       "OAX",
       "2014", "06", "17", "00",
       hush=True,
   )

   merged = spy.merge_profiles(data, later, weight_a=0.6)

The merged profile holds only the data arrays (``p``, ``z``, ``T``, ``Td``,
``u``, ``v``). Add ``site_info`` and plot titles before plotting, here copied
from the 18Z sounding:

.. code-block:: python

   merged["site_info"] = dict(data["site_info"], source="MERGED RAOB PROFILES")
   merged["titles"] = {
       "top_title": "MERGED VERTICAL PROFILE",
       "left_title": "60% 18Z OAX 06/16/2014 | 40% 00Z OAX 06/17/2014",
       "right_title": data["titles"]["right_title"],
   }

   spy.build_sounding(merged)

.. image:: /_static/gallery/profile_tools/merged.png
   :alt: A SounderPy sounding of two merged OAX profiles
   :width: 82%
   :align: center

See :ref:`The clean_data scheme <clean_data_scheme>` for the ``site_info`` keys.


***************************************************************

Fewer Levels and a Zoomed Skew-T
================================

Two plot settings change how much of the profile is used and shown, without
changing the data:

* ``max_levels`` thins a profile to that many evenly spaced levels (250 by
  default), which speeds up dense profiles such as high-resolution radiosonde
  data. ``max_levels=None`` uses every level.
* ``ylim`` zooms the skew-T to a pressure range, (bottom, top) in hPa. The
  rest of the figure still uses the full profile.

.. code-block:: python

   spy.build_sounding(data, max_levels=30)

   spy.build_sounding(data, ylim=(1000, 500))

.. list-table::
   :widths: 50 50
   :class: gallery-table

   * - .. image:: /_static/gallery/profile_tools/max_levels.png
          :alt: The OAX sounding thinned to 30 levels
          :width: 100%

       ``max_levels=30``

     - .. image:: /_static/gallery/profile_tools/zoom_ylim.png
          :alt: The OAX sounding with the skew-T zoomed to 1000-500 hPa
          :width: 100%

       ``ylim=(1000, 500)``


***************************************************************

Next Steps
==========

* :doc:`Helper Tools </reference/helpertools>`: every function in this tutorial
* :ref:`Parameters & References <params-severe>`: every parameter key, its
  units, and its source
* :doc:`Exporting Data </tutorials/exporting_data>`: save modified, smoothed,
  or merged profiles
* :doc:`Composite Soundings </tutorials/composite_soundings>`: compare profiles
  on one figure instead of merging them
* :doc:`Publication & Custom Figures </tutorials/custom_figures>`: compare an
  original and a modified profile side by side

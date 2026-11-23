.. _home-test:

=============
🎈 SounderPy
=============

.. container:: home-hero

   .. container:: home-hero-copy

      .. raw:: html

         <p class="home-kicker">VERTICAL PROFILE DATA RETRIEVAL &amp; ANALYSIS FOR PYTHON</p>
         <p class="home-lead">
           Retrieve, analyze, visualize, and export atmospheric vertical profiles
           from observations, forecasts, reanalyses, and numerical models.
         </p>

      .. container:: home-actions

         :doc:`🚀 Get Started </user_guide/getting_started>`   
         :doc:`👨‍🍳 Tutorials </tutorials/index>`   
         :doc:`📊 Plot Gallery </gallery/examplegallery>`   
         :doc:`📚 Reference </reference/gettingdata>`

      .. code-block:: bash

         pip install sounderpy

      .. code-block:: bash

         conda install conda-forge::sounderpy

   .. container:: home-hero-figure

      .. image:: /_static/gallery/map_radar/radar_single.png
         :alt: Example SounderPy vertical profile analysis figure
         :width: 100%

      .. image:: /_static/gallery/themes/hodograph_light.png
         :alt: Example SounderPy vertical profile analysis figure
         :width: 100%

Start Here
==========

.. container:: home-card-grid

   .. container:: home-card

      **🚀 Getting Started**

      Install SounderPy, retrieve a profile, and create your first sounding.

      :doc:`Start learning → </user_guide/getting_started>`

   .. container:: home-card

      **👨‍🍳 Tutorials**

      Step-by-step guides for plots, radar, parcels, hodographs, data sources,
      exports, and command-line workflows.

      :doc:`Browse tutorials → </tutorials/index>`

   .. container:: home-card

      **📊 Plot Gallery**

      Explore sounding, hodograph, composite, radar, parcel, and data-source
      examples.

      :doc:`View gallery → </gallery/examplegallery>`

   .. container:: home-card

      **🌎 Data Sources**

      Work with RAOB, ACARS, RAP/RUC, BUFKIT, ERA5, NCEP, WRF, CM1, and custom
      vertical profiles.

      :doc:`Explore data sources → </reference/gettingdata>`

   .. container:: home-card

      **📚 Reference**

      Find detailed retrieval, plotting, calculation, and helper-function
      documentation.

      :doc:`Open reference → </reference/plottingdata>`

   .. container:: home-card

      **💻 Command Line Interface**

      Retrieve, plot, export, and automate SounderPy directly from a terminal.

      :doc:`Use the CLI → </reference/cli>`


Get a Sounding in Seconds
=========================

SounderPy provides the same core workflow through both Python and the command
line.

.. container:: home-code-grid

   .. container:: home-code-panel

      **Python**

      .. code-block:: python

         import sounderpy as spy

         data = spy.get_obs_data(
             "OAX",
             "2014", "06", "16", "18")

         spy.build_sounding(data)

   .. container:: home-code-panel

      **Command Line**

      .. code-block:: bash

         sounderpy obs OAX 2014-06-16 18 \
             --plot sounding




One Workflow, Many Data Sources, Several Plotting Options
==========================================================

.. container:: source-strip

   .. container:: source-pill

      **Observations**

      RAOB, IGRAv2, ACARS

   .. container:: source-pill

      **Forecasts**

      BUFKIT: HRRR, RAP, GFS, etc

   .. container:: source-pill

      **Reanalysis**

      RAP, RUC, ERA5, NCEP-FNL

   .. container:: source-pill

      **Modeling Systems**

      CM1, WRF

   .. container:: source-pill

      **Custom Sources**

      SHARPPY + more

   .. container:: source-pill

      **Soundings**

      Full Skew-T, Log-P diagram and hodograph figure.

   .. container:: source-pill

      **Hodographs**

      Simple hodograph figure.

   .. container:: source-pill

      **Composite Soundings**

      Several profiles on one sounding.


:doc:`Learn how SounderPy handles data → </tutorials/data_sources>`


What's New in SounderPy 3.3.0
=============================

.. container:: release-panel

   **SounderPy 3.3.0 rebuilds the plotting module, adds winter and publication figures, and makes every
   figure easier to customize.**

   .. rubric:: 📈 Figure Updates

   * **Rebuilt plot module**: most figures are about 50% faster, with custom colors (``theme``), a custom watermark (``spy.set_watermark()``), and a ``ylim`` to zoom the skew-T.
   * **New** ``build_sounding()`` **kwargs**: ``max_levels`` (reduce the number of vertical levels of a profile), ``ylim`` (zoom the skew-T into a particular layer), ``theme`` (custom figure colors), ``sr_hodo`` (storm-relative hodograph), ``mode`` (``'severe'``, ``'winter'``, or ``'auto'``), ``show_obs`` (the nearest surface observation for ``mode='winter'`` figures), ``map_layers`` and ``map_label`` (adding your own data to the map inset), and ``return_fig`` (return the figure and its panels for customization). Many of these also work on ``build_hodograph()``, ``build_composite()``, and ``build_vad_hodograph()``. See :doc:`/reference/plottingdata`.
   * **Publication figures (experimental)**: ``spy.build_publication_sounding()`` draws journal-sized figures for AMS and AGU journals, and every ``build_*`` function takes ``return_fig=True`` to hand back its figure (see :doc:`/tutorials/custom_figures`).
   * **Your own data on the map inset** with ``map_layers`` and ``map_label`` (see :doc:`/tutorials/custom_map_data`).
   * **Rebuilt composite soundings** with figure improvements and parameter read-outs, and rebuilt VAD hodographs with map, radar, and boundary options.
   * **New parameter readouts** on the full sounding: composite indices in the thermodynamics box, EIL depth and heights, 0-10 km bulk shear and theta-e difference, PBL height and transport wind, improved layout for ML/MU parcel heights, updated watermark, new data source attribution line, and numerous other small improvements to the plot layout.

   .. rubric:: ❄️ NEW Winter Soundings

   * ``build_sounding(mode='winter')`` draws a winter-focused sounding built around precipitation type, snow growth, and winter hazards. Setting ``mode='auto'`` will generate a winter figure when the surface wet-bulb is below 35°F, and a new ``show_obs=True`` kwarg adds the latest (nearest) surface observation to the figure (See the :doc:`winter soundings tutorial </tutorials/winter_soundings>`!).
   * ``spy.winter_params(clean_data).calc()`` returns the dictionary of winter parameters.

   .. rubric:: 🧮 Calculations

   * NEW severe parameters: 0-10 km bulk shear and the low-to-midlevel theta-e difference (for MCSs; Cohen et al. 2007), PBL depth and transport wind (for air quality), EIL bottom height, top height, and depth, SHIP, MLBRN, TT, KI, and 700-500 hPa lapse rate.
   * UPDATED/CORRECTED parameters: Streamwise vorticity bug counted antistreamwise layers as streamwise; STP/SCP used EIL shear in kts (now m/s).
   * Every parameter, its units, and its source are listed on the :doc:`Parameters & References </reference/references>` documentation page.
   * Settings: ``modify_sfc`` adds ``'mix'``, ``'blend'``, and ``'replace'`` methods (see :ref:`modify_sfc`), and the NEW ``spy.modify_surface()`` returns the modified profile itself, to compare or save (see the :doc:`CM1 surface heating case study </case_studies/cm1_surface_heating>`); ``spy.smooth_profile()`` adds a Savitzky-Golay filter.

   .. rubric:: 🌦️ Data Sources

   * BUFKIT now serves the RRFS; NAM, NAMNEST, SREF, and HIRESW are no longer available as "latest runs" (discontinued; but NAM and NAMNEST are available as archived runs).
   * RRFS, HRRR, and RAP may now fall back to NOMADS NCEP's BUFR dataset when Penn State is behind or an archived run is missing.
   * NEW Station list: observed soundings use a packaged WMO/ICAO/IGRAv2 station list, with an ICAO fallback for the Wyoming archive. This list is a best-guess combination of 5 data sources to limit incorrect site names, IDs, and locations.
   * Dewpoints are now calculated (when possible) for older IGRA v2 soundings that don't ship dewpoint.
   * Thread-safe RAP/RUC GRIB parsing.

   .. rubric:: 💻 Command Line Interface Updates

   * NEW ``file`` command: read a local SHARPpy-format sounding (``--latlon LAT LON`` adds its location for the map).
   * NEW options:

     * ``--params``: print key calculated parameters (CAPE, SRH, shear, and precip type with ``--mode winter``);
     * ``--mode {severe,winter,auto}`` and ``--obs``: winter soundings and the nearest surface observation;
     * ``--storm-motion``: ``right_moving``, ``left_moving``, ``mean_wind``, or a direction and speed such as ``270 25``;
     * ``--modify-sfc``: change the surface, e.g. ``T=30 Td=21`` or ``T=-6 method=blend depth=300``;
     * ``--radar {mosaic,single,none}`` and ``--radar-time {sounding,now}``: the map's radar;
     * ``--parcels`` and ``--background-parcels``: custom parcel codes.

   * Options now match the Python functions, with clearer warnings. See :doc:`/reference/cli`.

   .. rubric:: 🐞 Bug Fixes

   * ACARS airport elevations corrected, CIN doesn't accumulate above the LFC, HGZ identification, extra moist adiabats on thin zoomed skew-Ts, and figures now close after ``save=True`` to preserve memory.

   .. rubric:: 📚 NEW HOME for Documentation!

   * Freshly redesigned docs now live at `sounderpy.app/docs <https://sounderpy.app/docs/>`_, with new :doc:`/reference/plotsettings`, :doc:`/reference/publicationfigures`, and :doc:`/reference/references` pages.
   * Also, each source module now includes a text header with its docs and references.

   .. important::
      **Changes to note:**

      * ``sounding_params(...).calc()`` and ``vad_params(...).calc()`` return a dictionary instead of a tuple.
      * Renamed, with the old names still working (with a ``FutureWarning``) until a later release: ``build_composite()``'s ``colors_to_use``, ``ls_to_use``, ``alphas_to_use``, and ``lw_to_use`` (now ``colors``, ``linestyles``, ``alphas``, and ``linewidths``), ``make_cm1_profile()`` / ``make_wrf_profile()`` (now ``spy.from_cm1()`` / ``spy.from_wrf()``), and the data modules (now in ``sounderpy.data``).
      * ``spy.smooth_profile()`` defaults to the Savitzky-Golay filter (passing ``sigma`` alone still uses the Gaussian filter).

      See :doc:`/reference/helpertools` and :doc:`/reference/plottingdata`.

   :doc:`See the plotting reference → </reference/plottingdata>`


Using SounderPy in Research?
============================

SounderPy is open-source scientific software developed for atmospheric vertical
profile retrieval, analysis, and visualization.

.. container:: research-links

   :doc:`About SounderPy </user_guide/about>`   
   `GitHub <https://github.com/kylejgillett/sounderpy>`_   
   `PyPI <https://pypi.org/project/sounderpy/>`_   
   `Operational SounderPy Site <https://sounderpy.app/>`_


☕ SounderPy is an open-source package developed on my own time. Would you like to support continued SounderPy development? Consider "`Buying me a coffee <https://www.buymeacoffee.com/kylejgillett>`_"! ☕

Directory
==================
.. toctree::
   :maxdepth: 5
   :caption: General:
   
   user_guide/getting_started
   user_guide/about
   user_guide/troubleshooting
   user_guide/community
   user_guide/additional_resources

.. toctree::
   :maxdepth: 5
   :caption: API Reference:
 
   reference/gettingdata
   reference/customdatasources
   reference/clean_data_scheme
   reference/plottingdata
   reference/plotsettings
   reference/publicationfigures
   reference/helpertools
   reference/cli
   reference/references

.. toctree::
   :maxdepth: 5
   :caption: Tutorials:

   tutorials/beginners_sounderpy_cookbook
   tutorials/index

.. toctree::
   :maxdepth: 5
   :caption: Examples:

   gallery/examplegallery
   case_studies/mound_city_20240828
   case_studies/cm1_surface_heating

* :ref:`genindex`
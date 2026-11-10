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

         conda install sounderpy

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

      RAOB, IRGAv2, ACARS

   .. container:: source-pill

      **Forecasts**

      BUFKIT: HRRR, RAP, GFS, etc

   .. container:: source-pill

      **Reanalysis**

      RAP, RUC, ERA5, NCEP-FNL

   .. container:: source-pill

      **Modeling Systems**

      CM1, WRF, MPAS

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

   **SounderPy 3.3.0 rebuilds the plotting module and makes every figure easier
   to customize.**

   * 📈 SounderPy's plot module is fully rebuilt and refactored allowing for increased figure customization and generation speed. Most figures are roughly 50% faster.
   * 🎨 The new figure kwarg ``theme`` lets you set your own colors while credit text and logo can be changed with ``spy.set_watermark()``.
   * 🔍 Full sounding and composite figures now include a ``ylim`` kwarg, allowing you to adjust the Skew-T's y-axis (pressure) limits.
   * 🧩 Composite soundings have been completely rebuilt. Composites include a table of various parameters, normalized kwarg names, and support for custom data sources.
   * 🗺️ VAD hodographs now include the map, radar, hodograph boundary, and DPI options of the other figures.
   * ❄️ New winter sounding: ``build_sounding(mode='winter')`` draws a sounding built around precipitation type, snow growth, and winter hazards. ``mode='auto'`` picks the winter sounding when the surface wet-bulb is below 35 degrees F, and ``show_obs=True`` adds the latest surface observation. See :doc:`/tutorials/winter_soundings`.
   * 🧮 ``spy.winter_params(clean_data).calc()`` returns every ``sounding_params`` value plus the winter values behind the winter sounding.
   * 💻 The command line tool adds ``--mode {severe,winter,auto}`` and ``--obs`` for sounding plots.
   * ✈️ ACARS profiles now use the airport elevation in meters (it was read in feet), which corrects the surface label and surface pressure.
   * 🌦️ The NAM and NAMNEST (NAM3km) BUFKIT forecasts are archive only after their NCEP decommissioning on Nov. 2nd, 2026.
   * 📚 Refactored and redesigned documentation to match new SounderPy web app: sounderpy.app

   .. important::
      ``sounding_params(...).calc()`` and ``vad_params(...).calc()`` now return a dictionary
      instead of a tuple, and ``build_composite()``'s ``colors_to_use``, ``ls_to_use``,
      ``alphas_to_use``, and ``lw_to_use`` are renamed to ``colors``, ``linestyles``,
      ``alphas``, and ``linewidths``. The old names still work in 3.3.0 but show a
      ``FutureWarning``, and will be removed in a later release. The data retrieval and file modules
      (such as ``sounderpy.bufkit_data``) moved into ``sounderpy.data`` (``sounderpy.data.bufkit_data``);
      ``spy.*`` functions are unchanged, and the old module paths still work with a ``FutureWarning``.
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


☕ SounderPy is a open-source package developed on my own time. Would you like to support continued SounderPy development? Consider "`Buying me a coffee <https://www.buymeacoffee.com/kylejgillett>`_"! ☕

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
   reference/helpertools
   reference/cli

.. toctree::
   :maxdepth: 5
   :caption: Tutorials:  

   tutorials/beginners_sounderpy_cookbook
   gallery/examplegallery
   tutorials/index
   case_studies/mound_city_20240828

* :ref:`genindex`
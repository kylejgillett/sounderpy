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

         :doc:`🚀 Get Started <getting_started>`   
         :doc:`👨‍🍳 Tutorials <tutorials/index>`   
         :doc:`📊 Plot Gallery <examplegallery>`   
         :doc:`📚 Reference <gettingdata>`

      .. code-block:: bash

         pip install sounderpy

      .. code-block:: bash

         conda install sounderpy

   .. container:: home-hero-figure

      .. image:: _static/gallery/map_radar/radar_single.png
         :alt: Example SounderPy vertical profile analysis figure
         :width: 100%

      .. image:: _static/gallery/themes/hodograph_light.png
         :alt: Example SounderPy vertical profile analysis figure
         :width: 100%

Start Here
==========

.. container:: home-card-grid

   .. container:: home-card

      **🚀 Getting Started**

      Install SounderPy, retrieve a profile, and create your first sounding.

      :doc:`Start learning → <getting_started>`

   .. container:: home-card

      **👨‍🍳 Tutorials**

      Step-by-step guides for plots, radar, parcels, hodographs, data sources,
      exports, and command-line workflows.

      :doc:`Browse tutorials → <tutorials/index>`

   .. container:: home-card

      **📊 Plot Gallery**

      Explore sounding, hodograph, composite, radar, parcel, and data-source
      examples.

      :doc:`View gallery → <examplegallery>`

   .. container:: home-card

      **🌎 Data Sources**

      Work with RAOB, ACARS, RAP/RUC, BUFKIT, ERA5, NCEP, WRF, CM1, and custom
      vertical profiles.

      :doc:`Explore data sources → <gettingdata>`

   .. container:: home-card

      **📚 Reference**

      Find detailed retrieval, plotting, calculation, and helper-function
      documentation.

      :doc:`Open reference → <plottingdata>`

   .. container:: home-card

      **💻 Command Line Interface**

      Retrieve, plot, export, and automate SounderPy directly from a terminal.

      :doc:`Use the CLI → <cli>`


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


:doc:`Learn how SounderPy handles data → <tutorials/data_sources>`


What's New in SounderPy 3.2.0
===========================

.. container:: release-panel

   **SounderPy 3.2.0 expands both how data are retrieved and how SounderPy can
   be used.**

   * 💻 New command line interface for retrieval, plotting, JSON, and export.
   * 🌦️ Updated RAP/RUC reanalysis data access with modern GRIB2 archives.
   * 📈 New and improved sounding plot features (such as storm-relative wind barbs and DPI settings)
   * 🎈 Restored access to UW RAOB database.
   * 📁 New file-io workflow for saving data to a file and for opening SHARPPY files.
   * 📚 Expanded task-oriented tutorials and visual plot gallery.
   * 🧰 Improved packaging, testing, and release infrastructure.

   :doc:`Try the CLI → <cli>`


Using SounderPy in Research?
============================

SounderPy is open-source scientific software developed for atmospheric vertical
profile retrieval, analysis, and visualization.

.. container:: research-links

   :doc:`About SounderPy <about>`   
   `GitHub <https://github.com/kylejgillett/sounderpy>`_   
   `PyPI <https://pypi.org/project/sounderpy/>`_   
   `Operational SounderPy Site <https://sounderpysoundings.anvil.app/>`_


☕ SounderPy is a open-source package developed on my own time. Would you like to support continued SounderPy development? Consider "`Buying me a coffee <https://www.buymeacoffee.com/kylejgillett>`_"! ☕

Directory
==================
  .. toctree::
     :maxdepth: 5
     :caption: General:
     
     getting_started
     about
     troubleshooting
     community

  .. toctree::
     :maxdepth: 5
     :caption: API Reference:
   
     gettingdata
     customdatasources
     clean_data_scheme
     plottingdata
     helpertools
     cli

  .. toctree::
     :maxdepth: 5
     :caption: Tutorials:  

     beginners_sounderpy_cookbook
     examplegallery
     tutorials/index
     case_studies/mound_city_20240828

* :ref:`genindex`
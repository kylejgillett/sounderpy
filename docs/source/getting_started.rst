.. _getting-started:

==================
🚀 Getting Started
==================

SounderPy is designed to make atmospheric vertical-profile retrieval, analysis, and visualization *straightforward*.

You can install SounderPy and create your first sounding in only a few minutes with only a few lines of code, in Python or on the command line!


Step 1: Install SounderPy
==========================

The easiest way to install SounderPy is with ``pip``:

.. code-block:: console

   pip install sounderpy

SounderPy is also available from ``conda-forge``:

.. code-block:: console

   conda install conda-forge::sounderpy

or:

.. code-block:: console

   mamba install conda-forge::sounderpy


You may also install a specific release version of SounderPy from GitHub:

.. code-block:: bash

   pip install "git+https://github.com/kylejgillett/sounderpy.git@v3.2.0"


or you may install the development version of SounderPy from Github:

.. code-block:: bash

   pip install "git+https://github.com/kylejgillett/sounderpy.git@main"

.. note::

   - Installing directly from GitHub requires ``git`` to be installed and available on your system.
   - The ``main`` branch may contain changes that have not yet been included in an official SounderPy release. For reproducible installations, use a tagged release such as ``v3.2.0``.

**********************************************************




Step 2: Import SounderPy
========================

A common convention is to import SounderPy as ``spy``:

.. code-block:: python

   import sounderpy as spy


**********************************************************

Step 3: Create Your First Sounding
====================================

Retrieve an observed sounding from Omaha, Nebraska (OAX) from
16 June 2014 at 18 UTC:

.. code-block:: python

   import sounderpy as spy

   data = spy.get_obs_data(
       "OAX",
       "2014", "06", "16", "18",
       hush=True,
   )

Then create a SounderPy sounding:

.. code-block:: python

   spy.build_sounding(data)

That's it!

SounderPy retrieves the profile, converts it into its standard
``clean_data`` structure, calculates the atmospheric parameters needed by the
plot, and creates the sounding.


Step 4: Examining the Data
============================

The object returned by SounderPy is a ``clean_data`` dictionary containing the
vertical profile and associated metadata.

Common fields include:

.. code-block:: text

   p          pressure
   z          height
   T          temperature
   Td         dewpoint
   u          u-component wind
   v          v-component wind
   site_info  profile metadata

Because supported data sources are converted into this same structure, the
same plotting and analysis functions can be used with observations, forecasts,
reanalyses, and supported model output.


Step 5: Where to Go Next
==========================

Choose the path that best matches what you want to do:

.. list-table::
   :widths: 30 70
   :class: home-card-table

   * - **👨‍🍳 Tutorials**
     - Learn SounderPy through complete, step-by-step workflows.

       :doc:`Open the tutorials <tutorials/index>`

   * - **🌎 Get Data**
     - Learn how to retrieve RAOB, ACARS, RAP/RUC, BUFKIT, ERA5, and other
       supported profiles.

       :doc:`Explore data sources <gettingdata>`

   * - **📊 Plot Gallery**
     - See examples of SounderPy soundings, hodographs, composites, parcels,
       radar, and plotting options.

       :doc:`Browse the gallery <examplegallery>`

   * - **💻 Command Line**
     - Retrieve, plot, and export profiles directly from a terminal.

       :doc:`Use the CLI <cli>`

   * - **📚 Reference**
     - Find detailed plotting, calculation, and helper-function documentation.

       :doc:`Open the plotting reference <plottingdata>`


A Complete First Example
========================

Putting everything together:

.. code-block:: python

   import sounderpy as spy

   data = spy.get_obs_data("OAX", "2014", "06", "16", "18", hush=True)

   spy.build_sounding(data, radar=None, map_zoom=0)


Ready for More?
===============

The :doc:`SounderPy Tutorials <tutorials/index>` are the best next step.
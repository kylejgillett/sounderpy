💡 Additional Resources
========================

SounderPy is designed to make working with vertical atmospheric profiles easier,
but users who are new to Python may find it helpful to first become familiar with
the basic scientific Python ecosystem. The resources below provide introductions
to Python, scientific computing, meteorological analysis, and other commonly used
atmospheric-science Python packages.

New to Python?
--------------

If you are completely new to Python, the following resources are good places to
start.

* `Python Downloads <https://www.python.org/downloads/>`_
  — Official Python installers and release information.

* `Project Pythia: Getting Started with Python <https://foundations.projectpythia.org/foundations/getting-started-python/>`_
  — A geoscience-focused introduction to installing, running, and managing Python.

* `Unidata: Introduction to Python for Atmospheric Science & Meteorology <https://unidata.github.io/python-training/python/intro-to-python/>`_
  — Introductory Python material written specifically for atmospheric-science users.

* `Official Python Tutorial <https://docs.python.org/3/tutorial/>`_
  — The official introduction to the Python language.

* `Python Virtual Environments <https://docs.python.org/3/library/venv.html>`_
  — An introduction to creating isolated Python environments with ``venv``.

Running and Writing Python
--------------------------

Python can be used from a terminal, an interactive notebook, or an integrated
development environment (IDE).

* `Project Pythia: Choosing a Python Platform <https://foundations.projectpythia.org/foundations/how-to-run-python/>`_
  — Overview of terminals, Jupyter notebooks, and IDEs.

* `Project Jupyter: Installing Jupyter <https://jupyter.org/install>`_
  — Instructions for installing JupyterLab and Jupyter Notebook.

* `Python in Visual Studio Code <https://code.visualstudio.com/docs/python/python-tutorial>`_
  — Guide to installing Python support, running scripts, using virtual
  environments, and debugging code in VS Code.

For users interested in scientific and atmospheric applications, Project Pythia
is an especially useful starting point because its examples are built around the
scientific Python ecosystem commonly used in the geosciences.

Scientific Python Fundamentals
------------------------------

Many atmospheric-science packages build upon a small group of widely used
scientific Python libraries. A good introduction to these packages is available
through the `Project Pythia Foundations <https://foundations.projectpythia.org/>`_
and its
`Core Scientific Python Packages <https://foundations.projectpythia.org/core/overview/>`_
section.

Common foundational packages include:

* `NumPy <https://numpy.org/>`_ — numerical arrays and numerical computation.
* `Matplotlib <https://matplotlib.org/>`_ — scientific plotting and visualization.
* `Pandas <https://pandas.pydata.org/>`_ — tabular and time-series data.
* `Xarray <https://docs.xarray.dev/en/latest/>`_ — labeled multidimensional and
  gridded data, including many NetCDF workflows.
* `SciPy <https://scipy.org/>`_ — scientific algorithms, interpolation,
  statistics, signal processing, and numerical methods.
* `Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ — geographic
  projections and meteorological mapping.

Meteorology-Focused Training
----------------------------

Several freely available training resources focus specifically on using Python
for meteorology and atmospheric science.

Project Pythia
~~~~~~~~~~~~~~

`Project Pythia <https://projectpythia.org/>`_ provides community-developed
training for Python-based computing in the geosciences.

Useful starting points include:

* `Pythia Foundations <https://foundations.projectpythia.org/>`_
* `Getting Started with Python <https://foundations.projectpythia.org/foundations/getting-started-python/>`_
* `Core Scientific Python Packages <https://foundations.projectpythia.org/core/overview/>`_

Pythia is particularly useful for users moving from basic Python into scientific
data analysis with NumPy, Matplotlib, Xarray, Pandas, and Cartopy.

Unidata Python Training
~~~~~~~~~~~~~~~~~~~~~~~

`Unidata Python Training <https://unidata.github.io/python-training/>`_ provides
Python lessons and examples specifically for atmospheric science and meteorology.

The
`Unidata Python Workshop <https://unidata.github.io/python-training/workshop/workshop-intro/>`_
includes examples involving upper-air data, Skew-T diagrams and hodographs,
weather-model output, satellite data, surface observations, Cartopy, MetPy,
Xarray, and Siphon.

MetPy
~~~~~

`MetPy <https://unidata.github.io/MetPy/latest/>`_ is a meteorological Python
library developed by Unidata. It provides tools for meteorological calculations,
units, thermodynamics, kinematics, plotting, upper-air analysis, and working with
Xarray datasets.

The `MetPy Tutorials <https://unidata.github.io/MetPy/latest/tutorials/index.html>`_
and example gallery are useful companions to SounderPy.

Other Atmospheric-Science Python Packages
-----------------------------------------

The following packages may also be useful depending on the type of atmospheric
data being analyzed. This is not intended to be an exhaustive list.

Data Access and Model Data
~~~~~~~~~~~~~~~~~~~~~~~~~~

* `Siphon <https://unidata.github.io/siphon/latest/>`_
  — Access to meteorological data services, including THREDDS-based sources.

* `Herbie <https://herbie.readthedocs.io/en/latest/>`_
  — Retrieval, subsetting, and reading of numerical weather prediction model data.

* `cfgrib <https://github.com/ecmwf/cfgrib>`_
  — GRIB and GRIB2 data access through the Xarray data model.

* `wrf-python <https://wrf-python.readthedocs.io/en/latest/>`_
  — Extraction, interpolation, diagnostics, and plotting utilities for WRF-ARW output.

Radar and Remote Sensing
~~~~~~~~~~~~~~~~~~~~~~~~

* `Py-ART <https://arm-doe.github.io/pyart/>`_
  — Reading, processing, analyzing, and visualizing weather-radar data.

* `ACT <https://arm-doe.github.io/ACT/>`_
  — Atmospheric data discovery, quality control, analysis, and visualization tools.

* `Satpy <https://satpy.readthedocs.io/en/latest/>`_
  — Reading, processing, resampling, and visualizing satellite remote-sensing data.

General Geoscience Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `GeoCAT-comp <https://geocat-comp.readthedocs.io/en/latest/>`_
  — Computational tools for geoscience data, including many workflows familiar
  to former NCL users.

* `Project Pythia <https://projectpythia.org/>`_
  — Community tutorials and cookbook-style workflows across the geosciences.

Severe Weather and Sounding Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `SHARPpy <https://sharppy.github.io/SHARPpy/>`_
  — Sounding and hodograph analysis tools commonly used in severe-weather meteorology.

Where Should I Start?
---------------------

For a new Python and SounderPy user, a reasonable learning path is:

#. Work through the
   `Project Pythia Getting Started with Python <https://foundations.projectpythia.org/foundations/getting-started-python/>`_
   material.
#. Become familiar with Jupyter or an editor such as VS Code.
#. Learn the basics of NumPy, Matplotlib, and Xarray.
#. Work through the
   `Unidata Python Training <https://unidata.github.io/python-training/>`_
   meteorology examples.
#. Explore the
   `MetPy Tutorials <https://unidata.github.io/MetPy/latest/tutorials/index.html>`_
   for meteorological calculations and plotting.
#. Return to the SounderPy documentation and tutorials to begin retrieving,
   analyzing, and plotting vertical atmospheric profiles.

You do not need to master every package listed above before using SounderPy.
These resources are intended as references that can be explored as your Python
and meteorological-analysis needs grow.

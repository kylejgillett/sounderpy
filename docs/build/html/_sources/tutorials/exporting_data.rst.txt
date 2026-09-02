.. _tutorial-exporting-data:

==================
💾 Exporting Data
==================

SounderPy ``clean_data`` profiles can be exported for use outside the plotting
API. The current file-writing helper supports:

* CSV;
* SHARPpy/NSHARP-style sounding files;
* CM1 ``input_sounding``-style files.

Exports are created with:

.. code-block:: python

   spy.to_file(
       file_type,
       clean_data,
       filename=...,
   )

***************************************************************
 

Retrieve an Example Profile
===========================

.. code-block:: python

   import sounderpy as spy

   data = spy.get_obs_data(
       "OAX",
       "2014", "06", "16", "18",
       hush=True,
   )

***************************************************************
 

CSV Export
==========

CSV is useful for spreadsheets, pandas, R, or general data exchange.

.. code-block:: python

   spy.to_file(
       "csv",
       data,
       filename="oax_20140616_18z.csv",
   )

The CSV contains the core sounding fields:

.. code-block:: text

   p
   z
   T
   Td
   u
   v

Units are removed when the values are written, so retain the SounderPy variable
definitions when using the file elsewhere:

.. list-table::
   :header-rows: 1
   :widths: 20 45 35

   * - Column
     - Variable
     - SounderPy unit
   * - ``p``
     - Pressure
     - hPa
   * - ``z``
     - Height
     - m
   * - ``T``
     - Temperature
     - degC
   * - ``Td``
     - Dewpoint
     - degC
   * - ``u``
     - U-component wind
     - knots
   * - ``v``
     - V-component wind
     - knots


***************************************************************
 
SHARPpy Export
==============

Export a SHARPpy-compatible sounding file with:

.. code-block:: python

   spy.to_file(
       "sharppy",
       data,
       filename="oax_20140616_18z.snd",
   )

The file contains pressure, height, temperature, dewpoint, wind direction, and
wind speed in the SHARPpy/NSHARP-style text format.

***************************************************************
 

Read a SHARPpy File Back into SounderPy
=======================================

SounderPy also provides ``from_sharppy()``:

.. code-block:: python

   data = spy.from_sharppy(
       "oax_20140616_18z.snd"
   )

The imported file is converted back into a SounderPy-style dictionary.

Depending on the source file, additional site metadata or plot titles may need
to be supplied before producing a fully annotated SounderPy figure.

***************************************************************
 

CM1 Export
==========

A SounderPy profile can be converted into CM1 ``input_sounding`` format:

.. code-block:: python

   spy.to_file(
       "cm1",
       data,
       filename="input_sounding",
   )

The exporter derives the thermodynamic quantities used by CM1 and writes the
profile in CM1 input format.

***************************************************************
 

AGL Height Handling for CM1
===========================

The ``convert_to_AGL`` argument controls whether exported CM1 vertical heights
are converted to above-ground-level values.

Default:

.. code-block:: python

   spy.to_file(
       "cm1",
       data,
       filename="input_sounding",
       convert_to_AGL=True,
   )

To retain the profile's existing height values:

.. code-block:: python

   spy.to_file(
       "cm1",
       data,
       filename="input_sounding",
       convert_to_AGL=False,
   )

***************************************************************
 

Default Filename
================

``filename`` is optional in the Python API. If omitted, SounderPy uses its
default output name:

.. code-block:: python

   spy.to_file(
       "csv",
       data,
   )

For reproducible workflows, explicitly providing a filename is recommended.


***************************************************************
 
CLI Export
==========

The CLI can infer output format from recognized filename extensions.

CSV:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --output oax.csv

SHARPpy:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --output oax.snd

CM1:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --output input_sounding.cm1

You can also specify the format explicitly:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --output profile.dat \
       --format csv


***************************************************************
 
Export and Plot in One CLI Command
==================================

The CLI retrieves the profile once and can use the result for several outputs:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --output oax.csv \
       --plot sounding \
       --map-zoom 0 \
       --plot-file oax.png

***************************************************************
 

Choosing an Export Format
=========================

Use **CSV** when:

* you want a simple tabular file;
* you plan to use pandas, spreadsheets, R, or another analysis tool.

Use **SHARPpy** when:

* you want a conventional sounding-text format;
* you want compatibility with SHARPpy-style workflows.

Use **CM1** when:

* you are preparing an idealized CM1 simulation;
* you want to turn a retrieved or observed SounderPy profile into model input.


***************************************************************
 
Next Steps
==========

Continue to :doc:`Command Line Workflows <cli_workflows>` for shell-based
retrieval, plotting, JSON, and batch examples.

See also:

* :doc:`Helper Tools </helpertools>`
* :doc:`Command Line Tool </cli>`

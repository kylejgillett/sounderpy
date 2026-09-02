.. _tutorial-cli-workflows:

=========================
💻 Command Line Workflows
=========================

SounderPy includes a command line interface (CLI) for common retrieval,
plotting, export, and machine-readable output workflows.

The CLI is useful when:

* you want a quick sounding without opening Python;
* you want to retrieve and save a profile in one command;
* you want to script repeated SounderPy operations in a shell;
* another program needs SounderPy data as JSON.

For the full command reference, see :doc:`Command Line Tool </cli>`.

***************************************************************
 
Check the Installation
======================

Display help:

.. code-block:: bash

   sounderpy --help

Check the installed version:

.. code-block:: bash

   sounderpy --version

The CLI can also be launched as a Python module:

.. code-block:: bash

   python -m sounderpy --help


***************************************************************
 
Observed Soundings
==================

Retrieve an observed RAOB/IGRA profile:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18

The default terminal output is a concise profile summary.

***************************************************************
 

Model/Reanalysis Profiles
=========================

Retrieve RAP/RUC:

.. code-block:: bash

   sounderpy model rap-ruc 44.58 -100.82 2024-08-28 18 \
       --box-size 0.25

***************************************************************
 

BUFKIT Forecasts
================

Most recent supported run:

.. code-block:: bash

   sounderpy bufkit gfs KMOP 6

Archived run:

.. code-block:: bash

   sounderpy bufkit gfs KMOP 6 \
       --run 2023-08-05 12


***************************************************************
 
ACARS Profiles
==============

List available profiles:

.. code-block:: bash

   sounderpy acars list 2024-05-21 18

Retrieve one profile returned by the list:

.. code-block:: bash

   sounderpy acars get 2024-05-21 18 PROFILE_ID


***************************************************************
 
Create a Sounding Plot
======================

Interactive:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot sounding

Save the plot:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot sounding \
       --map-zoom 0 \
       --plot-file oax_sounding.png

***************************************************************
 

Create a Hodograph
==================

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot hodograph \
       --map-zoom 0 \
       --plot-file oax_hodograph.png

Storm-relative:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot hodograph \
       --storm-relative \
       --map-zoom 0 \
       --plot-file oax_sr_hodograph.png


***************************************************************
 
Plot Styling
=============

Dark mode:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot sounding \
       --dark-mode \
       --map-zoom 0 \
       --plot-file oax_dark.png

Color-deficiency-friendly sounding:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot sounding \
       --color-blind \
       --map-zoom 0 \
       --plot-file oax_colorblind.png

Theta/theta-e:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot sounding \
       --show-theta \
       --map-zoom 0 \
       --plot-file oax_theta.png

Higher DPI:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot sounding \
       --map-zoom 0 \
       --dpi 200 \
       --plot-file oax_200dpi.png


***************************************************************
 
Export a Profile
================

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

***************************************************************
 

Combine Retrieval, Export, and Plotting
=======================================

One retrieval can produce several outputs:

.. code-block:: bash

   sounderpy model rap-ruc 44.58 -100.82 2024-08-28 18 \
       --box-size 0.25 \
       --output rap_profile.csv \
       --plot sounding \
       --map-zoom 0 \
       --plot-file rap_sounding.png


***************************************************************
 
JSON Output
===========

Write ``clean_data`` to standard output as JSON:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 --json

Redirect it to a file:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --json > oax.json

Validate/pretty-print the output with Python:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 --json | \
       python -m json.tool

SounderPy preserves quantity units in JSON using ``value`` and ``units``
fields.


***************************************************************
 
JSON and a Saved Plot
=====================

Machine-readable JSON can be combined with a saved plot:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --json \
       --plot sounding \
       --map-zoom 0 \
       --plot-file oax.png > oax.json

An interactive plot is intentionally not allowed with ``--json`` because
standard output must remain clean for machine use.

***************************************************************
 

Verbose Mode
============

The CLI hides most of the underlying SounderPy status output by default.

Show it with:

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 --verbose


***************************************************************
 
Simple Shell Batch Workflow
===========================

Because the CLI is a normal shell command, it can be used in loops.

For example:

.. code-block:: bash

   for station in OAX TOP DDC; do
       sounderpy obs "$station" 2024-05-21 00 \
           --output "${station}_20240521_00z.csv"
   done

Or save several hodographs:

.. code-block:: bash

   for station in OAX TOP DDC; do
       sounderpy obs "$station" 2024-05-21 00 \
           --plot hodograph \
           --map-zoom 0 \
           --plot-file "${station}_hodo.png"
   done

***************************************************************
 

Check Exit Status in Scripts
============================

The CLI returns a non-zero exit status when an error occurs, so shell scripts
can react to failures.

Example:

.. code-block:: bash

   if sounderpy obs OAX 2014-06-16 18 --output oax.csv; then
       echo "Sounding retrieved successfully"
   else
       echo "Sounding retrieval failed"
   fi

***************************************************************
 

Use ``python -m sounderpy``
===========================

Every CLI workflow can also use the module entry point:

.. code-block:: bash

   python -m sounderpy obs OAX 2014-06-16 18

This can be helpful when several Python environments are installed and you want
to explicitly use the SounderPy installation associated with a particular
Python interpreter.

***************************************************************
 

Next Steps
==========

You have now reached the end of the core tutorial series.

Useful reference pages:

* :doc:`Command Line Tool </cli>`
* :doc:`Getting Data </gettingdata>`
* :doc:`Plotting Data </plottingdata>`
* :doc:`Helper Tools </helpertools>`
* :doc:`Plot Gallery </examplegallery>`

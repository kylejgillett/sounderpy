.. _cli:

================================
💻 Command Line Interface (CLI)
================================

As of v3.2.0, SounderPy includes a command line interface (CLI) tool for retrieving, saving,
and plotting atmospheric vertical profile data directly from the command line.

After installing SounderPy, the CLI can be launched with either:

.. code-block:: bash

   sounderpy --help

or:

.. code-block:: bash

   python -m sounderpy --help


**********************************************************************

Checking the Installed Version
==============================

To print the installed SounderPy version:

.. code-block:: bash

   sounderpy --version


Gettings Data With The CLI
=============

The CLI data retrieval is organized around four primary commands:

.. code-block:: text

   sounderpy
   ├── obs
   ├── model
   ├── bufkit
   └── acars
       ├── list
       └── get

The commands provide access to:

- ``obs``: Observed RAOB and IGRAv2 vertical profiles.

- ``model``: Model and reanalysis vertical profiles, including RAP/RUC, ERA5, and NCEP datasets.

- ``bufkit``: Forecast vertical profiles from BUFKIT.

- ``acars``: Listing and retrieval of ACARS aircraft vertical profiles.

Help is available for the main CLI and for each command by using ``--help``



Observed Soundings
------------------

Observed RAOB or IGRAv2 profiles can be retrieved with:

.. code-block:: text

   sounderpy obs STATION YYYY-MM-DD HH

where ``STATION`` is a valid RAOB or IGRAv2 station identifier and ``HH`` is
the UTC hour from ``00`` through ``23``.

For example:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00

By default, SounderPy prints a concise summary of the returned profile,
including site information, valid time, number of vertical levels, and the
pressure and height ranges.


**********************************************************************

Model Reanalysis Data
--------------------------

Model reanalysis profiles are retrieved with:

.. code-block:: text

   sounderpy model MODEL LAT LON YYYY-MM-DD HH

For example, a RAP reanalysis profile can be requested with:

.. code-block:: bash

   sounderpy model rap-ruc 44.58 -100.82 2024-08-28 18




**Area Averaging**

The model retrieval command uses an area-average box size of 0.10 degrees by
default. A different box size can be supplied with ``--box-size``:

.. code-block:: bash

   sounderpy model rap-ruc 44.58 -100.82 2024-08-28 18 \
       --box-size 0.25

The box size must be greater than zero.




**Targeting a Dataset**

When supported by the underlying SounderPy retrieval function, a specific
dataset may be supplied with ``--dataset`` (such as with RAP & RUC data from NCEI):

.. code-block:: bash

   sounderpy model MODEL LAT LON YYYY-MM-DD HH \
       --dataset DATASET

If ``--dataset`` is not supplied, SounderPy uses its normal automatic dataset
selection logic.


**********************************************************************


BUFKIT Forecast Profiles
-------------------------

BUFKIT profiles use the following structure:

.. code-block:: text

   sounderpy bufkit MODEL STATION FHR

where ``FHR`` is the requested forecast hour.

For example:

.. code-block:: bash

   sounderpy bufkit gfs KMOP 6

If no model initialization time is given, SounderPy attempts to retrieve the
most recent available BUFKIT run.

**Archived BUFKIT Runs**

An archived initialization time can be selected with ``--run``:

.. code-block:: bash

   sounderpy bufkit gfs KMOP 6 \
       --run 2023-08-05 12

The two values following ``--run`` are the model initialization date
(``YYYY-MM-DD``) and UTC initialization hour (``HH``).


**********************************************************************

ACARS Aircraft Observations
-----------------------------

ACARS data use two CLI actions: ``list`` and ``get``.

**Listing Available Profiles**


To list available ACARS profile identifiers for a particular date and hour:

.. code-block:: bash

   sounderpy acars list 2024-05-21 18

Each available profile identifier is printed to the terminal.

The list may instead be returned as JSON:

.. code-block:: bash

   sounderpy acars list 2024-05-21 18 --json

**Retrieving a Profile**

After identifying an available ACARS profile, retrieve it with:

.. code-block:: text

   sounderpy acars get YYYY-MM-DD HH PROFILE_ID

For example:

.. code-block:: bash

   sounderpy acars get 2024-05-21 18 BNA_2320

Use a profile identifier returned by ``sounderpy acars list`` for the requested
date and hour.




**********************************************************************


Saving Data to a File
===========

Saving to JSON
---------------

Retrieved profile data can be written to standard output as JSON with
``--json``:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 --json

The JSON output preserves units for MetPy/Pint quantities. Quantity fields are
represented using separate ``value`` and ``units`` entries.

For example, a quantity may have the general form:

.. code-block:: json

   {
     "value": [
       1000.0,
       975.0,
       950.0
     ],
     "units": "hectopascal"
   }

JSON output can be redirected to a file:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 --json > dtx.json

It can also be piped into another command-line program:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 --json | python -m json.tool

When ``--json`` is used, SounderPy suppresses its normal status output so that
standard output contains valid JSON.

.. note::

   ``--json`` cannot be combined with an interactive plot. To produce JSON and
   a plot in the same command, save the plot using ``--plot-file``.

**********************************************************************


Saving to CSV, SHARPPY, CM1
----------------------------

Retrieved profiles can be saved to standard sounding file types from the CLI with ``-o`` or ``--output``:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --output dtx.csv

.. list-table:: Supported CLI export formats
   :header-rows: 1
   :widths: 25 35 40

   * - Format
     - ``--format`` value
     - Recognized filename extensions
   * - CSV
     - ``csv``
     - ``.csv``
   * - SHARPpy
     - ``sharppy``
     - ``.snd``, ``.sharp``, ``.sharppy``
   * - CM1
     - ``cm1``
     - ``.cm1``, ``.input``

For example:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --output dtx.snd

or:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --output input_sounding.cm1


**********************************************************************





Plotting from the CLI
=====================

SounderPy can generate a sounding or hodograph immediately after retrieving a
profile.

Sounding Plot
-------------

To open a sounding interactively:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot sounding

To save the sounding instead:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot sounding \
       --plot-file dtx_sounding.png

Hodograph
---------

To create a hodograph:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot hodograph

To save the hodograph:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot hodograph \
       --plot-file dtx_hodograph.png


**********************************************************************


Plot Options
-------------

Several SounderPy plotting options are available through the CLI.

**Dark Mode**

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot sounding \
       --dark-mode \
       --plot-file dtx_dark.png


**Color Blind Friendly Mode**

For sounding plots, ``--color-blind`` enables SounderPy's
color-deficiency-friendly temperature/dewpoint styling:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot sounding \
       --color-blind \
       --plot-file dtx_colorblind.png


**Storm-Relative Hodograph**

A hodograph can be transformed to storm-relative coordinates with
``--storm-relative``:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot hodograph \
       --storm-relative \
       --plot-file dtx_sr_hodograph.png

**Map Zoom**

The map inset zoom level can be changed with ``--map-zoom``. A value of zero
disables the map/radar inset:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot sounding \
       --map-zoom 0 \
       --plot-file dtx_sounding.png

Disabling the map inset can be useful when working offline or when a map/radar
inset is not needed.

**Plot Resolution**

The output resolution can be changed with ``--dpi``:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot sounding \
       --dpi 200 \
       --plot-file dtx_sounding.png

**********************************************************************




Combining Retrieval, Export, and Plotting
========================================

A primary advantage of the CLI is that one retrieved profile can be used for
multiple outputs without repeating the retrieval.

For example, the following command retrieves a RAP profile, exports the
profile to CSV, and saves a sounding image:

.. code-block:: bash

   sounderpy model rap-ruc 44.58 -100.82 2024-08-28 18 \
       --box-size 0.25 \
       --output rap_20240828_18z.csv \
       --plot sounding \
       --map-zoom 0 \
       --plot-file rap_20240828_18z.png



**********************************************************************



Verbose Output
==============

The CLI suppresses most of SounderPy's internal retrieval and plotting status
messages by default, producing cleaner command-line output.

To display the underlying SounderPy status messages, use ``-v`` or
``--verbose``:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 --verbose


**********************************************************************


Input Validation and Errors
===========================

The CLI validates several common inputs before beginning a retrieval.

Dates
-----

Dates must use the form:

.. code-block:: text

   YYYY-MM-DD

Hours
-----

Hours are interpreted as UTC and must be integers from 0 through 23. Both of the following are accepted:

.. code-block:: text

   0
   00

Latitude and Longitude
----------------------

For model retrieval:

* latitude must be between -90 and 90 degrees;
* longitude must be between -180 and 180 degrees.

Forecast Hours
--------------

BUFKIT forecast hours must be zero or greater.

Box Size
--------

The model ``--box-size`` value must be greater than zero.

Invalid command-line arguments produce a descriptive error and a non-zero
exit status.


Additional Examples
===================

Retrieve an observed profile and save it as CSV:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --output dtx.csv

Retrieve RAP reanalysis and save a sounding:

.. code-block:: bash

   sounderpy model rap-ruc 44.58 -100.82 2024-08-28 18 \
       --box-size 0.25 \
       --plot sounding \
       --map-zoom 0 \
       --plot-file rap_sounding.png

Retrieve an archived BUFKIT profile and export it:

.. code-block:: bash

   sounderpy bufkit gfs KMOP 6 \
       --run 2023-08-05 12 \
       --output kmop_gfs_f006.csv

Create a storm-relative hodograph from an observed sounding:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot hodograph \
       --storm-relative \
       --map-zoom 0 \
       --plot-file dtx_sr_hodograph.png

Retrieve JSON for use in another program:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 --json > dtx.json


Using the Python Module Entry Point
==================================

The CLI may also be executed through Python:

.. code-block:: bash

   python -m sounderpy obs DTX 2024-05-21 00

This form supports the same commands and options as the ``sounderpy`` console
command.

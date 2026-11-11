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


Getting Data With The CLI
=========================

The CLI data retrieval is organized around five primary commands:

.. code-block:: text

   sounderpy
   ├── obs
   ├── model
   ├── bufkit
   ├── acars
   │   ├── list
   │   └── get
   └── file

The commands provide access to:

- ``obs``: Observed RAOB and IGRAv2 vertical profiles.

- ``model``: Model and reanalysis vertical profiles, including RAP/RUC, ERA5, and NCEP datasets.

- ``bufkit``: Forecast vertical profiles from BUFKIT.

- ``acars``: Listing and retrieval of ACARS aircraft vertical profiles.

- ``file``: A local SHARPpy-format sounding file.

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

``MODEL`` is ``rap-ruc`` (or ``rap`` / ``ruc``), ``era5``, or ``ncep``. ERA5 comes from the Copernicus Climate
Data Store and needs a free CDS API key saved in ``~/.cdsapirc``.


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

where ``FHR`` is the requested forecast hour and ``MODEL`` is one of ``gfs``, ``rap``, ``hrrr``, or ``rrfs``.
NAM and NAM Nest (``nam``, ``namnest``) are discontinued, so they are only available for archived runs
(``--run``).

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

   sounderpy acars get 2024-05-21 18 BNA_1820

Use a profile identifier returned by ``sounderpy acars list`` for the requested
date and hour.


**********************************************************************

Local SHARPpy Files
-------------------

A SHARPpy-format sounding file, such as one saved earlier with ``--output profile.snd``, can be read with:

.. code-block:: text

   sounderpy file PATH

For example, to export it as CSV:

.. code-block:: bash

   sounderpy file oax.snd --output oax.csv

A SHARPpy file holds the station ID and valid time but not the location, so ``--plot`` and ``--params``
need ``--latlon LAT LON`` (the map inset, and the latitude used for temperature advection):

.. code-block:: bash

   sounderpy file oax.snd \
       --latlon 41.32 -96.37 \
       --plot sounding \
       --plot-file oax_sounding.png

The plot title uses the file name, station ID, and valid time, and the station elevation is taken from
the lowest height in the file.




**********************************************************************


Saving Data to a File
=====================

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
standard output contains valid JSON. With ``--json --verbose``, the status
output is written to standard error instead, so standard output stays valid JSON.

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

Opening a plot in a window needs a display and an interactive Matplotlib backend. Without one
(for example, on a server or over SSH), the CLI stops before downloading anything and asks for
``--plot-file`` instead.

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

The ``--plot-file`` extension sets the image format (``.png``, ``.pdf``, ``.svg``, ``.jpg``, and the
other formats Matplotlib writes). A file name without an extension is saved as ``.png``.

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

Several SounderPy plotting options are available through the CLI. They all need ``--plot``.
``--color-blind``, ``--show-theta``, ``--mode``, ``--obs``, ``--parcels``, and ``--background-parcels``
apply to ``--plot sounding`` only. ``--storm-motion``, ``--modify-sfc``, and ``--mode`` change the
calculation as well, so they also work with ``--params`` (see `Calculated Parameters`_).

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


**Theta and Theta-e**

For sounding plots, ``--show-theta`` shows theta and theta-e:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot sounding \
       --show-theta \
       --plot-file dtx_theta.png


**Storm-Relative Hodograph**

A hodograph can be transformed to storm-relative coordinates with
``--storm-relative``. It works with ``--plot hodograph`` and with the hodograph on
``--plot sounding``:

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

**Winter Sounding**

``--mode`` chooses the sounding plot: ``severe`` (the default), ``winter``, or
``auto``, which draws the winter sounding when the surface wet-bulb temperature
is below 35 degrees F. ``--obs`` adds the latest surface observation near the
profile to a winter sounding (from the Iowa Environmental Mesonet), so it needs
``--mode winter`` or ``--mode auto``:

.. code-block:: bash

   sounderpy bufkit hrrr KGFK 1 \
       --run 2024-01-24 12 \
       --plot sounding \
       --mode winter \
       --obs \
       --plot-file kgfk_winter.png

See :doc:`/tutorials/winter_soundings` for a guide to the winter sounding. The winter sounding
doesn't use ``--show-theta``, ``--storm-relative``, ``--parcels``, or ``--background-parcels``, so
they can't be combined with ``--mode winter``.

**Storm Motion**

``--storm-motion`` sets the storm motion for the plot and the parameters: ``right_moving`` (the
default), ``left_moving``, ``mean_wind``, or a direction in degrees and a speed in knots:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot hodograph \
       --storm-motion 250 30 \
       --plot-file dtx_hodograph.png

**Surface Modification**

``--modify-sfc`` changes the surface with ``KEY=VALUE`` pairs: ``T`` and ``Td`` (°C), ``ws`` (kt)
with ``wd`` (degrees), ``method`` (``mix``, ``replace``, or ``blend``), and ``depth`` (m, for
``blend``). See :ref:`modify_sfc` for what each method does. The change applies to the plot and
``--params``; a file saved with ``--output`` keeps the original surface.

.. code-block:: bash

   sounderpy bufkit hrrr KGFK 1 \
       --run 2024-01-24 12 \
       --plot sounding \
       --mode winter \
       --modify-sfc T=-6 Td=-7 method=blend depth=300 \
       --plot-file kgfk_blend.png

**Radar**

``--radar`` chooses the map's radar: ``mosaic`` (the default; the CONUS mosaic is only available for
recent dates), ``single`` (the nearest NEXRAD site), or ``none``. ``--radar-time`` draws the radar at
the profile's valid time (``sounding``, the default) or ``now``. Both need the map, so they can't be
combined with ``--map-zoom 0``.

.. code-block:: bash

   sounderpy obs OAX 2014-06-16 18 \
       --plot sounding \
       --radar single \
       --plot-file oax_radar.png

**Special Parcels**

``--parcels`` draws special parcels from the ``ecape_parcel`` library on a severe sounding, by code
(``mu``, ``ml``, or ``sb``; ``ps`` pseudoadiabatic or ``ia`` irreversible adiabatic; ``cape`` or
``ecape``), and ``--background-parcels`` draws more of them faded behind. ``--parcels simple`` draws
only the standard parcels and skips the default MUECAPE trace, which is faster. See
:ref:`parcels_logic` for the parcel codes.

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 \
       --plot sounding \
       --parcels mu_ia_ecape \
       --background-parcels mu_ps_cape sb_ps_cape \
       --plot-file dtx_parcels.png

**********************************************************************




Combining Retrieval, Export, and Plotting
=========================================

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



Calculated Parameters
=====================

``--params`` prints key parameters under the profile summary: CAPE, CIN, LCL, and ECAPE for each
parcel, DCAPE, PWAT, lapse rates, freezing levels, the storm motion, SRH, bulk shear, and composite
indices. With ``--mode winter`` (or ``--mode auto`` when the surface wet-bulb is below 35 °F), it
adds the precipitation type, the phase from source to surface, the snow ratio, the DGZ, and the
freezing drizzle and freezing fog levels.

.. code-block:: bash

   sounderpy bufkit hrrr KGFK 1 \
       --run 2024-01-24 12 \
       --params \
       --mode auto

``--storm-motion`` and ``--modify-sfc`` change the parameters the same way they change the plot.
With ``--json``, every calculated value is added to the JSON as ``"params"``, next to the profile:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 --params --json | python -c "import json, sys; print(json.load(sys.stdin)['params']['thermo']['mlcape'])"


**********************************************************************


Verbose Output
==============

The CLI hides SounderPy's internal retrieval and plotting status messages and
warnings by default, producing cleaner command-line output.
SounderPy's own warnings about the figure or data, such as a radar mosaic that
isn't available for the profile's date, are still shown, one line each:

.. code-block:: text

   sounderpy: warning: Radar mosaic could not be plotted: no radar mosaic for 08/05/2023 18:00Z: the mosaic archive only keeps about the last month. Use the nearest NEXRAD site for older dates (radar='single', or --radar single in the CLI)

To display the underlying SounderPy status messages and warnings, and the full
Python traceback when something fails, use ``-v`` or ``--verbose``:

.. code-block:: bash

   sounderpy obs DTX 2024-05-21 00 --verbose


**********************************************************************


Input Validation and Errors
===========================

The CLI validates its inputs before beginning a retrieval, so a mistake is reported
without waiting for a download. Besides the values below, it checks that:

* ``--output`` has a recognized extension, or ``--format`` is given (and ``--format`` is only used with ``--output``);
* plot options are only used with ``--plot``, and sounding-only options only with ``--plot sounding``;
* ``--storm-motion``, ``--modify-sfc``, and ``--mode`` are used with ``--plot`` or ``--params``, and their values are valid;
* ``--parcels`` codes are valid, and severe-only options aren't combined with ``--mode winter``;
* ``--radar`` and ``--radar-time`` aren't combined with ``--map-zoom 0``;
* ``file`` has ``--latlon`` for ``--plot`` and ``--params``;
* ``--plot-file`` has an extension Matplotlib can write;
* ``--json`` is combined with a saved plot, not a plot window.

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

Map Zoom and DPI
----------------

``--map-zoom`` must be zero or greater, and ``--dpi`` greater than zero.

Box Size
--------

The model ``--box-size`` value must be greater than zero.

Invalid command-line arguments produce a descriptive error and a non-zero
exit status:

* ``0``: success;
* ``1``: the retrieval, export, or plot failed (for example, no data for that station and time);
* ``2``: invalid arguments or options;
* ``130``: interrupted with Ctrl+C.


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
===================================

The CLI may also be executed through Python:

.. code-block:: bash

   python -m sounderpy obs DTX 2024-05-21 00

This form supports the same commands and options as the ``sounderpy`` console
command.

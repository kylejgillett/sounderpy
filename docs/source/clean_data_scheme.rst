.. _clean_data_scheme:

🧼 The ``clean_data`` Schema
==============================

SounderPy converts supported atmospheric-profile data sources into a common dictionary structure called ``clean_data``. This common structure is the interface between SounderPy's data-retrieval, analysis, plotting, and export tools.

In other words:

.. code-block:: text

   RAOB ────────┐
   ACARS ───────┤
   BUFKIT ──────┤
   RAP / RUC ───┤
   ERA5 ────────┼──> clean_data ──> calculations
   NCEP ────────┤                  ├─> sounding
   WRF ─────────┤                  ├─> hodograph
   CM1 ─────────┤                  ├─> composite
   custom data ─┘                  └─> export

Once a profile has been converted into ``clean_data``, the same SounderPy functions can be used regardless of the original data source.

*************************************


Core Schema
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A valid SounderPy profile contains the following core fields:

.. list-table::
   :header-rows: 1
   :widths: 15 25 20 40

   * - Key
     - Variable
     - Standard unit
     - Description

   * - ``p``
     - Pressure
     - hPa
     - Atmospheric pressure at each profile level.

   * - ``z``
     - Height
     - m
     - Vertical coordinate associated with each profile level.

   * - ``T``
     - Temperature
     - degC
     - Environmental air temperature.

   * - ``Td``
     - Dewpoint
     - degC
     - Environmental dewpoint temperature.

   * - ``u``
     - U-component wind
     - knots
     - Zonal wind component.

   * - ``v``
     - V-component wind
     - knots
     - Meridional wind component.

   * - ``site_info``
     - Profile metadata
     - —
     - Dictionary describing the source, location, and valid time of the
       profile.

*************************************



Profile Array Requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The six atmospheric-profile variables:

.. code-block:: text

   p
   z
   T
   Td
   u
   v

must follow several rules.

**1. Each variable must be one-dimensional.**

For example:

.. code-block:: python

   data["T"].shape

might return:

.. code-block:: text

   (73,)


**2. All profile arrays must have the same length.**

A single index must represent the same atmospheric level across all variables:

.. code-block:: python

   i = 10

   print(data["z"][i])
   print(data["p"][i])
   print(data["T"][i])
   print(data["Td"][i])
   print(data["u"][i])
   print(data["v"][i])


**3. Profile arrays must carry physical units.**

SounderPy uses Pint/MetPy quantities rather than unitless NumPy arrays.

For example:

.. code-block:: python

   print(data["p"][0])
   print(data["T"][0])
   print(data["u"][0])

may produce values similar to:

.. code-block:: text

   970 hectopascal
   24.3 degree_Celsius
   12.5 knot


**4. The vertical profile should be ordered from the surface upward.**

Generally:

* ``z`` increases monotonically with array index;
* ``p`` decreases monotonically with array index.

Conceptually:

.. code-block:: text

   index       z          p
   -----     ------     ------
     0         0 m       970 hPa
     1       250 m       942 hPa
     2       500 m       915 hPa
     ...       ...         ...
     n      15000 m       120 hPa


*************************************




Units
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SounderPy normally uses the following units internally:

.. list-table::
   :header-rows: 1
   :widths: 20 40

   * - Variable
     - Unit

   * - ``p``
     - hectopascals (``hPa``)

   * - ``z``
     - meters (``m``)

   * - ``T``
     - degrees Celsius (``degC``)

   * - ``Td``
     - degrees Celsius (``degC``)

   * - ``u``
     - knots (``kt``)

   * - ``v``
     - knots (``kt``)

These values should remain unit-aware when constructing custom SounderPy
profiles.

For example:

.. code-block:: python

   from metpy.units import units

   pressure = pressure_values * units.hPa
   height = height_values * units.m
   temperature = temperature_values * units.degC
   dewpoint = dewpoint_values * units.degC
   u_wind = u_values * units.kt
   v_wind = v_values * units.kt

*************************************





``site_info`` Metadata
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Profile metadata are stored in:

.. code-block:: python

   data["site_info"]

The exact metadata available can vary with the original data source.

Common SounderPy metadata fields include:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Key
     - Meaning

   * - ``site-id``
     - Station, airport, or profile identifier.

   * - ``site-name``
     - Human-readable site name.

   * - ``site-lctn``
     - Location description.

   * - ``site-latlon``
     - ``[latitude, longitude]`` pair.

   * - ``site-elv``
     - Surface/site elevation.

   * - ``source``
     - Description of the original data source.

   * - ``model``
     - Model identifier when applicable.

   * - ``fcst-hour``
     - Forecast hour when applicable.

   * - ``run-time``
     - Model initialization time when applicable.

   * - ``valid-time``
     - Valid time of the atmospheric profile.

Some data sources may contain additional metadata such as ``box_area``.

For example, a model profile may contain:

.. code-block:: python

   data["site_info"] = {
       "site-id": "no-site-id",
       "site-name": "no-site-name",
       "site-lctn": "no-site-location",
       "site-latlon": [44.58, -100.82],
       "site-elv": 520.0,
       "source": "MODEL REANALYSIS",
       "model": "RAP",
       "fcst-hour": "F00",
       "run-time": ["2024", "08", "28", "18"],
       "valid-time": ["2024", "08", "28", "18"],
       "box_area": "0.25° x 0.25° BOX AVG",
   }

*************************************



Plot Titles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Many SounderPy-generated profiles also contain a ``titles`` dictionary used for figure annotation.

.. code-block:: python

   data["titles"]

A typical example is:

.. code-block:: python

   {
       "top_title": "MODEL REANALYSIS VERTICAL PROFILE | 18Z RAP F00",
       "left_title": "18Z RAP F00 | VALID: 08/28/2024 18Z",
       "right_title": "44.58, -100.82 | 0.25° x 0.25° BOX AVG",
   }



*************************************



Optional Profile Variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some SounderPy data sources contain additional atmospheric variables.

Examples include:

.. list-table::
   :header-rows: 1
   :widths: 20 30 25

   * - Key
     - Variable
     - Typical unit

   * - ``rh``
     - Relative humidity
     - percent

   * - ``omega``
     - Pressure vertical velocity
     - Pa s^-1

These fields are **not part of the minimum ``clean_data`` contract** required for general SounderPy plotting and analysis.

Additional source-specific fields may also be present. For example ``omega`` is available for model forecasts and reanalysis but not for observations.

*************************************



Example ``clean_data`` Dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simplified SounderPy profile can be represented as:

.. code-block:: python

   from metpy.units import units

   clean_data = {
       "p": [1000, 950, 900, 850] * units.hPa,
       "z": [0, 500, 1000, 1500] * units.m,
       "T": [24, 21, 18, 14] * units.degC,
       "Td": [20, 18, 14, 10] * units.degC,
       "u": [10, 15, 20, 25] * units.kt,
       "v": [5, 10, 15, 20] * units.kt,

       "site_info": {
           "site-id": "EXAMPLE",
           "site-name": "Example Profile",
           "site-lctn": "North Dakota",
           "site-latlon": [47.9, -97.0],
           "site-elv": 250,
           "source": "CUSTOM PROFILE",
           "model": "no-model",
           "fcst-hour": "no-fcst-hour",
           "run-time": ["none", "none", "none", "none"],
           "valid-time": ["2026", "09", "01", "18"],
       },

       "titles": {
           "top_title": "CUSTOM VERTICAL PROFILE",
           "left_title": "VALID: 09/01/2026 18Z",
           "right_title": "47.9, -97.0",
       },
   }


*************************************



Inspecting a Profile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can inspect the structure of any retrieved SounderPy profile with:

.. code-block:: python

   print(data.keys())

and:

.. code-block:: python

   print(data["site_info"])

To inspect units:

.. code-block:: python

   for key in ["p", "z", "T", "Td", "u", "v"]:
       print(key, data[key].units)


*************************************


Why the Schema Matters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``clean_data`` schema allows SounderPy's major tools to operate independently of the original data source.

For example:

.. code-block:: python

   spy.build_sounding(data)

   spy.build_hodograph(data)

   general, thermo, kinem, intrp = spy.sounding_params(data).calc()

   spy.to_file("csv", data)

all consume the same general profile structure.

This means data retrieval and atmospheric analysis remain separate:

.. code-block:: text

   retrieve / import
          ↓
      clean_data
          ↓
   ┌──────┼─────────┐
   ↓      ↓         ↓
   plot  calculate  export


*************************************




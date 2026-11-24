.. _readout_file:

===================
📄 The Readout File
===================

``readout=True`` saves ``NAME.readout.json`` next to ``NAME.png``. ``sounderpy-readout.js`` reads it, but it's plain
JSON, so any tool can map the figure's pixels to data with it. **Experimental:** the format may change between
releases; ``version`` says which one a file uses.

***************************************************************

Overview
========

About 30 KB for a full sounding (a few KB compressed). Pixels are the saved image's, with the origin at the top left.

.. code-block:: javascript

   {
     "format": "sounderpy-readout", "version": 1, "sounderpy": "3.3.0",
     "figure": "sounding",                      // "sounding", "winter", or "hodograph"
     "image": {"width": 2454, "height": 1533},
     "dark": false,                             // the figure's dark mode
     "site": {"id": "KOKC", "name": "...", "source": "...", "valid": ["2026", "11", "23", "21"], "latlon": [35.4, -97.6]},
     "edited": null,                            // what was changed, for a redrawn edited profile
     "units": {"p": "hPa", "z": "m", "T": "degC", "u": "kt", ...},
     "panels": {"skewt": {...}, "hodograph": {...}},
     "profile": {"p": [...], "T": [...], ...},
     "hodograph_trace": {"z_agl": [...], "u": [...], "v": [...]},
     "storm_motion": [u, v],
     "params": [{...}, ...],
     "height_panels": [{...}, ...]
   }

A hodograph figure has no ``skewt`` panel. ``hodograph_trace``, ``storm_motion``, ``params``, and ``height_panels``
are left out when the figure doesn't have them.

***************************************************************

Panels: pixels to data
======================

.. code-block:: javascript

   "panels": {
     "skewt":     {"bbox": [x0, y0, x1, y1], "transform": [[a, b, c], [d, e, f]],
                   "coords": ["T_degC", "log10_p_hPa"], "p_limits": [1050, 100]},
     "hodograph": {"bbox": [x0, y0, x1, y1], "transform": [[a, b, c], [d, e, f]],
                   "coords": ["u_kt", "v_kt"], "offset": [0, 0], "max_height": 12000}
   }

``bbox`` is the panel's area in pixels. ``transform`` maps the panel's data to pixels exactly (a skew-T is straight
lines in temperature and log pressure):

- **Skew-T:** ``x = a·T + b·log10(p) + c`` and ``y = d·T + e·log10(p) + f``, with T in °C and p in hPa.
- **Hodograph:** the same in ``(u - offset[0], v - offset[1])``, in kt. ``offset`` is the storm motion on a
  storm-relative hodograph, otherwise ``[0, 0]``.

To read a pixel, solve for the data:

.. code-block:: python

   import json

   ro = json.load(open('kokc.readout.json'))
   (a, b, c), (d, e, f) = ro['panels']['skewt']['transform']

   def skewt_at(x, y):
       det = a * e - b * d
       T = (e * (x - c) - b * (y - f)) / det
       logp = (a * (y - f) - d * (x - c)) / det
       return T, 10 ** logp          # °C, hPa

On a page that shows the image at another size, scale the pointer's position by ``image.width`` over the shown
width first.

***************************************************************

Profile
=======

Each list has one value per plotted level, surface first; ``null`` where a value is missing.

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 30 20 50

   * - Key
     - Units
     - Value
   * - ``p``, ``z``, ``z_agl``
     - hPa, m, m
     - Pressure, height above sea level, height above ground.
   * - ``T``, ``Td``, ``Tw``, ``Tv``
     - °C
     - Temperature, dewpoint, wet-bulb, virtual temperature.
   * - ``rh``, ``mr``
     - %, g/kg
     - Relative humidity, mixing ratio.
   * - ``theta``, ``thetae``
     - K
     - Potential and equivalent potential temperature.
   * - ``u``, ``v``, ``wdir``, ``wspd``
     - kt, °
     - Wind components, the direction it comes from, and speed.
   * - ``omega``
     - Pa/s
     - Vertical velocity, when the data has it.

``hodograph_trace`` is the wind every 100 m above ground, up to ``max_height``, as the figure's hodograph draws it.
``storm_motion`` is the ground-relative storm motion ``[u, v]`` (kt) that the figure's SRH and storm-relative values
use.

***************************************************************

Params
======

Each value in the figure's tables:

.. code-block:: javascript

   {"key": "kinem.srh_0_to_3000", "name": "0-3 km storm-relative helicity", "value": 100.2, "units": "m²/s²",
    "about": "Streamwise vorticity in the storm's inflow over the layer, ...",
    "bbox": [x0, y0, x1, y1], "layer": {"z": [0, 3000]}}

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 20 80

   * - Key
     - Value
   * - ``key``
     - Where it comes from in SounderPy's parameters, e.g. ``thermo.mlcape``.
   * - ``name``, ``units``, ``about``
     - A readable name, its units, and a sentence about it.
   * - ``value``
     - A number; a wind ``[u, v]`` in kt with ``"kind": "wind"``; or text, such as a precipitation type, with
       ``"kind": "text"``.
   * - ``bbox``
     - The value's text on the figure, in pixels.
   * - ``layer``
     - What it describes: ``{"z": [bottom, top]}`` (m above ground), ``{"p": [bottom, top]}`` (hPa),
       ``{"z_at": height}`` or a list of heights (a single height such as an LCL or a freezing level), or ``null``.
   * - ``extra``
     - Optional: more ``[label, value]`` rows, such as the checks behind a winter value.

***************************************************************

Height panels
=============

Each 0-3 km panel beside the skew-T:

.. code-block:: javascript

   {"title": "Storm-relative wind", "bbox": [x0, y0, x1, y1], "transform": [[a, b, c], [d, e, f]],
    "coords": ["x", "z_agl_m"], "z_limits": [0, 3000],
    "z": [0, 100, ...], "series": [{"label": "SR wind", "units": "kt", "values": [...]}]}

``transform`` maps ``(x, height m AGL)`` to pixels; ``series`` has the plotted values at each ``z``.

***************************************************************

Versions
========

``version`` goes up when the format changes in a way an older script can't read. ``sounderpy-readout.js`` warns
in the console when a file's version is newer than it knows; copy the script from the same SounderPy that wrote the
file (``spy.readout_js_path()``).

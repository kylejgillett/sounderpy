.. _webreadouts:

=========================
🖱️ Web Readouts
=========================

Put a SounderPy figure on a web page and let readers hover it for the data, or drag the profile to edit it.
SounderPy saves a small readout file next to the figure, and ``sounderpy-readout.js`` (shipped with SounderPy, no
dependencies) lays a layer over the ``<img>`` that reads it.

- **Hover** the skew-t or hodograph: pressure, height, temperature, dewpoint, wet-bulb, RH, mixing ratio, θe, wind,
  and omega at that level, with the level marked on the skew-t, hodograph, and 0-3 km panels.
- **Hover** a 0-3 km panel: its values at that height, with the same level readout.
- **Hover** a value in the parameter tables: its full name, exact value, and what it measures, with its layer
  shaded on the skew-t, hodograph, and 0-3 km panels (e.g. 0-3 km SRH, the 700-500 hPa lapse rate, or a winter
  sounding's melting layer, DGZ, and freezing levels).
- **Edit** (optional): drag temperature and dewpoint points along their pressure level, or hodograph points anywhere,
  then send the edited profile to your server to redraw it with SounderPy.

Works with ``build_sounding()`` (severe and winter) and ``build_hodograph()``, at any size the image is shown,
with mouse or touch.

***************************************************************

1. Save the figure with a readout
=================================

.. code-block:: python

   import sounderpy as spy

   clean_data = spy.get_bufkit_data('hrrr', 'kokc', 3)
   spy.build_sounding(clean_data, save=True, filename='kokc.png', readout=True)
   # writes kokc.png and kokc.readout.json

From the command line, add ``--readout`` to a ``--plot-file`` (see :doc:`/reference/cli`).

2. Add the script to the page
=============================

Copy the script next to your page:

.. code-block:: python

   import shutil
   shutil.copy(spy.readout_js_path(), 'site/')

and attach it to the image:

.. code-block:: html

   <img id="sounding" src="kokc.png" alt="KOKC sounding">

   <script type="module">
     import { SounderPyReadout } from "./sounderpy-readout.js";

     const readout = await SounderPyReadout.attach(document.getElementById("sounding"), "kokc.readout.json");
   </script>

The image can be scaled or resized freely, including inside a zoom-and-pan viewer that zooms it with a CSS transform:
the readout follows the zoom, and dragging a point doesn't pan the view. Keep the PNG and its readout file together: a
readout only fits the image it was saved with.

Options
-------

.. code-block:: javascript

   SounderPyReadout.attach(img, "kokc.readout.json", {
     fields: ["p", "z_agl", "T", "Td", "Tw", "rh", "mr", "thetae", "wind", "omega"],  // rows, in order
     temperature: "both",   // "C", "F", or "both"
     height: "m",           // "m" or "ft"
     theme: "auto",         // "auto" follows the figure's dark mode; or "light", "dark"
     editable: false,       // start in edit mode
     tapReadout: true,      // touch screens: a tap shows the readout (false leaves taps to the page unless editing)
     parts: ["skewt", "hodograph", "panels", "params"],  // what responds to the pointer
     onHover: (level) => {},       // the level under the pointer, or null
     onEdit: (profile) => {},      // while a point is dragged
     onEditEnd: (profile) => {},   // when a drag ends
   });

Other ``fields``: ``z`` (height MSL), ``Tv``, and ``theta``. The tooltip's colors are CSS variables on
``.spy-readout-tip`` (``--spy-readout-bg``, ``--spy-readout-text``, ``--spy-readout-muted``, ``--spy-readout-border``).

Methods
-------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Method
     - Does
   * - ``update(readout)``
     - Swaps in a new readout file (after a redraw); clears edits.
   * - ``setEditing(on)``
     - Turns editing on or off.
   * - ``setEnabled(on)``
     - Shows or hides the layer (e.g. a readout mode that can be switched off).
   * - ``getProfile()``
     - The edited profile: ``p`` (hPa), ``z`` (m), ``T`` and ``Td`` (°C), ``u`` and ``v`` (kt), one value per level.
   * - ``isEdited()``, ``undo()``, ``resetEdits()``
     - Whether anything changed; undo the last drag; go back to the figure's profile.
   * - ``levelAt(p)``
     - The readout values at a pressure (hPa).
   * - ``destroy()``
     - Removes the layer.

***************************************************************

3. Redraw an edited profile
===========================

While editing, the readout's values follow the edited lines. To recompute everything (CAPE, SRH, and the rest), send
``getProfile()`` to a server and redraw with :py:func:`spy.apply_profile_edits`:

.. code-block:: javascript

   const res = await fetch("/redraw", { method: "POST", body: JSON.stringify(readout.getProfile()) });

.. code-block:: python

   edited = spy.apply_profile_edits(clean_data, profile)   # clean_data: the profile the figure was drawn from
   spy.build_sounding(edited, save=True, filename='kokc_edited.png', readout=True)

Send the new image and readout file back, then call ``readout.update(newReadout)``. The redrawn figure says
**EDITED PROFILE** and lists what changed, so it can't pass as the original data.

.. py:function:: spy.apply_profile_edits(clean_data, edits)

   Rebuild ``clean_data`` from a profile edited in the browser.

   :param clean_data: the profile the readout was made from.
   :type clean_data: dict, required
   :param edits: the edited profile from ``getProfile()``: lists of ``p`` (hPa), ``z`` (m), ``T`` and ``Td`` (°C), and ``u`` and ``v`` (kt), one value per plotted level. A dewpoint above the temperature is lowered to it.
   :type edits: dict, required
   :return: the edited profile, with the original's site information and titles, other per-level values (such as omega) moved to the edited levels, and ``site_info['edited']`` describing the change.
   :rtype: dict

.. py:function:: spy.readout_js_path()

   :return: the path of the installed ``sounderpy-readout.js``.
   :rtype: pathlib.Path

***************************************************************

The readout file
================

JSON, about 30 KB for a full sounding. ``panels`` maps data to image pixels (origin top left) with a 2×3 affine
matrix: the skew-t in (temperature °C, log10 pressure hPa), the hodograph in (u, v) kt, less ``offset`` (the storm
motion of a storm-relative hodograph). ``profile`` has the plotted levels; ``hodograph_trace`` the wind every 100 m up
to 12 km. ``params`` has each value in the tables, with its ``layer``: ``{"z": [bottom, top]}`` (m AGL), ``{"p": [bottom, top]}``
(hPa), ``{"z_at": height}`` (or a list of heights) for single heights such as an LCL or freezing levels, or ``null``.
A value is a number, a wind ``[u, v]`` in kt (``"kind": "wind"``), or text such as a precipitation type
(``"kind": "text"``); ``extra`` lists more rows, such as a freezing fog check and whether it's met. ``height_panels`` has each 0-3 km panel,
mapped in (x, height m AGL).

.. code-block:: javascript

   {
     "format": "sounderpy-readout", "version": 1, "sounderpy": "3.3.0",
     "figure": "sounding",                      // "sounding", "winter", or "hodograph"
     "image": {"width": 2454, "height": 1533},
     "dark": false,
     "site": {"id": "KOKC", "name": "...", "source": "...", "valid": ["2026", "11", "23", "21"], "latlon": [35.4, -97.6]},
     "edited": null,
     "panels": {
       "skewt":     {"bbox": [x0, y0, x1, y1], "transform": [[a, b, c], [d, e, f]], "coords": ["T_degC", "log10_p_hPa"], "p_limits": [1050, 100]},
       "hodograph": {"bbox": [x0, y0, x1, y1], "transform": [[a, b, c], [d, e, f]], "coords": ["u_kt", "v_kt"], "offset": [0, 0], "max_height": 12000}
     },
     "units": {"p": "hPa", "T": "degC", ...},
     "profile": {"p": [...], "z": [...], "z_agl": [...], "T": [...], "Td": [...], "Tw": [...], "Tv": [...], "rh": [...],
                 "mr": [...], "theta": [...], "thetae": [...], "u": [...], "v": [...], "wdir": [...], "wspd": [...], "omega": [...]},
     "hodograph_trace": {"z_agl": [...], "u": [...], "v": [...]},
     "params": [
       {"key": "kinem.srh_0_to_3000", "name": "0-3 km storm-relative helicity", "value": 104.5, "units": "m²/s²",
        "about": "...", "bbox": [x0, y0, x1, y1], "layer": {"z": [0, 3000]}}
     ],
     "height_panels": [
       {"title": "Storm-relative wind", "bbox": [...], "transform": [[a, b, c], [d, e, f]], "z_limits": [0, 3000],
        "z": [0, 100, ...], "series": [{"label": "SR wind", "units": "kt", "values": [...]}]}
     ]
   }

A pixel ``(x, y)`` on the skew-t is the temperature and pressure that solve ``x = a·T + b·log10(p) + c`` and
``y = d·T + e·log10(p) + f``.

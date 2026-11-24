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
- **Layer**: drag a height range on the skew-t or a 0-3 km panel for its depth, lapse rate, mean RH and mixing ratio,
  precipitable water, θe range, mean wind, bulk shear, and SRH.
- **Measure**: drag between two points on the skew-t (ΔT, Δz, lapse rate, Δθ) or the hodograph (the wind difference).
- **Storm motion**: hover (or drag) on the hodograph to try a storm motion, with live SRH and storm-relative wind
  and the SRH area shaded; a click picks it, and your server can redraw with ``storm_motion=[direction, speed]``.
- Hovering also shows the 0-to-that-height bulk shear and SRH, with the figure's storm motion.
- **Parcel**: hover a level on the skew-t to lift a parcel from it: an approximate trace with CAPE and CIN shaded,
  and CAPE, CIN, LCL, LFC, EL, and lifted index (within a few percent of SounderPy's parcels). Once a parcel is pinned,
  hovering a height reads its CAPE from the LFC up to there, its share of the total, and the parcel's buoyancy.
  ``liftParcel("sb" | "ml" | "mu")`` lifts SounderPy's surface-based, 100 hPa mixed-layer, or most-unstable (lowest
  300 hPa) parcel.
- **Draw**: pen, arrows, and text on the plot; ``exportImage()`` gives the figure with the drawings on it.
- Optionally, a **panel** beside the image holds the full readout, tool results, and a standard levels table
  (surface, 925, 850, 700, 500, 300, 250 hPa), with a short readout at the pointer.

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
     fields: ["p", "z_agl", "T", "Td", "Tw", "rh", "mr", "thetae", "wind", "cumulative", "omega"],  // rows, in order
     temperature: "both",   // "C", "F", or "both"
     height: "m",           // "m" or "ft"
     wind: "kt",            // "kt", "m/s", or "mph"
     theme: "auto",         // "auto" follows the figure's dark mode; or "light", "dark"
     editable: false,       // start in edit mode
     tool: "inspect",       // what dragging does: "inspect", "edit", "layer", "measure", "storm", "parcel", or "draw"
     panel: null,           // an element for the full readout; the pointer then shows only tipFields
     tipFields: ["p", "z_agl", "T", "Td", "wind"],
     levels: [925, 850, 700, 500, 300, 250],  // the panel's levels table (after the surface), or false
     tipFooter: "",         // a line at the bottom of the pointer readout, e.g. "Click for more"
     cursor: true,          // add the temperature (skew-t) or wind (hodograph) under the pointer
     tapReadout: true,      // touch screens: a tap shows the readout (false leaves taps to the page unless a tool is on)
     parts: ["skewt", "hodograph", "panels", "params"],  // what responds to the pointer
     onHover: (level) => {},       // the level under the pointer, or null
     onEdit: (profile) => {},      // while a point is dragged
     onEditEnd: (profile) => {},   // when a drag ends
     onMark: (mark) => {},         // a layer or measurement's results (or null when cleared)
     onStormMotion: ({u, v, direction, speed}) => {},  // a picked storm motion: direction it comes from (°), speed (kt)
     onDraw: (drawings) => {},     // drawings added, undone, or cleared
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
   * - ``setTool(tool)``
     - ``"inspect"``, ``"edit"``, ``"layer"``, ``"measure"``, ``"storm"``, ``"parcel"``, or ``"draw"`` (``setEditing(on)`` is the same as ``"edit"`` / ``"inspect"``).
   * - ``setPanel(element)``
     - Puts the full readout in an element beside the image (``null`` for the pointer only).
   * - ``setUnits({temperature, height, wind})``
     - Changes the readout's units (the figure keeps its own).
   * - ``setDrawStyle({mode, color, width})``
     - How the draw tool draws: ``"pen"``, ``"arrow"``, or ``"text"``, a color, and a line width (image pixels).
   * - ``hasDrawing()``, ``undoDrawing()``, ``clearDrawing()``
     - Whether anything is drawn; remove the last drawing; remove them all.
   * - ``exportImage()``
     - A Promise of the figure with its drawings as a PNG Blob, at full size (the image must be same-origin or CORS).
   * - ``liftParcel(preset)``
     - Lifts and pins a ``"sb"``, ``"ml"``, or ``"mu"`` parcel (switching to the parcel tool).
   * - ``getMark()``, ``clearMark()``
     - The current layer or measurement (its heights or points, and ``rows`` of results), or null; remove it.
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

A picked storm motion works the same way: redraw with ``spy.build_sounding(clean_data, storm_motion=[direction, speed],
...)``, where ``direction`` is where the storm moves from, as ``onStormMotion`` gives it. Send the new image and readout
file back, then call ``readout.update(newReadout)``. The redrawn figure says
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
     "storm_motion": [u, v],                    // kt, ground relative; what the figure's SRH uses
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

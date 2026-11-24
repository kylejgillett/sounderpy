.. _readout_reference:

=====================
🧩 Readout Reference
=====================

Everything ``sounderpy-readout.js`` and SounderPy's readout functions take. **Experimental:** names and options may
change between releases. For how to use them, see :doc:`/developers/interactive_figures`.

***************************************************************

Python
======

``readout=True`` on :py:func:`spy.build_sounding` or :py:func:`spy.build_hodograph` (with ``save=True``) writes the
readout file next to the figure: ``kokc.png`` gets ``kokc.readout.json``. From the command line, ``--readout`` with
``--plot-file`` does the same (see :doc:`/reference/cli`).

.. py:function:: spy.readout_js_path()

   :return: the path of the installed ``sounderpy-readout.js``, to copy next to your page.
   :rtype: pathlib.Path

.. py:function:: spy.apply_profile_edits(clean_data, edits)

   Rebuild ``clean_data`` from a profile edited in the browser, to redraw with SounderPy.

   :param clean_data: the profile the readout was made from.
   :type clean_data: dict, required
   :param edits: the edited profile from ``getProfile()``: lists of ``p`` (hPa), ``z`` (m), ``T`` and ``Td`` (°C), and ``u`` and ``v`` (kt), one value per plotted level. A dewpoint above the temperature is lowered to it.
   :type edits: dict, required
   :return: the edited profile, with the original's site information and titles, other per-level values (such as omega) moved to the edited levels, and ``site_info['edited']`` describing the change (shown on the figure as **EDITED PROFILE**).
   :rtype: dict

A picked storm motion needs no helper: redraw with ``storm_motion=[direction, speed]``, as ``onStormMotion`` gives it.

***************************************************************

The script
==========

.. code-block:: javascript

   import { SounderPyReadout } from "./sounderpy-readout.js";

   // fetch the readout file and attach
   const readout = await SounderPyReadout.attach(img, "kokc.readout.json", options);

   // or, with the file already loaded
   const readout2 = new SounderPyReadout(img, readoutJson, options);

It's an ES module; ``window.SounderPyReadout`` is also set. ``img`` is the figure's ``<img>``: the script adds its
layer and tooltips next to it, in the same parent element.

Options
-------

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 22 26 52

   * - Option
     - Default
     - What it does
   * - ``fields``
     - ``["p", "z_agl", "T", "Td", "Tw", "rh", "mr", "thetae", "wind", "cumulative", "omega"]``
     - The readout's rows, in order (see `Readout fields`_).
   * - ``temperature``
     - ``"both"``
     - ``"C"``, ``"F"``, or ``"both"``.
   * - ``height``
     - ``"m"``
     - ``"m"`` or ``"ft"``.
   * - ``wind``
     - ``"kt"``
     - ``"kt"``, ``"m/s"``, or ``"mph"``.
   * - ``theme``
     - ``"auto"``
     - Tooltip colors: ``"auto"`` follows the figure's dark mode; or ``"light"``, ``"dark"``.
   * - ``tool``
     - ``"inspect"``
     - What dragging does (see `Tools`_).
   * - ``editable``
     - ``false``
     - Start with the edit tool (the same as ``tool: "edit"``).
   * - ``parts``
     - ``["skewt", "hodograph", "panels", "params"]``
     - What responds to the pointer: the skew-T, hodograph, 0-3 km panels, and table values.
   * - ``panel``
     - none
     - An element for the full readout, tool results, and the levels table.
   * - ``panelTip``
     - ``true``
     - With a panel, also show a short readout at the pointer.
   * - ``tipFields``
     - ``["p", "z_agl", "T", "Td", "wind"]``
     - The pointer readout's rows when there's a panel.
   * - ``panelHint``
     - ``"Hover the plot to read it here."``
     - What the panel says before anything is read.
   * - ``levels``
     - ``[925, 850, 700, 500, 300, 250]``
     - The panel's levels table (after the surface), or ``false`` for none.
   * - ``cursor``
     - ``true``
     - Add the temperature (skew-T) or wind (hodograph) under the pointer.
   * - ``tapReadout``
     - ``true``
     - On touch screens, a tap shows the readout; ``false`` leaves taps to the page unless a tool is on.
   * - ``tipFooter``
     - none
     - A line at the bottom of the pointer readout, e.g. ``"Click for more"``.
   * - ``clouds``
     - ``false``
     - Shade cloud layers on the skew-T (RH ≥ 90%, over ice below 0 °C, at least 200 m deep), and list them in the
       panel.
   * - ``keys``
     - ``true``
     - Up and down arrow keys step the readout (or an unpinned parcel or precipitation source) one level while the
       pointer is on the plot; ``false`` leaves arrow keys to the page.

Callbacks
---------

Also options:

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 34 66

   * - Callback
     - Called with
   * - ``onHover(level)``
     - The readout values under the pointer, or ``null``.
   * - ``onEdit(profile)``
     - The edited profile, while a point is dragged (as ``getProfile()``).
   * - ``onEditEnd(profile)``
     - The edited profile, when a drag ends (and after ``undo()`` or ``resetEdits()``).
   * - ``onMark(mark)``
     - A layer, measurement, parcel, or storm motion's results (as ``getMark()``), or ``null`` when cleared.
   * - ``onStormMotion(motion)``
     - A picked storm motion: ``{u, v, direction, speed}``, with ``direction`` where it comes from (°) and ``speed`` in kt, whole numbers, as ``storm_motion=[direction, speed]`` takes them.
   * - ``onDraw(drawings)``
     - The drawings, when one is added, undone, or cleared.

Methods
-------

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 34 66

   * - Method
     - What it does
   * - ``update(readout)``
     - Use a new readout file, e.g. after a redraw (set the ``<img>``'s new ``src`` too). Clears edits and marks.
   * - ``setTool(tool)``
     - Change what dragging does (see `Tools`_); clears a mark. ``setEditing(on)`` is ``setTool("edit")`` or ``setTool("inspect")``.
   * - ``getMark()``, ``clearMark()``
     - The current layer, measurement, parcel, or storm motion and its ``rows`` of results, or ``null``; remove it.
   * - ``liftParcel(preset)``
     - Lift and pin a ``"sb"``, ``"ml"`` (lowest 100 hPa), or ``"mu"`` (lowest 300 hPa) parcel; returns its results.
   * - ``setPanel(element, {tip})``
     - Move the full readout to an element (``null`` for the pointer only); ``tip: false`` drops the pointer readout.
   * - ``setUnits({temperature, height, wind})``
     - Change the readout's units.
   * - ``setClouds(on)``, ``getClouds()``
     - Shade cloud layers or not; the cloud layers: ``bottom`` and ``top`` (m AGL), ``bottom_p`` and ``top_p`` (hPa).
   * - ``setEnabled(on)``
     - Show or hide the readout layer.
   * - ``getProfile()``
     - The (edited) profile: ``p`` (hPa), ``z`` (m MSL), ``T`` and ``Td`` (°C), ``u`` and ``v`` (kt), one value per level.
   * - ``isEdited()``, ``undo()``, ``resetEdits()``
     - Whether the profile was edited; undo the last drag; go back to the figure's profile.
   * - ``levelAt(p)``
     - The readout values at a pressure (hPa).
   * - ``setDrawStyle({mode, color, width})``
     - How the draw tool draws: ``"pen"``, ``"arrow"``, or ``"text"``; a CSS color; a line width in image pixels.
   * - ``hasDrawing()``, ``undoDrawing()``, ``clearDrawing()``
     - Whether anything is drawn; remove the last drawing; remove them all.
   * - ``exportImage(type)``
     - A Promise of the figure with its drawings, as an image Blob at full size (``"image/png"`` by default). The image must be same-origin or served with CORS.
   * - ``destroy()``
     - Remove the layer, tooltips, and listeners.

Tools
-----

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 16 84

   * - Tool
     - Dragging (or a click)
   * - ``"inspect"``
     - Nothing: the pointer reads the plot.
   * - ``"layer"``
     - Up or down on the skew-T or a 0-3 km panel: the layer's depth (also in dam), lapse rate (virtual temperature), mean RH and mixing ratio, precipitable water, θe range, mean wind, bulk shear, and SRH. On the skew-T, an end near 1000, 925, 850, 700, 500, 300, 250, or 200 hPa snaps to it (and is tagged); with both ends snapped, the results include that thickness. A click clears it.
   * - ``"measure"``
     - Between two points on the skew-T (ΔT, Δz, lapse rate, Δθ, Δp) or the hodograph (the two winds and their difference). A click clears it.
   * - ``"storm"``
     - Hover the hodograph to try a storm motion (SRH, storm-relative wind, the SRH area shaded); a click picks it (``onStormMotion``).
   * - ``"parcel"``
     - Hover a level on the skew-T to lift a parcel (trace, CAPE, CIN, LCL, LFC, EL, lifted index); a click pins it, and hovering then reads CAPE from the LFC to the pointer.
   * - ``"ptype"``
     - Hover a level on the skew-T to follow precipitation down from it: whether it's saturated, what it starts as (from the cloud top's temperature), each melting and refreezing layer's energy and what it does, and the type at the surface; a click pins it. The same method as SounderPy's winter soundings.
   * - ``"draw"``
     - Draw with ``setDrawStyle``'s mode: a pen stroke, an arrow, or (click) a text label.
   * - ``"edit"``
     - Move temperature and dewpoint points along their pressure level, or hodograph points anywhere (``getProfile()``).

On touch screens, these tools take one-finger drags on the plot, so a page that also pans the image should pan with
two fingers.

Readout fields
--------------

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 18 82

   * - Field
     - Row
   * - ``p``, ``z_agl``, ``z``
     - Pressure; height above ground; height above sea level.
   * - ``T``, ``Td``, ``Tw``, ``Tv``
     - Temperature, dewpoint, wet-bulb, and virtual temperature.
   * - ``rh``, ``mr``
     - Relative humidity and mixing ratio.
   * - ``theta``, ``thetae``
     - Potential and equivalent potential temperature.
   * - ``wind``
     - Wind direction (where it comes from) and speed.
   * - ``cumulative``
     - Bulk shear and SRH from the surface to that height (within the hodograph's 12 km).
   * - ``omega``
     - Vertical velocity, when the data has it.

Styling
-------

Set these CSS variables on the page (or any parent of the image) to restyle the readouts:

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 34 66

   * - Variable
     - What it styles
   * - ``--spy-readout-font``
     - The readouts' font.
   * - ``--spy-readout-label``
     - Field and parameter names.
   * - ``--spy-readout-value``
     - Values.
   * - ``--spy-readout-bg``, ``--spy-readout-text``, ``--spy-readout-muted``, ``--spy-readout-border``
     - The pointer readout and result cards (set them on ``.spy-readout-tip`` and ``.spy-readout-tip[data-dark]``).

The elements are ``.spy-readout-layer`` (the canvas over the image), ``.spy-readout-tip`` (the pointer readout),
``.spy-readout-card`` (a tool's results without a panel), and ``.spy-readout-panel`` (added to your panel element).

.. _interactive_figures:

=================================
🖱️ Interactive Figures on the Web
=================================

Turn a SounderPy PNG into an interactive sounding on your own web page: readers hover it for the data, read the
parameter tables, measure layers, try storm motions, lift parcels, draw on it, and edit the profile. SounderPy saves a
small readout file next to the figure, and ``sounderpy-readout.js`` (shipped with SounderPy, no dependencies) reads it.
The figure is still SounderPy's own PNG; nothing is redrawn in the browser.

**Experimental:** the readout file and the script's options may change between releases. The
`SounderPy web app <https://sounderpy.app/>`_ uses everything on this page in its Interactive Sounding Viewer.

This page shows how to set it up. For every option, see :doc:`/developers/readout_reference`; for the JSON itself,
:doc:`/developers/readout_file`.

***************************************************************

Quickstart
==========

**1. Save a figure with its readout, and copy the script**

.. code-block:: python

   import shutil
   from pathlib import Path

   import sounderpy as spy

   site = Path('site')
   site.mkdir(exist_ok=True)

   clean_data = spy.get_bufkit_data('hrrr', 'kokc', 3)
   spy.build_sounding(clean_data, save=True, filename=str(site / 'kokc.png'), readout=True)
   # site/kokc.png and site/kokc.readout.json

   shutil.copy(spy.readout_js_path(), site)   # site/sounderpy-readout.js

From the command line, add ``--readout`` to a ``--plot-file`` (see :doc:`/reference/cli`). ``build_hodograph()`` and
winter soundings (``mode='winter'``) work the same way.

**2. Make the page:** ``site/index.html``

.. code-block:: html

   <!doctype html>
   <html lang="en">
   <head>
     <meta charset="utf-8">
     <meta name="viewport" content="width=device-width, initial-scale=1">
     <title>KOKC sounding</title>
     <style>
       .plot { max-width: 1200px; margin: 0 auto; }
       .plot img { display: block; width: 100%; height: auto; }
     </style>
   </head>
   <body>
     <div class="plot"><img id="sounding" src="kokc.png" alt="KOKC sounding"></div>

     <script type="module">
       import { SounderPyReadout } from "./sounderpy-readout.js";

       const img = document.getElementById("sounding");
       const readout = await SounderPyReadout.attach(img, "kokc.readout.json");
     </script>
   </body>
   </html>

**3. Open it over HTTP**

.. code-block:: bash

   cd site && python -m http.server 8000

Then open http://localhost:8000 and hover the skew-T, the hodograph, a 0-3 km panel, or a value in the tables.
Browsers don't load module scripts from ``file://`` pages, so open it through a server, not by double-clicking.

***************************************************************

Recipes
=======

Each recipe adds to the quickstart page.

A readout panel beside the plot
-------------------------------

Give the full readout its own panel, with a short one at the pointer. The panel also holds tool results and a table
of standard levels (surface, 925, 850, 700, 500, 300, and 250 hPa).

.. code-block:: html

   <style>
     .layout { display: grid; grid-template-columns: 280px 1fr; gap: 12px; }
     .readout { height: 80vh; overflow: auto; padding: 12px; border: 1px solid #ccc; border-radius: 8px; }
   </style>

   <div class="layout">
     <aside id="readout" class="readout"></aside>
     <div class="plot"><img id="sounding" src="kokc.png" alt="KOKC sounding"></div>
   </div>

.. code-block:: javascript

   const readout = await SounderPyReadout.attach(img, "kokc.readout.json", {
     panel: document.getElementById("readout"),
   });

The pointer shows ``tipFields`` (pressure, height, T, Td, and wind by default); the panel shows everything, with
standard thicknesses (1000–500, 1000–850, and 850–700 hPa) under the levels table. With the pointer on the plot, the
up and down arrow keys step the readout one level at a time.

On phones
---------

On a small screen, put the panel under the plot and leave the plot clear of tooltips, so a finger can move over it
and read the panel:

.. code-block:: css

   @media (max-width: 700px) {
     .layout { grid-template-columns: 1fr; }
     .readout { order: 2; height: 30vh; }
   }

.. code-block:: javascript

   const phone = window.matchMedia("(max-width: 700px)");
   const place = () => readout.setPanel(document.getElementById("readout"), { tip: !phone.matches });
   phone.addEventListener("change", place);
   place();

A touch reads the plot. Dragging keeps reading sideways, but an up-and-down drag scrolls the page; to read in
every direction, give the plot's container ``touch-action: none`` (the page then scrolls from outside the plot). A
fixed panel height keeps the plot from moving as the readout changes.

Tool buttons
------------

The tools change what dragging does. Results (a layer's lapse rate and SRH, a parcel's CAPE, and so on) go to the
panel, or a card by the pointer without one, and to ``onMark``:

.. code-block:: html

   <div id="tools">
     <button data-tool="inspect">Inspect</button>
     <button data-tool="layer">Layer</button>
     <button data-tool="measure">Measure</button>
     <button data-tool="storm">Storm motion</button>
     <button data-tool="parcel">Parcel</button>
     <button data-tool="ptype">Precip type</button>
     <button data-tool="draw">Draw</button>
     <button id="clear">Clear</button>
     <label><input type="checkbox" id="clouds"> Clouds</label>
   </div>

.. code-block:: javascript

   const readout = await SounderPyReadout.attach(img, "kokc.readout.json", {
     onMark: (mark) => console.log(mark?.rows),   // [label, value] pairs, e.g. ["SRH", "104 m²/s²"]
   });

   document.getElementById("tools").addEventListener("click", (e) => {
     if (e.target.dataset.tool) readout.setTool(e.target.dataset.tool);
   });
   document.getElementById("clear").onclick = () => readout.clearMark();

   readout.liftParcel("ml");   // lift and pin SounderPy's mixed-layer parcel (or "sb", "mu")

   // shade cloud layers (RH ≥ 90%) on the skew-T
   document.getElementById("clouds").onchange = (e) => readout.setClouds(e.target.checked);

- **Layer:** drag up or down on the skew-T or a 0-3 km panel. Its ends snap to standard levels, so dragging from
  1000 to 500 hPa gives that thickness.
- **Measure:** drag between two points on the skew-T or hodograph.
- **Storm motion:** hover the hodograph to preview SRH and storm-relative wind; a click calls ``onStormMotion``
  (see the server recipe below).
- **Parcel:** hover a level to lift a parcel from it; a click pins it, and hovering then reads CAPE from the LFC up
  to the pointer.
- **Precip type:** hover a level to follow precipitation down from it, through melting and refreezing layers, to
  its type at the surface; a click pins it.
- **Draw:** pen, arrows, and text; ``setDrawStyle({ mode: "arrow", color: "#2563eb" })`` changes how it draws.

The browser's parcels, precipitation types, SRH, and lapse rates are quick approximations; the figure's tables are
SounderPy's exact values.

Units and colors
----------------

.. code-block:: javascript

   readout.setUnits({ temperature: "F", height: "ft", wind: "mph" });   // or set them as options

The readouts take their font and colors from CSS variables, so they can match your site:

.. code-block:: css

   body {
     --spy-readout-font: "Inter", system-ui, sans-serif;
     --spy-readout-label: #c89b3c;   /* parameter names */
     --spy-readout-value: #fff;      /* values */
   }
   .spy-readout-tip, .spy-readout-tip[data-dark] {
     --spy-readout-bg: rgb(30 30 30 / 95%);
     --spy-readout-text: #eee;
     --spy-readout-border: #555;
   }

Redraw on a server
------------------

Edits and picked storm motions are recalculated by SounderPy, so they need a small server. This one, with
`FastAPI <https://fastapi.tiangolo.com/>`_, serves the page and redraws on request:

.. code-block:: python

   # server.py: run with  uvicorn server:app --port 8000
   import base64
   import json
   import tempfile
   import threading
   from functools import lru_cache
   from pathlib import Path
   from typing import Optional

   import matplotlib
   matplotlib.use('Agg')
   from fastapi import FastAPI
   from fastapi.staticfiles import StaticFiles
   from pydantic import BaseModel

   import sounderpy as spy
   from sounderpy.readout import readout_filename

   app = FastAPI()
   PLOT_LOCK = threading.Lock()   # Matplotlib draws one figure at a time


   @lru_cache(maxsize=32)
   def profile(model, site, hour):
       return spy.get_bufkit_data(model, site, hour)


   class Redraw(BaseModel):
       model: str
       site: str
       hour: int
       edits: Optional[dict] = None                # readout.getProfile()
       storm_motion: Optional[list[float]] = None  # [direction, speed], from onStormMotion


   @app.post('/redraw')
   def redraw(req: Redraw):
       data = profile(req.model, req.site, req.hour)
       if req.edits:
           data = spy.apply_profile_edits(data, req.edits)
       with tempfile.TemporaryDirectory() as tmp, PLOT_LOCK:
           png = Path(tmp) / 'plot.png'
           spy.build_sounding(data, save=True, filename=str(png), readout=True, map_zoom=0,
                              storm_motion=req.storm_motion or 'right_moving')
           return {'image': 'data:image/png;base64,' + base64.b64encode(png.read_bytes()).decode(),
                   'readout': json.loads(readout_filename(png).read_text())}


   app.mount('/', StaticFiles(directory='site', html=True), name='site')

And on the page:

.. code-block:: javascript

   async function redraw(change) {
     const res = await fetch("/redraw", {
       method: "POST",
       headers: { "Content-Type": "application/json" },
       body: JSON.stringify({ model: "hrrr", site: "kokc", hour: 3, ...change }),
     });
     const { image, readout: file } = await res.json();
     img.src = image;
     readout.update(file);
   }

   const readout = await SounderPyReadout.attach(img, "kokc.readout.json", {
     onStormMotion: ({ direction, speed }) => redraw({ storm_motion: [direction, speed] }),
   });

   // with the Edit tool on (readout.setTool("edit")), drag points, then:
   document.getElementById("redraw").onclick = () => redraw({ edits: readout.getProfile() });

A redrawn edited profile says **EDITED PROFILE** on the figure, with what changed, so it can't pass as the
original data. A picked storm motion shows as ``SM: USER DEFINED``.

A whole model run
-----------------

For a model site, write every forecast hour's figure and readout once, from one BUFKIT download, and serve them as
static files:

.. code-block:: python

   from pathlib import Path

   import sounderpy as spy

   out = Path('site/hrrr/kokc')
   out.mkdir(parents=True, exist_ok=True)

   for data in spy.get_bufkit_data('hrrr', 'kokc', 0, all_fcst_hours=True):
       hour = data['site_info']['fcst-hour']   # 'F00', 'F01', ...
       spy.build_sounding(data, save=True, filename=str(out / f'{hour}.png'), readout=True, map_zoom=0)

A readout file is about 30 KB (a few KB compressed), and everything but redraws runs in the reader's browser.

Save the drawing
----------------

.. code-block:: javascript

   const blob = await readout.exportImage();   // the figure with its drawings, full size
   const link = Object.assign(document.createElement("a"), {
     href: URL.createObjectURL(blob), download: "sounding_annotated.png",
   });
   link.click();

***************************************************************

Full example
============

One page with everything above that runs in the browser: a readout panel (under the plot on phones), tool buttons,
SounderPy's parcels, a units switch, and saving the drawing. Copy the two files, run them, and change what you need.

.. code-block:: text

   build_site.py
   site/
     index.html               (below)
     kokc.png                 (written by build_site.py)
     kokc.readout.json        (written by build_site.py)
     sounderpy-readout.js     (copied by build_site.py)

**build_site.py**

.. code-block:: python

   import shutil
   from pathlib import Path

   import sounderpy as spy

   site = Path('site')
   site.mkdir(exist_ok=True)

   clean_data = spy.get_bufkit_data('hrrr', 'kokc', 3)
   spy.build_sounding(clean_data, save=True, filename=str(site / 'kokc.png'), readout=True)
   shutil.copy(spy.readout_js_path(), site)

**site/index.html**

.. code-block:: html

   <!doctype html>
   <html lang="en">
   <head>
     <meta charset="utf-8">
     <meta name="viewport" content="width=device-width, initial-scale=1">
     <title>KOKC sounding</title>
     <style>
       body { margin: 0; font: 15px/1.4 system-ui, sans-serif; background: #f4f4f5; color: #111; }
       h1 { margin: 0; padding: 12px 16px; font-size: 20px; }

       .toolbar { display: flex; flex-wrap: wrap; gap: 6px; padding: 0 16px 12px; }
       .toolbar button, .toolbar select {
         padding: 6px 10px; border: 1px solid #bbb; border-radius: 6px; background: #fff; font: inherit; cursor: pointer;
       }
       .toolbar button[aria-pressed="true"] { background: #111; color: #fff; border-color: #111; }
       .toolbar label { display: flex; align-items: center; gap: 4px; }

       .layout { display: grid; grid-template-columns: 300px 1fr; gap: 12px; padding: 0 16px 16px; }
       .readout { height: 80vh; overflow: auto; padding: 12px; background: #fff; border: 1px solid #ddd; border-radius: 8px; }
       .plot img { display: block; width: 100%; height: auto; }

       /* phones: the panel goes under the plot */
       @media (max-width: 760px) {
         .layout { grid-template-columns: 1fr; }
         .readout { order: 2; height: 35vh; }
       }
     </style>
   </head>
   <body>
     <h1>KOKC sounding</h1>

     <div class="toolbar">
       <button data-tool="inspect" aria-pressed="true">Inspect</button>
       <button data-tool="layer">Layer</button>
       <button data-tool="measure">Measure</button>
       <button data-tool="storm">Storm motion</button>
       <button data-tool="parcel">Parcel</button>
       <button data-tool="ptype">Precip type</button>
       <button data-tool="draw">Draw</button>
       <select id="parcel" aria-label="Lift a parcel">
         <option value="">Lift a parcel…</option>
         <option value="sb">Surface-based</option>
         <option value="ml">Mixed-layer</option>
         <option value="mu">Most unstable</option>
       </select>
       <select id="units" aria-label="Units">
         <option value="metric">°C, m, kt</option>
         <option value="imperial">°F, ft, mph</option>
       </select>
       <label><input type="checkbox" id="clouds"> Clouds</label>
       <button id="clear">Clear</button>
       <button id="save">Save image</button>
     </div>

     <div class="layout">
       <aside id="readout" class="readout"></aside>
       <div class="plot"><img id="sounding" src="kokc.png" alt="KOKC sounding"></div>
     </div>

     <script type="module">
       import { SounderPyReadout } from "./sounderpy-readout.js";

       const img = document.getElementById("sounding");
       const panel = document.getElementById("readout");
       const toolButtons = [...document.querySelectorAll("[data-tool]")];

       const readout = await SounderPyReadout.attach(img, "kokc.readout.json", {
         panel,
         temperature: "C",
         onStormMotion: ({ direction, speed }) => console.log(`storm motion: ${direction}° at ${speed} kt`),
       });

       // tool buttons: one is pressed at a time
       const press = (tool) => toolButtons.forEach((b) => b.setAttribute("aria-pressed", b.dataset.tool === tool));
       toolButtons.forEach((b) => b.addEventListener("click", () => {
         readout.setTool(b.dataset.tool);
         press(b.dataset.tool);
       }));

       // SounderPy's parcels (this switches to the parcel tool)
       document.getElementById("parcel").addEventListener("change", (e) => {
         if (!e.target.value) return;
         readout.liftParcel(e.target.value);
         press("parcel");
         e.target.value = "";
       });

       document.getElementById("units").addEventListener("change", (e) => {
         readout.setUnits(e.target.value === "imperial"
           ? { temperature: "F", height: "ft", wind: "mph" }
           : { temperature: "C", height: "m", wind: "kt" });
       });

       document.getElementById("clouds").addEventListener("change", (e) => readout.setClouds(e.target.checked));

       document.getElementById("clear").addEventListener("click", () => {
         readout.clearMark();
         readout.clearDrawing();
       });

       document.getElementById("save").addEventListener("click", async () => {
         const blob = await readout.exportImage();
         const link = Object.assign(document.createElement("a"), {
           href: URL.createObjectURL(blob), download: "kokc_annotated.png",
         });
         link.click();
       });

       // phones: no readout at the finger, so it doesn't cover the plot
       const phone = window.matchMedia("(max-width: 760px)");
       const place = () => readout.setPanel(panel, { tip: !phone.matches });
       phone.addEventListener("change", place);
       place();
     </script>
   </body>
   </html>

**Run it**

.. code-block:: bash

   python build_site.py
   cd site && python -m http.server 8000

Then open http://localhost:8000. To publish it, upload the ``site`` folder to any web host; it needs no server code.
To also redraw edited profiles and picked storm motions with SounderPy, add the server from
`Redraw on a server`_ and an **Edit** button (``data-tool="edit"``).

***************************************************************

Common problems
===============

.. list-table::
   :class: param-table
   :header-rows: 1
   :widths: 40 60

   * - Problem
     - Fix
   * - Nothing happens, and the console mentions modules or CORS
     - Open the page over HTTP (``python -m http.server``), not as a ``file://`` page.
   * - The readout is offset from the plot
     - The readout file isn't from this PNG: save both together, and don't crop the image. Showing it at any size is
       fine.
   * - ``exportImage()`` fails
     - The image must come from the page's own site, or be served with CORS headers.
   * - "file version … is newer than this script" in the console
     - Copy the script from the same SounderPy that wrote the file (``spy.readout_js_path()``).
   * - No readout file was written
     - Readouts need ``save=True``. They work with ``build_sounding()`` (severe and winter) and ``build_hodograph()``;
       composites, VAD hodographs, and publication figures don't write them.

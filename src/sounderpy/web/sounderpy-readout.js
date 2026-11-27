/*!
 * SounderPy readout — hover readouts and profile editing for SounderPy figures on the web.
 *
 * EXPERIMENTAL: the readout file format and this script's options may change between releases.
 *
 * Reads the readout file SounderPy writes next to a figure (build_sounding(..., save=True, readout=True)
 * writes kokc.png and kokc.readout.json) and lays a transparent layer over the <img>:
 *   - hover the skew-t or hodograph to read the profile at that level, linked across both panels
 *   - optionally drag temperature, dewpoint, and wind points to edit the profile, then send
 *     getProfile() to a server that redraws it with sounderpy.apply_profile_edits()
 *   - measure a layer (drag a height range: lapse rate, moisture, shear, SRH, ...) or the difference between two
 *     points on the skew-t or hodograph
 *   - try a storm motion on the hodograph with live SRH, then send it to a server to redraw with it
 *   - lift a parcel from any level for an approximate trace, CAPE, CIN, LCL, LFC, and EL
 *   - follow precipitation down from any level for its type at the surface (melting and refreezing layers)
 *   - shade cloud layers, and read thicknesses; arrow keys step the readout up and down the profile
 *   - draw on the plot (pen, arrows, text) and export the annotated image
 *   - optionally put the full readout in a panel beside the image, with a short one at the pointer
 *
 * No dependencies. Works with a responsive or scaled image, in light or dark mode, with mouse or touch.
 *
 *   import { SounderPyReadout } from "./sounderpy-readout.js";
 *   const readout = await SounderPyReadout.attach(img, "kokc.readout.json");
 *
 * Docs: https://sounderpy.app/docs/developers/interactive_figures.html
 * Part of SounderPy (https://github.com/kylejgillett/sounderpy), MIT License, (C) Kyle J Gillett 2026
 */

const FORMAT = "sounderpy-readout";
const VERSION = 1;

const DEFAULT_FIELDS = ["p", "z_agl", "T", "Td", "Tw", "rh", "mr", "thetae", "wind", "cumulative", "omega"];
const TOOLS = ["inspect", "edit", "layer", "measure", "storm", "parcel", "ptype", "draw"];
// tools that mark the plot (their results replace the pointer readout)
const MARK_TOOLS = ["layer", "measure", "storm", "parcel", "ptype"];
// the layer tool's ends snap to these levels (hPa), for thicknesses
const SNAP_LEVELS = [1000, 925, 850, 700, 500, 300, 250, 200];
// the panel's thicknesses (hPa)
const THICKNESSES = [[1000, 500], [1000, 850], [850, 700]];
// precipitation types: SounderPy's codes, names, and the winter figure's colors
const PTYPES = {
  SN: ["Snow", "#1f6fd1"], PL: ["Sleet", "#9b4fd1"], FZRA: ["Freezing rain", "#d1328a"], "FZRA/PL": ["Freezing rain / sleet", "#c0409f"],
  FZDZ: ["Freezing drizzle", "#e070b0"], RA: ["Rain", "#2e9b47"], "RA/SN": ["Rain / snow", "#3aa3a3"], "RA/PL": ["Rain / sleet", "#6c7fbf"],
  DZ: ["Drizzle", "#5cae6c"],
};
const LEVELS = [925, 850, 700, 500, 300, 250]; // the standard levels table, below the surface
const TIP_FIELDS = ["p", "z_agl", "T", "Td", "wind"]; // the pointer readout when a panel has the rest
const KT = 0.514444; // m/s per kt
const MARK = "#0aa5a5"; // layer and measure marks

const LABELS = {
  p: "Pressure", z_agl: "Height", z: "Height MSL", T: "Temp", Td: "Dewpoint", Tw: "Wet-bulb", Tv: "Virtual T",
  rh: "RH", mr: "Mix ratio", theta: "θ", thetae: "θe", wind: "Wind", omega: "ω",
};

const STYLE_ID = "spy-readout-style";
// Pages can restyle the readouts with CSS variables on any parent: --spy-readout-font, --spy-readout-label (parameter
// names), --spy-readout-value (values), and on .spy-readout-tip --spy-readout-bg, --spy-readout-text,
// --spy-readout-muted, and --spy-readout-border.
const STYLE = `
.spy-readout-layer { position: absolute; z-index: 2; cursor: crosshair; touch-action: pan-y pinch-zoom; }
.spy-readout-layer:is([data-editing], [data-marking]) { touch-action: none; }
.spy-readout-tip {
  --spy-readout-bg: rgba(255, 255, 255, 0.96); --spy-readout-text: #111; --spy-readout-muted: #555;
  --spy-readout-border: rgba(0, 0, 0, 0.35); --spy-readout-edited: #c2187a;
  position: absolute; z-index: 3; display: none; pointer-events: none; box-sizing: border-box;
  min-width: min(190px, calc(100vw - 16px)); max-width: calc(100vw - 16px); padding: 6px 10px; border: 1px solid var(--spy-readout-border); border-radius: 6px;
  background: var(--spy-readout-bg); color: var(--spy-readout-text);
  font: 12px/1.45 var(--spy-readout-font, ui-monospace, SFMono-Regular, Menlo, Consolas, monospace); font-variant-numeric: tabular-nums;
  box-shadow: 0 2px 10px rgba(0, 0, 0, 0.25); text-align: left; white-space: nowrap;
}
.spy-readout-tip[data-dark] {
  --spy-readout-bg: rgba(24, 26, 30, 0.95); --spy-readout-text: #eee; --spy-readout-muted: #aaa;
  --spy-readout-border: rgba(255, 255, 255, 0.3); --spy-readout-edited: #ff6ec7;
}
.spy-readout-tip div { display: flex; justify-content: space-between; gap: 14px; }
.spy-readout-tip span { color: var(--spy-readout-label, var(--spy-readout-muted)); }
.spy-readout-tip b { color: var(--spy-readout-value, inherit); font-weight: 600; }
.spy-readout-tip .spy-readout-head { margin: 0 0 2px; color: var(--spy-readout-edited); font-weight: 600; }
.spy-readout-tip .spy-readout-title { margin: 0 0 2px; font-weight: 700; }
.spy-readout-tip .spy-readout-foot { margin: 5px 0 0; padding-top: 4px; border-top: 1px solid var(--spy-readout-border); color: var(--spy-readout-muted); font: 600 11.5px/1.4 var(--spy-readout-font, system-ui, sans-serif); white-space: normal; }
.spy-readout-tip .spy-readout-about { max-width: 280px; margin: 4px 0 0; color: var(--spy-readout-muted); white-space: normal; font-family: var(--spy-readout-font, system-ui, sans-serif); }
.spy-readout-tip hr { margin: 4px 0; border: 0; border-top: 1px solid var(--spy-readout-border); }
.spy-readout-card { --spy-readout-edited: #0a8f8f; border-color: #0aa5a5; border-width: 2px; }
.spy-readout-card[data-dark] { --spy-readout-edited: #3fd6d6; }
.spy-readout-panel { font: 13px/1.5 var(--spy-readout-font, ui-monospace, SFMono-Regular, Menlo, Consolas, monospace); font-variant-numeric: tabular-nums; }
.spy-readout-panel div { display: flex; justify-content: space-between; gap: 12px; }
.spy-readout-panel span { flex: none; color: var(--spy-readout-label, color-mix(in srgb, currentColor 70%, transparent)); white-space: nowrap; }
.spy-readout-panel b { color: var(--spy-readout-value, inherit); font-weight: 600; text-align: right; }
.spy-readout-panel section[hidden] { display: none; }
.spy-readout-panel section:not([hidden]) ~ section { margin-top: 12px; padding-top: 12px; border-top: 1px solid rgba(128, 128, 128, .4); }
.spy-readout-panel p { margin: 0; }
.spy-readout-panel .spy-readout-title { margin: 0 0 4px; font-weight: 700; }
.spy-readout-panel .spy-readout-sub { margin-top: 10px; }
.spy-readout-panel .spy-readout-head { margin: 0 0 4px; color: #0aa5a5; font-weight: 600; }
.spy-readout-panel .spy-readout-about { margin-top: 6px; opacity: .75; white-space: normal; font-family: var(--spy-readout-font, system-ui, sans-serif); }
.spy-readout-panel .spy-readout-empty { opacity: .6; font-family: var(--spy-readout-font, system-ui, sans-serif); }
.spy-readout-panel hr { margin: 6px 0; border: 0; border-top: 1px solid rgba(128, 128, 128, .4); }
.spy-readout-panel table { width: 100%; border-collapse: collapse; font-size: 12px; font-variant-numeric: tabular-nums; }
.spy-readout-panel th { padding: 0 0 2px; color: var(--spy-readout-label, color-mix(in srgb, currentColor 70%, transparent)); font-weight: 600; text-align: right; }
.spy-readout-panel td { padding: 1px 0; color: var(--spy-readout-value, inherit); text-align: right; }
.spy-readout-panel :is(th, td):first-child { text-align: left; }
.spy-readout-layer[data-readout-tool="draw"] { cursor: crosshair; }
.spy-readout-text-input { position: absolute; z-index: 4; box-sizing: border-box; min-height: 30px; padding: 4px 6px; border: 2px solid #e11d48; border-radius: 4px; background: #fff; color: #111; font: 600 15px/1.3 system-ui, sans-serif; resize: none; overflow: hidden; white-space: pre-wrap; overflow-wrap: anywhere; }
`;

// ---------------------------------------------------------------------------------------------------
// thermodynamics for live values while editing (the server's redraw replaces them with SounderPy's)
// ---------------------------------------------------------------------------------------------------
const es = (T) => 6.112 * Math.exp((17.67 * T) / (T + 243.5)); // hPa, Bolton (1980)
const esIce = (T) => 6.112 * Math.exp((22.46 * T) / (272.62 + T)); // hPa, over ice
const mixingRatio = (p, Td) => (622 * es(Td)) / (p - es(Td)); // g/kg
// relative humidity (%) for cloud: over ice below 0 °C, as SounderPy's winter soundings find clouds
const rhCloud = (T, Td) => (T == null || Td == null ? null : (100 * es(Td)) / (T < 0 ? esIce(T) : es(T)));
function wetBulb(p, T, Td) {
  // psychrometric equation, solved by bisection
  const e = es(Td);
  let lo = Td, hi = T;
  for (let i = 0; i < 40; i++) {
    const tw = (lo + hi) / 2;
    const f = es(tw) - 0.00066 * (1 + 0.00115 * tw) * p * (T - tw) - e;
    if (f > 0) hi = tw; else lo = tw;
  }
  return (lo + hi) / 2;
}
function thetaE(p, T, Td) {
  // Bolton (1980) equation 43
  const r = mixingRatio(p, Td), tk = T + 273.15, tdk = Td + 273.15;
  const tl = 1 / (1 / (tdk - 56) + Math.log(tk / tdk) / 800) + 56;
  const th = tk * Math.pow(1000 / p, 0.2854 * (1 - 0.00028 * r));
  return th * Math.exp((3.376 / tl - 0.00254) * r * (1 + 0.00081 * r));
}

// a number for display, without a "-0"
const num = (v, digits = 0) => v.toFixed(digits).replace(/^-(0(\.0+)?)$/, "$1");
// parcels: an approximate surface-style lift (Bolton 1980 LCL, a dry adiabat, then a pseudoadiabat), in virtual
// temperature, for a quick look; SounderPy computes the figure's parcels exactly
const RD = 287.04, CP = 1005.7, LV = 2.501e6, EPS = 0.622, KAPPA = RD / CP;
const satMixing = (p, T) => { const e = es(T); return (EPS * e) / (p - e); }; // kg/kg; p hPa, T °C
const virtualK = (Tk, r) => (Tk * (1 + r / EPS)) / (1 + r);
function moistLapse(p, Tk) {
  // dT/dp on a pseudoadiabat (K/hPa)
  const rs = satMixing(p, Tk - 273.15);
  return ((RD * Tk + LV * rs) / (CP + (LV * LV * rs * EPS) / (RD * Tk * Tk))) / p;
}

const fmt = (v, digits = 0, unit = "") => (v == null || !Number.isFinite(v) ? "–" : num(v, digits) + unit);
const windDir = (u, v) => ((Math.atan2(-u, -v) * 180) / Math.PI + 360) % 360;
const COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
const compass = (deg) => COMPASS[Math.round(deg / 22.5) % 16];

export class SounderPyReadout {
  /**
   * Fetch (if given a URL) the readout file and attach to the image.
   * @param {HTMLImageElement} img the SounderPy figure
   * @param {object|string} readout the readout JSON, or its URL
   * @param {object} [options] see the constructor
   */
  static async attach(img, readout, options = {}) {
    if (typeof readout === "string") {
      const res = await fetch(readout);
      if (!res.ok) throw new Error(`Couldn't load the readout file (${res.status}).`);
      readout = await res.json();
    }
    return new SounderPyReadout(img, readout, options);
  }

  /**
   * @param {HTMLImageElement} img the SounderPy figure
   * @param {object} readout the readout JSON
   * @param {object} [options]
   * @param {string[]} [options.fields] readout rows, from: p, z_agl, z, T, Td, Tw, Tv, rh, mr, theta, thetae, wind, cumulative, omega
   * @param {"C"|"F"|"both"} [options.temperature="both"] temperature units
   * @param {"m"|"ft"} [options.height="m"] height units
   * @param {"kt"|"m/s"|"mph"} [options.wind="kt"] wind speed units
   * @param {"auto"|"light"|"dark"} [options.theme="auto"] tooltip colors; auto follows the figure
   * @param {boolean} [options.editable=false] start with editing on (the same as tool: "edit")
   * @param {"inspect"|"edit"|"layer"|"measure"|"storm"|"parcel"|"ptype"|"draw"} [options.tool="inspect"] what dragging
   *   does: nothing (inspect), move profile points (edit), measure a height range (layer), measure between two points
   *   (measure), pick a storm motion on the hodograph (storm), lift a parcel (parcel), follow precipitation down from
   *   a level (ptype), or draw on the plot (draw)
   * @param {boolean} [options.clouds=false] shade cloud layers (RH ≥ 90%, over ice below 0 °C) on the skew-t
   * @param {boolean} [options.keys=true] up and down arrow keys step the readout (or a parcel or precipitation
   *   source) through the profile's levels while the pointer is on the plot
   * @param {HTMLElement} [options.panel] an element beside the image for the full readout and tool results; the
   *   pointer then shows only `tipFields`
   * @param {string[]} [options.tipFields] the pointer readout's rows when there's a panel (p, height, T, Td, wind)
   * @param {boolean} [options.panelTip=true] with a panel, also show the short readout at the pointer; false keeps the
   *   plot clear (e.g. a phone, with the panel under the plot)
   * @param {string} [options.panelHint] what the panel says before anything is hovered
   * @param {string} [options.tipFooter] a line at the bottom of the pointer readout, e.g. "Click for more"
   * @param {number[]|false} [options.levels] pressures (hPa) for the panel's standard levels table, after the surface;
   *   false for none. Default is 925, 850, 700, 500, 300, and 250
   * @param {(drawings: object[]) => void} [options.onDraw] called when drawings are added, undone, or cleared
   * @param {boolean} [options.tapReadout=true] on touch screens, a tap shows the readout; false leaves taps to the
   *   page (e.g. a link around the image) unless editing
   * @param {number} [options.touchOffset=0] on touch screens, read this many screen pixels above a dragging finger,
   *   which would otherwise cover the point; the offset grows from zero as the finger moves, so a plain tap reads
   *   where it lands. Not used while drawing or editing
   * @param {boolean} [options.cursor=true] add the temperature (skew-t) or wind (hodograph) under the pointer
   * @param {string[]} [options.parts] what responds to the pointer: "skewt", "hodograph", "panels" (the 0-3 km panels),
   *   and "params" (values in the tables). Default is all four
   * @param {(level: object|null) => void} [options.onHover] called with the level under the pointer
   * @param {(profile: object) => void} [options.onEdit] called while a point is dragged
   * @param {(profile: object) => void} [options.onEditEnd] called when a drag finishes
   * @param {(mark: object|null) => void} [options.onMark] called with a layer or measurement's results (see getMark)
   * @param {(motion: {u: number, v: number, direction: number, speed: number}) => void} [options.onStormMotion] called
   *   when a storm motion is picked: u and v in kt, and the direction it comes from (degrees) and speed (kt), as
   *   SounderPy's storm_motion=[direction, speed] takes them
   */
  constructor(img, readout, options = {}) {
    if (!(img instanceof HTMLImageElement)) throw new TypeError("SounderPyReadout needs an <img> element.");
    this.img = img;
    this.options = { fields: DEFAULT_FIELDS, temperature: "both", height: "m", wind: "kt", theme: "auto", editable: false, tapReadout: true, touchOffset: 0, cursor: true, clouds: false, keys: true, panelTip: true, tipFields: TIP_FIELDS, panelHint: "Hover the plot to read it here.", levels: LEVELS,
      parts: ["skewt", "hodograph", "panels", "params"], ...options };
    this.tool = "inspect";
    this.mark = null;
    this.drawings = [];
    // fingers on the layer: with two or more, the page is pinching or panning, so the readout waits
    this._touches = new Set();
    this.drawStyle = { mode: "pen", color: "#e11d48", width: 5 };
    this.hover = null;
    this.drag = null;
    this.history = [];
    injectStyle();

    const parent = img.parentElement;
    if (getComputedStyle(parent).position === "static") parent.style.position = "relative";
    this.layer = document.createElement("canvas");
    this.layer.className = "spy-readout-layer";
    this.layer.setAttribute("aria-hidden", "true");
    this.tip = document.createElement("div");
    this.tip.className = "spy-readout-tip";
    this.tip.setAttribute("role", "status");
    // a layer or measurement's results, pinned beside it
    this.card = document.createElement("div");
    this.card.className = "spy-readout-tip spy-readout-card";
    this.card.setAttribute("role", "status");
    img.after(this.layer, this.tip, this.card);
    this.ctx = this.layer.getContext("2d");

    this._listeners = [];
    const on = (target, type, fn, opts) => { target.addEventListener(type, fn, opts); this._listeners.push([target, type, fn, opts]); };
    on(this.layer, "pointermove", (e) => this._onMove(e));
    on(this.layer, "pointerdown", (e) => this._onDown(e));
    on(this.layer, "pointerup", (e) => { this._touches.delete(e.pointerId); this._onUp(e); });
    on(this.layer, "pointercancel", (e) => { this._touches.delete(e.pointerId); this._onUp(e); });
    on(this.layer, "pointerleave", (e) => { if (!this.drag && e.pointerType === "mouse") this._setHover(null); });
    // with a tool on, clicks belong to the tool (not, say, a link or button around the image)
    on(this.layer, "click", (e) => { if (this.tool !== "inspect") { e.preventDefault(); e.stopPropagation(); } });
    // capture phase: the readout's arrow keys don't also pan a zoomable view around the image
    on(document, "keydown", (e) => this._onKey(e), true);
    on(img, "load", () => this._layout());
    this._resize = new ResizeObserver(() => this._layout());
    this._resize.observe(img);
    // an image zoomed or panned with a CSS transform (e.g. in a zoomable viewer) takes the layer with it
    this._restyle = new MutationObserver(() => this._followTransform());
    this._restyle.observe(img, { attributes: true, attributeFilter: ["style", "class"] });

    this.update(readout);
    if (this.options.panel) this.setPanel(this.options.panel);
    if (this.options.editable) this.setTool("edit");
    else if (this.options.tool) this.setTool(this.options.tool);
  }

  // ------------------------------------------------------------------------------------------------
  // public API
  // ------------------------------------------------------------------------------------------------

  /** Swap in a new readout, e.g. after the server redraws the figure. Clears edits. */
  update(readout) {
    if (!readout || readout.format !== FORMAT) throw new Error("That isn't a SounderPy readout file.");
    if (readout.version > VERSION) console.warn(`SounderPy readout: file version ${readout.version} is newer than this script (${VERSION}); update sounderpy-readout.js.`);
    this.ro = readout;
    const P = readout.profile;
    this.n = P.p.length;
    this.lp = P.p.map((p) => Math.log10(p));
    this.original = { T: [...P.T], Td: [...P.Td], u: [...P.u], v: [...P.v] };
    this.edit = structuredClone(this.original);
    this.history = [];
    this.hover = null;
    const dark = this.options.theme === "dark" || (this.options.theme === "auto" && readout.dark);
    this.tip.toggleAttribute("data-dark", dark);
    this.card.toggleAttribute("data-dark", dark);
    this.mark = null;
    // drawings stay on a redrawn figure of the same size (e.g. light to dark); a new layout clears them
    if (this._size && (this._size[0] !== readout.image.width || this._size[1] !== readout.image.height)) this.drawings = [];
    this._size = [readout.image.width, readout.image.height];
    this._layout();
    this._renderLevels();
  }

  /** True while dragging edits the profile. */
  get editing() { return this.tool === "edit"; }

  /** Turn profile editing on or off (the same as setTool("edit") and setTool("inspect")). */
  setEditing(on) {
    this.setTool(on ? "edit" : "inspect");
  }

  /**
   * What dragging does: "inspect" (nothing; hover reads the plot), "edit" (move profile points), "layer" (drag a
   * height range on the skew-t or a 0-3 km panel for its stats), "measure" (drag between two points on the
   * skew-t or hodograph), or "storm" (hover the hodograph to try a storm motion, click to pick it; see
   * onStormMotion). Changing tools clears a layer, measurement, or storm motion.
   */
  setTool(tool) {
    if (!TOOLS.includes(tool)) throw new Error(`Unknown tool "${tool}"; use one of ${TOOLS.join(", ")}.`);
    this.tool = tool;
    this.layer.toggleAttribute("data-editing", tool === "edit");
    this.layer.toggleAttribute("data-marking", [...MARK_TOOLS, "draw"].includes(tool));
    this.layer.dataset.readoutTool = tool;
    this.clearMark();
  }

  /**
   * The current layer or measurement and its results, or null. A layer has z (m AGL) and p (hPa) ranges and its
   * stats; a measurement has its two points and their differences. Each has `rows`: [label, value] pairs.
   */
  getMark() {
    return this.mark ? this._markResults(this.mark) : null;
  }

  /**
   * Put the full readout and tool results in an element beside the image (or null to go back to the pointer only).
   * The element's size and colors are the page's; the readout fills it with sections of rows. `{ tip: false }` drops
   * the short readout at the pointer while the panel is in use.
   */
  setPanel(el, { tip } = {}) {
    if (tip !== undefined) this.options.panelTip = tip;
    this.panel?.classList.remove("spy-readout-panel");
    this.panel = el || null;
    if (this.panel) {
      this.panel.classList.add("spy-readout-panel");
      this.panel.innerHTML = `<section data-readout-mark hidden></section><section data-readout-cursor><p class="spy-readout-empty"></p></section><section data-readout-levels hidden></section>`;
      this.panel.querySelector(".spy-readout-empty").textContent = this.options.panelHint;
      this.card.style.display = "none";
    }
    if (this.ro) { this._draw(); this._renderLevels(); }
  }

  // ------------------------------------------------------------------------------------------------
  // drawing on the plot
  // ------------------------------------------------------------------------------------------------

  /** How the draw tool draws: { mode: "pen" | "arrow" | "text", color: "#rrggbb", width: image pixels }. */
  setDrawStyle(style = {}) {
    Object.assign(this.drawStyle, style);
  }

  /** True when there's something drawn on the plot. */
  hasDrawing() { return this.drawings.length > 0; }

  /** Remove the last drawing. */
  undoDrawing() {
    this.drawings.pop();
    this._draw();
    this.options.onDraw?.(this.drawings);
  }

  /** Remove every drawing. */
  clearDrawing() {
    this.drawings = [];
    this._draw();
    this.options.onDraw?.(this.drawings);
  }

  /**
   * The figure with the drawings on it, as an image Blob at the image's full size. The image must be same-origin
   * (or served with CORS) for the browser to allow this.
   */
  async exportImage(type = "image/png") {
    const img = this.img, w = img.naturalWidth || this.ro.image.width, h = img.naturalHeight || this.ro.image.height;
    const canvas = document.createElement("canvas");
    canvas.width = w; canvas.height = h;
    const ctx = canvas.getContext("2d");
    ctx.drawImage(img, 0, 0, w, h);
    ctx.scale(w / this.ro.image.width, h / this.ro.image.height);
    this._drawDrawings(ctx);
    return new Promise((resolve, reject) => canvas.toBlob((b) => (b ? resolve(b) : reject(new Error("Couldn't export the image."))), type));
  }

  /**
   * Lift and pin a preset parcel like SounderPy's: "sb" (surface-based), "ml" (lowest 100 hPa mixed), or "mu" (most
   * unstable in the lowest 300 hPa). Switches to the parcel tool; returns the parcel's results (see getMark).
   */
  liftParcel(preset) {
    if (!["sb", "ml", "mu"].includes(preset)) throw new Error(`Unknown parcel "${preset}"; use "sb", "ml", or "mu".`);
    if (this.tool !== "parcel") this.setTool("parcel");
    const S = this.ro.panels.skewt;
    this.mark = { kind: "parcel", preset, i: 0, picked: true, end: [S.bbox[2] - 20, S.bbox[1] + 20] };
    this._draw();
    const res = this.getMark();
    this.options.onMark?.(res);
    return res;
  }

  /** Remove a layer or measurement. */
  clearMark() {
    this.mark = null;
    this.card.style.display = "none";
    const section = this.panel?.querySelector("[data-readout-mark]");
    if (section) section.hidden = true;
    this._draw();
    this.options.onMark?.(null);
  }

  /**
   * Change the readout's units: { temperature: "C" | "F" | "both", height: "m" | "ft", wind: "kt" | "m/s" | "mph" }.
   * The figure keeps its own units; only the readouts change.
   */
  setUnits(units = {}) {
    for (const key of ["temperature", "height", "wind"]) if (units[key]) this.options[key] = units[key];
    if (this.ro) this._draw();
    this._renderLevels();
  }

  /** Shade cloud layers (RH ≥ 90%, over ice below 0 °C) on the skew-t, and list them in the panel. */
  setClouds(on) {
    this.options.clouds = Boolean(on);
    if (this.ro) this._draw();
    this._renderLevels();
  }

  /** The cloud layers: [{bottom, top} m AGL, {bottom_p, top_p} hPa], lowest first. */
  getClouds() {
    return this._clouds().map(([a, b]) => ({ bottom: a, top: b, bottom_p: 10 ** this._lpAtHeight(a), top_p: 10 ** this._lpAtHeight(b) }));
  }

  /** True when the profile has been edited. */
  isEdited() {
    const { T, Td, u, v } = this.original;
    const differs = (a, b) => a.some((x, i) => (x == null) !== (b[i] == null) || Math.abs(x - b[i]) > 0.05);
    return differs(T, this.edit.T) || differs(Td, this.edit.Td) || differs(u, this.edit.u) || differs(v, this.edit.v);
  }

  /**
   * The (edited) profile, one value per plotted level: p (hPa), z (m MSL), T and Td (°C), u and v (kt).
   * Send it to a server and pass it to sounderpy.apply_profile_edits(clean_data, profile) to redraw.
   */
  getProfile() {
    const r = (a) => a.map((x) => (x == null ? null : Math.round(x * 100) / 100));
    return { p: [...this.ro.profile.p], z: [...this.ro.profile.z], T: r(this.edit.T), Td: r(this.edit.Td), u: r(this.edit.u), v: r(this.edit.v) };
  }

  /** Undo the last drag. */
  undo() {
    if (!this.history.length) return;
    this.edit = this.history.pop();
    this._draw();
    this._renderLevels();
    this.options.onEditEnd?.(this.getProfile());
  }

  /** Put the profile back to the figure's. */
  resetEdits() {
    this.edit = structuredClone(this.original);
    this.history = [];
    this._draw();
    this._renderLevels();
    this.options.onEditEnd?.(this.getProfile());
  }

  /** The readout values at a pressure (hPa), or null outside the profile. */
  levelAt(p) {
    return this._level(Math.log10(p));
  }

  /** Show or hide the layer, e.g. for a readout mode that can be switched off. */
  setEnabled(on) {
    this.layer.hidden = !on;
    if (!on) { this._setHover(null); this.card.style.display = "none"; }
    else if (this.mark) this._draw();
  }

  /** Remove the layer and listeners. */
  destroy() {
    this._resize.disconnect();
    this._restyle.disconnect();
    cancelAnimationFrame(this._frame);
    for (const [target, type, fn, opts] of this._listeners) target.removeEventListener(type, fn, opts);
    this.layer.remove();
    this.tip.remove();
    this.card.remove();
  }

  // ------------------------------------------------------------------------------------------------
  // geometry: readout panels map data to image pixels with a 2x3 affine matrix
  // ------------------------------------------------------------------------------------------------
  // image pixels -> the image's layout pixels, and -> screen pixels (including any CSS zoom)
  get _scale() { return this.img.offsetWidth / this.ro.image.width; }
  get _screenScale() { return this.img.getBoundingClientRect().width / this.ro.image.width || this._scale; }
  _wants(part) { return this.options.parts.includes(part); }
  _fwd(M, a, b) { return [M[0][0] * a + M[0][1] * b + M[0][2], M[1][0] * a + M[1][1] * b + M[1][2]]; }
  _inv(M, x, y) {
    const [a, b, c] = M[0], [d, e, f] = M[1], det = a * e - b * d;
    x -= c; y -= f;
    return [(e * x - b * y) / det, (a * y - d * x) / det];
  }
  _in(bbox, x, y) { return x >= bbox[0] && x <= bbox[2] && y >= bbox[1] && y <= bbox[3]; }
  _skewPx(T, i) { return this._fwd(this.ro.panels.skewt.transform, T, this.lp[i]); }
  // the hodograph is drawn shifted by the storm motion when storm relative
  _hodoPx(u, v) {
    const H = this.ro.panels.hodograph;
    return this._fwd(H.transform, u - H.offset[0], v - H.offset[1]);
  }
  _hodoUV(x, y) {
    const H = this.ro.panels.hodograph, [u, v] = this._inv(H.transform, x, y);
    return [u + H.offset[0], v + H.offset[1]];
  }
  _eventPx(e) {
    // works when the image (or a parent) is scaled with CSS, since the rect includes transforms
    const rect = this.img.getBoundingClientRect();
    let y = e.clientY;
    // a finger covers what it points at: once it drags, read a little above it (options.touchOffset)
    const lift = this.options.touchOffset;
    if (lift && e.pointerType === "touch" && this._touchStart && this.tool !== "draw" && this.tool !== "edit") {
      y -= Math.min(lift, Math.hypot(e.clientX - this._touchStart[0], e.clientY - this._touchStart[1]));
    }
    return [((e.clientX - rect.left) / rect.width) * this.ro.image.width, ((y - rect.top) / rect.height) * this.ro.image.height];
  }
  _hodoLevels() {
    const H = this.ro.panels.hodograph, z = this.ro.profile.z_agl, out = [];
    for (let i = 0; i < this.n; i++) if (z[i] != null && z[i] <= H.max_height && this.edit.u[i] != null) out.push(i);
    return out;
  }

  _layout() {
    const img = this.img, w = img.offsetWidth, h = img.offsetHeight;
    Object.assign(this.layer.style, { left: img.offsetLeft + "px", top: img.offsetTop + "px", width: w + "px", height: h + "px" });
    this._copyTransform();
    // drawn at the zoom it's shown at (and the screen's pixel density), up to about 16 megapixels
    const zoom = w ? this.img.getBoundingClientRect().width / w : 1;
    const res = Math.min((window.devicePixelRatio || 1) * Math.max(zoom, 1), Math.sqrt(16e6 / Math.max(w * h, 1)));
    this._res = res;
    this.layer.width = Math.max(1, Math.round(w * res));
    this.layer.height = Math.max(1, Math.round(h * res));
    this._draw();
  }

  _copyTransform() {
    const style = getComputedStyle(this.img);
    this.layer.style.transform = style.transform === "none" ? "" : style.transform;
    this.layer.style.transformOrigin = style.transformOrigin;
  }

  _followTransform() {
    this._copyTransform();
    if (this.mark) this._placeCard();
    // sharpen once the zoom settles
    cancelAnimationFrame(this._frame);
    this._frame = requestAnimationFrame(() => {
      const w = this.img.offsetWidth, zoom = w ? this.img.getBoundingClientRect().width / w : 1;
      const want = Math.min((window.devicePixelRatio || 1) * Math.max(zoom, 1), Math.sqrt(16e6 / Math.max(w * this.img.offsetHeight, 1)));
      if (Math.abs(want - (this._res || 0)) / want > 0.15) this._layout();
    });
  }

  // ------------------------------------------------------------------------------------------------
  // values
  // ------------------------------------------------------------------------------------------------
  _interp(arr, lp) {
    const L = this.lp;
    for (let i = 0; i < this.n - 1; i++) {
      if ((L[i] - lp) * (L[i + 1] - lp) <= 0) {
        if (arr[i] == null || arr[i + 1] == null) return arr[i] ?? arr[i + 1] ?? null;
        const t = L[i + 1] === L[i] ? 0 : (lp - L[i]) / (L[i + 1] - L[i]);
        return arr[i] + t * (arr[i + 1] - arr[i]);
      }
    }
    return null;
  }

  _zAtLp(lp) { return this._interp(this.ro.profile.z_agl, lp); }
  _panelY(panel, z) { return this._fwd(panel.transform, 0, z)[1]; }
  _panelZ(panel, x, y) { return this._inv(panel.transform, x, y)[1]; }

  // a pointer just past the lowest or highest level (within a few pixels) reads that level
  _snapToProfile(lp) {
    const lo = Math.min(this.lp[0], this.lp[this.n - 1]), hi = Math.max(this.lp[0], this.lp[this.n - 1]);
    const slack = 6 / Math.abs(this.ro.panels.skewt.transform[1][1]);
    if (lp > hi && lp - hi <= slack) return hi;
    if (lp < lo && lo - lp <= slack) return lo;
    return lp;
  }

  _lpAtHeight(z) {
    const Z = this.ro.profile.z_agl;
    for (let i = 0; i < this.n - 1; i++) {
      if (Z[i] <= z && z <= Z[i + 1]) return this.lp[i] + ((z - Z[i]) / (Z[i + 1] - Z[i] || 1)) * (this.lp[i + 1] - this.lp[i]);
    }
    return null;
  }

  _level(lp) {
    if (lp == null) return null;
    const P = this.ro.profile, E = this.edit, f = (a) => (a ? this._interp(a, lp) : null);
    const p = 10 ** lp, T = f(E.T);
    if (T == null) return null;
    const Td = f(E.Td), u = f(E.u), v = f(E.v);
    const edited = this.isEdited();
    const lvl = { p, z: f(P.z), z_agl: f(P.z_agl), T, Td, u, v, omega: f(P.omega), edited };
    lvl.wspd = u == null ? null : Math.hypot(u, v);
    lvl.wdir = u == null ? null : windDir(u, v);
    if (!edited) {
      for (const key of ["Tw", "Tv", "rh", "mr", "theta", "thetae"]) lvl[key] = f(P[key]);
    } else if (Td != null) {
      // computed here from the edited T and Td
      lvl.rh = (100 * es(Td)) / es(T);
      lvl.mr = mixingRatio(p, Td);
      lvl.Tw = wetBulb(p, T, Td);
      lvl.Tv = (T + 273.15) * (1 + 0.61 * lvl.mr / 1000) - 273.15;
      lvl.theta = (T + 273.15) * Math.pow(1000 / p, 0.2854);
      lvl.thetae = thetaE(p, T, Td);
    }
    return lvl;
  }

  // ------------------------------------------------------------------------------------------------
  // pointer
  // ------------------------------------------------------------------------------------------------
  _onMove(e) {
    if (this._touches.size > 1) return;
    this._client = [e.clientX, e.clientY];
    this._touch = e.pointerType !== "mouse";
    const [x, y] = this._eventPx(e);
    if (this.drag) { this._dragTo(x, y); return; }
    if (this.marking) { this._markTo(x, y); return; }
    if (this.drawing) { this._drawTo(x, y); return; }
    // the parcel tool previews a parcel from the level under the pointer, until one is picked
    // with a parcel pinned, the pointer's height on the skew-t reads its CAPE from the LFC up to there
    if (this.tool === "parcel" && this.mark?.kind === "parcel" && this.mark.picked) {
      const S = this.ro.panels.skewt, inside = S && this._in(S.bbox, x, y);
      this.mark.cursor = inside ? this._inv(S.transform, x, y)[1] : null;
      this._draw();
      return;
    }
    if ((this.tool === "parcel" || this.tool === "ptype") && !this.mark?.picked) {
      const S = this.ro.panels.skewt;
      if (S && this._in(S.bbox, x, y)) { this.mark = this._markStart(x, y); this.hover = null; this._draw(); return; }
      if (this.mark) { this.mark = null; this._draw(); }
    }
    if (this.tool === "draw") return;
    // the storm tool previews a storm motion wherever the pointer is on the hodograph, until one is picked
    if (this.tool === "storm" && !this.mark?.picked) {
      const H = this.ro.panels.hodograph;
      if (H && this._in(H.bbox, x, y)) { this.mark = this._markStart(x, y); this.hover = null; this._draw(); return; }
      if (this.mark) { this.mark = null; this._draw(); }
    }
    this._setHover(this._hoverAt(x, y));
    this.layer.style.cursor = this.editing && this._hit(x, y) ? "grab" : "crosshair";
  }

  _hoverAt(x, y) {
    const S = this.ro.panels.skewt, H = this.ro.panels.hodograph;
    // a value in the figure's tables
    const pad = 3 / Math.max(this._screenScale, 0.05);
    for (const param of this._wants("params") ? this.ro.params || [] : []) {
      const b = param.bbox;
      if (x >= b[0] - pad && x <= b[2] + pad && y >= b[1] - pad && y <= b[3] + pad) return { panel: "param", param, x, y };
    }
    // a 0-3 km height panel
    for (const hp of this._wants("panels") ? this.ro.height_panels || [] : []) {
      if (this._in(hp.bbox, x, y)) {
        const z = Math.max(0, this._panelZ(hp, x, y));
        return { panel: "height", hp, z, lp: this._lpAtHeight(z), x, y };
      }
    }
    if (S && this._wants("skewt") && this._in(S.bbox, x, y)) {
      const [Tc, lp] = this._inv(S.transform, x, y);
      return { panel: "skewt", lp: this._snapToProfile(lp), x, y, cursorT: Tc };
    }
    if (H && this._wants("hodograph") && this._in(H.bbox, x, y)) {
      // the nearest point on the drawn trace, every 100 m from the readout file
      const tr = this.ro.hodograph_trace;
      let best = null, bd = Infinity;
      if (tr && !this.isEdited()) {
        for (let i = 0; i < tr.z_agl.length; i++) {
          const [hx, hy] = this._hodoPx(tr.u[i], tr.v[i]), d = Math.hypot(hx - x, hy - y);
          if (d < bd) { bd = d; best = this._lpAtHeight(tr.z_agl[i]); }
        }
      } else {
        for (const i of this._hodoLevels()) {
          const [hx, hy] = this._hodoPx(this.edit.u[i], this.edit.v[i]), d = Math.hypot(hx - x, hy - y);
          if (d < bd) { bd = d; best = this.lp[i]; }
        }
      }
      return best == null ? null : { panel: "hodograph", lp: best, x, y, cursorUV: this._hodoUV(x, y) };
    }
    return null;
  }

  _hit(x, y) {
    const tol = 14 / Math.max(this._screenScale, 0.05);
    const S = this.ro.panels.skewt, H = this.ro.panels.hodograph;
    let best = null, bd = tol;
    if (S && this._in(S.bbox, x, y)) {
      for (let i = 0; i < this.n; i++) {
        for (const field of ["T", "Td"]) {
          if (this.edit[field][i] == null) continue;
          const [px, py] = this._skewPx(this.edit[field][i], i);
          if (!this._in(S.bbox, px, py)) continue;
          const d = Math.hypot(px - x, py - y);
          if (d < bd) { bd = d; best = { field, i }; }
        }
      }
    } else if (H && this._in(H.bbox, x, y)) {
      for (const i of this._hodoLevels()) {
        const [px, py] = this._hodoPx(this.edit.u[i], this.edit.v[i]), d = Math.hypot(px - x, py - y);
        if (d < bd) { bd = d; best = { field: "wind", i }; }
      }
    }
    return best;
  }

  _onDown(e) {
    if (e.pointerType === "touch") {
      this._touches.add(e.pointerId);
      if (this._touches.size === 1) this._touchStart = [e.clientX, e.clientY]; // where touchOffset measures the drag from
    }
    if (this._touches.size > 1) { this.tip.style.display = "none"; return; }
    this._client = [e.clientX, e.clientY];
    this._touch = e.pointerType !== "mouse";
    const [x, y] = this._eventPx(e);
    if (this.editing) {
      const hit = this._hit(x, y);
      if (hit) {
        // the point moves, not the page or a zoomable view around the image
        e.preventDefault();
        e.stopPropagation();
        this.history.push(structuredClone(this.edit));
        this.drag = hit;
        this.layer.setPointerCapture(e.pointerId);
        this.layer.style.cursor = "grabbing";
        return;
      }
    }
    if (this.tool === "draw") {
      // drawing, not panning a zoomable view around the image
      e.preventDefault();
      e.stopPropagation();
      this._drawStart(x, y, e);
      return;
    }
    if (MARK_TOOLS.includes(this.tool)) {
      const start = this._markStart(x, y);
      if (start) {
        // the drag measures, not pans a zoomable view around the image
        e.preventDefault();
        e.stopPropagation();
        this.mark = start;
        this.hover = null;
        this.marking = { x0: e.clientX, y0: e.clientY };
        this.layer.setPointerCapture(e.pointerId);
        this._draw();
        return;
      }
    }
    // a tap shows the readout on touch screens
    if (e.pointerType !== "mouse" && (this.options.tapReadout || this.tool !== "inspect")) this._setHover(this._hoverAt(x, y));
  }

  _onUp(e) {
    if (this.drawing) { this._drawEnd(e); return; }
    if (this.marking) {
      const moved = Math.hypot(e.clientX - this.marking.x0, e.clientY - this.marking.y0);
      this.marking = null;
      if (this.layer.hasPointerCapture?.(e.pointerId)) this.layer.releasePointerCapture(e.pointerId);
      if (this.mark?.kind === "parcel" || this.mark?.kind === "ptype") {
        // a click (or a drag that ends) on the skew-t pins the parcel or precipitation source
        this.mark.picked = true;
        this._draw();
        this.options.onMark?.(this._markResults(this.mark));
        return;
      }
      if (this.mark?.kind === "storm") {
        // a click (or a drag that ends) on the hodograph picks the storm motion
        this.mark.picked = true;
        this._draw();
        const res = this._markResults(this.mark);
        this.options.onMark?.(res);
        this.options.onStormMotion?.({ u: this.mark.u, v: this.mark.v, direction: Math.round(res.direction) % 360, speed: Math.round(res.speed) });
        return;
      }
      // a click without a drag clears the mark
      if (moved < 4) this.clearMark();
      else { this._draw(); this.options.onMark?.(this.getMark()); }
      return;
    }
    if (!this.drag) return;
    this.drag = null;
    if (this.layer.hasPointerCapture?.(e.pointerId)) this.layer.releasePointerCapture(e.pointerId);
    this.layer.style.cursor = "crosshair";
    if (this.history.length && JSON.stringify(this.history.at(-1)) === JSON.stringify(this.edit)) this.history.pop();
    this.options.onEditEnd?.(this.getProfile());
  }

  _dragTo(x, y) {
    const { field, i } = this.drag;
    if (field === "wind") {
      const [u, v] = this._hodoUV(x, y);
      this.edit.u[i] = u;
      this.edit.v[i] = v;
    } else {
      // along the level's isobar: solve x = a*T + b*log10(p) + c for T
      const [a, b, c] = this.ro.panels.skewt.transform[0];
      const T = (x - b * this.lp[i] - c) / a;
      if (field === "T") {
        this.edit.T[i] = T;
        if (this.edit.Td[i] != null) this.edit.Td[i] = Math.min(this.edit.Td[i], T);
      } else {
        this.edit.Td[i] = Math.min(T, this.edit.T[i]);
      }
    }
    this.hover = { panel: field === "wind" ? "hodograph" : "skewt", lp: this.lp[i], x, y };
    this._draw();
    this._renderLevels();
    this.options.onEdit?.(this.getProfile());
  }

  // ------------------------------------------------------------------------------------------------
  // layer and measure tools
  // ------------------------------------------------------------------------------------------------
  // where a layer or measurement starts, or null when the pointer isn't on a panel the tool works on
  _markStart(x, y) {
    const S = this.ro.panels.skewt, H = this.ro.panels.hodograph;
    if (this.tool === "layer") {
      for (const hp of this.ro.height_panels || []) {
        if (this._in(hp.bbox, x, y)) { const z = this._clampZ(this._panelZ(hp, x, y)); return { kind: "layer", panel: hp, z0: z, z1: z, end: [x, y] }; }
      }
      if (S && this._in(S.bbox, x, y)) {
        const z = this._zFromSkew(x, y), snap = this._snapP(y);
        return z == null ? null : { kind: "layer", panel: "skewt", z0: z, z1: z, snap0: snap, snap1: snap, end: [x, y] };
      }
    } else if (this.tool === "measure") {
      if (S && this._in(S.bbox, x, y)) return { kind: "measure", panel: "skewt", a: [x, y], b: [x, y], end: [x, y] };
      if (H && this._in(H.bbox, x, y)) return { kind: "measure", panel: "hodograph", a: [x, y], b: [x, y], end: [x, y] };
    } else if (this.tool === "storm") {
      if (H && this._in(H.bbox, x, y)) { const [u, v] = this._snapMotion(...this._hodoUV(x, y)); return { kind: "storm", u, v, end: [x, y] }; }
    } else if (this.tool === "ptype") {
      if (S && this._in(S.bbox, x, y)) return { kind: "ptype", lp: this._clampLp(this._inv(S.transform, x, y)[1]), end: [x, y] };
    } else if (this.tool === "parcel") {
      if (S && this._in(S.bbox, x, y)) {
        const i = this._levelNear(this._inv(S.transform, x, y)[1]);
        return i == null ? null : { kind: "parcel", i, end: [x, y] };
      }
    }
    return null;
  }

  _markTo(x, y) {
    const m = this.mark;
    if (m.kind === "layer") {
      const z = m.panel === "skewt" ? this._zFromSkew(x, y) : this._clampZ(this._panelZ(m.panel, x, y));
      if (z != null) m.z1 = z;
      if (m.panel === "skewt") m.snap1 = this._snapP(y);
    } else if (m.kind === "parcel") {
      const i = this._levelNear(this._inv(this.ro.panels.skewt.transform, x, y)[1]);
      if (i != null) m.i = i;
    } else if (m.kind === "ptype") {
      m.lp = this._clampLp(this._inv(this.ro.panels.skewt.transform, x, y)[1]);
    } else if (m.kind === "storm") {
      const H = this.ro.panels.hodograph;
      [m.u, m.v] = this._snapMotion(...this._hodoUV(Math.min(Math.max(x, H.bbox[0]), H.bbox[2]), Math.min(Math.max(y, H.bbox[1]), H.bbox[3])));
    } else {
      const bbox = this.ro.panels[m.panel].bbox;
      m.b = [Math.min(Math.max(x, bbox[0]), bbox[2]), Math.min(Math.max(y, bbox[1]), bbox[3])];
    }
    m.end = [x, y];
    this._draw();
  }

  // ---- the draw tool ----
  _drawStart(x, y, e) {
    const { mode, color, width } = this.drawStyle;
    this.layer.setPointerCapture(e.pointerId);
    this.drawing = { start: [e.clientX, e.clientY], item: mode === "pen" ? { type: "pen", color, width, pts: [[x, y]] }
      : mode === "arrow" ? { type: "arrow", color, width, a: [x, y], b: [x, y] } : null, text: mode === "text" ? [x, y] : null };
  }

  _drawTo(x, y) {
    const item = this.drawing.item;
    if (item?.type === "pen") item.pts.push([x, y]);
    else if (item?.type === "arrow") item.b = [x, y];
    if (item && !this.drawings.includes(item)) this.drawings.push(item);
    this._draw();
  }

  _drawEnd(e) {
    const { item, text, start } = this.drawing;
    this.drawing = null;
    if (this.layer.hasPointerCapture?.(e.pointerId)) this.layer.releasePointerCapture(e.pointerId);
    if (text) { this._textInput(text, [e.clientX, e.clientY]); return; }
    const moved = Math.hypot(e.clientX - start[0], e.clientY - start[1]);
    // a tap without a drag draws nothing
    if (moved < 3 && item) this.drawings = this.drawings.filter((d) => d !== item);
    this._draw();
    this.options.onDraw?.(this.drawings);
  }

  // a small box to type a label where the plot was clicked: it wraps as you type, Enter adds the label (Shift+Enter
  // starts a new line), and Escape cancels. It stays on screen and doesn't scroll the view under it.
  _textInput(at, client) {
    const input = document.createElement("textarea");
    input.className = "spy-readout-text-input";
    input.rows = 1;
    input.placeholder = "Label, then Enter";
    const parent = this.layer.offsetParent || this.img.parentElement, box = parent.getBoundingClientRect();
    // inside the part of the image that's on screen, so the label stays on the figure
    const img = this.img.getBoundingClientRect();
    const area = { left: Math.max(box.left, img.left), right: Math.min(box.right, img.right), top: Math.max(box.top, img.top), bottom: Math.min(box.bottom, img.bottom) };
    const width = Math.max(120, Math.min(260, area.right - area.left - 16));
    const left = Math.min(Math.max(client[0], area.left + 8), area.right - width - 8) - box.left;
    const top = Math.min(Math.max(client[1] - 14, area.top + 8), area.bottom - 60) - box.top;
    // the label's size on the plot, so what's typed looks like what's drawn
    // at least 12 px on screen while typing; drawn at that size too, so the label matches the box
    const screenSize = Math.max(12, Math.max(24, this.drawStyle.width * 6) * this._screenScale);
    const size = screenSize / Math.max(this._screenScale, 0.05);
    Object.assign(input.style, { left: `${left}px`, top: `${top}px`, width: `${width}px`, borderColor: this.drawStyle.color,
      fontSize: `${screenSize}px`, lineHeight: "1.2", color: this.drawStyle.color });
    parent.append(input);
    input.focus({ preventScroll: true });
    const grow = () => { input.style.height = "auto"; input.style.height = `${input.scrollHeight + 2}px`; };
    input.addEventListener("input", grow);
    let done = false;
    const finish = (keep) => {
      if (done) return;
      done = true;
      const text = input.value.trim();
      // drawn where the box is (it may have moved to stay on screen), wrapped at the box's width, in image pixels
      const scale = Math.max(this._screenScale, 0.05), box = input.getBoundingClientRect(), img = this.img.getBoundingClientRect();
      const maxWidth = (input.clientWidth - 12) / scale;
      const x = ((box.left + 8 - img.left) / img.width) * this.ro.image.width;
      const y = ((box.top + 6 - img.top) / img.height) * this.ro.image.height + size * 0.9;
      input.remove();
      if (keep && text) {
        this.drawings.push({ type: "text", color: this.drawStyle.color, size, at: [x, y], text, lines: this._wrapText(text, size, maxWidth) });
        this._draw();
        this.options.onDraw?.(this.drawings);
      }
    };
    input.addEventListener("keydown", (ev) => {
      ev.stopPropagation(); // typing doesn't zoom or pan a viewer
      if (ev.key === "Enter" && !ev.shiftKey) { ev.preventDefault(); finish(true); }
      if (ev.key === "Escape") { ev.preventDefault(); finish(false); }
    });
    input.addEventListener("blur", () => finish(true));
  }

  // lines of text that fit maxWidth (image pixels) at a font size, keeping the typed line breaks
  _wrapText(text, size, maxWidth) {
    this.ctx.save();
    this.ctx.font = `bold ${size}px system-ui, sans-serif`;
    const lines = [];
    for (const paragraph of text.split("\n")) {
      let line = "";
      for (const word of paragraph.split(/\s+/)) {
        const next = line ? `${line} ${word}` : word;
        if (line && this.ctx.measureText(next).width > maxWidth) { lines.push(line); line = word; } else line = next;
      }
      lines.push(line);
    }
    this.ctx.restore();
    return lines;
  }

  _drawDrawings(ctx) {
    for (const d of this.drawings) {
      ctx.save();
      ctx.strokeStyle = d.color; ctx.fillStyle = d.color; ctx.lineWidth = d.width ?? 5; ctx.lineCap = "round"; ctx.lineJoin = "round";
      if (d.type === "pen" && d.pts.length > 1) {
        ctx.beginPath(); ctx.moveTo(...d.pts[0]); d.pts.slice(1).forEach((pt) => ctx.lineTo(...pt)); ctx.stroke();
      } else if (d.type === "arrow") {
        const [ax, ay] = d.a, [bx, by] = d.b, ang = Math.atan2(by - ay, bx - ax), head = 10 + d.width * 3;
        ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx - Math.cos(ang) * head * 0.6, by - Math.sin(ang) * head * 0.6); ctx.stroke();
        ctx.beginPath(); ctx.moveTo(bx, by);
        ctx.lineTo(bx - head * Math.cos(ang - 0.45), by - head * Math.sin(ang - 0.45));
        ctx.lineTo(bx - head * Math.cos(ang + 0.45), by - head * Math.sin(ang + 0.45));
        ctx.closePath(); ctx.fill();
      } else if (d.type === "text") {
        // a light halo keeps labels readable on any part of the figure
        ctx.font = `bold ${d.size}px system-ui, sans-serif`;
        ctx.lineWidth = d.size / 5; ctx.strokeStyle = "rgba(255,255,255,0.9)";
        (d.lines ?? [d.text]).forEach((line, k) => {
          const y = d.at[1] + k * d.size * 1.2;
          ctx.strokeText(line, d.at[0], y); ctx.fillText(line, d.at[0], y);
        });
      }
      ctx.restore();
    }
  }

  // ---- the panel's standard levels table: the surface and LEVELS, in the readout's units ----
  _renderLevels() {
    const section = this.panel?.querySelector("[data-readout-levels]");
    if (!section) return;
    const levels = this.options.levels;
    if (!levels || !this.ro) { section.hidden = true; return; }
    const psfc = this.ro.profile.p[0], ptop = this.ro.profile.p[this.n - 1];
    const tempUnit = this.options.temperature === "F" ? "°F" : "°C", ft = this.options.height === "ft";
    const t = (c) => (c == null ? "–" : num(this.options.temperature === "F" ? (c * 9) / 5 + 32 : c, this.options.temperature === "F" ? 0 : 1));
    const spd = (kt) => (this.options.wind === "m/s" ? kt * KT : this.options.wind === "mph" ? kt * 1.150779 : kt);
    const rows = [["SFC", this.lp[0]], ...levels.filter((p) => p < psfc && p >= ptop).map((p) => [String(p), Math.log10(p)])].map(([label, lp]) => {
      const l = this._level(lp);
      if (!l) return "";
      const z = l.z_agl == null ? "–" : num(ft ? l.z_agl * 3.28084 : l.z_agl);
      const wind = l.u == null ? "–" : Math.hypot(l.u, l.v) < 0.5 ? "calm" : `${num(windDir(l.u, l.v) / 5) * 5}/${num(spd(Math.hypot(l.u, l.v)))}`;
      return `<tr><td>${label}</td><td>${z}</td><td>${t(l.T)}</td><td>${t(l.Td)}</td><td>${wind}</td></tr>`;
    }).join("");
    // thicknesses (dam) between standard levels in the profile
    const zMsl = (p) => {
      if (p <= psfc) return this._interp(this.ro.profile.z, Math.log10(p));
      // below ground: the hypsometric equation with the surface's virtual temperature
      const sfc = this._level(this.lp[0]), Tv = (sfc?.Tv ?? sfc?.T) + 273.15;
      return this.ro.profile.z[0] - ((RD * Tv) / 9.81) * Math.log(p / psfc);
    };
    const thick = THICKNESSES.filter(([a, b]) => a <= psfc + 50 && b >= ptop).map(([a, b]) => {
      const za = zMsl(a), zb = zMsl(b);
      return za == null || zb == null ? "" : `<tr><td>${a}–${b}</td><td>${num((zb - za) / 10)}</td></tr>`;
    }).join("");
    let clouds = "";
    if (this.options.clouds) {
      const list = this.getClouds().map((c) => `<tr><td>${num(c.bottom_p)}–${num(c.top_p)}</td><td>${this._heightRange(c.bottom, c.top)}</td></tr>`).join("");
      clouds = `<p class="spy-readout-title spy-readout-sub">Clouds</p>${list ? `<table><thead><tr><th>hPa</th><th>AGL</th></tr></thead><tbody>${list}</tbody></table>`
        : '<p class="spy-readout-empty">None (RH under 90%)</p>'}`;
    }
    section.innerHTML = `<p class="spy-readout-title">Levels</p><table><thead><tr><th>hPa</th><th>${ft ? "ft" : "m"} AGL</th>`
      + `<th>T ${tempUnit}</th><th>Td</th><th>Wind ${this.options.wind}</th></tr></thead><tbody>${rows}</tbody></table>`
      + (thick ? `<p class="spy-readout-title spy-readout-sub">Thickness</p><table><thead><tr><th>hPa</th><th>dam</th></tr></thead><tbody>${thick}</tbody></table>` : "")
      + clouds;
    section.hidden = false;
  }

  // the profile level nearest log10(p), among levels with a temperature and dewpoint
  _levelNear(lp) {
    let best = null, bd = Infinity;
    for (let i = 0; i < this.n; i++) {
      if (this.edit.T[i] == null || this.edit.Td[i] == null) continue;
      const d = Math.abs(this.lp[i] - lp);
      if (d < bd) { bd = d; best = i; }
    }
    return best;
  }

  // the environment's virtual temperature (K) at log10(p), from the (edited) temperature and dewpoint
  _envTv(lp) {
    const T = this._interp(this.edit.T, lp), Td = this._interp(this.edit.Td, lp);
    if (T == null) return null;
    const p = 10 ** lp;
    return virtualK(T + 273.15, Td == null ? 0 : satMixing(p, Td));
  }

  // a parcel mark's parcel, recomputed only when its source or the profile changes
  _parcelFor(m) {
    const key = `${m.preset ?? m.i}:${JSON.stringify(this.edit.T)}:${JSON.stringify(this.edit.Td)}`;
    if (this._parcelCache?.key !== key) this._parcelCache = { key, r: this._parcel(this._parcelSource(m)) };
    return this._parcelCache.r;
  }

  // where a parcel starts: a level, or a preset like SounderPy's parcels
  //   sb: the surface; ml: the lowest 100 hPa mixed (mean potential temperature and mixing ratio, by pressure),
  //   from the surface; mu: the highest theta-e in the lowest 300 hPa
  _parcelSource(m) {
    const P = this.ro.profile.p, ok = (k) => this.edit.T[k] != null && this.edit.Td[k] != null;
    const sfc = P.findIndex((_, k) => ok(k));
    const at = (k, label) => ({ p0: P[k], T0: this.edit.T[k], Td0: Math.min(this.edit.Td[k], this.edit.T[k]), label });
    if (m.preset === "sb") return at(sfc, "Surface-based");
    if (m.preset === "mu") {
      let best = sfc, top = -Infinity;
      for (let k = sfc; k < this.n && P[k] >= P[sfc] - 300; k++) {
        if (!ok(k)) continue;
        const te = thetaE(P[k], this.edit.T[k], this.edit.Td[k]);
        if (te > top) { top = te; best = k; }
      }
      return at(best, "Most unstable");
    }
    if (m.preset === "ml") {
      // pressure-weighted means every 1 hPa through the lowest 100 hPa
      const p0 = P[sfc];
      let theta = 0, r = 0, weight = 0;
      for (let p = p0; p >= p0 - 100; p -= 1) {
        const lp = Math.log10(p), T = this._interp(this.edit.T, lp), Td = this._interp(this.edit.Td, lp);
        if (T == null || Td == null) continue;
        theta += (T + 273.15) * Math.pow(1000 / p, KAPPA); r += satMixing(p, Td); weight += 1;
      }
      theta /= weight; r /= weight;
      const T0 = theta * Math.pow(p0 / 1000, KAPPA) - 273.15;
      const e = (r * p0) / (EPS + r), x = Math.log(e / 6.112);
      return { p0, T0, Td0: Math.min((243.5 * x) / (17.67 - x), T0), label: "Mixed layer (100 hPa)" };
    }
    return at(m.i, null);
  }

  // CAPE and CIN of a parcel from the start up to log10(p) = lpTop (CAPE counts from the LFC), and its buoyancy there
  _parcelTo(r, lpTop) {
    const pts = r.trace, lo = pts[pts.length - 1].lp;
    const top = Math.max(Math.min(lpTop, pts[0].lp), lo);
    let cape = 0, cin = 0, b = pts[0].b;
    for (let k = 0; k < pts.length - 1; k++) {
      const a = pts[k];
      if (a.lp <= top) break;
      let c = pts[k + 1];
      // the last step ends at the cursor
      if (c.lp < top) { const f = (a.lp - top) / (a.lp - c.lp); c = { lp: top, p: 10 ** top, b: a.b + f * (c.b - a.b) }; }
      const area = RD * ((a.b + c.b) / 2) * Math.log(a.p / c.p);
      if (r.lfc != null && c.lp <= r.lfc && (r.el == null || a.lp >= r.el) && area > 0) cape += area;
      if (r.lfc != null && a.lp > r.lfc && area < 0) cin += area;
      b = c.b;
    }
    return { lp: top, cape, cin, b };
  }

  // a parcel lifted from level i: its trace and CAPE, CIN, LCL, LFC, EL, and lifted index
  _parcel({ p0, T0, Td0 }) {
    const tk0 = T0 + 273.15, tdk0 = Td0 + 273.15, r0 = satMixing(p0, Td0);
    const tlcl = 1 / (1 / (tdk0 - 56) + Math.log(tk0 / tdk0) / 800) + 56;
    const plcl = p0 * Math.pow(tlcl / tk0, 1 / KAPPA);
    const ptop = Math.min(...this.ro.profile.p.filter((v, k) => this.edit.T[k] != null));
    // pressures every 2 hPa from the start to the profile's top, with the LCL exactly
    const ps = [];
    for (let q = p0; q > ptop; q -= 2) ps.push(q);
    ps.push(ptop);
    if (plcl < p0 && plcl > ptop) { ps.push(plcl); ps.sort((a, b) => b - a); }
    // dry adiabat (constant mixing ratio) to the LCL, then RK4 steps along the pseudoadiabat
    const trace = [];
    let Tk = tk0;
    for (let k = 0; k < ps.length; k++) {
      const p = ps[k];
      let Tv;
      if (p >= plcl) {
        Tk = tk0 * Math.pow(p / p0, KAPPA);
        Tv = virtualK(Tk, r0);
      } else {
        const prev = ps[k - 1], h = p - prev;
        const k1 = moistLapse(prev, Tk), k2 = moistLapse(prev + h / 2, Tk + (h / 2) * k1);
        const k3 = moistLapse(prev + h / 2, Tk + (h / 2) * k2), k4 = moistLapse(p, Tk + h * k3);
        Tk += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        Tv = virtualK(Tk, satMixing(p, Tk - 273.15));
      }
      trace.push({ p, lp: Math.log10(p), Tv });
    }
    // buoyancy (parcel minus environment virtual temperature) through the trace
    for (const t of trace) { const env = this._envTv(t.lp); t.env = env; t.b = env == null ? null : t.Tv - env; }
    const pts = trace.filter((t) => t.b != null);
    const cross = (a, b) => a.lp + (0 - a.b) * (b.lp - a.lp) / (b.b - a.b);
    let lfc = null, el = null, cape = 0, cin = 0;
    for (let k = 0; k < pts.length - 1; k++) {
      const a = pts[k], b = pts[k + 1];
      if (a.p > plcl + 1e-6 && b.p > plcl + 1e-6) continue; // the LFC is at or above the LCL
      if (lfc == null && a.b <= 0 && b.b > 0) lfc = cross(a, b);
      if (lfc == null && k === 0 && a.b > 0 && a.p <= plcl) lfc = a.lp;
      if (a.b > 0 && b.b <= 0) el = cross(a, b);
    }
    if (lfc == null && pts.some((t) => t.p <= plcl && t.b > 0)) lfc = pts.find((t) => t.p <= plcl && t.b > 0).lp;
    for (let k = 0; k < pts.length - 1; k++) {
      const a = pts[k], b = pts[k + 1], dlnp = Math.log(a.p / b.p), area = RD * ((a.b + b.b) / 2) * dlnp;
      if (lfc != null && b.lp <= lfc && (el == null || a.lp >= el) && area > 0) cape += area;
      if (lfc != null && a.lp > lfc && area < 0) cin += area;
    }
    if (lfc == null) cin = 0;
    const at500 = pts.find((t) => t.p <= 500);
    const li = at500 && at500.env != null ? at500.env - at500.Tv : null;
    return { p0, T0, Td0, trace: pts, lcl: { p: plcl, lp: Math.log10(plcl) }, lfc, el, cape, cin, li };
  }

  // a storm motion at a whole degree and knot, so what's previewed is what a redraw uses
  _snapMotion(u, v) {
    const dir = Math.round(windDir(u, v)) % 360, speed = Math.round(Math.hypot(u, v));
    const rad = (dir * Math.PI) / 180;
    return [-speed * Math.sin(rad), -speed * Math.cos(rad)];
  }

  // ---- precipitation type: follow precipitation down from a level, as SounderPy's winter soundings do ----
  // It starts as snow, drizzle, or rain from the temperature at the top of the cloud the level is in (the level's
  // own when it isn't in cloud): below -10 °C snow, 0 to -10 °C drizzle (little ice), above 0 °C rain (Huffman &
  // Norman 1988; Rauber et al. 2000). Then each melting and refreezing layer below changes it by its energy, from the
  // wet-bulb profile (Bourgouin 2000; the revised method of Birk et al. 2021); liquid on a surface at or below 0 °C
  // freezes. An approximation of SounderPy's: the figure's precipitation type is the exact value.
  _ptypeFor(m) {
    const key = `${m.lp}:${JSON.stringify(this.edit.T)}:${JSON.stringify(this.edit.Td)}`;
    if (this._ptypeCache?.key !== key) this._ptypeCache = { key, r: this._ptype(m.lp) };
    return this._ptypeCache.r;
  }

  _ptype(lp0) {
    const Z = this.ro.profile.z_agl, src = this._level(lp0);
    if (!src || src.z_agl == null) return null;
    const zs = src.z_agl, rh = rhCloud(src.T, src.Td), saturated = rh != null && rh >= 90;
    // in cloud, the cloud's top temperature decides whether it makes ice
    let startT = src.T, where = "here";
    if (saturated) {
      let top = src;
      for (let z = zs + 25; z <= Z[this.n - 1]; z += 25) {
        const l = this._level(this._lpAtHeight(z));
        if (!l || !(rhCloud(l.T, l.Td) >= 90)) break;
        top = l;
      }
      startT = top.T;
      where = top === src ? "here" : "cloud top";
    }
    let phase = startT > 0 ? "RA" : startT > -10 ? "DZ" : "SN";
    const trace = [{ z: zs, phase }], fromIce = phase === "SN";

    // the wet-bulb profile below the source, at the profile's levels, with its exact 0 °C crossings
    const pts = [];
    for (let i = 0; i < this.n; i++) {
      if (Z[i] == null || Z[i] >= zs) continue;
      const tw = this._level(this.lp[i])?.Tw;
      if (tw != null) pts.push([Z[i], tw]);
    }
    if (src.Tw != null) pts.push([zs, src.Tw]);
    const z = [], tw = [];
    pts.forEach(([zb, tb], k) => {
      const [za, ta] = pts[k - 1] ?? [];
      if (k && ta * tb < 0) { z.push(za + ((zb - za) * ta) / (ta - tb)); tw.push(0); }
      z.push(zb); tw.push(tb);
    });
    // melting (> 0 °C) and refreezing (< 0 °C) layers, each with its energy g ∫ (Tw - T0) / T0 dz (J/kg)
    let layers = [];
    for (let k = 0; k < z.length - 1; k++) {
      const mean = (tw[k] + tw[k + 1]) / 2;
      if (mean === 0) continue;
      const kind = mean > 0 ? "warm" : "cold", energy = (9.81 * mean * (z[k + 1] - z[k])) / 273.15, last = layers[layers.length - 1];
      if (last?.kind === kind) { last.top = z[k + 1]; last.energy += energy; }
      else layers.push({ kind, bottom: z[k], top: z[k + 1], energy });
    }
    // top first, without slivers (under 50 m with under 1 J/kg)
    layers = layers.filter((L) => L.top - L.bottom >= 50 || Math.abs(L.energy) >= 1).reverse();

    // melting: ice stays frozen, partly melts, or melts (Bourgouin 2000 at the surface; Birk et al. 2021 aloft);
    // refreezing: melted ice refreezes to sleet, a mix, or stays liquid (Birk et al. 2021)
    let melt = 0;
    for (const L of layers) {
      const before = phase;
      if (L.kind === "warm" && (phase === "SN" || phase === "PL")) {
        melt = L.energy;
        const [partly, fully] = L.bottom <= 1 ? [5.6, 13.2] : [9.4, 15];
        phase = melt < partly ? phase : melt <= fully ? `RA/${phase}` : "RA";
      } else if (L.kind === "cold") {
        if (phase.startsWith("RA/")) phase = "PL";
        else if (phase === "RA" && fromIce && melt > 0) {
          if (-L.energy > 208 + 0.08 * melt) phase = "PL";
          else if (-L.energy >= 117 + 0.08 * melt) phase = "FZRA/PL";
        }
      }
      L.phase = phase;
      L.changed = phase !== before;
      if (L.changed) trace.push({ z: L.bottom, phase });
    }
    const sfcT = this.edit.T.find((t) => t != null);
    if ((phase === "RA" || phase === "DZ") && sfcT <= 0) { phase = `FZ${phase}`; trace.push({ z: 0, phase }); }
    // the layers that matter: melting layers, and refreezing layers under one (above any melting, the cold air only
    // keeps snow snow)
    const firstWarm = layers.findIndex((L) => L.kind === "warm");
    const shown = firstWarm < 0 ? [] : layers.slice(firstWarm);
    return { lp: lp0, z: zs, src, rh, saturated, startT, where, layers, shown, trace, surface: phase };
  }

  // ---- cloud layers: RH ≥ 90% (over ice below 0 °C), at least 200 m deep, as SounderPy's winter soundings find them ----
  _clouds() {
    const key = `${JSON.stringify(this.edit.T)}:${JSON.stringify(this.edit.Td)}`;
    if (this._cloudCache?.key === key) return this._cloudCache.layers;
    const top = this.ro.profile.z_agl[this.n - 1] ?? 0, layers = [];
    let start = null, last = null;
    for (let z = 0; z <= top; z += 25) {
      const l = this._level(this._lpAtHeight(z)), wet = Boolean(l) && rhCloud(l.T, l.Td) >= 90;
      if (wet) { start ??= z; last = z; }
      if ((!wet || z + 25 > top) && start != null) { if (last - start >= 200) layers.push([start, last]); start = null; }
    }
    this._cloudCache = { key, layers };
    return layers;
  }

  _drawClouds(ctx) {
    const S = this.ro.panels.skewt;
    if (!S) return;
    this._clip(ctx, S.bbox, () => {
      for (const [a, b] of this._clouds()) {
        const y0 = this._fwd(S.transform, 0, this._lpAtHeight(a))[1], y1 = this._fwd(S.transform, 0, this._lpAtHeight(b))[1];
        ctx.fillStyle = this.ro.dark ? "rgba(200, 200, 210, 0.16)" : "rgba(110, 120, 135, 0.16)";
        ctx.fillRect(S.bbox[0], y1, S.bbox[2] - S.bbox[0], y0 - y1);
        ctx.strokeStyle = this.ro.dark ? "rgba(200, 200, 210, 0.6)" : "rgba(90, 100, 115, 0.6)"; ctx.lineWidth = 1.5; ctx.setLineDash([4, 4]);
        for (const y of [y0, y1]) { ctx.beginPath(); ctx.moveTo(S.bbox[0], y); ctx.lineTo(S.bbox[2], y); ctx.stroke(); }
        ctx.setLineDash([]);
        ctx.font = "bold 18px system-ui, sans-serif"; ctx.fillStyle = ctx.strokeStyle;
        ctx.fillText("CLOUD", S.bbox[0] + 70, Math.min(y1 + 20, y0 - 4));
      }
    });
  }

  // ---- arrow keys: step the readout (or an unpinned parcel or precipitation source) a level up or down ----
  _onKey(e) {
    if (!this.options.keys || this.layer.hidden || !this.ro || (e.key !== "ArrowUp" && e.key !== "ArrowDown")) return;
    if (e.ctrlKey || e.metaKey || e.altKey || e.target.closest?.("input, textarea, select, [contenteditable]")) return;
    const up = e.key === "ArrowUp";
    // the next profile level above (lower pressure) or below
    const step = (lp) => {
      const next = up ? this.lp.findIndex((l) => l < lp - 1e-6) : this.lp.findLastIndex((l) => l > lp + 1e-6);
      return next < 0 ? null : next;
    };
    const m = this.mark;
    if (m && !m.picked && !this.marking && (m.kind === "parcel" || m.kind === "ptype")) {
      let i = step(m.kind === "parcel" ? this.lp[m.i] : m.lp);
      if (m.kind === "parcel") while (i != null && i >= 0 && i < this.n && (this.edit.T[i] == null || this.edit.Td[i] == null)) i += up ? 1 : -1;
      if (i == null || i < 0 || i >= this.n) return;
      if (m.kind === "parcel") m.i = i; else m.lp = this.lp[i];
    } else if (this.hover?.lp != null && !this.drag && !this.marking) {
      const i = step(this.hover.lp);
      if (i == null) return;
      const hover = { ...this.hover, lp: this.lp[i], cursorT: null, cursorUV: null };
      if (hover.panel === "height") hover.z = this.ro.profile.z_agl[i];
      this.hover = hover;
      this.options.onHover?.(this._level(hover.lp));
    } else return;
    // the readout's, not the page's (or a zoomable view's) arrow keys
    e.preventDefault();
    e.stopPropagation();
    this._draw();
  }

  _clampZ(z) {
    const Z = this.ro.profile.z_agl;
    return Math.min(Math.max(z, 0), Z[this.n - 1]);
  }

  // the height (m AGL) of the skew-t row at y, kept within the profile; within a few pixels of a standard level
  // (1000, 850, 700, 500 hPa, ...) it snaps to it, so a layer's depth is that thickness
  _zFromSkew(x, y) {
    const snap = this._snapP(y);
    return this._zAtLp(this._clampLp(snap ? Math.log10(snap) : this._inv(this.ro.panels.skewt.transform, x, y)[1]));
  }

  // the standard level (hPa) nearest y on the skew-t, within about 12 screen pixels and inside the profile, or null
  _snapP(y) {
    const S = this.ro.panels.skewt, tol = 12 / Math.max(this._screenScale, 0.05);
    const psfc = this.ro.profile.p[0], ptop = this.ro.profile.p[this.n - 1];
    let best = null, bd = tol;
    for (const p of SNAP_LEVELS) {
      const d = Math.abs(this._fwd(S.transform, 0, Math.log10(p))[1] - y);
      if (p <= psfc && p >= ptop && d <= bd) { best = p; bd = d; }
    }
    return best;
  }

  // log10(p) kept within the profile
  _clampLp(lp) {
    const lo = Math.min(this.lp[0], this.lp[this.n - 1]), hi = Math.max(this.lp[0], this.lp[this.n - 1]);
    return Math.min(Math.max(lp, lo), hi);
  }

  // the wind (kt) at a height: from the 100 m hodograph trace, or the (edited) profile levels
  _windAt(z) {
    const tr = this.ro.hodograph_trace;
    if (tr && !this.isEdited() && z <= tr.z_agl[tr.z_agl.length - 1]) {
      const Z = tr.z_agl;
      for (let i = 0; i < Z.length - 1; i++) {
        if (Z[i] <= z && z <= Z[i + 1]) {
          const t = (z - Z[i]) / (Z[i + 1] - Z[i] || 1);
          return [tr.u[i] + t * (tr.u[i + 1] - tr.u[i]), tr.v[i] + t * (tr.v[i + 1] - tr.v[i])];
        }
      }
    }
    const lp = this._lpAtHeight(z);
    if (lp == null) return null;
    const u = this._interp(this.edit.u, lp), v = this._interp(this.edit.v, lp);
    return u == null || v == null ? null : [u, v];
  }

  // heights every 100 m from z0 to z1 (with both ends), the hodograph trace's spacing
  _heights(z0, z1, step = 100) {
    const out = [z0];
    for (let z = Math.floor(z0 / step) * step + step; z < z1; z += step) out.push(z);
    out.push(z1);
    return out;
  }

  // bulk shear (kt) and storm-relative helicity (m²/s², with the figure's storm motion) from z0 to z1
  _kinematics(z0, z1, sm = this.ro.storm_motion) {
    const winds = this._heights(z0, z1).map((z) => this._windAt(z));
    if (winds.some((w) => !w)) return null;
    const [a, b] = [winds[0], winds[winds.length - 1]];
    const out = { shear: Math.hypot(b[0] - a[0], b[1] - a[1]) };
    out.mean = [winds.reduce((s, w) => s + w[0], 0) / winds.length, winds.reduce((s, w) => s + w[1], 0) / winds.length];
    if (sm) {
      // storm-relative wind speed: the mean of the layer's bottom and top, as SounderPy's tables give it
      out.srw = (Math.hypot(a[0] - sm[0], a[1] - sm[1]) + Math.hypot(b[0] - sm[0], b[1] - sm[1])) / 2;
      let srh = 0;
      for (let i = 0; i < winds.length - 1; i++) {
        const [u1, v1] = winds[i], [u2, v2] = winds[i + 1];
        srh += ((u2 - sm[0]) * (v1 - sm[1]) - (u1 - sm[0]) * (v2 - sm[1])) * KT * KT;
      }
      out.srh = srh;
    }
    return out;
  }

  _markResults(m) {
    if (m.kind === "ptype") {
      const r = this._ptypeFor(m);
      if (!r) return null;
      const name = (code) => PTYPES[code]?.[0] ?? code;
      const rows = [
        ["From", `${fmt(r.src.p, 0)} hPa · ${this._height(r.z, " AGL")}`],
        ["Saturated", r.rh == null ? "–" : `${r.saturated ? "Yes" : "No"} · RH ${fmt(Math.min(r.rh, 100), 0)}%${r.src.T < 0 ? " over ice" : ""}`],
        ["Starts as", name(r.trace[0].phase)],
        [r.where === "cloud top" ? "Cloud top" : "Temp here", this._temp(r.startT)],
        ...r.shown.map((L) => [`${L.kind === "warm" ? "Melting" : "Refreezing"} ${this._heightRange(L.bottom, L.top)}`,
          `${fmt(L.energy, 0, " J/kg")}${L.changed ? ` → ${name(L.phase)}` : ""}`]),
        ...(r.shown.length ? [] : [["Melting layers", "None"]]),
        ["At the surface", name(r.surface)],
      ];
      return { kind: "ptype", picked: Boolean(m.picked), p: r.src.p, z: r.z, saturated: r.saturated, start: r.trace[0].phase,
        surface: r.surface, layers: r.layers.map(({ kind, bottom, top, energy, phase }) => ({ kind, bottom, top, energy, phase })), rows };
    }
    if (m.kind === "parcel") {
      const r = this._parcelFor(m), source = this._parcelSource(m);
      const z = (lp) => (lp == null ? null : this._zAtLp(Math.min(Math.max(lp, Math.min(...this.lp)), Math.max(...this.lp))));
      const level = (lp) => (lp == null ? "–" : `${fmt(10 ** lp, 0)} hPa · ${this._height(z(lp), " AGL")}`);
      return {
        kind: "parcel", i: m.i, preset: m.preset ?? null, picked: Boolean(m.picked), cape: r.cape, cin: r.cin, trace: r.trace,
        lcl: 10 ** r.lcl.lp, lfc: r.lfc == null ? null : 10 ** r.lfc, el: r.el == null ? null : 10 ** r.el,
        rows: [
          ...(source.label ? [["Parcel", source.label], ["Starts at", level(Math.log10(source.p0))]] : [["Lifted from", level(Math.log10(source.p0))]]),
          ["Parcel T / Td", `${this._temp(r.T0)} / ${this._temp(r.Td0)}`],
          ["CAPE", fmt(r.cape, 0, " J/kg")],
          ["CIN", fmt(r.cin, 0, " J/kg")],
          ["LCL", level(r.lcl.lp)],
          ["LFC", level(r.lfc)],
          ["EL", level(r.el)],
          ["Lifted index", r.li == null ? "–" : fmt(r.li, 1, " °C")],
          // the pointer's height: CAPE from the LFC up to it
          ...(m.picked && m.cursor != null ? (() => {
            const c = this._parcelTo(r, m.cursor);
            const below = r.lfc == null || c.lp > r.lfc;
            return [
              ["At cursor", level(c.lp)],
              ["CAPE to cursor", below ? "0 J/kg (below LFC)" : fmt(c.cape, 0, " J/kg")],
              ["Share of CAPE", r.cape > 0 ? fmt((100 * c.cape) / r.cape, 0, "%") : "–"],
              ["Buoyancy (Tv)", fmt(c.b, 1, " °C")],
            ];
          })() : []),
        ],
      };
    }
    if (m.kind === "storm") {
      // storm motions read like the figure's: the direction the storm moves toward
      const motion = ([u, v]) => {
        const toward = (windDir(u, v) + 180) % 360;
        return `${compass(toward)} (${fmt(toward, 0, "°")}) @ ${this._speed(Math.hypot(u, v))}`;
      };
      const sm = [m.u, m.v];
      const srh = (z0, z1) => { const k = this._kinematics(z0, z1, sm); return k?.srh == null ? "–" : fmt(k.srh, 0, " m²/s²"); };
      const srw = (z0, z1) => { const k = this._kinematics(z0, z1, sm); return k?.srw == null ? "–" : this._speed(k.srw); };
      const rows = [["Moving toward", motion(sm)], ["0–500 m SRH", srh(0, 500)], ["0–1 km SRH", srh(0, 1000)], ["0–3 km SRH", srh(0, 3000)]];
      const eil = (this.ro.params || []).find((q) => q.key === "kinem.srh_eil")?.layer?.z;
      if (eil) rows.push(["Effective SRH", srh(eil[0], eil[1])]);
      rows.push(["0–1 km SR wind", srw(0, 1000)], ["4–6 km SR wind", srw(4000, 6000)]);
      if (this.ro.storm_motion) rows.push(["Figure's motion", motion(this.ro.storm_motion)]);
      return { kind: "storm", u: m.u, v: m.v, direction: windDir(m.u, m.v), speed: Math.hypot(m.u, m.v), picked: Boolean(m.picked), rows };
    }
    if (m.kind === "layer") {
      const z0 = Math.min(m.z0, m.z1), z1 = Math.max(m.z0, m.z1);
      const lp0 = this._lpAtHeight(z0), lp1 = this._lpAtHeight(z1);
      if (lp0 == null || lp1 == null) return null;
      const bottom = this._level(lp0), top = this._level(lp1);
      const depth = z1 - z0;
      // values every ~1% of the layer for means and extremes
      const samples = Array.from({ length: 41 }, (_, i) => this._level(this._lpAtHeight(z0 + (depth * i) / 40))).filter(Boolean);
      const mean = (k) => { const v = samples.map((l) => l[k]).filter((x) => x != null && Number.isFinite(x)); return v.length ? v.reduce((s, x) => s + x, 0) / v.length : null; };
      const thetae = samples.map((l) => l.thetae).filter((x) => x != null && Number.isFinite(x));
      // precipitable water in the layer: sum of q dp / g, in inches
      let pw = 0;
      for (let i = 0; i < samples.length - 1; i++) {
        const q = (l) => (l.mr ?? 0) / 1000 / (1 + (l.mr ?? 0) / 1000);
        pw += ((q(samples[i]) + q(samples[i + 1])) / 2) * (samples[i].p - samples[i + 1].p) * 100 / 9.81;
      }
      const kin = depth > 0 ? this._kinematics(z0, z1) : null;
      const rows = [
        ["Height", this._heightRange(z0, z1, " AGL")],
        ["Pressure", `${fmt(bottom.p, 0)}–${fmt(top.p, 0)} hPa`],
        ["Depth", `${this._height(depth)} (${fmt(depth / 10, 0)} dam)`],
        ...(m.snap0 && m.snap1 && m.snap0 !== m.snap1
          ? [["Thickness", `${Math.max(m.snap0, m.snap1)}–${Math.min(m.snap0, m.snap1)} hPa: ${fmt(depth / 10, 0)} dam`]] : []),
        // in virtual temperature, as SounderPy's (SHARPpy's) lapse rates are
        ["Lapse rate", depth > 0 ? fmt(-((top.Tv ?? top.T) - (bottom.Tv ?? bottom.T)) / (depth / 1000), 1, " °C/km") : "–"],
        ["ΔT", fmt(top.T - bottom.T, 1, " °C")],
        ["Mean RH", fmt(mean("rh"), 0, "%")],
        ["Mean mix ratio", fmt(mean("mr"), 1, " g/kg")],
        ["PW in layer", fmt(pw / 25.4, 2, " in")],
        ["θe min / max", thetae.length ? `${fmt(Math.min(...thetae), 1)} / ${fmt(Math.max(...thetae), 1)} K` : "–"],
      ];
      if (kin) {
        rows.push(["Mean wind", this._windText(...kin.mean)]);
        rows.push(["Bulk shear", this._speed(kin.shear)]);
        if (kin.srh != null) rows.push(["SRH", fmt(kin.srh, 0, " m²/s²")]);
      }
      return { kind: "layer", z: [z0, z1], p: [bottom.p, top.p], rows };
    }
    if (m.panel === "skewt") {
      const S = this.ro.panels.skewt;
      const pt = ([x, y]) => {
        const [T, lp] = this._inv(S.transform, x, y), p = 10 ** lp;
        return { T, p, z: this._zAtLp(lp), theta: (T + 273.15) * Math.pow(1000 / p, 0.2854) - 273.15 };
      };
      const a = pt(m.a), b = pt(m.b);
      const dz = a.z != null && b.z != null ? b.z - a.z : null;
      const rows = [
        ["A", `${fmt(a.p, 0)} hPa, ${this._temp(a.T)}`],
        ["B", `${fmt(b.p, 0)} hPa, ${this._temp(b.T)}`],
        ["ΔT", fmt(b.T - a.T, 1, " °C")],
        ["Δz", dz == null ? "–" : this._height(dz)],
        ["Lapse rate", dz ? fmt(-(b.T - a.T) / (dz / 1000), 1, " °C/km") : "–"],
        ["Δθ", fmt(b.theta - a.theta, 1, " K")],
        ["Δp", fmt(b.p - a.p, 0, " hPa")],
      ];
      return { kind: "measure", panel: "skewt", a, b, rows };
    }
    const a = this._hodoUV(...m.a), b = this._hodoUV(...m.b);
    const d = [b[0] - a[0], b[1] - a[1]];
    const wind = ([u, v]) => this._windText(u, v);
    return {
      kind: "measure", panel: "hodograph", a, b,
      rows: [["A", wind(a)], ["B", wind(b)], ["Difference", this._speed(Math.hypot(...d))],
             ["Toward", `${fmt((Math.atan2(d[0], d[1]) * 180 / Math.PI + 360) % 360, 0, "°")}`]],
    };
  }

  _drawMark(ctx) {
    const m = this.mark;
    if (m.kind === "ptype") {
      const S = this.ro.panels.skewt, r = this._ptypeFor(m);
      if (!r) { this._placeCard(); return; }
      const at = (T, z) => { const lp = this._lpAtHeight(z); return lp == null ? null : this._fwd(S.transform, T, lp); };
      this._clip(ctx, S.bbox, () => {
        // melting (warm) and refreezing (cold) energy: the area between the wet-bulb profile and 0 °C
        for (const L of r.shown) {
          const zs = this._heights(L.bottom, L.top, 20);
          const tw = zs.map((z) => at(this._level(this._lpAtHeight(z))?.Tw ?? 0, z)).filter(Boolean);
          const zero = zs.map((z) => at(0, z)).filter(Boolean).reverse();
          if (tw.length < 2) continue;
          ctx.fillStyle = L.kind === "warm" ? "rgba(230, 60, 40, 0.32)" : "rgba(40, 110, 230, 0.3)";
          ctx.beginPath(); ctx.moveTo(...tw[0]); [...tw.slice(1), ...zero].forEach((pt) => ctx.lineTo(...pt)); ctx.closePath(); ctx.fill();
        }
        // the source level
        const [, y] = this._fwd(S.transform, 0, r.lp);
        ctx.strokeStyle = "#1f6fd1"; ctx.lineWidth = 2; ctx.setLineDash([8, 6]);
        ctx.beginPath(); ctx.moveTo(S.bbox[0], y); ctx.lineTo(S.bbox[2], y); ctx.stroke(); ctx.setLineDash([]);
        this._dot(ctx, this._fwd(S.transform, r.src.T, r.lp), "#1f6fd1", 8);
        // the phase at the source and after each change, left of 0 °C (a winter figure labels its own on the right)
        ctx.font = "bold 22px system-ui, sans-serif"; ctx.lineWidth = 5; ctx.lineJoin = "round"; ctx.textAlign = "right";
        for (const step of r.trace) {
          const pt = at(0, Math.max(step.z, 25));
          if (!pt) continue;
          const label = (PTYPES[step.phase]?.[0] ?? step.phase).toUpperCase(), x = pt[0] - 14, ty = pt[1] - 6;
          ctx.strokeStyle = this.ro.dark ? "rgba(0,0,0,0.8)" : "rgba(255,255,255,0.9)"; ctx.strokeText(label, x, ty);
          ctx.fillStyle = PTYPES[step.phase]?.[1] ?? "#808080"; ctx.fillText(label, x, ty);
        }
        ctx.textAlign = "start";
      });
      this._placeCard();
      return;
    }
    if (m.kind === "parcel") {
      const S = this.ro.panels.skewt, r = this._parcelFor(m);
      const px = (Tv, lp) => this._fwd(S.transform, Tv - 273.15, lp);
      const cursor = m.picked && m.cursor != null ? Math.max(m.cursor, r.trace[r.trace.length - 1].lp) : null;
      this._clip(ctx, S.bbox, () => {
        // CAPE (warm) and CIN (cool) between the parcel and the environment
        for (let k = 0; k < r.trace.length - 1; k++) {
          const a = r.trace[k], b = r.trace[k + 1];
          if (a.env == null || b.env == null) continue;
          const sign = (a.b + b.b) / 2;
          if (Math.abs(sign) < 0.02) continue;
          // CAPE between the LFC and EL, CIN below the LFC (nothing to shade without an LFC)
          if (r.lfc == null) continue;
          if (sign > 0 && (b.lp > r.lfc || (r.el != null && a.lp < r.el))) continue;
          if (sign < 0 && a.lp <= r.lfc) continue;
          // CAPE from the LFC up to the pointer is darker
          const upTo = sign > 0 && cursor != null && a.lp >= cursor;
          ctx.fillStyle = sign > 0 ? (upTo ? "rgba(230, 80, 40, 0.55)" : "rgba(230, 80, 40, 0.28)") : "rgba(40, 110, 230, 0.22)";
          ctx.beginPath(); ctx.moveTo(...px(a.Tv, a.lp)); ctx.lineTo(...px(b.Tv, b.lp)); ctx.lineTo(...px(b.env, b.lp)); ctx.lineTo(...px(a.env, a.lp)); ctx.closePath(); ctx.fill();
        }
        this._line(ctx, r.trace.map((t) => px(t.Tv, t.lp)), "#8b3fd9", [10, 6], 3.5);
        if (cursor != null) {
          const y = this._fwd(S.transform, 0, cursor)[1];
          ctx.strokeStyle = this.ro.dark ? "rgba(255,255,255,0.8)" : "rgba(0,0,0,0.7)"; ctx.lineWidth = 1.5; ctx.setLineDash([8, 6]);
          ctx.beginPath(); ctx.moveTo(S.bbox[0], y); ctx.lineTo(S.bbox[2], y); ctx.stroke(); ctx.setLineDash([]);
        }
        this._dot(ctx, px(r.trace[0].Tv, r.trace[0].lp), "#8b3fd9", 8);
        // LCL, LFC, and EL ticks with labels
        ctx.font = "bold 22px system-ui, sans-serif"; ctx.fillStyle = "#8b3fd9"; ctx.strokeStyle = "#8b3fd9"; ctx.lineWidth = 3;
        for (const [lp, label] of [[r.lcl.lp, "LCL"], [r.lfc, "LFC"], [r.el, "EL"]]) {
          if (lp == null) continue;
          const t = r.trace.reduce((best, q) => (Math.abs(q.lp - lp) < Math.abs(best.lp - lp) ? q : best), r.trace[0]);
          const [x, y] = px(t.Tv, lp);
          ctx.beginPath(); ctx.moveTo(x + 6, y); ctx.lineTo(x + 26, y); ctx.stroke();
          ctx.fillText(label, x + 30, y + 7);
        }
      });
      this._placeCard();
      return;
    }
    if (m.kind === "storm") {
      const H = this.ro.panels.hodograph, tr = this.ro.hodograph_trace;
      this._clip(ctx, H.bbox, () => {
        const c = this._hodoPx(m.u, m.v);
        // the 0-3 km SRH area: swept between the storm motion and the hodograph
        if (tr) {
          const pts = [];
          for (let i = 0; i < tr.z_agl.length && tr.z_agl[i] <= 3000; i++) pts.push(this._hodoPx(tr.u[i], tr.v[i]));
          if (pts.length > 1) {
            ctx.fillStyle = this.ro.dark ? "rgba(63, 214, 214, 0.25)" : "rgba(10, 165, 165, 0.2)";
            ctx.beginPath(); ctx.moveTo(...c); pts.forEach((pt) => ctx.lineTo(...pt)); ctx.closePath(); ctx.fill();
            ctx.strokeStyle = MARK; ctx.lineWidth = 2; ctx.setLineDash([6, 4]);
            for (const pt of [pts[0], pts[pts.length - 1]]) { ctx.beginPath(); ctx.moveTo(...c); ctx.lineTo(...pt); ctx.stroke(); }
            ctx.setLineDash([]);
          }
        }
        // the figure's storm motion, for comparison
        if (this.ro.storm_motion) {
          const f = this._hodoPx(...this.ro.storm_motion);
          ctx.strokeStyle = this.ro.dark ? "rgba(255,255,255,0.7)" : "rgba(0,0,0,0.6)"; ctx.lineWidth = 2;
          ctx.beginPath(); ctx.arc(...f, 8, 0, 2 * Math.PI); ctx.stroke();
        }
        // the storm motion being tried: a ring with a cross, filled once picked
        ctx.strokeStyle = MARK; ctx.lineWidth = 3;
        ctx.beginPath(); ctx.arc(...c, 11, 0, 2 * Math.PI); ctx.stroke();
        if (m.picked) { ctx.fillStyle = MARK; ctx.fill(); }
        ctx.beginPath(); ctx.moveTo(c[0] - 18, c[1]); ctx.lineTo(c[0] + 18, c[1]); ctx.moveTo(c[0], c[1] - 18); ctx.lineTo(c[0], c[1] + 18); ctx.stroke();
      });
      this._placeCard();
      return;
    }
    if (m.kind === "layer") {
      const z0 = Math.min(m.z0, m.z1), z1 = Math.max(m.z0, m.z1);
      if (z1 - z0 < 1) return;
      this._drawLayer(ctx, { z: [z0, z1] }, MARK);
      // a standard level an end snapped to, tagged at the skew-t's left edge
      const S = this.ro.panels.skewt;
      if (m.panel === "skewt" && S) {
        this._clip(ctx, S.bbox, () => {
          ctx.font = "bold 20px system-ui, sans-serif";
          for (const p of new Set([m.snap0, m.snap1].filter(Boolean))) {
            const y = this._fwd(S.transform, 0, Math.log10(p))[1], text = `${p} hPa`, w = ctx.measureText(text).width + 14;
            ctx.fillStyle = MARK; ctx.fillRect(S.bbox[0] + 4, y - 13, w, 26);
            ctx.fillStyle = "#fff"; ctx.fillText(text, S.bbox[0] + 11, y + 7);
          }
        });
      }
    } else {
      const bbox = this.ro.panels[m.panel].bbox;
      this._clip(ctx, bbox, () => {
        ctx.strokeStyle = MARK; ctx.lineWidth = 3; ctx.setLineDash([]);
        ctx.beginPath(); ctx.moveTo(...m.a); ctx.lineTo(...m.b); ctx.stroke();
        for (const [pt, label] of [[m.a, "A"], [m.b, "B"]]) {
          this._dot(ctx, pt, MARK, 7);
          ctx.fillStyle = MARK; ctx.font = "bold 22px system-ui, sans-serif";
          ctx.fillText(label, pt[0] + 10, pt[1] - 10);
        }
      });
    }
    this._placeCard();
  }

  _placeCard() {
    const res = this.mark && this._markResults(this.mark);
    if (!res) { this.card.style.display = "none"; return; }
    const title = res.kind === "layer" ? "Layer" : res.kind === "storm" ? "Storm motion" : res.kind === "ptype" ? "Precip type (approximate)"
      : res.kind === "parcel" ? `${res.preset ? `${res.preset.toUpperCase()} parcel` : "Parcel"} (approximate)`
      : res.panel === "skewt" ? "Measure (skew-T)" : "Measure (hodograph)";
    const head = res.kind === "storm" ? (res.picked ? "Picked" : `${this._touch ? "Let go" : "Click"} to use this storm motion`)
      : res.kind === "parcel" ? (res.picked ? "Pinned" : `${this._touch ? "Let go" : "Click"} to pin this parcel`)
      : res.kind === "ptype" ? (res.picked ? "Pinned" : `If precipitation fell from here · ${this._touch ? "let go" : "click"} to pin`)
      : this.isEdited() && res.kind === "layer" ? "Edited profile" : "";
    const html = `<p class="spy-readout-title">${title}</p>${head ? `<p class="spy-readout-head">${head}</p>` : ""}`
      + res.rows.map(([k, v]) => `<div><span>${k}</span><b>${v}</b></div>`).join("");
    if (this.panel) {
      const section = this.panel.querySelector("[data-readout-mark]");
      section.innerHTML = html;
      section.hidden = false;
      this.card.style.display = "none";
      return;
    }
    this.card.innerHTML = html;
    // beside where the drag ended
    const rect = this.img.getBoundingClientRect();
    const [x, y] = this.mark.end;
    this._place(this.card, [rect.left + (x * rect.width) / this.ro.image.width, rect.top + (y * rect.height) / this.ro.image.height]);
  }

  _setHover(hover) {
    this.hover = hover;
    this._draw();
    this.options.onHover?.(hover ? this._level(hover.lp) : null);
  }

  // ------------------------------------------------------------------------------------------------
  // drawing, in image pixels
  // ------------------------------------------------------------------------------------------------
  _draw() {
    if (!this.ro) return;
    const ctx = this.ctx, s = this._scale * (this._res || window.devicePixelRatio || 1);
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, this.layer.width, this.layer.height);
    if (!Number.isFinite(s) || s <= 0) return;
    ctx.setTransform(s, 0, 0, s, 0, 0);

    if (this.options.clouds) this._drawClouds(ctx);
    if (this.editing || this.isEdited()) this._drawEdits(ctx);
    if (this.mark) this._drawMark(ctx);
    this._drawDrawings(ctx);
    // with a layer or measurement showing, its card is the readout
    if (this.marking || this.tool === "draw" || (this.mark && MARK_TOOLS.includes(this.tool))) { this.tip.style.display = "none"; return; }

    if (this.hover?.panel === "param") { this._drawParam(ctx, this.hover.param); return; }

    const lvl = this.hover ? this._level(this.hover.lp) : null;
    const z = this.hover?.panel === "height" ? this.hover.z : lvl?.z_agl;
    if (z != null && this._wants("panels")) this._drawHeightLines(ctx, z);
    if (!lvl && this.hover?.panel === "height") { this._showHeightOnlyTip(); return; }
    if (!lvl) { this.tip.style.display = "none"; return; }
    const S = this.ro.panels.skewt, H = this.ro.panels.hodograph;
    const ink = this.ro.dark ? "rgba(255,255,255,0.75)" : "rgba(0,0,0,0.7)";

    if (S) {
      const [, y] = this._fwd(S.transform, 0, this.hover.lp);
      this._clip(ctx, S.bbox, () => {
        ctx.strokeStyle = ink; ctx.lineWidth = 1.5; ctx.setLineDash([8, 6]);
        ctx.beginPath(); ctx.moveTo(S.bbox[0], y); ctx.lineTo(S.bbox[2], y); ctx.stroke();
        ctx.setLineDash([]);
        this._dot(ctx, this._fwd(S.transform, lvl.T, this.hover.lp), "#e01b24");
        if (lvl.Td != null) this._dot(ctx, this._fwd(S.transform, lvl.Td, this.hover.lp), "#1a9e2c");
      });
    }
    if (H && lvl.u != null && lvl.z_agl != null && lvl.z_agl <= H.max_height) {
      this._clip(ctx, H.bbox, () => this._dot(ctx, this._hodoPx(lvl.u, lvl.v), "#7b2cbf", 9));
    }
    this._showTip(lvl);
  }

  _drawHeightLines(ctx, z) {
    const ink = this.ro.dark ? "rgba(255,255,255,0.75)" : "rgba(0,0,0,0.7)";
    for (const hp of this.ro.height_panels || []) {
      if (z < hp.z_limits[0] || z > hp.z_limits[1]) continue;
      const y = this._panelY(hp, z);
      this._clip(ctx, hp.bbox, () => {
        ctx.strokeStyle = ink; ctx.lineWidth = 1.5; ctx.setLineDash([6, 4]);
        ctx.beginPath(); ctx.moveTo(hp.bbox[0], y); ctx.lineTo(hp.bbox[2], y); ctx.stroke();
        ctx.setLineDash([]);
      });
    }
  }

  // a panel's values at a height, interpolated between its 100 m values
  _panelValues(hp, z) {
    const Z = hp.z;
    for (let i = 0; i < Z.length - 1; i++) {
      if (Z[i] <= z && z <= Z[i + 1]) {
        const t = (z - Z[i]) / (Z[i + 1] - Z[i] || 1);
        return hp.series.map((s) => [s, s.values[i] == null || s.values[i + 1] == null ? null : s.values[i] + t * (s.values[i + 1] - s.values[i])]);
      }
    }
    return hp.series.map((s) => [s, null]);
  }

  _paramValue(param) {
    const v = param.value;
    if (v == null) return "–";
    if (param.kind === "text") return String(v);
    if (param.kind === "wind") return this._windText(v[0], v[1]);
    // speeds, heights, and temperatures in the readout's units
    if (param.units === "kt") return this._speed(v, Math.abs(v) >= 10 ? 0 : 1);
    if (param.units === "m") return this._height(v);
    if (param.units === "°C") return this._temp(v);
    if (param.units === "°F") return this._temp(((v - 32) * 5) / 9);
    const a = Math.abs(v);
    const digits = param.units === "/s" ? 4 : a >= 100 ? 0 : a >= 10 ? 1 : 2;
    return fmt(v, digits, param.units ? ` ${param.units}` : "");
  }

  // a parameter's layer, shaded on the skew-t, hodograph, and height panels
  _drawLayer(ctx, layer, color = "#e08a00") {
    // several single heights, such as every freezing level
    if (Array.isArray(layer.z_at)) { for (const z_at of layer.z_at) this._drawLayer(ctx, { z_at }, color); return; }
    const S = this.ro.panels.skewt, H = this.ro.panels.hodograph;
    const band = color === "#e08a00" ? (this.ro.dark ? "rgba(255, 214, 102, 0.22)" : "rgba(255, 170, 0, 0.18)")
      : (this.ro.dark ? "rgba(63, 214, 214, 0.22)" : "rgba(10, 165, 165, 0.16)");
    let zs = null, lps = null;
    if (layer.z) { zs = layer.z; lps = zs.map((z) => this._lpAtHeight(z)); }
    else if (layer.p) { lps = layer.p.map((p) => Math.log10(p)); zs = lps.map((lp) => this._zAtLp(lp)); }
    else if (layer.z_at != null) { zs = [layer.z_at, layer.z_at]; lps = zs.map((z) => this._lpAtHeight(z)); }
    if (!zs) return;
    const line = zs[0] === zs[1];
    const shade = (bbox, y0, y1) => this._clip(ctx, bbox, () => {
      ctx.fillStyle = band; ctx.strokeStyle = color; ctx.lineWidth = 3;
      if (line) { ctx.beginPath(); ctx.moveTo(bbox[0], y0); ctx.lineTo(bbox[2], y0); ctx.stroke(); }
      else {
        ctx.fillRect(bbox[0], Math.min(y0, y1), bbox[2] - bbox[0], Math.abs(y1 - y0));
        // edges, so a thin layer still shows
        ctx.lineWidth = 1.5;
        for (const yy of [y0, y1]) { ctx.beginPath(); ctx.moveTo(bbox[0], yy); ctx.lineTo(bbox[2], yy); ctx.stroke(); }
      }
    });
    if (S && lps[0] != null && lps[1] != null) shade(S.bbox, ...lps.map((lp) => this._fwd(S.transform, 0, lp)[1]));
    const tr = this.ro.hodograph_trace;
    if (H && tr && zs[0] != null && zs[1] != null) {
      const pts = [];
      for (let i = 0; i < tr.z_agl.length; i++) if (tr.z_agl[i] >= zs[0] - 1 && tr.z_agl[i] <= zs[1] + 1) pts.push(this._hodoPx(tr.u[i], tr.v[i]));
      this._clip(ctx, H.bbox, () => {
        if (pts.length > 1) { ctx.globalAlpha = 0.4; this._line(ctx, pts, color, [], 16); ctx.globalAlpha = 1; }
        else if (pts.length === 1) this._dot(ctx, pts[0], color, 8);
      });
    }
    for (const hp of this.ro.height_panels || []) {
      if (zs[0] == null || zs[1] == null || zs[0] > hp.z_limits[1]) continue;
      shade(hp.bbox, this._panelY(hp, zs[0]), this._panelY(hp, Math.min(zs[1], hp.z_limits[1])));
    }
  }

  _drawParam(ctx, param) {
    if (param.layer) this._drawLayer(ctx, param.layer);
    const b = param.bbox;
    ctx.strokeStyle = "#e08a00"; ctx.lineWidth = 2.5;
    ctx.strokeRect(b[0] - 4, b[1] - 3, b[2] - b[0] + 8, b[3] - b[1] + 6);
    const esc = (t) => String(t).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
    const stale = this.isEdited() ? `<p class="spy-readout-head">From the figure; redraw to update</p>` : "";
    const L = param.layer;
    // a layer, or one level ("At 700 hPa")
    const span = (a, b, unit) => (Math.round(a) === Math.round(b) ? `${fmt(a, 0)} ${unit}` : `${fmt(a, 0)}–${fmt(b, 0)} ${unit}`);
    const where = L?.z ? (Math.round(L.z[0]) === Math.round(L.z[1]) ? this._height(L.z[0], " AGL") : this._heightRange(L.z[0], L.z[1], " AGL"))
      : L?.p ? span(L.p[0], L.p[1], "hPa") : "";
    const whereLabel = (L?.z && Math.round(L.z[0]) === Math.round(L.z[1])) || (L?.p && Math.round(L.p[0]) === Math.round(L.p[1])) ? "At" : "Layer";
    const title = `<p class="spy-readout-title">${esc(param.name)}</p>`, value = `<div><span>Value</span><b>${esc(this._paramValue(param))}</b></div>`;
    const full = stale + title + value
      + (where ? `<div><span>${whereLabel}</span><b>${where}</b></div>` : "")
      + (param.extra || []).map(([k, v]) => `<div><span>${esc(k)}</span><b>${esc(v)}</b></div>`).join("")
      + (param.about ? `<p class="spy-readout-about">${esc(param.about)}</p>` : "");
    this._show(full, title + value);
  }

  _showHeightOnlyTip() {
    this._show(this._heightRows(), this._heightRows());
  }

  // the hovered height panel's values, above the level readout
  _heightRows() {
    if (this.hover?.panel !== "height") return "";
    const { hp, z } = this.hover;
    const rows = this._panelValues(hp, z).map(([s, v]) => {
      const a = Math.abs(v ?? 0), digits = s.units === "/s" ? 4 : a >= 100 ? 0 : 1;
      return `<div><span>${s.label}</span><b>${fmt(v, digits, s.units === "%" ? "%" : ` ${s.units}`)}</b></div>`;
    }).join("");
    return `<p class="spy-readout-title">${hp.title} · ${this._height(z, " AGL")}</p>${rows}<hr>`;
  }

  _drawEdits(ctx) {
    const S = this.ro.panels.skewt, H = this.ro.panels.hodograph, edited = this.isEdited();
    const dash = edited ? [] : [5, 5];
    // while editing, fade the figure's own traces so the editable ones stand out
    if (this.editing) {
      ctx.fillStyle = this.ro.dark ? "rgba(0, 0, 0, 0.5)" : "rgba(255, 255, 255, 0.55)";
      for (const P of [S, H]) if (P) ctx.fillRect(P.bbox[0], P.bbox[1], P.bbox[2] - P.bbox[0], P.bbox[3] - P.bbox[1]);
    }
    if (S) {
      this._clip(ctx, S.bbox, () => {
        for (const [field, color] of [["T", "#ff1f8f"], ["Td", "#00a6ff"]]) {
          const pts = this.edit[field].map((t, i) => (t == null ? null : this._skewPx(t, i)));
          this._line(ctx, pts, color, dash);
          if (this.editing) pts.forEach((pt) => pt && this._dot(ctx, pt, color, 5));
        }
      });
    }
    if (H) {
      this._clip(ctx, H.bbox, () => {
        const pts = this._hodoLevels().map((i) => this._hodoPx(this.edit.u[i], this.edit.v[i]));
        this._line(ctx, pts, "#ff1f8f", dash);
        if (this.editing) pts.forEach((pt) => this._dot(ctx, pt, "#ff1f8f", 5));
      });
    }
  }

  _clip(ctx, b, draw) {
    ctx.save();
    ctx.beginPath(); ctx.rect(b[0], b[1], b[2] - b[0], b[3] - b[1]); ctx.clip();
    draw();
    ctx.restore();
  }

  _line(ctx, pts, color, dash, width = 3.5) {
    ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash); ctx.lineJoin = "round"; ctx.lineCap = "round";
    ctx.beginPath();
    let pen = false;
    for (const pt of pts) {
      if (!pt) { pen = false; continue; }
      pen ? ctx.lineTo(pt[0], pt[1]) : ctx.moveTo(pt[0], pt[1]);
      pen = true;
    }
    ctx.stroke();
    ctx.setLineDash([]);
  }

  _dot(ctx, [x, y], color, r = 7) {
    ctx.beginPath(); ctx.arc(x, y, r, 0, 2 * Math.PI);
    ctx.fillStyle = color; ctx.fill();
    ctx.lineWidth = 2; ctx.strokeStyle = this.ro.dark ? "#111" : "#fff"; ctx.stroke();
  }

  _temp(c) {
    if (c == null || !Number.isFinite(c)) return "–";
    const f = (c * 9) / 5 + 32, mode = this.options.temperature;
    if (mode === "F") return `${num(f)} °F`;
    if (mode === "C") return `${num(c, 1)} °C`;
    return `${num(c, 1)} °C (${num(f)} °F)`;
  }

  // a wind speed (given in kt) in the readout's units
  _speed(kt, digits = 0) {
    if (kt == null || !Number.isFinite(kt)) return "–";
    if (this.options.wind === "m/s") return `${num(kt * KT, digits)} m/s`;
    if (this.options.wind === "mph") return `${num(kt * 1.150779, digits)} mph`;
    return `${num(kt, digits)} kt`;
  }

  // a wind (u, v in kt): the direction it comes from, and its speed
  _windText(u, v) {
    if (u == null || v == null) return "–";
    return Math.hypot(u, v) < 0.5 ? "calm" : `${fmt(windDir(u, v), 0, "°")} @ ${this._speed(Math.hypot(u, v))}`;
  }

  // a height range (m) in one unit, e.g. "0–3000 m AGL"
  _heightRange(a, b, suffix = "") {
    if (a == null || b == null) return "–";
    const ft = this.options.height === "ft", f = ft ? 3.28084 : 1;
    return `${num(a * f)}–${num(b * f)} ${ft ? "ft" : "m"}${suffix}`;
  }

  // "0–2.7 km" (or kft) for cumulative rows
  _depthLabel(z) {
    return this.options.height === "ft" ? `0–${num((z * 3.28084) / 1000, 1)} kft` : `0–${num(z / 1000, 1)} km`;
  }

  _height(m, suffix = "") {
    if (m == null || !Number.isFinite(m)) return "–";
    return this.options.height === "ft" ? `${num(m * 3.28084)} ft${suffix}` : `${num(m)} m${suffix}`;
  }

  // [label, value] rows for a level, for the given fields
  _levelRows(l, fields, cursor) {
    const rows = [];
    for (const key of fields) {
      let value;
      if (key === "p") value = fmt(l.p, 0, " hPa");
      else if (key === "z_agl") value = this._height(l.z_agl, " AGL");
      else if (key === "z") value = this._height(l.z, " MSL");
      else if (["T", "Td", "Tw", "Tv"].includes(key)) value = this._temp(l[key]);
      else if (key === "rh") value = fmt(l.rh, 0, "%");
      else if (key === "mr") value = fmt(l.mr, 1, " g/kg");
      else if (key === "theta" || key === "thetae") value = fmt(l[key], 1, " K");
      else if (key === "wind") value = l.u == null ? "–" : this._windText(l.u, l.v);
      else if (key === "omega") { if (l.omega == null) continue; value = fmt(l.omega, 2, " Pa/s"); }
      else if (key === "cumulative") {
        // shear and SRH from the surface up to this level, within the hodograph's height
        const top = this.ro.panels.hodograph?.max_height ?? 12000;
        if (l.z_agl == null || l.z_agl < 50 || l.z_agl > top) continue;
        const kin = this._kinematics(0, l.z_agl);
        if (!kin) continue;
        const depth = this._depthLabel(l.z_agl);
        rows.push([`${depth} shear`, this._speed(kin.shear)]);
        if (kin.srh != null) rows.push([`${depth} SRH`, fmt(kin.srh, 0, " m²/s²")]);
        continue;
      }
      else continue;
      rows.push([LABELS[key] ?? key, value]);
    }
    if (cursor && this.hover.cursorT != null) rows.push(["Cursor T", this._temp(this.hover.cursorT)]);
    if (cursor && this.hover.cursorUV) {
      const [u, v] = this.hover.cursorUV;
      rows.push(["Cursor", this._windText(u, v)]);
    }
    return rows.map(([k, v]) => `<div><span>${k}</span><b>${v}</b></div>`).join("");
  }

  // the full readout to the panel (if there is one) and a short one at the pointer, or the full one at the pointer
  _show(full, short) {
    const footer = this._footer();
    const section = this.panel?.querySelector("[data-readout-cursor]");
    if (!section) { this.tip.innerHTML = full + footer; this._placeTip(); return; }
    section.innerHTML = full;
    if (short && this.options.panelTip) { this.tip.innerHTML = short + footer; this._placeTip(); } else this.tip.style.display = "none";
  }

  _footer() {
    const text = this.options.tipFooter;
    if (!text) return "";
    const p = document.createElement("p");
    p.className = "spy-readout-foot";
    p.textContent = text;
    return p.outerHTML;
  }

  _showTip(l) {
    const head = l.edited ? `<p class="spy-readout-head">Edited profile</p>` : "";
    const full = head + this._heightRows() + this._levelRows(l, this.options.fields, this.options.cursor);
    this._show(full, this.panel ? this._levelRows(l, this.options.tipFields, false) : null);
  }

  _placeTip() {
    this._place(this.tip, this._client);
  }

  // beside a screen point, inside the part of the image's container that's on screen, with a margin
  _place(box, client) {
    box.style.display = "block";
    const parent = box.offsetParent || this.img.parentElement;
    const rect = parent.getBoundingClientRect(), margin = 8, gap = 16;
    const area = {
      left: Math.max(rect.left, 0) + margin, top: Math.max(rect.top, 0) + margin,
      right: Math.min(rect.right, window.innerWidth) - margin, bottom: Math.min(rect.bottom, window.innerHeight) - margin,
    };
    const [cx, cy] = client ?? [area.left, area.top];
    const w = box.offsetWidth, h = box.offsetHeight;
    // right of and below the point, or flipped to the left and above when it won't fit
    let x = cx + gap, y = cy + gap;
    if (x + w > area.right) x = cx - gap - w;
    if (x < area.left) x = Math.min(Math.max(area.left, cx - w / 2), area.right - w);
    if (y + h > area.bottom) y = cy - gap - h;
    if (y < area.top) y = Math.min(Math.max(area.top, area.bottom - h), cy + gap);
    box.style.left = `${Math.max(area.left, x) - rect.left - parent.clientLeft}px`;
    box.style.top = `${Math.max(area.top, y) - rect.top - parent.clientTop}px`;
  }
}

function injectStyle() {
  if (document.getElementById(STYLE_ID)) return;
  const style = document.createElement("style");
  style.id = STYLE_ID;
  style.textContent = STYLE;
  document.head.append(style);
}

// for pages that load this as a plain <script type="module"> and use the global
if (typeof window !== "undefined") window.SounderPyReadout ??= SounderPyReadout;

export default SounderPyReadout;

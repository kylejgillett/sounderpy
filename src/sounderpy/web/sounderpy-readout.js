/*!
 * SounderPy readout — hover readouts and profile editing for SounderPy figures on the web.
 *
 * Reads the readout file SounderPy writes next to a figure (build_sounding(..., save=True, readout=True)
 * writes kokc.png and kokc.readout.json) and lays a transparent layer over the <img>:
 *   - hover the skew-t or hodograph to read the profile at that level, linked across both panels
 *   - optionally drag temperature, dewpoint, and wind points to edit the profile, then send
 *     getProfile() to a server that redraws it with sounderpy.apply_profile_edits()
 *
 * No dependencies. Works with a responsive or scaled image, in light or dark mode, with mouse or touch.
 *
 *   import { SounderPyReadout } from "./sounderpy-readout.js";
 *   const readout = await SounderPyReadout.attach(img, "kokc.readout.json");
 *
 * Docs: https://sounderpy.app/docs/reference/webreadouts.html
 * Part of SounderPy (https://github.com/kylejgillett/sounderpy), MIT License, (C) Kyle J Gillett 2026
 */

const FORMAT = "sounderpy-readout";
const VERSION = 1;

const DEFAULT_FIELDS = ["p", "z_agl", "T", "Td", "Tw", "rh", "mr", "thetae", "wind", "omega"];

const LABELS = {
  p: "Pressure", z_agl: "Height", z: "Height MSL", T: "Temp", Td: "Dewpoint", Tw: "Wet-bulb", Tv: "Virtual T",
  rh: "RH", mr: "Mix ratio", theta: "θ", thetae: "θe", wind: "Wind", omega: "ω",
};

const STYLE_ID = "spy-readout-style";
const STYLE = `
.spy-readout-layer { position: absolute; z-index: 2; cursor: crosshair; touch-action: pan-y pinch-zoom; }
.spy-readout-layer[data-editing] { touch-action: none; }
.spy-readout-tip {
  --spy-readout-bg: rgba(255, 255, 255, 0.96); --spy-readout-text: #111; --spy-readout-muted: #555;
  --spy-readout-border: rgba(0, 0, 0, 0.35); --spy-readout-edited: #c2187a;
  position: absolute; z-index: 3; display: none; pointer-events: none; box-sizing: border-box;
  min-width: min(190px, calc(100vw - 16px)); max-width: calc(100vw - 16px); padding: 6px 10px; border: 1px solid var(--spy-readout-border); border-radius: 6px;
  background: var(--spy-readout-bg); color: var(--spy-readout-text);
  font: 12px/1.45 ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
  box-shadow: 0 2px 10px rgba(0, 0, 0, 0.25); text-align: left; white-space: nowrap;
}
.spy-readout-tip[data-dark] {
  --spy-readout-bg: rgba(24, 26, 30, 0.95); --spy-readout-text: #eee; --spy-readout-muted: #aaa;
  --spy-readout-border: rgba(255, 255, 255, 0.3); --spy-readout-edited: #ff6ec7;
}
.spy-readout-tip div { display: flex; justify-content: space-between; gap: 14px; }
.spy-readout-tip span { color: var(--spy-readout-muted); }
.spy-readout-tip b { font-weight: 600; }
.spy-readout-tip .spy-readout-head { margin: 0 0 2px; color: var(--spy-readout-edited); font-weight: 600; }
.spy-readout-tip .spy-readout-title { margin: 0 0 2px; font-weight: 700; }
.spy-readout-tip .spy-readout-about { max-width: 280px; margin: 4px 0 0; color: var(--spy-readout-muted); white-space: normal; font-family: system-ui, sans-serif; }
.spy-readout-tip hr { margin: 4px 0; border: 0; border-top: 1px solid var(--spy-readout-border); }
`;

// ---------------------------------------------------------------------------------------------------
// thermodynamics for live values while editing (the server's redraw replaces them with SounderPy's)
// ---------------------------------------------------------------------------------------------------
const es = (T) => 6.112 * Math.exp((17.67 * T) / (T + 243.5)); // hPa, Bolton (1980)
const mixingRatio = (p, Td) => (622 * es(Td)) / (p - es(Td)); // g/kg
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
const fmt = (v, digits = 0, unit = "") => (v == null || !Number.isFinite(v) ? "–" : num(v, digits) + unit);
const windDir = (u, v) => ((Math.atan2(-u, -v) * 180) / Math.PI + 360) % 360;

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
   * @param {string[]} [options.fields] readout rows, from: p, z_agl, z, T, Td, Tw, Tv, rh, mr, theta, thetae, wind, omega
   * @param {"C"|"F"|"both"} [options.temperature="both"] temperature units
   * @param {"m"|"ft"} [options.height="m"] height units
   * @param {"auto"|"light"|"dark"} [options.theme="auto"] tooltip colors; auto follows the figure
   * @param {boolean} [options.editable=false] start with editing on
   * @param {boolean} [options.tapReadout=true] on touch screens, a tap shows the readout; false leaves taps to the
   *   page (e.g. a link around the image) unless editing
   * @param {boolean} [options.cursor=true] add the temperature (skew-t) or wind (hodograph) under the pointer
   * @param {string[]} [options.parts] what responds to the pointer: "skewt", "hodograph", "panels" (the 0-3 km panels),
   *   and "params" (values in the tables). Default is all four
   * @param {(level: object|null) => void} [options.onHover] called with the level under the pointer
   * @param {(profile: object) => void} [options.onEdit] called while a point is dragged
   * @param {(profile: object) => void} [options.onEditEnd] called when a drag finishes
   */
  constructor(img, readout, options = {}) {
    if (!(img instanceof HTMLImageElement)) throw new TypeError("SounderPyReadout needs an <img> element.");
    this.img = img;
    this.options = { fields: DEFAULT_FIELDS, temperature: "both", height: "m", theme: "auto", editable: false, tapReadout: true, cursor: true,
      parts: ["skewt", "hodograph", "panels", "params"], ...options };
    this.editing = false;
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
    img.after(this.layer, this.tip);
    this.ctx = this.layer.getContext("2d");

    this._listeners = [];
    const on = (target, type, fn, opts) => { target.addEventListener(type, fn, opts); this._listeners.push([target, type, fn, opts]); };
    on(this.layer, "pointermove", (e) => this._onMove(e));
    on(this.layer, "pointerdown", (e) => this._onDown(e));
    on(this.layer, "pointerup", (e) => this._onUp(e));
    on(this.layer, "pointercancel", (e) => this._onUp(e));
    on(this.layer, "pointerleave", (e) => { if (!this.drag && e.pointerType === "mouse") this._setHover(null); });
    // while editing, clicks belong to the editor (not, say, a link or button around the image)
    on(this.layer, "click", (e) => { if (this.editing) { e.preventDefault(); e.stopPropagation(); } });
    on(img, "load", () => this._layout());
    this._resize = new ResizeObserver(() => this._layout());
    this._resize.observe(img);
    // an image zoomed or panned with a CSS transform (e.g. in a zoomable viewer) takes the layer with it
    this._restyle = new MutationObserver(() => this._followTransform());
    this._restyle.observe(img, { attributes: true, attributeFilter: ["style", "class"] });

    this.update(readout);
    if (this.options.editable) this.setEditing(true);
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
    this._layout();
  }

  /** Turn profile editing on or off. */
  setEditing(on) {
    this.editing = Boolean(on);
    this.layer.toggleAttribute("data-editing", this.editing);
    this._draw();
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
    this.options.onEditEnd?.(this.getProfile());
  }

  /** Put the profile back to the figure's. */
  resetEdits() {
    this.edit = structuredClone(this.original);
    this.history = [];
    this._draw();
    this.options.onEditEnd?.(this.getProfile());
  }

  /** The readout values at a pressure (hPa), or null outside the profile. */
  levelAt(p) {
    return this._level(Math.log10(p));
  }

  /** Show or hide the layer, e.g. for a readout mode that can be switched off. */
  setEnabled(on) {
    this.layer.hidden = !on;
    if (!on) this._setHover(null);
  }

  /** Remove the layer and listeners. */
  destroy() {
    this._resize.disconnect();
    this._restyle.disconnect();
    cancelAnimationFrame(this._frame);
    for (const [target, type, fn, opts] of this._listeners) target.removeEventListener(type, fn, opts);
    this.layer.remove();
    this.tip.remove();
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
    return [((e.clientX - rect.left) / rect.width) * this.ro.image.width, ((e.clientY - rect.top) / rect.height) * this.ro.image.height];
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
    this._client = [e.clientX, e.clientY];
    const [x, y] = this._eventPx(e);
    if (this.drag) { this._dragTo(x, y); return; }
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
    this._client = [e.clientX, e.clientY];
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
    // a tap shows the readout on touch screens
    if (e.pointerType !== "mouse" && (this.options.tapReadout || this.editing)) this._setHover(this._hoverAt(x, y));
  }

  _onUp(e) {
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
    this.options.onEdit?.(this.getProfile());
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

    if (this.editing || this.isEdited()) this._drawEdits(ctx);

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
    if (param.kind === "wind") {
      const [u, w] = v;
      return `${fmt(windDir(u, w), 0, "°")} @ ${fmt(Math.hypot(u, w), 0, " kt")}`;
    }
    const a = Math.abs(v);
    const digits = param.units === "/s" ? 4 : a >= 100 ? 0 : a >= 10 ? 1 : 2;
    return fmt(v, digits, param.units ? ` ${param.units}` : "");
  }

  // a parameter's layer, shaded on the skew-t, hodograph, and height panels
  _drawLayer(ctx, layer) {
    // several single heights, such as every freezing level
    if (Array.isArray(layer.z_at)) { for (const z_at of layer.z_at) this._drawLayer(ctx, { z_at }); return; }
    const S = this.ro.panels.skewt, H = this.ro.panels.hodograph;
    const band = this.ro.dark ? "rgba(255, 214, 102, 0.22)" : "rgba(255, 170, 0, 0.18)";
    let zs = null, lps = null;
    if (layer.z) { zs = layer.z; lps = zs.map((z) => this._lpAtHeight(z)); }
    else if (layer.p) { lps = layer.p.map((p) => Math.log10(p)); zs = lps.map((lp) => this._zAtLp(lp)); }
    else if (layer.z_at != null) { zs = [layer.z_at, layer.z_at]; lps = zs.map((z) => this._lpAtHeight(z)); }
    if (!zs) return;
    const line = zs[0] === zs[1];
    const shade = (bbox, y0, y1) => this._clip(ctx, bbox, () => {
      ctx.fillStyle = band; ctx.strokeStyle = "#e08a00"; ctx.lineWidth = 3;
      if (line) { ctx.beginPath(); ctx.moveTo(bbox[0], y0); ctx.lineTo(bbox[2], y0); ctx.stroke(); }
      else ctx.fillRect(bbox[0], Math.min(y0, y1), bbox[2] - bbox[0], Math.abs(y1 - y0));
    });
    if (S && lps[0] != null && lps[1] != null) shade(S.bbox, ...lps.map((lp) => this._fwd(S.transform, 0, lp)[1]));
    const tr = this.ro.hodograph_trace;
    if (H && tr && zs[0] != null && zs[1] != null) {
      const pts = [];
      for (let i = 0; i < tr.z_agl.length; i++) if (tr.z_agl[i] >= zs[0] - 1 && tr.z_agl[i] <= zs[1] + 1) pts.push(this._hodoPx(tr.u[i], tr.v[i]));
      this._clip(ctx, H.bbox, () => {
        if (pts.length > 1) { ctx.globalAlpha = 0.4; this._line(ctx, pts, "#e08a00", [], 16); ctx.globalAlpha = 1; }
        else if (pts.length === 1) this._dot(ctx, pts[0], "#e08a00", 8);
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
    const where = L?.z ? span(L.z[0], L.z[1], "m AGL") : L?.p ? span(L.p[0], L.p[1], "hPa") : "";
    const whereLabel = (L?.z && Math.round(L.z[0]) === Math.round(L.z[1])) || (L?.p && Math.round(L.p[0]) === Math.round(L.p[1])) ? "At" : "Layer";
    this.tip.innerHTML = stale + `<p class="spy-readout-title">${esc(param.name)}</p>`
      + `<div><span>Value</span><b>${esc(this._paramValue(param))}</b></div>`
      + (where ? `<div><span>${whereLabel}</span><b>${where}</b></div>` : "")
      + (param.extra || []).map(([k, v]) => `<div><span>${esc(k)}</span><b>${esc(v)}</b></div>`).join("")
      + (param.about ? `<p class="spy-readout-about">${esc(param.about)}</p>` : "");
    this._placeTip();
  }

  _showHeightOnlyTip() {
    this.tip.innerHTML = this._heightRows();
    this._placeTip();
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

  _height(m, suffix = "") {
    if (m == null || !Number.isFinite(m)) return "–";
    return this.options.height === "ft" ? `${num(m * 3.28084)} ft${suffix}` : `${num(m)} m${suffix}`;
  }

  _showTip(l) {
    const rows = [];
    for (const key of this.options.fields) {
      let value;
      if (key === "p") value = fmt(l.p, 0, " hPa");
      else if (key === "z_agl") value = this._height(l.z_agl, " AGL");
      else if (key === "z") value = this._height(l.z, " MSL");
      else if (["T", "Td", "Tw", "Tv"].includes(key)) value = this._temp(l[key]);
      else if (key === "rh") value = fmt(l.rh, 0, "%");
      else if (key === "mr") value = fmt(l.mr, 1, " g/kg");
      else if (key === "theta" || key === "thetae") value = fmt(l[key], 1, " K");
      else if (key === "wind") value = l.u == null ? "–" : `${fmt(l.wdir, 0, "°")} @ ${fmt(l.wspd, 0, " kt")}`;
      else if (key === "omega") { if (l.omega == null) continue; value = fmt(l.omega, 2, " Pa/s"); }
      else continue;
      rows.push([LABELS[key] ?? key, value]);
    }
    if (this.options.cursor && this.hover.cursorT != null) rows.push(["Cursor T", this._temp(this.hover.cursorT)]);
    if (this.options.cursor && this.hover.cursorUV) {
      const [u, v] = this.hover.cursorUV;
      rows.push(["Cursor", `${fmt(windDir(u, v), 0, "°")} @ ${fmt(Math.hypot(u, v), 0, " kt")}`]);
    }
    const head = l.edited ? `<p class="spy-readout-head">Edited profile</p>` : "";
    this.tip.innerHTML = head + this._heightRows() + rows.map(([k, v]) => `<div><span>${k}</span><b>${v}</b></div>`).join("");
    this._placeTip();
  }

  _placeTip() {
    const tip = this.tip;
    tip.style.display = "block";
    // beside the pointer, inside the part of the image's container that's on screen, with a margin
    const parent = tip.offsetParent || this.img.parentElement;
    const box = parent.getBoundingClientRect(), margin = 8, gap = 16;
    const area = {
      left: Math.max(box.left, 0) + margin, top: Math.max(box.top, 0) + margin,
      right: Math.min(box.right, window.innerWidth) - margin, bottom: Math.min(box.bottom, window.innerHeight) - margin,
    };
    const [cx, cy] = this._client ?? [area.left, area.top];
    const w = tip.offsetWidth, h = tip.offsetHeight;
    // right of and below the pointer, or flipped to the left and above when it won't fit
    let x = cx + gap, y = cy + gap;
    if (x + w > area.right) x = cx - gap - w;
    if (x < area.left) x = Math.min(Math.max(area.left, cx - w / 2), area.right - w);
    if (y + h > area.bottom) y = cy - gap - h;
    if (y < area.top) y = Math.min(Math.max(area.top, area.bottom - h), cy + gap);
    tip.style.left = `${Math.max(area.left, x) - box.left - parent.clientLeft}px`;
    tip.style.top = `${Math.max(area.top, y) - box.top - parent.clientTop}px`;
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

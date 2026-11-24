// Tests for the browser script, sounderpy-readout.js, run with Node's test runner (node --test).
//
// tests/test_readout_js.py runs this: it draws a sounding with readout=True and passes the script (as an .mjs copy) and
// the readout file in READOUT_JS and READOUT_JSON. The checks compare the script's own math with what SounderPy computed
// for the figure, and need no browser: the readout's calculations don't touch the page.

import { test } from "node:test";
import assert from "node:assert/strict";
import { close, load, view as viewOf } from "./readout_view.mjs";

const READOUT = load(process.env.READOUT_JSON);
const FIGURE = Object.fromEntries(READOUT.params.map((p) => [p.key, p.value]));
const view = (options) => viewOf(READOUT, options);

const number = (text) => Number(String(text).replace(/[^0-9.+-]/g, ""));
const row = (rows, label) => rows.find(([k]) => k === label)?.[1];

test("the skew-t and hodograph pixel mappings invert exactly", () => {
  const v = view(), S = READOUT.panels.skewt;
  for (let i = 0; i < v.n; i++) {
    const [x, y] = v._skewPx(READOUT.profile.T[i], i);
    const [T, lp] = v._inv(S.transform, x, y);
    close(T, READOUT.profile.T[i], 1e-6, `T at level ${i}`);
    close(lp, v.lp[i], 1e-9, `log p at level ${i}`);
    const [u, w] = v._hodoUV(...v._hodoPx(READOUT.profile.u[i], READOUT.profile.v[i]));
    close(u, READOUT.profile.u[i], 1e-6, `u at level ${i}`);
    close(w, READOUT.profile.v[i], 1e-6, `v at level ${i}`);
  }
});

test("the readout at a level is the profile at that level", () => {
  const v = view(), P = READOUT.profile;
  for (let i = 0; i < v.n; i++) {
    const l = v._level(v.lp[i]);
    close(l.T, P.T[i], 1e-9, `T at ${P.p[i]} hPa`);
    close(l.Td, P.Td[i], 1e-9, `Td at ${P.p[i]} hPa`);
    close(l.z_agl, P.z_agl[i], 1e-6, `height at ${P.p[i]} hPa`);
    close(l.wspd, P.wspd[i], 0.05, `wind speed at ${P.p[i]} hPa`);
    if (P.wspd[i] > 1) close(Math.abs(((l.wdir - P.wdir[i] + 540) % 360) - 180), 0, 0.5, `wind direction at ${P.p[i]} hPa`);
  }
});

test("between levels, values are linear in log pressure", () => {
  const v = view(), P = READOUT.profile;
  for (let i = 0; i < v.n - 1; i++) {
    const mid = v._level((v.lp[i] + v.lp[i + 1]) / 2);
    close(mid.T, (P.T[i] + P.T[i + 1]) / 2, 1e-9, `T between levels ${i} and ${i + 1}`);
  }
});

test("bulk shear and SRH match the figure's", () => {
  const v = view();
  for (const [z0, z1, key] of [[0, 500, "0_to_500"], [0, 1000, "0_to_1000"], [1000, 3000, "1_to_3000"],
                               [3000, 6000, "3_to_6000"], [0, 3000, "0_to_3000"], [0, 6000, "0_to_6000"]]) {
    const k = v._kinematics(z0, z1);
    const srh = FIGURE[`kinem.srh_${key}`], shear = FIGURE[`kinem.shear_${key}`];
    if (srh != null) close(k.srh, srh, Math.max(2, Math.abs(srh) * 0.03), `${key} SRH`);
    if (shear != null) close(k.shear, shear, 0.6, `${key} bulk shear`);
  }
});

test("storm-relative wind matches the figure's", () => {
  const v = view();
  for (const [z0, z1, key] of [[0, 500, "0_to_500"], [0, 1000, "0_to_1000"], [1000, 3000, "1_to_3000"], [3000, 6000, "3_to_6000"]]) {
    const expected = FIGURE[`kinem.srw_${key}`];
    if (expected != null) close(v._kinematics(z0, z1).srw, expected, 1, `${key} SR wind`);
  }
});

test("the storm motion tool, at the figure's storm motion, gives the figure's SRH", () => {
  const v = view(), [u, w] = READOUT.storm_motion;
  const rows = v._markResults({ kind: "storm", u, v: w }).rows;
  close(number(row(rows, "0–1 km SRH")), FIGURE["kinem.srh_0_to_1000"], 2, "0-1 km SRH");
  close(number(row(rows, "0–3 km SRH")), FIGURE["kinem.srh_0_to_3000"], Math.max(2, Math.abs(FIGURE["kinem.srh_0_to_3000"]) * 0.03), "0-3 km SRH");
});

test("storm motions snap to whole degrees and knots, and read by where the storm moves", () => {
  const v = view();
  const [u, w] = v._snapMotion(12.34, -7.89);
  const from = (Math.atan2(-u, -w) * 180) / Math.PI;
  close(Math.hypot(u, w), Math.round(Math.hypot(12.34, -7.89)), 1e-9, "speed");
  close(Math.abs(from - Math.round(from)), 0, 1e-9, "direction is a whole degree");
  // a storm moving east at 20 kt comes from 270° and moves toward E (90°)
  const res = v._markResults({ kind: "storm", u: 20, v: 0 });
  close(res.direction, 270, 1e-9, "direction it comes from");
  assert.match(row(res.rows, "Moving toward"), /^E \(90°\) @ 20 kt$/);
});

test("SB, ML, and MU parcels are close to SounderPy's", () => {
  const v = view();
  for (const [preset, tolerance] of [["sb", 0.03], ["ml", 0.08], ["mu", 0.03]]) {
    const r = v._parcelFor({ kind: "parcel", preset, i: 0 });
    const cape = FIGURE[`thermo.${preset}cape`], lcl = FIGURE[`thermo.${preset}_lcl_z`];
    close(r.cape, cape, Math.max(25, cape * tolerance), `${preset.toUpperCase()}CAPE`);
    close(v._zAtLp(r.lcl.lp), lcl, 25, `${preset.toUpperCase()} LCL height`);
  }
});

test("CAPE up to the EL is the full CAPE, and there's none below the LFC", () => {
  const v = view(), r = v._parcelFor({ kind: "parcel", preset: "sb", i: 0 });
  assert.ok(r.cape > 0 && r.lfc != null && r.el != null, "the test profile has a buoyant surface parcel");
  close(v._parcelTo(r, r.el).cape, r.cape, 0.5, "CAPE to the EL");
  close(v._parcelTo(r, r.lfc + 0.01).cape, 0, 1e-9, "CAPE below the LFC");
});

test("the layer tool's 0-3 km lapse rate (virtual temperature) and SRH match the figure's", () => {
  const v = view(), rows = v._markResults({ kind: "layer", z0: 0, z1: 3000 }).rows;
  close(number(row(rows, "Lapse rate")), FIGURE["thermo.lr_03km"], 0.1, "0-3 km lapse rate");
  close(number(row(rows, "SRH")), FIGURE["kinem.srh_0_to_3000"], Math.max(3, Math.abs(FIGURE["kinem.srh_0_to_3000"]) * 0.03), "0-3 km SRH");
});

test("edits change the profile that's sent back, and the readout follows them", () => {
  const v = view();
  assert.equal(v.isEdited(), false);
  v.edit.T[0] += 5.004;
  v.edit.Td[0] = v.edit.T[0];
  assert.equal(v.isEdited(), true);
  const sent = v.getProfile();
  assert.equal(sent.T[0], Math.round((READOUT.profile.T[0] + 5.004) * 100) / 100);
  assert.equal(sent.p.length, v.n);
  const l = v._level(v.lp[0]);
  close(l.rh, 100, 1e-9, "RH of a saturated edited level");
});

test("units: speeds, heights, and temperatures", () => {
  assert.equal(view({ wind: "m/s" })._speed(10), "5 m/s");
  assert.equal(view({ wind: "mph" })._speed(10), "12 mph");
  assert.equal(view({ height: "ft" })._height(1000), "3281 ft");
  assert.equal(view({ temperature: "F" })._temp(0), "32 °F");
  assert.equal(view({ temperature: "both" })._temp(10), "10.0 °C (50 °F)");
  // no "-0.0" for a value that rounds to zero
  assert.equal(view({ temperature: "C" })._temp(-0.01), "0.0 °C");
});

test("the layer tool's ends snap to standard levels, so its depth is a thickness", () => {
  const v = view(), S = READOUT.panels.skewt;
  v.img = { getBoundingClientRect: () => ({ width: READOUT.image.width }) }; // shown at full size
  const y850 = v._fwd(S.transform, 0, Math.log10(850))[1];
  const z850 = v._zAtLp(Math.log10(850));
  close(v._zFromSkew(0, y850 + 4), z850, 1e-9, "4 px from 850 hPa snaps to it");
  assert.ok(Math.abs(v._zFromSkew(0, y850 + 30) - z850) > 1, "30 px away doesn't snap");
  const rows = v._markResults({ kind: "layer", z0: 0, z1: 3000 }).rows;
  assert.match(row(rows, "Depth"), /^3000 m \(300 dam\)$/);
  // on a small figure, a point between two levels snaps to the nearer one
  v.img = { getBoundingClientRect: () => ({ width: READOUT.image.width * 0.3 }) };
  const y925 = v._fwd(S.transform, 0, Math.log10(925))[1];
  assert.equal(v._snapP(y850 + 0.3 * (y925 - y850)), 850);
  assert.equal(v._snapP(y850 + 0.7 * (y925 - y850)), 925);
  // both ends on standard levels: the layer's thickness
  const z500 = v._zAtLp(Math.log10(500));
  const thick = row(v._markResults({ kind: "layer", z0: z850, z1: z500, snap0: 850, snap1: 500 }).rows, "Thickness");
  assert.equal(thick, `850–500 hPa: ${Math.round((z500 - z850) / 10)} dam`);
});

test("arrow keys step the readout, and an unpinned parcel, one level at a time", () => {
  const v = view();
  v.layer = { hidden: false };
  v._draw = () => {};
  const key = (k) => { let stopped = false; v._onKey({ key: k, target: {}, preventDefault() {}, stopPropagation() { stopped = true; } }); return stopped; };
  v.hover = { panel: "skewt", lp: v.lp[3] };
  assert.equal(key("ArrowUp"), true, "the readout handles the key");
  assert.equal(v.hover.lp, v.lp[4]);
  key("ArrowDown"); key("ArrowDown");
  assert.equal(v.hover.lp, v.lp[2]);
  v.hover = { panel: "skewt", lp: v.lp[0] };
  assert.equal(key("ArrowDown"), false, "at the surface, the page keeps the key");
  v.mark = { kind: "parcel", i: 2 };
  key("ArrowUp");
  assert.equal(v.mark.i, 3);
  v.mark.picked = true;
  v.hover = null;
  assert.equal(key("ArrowUp"), false, "a pinned parcel stays put");
  v.options.keys = false;
  v.mark = null; v.hover = { panel: "skewt", lp: v.lp[3] };
  assert.equal(key("ArrowUp"), false, "keys: false leaves arrow keys to the page");
});

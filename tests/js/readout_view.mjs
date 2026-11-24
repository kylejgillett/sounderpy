// Shared by the browser script's tests: the script, a readout file, and a readout on no page.

import { readFileSync } from "node:fs";
import { pathToFileURL } from "node:url";

export const { SounderPyReadout } = await import(pathToFileURL(process.env.READOUT_JS).href);

export const load = (path) => JSON.parse(readFileSync(path, "utf8"));

// a readout on no page: the state update() sets, without the layer and tooltip elements
export function view(readout, options = {}) {
  const v = Object.create(SounderPyReadout.prototype);
  const P = readout.profile;
  v.ro = readout;
  v.n = P.p.length;
  v.lp = P.p.map((p) => Math.log10(p));
  v.original = { T: [...P.T], Td: [...P.Td], u: [...P.u], v: [...P.v] };
  v.edit = structuredClone(v.original);
  v.history = [];
  v.hover = null;
  v.mark = null;
  v.drawings = [];
  v._touches = new Set();
  v.options = { fields: ["p", "z_agl", "T", "Td", "wind"], temperature: "C", height: "m", wind: "kt", cursor: true, keys: true, ...options };
  return v;
}

export const close = (actual, expected, tolerance, what) => {
  if (!(Math.abs(actual - expected) <= tolerance)) throw new Error(`${what}: ${actual} vs ${expected} (± ${tolerance})`);
};

// The browser script's precipitation type and cloud layers, against SounderPy's winter calculations.
//
// tests/test_readout_js.py runs this once per idealized winter case: READOUT_JSON is a winter figure's readout file and
// READOUT_EXPECTED what calc_winter computed for it (the precipitation source, surface type, melting and refreezing
// layers, and cloud layers).

import { test } from "node:test";
import assert from "node:assert/strict";
import { close, load, view } from "./readout_view.mjs";

const READOUT = load(process.env.READOUT_JSON);
const EXPECTED = load(process.env.READOUT_EXPECTED);

test("precipitation from SounderPy's source reaches the surface as the same type", () => {
  const r = view(READOUT)._ptype(Math.log10(EXPECTED.source_p));
  assert.equal(r.surface, EXPECTED.ptype);
});

test("the melting and refreezing layers and their energies match SounderPy's", () => {
  const r = view(READOUT)._ptype(Math.log10(EXPECTED.source_p));
  assert.deepEqual(r.layers.map((L) => L.kind), EXPECTED.layers.map((L) => L.kind));
  r.layers.forEach((L, k) => {
    const want = EXPECTED.layers[k];
    close(L.energy, want.energy, Math.max(2, Math.abs(want.energy) * 0.1), `layer ${k} energy`);
    close(L.bottom, want.bottom, 60, `layer ${k} base`);
    close(L.top, want.top, 60, `layer ${k} top`);
  });
});

test("cloud layers match SounderPy's (below 8 km, where it looks for them)", () => {
  const clouds = view(READOUT)._clouds().filter(([a]) => a <= 8000);
  assert.equal(clouds.length, EXPECTED.clouds.length);
  clouds.forEach(([a, b], k) => {
    close(a, EXPECTED.clouds[k].bottom, 150, `cloud ${k} base`);
    close(Math.min(b, 8000), Math.min(EXPECTED.clouds[k].top, 8000), 150, `cloud ${k} top`);
  });
});

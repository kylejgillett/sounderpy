"""The browser script, sounderpy-readout.js, tested against a readout file SounderPy just made.

Draws the synthetic sounding with ``readout=True`` and runs ``tests/js/test_readout.mjs`` with Node's test runner, which
checks the script's math (pixel mappings, interpolation, shear, SRH, storm-relative wind, parcels, layer stats, edits,
units, snapping, and arrow keys) against the figure's values. Then, for each idealized winter case,
``tests/js/test_readout_winter.mjs`` checks the script's precipitation type and cloud layers against SounderPy's winter
calculations. Skipped without Node.js 18 or newer.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

import sounderpy as spy
from sounderpy.readout import readout_filename

from helpers import WINTER_CASES, make_winter_profile


NODE = shutil.which("node")
JS = Path(__file__).parent / "js"


def node_version():
    if NODE is None:
        return None
    match = re.match(r"v(\d+)", subprocess.run([NODE, "--version"], capture_output=True, text=True).stdout)
    return int(match.group(1)) if match else None


needs_node = pytest.mark.skipif((node_version() or 0) < 18, reason="needs Node.js 18 or newer for the browser script's tests")


def run_node(tests, tmp_path, readout, **env):
    # an .mjs copy, so Node loads it as a module whatever its version
    script = tmp_path / "sounderpy-readout.mjs"
    shutil.copy(spy.readout_js_path(), script)
    env = {**os.environ, "READOUT_JS": str(script), "READOUT_JSON": str(readout), **env}
    result = subprocess.run([NODE, "--test", str(JS / tests)], env=env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.slow
@needs_node
def test_readout_script(clean_profile, tmp_path):
    figure = tmp_path / "sounding.png"
    spy.build_sounding(clean_profile, map_zoom=0, dpi=50, special_parcels="simple", save=True, filename=str(figure),
                       readout=True)
    run_node("test_readout.mjs", tmp_path, readout_filename(figure))


@pytest.mark.slow
@needs_node
@pytest.mark.parametrize("case", list(WINTER_CASES))
def test_readout_script_winter(case, winter_params_for, tmp_path):
    figure = tmp_path / f"winter_{case}.png"
    spy.build_sounding(make_winter_profile(case), mode="winter", map_zoom=0, dpi=50, save=True, filename=str(figure),
                       readout=True)
    winter = winter_params_for(case)["winter"]
    expected = {
        "source_p": float(winter["source_p"]),
        "ptype": winter["ptype"],
        "layers": [{key: (str(lyr[key]) if key == "kind" else float(lyr[key])) for key in ("kind", "bottom", "top", "energy")}
                   for lyr in winter["layers"]],
        "clouds": [{"bottom": float(c["bottom"]), "top": float(c["top"])} for c in winter["clouds"]],
    }
    (tmp_path / "expected.json").write_text(json.dumps(expected))
    run_node("test_readout_winter.mjs", tmp_path, readout_filename(figure), READOUT_EXPECTED=str(tmp_path / "expected.json"))

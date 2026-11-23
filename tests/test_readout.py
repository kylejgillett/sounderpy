"""Web readouts: the readout file saved with a figure, its pixel mapping, and editing a profile from the browser.

The renders check the mapping against the image itself: a point mapped from the profile's temperature lands on the
red temperature trace in the saved PNG.
"""

from __future__ import annotations

import copy
import json
import math
import warnings

import numpy as np
import pytest
from PIL import Image

import sounderpy as spy
from sounderpy import cli
from sounderpy.readout import READOUT_FORMAT, READOUT_VERSION, readout_filename

from helpers import make_winter_profile


def save_with_readout(function, data, path, **kwargs):
    kwargs.setdefault("map_zoom", 0)
    kwargs.setdefault("dpi", 50)
    function(data, save=True, filename=str(path), readout=True, **kwargs)
    readout = json.loads(readout_filename(path).read_text())
    image = Image.open(path).convert("RGB")
    return readout, np.asarray(image, dtype=int)


def to_pixel(matrix, a, b):
    return (matrix[0][0] * a + matrix[0][1] * b + matrix[0][2], matrix[1][0] * a + matrix[1][1] * b + matrix[1][2])


def near_color(pixels, x, y, test, radius=3):
    """Whether any pixel within ``radius`` of (x, y) passes ``test``."""
    x, y = int(round(x)), int(round(y))
    patch = pixels[max(y - radius, 0):y + radius + 1, max(x - radius, 0):x + radius + 1].reshape(-1, 3)
    return any(test(*rgb) for rgb in patch)


def is_red(r, g, b):
    return r > 180 and g < 90 and b < 90


def check_readout_file(readout, pixels, kind):
    assert readout["format"] == READOUT_FORMAT
    assert readout["version"] == READOUT_VERSION
    assert readout["figure"] == kind
    assert (readout["image"]["height"], readout["image"]["width"]) == pixels.shape[:2]
    profile = readout["profile"]
    assert len({len(values) for values in profile.values()}) == 1
    for name, panel in readout["panels"].items():
        x0, y0, x1, y1 = panel["bbox"]
        assert 0 <= x0 < x1 <= pixels.shape[1] and 0 <= y0 < y1 <= pixels.shape[0], name


def check_temperature_trace(readout, pixels):
    """Temperatures mapped through the readout land on the red trace, at levels inside the skew-t."""
    skewt = readout["panels"]["skewt"]
    x0, y0, x1, y1 = skewt["bbox"]
    profile = readout["profile"]
    hits = checked = 0
    for p, T in zip(profile["p"], profile["T"]):
        x, y = to_pixel(skewt["transform"], T, math.log10(p))
        if x0 + 10 < x < x1 - 10 and y0 + 10 < y < y1 - 10:
            checked += 1
            hits += near_color(pixels, x, y, is_red)
    assert checked >= 5
    # a level can sit under a label or marker drawn on the trace
    assert hits >= 0.9 * checked


# -----------------------------------------------------------------------------
# READOUT FILES
# -----------------------------------------------------------------------------

@pytest.mark.slow
def test_sounding_readout_lines_up_with_the_image(clean_profile, tmp_path):
    readout, pixels = save_with_readout(spy.build_sounding, clean_profile, tmp_path / "sounding.png", special_parcels="simple")
    check_readout_file(readout, pixels, "sounding")
    check_temperature_trace(readout, pixels)
    assert readout["dark"] is False
    assert readout["panels"]["hodograph"]["offset"] == [0.0, 0.0]
    assert readout["hodograph_trace"]["z_agl"][-1] <= 12000


@pytest.mark.slow
def test_zoomed_dark_storm_relative_readout(clean_profile, tmp_path):
    readout, pixels = save_with_readout(spy.build_sounding, clean_profile, tmp_path / "zoomed.png", special_parcels="simple",
                                        ylim=(1050, 500), dark_mode=True, sr_hodo=True, dpi=60)
    check_readout_file(readout, pixels, "sounding")
    check_temperature_trace(readout, pixels)
    assert readout["dark"] is True
    assert readout["panels"]["skewt"]["p_limits"][1] == pytest.approx(500)
    # the storm-relative hodograph is shifted by the storm motion
    assert readout["panels"]["hodograph"]["offset"] != [0.0, 0.0]


@pytest.mark.slow
def test_winter_readout(tmp_path):
    readout, pixels = save_with_readout(spy.build_sounding, make_winter_profile("fzra"), tmp_path / "winter.png", mode="winter")
    check_readout_file(readout, pixels, "winter")
    assert set(readout["panels"]) == {"skewt", "hodograph"}


@pytest.mark.slow
def test_hodograph_readout(clean_profile, tmp_path):
    readout, pixels = save_with_readout(spy.build_hodograph, clean_profile, tmp_path / "hodograph.png", sr_hodo=True)
    check_readout_file(readout, pixels, "hodograph")
    assert set(readout["panels"]) == {"hodograph"}
    # winds mapped through the readout land on the colored trace, not the background
    hodo = readout["panels"]["hodograph"]
    trace = readout["hodograph_trace"]
    background = pixels[int(hodo["bbox"][1]) + 5, int(hodo["bbox"][0]) + 5]
    off_background = lambda r, g, b: abs(r - background[0]) + abs(g - background[1]) + abs(b - background[2]) > 120
    points = [to_pixel(hodo["transform"], u - hodo["offset"][0], v - hodo["offset"][1]) for u, v in zip(trace["u"], trace["v"])]
    assert sum(near_color(pixels, x, y, off_background) for x, y in points) >= 0.9 * len(points)


def check_param_boxes(params, pixels):
    """Each tagged value's box holds its drawn text, and has a name and description from the catalog."""
    filled = 0
    for param in params:
        assert param["name"] != param["key"].split(".", 1)[1] and param["about"], param["key"]
        x0, y0, x1, y1 = (int(round(v)) for v in param["bbox"])
        patch = pixels[y0:y1 + 1, x0:x1 + 1].reshape(-1, 3)
        filled += (np.abs(patch - patch.mean(axis=0)).sum(axis=1) > 150).any()
    assert filled >= 0.95 * len(params)


def background_test(pixels, bbox):
    """A test for pixels that differ from the background just inside a panel's corner."""
    background = pixels[int(bbox[1]) + 3, int(bbox[0]) + 3]
    return lambda r, g, b: abs(r - background[0]) + abs(g - background[1]) + abs(b - background[2]) > 120


@pytest.mark.slow
def test_table_values_and_height_panels_line_up(clean_profile, tmp_path):
    readout, pixels = save_with_readout(spy.build_sounding, clean_profile, tmp_path / "tables.png", special_parcels="simple", dpi=80)

    params = {param["key"]: param for param in readout["params"]}
    for key in ("thermo.mlcape", "kinem.srh_0_to_3000", "general.pwat", "kinem.eil_stp", "kinem.mean_wind_pbl"):
        assert key in params
    assert params["kinem.srh_0_to_3000"]["layer"] == {"z": [0, 3000]}
    assert params["thermo.lr_75"]["layer"] == {"p": [700, 500]}
    assert params["kinem.srh_0_to_3000"]["value"] == pytest.approx(float(spy.sounding_params(clean_profile).calc()["kinem"]["srh_0_to_3000"]), rel=0.05)

    check_param_boxes(readout["params"], pixels)

    # each panel's first series, mapped through the readout, lands on its drawn line
    titles = [panel["title"] for panel in readout["height_panels"]]
    assert titles == ["Streamwiseness of ω", "Total & streamwise ω", "Storm-relative wind", "Stepwise CIN & CAPE"]
    panel = readout["height_panels"][2]
    drawn = background_test(pixels, panel["bbox"])
    points = [to_pixel(panel["transform"], v, z) for z, v in zip(panel["z"], panel["series"][0]["values"]) if v is not None and 50 < z < 2950]
    assert sum(near_color(pixels, x, y, drawn) for x, y in points) >= 0.9 * len(points)


@pytest.mark.slow
def test_winter_tables_and_height_panels(tmp_path):
    readout, pixels = save_with_readout(spy.build_sounding, make_winter_profile("fzra"), tmp_path / "winter_panels.png", mode="winter", dpi=80)
    assert [panel["title"] for panel in readout["height_panels"]] == ["T & Tw", "RHi & RH", "Wind", "θe"]

    check_param_boxes(readout["params"], pixels)
    by_key = {}
    for param in readout["params"]:
        by_key.setdefault(param["key"], []).append(param)
    # text values, winds, layers, single heights, and checklists
    assert by_key["winter.ptype"][0]["kind"] == "text" and by_key["winter.ptype"][0]["value"] == "FREEZING RAIN"
    assert by_key["winter.sfc_wind_uv"][0]["kind"] == "wind"
    melt = next(p for p in by_key["winter.layer.energy"] if p["value"] > 0)
    assert melt["name"] == "Melting layer: energy" and len(melt["layer"]["z"]) == 2
    assert len(by_key["winter.frz_pts_z"][0]["layer"]["z_at"]) == 2
    assert by_key["winter.source_p"][0]["layer"]["z_at"] > 0
    assert len(by_key["winter.fzfg_check"]) == 4 and all(p["extra"][-1][0] == "Check" for p in by_key["winter.fzfg_check"])
    # no CAPE, no LFC or EL: the figure says NA and so does the readout
    sb_el = by_key["thermo.sb_el_z"][0]
    assert sb_el["value"] is None and sb_el["layer"] is None


@pytest.mark.slow
def test_readout_needs_a_saved_figure(clean_profile):
    with pytest.warns(Warning, match="readout=True writes its file next to a saved figure"):
        fig, _ = spy.build_hodograph(clean_profile, map_zoom=0, dpi=50, readout=True, return_fig=True)


def test_readout_filename():
    assert readout_filename("plots/kokc.png").as_posix() == "plots/kokc.readout.json"
    assert readout_filename("kokc").as_posix() == "kokc.readout.json"
    assert readout_filename("kokc.v2.svg").as_posix() == "kokc.v2.readout.json"


def test_readout_script_is_installed():
    script = spy.readout_js_path()
    assert script.is_file()
    assert "export class SounderPyReadout" in script.read_text(encoding="utf-8")


# -----------------------------------------------------------------------------
# EDITED PROFILES
# -----------------------------------------------------------------------------

def browser_profile(data):
    """The profile as the readout script returns it from getProfile()."""
    as_list = lambda values, unit: [round(float(v), 2) for v in values.to(unit).m]
    return {"p": as_list(data["p"], "hPa"), "z": as_list(data["z"], "m"), "T": as_list(data["T"], "degC"),
            "Td": as_list(data["Td"], "degC"), "u": as_list(data["u"], "kt"), "v": as_list(data["v"], "kt")}


def test_apply_profile_edits(clean_profile):
    original = copy.deepcopy(clean_profile)
    edits = browser_profile(clean_profile)
    edits["T"][0] += 4
    edits["Td"][1] = edits["T"][1] + 5  # above the temperature, so it's capped
    edits["u"][2] += 10

    edited = spy.apply_profile_edits(clean_profile, edits)

    assert edited["T"][0].m == pytest.approx(edits["T"][0])
    assert edited["Td"][1].m == pytest.approx(edited["T"][1].m)
    assert edited["u"].to("kt").m[2] == pytest.approx(edits["u"][2])
    assert edited["site_info"]["edited"].startswith("T AT 1 LEVEL, TD AT 1 LEVEL, WIND AT 1 LEVEL (")
    assert edited["titles"] == clean_profile["titles"]
    # the original isn't changed
    np.testing.assert_array_equal(clean_profile["T"].m, original["T"].m)
    assert "edited" not in clean_profile["site_info"]


def test_unchanged_profile_is_not_marked_edited(clean_profile):
    edited = spy.apply_profile_edits(clean_profile, browser_profile(clean_profile))
    assert "edited" not in edited["site_info"]


def test_edits_on_thinned_levels_carry_other_values_over(clean_profile):
    # a readout of a thinned profile has fewer levels than the original
    data = copy.deepcopy(clean_profile)
    data["omega"] = np.linspace(-1, 1, len(data["p"])) * spy.readout.units("Pa/s")
    keep = slice(None, None, 2)
    thinned = {key: data[key][keep] for key in ("p", "z", "T", "Td", "u", "v")}
    edited = spy.apply_profile_edits(data, browser_profile(thinned))
    assert len(edited["p"]) == len(thinned["p"])
    assert len(edited["omega"]) == len(thinned["p"])
    np.testing.assert_allclose(edited["omega"].m, data["omega"].m[keep], atol=1e-6)


@pytest.mark.parametrize("change, message", [
    (lambda e: e.update(p=e["p"][::-1]), "pressures must be at least 3 values that decrease"),
    (lambda e: e.update(T=e["T"][:-1]), "'T' has"),
    (lambda e: e.pop("Td"), "needs 'Td'"),
    (lambda e: e.update(z=e["z"][::-1]), "heights must increase"),
])
def test_bad_edited_profiles_are_rejected(change, message, clean_profile):
    edits = browser_profile(clean_profile)
    change(edits)
    with pytest.raises(ValueError, match=message):
        spy.apply_profile_edits(clean_profile, edits)


@pytest.mark.slow
@pytest.mark.parametrize("function, kwargs", [(spy.build_sounding, {"special_parcels": "simple"}), (spy.build_hodograph, {})])
def test_edited_profiles_are_labeled(function, kwargs, clean_profile):
    edits = browser_profile(clean_profile)
    edits["T"][0] += 3
    edits["u"][0] += 5
    edited = spy.apply_profile_edits(clean_profile, edits)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig, axes = function(edited, map_zoom=0, dpi=50, return_fig=True, **kwargs)
    panel = axes["skewt" if "skewt" in axes else "hodograph"].ax
    assert any(text.get_text().startswith("EDITED PROFILE\n") for text in panel.texts)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def test_cli_readout_needs_a_plot_file(capsys):
    assert cli.main(["bufkit", "hrrr", "KOKC", "3", "--plot", "sounding", "--readout"]) != 0
    assert "--readout writes its file next to the saved plot" in capsys.readouterr().err

"""Winter figure components (``plot/components_winter.py``), drawn one at a time from the idealized winter cases."""

from __future__ import annotations

import copy
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pytest
from metpy.units import units

from sounderpy.calc_winter import winter_params
from sounderpy.plot import components as c
from sounderpy.plot import components_winter as w
from sounderpy.plot.style import HAZARDS_LAYOUT, PANEL_TOP, PTYPE_NAMES, get_winter_theme

from helpers import WINTER_CASES, make_winter_profile, texts_on


THEME = get_winter_theme()

OBS = {"station": "KSYN", "name": "SYNTHETIC", "distance_km": 3.2, "time": datetime(2026, 1, 1, 0, 0, tzinfo=timezone.utc),
       "offset_min": 0.0, "tmpf": 27.0, "dwpf": 25.0, "drct": 20.0, "sknt": 12.0, "gust": 22.0, "vsby": 1.5,
       "wxcodes": "-FZRA BR", "mslp": 1018.2, "skyc1": "BKN", "skyc2": "OVC", "skyc3": None}


def new_ax():
    return plt.figure().add_subplot()


def winter_skewt(case):
    data = make_winter_profile(case)
    return c._setup_skewt(plt.figure(), data["T"], THEME), data


# -----------------------------------------------------------------------------
# TEXT HELPERS
# -----------------------------------------------------------------------------

def test_value_formatting():
    assert w._fmt(12.0, " m") == "12 m"
    assert w._fmt(1.234, decimals=2) == "1.23"
    assert w._fmt(ma.masked, " m") == "NA"
    assert w._fmt_signed(6, "°C") == "+6.0°C"
    assert w._fmt_signed(-2.34, "°C") == "-2.3°C"
    assert w._fmt_signed(np.nan) == "NA" and w._fmt_signed(ma.masked) == "NA"


@pytest.mark.parametrize(("u", "v", "expected"), [(0.0, -10.0, "N @ 10 kt"), (10.0, 0.0, "W @ 10 kt"), (ma.masked, 5.0, "NA")])
def test_wind_text_names_where_the_wind_blows_from(u, v, expected):
    assert w._wind_text(u, v) == expected


def test_snow_ratio_text():
    assert w._slr_text({"slr": 12.4}) == "12:1"
    assert w._slr_text({"slr": ma.masked}) == "NA"


def test_pressure_at_height():
    z, p = np.array([0.0, 1000, 2000]), np.array([1000.0, 900, 800])
    assert w._p_at(z, p, 500) == pytest.approx(950)


def test_fit_width_shrinks_only_wide_text():
    ax = new_ax()
    short = w._fit_width(ax.text(0, 0, "x", fontsize=10), ax, 0.5)
    wide = w._fit_width(ax.text(0, 0, "a much wider line of text " * 4, fontsize=10), ax, 0.5)
    assert short.get_fontsize() == 10 and wide.get_fontsize() < 10


# -----------------------------------------------------------------------------
# SKEW-T
# -----------------------------------------------------------------------------

def test_melt_and_refreeze_layers_are_shaded_and_labeled(winter_params_for):
    skew, data = winter_skewt("fzra")
    winter = winter_params_for("fzra")["winter"]
    w._plot_melt_refreeze(skew, data, winter, THEME)
    labels = texts_on(skew.ax)
    # the warm nose and the surface cold layer under it; the cold layer above the warm nose is where snow forms
    assert "+85 J/kg" in labels and "-68 J/kg" in labels
    assert not any("-905" in label for label in labels)
    assert len(skew.ax.collections) == 2  # melt and refreeze fills


def test_no_layers_draws_nothing(winter_params_for):
    skew, data = winter_skewt("fzdz")
    winter = dict(winter_params_for("fzdz")["winter"], layers=[])
    w._plot_melt_refreeze(skew, data, winter, THEME)
    assert not skew.ax.texts and not skew.ax.collections


@pytest.mark.parametrize("case", WINTER_CASES)
def test_every_freezing_level_is_labeled(case, winter_params_for):
    skew, data = winter_skewt(case)
    winter = winter_params_for(case)["winter"]
    w._plot_freezing_levels(skew, data, winter, THEME)
    labels = texts_on(skew.ax)
    assert labels.count("←FRZ") == len(winter["frz_pts_z"])
    assert labels.count("←WB0") == len(winter["wb_frz_pts_z"])


@pytest.mark.parametrize(("case", "expected"), [("fzra", ["← SN", "← RA", "← FZRA"]), ("pl", ["← SN", "← RA", "← PL"]),
                                                ("fzdz", ["← DZ", "← FZDZ"])])
def test_phase_changes_are_labeled_on_the_way_down(case, expected, winter_params_for):
    skew, data = winter_skewt(case)
    w._plot_phase_annotations(skew, data, winter_params_for(case)["winter"], THEME)
    assert texts_on(skew.ax) == expected


def test_snow_label_carries_its_ratio(winter_params_for):
    skew, data = winter_skewt("rasn")
    winter = winter_params_for("rasn")["winter"]
    w._plot_phase_annotations(skew, data, winter, THEME)
    assert f"← SN {w._slr_text(winter)}" in texts_on(skew.ax)


def test_layer_bars_and_markers(winter_params_for):
    skew, data = winter_skewt("fzra")
    snd = winter_params_for("fzra")
    winter = snd["winter"]
    n_lines = len(skew.ax.lines)
    w._plot_dgz_bar(skew, winter, snd, THEME)
    w._plot_pbl(skew, data, snd, THEME)
    w._plot_clouds(skew, winter, THEME)
    w._plot_icing_layers(skew, winter, THEME)
    w._plot_warm_nose(skew, winter, THEME)
    labels = texts_on(skew.ax)
    assert {"DGZ", "PBL"} <= set(labels) and "ML" not in labels
    assert any(label.startswith("CLD TOP\n") for label in labels)
    assert len(skew.ax.lines) > n_lines
    assert "ICING LAYER" in [line.get_label() for line in skew.ax.lines]
    assert labels.count("ICING") == len(winter["icing_layers"])
    assert labels.count("CLOUD") == len(winter["clouds"])


def test_no_warm_nose_or_dgz_draws_nothing(winter_params_for):
    skew, data = winter_skewt("rasn")
    winter = dict(winter_params_for("rasn")["winter"], dgz_depth=ma.masked)
    n_lines = len(skew.ax.lines)
    w._plot_warm_nose(skew, winter, THEME)  # rasn has no warm nose
    w._plot_dgz_bar(skew, winter, None, THEME)
    assert len(skew.ax.lines) == n_lines and not skew.ax.texts


def test_rh_ice_labels_below_freezing(winter_params_for):
    skew, _ = winter_skewt("fzra")
    snd = winter_params_for("fzra")
    w._plot_rh_ice_annotations(skew, snd, snd["winter"], THEME)
    labels = texts_on(skew.ax)
    assert labels and all(label.startswith(("RHᵢ\n", "RH\n")) and label.endswith("% →") for label in labels)


# -----------------------------------------------------------------------------
# HODOGRAPH
# -----------------------------------------------------------------------------

def hodograph():
    return c._setup_hodograph(plt.figure(figsize=(10, 10)).add_axes((0.1, 0.1, 0.6, 0.6)), THEME, (-60, 60), (-60, 60))


def test_mean_wind_vectors(winter_params_for):
    hodo = hodograph()
    w._plot_mean_wind_vectors(hodo, winter_params_for("fzra")["winter"], THEME)
    assert [text.get_text() for text in hodo.ax.texts if text.get_text()] == ["PBL", "DGZ"]


def test_winds_box_lists_each_level(winter_params_for):
    ax = new_ax()
    winter = winter_params_for("fzra")["winter"]
    w._plot_winds_box(ax, winter, THEME)
    labels = texts_on(ax)
    assert {"WINDS", "10 m:", "PBL mean:", "DGZ mean:", "Max 0-2 km:"} <= set(labels)
    assert [f"{level} hPa:" for level in winter["level_winds"]] == [label for label in labels if label.endswith("hPa:")]


def test_station_model_draws_the_observation():
    ax = w._plot_station_obs(hodograph(), copy.deepcopy(OBS), THEME)
    labels = texts_on(ax)
    assert "Latest Observation\nKSYN (0000Z) · 3 km" in labels
    assert {"27", "25", "182", "G22", "1.5"} <= set(labels)  # T, Td, coded SLP (1018.2 hPa), gust, visibility


def test_station_model_with_missing_values():
    obs = dict(OBS, tmpf=None, dwpf=None, mslp=None, vsby=None, gust=None, wxcodes=None, sknt=None)
    labels = texts_on(w._plot_station_obs(hodograph(), obs, THEME))
    assert labels.count("NA") == 4
    assert not any(label.startswith("G") for label in labels)


# -----------------------------------------------------------------------------
# SIDE PANELS
# -----------------------------------------------------------------------------

def test_temperature_panel_fits_the_data(winter_params_for):
    ax = new_ax()
    winter = winter_params_for("fzra")["winter"]
    w._plot_temperature_panel(ax, winter, THEME)
    low, high = ax.get_xlim()
    sel = winter["z"] <= PANEL_TOP
    assert low <= min(np.nanmin(winter["T"][sel]), np.nanmin(winter["Tw"][sel])) and high >= np.nanmax(winter["T"][sel])
    assert ax.get_ylim() == (0, PANEL_TOP)


def test_wind_panel_has_at_most_two_ticks(winter_params_for):
    ax = new_ax()
    w._plot_wind_panel(ax, winter_params_for("pl")["winter"], THEME)
    assert 1 <= len(ax.get_xticks()) <= 2


@pytest.mark.parametrize("draw", [w._plot_temperature_panel, w._plot_rh_panel, w._plot_wind_panel, w._plot_thetae_panel])
@pytest.mark.parametrize("case", WINTER_CASES)
def test_side_panels_draw(draw, case, winter_params_for):
    ax = new_ax()
    draw(ax, winter_params_for(case)["winter"], THEME)
    assert ax.lines


@pytest.mark.parametrize(("depth", "drawn"), [(ma.masked, False), (0.0, False), (500.0, True), (150.0, True), (PANEL_TOP + 1.0, False)])
def test_pbl_top_mark(depth, drawn):
    ax = new_ax()
    w._mark_pbl_top(ax, {"pbl_depth": depth}, THEME)
    assert (texts_on(ax) == ["PBL"]) is drawn


# -----------------------------------------------------------------------------
# TABLES
# -----------------------------------------------------------------------------

@pytest.mark.parametrize(("case", "rows"), [("fzra", ["Cold aloft", "Melt", "Sfc refrz"]), ("pl", ["Cold aloft", "Melt", "Sfc refrz"]),
                                            ("rasn", ["Cold aloft", "Sfc melt"]), ("fzdz", ["Sfc cold"])])
def test_precip_table(case, rows, winter_params_for):
    ax = new_ax()
    snd = winter_params_for(case)
    w._plot_precip_table(ax, snd["winter"], snd, THEME)
    labels = texts_on(ax)
    assert PTYPE_NAMES.get(snd["winter"]["ptype"], snd["winter"]["ptype"]) in labels
    assert snd["winter"]["sharppy_ptype"].upper() in labels
    assert [label for label in labels if label in ("Cold aloft", "Melt", "Sfc melt", "Refreeze", "Sfc refrz", "Sfc cold")] == rows
    assert {"SB:", "MU:", "CAPE", "EL"} <= set(labels)


def test_precip_table_without_a_saturated_source():
    data = make_winter_profile("fzra")
    data["Td"] = data["T"] - 25 * units.delta_degC
    snd = winter_params(data).calc()
    ax = new_ax()
    w._plot_precip_table(ax, snd["winter"], snd, THEME)
    labels = texts_on(ax)
    assert "none saturated" in labels
    assert "no saturated source (no precip expected)" in labels


def test_snow_growth_table(winter_params_for):
    for case, slr in (("rasn", True), ("fzra", False)):
        ax = new_ax()
        snd = winter_params_for(case)
        w._plot_snow_growth_table(ax, snd["winter"], snd, THEME)
        labels = texts_on(ax)
        assert "SNOW GROWTH (DGZ)" in labels and "Kuchera SLR:" in labels
        assert (w._slr_text(snd["winter"]) in labels) and (("NA" != w._slr_text(snd["winter"])) is slr)


@pytest.mark.parametrize("case", WINTER_CASES)
def test_hazards_table(case, winter_params_for):
    ax = new_ax()
    winter = winter_params_for(case)["winter"]
    w._plot_hazards_table(ax, winter, THEME)
    labels = texts_on(ax)
    assert {"COLD", "BLOWING SNOW", "OTHER", "FREEZING FOG", "FREEZING DRIZZLE"} <= set(labels)
    assert winter["fzdz"]["level"] in labels and winter["fzfg"]["level"] in labels


def test_freezing_drizzle_checks_include_lift_only_with_omega(winter_params_for):
    with_omega = w._fzdz_lines(winter_params_for("fzdz")["winter"])
    data = make_winter_profile("fzdz")
    data.pop("omega")
    without_omega = w._fzdz_lines(winter_params(data).calc()["winter"])
    assert with_omega[0] == ("title", "FREEZING DRIZZLE", "HIGH")
    assert [line[2] for line in with_omega if line[0] == "check"][-1] == "Lift in cloud"
    assert sum(line[0] == "check" for line in with_omega) == 5
    assert sum(line[0] == "check" for line in without_omega) == 4


def test_freezing_fog_checks(winter_params_for):
    lines = w._fzfg_lines(winter_params_for("fzra")["winter"])
    assert lines[0] == ("title", "FREEZING FOG", "MOD")
    assert [line[2] for line in lines if line[0] == "check"] == ["Sfc T < 0°C", "Sfc RH ≥ 95%", "10 m wind ≤ 6 kt", "Sfc inversion"]


def test_even_column_spreads_lines_from_top_to_bottom():
    ax = new_ax()
    top, bottom = HAZARDS_LAYOUT["top"], HAZARDS_LAYOUT["bottom"]
    # 5 slots: title, row, gap, and a two-line note
    w._draw_even_column(ax, [("title", "A", "HIGH"), ("row", "b", "c", "k"), ("gap",), ("note", "n\nn", 2)], THEME, 0.1, 0.9)
    y = {text.get_text(): text.get_position()[1] for text in ax.texts}
    step = (top - bottom) / 4
    assert y["A"] == pytest.approx(top) and y["HIGH"] == pytest.approx(top)
    assert y["b"] == pytest.approx(top - step)
    assert y["n\nn"] == pytest.approx(top - 3 * step + 0.03)

"""Severe figure components (``plot/components.py``), tested one helper at a time on small figures."""

from __future__ import annotations

import copy

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pytest
from metpy.units import units

from sounderpy.plot import components as c
from sounderpy.plot.style import (MISSING_VALUE, WATERMARK, WATERMARK_DEFAULTS, WATERMARK_TEXT, _get_plot_theme,
                                  set_watermark)

from helpers import make_clean_profile, make_winter_profile, texts_on


THEME = _get_plot_theme()


def axis_position(t, p, temperature_limits, ylim):
    """Where a temperature at a pressure sits across a zoomed skew-t, as a fraction of its width."""
    left, right = temperature_limits
    return (t - left) / (right - left) + c._zoom_lean(p, ylim)


def zoomed_skewt(data, ylim):
    return c._setup_zoomed_skewt(plt.figure(figsize=(8, 8)), [data], THEME, ylim, rect=(0.1, 0.1, 0.8, 0.8))


@pytest.fixture
def watermark():
    """Restore the default watermark after a test changes it (it is shared by the whole process)."""
    yield WATERMARK
    set_watermark()


# -----------------------------------------------------------------------------
# DATA PREP: thinning to max_levels
# -----------------------------------------------------------------------------

def test_thin_sounding_keeps_surface_and_top(clean_profile):
    thinned = c._thin_sounding(clean_profile, 10)
    for key in ("p", "z", "T", "Td", "u", "v"):
        assert len(thinned[key]) == 10
        assert thinned[key][0] == clean_profile[key][0]
        assert thinned[key][-1] == clean_profile[key][-1]
    assert thinned["site_info"] is clean_profile["site_info"]


def test_thin_sounding_thins_every_level_array():
    data = make_winter_profile("fzra")
    thinned = c._thin_sounding(data, 12)
    assert len(thinned["omega"]) == len(thinned["p"]) == 12


@pytest.mark.parametrize("max_levels", [None, 73, 500])
def test_thin_sounding_leaves_short_profiles_alone(clean_profile, max_levels):
    assert c._thin_sounding(clean_profile, max_levels) is clean_profile


def test_prepare_sounding_data_copies_thins_and_modifies(clean_profile):
    original = copy.deepcopy(clean_profile)
    data = c._prepare_sounding_data(clean_profile, {"T": 35}, max_levels=20)
    assert len(data["p"]) == 20
    assert data["T"][0].m == pytest.approx(35)
    np.testing.assert_array_equal(clean_profile["T"].m, original["T"].m)


# -----------------------------------------------------------------------------
# PRESSURE RANGE (ptop / ylim)
# -----------------------------------------------------------------------------

@pytest.mark.parametrize(("ylim", "expected"), [(None, None), ((1000, 700), (1000, 700)), ([1000, 700], (1000, 700)),
                                                ((np.float64(950), np.int64(500)), (950, 500))])
def test_validate_pressure_limits_accepts(ylim, expected):
    assert c._validate_pressure_limits(ylim) == expected


@pytest.mark.parametrize("ylim", [(700, 1000), (1000, 1000), (1000, 0), (1000, -5), (1000,), (1000, 700, 500),
                                  "1000,700", (True, 700), (1000, "700"), 1000])
def test_validate_pressure_limits_rejects(ylim):
    with pytest.raises(ValueError, match="'ylim' must be"):
        c._validate_pressure_limits(ylim)


def test_pressure_limits_from_ptop():
    assert c._get_pressure_limits(None) is None
    assert c._get_pressure_limits(700) == (1050, 700)
    for ptop in (0, 1050, 1100, True, "700"):
        with pytest.raises(ValueError, match="'ptop' must be"):
            c._get_pressure_limits(ptop)


@pytest.mark.parametrize(("ylim", "step"), [((1050, 100), 100), ((1000, 600), 50), ((1000, 900), 25), ((1000, 950), 10)])
def test_pressure_ticks_suit_the_range(ylim, step):
    ticks = c._get_pressure_ticks(ylim)
    assert all(ylim[1] < tick <= ylim[0] for tick in ticks)
    assert set(np.diff(ticks)) == {step}


def test_temperature_ticks_stay_in_view_and_widen_with_span():
    for limits in ((-30, 10), (-10, 16), (-60, 20)):
        ticks = c._get_temperature_ticks(limits)
        assert len(ticks) >= 3
        assert np.all((ticks - limits[0]) / (limits[1] - limits[0]) > 0.07)  # the left edge is for pressure labels
        assert np.all(ticks <= limits[1])
        assert len(set(np.diff(ticks))) == 1
    assert np.diff(c._get_temperature_ticks((-30, 10)))[0] < np.diff(c._get_temperature_ticks((-60, 60)))[0]


@pytest.mark.parametrize(("case", "ylim"), [("clean", (1000, 700)), ("clean", (1000, 300)), ("fzra", (1000, 700)),
                                            ("pl", (1000, 600)), ("rasn", (1000, 850))])
def test_zoomed_temperature_limits_center_the_lowest_3_km(case, ylim):
    data = make_clean_profile() if case == "clean" else make_winter_profile(case)
    limits = c._get_zoomed_temperature_limits([data], ylim)
    assert limits[1] - limits[0] <= c.ZOOM_MAX_SPAN

    # the lowest ~3 km of the levels in view
    p = data["p"].m
    in_view = (p <= ylim[0]) & (p >= ylim[1])
    low = in_view & (p >= p[in_view].max() * np.exp(-3000 / 8000))
    positions = np.concatenate([axis_position(data[key].m[low], p[low], limits, ylim) for key in ("T", "Td")])
    assert np.all((positions >= 0) & (positions <= 1))
    assert (positions.min() + positions.max()) / 2 == pytest.approx(0.5, abs=0.05)


def test_zoomed_sleet_profile_does_not_raise():
    # regression: the sleet profile warms with height in (1000, 600); no tick step once fit, raising StopIteration
    limits = c._get_zoomed_temperature_limits([make_winter_profile("pl")], (1000, 600))
    assert np.isfinite(limits).all()
    assert len(c._get_temperature_ticks(limits)) > 0
    zoomed_skewt(make_winter_profile("pl"), (1000, 600))


def test_zoomed_temperature_limits_without_data_in_view():
    # a profile entirely below the view falls back to the full skew-t's range
    data = {"p": np.array([1000.0, 950, 900]) * units.hPa, "T": np.array([10.0, 8, 6]) * units.degC,
            "Td": np.array([5.0, 4, 3]) * units.degC}
    assert c._get_zoomed_temperature_limits([data], (500, 300)) == (-32, 52)


def test_zoomed_skewt_axes_match_its_limits():
    data = make_winter_profile("fzra")
    skew = zoomed_skewt(data, (1000, 700))
    left, right = skew.ax.get_xlim()
    assert (left, right) == pytest.approx(c._get_zoomed_temperature_limits([data], (1000, 700), right_margin=0.04))
    assert list(skew.ax.get_yticks()) == c._get_pressure_ticks((1000, 700))


THIN_VIEW_BUG = pytest.mark.xfail(strict=True, reason="PLOT_ISSUES.md: on a thin or inverted zoomed view the x-axis "
                                   "reaches ZOOM_MAX_SPAN and the moist adiabats' wet-bulb potential temperatures at "
                                   "the hot corner blow up, giving hundreds of lines")


@pytest.mark.parametrize(("case", "ylim", "has_mixing_lines"), [
    ("fzra", (1000, 700), True),
    ("fzra", (500, 300), False),
    pytest.param("fzra", (1000, 900), True, marks=THIN_VIEW_BUG),
    pytest.param("pl", (1000, 700), True, marks=THIN_VIEW_BUG),
])
def test_zoomed_background_lines_are_spaced_to_the_view(case, ylim, has_mixing_lines):
    lines = c._get_zoomed_background_lines(zoomed_skewt(make_winter_profile(case), ylim))
    for kind in ("dry", "moist"):
        assert 3 <= len(lines[kind]["t0"]) <= 16
        assert lines[kind]["pressure"].m.max() == pytest.approx(ylim[0])
    # mixing lines only below 600 hPa, as on the full skew-t
    assert (lines["mixing"]["mixing_ratio"].size > 0) == has_mixing_lines
    assert lines["mixing"]["mixing_ratio"].size <= 12


@pytest.mark.parametrize(("sfc_T", "bounds"), [(-25, (-52, 32)), (-15, (-42, 42)), (5, (-32, 52)), (25, (-22, 62))])
def test_full_skewt_temperature_bounds_follow_the_surface(sfc_T, bounds):
    assert c._get_skewt_temperature_bounds(sfc_T * units.degC) == bounds


def test_height_markers_suit_the_view(severe_params):
    assert c._get_height_markers(severe_params, (1000, 100)) == [2000, 4000, 6000, 8000, 10000, 12000, 14000, 16000]
    markers = c._get_height_markers(severe_params, (1000, 850))
    assert markers and np.all(np.diff(markers) == 500)
    assert c._get_height_markers(severe_params, (1100, 1050)) == []


# -----------------------------------------------------------------------------
# SURFACE MODIFICATION LABEL
# -----------------------------------------------------------------------------

@pytest.mark.parametrize(("modify_sfc", "expected"), [
    ({"T": 35, "Td": 20}, "T 35°C, Td 20°C · mixed to 1.0 km"),
    ({"T": 31}, "T 31°C · mixed to 250 m"),
    ({"T": 20}, "T 20°C · surface only"),
    ({"T": 35, "method": "replace"}, "T 35°C · surface only"),
    ({"T": 35, "method": "blend", "depth": 250}, "T 35°C · blended to 250 m"),
    ({"T": 35, "method": "blend"}, "T 35°C · blended to 500 m"),
    ({"ws": 10, "wd": 270}, "10 kt @ 270° · surface only"),
    ({"ws": 10}, "surface only"),
])
def test_describe_surface_modification(clean_profile, modify_sfc, expected):
    assert c._describe_surface_modification(clean_profile, modify_sfc) == expected


def test_modified_surface_is_marked_on_the_skewt(clean_profile):
    modify_sfc = {"T": 35}
    data = c._prepare_sounding_data(clean_profile, modify_sfc)
    skew = c._setup_skewt(plt.figure(), data["T"], THEME)
    c._plot_skewt_surface_annotations(skew, data, THEME, modify_sfc=modify_sfc, original=clean_profile, temperature_unit="degC")
    texts = texts_on(skew.ax)
    assert "35°C*" in texts and "25°C" in texts
    assert any(text.startswith("* MODIFIED SURFACE\nT 35°C") for text in texts)


# -----------------------------------------------------------------------------
# WATERMARK
# -----------------------------------------------------------------------------

def test_set_watermark_changes_and_restores_settings(watermark):
    set_watermark(title="MY LAB", alpha=0.2)
    assert watermark["title"] == "MY LAB" and watermark["alpha"] == 0.2
    set_watermark()
    assert watermark == WATERMARK_DEFAULTS


def test_set_watermark_rejects_unknown_settings(watermark):
    with pytest.raises(ValueError, match="Invalid watermark setting"):
        set_watermark(font="Comic Sans")
    assert watermark == WATERMARK_DEFAULTS


def test_credits_follow_the_watermark(watermark):
    fig = plt.figure()
    c._plot_credits(fig, THEME)
    assert [text.get_text() for text in fig.texts] == list(WATERMARK_TEXT["full"])

    set_watermark(title="MY LAB", subtitle="", color="red")
    fig = plt.figure()
    c._plot_credits(fig, THEME)
    assert fig.texts[0].get_text() == "MY LAB" and fig.texts[0].get_color() == "red"

    set_watermark(show=False)
    fig = plt.figure()
    c._plot_credits(fig, THEME)
    ax = c._plot_logo(fig.add_subplot())
    assert not fig.texts and not ax.images


def test_long_watermark_text_is_shrunk_to_fit(watermark):
    set_watermark(title="A VERY LONG CUSTOM WATERMARK TITLE " * 6)
    fig = plt.figure()
    c._plot_credits(fig, THEME)
    assert fig.texts[0].get_fontsize() < c.FONT_SIZES["watermark"]


# -----------------------------------------------------------------------------
# MISSING VALUES
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("value", [np.nan, ma.masked, np.nan * units.m, float("nan")])
def test_format_value_shows_missing_values(value):
    assert c._format_value(value, " m") == MISSING_VALUE


def test_format_value_rounds():
    assert c._format_value(12.345, " m", decimals=1) == "12.3 m"
    assert c._format_value(1234.0 * units.m, " m") == "1234 m"


def test_missing_profile_details():
    assert c._get_composite_parameter(None, "thermo", "mlcape", " {jkg}", None) == MISSING_VALUE
    assert c._get_profile_location({}) == MISSING_VALUE
    assert c._format_profile_time(None) == MISSING_VALUE
    assert c._format_profile_time(["2024", "05"]) == MISSING_VALUE
    assert not c._is_known("no-site-id") and not c._is_known("none") and c._is_known("KOUN")


@pytest.mark.parametrize(("eil", "expected"), [(None, False), ((ma.masked, ma.masked), False), ((900.0, ma.masked), False),
                                               ((900.0, 800.0), True), ((800.0, 900.0), False)])
def test_has_eil(eil, expected):
    assert c._has_eil({} if eil is None else {"eil": eil}) is expected


def test_hodograph_data_drops_missing_levels():
    snd = {"intrp": {"zINTRP": np.array([0.0, 100, 200, 300, 13000]), "uINTRP": np.array([1.0, np.nan, 3, 4, 5]),
                     "vINTRP": np.array([1.0, 2, np.inf, 4, 5])}}
    data = c._prepare_hodograph_data(snd)
    np.testing.assert_array_equal(data["z"], [0, 300])


def test_tables_show_missing_values(severe_params):
    snd = copy.deepcopy(severe_params)
    snd["thermo"]["mucape"] = ma.masked
    snd["thermo"]["dcape"] = np.nan
    snd["kinem"]["srh_0_to_1000"] = ma.masked
    fig, (thermo, kinem) = plt.subplots(2)
    c._plot_thermo_table(thermo, snd, THEME)
    c._plot_kinematics_table(kinem, snd, THEME)
    assert texts_on(thermo).count(MISSING_VALUE) >= 2
    assert MISSING_VALUE in texts_on(kinem)


def test_panels_warn_without_a_storm_motion(severe_params):
    snd = copy.deepcopy(severe_params)
    snd["kinem"]["sm_u"] = snd["kinem"]["sm_v"] = ma.masked
    ax = plt.figure().add_subplot()
    with pytest.warns(Warning, match="Streamwiseness could not be plotted"):
        c._plot_streamwiseness_panel(ax, snd, THEME)
    hodo = c._setup_hodograph(plt.figure().add_subplot(), THEME, (-50, 50), (-50, 50))
    n_patches = len(hodo.ax.patches)
    c._plot_hodograph_srh(hodo, snd, THEME)
    assert len(hodo.ax.patches) == n_patches  # no SRH fill


# -----------------------------------------------------------------------------
# PANELS AND TABLES (each draws on its own axis)
# -----------------------------------------------------------------------------

@pytest.mark.parametrize("draw", [c._plot_thermo_table, c._plot_moisture_table, c._plot_kinematics_table, c._plot_indices_box,
                                  c._plot_wind_parameter_box, c._plot_hodograph_parameter_box, c._plot_streamwiseness_panel,
                                  c._plot_vorticity_panel, c._plot_sr_wind_panel, c._plot_theta_panel,
                                  c._plot_wind_component_panel, c._plot_cape_cin_panel, c._plot_temperature_advection])
def test_panels_and_tables_draw(draw, severe_params):
    ax = plt.figure().add_subplot()
    assert draw(ax, severe_params, THEME) is ax
    assert ax.texts or ax.lines or ax.collections or ax.patches


def test_parcel_table_rows():
    snd = {"thermo": {f"{key}{name}": 100.0 for key in ("sb", "ml", "mu")
                      for name in ("_ecape", "cape", "6cape", "3cape", "cin", "_lcl_z")}}
    table = c._get_parcel_table(snd, THEME)
    assert [row["label"] for row in table["rows"]] == ["SB:", "ML:", "MU:"]
    assert all(len(row["values"]) == len(table["headers"]) for row in table["rows"])


# -----------------------------------------------------------------------------
# HODOGRAPH
# -----------------------------------------------------------------------------

def test_hodograph_bounds_include_the_origin():
    bounds = c._get_hodograph_bounds({"u": np.array([10.0, 30]), "v": np.array([5.0, 20])})
    assert bounds == ((0 - 45, 30 + 30), (0 - 45, 20 + 30))
    assert c._get_hodograph_bounds({"u": np.array([]), "v": np.array([])}) == ((-100, 100), (-100, 100))
    assert c._get_hodograph_bounds({"u": np.array([1.0]), "v": np.array([1.0])}, dynamic=False) == ((-100, 100), (-100, 100))


def test_storm_motion_text():
    assert c._format_motion("RM", np.nan, 5) == "RM: ---"
    assert c._format_motion("RM", 0.0, 10.0).endswith("@ 10 kt")
    assert c._storm_motion_name("right_moving") == "RIGHT MOVING"
    assert c._storm_motion_name([270, 30]) == "USER DEFINED"


def test_make_storm_relative_shifts_winds_and_motions(severe_params):
    relative = c._make_storm_relative(severe_params)
    kinem = severe_params["kinem"]
    assert (relative["kinem"]["sm_u"], relative["kinem"]["sm_v"]) == (0.0, 0.0)
    np.testing.assert_allclose(relative["intrp"]["uINTRP"], severe_params["intrp"]["uINTRP"] - kinem["sm_u"])
    np.testing.assert_allclose(relative["kinem"]["sm_rm"], np.subtract(kinem["sm_rm"], (kinem["sm_u"], kinem["sm_v"])))
    assert severe_params["kinem"]["sm_u"] == kinem["sm_u"]  # the input is not changed


def test_hodograph_srh_fill_layer(severe_params):
    hodo = c._setup_hodograph(plt.figure().add_subplot(), THEME, (-50, 50), (-50, 50))
    c._plot_hodograph_srh(hodo, severe_params, THEME, layer=(0, 1000))
    assert hodo.ax.patches[-1].get_label() == "0-1 km SRH"
    with pytest.raises(ValueError, match="'layer' must be"):
        c._plot_hodograph_srh(hodo, severe_params, THEME, layer="deep")


def test_hodograph_boundary_needs_matching_angles_and_colors(severe_params):
    hodo = c._setup_hodograph(plt.figure().add_subplot(), THEME, (-50, 50), (-50, 50))
    with pytest.raises(ValueError, match="same number of angles and colors"):
        c._plot_hodograph_boundary(hodo, {"angle": [45, 90], "color": ["k"]}, severe_params, THEME)


# -----------------------------------------------------------------------------
# PARCELS, LABELS, AND FIGURE UTILITIES
# -----------------------------------------------------------------------------

def test_special_parcel_codes():
    assert c._validate_special_parcels([["mu_ps_ecape"], ["sb_ia_cape"]]) == ["mu_ps_ecape", "sb_ia_cape"]
    assert c._format_special_parcel_name("mu_ps_ecape") == "MU PSEUDOADIABATIC ECAPE"
    with pytest.raises(ValueError, match="Invalid parcel code"):
        c._validate_special_parcels([["mu_ps_ecape", "xx_ps_cape"], []])


def test_spread_label_pressures_keeps_labels_apart():
    skew = c._setup_skewt(plt.figure(), np.array([20.0]) * units.degC, THEME)
    pressures = c._spread_label_pressures(skew.ax, [850.0, 845.0, 840.0, 500.0], gap=20)
    ys = sorted(skew.ax.transData.transform((0, p))[1] for p in pressures)
    assert np.all(np.diff(ys) >= 20 - 1e-6)
    assert pressures[3] == pytest.approx(500.0)  # a label with room stays at its own pressure


def test_scale_figure_scales_text_and_rejects_nonpositive():
    fig = plt.figure()
    text = fig.text(0.5, 0.5, "x", fontsize=10)
    c._scale_figure(fig, font_scale=1.5)
    assert text.get_fontsize() == pytest.approx(15)
    with pytest.raises(ValueError, match="must be above 0"):
        c._scale_figure(fig, font_scale=0)


# -----------------------------------------------------------------------------
# COMPOSITES
# -----------------------------------------------------------------------------

@pytest.mark.parametrize(("index", "letter"), [(0, "A"), (25, "Z"), (26, "AA"), (27, "AB"), (701, "ZZ"), (702, "AAA")])
def test_profile_letters(index, letter):
    assert c._get_profile_letter(index) == letter


def test_profile_styles_and_colors():
    assert c._get_profile_styles(None, "-", 3) == ["-", "-", "-"]
    assert c._get_profile_styles("none", "-", 2) == ["-", "-"]
    assert c._get_profile_styles(("--", ":"), "-", 2) == ["--", ":"]
    assert c._get_profile_colors(["red", "blue"], "viridis", 2) == ["red", "blue"]
    assert len(c._get_profile_colors(None, "viridis", 4)) == 4
    assert c._get_profile_colors("none", "viridis", 2, dark_mode=True) != c._get_profile_colors("none", "viridis", 2)


@pytest.mark.parametrize(("site", "expected"), [({"source": "RAOB OBSERVED PROFILE"}, "RAOB"),
                                                ({"source": "ACARS OBSERVED PROFILE"}, "ACARS"),
                                                ({"source": "BUFKIT", "model": "hrrr", "fcst-hour": "F03"}, "HRRR F03"),
                                                ({"source": "BUFKIT", "model": "hrrr", "fcst-hour": "no-fcst-hour"}, "HRRR"),
                                                ({"source": "none"}, "PROFILE")])
def test_profile_type(site, expected):
    assert c._get_profile_type(site) == expected


def test_profile_location():
    site = {"site-id": "KOUN", "site-name": "NORMAN", "site-latlon": [35.2, -97.4]}
    assert c._get_profile_location(site) == "KOUN, NORMAN"
    assert c._get_profile_location(site, short=True) == "KOUN"
    assert c._get_profile_location({"site-id": "no-site-id", "site-latlon": [35.2, -97.4]}) == "35.2, -97.4"


def test_composite_titles_sort_times_of_either_type():
    first, second = make_clean_profile(), make_clean_profile()
    first["site_info"]["valid-time"] = [2024, 5, 21, 6]  # ints, as some sources give them
    second["site_info"]["valid-time"] = ["2024", "05", "21", "00"]
    titles = c._get_composite_titles([first, second], style="hodograph", sr_hodo=True)
    assert titles["top_title"] == "COMPOSITE STORM RELATIVE HODOGRAPH"
    assert titles["left_title"] == "2 PROFILES | VALID: 00Z | 05-21-2024 - 6Z | 5-21-2024"
    assert titles["right_title"] == "TEST, Synthetic SounderPy Test Profile"


def test_profile_winds_are_interpolated_to_regular_heights(clean_profile):
    winds = c._get_profile_winds(clean_profile, step=100, max_height=9000)
    np.testing.assert_array_equal(winds["z"], np.arange(0, 9001, 100))
    assert winds["u"][0] == pytest.approx(clean_profile["u"][0].m)


def test_composite_storm_relative_winds():
    winds = [{"z": np.array([0, 100]), "u": np.array([10.0, 20.0]), "v": np.array([5.0, 5.0])}]
    params = [{"kinem": {"sm_u": 4.0, "sm_v": -1.0}}]
    shifted = c._make_winds_storm_relative(winds, params)[0]
    np.testing.assert_allclose(shifted["u"], [6.0, 16.0])
    np.testing.assert_allclose(shifted["v"], [6.0, 6.0])

    missing = [None, {"kinem": {"sm_u": ma.masked, "sm_v": ma.masked}}, {"kinem": {"sm_u": np.nan}}, params[0]]
    assert c._get_missing_storm_motions(missing) == ["A", "B", "C"]


def test_composite_table_without_parameters(clean_profile, second_clean_profile):
    ax = plt.figure().add_subplot()
    c._plot_composite_table(ax, [clean_profile, second_clean_profile], None, ["red", "blue"], THEME)
    texts = texts_on(ax)
    assert {"PROFILE", "LOCATION", "VALID", "A", "B"} <= set(texts)
    assert "MLCAPE" not in texts



def test_indices_box_shows_the_pbl(severe_params):
    snd = copy.deepcopy(severe_params)
    snd["thermo"]["pbl_depth"], snd["kinem"]["mean_wind_pbl"] = 1078.0, (7.53, 5.74)  # blowing toward the NE
    ax = plt.figure().add_subplot()
    c._plot_indices_box(ax, snd, THEME)
    texts = texts_on(ax)
    assert "  PBL depth," in texts and "1078m" in texts and "SW 9kt" in texts


def test_indices_box_without_a_pbl(severe_params):
    snd = copy.deepcopy(severe_params)
    snd["thermo"]["pbl_depth"], snd["kinem"]["mean_wind_pbl"] = np.ma.masked, (np.ma.masked, np.ma.masked)
    ax = plt.figure().add_subplot()
    c._plot_indices_box(ax, snd, THEME)
    assert texts_on(ax).count(MISSING_VALUE) == 2

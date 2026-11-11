"""End-to-end figures: build_sounding (severe and winter), build_hodograph, build_composite, build_vad_hodograph.

Each render saves a small PNG offline (``map_zoom=0`` turns off the map and radar). Renders are marked
``slow``; skip them with ``pytest -m "not slow"``. Only one severe sounding draws the default MUECAPE parcel
(about 15 s in the ecape_parcel library); the others pass ``special_parcels="simple"``.
"""

from __future__ import annotations

import copy
import warnings
from datetime import datetime, timezone

import numpy as np
import pytest

import sounderpy as spy
from sounderpy.plot import _resolve_sounding_mode, operational

from helpers import WINTER_CASES, assert_png, make_winter_profile


# a fixed surface observation, so the winter station model can be drawn offline
FIXED_OBS = {"station": "KSYN", "name": "SYNTHETIC", "distance_km": 3.2,
             "time": datetime(2026, 1, 1, 0, 0, tzinfo=timezone.utc), "offset_min": 0.0,
             "tmpf": 27.0, "dwpf": 25.0, "drct": 20.0, "sknt": 12.0, "gust": 22.0, "vsby": 1.5,
             "wxcodes": "-FZRA BR", "mslp": 1018.2, "skyc1": "BKN", "skyc2": "OVC", "skyc3": None,
             "timing": {"find_station": 0.0, "fetch": 0.0}}


def build(function, data, path, **kwargs):
    """Save a figure to ``path`` with the offline test defaults, check the PNG, and return the function's result."""
    if function is not spy.build_composite:
        kwargs.setdefault("map_zoom", 0)
    if function in (spy.build_sounding, spy.build_hodograph):
        kwargs.setdefault("dpi", 50)
    result = function(data, save=True, filename=str(path), **kwargs)
    assert_png(path)
    return result


# -----------------------------------------------------------------------------
# SEVERE SOUNDING AND HODOGRAPH
# -----------------------------------------------------------------------------

@pytest.mark.slow
def test_build_sounding(clean_profile, tmp_path):
    # the public wrapper saves or shows the figure and returns nothing
    assert build(spy.build_sounding, clean_profile, tmp_path / "sounding.png") is None


@pytest.mark.slow
def test_build_sounding_dark_colorblind(clean_profile, tmp_path):
    build(spy.build_sounding, clean_profile, tmp_path / "sounding_dark_cb.png", color_blind=True, dark_mode=True,
          special_parcels="simple")


@pytest.mark.slow
def test_build_sounding_shared_options(clean_profile, tmp_path):
    build(spy.build_sounding, clean_profile, tmp_path / "sounding_options.png", special_parcels="simple",
          modify_sfc={"T": 33, "Td": 22, "ws": 20, "wd": 180}, max_levels=40, ylim=(1000, 500), show_theta=True,
          theme={"temperature": "black", "parcel": {"mu": "navy"}})


@pytest.mark.slow
@pytest.mark.parametrize("problem", ["shallow", "missing_levels"])
def test_build_sounding_with_incomplete_data(problem, clean_profile, tmp_path):
    if problem == "shallow":  # 0-4 km, so no 0-6 km storm motion
        data = {key: value[:17] if hasattr(value, "shape") else value for key, value in clean_profile.items()}
    else:
        data = clean_profile
        for key in ("T", "Td", "u", "v"):
            values = data[key].m.copy()
            values[[5, 6, 20]] = np.nan
            data[key] = values * data[key].units
    build(spy.build_sounding, data, tmp_path / f"sounding_{problem}.png", special_parcels="simple")


@pytest.mark.slow
def test_build_hodograph(clean_profile, tmp_path):
    assert build(spy.build_hodograph, clean_profile, tmp_path / "hodograph.png") is None


@pytest.mark.slow
def test_build_hodograph_storm_relative(clean_profile, tmp_path):
    build(spy.build_hodograph, clean_profile, tmp_path / "hodograph_sr.png", storm_motion="right_moving", sr_hodo=True)


@pytest.mark.slow
def test_build_vad_hodograph(vad_profile, tmp_path):
    assert build(spy.build_vad_hodograph, vad_profile, tmp_path / "vad_hodograph.png") is None


# -----------------------------------------------------------------------------
# COMPOSITES
# -----------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.parametrize("style", ["full", "sounding", "hodograph"])
def test_build_composite(style, clean_profile, second_clean_profile, tmp_path):
    result = build(spy.build_composite, [clean_profile, second_clean_profile], tmp_path / f"composite_{style}.png",
                   style=style, shade_between=True)
    assert result is None


@pytest.mark.slow
@pytest.mark.parametrize("style", ["full", "hodograph"])
def test_build_composite_storm_relative(style, clean_profile, second_clean_profile, tmp_path):
    build(spy.build_composite, [clean_profile, second_clean_profile], tmp_path / f"composite_sr_{style}.png",
          style=style, sr_hodo=True, show_params=False)


def test_build_composite_rejects_unknown_style(clean_profile):
    with pytest.raises(ValueError, match="style"):
        spy.build_composite([clean_profile], style="skewt")


# -----------------------------------------------------------------------------
# MODE (severe or winter)
# -----------------------------------------------------------------------------

def test_auto_mode_uses_surface_wet_bulb(clean_profile, capsys):
    # a warm, moist surface (29 deg C) is severe; every idealized winter profile is below 35 deg F wet-bulb
    assert _resolve_sounding_mode("auto", clean_profile) == "severe"
    for case in WINTER_CASES:
        assert _resolve_sounding_mode("auto", make_winter_profile(case)) == "winter"
    assert "AUTO MODE" in capsys.readouterr().out


def test_explicit_modes_skip_the_auto_check(clean_profile):
    assert _resolve_sounding_mode("severe", make_winter_profile("fzra")) == "severe"
    assert _resolve_sounding_mode("winter", clean_profile) == "winter"


def test_invalid_mode_is_rejected(clean_profile):
    with pytest.raises(ValueError, match="'mode' must be"):
        spy.build_sounding(clean_profile, mode="snow")


# -----------------------------------------------------------------------------
# WINTER SOUNDING
# -----------------------------------------------------------------------------

@pytest.mark.slow
def test_build_winter_sounding(winter_profile, tmp_path):
    case, data = winter_profile
    build(spy.build_sounding, data, tmp_path / f"winter_{case}.png", mode="winter")


@pytest.mark.slow
def test_build_winter_sounding_dark_colorblind(tmp_path):
    build(spy.build_sounding, make_winter_profile("fzra"), tmp_path / "winter_dark_cb.png", mode="winter",
          color_blind=True, dark_mode=True)


@pytest.mark.slow
def test_build_winter_sounding_shared_options(tmp_path):
    # the sleet profile warms with height in (1000, 600), which once made the zoomed skew-t raise StopIteration
    build(spy.build_sounding, make_winter_profile("pl"), tmp_path / "winter_options.png", mode="winter",
          modify_sfc={"T": -6, "Td": -7, "ws": 25, "wd": 20}, max_levels=20, ylim=(1000, 600),
          theme={"melt": "orange", "border": "navy", "parcel": {"mu": "black"}})


def test_winter_theme_rejects_unknown_settings(tmp_path):
    with pytest.raises(ValueError, match="Invalid theme setting"):
        spy.build_sounding(make_winter_profile("fzra"), mode="winter", map_zoom=0, save=True,
                           filename=str(tmp_path / "x.png"), theme={"not_a_color": "k"})


@pytest.mark.slow
def test_winter_sounding_warns_about_severe_only_options(tmp_path):
    with pytest.warns(Warning, match="winter sounding doesn't use show_theta, sr_hodo"):
        build(spy.build_sounding, make_winter_profile("fzra"), tmp_path / "winter.png", mode="winter",
              show_theta=True, sr_hodo=True)


@pytest.mark.slow
def test_winter_sounding_without_severe_options_does_not_warn(tmp_path):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        build(spy.build_sounding, make_winter_profile("fzra"), tmp_path / "winter.png", mode="winter")
    assert not [w for w in caught if "winter sounding doesn't use" in str(w.message)]


@pytest.mark.slow
def test_severe_sounding_warns_about_show_obs(clean_profile, tmp_path):
    with pytest.warns(Warning, match="show_obs is only used by the winter sounding"):
        build(spy.build_sounding, clean_profile, tmp_path / "severe.png", special_parcels="simple", show_obs=True)


# -----------------------------------------------------------------------------
# WINTER STATION MODEL (show_obs, with the observation mocked)
# -----------------------------------------------------------------------------

@pytest.mark.slow
def test_winter_sounding_draws_station_model(monkeypatch, tmp_path):
    requested = []

    def fake_sounding_obs(clean_data):
        requested.append(clean_data["site_info"]["site-id"])
        return copy.deepcopy(FIXED_OBS)

    monkeypatch.setattr(operational, "sounding_obs", fake_sounding_obs)
    build(spy.build_sounding, make_winter_profile("fzra"), tmp_path / "winter_obs.png", mode="winter", show_obs=True)
    assert requested == ["SYN"]


@pytest.mark.slow
def test_station_model_with_missing_values(monkeypatch, tmp_path):
    # a report with only wind, visibility and weather (KFET, 2024-01-12 1155Z): T, Td and SLP read NA
    obs = dict(copy.deepcopy(FIXED_OBS), tmpf=None, dwpf=None, mslp=None, vsby=None, skyc1=None, skyc2=None)
    monkeypatch.setattr(operational, "sounding_obs", lambda clean_data: obs)
    build(spy.build_sounding, make_winter_profile("fzra"), tmp_path / "winter_obs_missing.png", mode="winter", show_obs=True)


@pytest.mark.slow
def test_winter_sounding_without_obs_found(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(operational, "sounding_obs", lambda clean_data: None)
    build(spy.build_sounding, make_winter_profile("fzra"), tmp_path / "winter_no_obs.png", mode="winter", show_obs=True)
    assert "SFC OB: none found" in capsys.readouterr().out


# -----------------------------------------------------------------------------
# SPARSE HISTORICAL SOUNDINGS
# -----------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.parametrize("moisture, winds", [("none", False), ("low", True)])
def test_build_sounding_sparse_old_profile(moisture, winds, tmp_path):
    # a 1958/1965-style sounding: mandatory levels only, no (or only low-level) dewpoints, maybe no winds.
    # The default parcels (MUECAPE trace included) must fall back to missing values, not crash
    from helpers import make_sparse_old_profile
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        build(spy.build_sounding, make_sparse_old_profile(moisture, winds), tmp_path / f"old_{moisture}.png")
        build(spy.build_hodograph, make_sparse_old_profile(moisture, winds), tmp_path / f"old_{moisture}_hodo.png")

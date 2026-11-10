"""Tests for the winter sounding: winter calc, mode selection, figure, and surface observation.

The calc and figure tests use the idealized profiles in ``conftest.py`` (``WINTER_CASES``) and run
offline. The live surface-observation test contacts the Iowa Environmental Mesonet and only runs with::

    SOUNDERPY_RUN_NETWORK_TESTS=1 pytest -m network -v
"""

from __future__ import annotations

import copy
import os
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

import sounderpy as spy
from sounderpy.calc_winter import (MELT_ALOFT, REFREEZE_MIX, REFREEZE_PL, REFREEZE_SLOPE, dendritic_growth_zone,
                                   mixed_layer, precip_phase_trace, winter_params)
from sounderpy.data import sfc_obs
from sounderpy.plot import _resolve_sounding_mode
from sounderpy.plot import operational

from conftest import make_winter_profile


RUN_NETWORK = os.getenv("SOUNDERPY_RUN_NETWORK_TESTS", "0").lower() in {"1", "true", "yes", "on"}

network_only = pytest.mark.skipif(
    not RUN_NETWORK,
    reason="Set SOUNDERPY_RUN_NETWORK_TESTS=1 to run live retrieval tests.",
)


# precip type (Bourgouin), SPC best guess, the phase trace from source to surface, and freezing drizzle level
EXPECTED = {
    "fzra": ("FZRA", "Freezing Rain", ["SN", "RA", "FZRA"], "LOW"),
    "pl": ("PL", "Sleet", ["SN", "RA", "PL"], "LOW"),
    "rasn": ("RA/SN", "Snow", ["SN", "RA/SN"], "NONE"),
    "fzdz": ("FZDZ", "Freezing Rain/Drizzle", ["DZ", "FZDZ"], "HIGH"),
}

# a fixed surface observation, so the station model can be drawn offline
FIXED_OBS = {"station": "KSYN", "name": "SYNTHETIC", "distance_km": 3.2,
             "time": datetime(2026, 1, 1, 0, 0, tzinfo=timezone.utc), "offset_min": 0.0,
             "tmpf": 27.0, "dwpf": 25.0, "drct": 20.0, "sknt": 12.0, "gust": 22.0, "vsby": 1.5,
             "wxcodes": "-FZRA BR", "mslp": 1018.2, "skyc1": "BKN", "skyc2": "OVC", "skyc3": None,
             "timing": {"find_station": 0.0, "fetch": 0.0}}


def assert_nonempty_png(path: Path):
    assert path.exists(), f"expected plot was not created: {path}"
    assert path.stat().st_size > 1_000, f"plot file is unexpectedly small: {path.stat().st_size} bytes"
    with path.open("rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n"


# -----------------------------------------------------------------------------
# WINTER CALC
# -----------------------------------------------------------------------------

def test_winter_calc_precip_type(winter_profile):
    case, data = winter_profile
    ptype, sharppy_ptype, trace, fzdz_level = EXPECTED[case]

    snd = winter_params(data).calc()
    winter = snd["winter"]

    # the standard calc is still there, with winter values added
    assert {"general", "thermo", "kinem", "intrp", "winter"} <= set(snd)
    assert winter["ptype"] == ptype
    assert winter["sharppy_ptype"] == sharppy_ptype
    assert [step["phase"] for step in winter["phase_trace"]] == trace
    assert winter["fzdz"]["level"] == fzdz_level
    assert winter["has_omega"]


def test_winter_params_is_public_and_takes_the_sounding_params_options():
    data = make_winter_profile("fzra")
    snd = spy.winter_params(data, storm_motion="left_moving", include_all_parcels=True).calc()
    assert {"general", "thermo", "kinem", "intrp", "winter"} <= set(snd)

    # a modified surface reaches the winter values: a warmer surface turns freezing rain to rain
    warm = spy.winter_params(make_winter_profile("fzra"), modify_sfc={"T": 4, "Td": 3}).calc()["winter"]
    assert warm["sfc_T_F"] > 32
    assert warm["ptype"] != "FZRA"


def test_winter_calc_energies_follow_bourgouin():
    # freezing rain: the warm nose melts snow completely, the surface layer is too shallow to refreeze it
    winter = winter_params(make_winter_profile("fzra")).calc()["winter"]
    melt = [lyr for lyr in winter["layers"] if lyr["kind"] == "warm"]
    refreeze = [lyr for lyr in winter["layers"] if lyr["kind"] == "cold" and lyr["surface"]]
    assert len(melt) == 1 and melt[0]["energy"] > MELT_ALOFT[1]
    assert len(refreeze) == 1 and -refreeze[0]["energy"] < REFREEZE_MIX + REFREEZE_SLOPE * melt[0]["energy"]

    # sleet: a deep cold layer refreezes the melted drops
    winter = winter_params(make_winter_profile("pl")).calc()["winter"]
    melt = [lyr for lyr in winter["layers"] if lyr["kind"] == "warm"][0]
    refreeze = [lyr for lyr in winter["layers"] if lyr["kind"] == "cold" and lyr["surface"]][0]
    assert -refreeze["energy"] > REFREEZE_PL + REFREEZE_SLOPE * melt["energy"]


def test_precip_phase_trace_thresholds():
    # (kind, energy) layers top first, as melting_refreezing_layers returns them; the last one is at the surface.
    # the source is snow (-20 °C)
    def surface_type(*layers):
        layers = [{"kind": kind, "energy": energy, "bottom": 1000.0 - 100 * i, "surface": i == len(layers) - 1}
                  for i, (kind, energy) in enumerate(layers)]
        return precip_phase_trace(layers, -20.0, 3000.0, -20.0, 2.0 if layers[-1]["kind"] == "warm" else -2.0)[1]

    # warm layer at the surface (Bourgouin 2000): snow below 5.6 J/kg, a mix to 13.2, rain above
    assert surface_type(("warm", 4.0)) == "SN"
    assert surface_type(("warm", 9.0)) == "RA/SN"
    assert surface_type(("warm", 20.0)) == "RA"
    # warm layer aloft (Birk et al. 2021): frozen below 9.4 J/kg, partly melted (refreezes to sleet) to 15
    assert surface_type(("warm", 9.0), ("cold", -50)) == "SN"
    assert surface_type(("warm", 12.0), ("cold", -50)) == "PL"
    # fully melted (ME 50): liquid below RE 121, a mix to 212, ice pellets above
    assert surface_type(("warm", 50.0), ("cold", -100)) == "FZRA"
    assert surface_type(("warm", 50.0), ("cold", -150)) == "FZRA/PL"
    assert surface_type(("warm", 50.0), ("cold", -250)) == "PL"


def test_mixed_layer_top_from_theta_v():
    # constant theta-v to 1000 m, stable (+5 K/km) above: the top is just above 1000 m, where theta-v passes the
    # layer mean + 0.5 K
    z = np.arange(0, 3001, 50.0)
    p = 1000 * np.exp(-z / 8000)
    thetav = np.where(z <= 1000, 270.0, 270.0 + 5 * (z - 1000) / 1000)
    assert mixed_layer(p, z, thetav) == pytest.approx(1100, abs=10)

    # a surface inversion (+10 K/km from the ground) leaves only a shallow layer
    assert mixed_layer(p, z, 265 + 10 * z / 1000) < 100

    # never warms past the mean: no top
    assert np.ma.is_masked(mixed_layer(p, z, np.full_like(z, 270.0)))


def test_mixed_layer_values_on_the_winter_dict():
    winter = winter_params(make_winter_profile("rasn")).calc()["winter"]
    assert winter["ml_depth"] > 0
    assert np.isfinite(winter["ml_lapse_rate"])
    assert not any(np.ma.is_masked(c) for c in winter["mean_wind_ml"])
    assert not {"blsn_ratio", "dT_0_to_1000", "mean_wind_pbl"} & set(winter)


def test_dgz_is_the_layer_in_cloud():
    # a surface at -14 C cooling to -19 C, an inversion back to -8 C, then cooling through -12 to -17 C in cloud
    # aloft (like the OAX RAOB, 2024-01-12 12Z). The first -12 and -17 C levels give an inverted layer
    z = np.array([0, 300, 600, 1000, 1500, 2000, 2500, 3000, 4000])
    T = np.array([-14, -17.5, -19, -8, -9, -11, -14, -18, -25.0])
    p = 1000 * np.exp(-z / 8000)
    cloud_z = np.arange(0, 4001, 100.0)
    rh_ice = np.where((cloud_z >= 1800) & (cloud_z <= 3500), 100.0, 70.0)

    (bottom, top), (p_bottom, p_top) = dendritic_growth_zone(z, T, p, cloud_z, rh_ice)
    assert 2000 < bottom < top < 3000 and p_bottom > p_top

    # with no cloud anywhere, the lowest layer
    (bottom, top), _ = dendritic_growth_zone(z, T, p, cloud_z, np.full_like(cloud_z, 50.0))
    assert bottom < 300

    # never between -12 and -17 C: none
    assert dendritic_growth_zone(z, T + 30, p, cloud_z, rh_ice) is None


def test_snow_ratio_only_when_snow_reaches_the_surface():
    assert np.ma.is_masked(winter_params(make_winter_profile("fzra")).calc()["winter"]["slr"])
    data = make_winter_profile("rasn")
    slr = winter_params(data).calc()["winter"]["slr"]
    assert not np.ma.is_masked(slr) and slr > 0

    # Kuchera: 12 + 2 (271.16 - Tmax) above 271.16 K, 12 + (271.16 - Tmax) below it, no cap
    Tmax = float(np.max(data["T"].m[data["p"].m >= 500])) + 273.15
    expected = 12 + 2 * (271.16 - Tmax) if Tmax > 271.16 else 12 + (271.16 - Tmax)
    assert slr == pytest.approx(max(expected, 0.0))


# -----------------------------------------------------------------------------
# MODE
# -----------------------------------------------------------------------------

def test_auto_mode_uses_surface_wet_bulb(clean_profile, capsys):
    # a warm, moist surface (29 deg C) is severe; every idealized winter profile is below 35 deg F wet-bulb
    assert _resolve_sounding_mode("auto", clean_profile) == "severe"
    for case in EXPECTED:
        assert _resolve_sounding_mode("auto", make_winter_profile(case)) == "winter"
    assert "AUTO MODE" in capsys.readouterr().out


def test_explicit_modes_skip_the_auto_check(clean_profile):
    assert _resolve_sounding_mode("severe", make_winter_profile("fzra")) == "severe"
    assert _resolve_sounding_mode("winter", clean_profile) == "winter"


def test_invalid_mode_is_rejected(clean_profile):
    with pytest.raises(ValueError, match="'mode' must be"):
        spy.build_sounding(clean_profile, mode="snow")


# -----------------------------------------------------------------------------
# WINTER FIGURE
# -----------------------------------------------------------------------------

@pytest.mark.slow
def test_build_winter_sounding_saves_rendered_png(winter_profile, offline_plot_assets, tmp_path):
    case, data = winter_profile
    output = tmp_path / f"winter_{case}.png"

    spy.build_sounding(data, mode="winter", map_zoom=0, dpi=50, save=True, filename=str(output))

    assert_nonempty_png(output)


@pytest.mark.slow
def test_build_winter_sounding_dark_colorblind_variant(offline_plot_assets, tmp_path):
    output = tmp_path / "winter_dark_cb.png"

    spy.build_sounding(make_winter_profile("fzra"), mode="winter", color_blind=True, dark_mode=True,
                       map_zoom=0, dpi=50, save=True, filename=str(output))

    assert_nonempty_png(output)


@pytest.mark.slow
def test_winter_sounding_supports_shared_options(offline_plot_assets, tmp_path):
    output = tmp_path / "winter_options.png"

    # the sleet profile warms with height in (1000, 600), which once made the zoomed skew-t raise StopIteration
    spy.build_sounding(make_winter_profile("pl"), mode="winter", map_zoom=0, dpi=50, save=True, filename=str(output),
                       modify_sfc={"T": -6, "Td": -7, "ws": 25, "wd": 20}, max_levels=20, ylim=(1000, 600),
                       theme={"melt": "orange", "border": "navy", "parcel": {"mu": "black"}})

    assert_nonempty_png(output)


def test_winter_theme_rejects_unknown_settings(tmp_path):
    with pytest.raises(ValueError, match="Invalid theme setting"):
        spy.build_sounding(make_winter_profile("fzra"), mode="winter", map_zoom=0, save=True,
                           filename=str(tmp_path / "x.png"), theme={"not_a_color": "k"})


@pytest.mark.slow
def test_winter_sounding_warns_about_severe_only_options(offline_plot_assets, tmp_path):
    with pytest.warns(Warning, match="winter sounding doesn't use show_theta, sr_hodo"):
        spy.build_sounding(make_winter_profile("fzra"), mode="winter", map_zoom=0, dpi=50, save=True,
                           filename=str(tmp_path / "winter.png"), show_theta=True, sr_hodo=True)


@pytest.mark.slow
def test_winter_sounding_without_severe_options_does_not_warn(offline_plot_assets, tmp_path):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        spy.build_sounding(make_winter_profile("fzra"), mode="winter", map_zoom=0, dpi=50, save=True,
                           filename=str(tmp_path / "winter.png"))
    assert not [w for w in caught if "winter sounding doesn't use" in str(w.message)]


@pytest.mark.slow
def test_severe_sounding_warns_about_show_obs(clean_profile, offline_plot_assets, tmp_path):
    with pytest.warns(Warning, match="show_obs is only used by the winter sounding"):
        spy.build_sounding(clean_profile, map_zoom=0, dpi=50, save=True, filename=str(tmp_path / "severe.png"),
                           show_obs=True)


# -----------------------------------------------------------------------------
# SURFACE OBSERVATION
# -----------------------------------------------------------------------------

@pytest.mark.slow
def test_winter_sounding_draws_station_model(monkeypatch, offline_plot_assets, tmp_path):
    requested = []

    def fake_sounding_obs(clean_data):
        requested.append(clean_data["site_info"]["site-id"])
        return copy.deepcopy(FIXED_OBS)

    monkeypatch.setattr(operational, "sounding_obs", fake_sounding_obs)
    output = tmp_path / "winter_obs.png"

    spy.build_sounding(make_winter_profile("fzra"), mode="winter", show_obs=True, map_zoom=0, dpi=50, save=True,
                       filename=str(output))

    assert requested == ["SYN"]
    assert_nonempty_png(output)


@pytest.mark.slow
def test_station_model_with_missing_values(monkeypatch, offline_plot_assets, tmp_path):
    # a report with only wind, visibility and weather (KFET, 2024-01-12 1155Z): T, Td and SLP read NA
    obs = dict(copy.deepcopy(FIXED_OBS), tmpf=None, dwpf=None, mslp=None, vsby=None, skyc1=None, skyc2=None)
    monkeypatch.setattr(operational, "sounding_obs", lambda clean_data: obs)
    output = tmp_path / "winter_obs_missing.png"

    spy.build_sounding(make_winter_profile("fzra"), mode="winter", show_obs=True, map_zoom=0, dpi=50, save=True,
                       filename=str(output))

    assert_nonempty_png(output)


@pytest.mark.slow
def test_winter_sounding_without_obs_found(monkeypatch, offline_plot_assets, tmp_path, capsys):
    monkeypatch.setattr(operational, "sounding_obs", lambda clean_data: None)
    output = tmp_path / "winter_no_obs.png"

    spy.build_sounding(make_winter_profile("fzra"), mode="winter", show_obs=True, map_zoom=0, dpi=50, save=True,
                       filename=str(output))

    assert "SFC OB: none found" in capsys.readouterr().out
    assert_nonempty_png(output)


def test_sounding_valid_time_handles_hours_and_minutes():
    data = make_winter_profile("fzra")
    assert sfc_obs.sounding_valid_time(data) == datetime(2026, 1, 1, 0, 0, tzinfo=timezone.utc)

    # ACARS gives hour and minute
    data["site_info"]["valid-time"] = ["2024", "01", "15", "14:18"]
    assert sfc_obs.sounding_valid_time(data) == datetime(2024, 1, 15, 14, 18, tzinfo=timezone.utc)


def test_future_valid_time_has_no_observation(monkeypatch):
    def should_not_fetch(*args, **kwargs):
        raise AssertionError("a future valid time should not be fetched")

    monkeypatch.setattr(sfc_obs, "fetch_obs", should_not_fetch)
    monkeypatch.setattr(sfc_obs, "nearest_stations", lambda *args, **kwargs: [("KSYN", "SYNTHETIC", 1.0)])
    data = make_winter_profile("fzra")
    data["site_info"]["valid-time"] = [str(datetime.now(timezone.utc).year + 1), "01", "01", "00"]

    assert sfc_obs.sounding_obs(data) is None


@pytest.mark.network
@network_only
def test_live_surface_observation():
    # Grand Forks, ND, during freezing drizzle (HRRR BUFKIT 2024-01-24 12Z F01, valid 13Z)
    data = make_winter_profile("fzdz")
    data["site_info"]["site-latlon"] = [47.95, -97.18]
    data["site_info"]["valid-time"] = ["2024", "01", "24", "13"]

    obs = sfc_obs.sounding_obs(data)

    assert obs is not None
    assert obs["station"] == "KGFK"
    assert abs(obs["offset_min"]) <= 60
    assert "FZDZ" in obs["wxcodes"]

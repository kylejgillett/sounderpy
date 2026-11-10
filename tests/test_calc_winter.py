"""Winter calc (``calc_winter.py``): precipitation type, Bourgouin energies, mixed layer, DGZ, and snow ratio.

Uses the idealized profiles in ``helpers.py`` (``WINTER_CASES``).
"""

from __future__ import annotations

import numpy as np
import pytest

import sounderpy as spy
from sounderpy.calc_winter import (MELT_ALOFT, REFREEZE_MIX, REFREEZE_PL, REFREEZE_SLOPE, dendritic_growth_zone,
                                   mixed_layer, precip_phase_trace, winter_params)

from helpers import make_winter_profile


# precip type (Bourgouin), SPC best guess, the phase trace from source to surface, and freezing drizzle level
EXPECTED = {
    "fzra": ("FZRA", "Freezing Rain", ["SN", "RA", "FZRA"], "LOW"),
    "pl": ("PL", "Sleet", ["SN", "RA", "PL"], "LOW"),
    "rasn": ("RA/SN", "Snow", ["SN", "RA/SN"], "NONE"),
    "fzdz": ("FZDZ", "Freezing Rain/Drizzle", ["DZ", "FZDZ"], "HIGH"),
}


@pytest.mark.parametrize("case", EXPECTED)
def test_winter_calc_precip_type(case, winter_params_for):
    ptype, sharppy_ptype, trace, fzdz_level = EXPECTED[case]

    snd = winter_params_for(case)
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


def test_mixed_layer_values_on_the_winter_dict(winter_params_for):
    winter = winter_params_for("rasn")["winter"]
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

"""Calculations behind the figures: surface modification (``utils.modify_surface``) and ``sounding_params``' parcel_depth."""

from __future__ import annotations

import copy

import metpy.calc as mpcalc
import numpy as np
import pytest
from metpy.units import units

from sounderpy.calc import sounding_params
from sounderpy.utils import modify_surface


def changed_levels(before, after, key="T"):
    """Indices of the levels where a variable changed."""
    return np.flatnonzero(~np.isclose(before[key].to("degC").m, after[key].to("degC").m))


# -----------------------------------------------------------------------------
# modify_surface
# -----------------------------------------------------------------------------

def test_mix_raises_a_warmer_surface_up_its_dry_adiabat(clean_profile):
    # the clean profile's surface is 29 deg C
    modified = modify_surface(clean_profile, {"T": 35, "Td": 20})
    changed = changed_levels(clean_profile, modified)
    assert changed[0] == 0 and changed.size > 1
    assert np.array_equal(changed, np.arange(changed.size))  # one layer from the ground up

    # through the mixed layer: one potential temperature and one mixing ratio (the surface's)
    mixed = slice(0, changed.size)
    p = modified["p"][mixed]
    theta = mpcalc.potential_temperature(p, modified["T"][mixed]).m
    np.testing.assert_allclose(theta, theta[0], atol=0.01)
    w = mpcalc.saturation_mixing_ratio(p, modified["Td"][mixed]).m
    np.testing.assert_allclose(w, w[0], rtol=1e-3)

    # above it, nothing changes
    above = slice(changed.size, None)
    np.testing.assert_allclose(modified["T"][above].m, clean_profile["T"][above].m)


@pytest.mark.parametrize("modify_sfc", [{"T": 20}, {"Td": 10}, {"T": 35, "method": "replace"}])
def test_cooler_dewpoint_only_and_replace_change_only_the_surface(clean_profile, modify_sfc):
    modified = modify_surface(clean_profile, modify_sfc)
    for key in ("T", "Td"):
        assert set(changed_levels(clean_profile, modified, key)) <= {0}
    if "T" in modify_sfc:
        assert modified["T"][0].m == pytest.approx(modify_sfc["T"])


def test_blend_tapers_the_change_to_zero_at_depth(clean_profile):
    # levels every 250 m: the +6 deg C change is +3 at 250 m and gone at 500 m
    modified = modify_surface(clean_profile, {"T": 35, "method": "blend", "depth": 500})
    change = modified["T"].m - clean_profile["T"].m
    np.testing.assert_allclose(change[:3], [6, 3, 0], atol=1e-6)
    np.testing.assert_allclose(change[3:], 0)


def test_blend_depth_defaults_to_500_m(clean_profile):
    default = modify_surface(clean_profile, {"T": 35, "method": "blend"})
    explicit = modify_surface(clean_profile, {"T": 35, "method": "blend", "depth": 500})
    np.testing.assert_allclose(default["T"].m, explicit["T"].m)


@pytest.mark.parametrize("method", ["mix", "replace", "blend"])
def test_dewpoint_is_capped_at_the_temperature(clean_profile, method):
    modified = modify_surface(clean_profile, {"T": 25, "Td": 30, "method": method})
    assert modified["Td"][0].m == pytest.approx(25)
    assert np.all(modified["Td"].m <= modified["T"].m + 1e-9)


def test_surface_wind_needs_speed_and_direction(clean_profile):
    modified = modify_surface(clean_profile, {"ws": 20, "wd": 270})
    assert modified["u"][0].to("kt").m == pytest.approx(20)
    assert modified["v"][0].to("kt").m == pytest.approx(0, abs=1e-9)
    np.testing.assert_allclose(modified["u"][1:].m, clean_profile["u"][1:].m)

    # a speed without a direction is ignored
    unchanged = modify_surface(clean_profile, {"ws": 20})
    np.testing.assert_allclose(unchanged["u"].m, clean_profile["u"].m)


def test_modify_surface_keeps_the_input_and_its_units(clean_profile):
    original = copy.deepcopy(clean_profile)
    data = dict(clean_profile, T=clean_profile["T"].to("degF"))
    modified = modify_surface(data, {"T": 35})
    assert modified["T"].units == units.degF
    assert modified["T"][0].to("degC").m == pytest.approx(35)
    np.testing.assert_allclose(clean_profile["T"].m, original["T"].m)


@pytest.mark.parametrize(("modify_sfc", "message"), [
    ({"T": 30, "method": "mixed"}, "'method' must be 'mix', 'replace' or 'blend'"),
    ({"T": 30, "method": "blend", "depth": 0}, "'depth' must be a height in m above 0"),
    ({"T": 30, "method": "blend", "depth": -100}, "'depth' must be a height in m above 0"),
])
def test_modify_surface_rejects_bad_method_or_depth(clean_profile, modify_sfc, message):
    with pytest.raises(ValueError, match=message):
        modify_surface(clean_profile, modify_sfc)


# -----------------------------------------------------------------------------
# parcel_depth
# -----------------------------------------------------------------------------

def test_parcel_depth_lifts_only_the_lowest_levels(clean_profile, severe_params):
    # severe_params is calculated with parcel_depth=3000, as the severe figure's CAPE/CIN panel needs
    every_level = sounding_params(clean_profile, include_all_parcels=True).calc()["thermo"]
    limited = severe_params["thermo"]

    z = clean_profile["z"].m - clean_profile["z"].m[0]
    lifted = np.flatnonzero(z >= 3000)[0] + 2  # two levels past the depth, for interpolation
    for key in ("cape_profile", "cin_profile", "3cape_profile"):
        assert np.isfinite(limited[key][:lifted]).all()
        assert np.isnan(limited[key][lifted:]).all()
        np.testing.assert_allclose(limited[key][:lifted], every_level[key][:lifted])


# -----------------------------------------------------------------------------
# MISSING DATA
# -----------------------------------------------------------------------------

def test_shallow_profile_masks_deep_layer_values(clean_profile):
    # 0-4 km: no 0-6 km layer, so no Bunkers motion, and the values that need it are masked rather than raising
    shallow = {key: value[:17] if hasattr(value, "shape") else value for key, value in clean_profile.items()}
    snd = sounding_params(shallow).calc()
    assert np.ma.is_masked(snd["kinem"]["sm_u"])
    assert np.ma.is_masked(snd["kinem"]["srh_0_to_3000"])
    assert snd["thermo"]["mlcape"] > 0


def test_missing_levels_are_skipped(clean_profile):
    for key in ("T", "Td", "u", "v"):
        values = clean_profile[key].m.copy()
        values[[5, 6, 20]] = np.nan
        clean_profile[key] = values * clean_profile[key].units
    snd = sounding_params(clean_profile).calc()
    assert np.isfinite(snd["kinem"]["srh_0_to_3000"]) and snd["thermo"]["mlcape"] > 0

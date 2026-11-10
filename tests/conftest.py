"""Shared pytest fixtures for SounderPy.

The plotting fixtures are deliberately synthetic and offline.  They provide a
meteorologically reasonable deep profile without depending on NOAA/OU/PSU/IEM
servers, so plotting regressions can be tested independently of data-service
availability.
"""

from __future__ import annotations

import copy

import matplotlib

# Must be selected before importing pyplot/SounderPy plotting code.  This makes
# the suite safe on GitHub Actions, servers, WSL, and other headless systems.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from metpy.units import units



@pytest.fixture(autouse=True)
def close_matplotlib_figures():
    """Prevent figures from leaking between tests."""
    yield
    plt.close("all")


@pytest.fixture
def offline_plot_assets():
    """Mark a test as needing no network.

    The SounderPy logo is bundled with the package.  Plot tests also pass
    ``map_zoom=0`` for sounding/hodograph plots so radar and map retrieval code
    is not exercised here.  Live data-service behavior belongs in the explicitly
    enabled network tests.
    """


@pytest.fixture
def clean_profile():
    """Return a deep, physically plausible SounderPy ``clean_data`` dict.

    Heights are AGL and extend above 100 hPa.  The profile intentionally has
    temperature/dewpoint/wind structure across the full troposphere so that
    SounderPy's thermodynamic and kinematic plotting code has enough depth for
    parcel calculations, Bunkers motion, shear, SRH, lapse rates, etc.
    """
    z_m = np.arange(0.0, 18000.0 + 250.0, 250.0)
    z_km = z_m / 1000.0

    # Smooth pressure decrease to ~90 hPa near 18 km.
    p_hpa = 970.0 * np.exp(-z_m / 7600.0)

    # Warm/moist low levels with a realistic environmental lapse rate.
    t_c = 29.0 - 6.2 * z_km
    td_c = t_c - (4.0 + 0.55 * z_km)

    # Curved/deepening wind profile in knots.
    u_kt = 5.0 + 2.25 * z_km
    v_kt = 2.0 + 1.10 * z_km + 3.0 * np.sin(z_km / 3.0)

    profile = {
        "p": p_hpa * units.hPa,
        "z": z_m * units.m,
        "T": t_c * units.degC,
        "Td": td_c * units.degC,
        "u": u_kt * units.kts,
        "v": v_kt * units.kts,
        "site_info": {
            "site-id": "TEST",
            "site-name": "Synthetic SounderPy Test Profile",
            "site-lctn": "North Dakota",
            "site-latlon": [47.9, -97.0],
            "site-elv": 250,
            "source": "PYTEST SYNTHETIC PROFILE",
            "model": "TEST",
            "fcst-hour": "F00",
            "run-time": ["2024", "05", "21", "00"],
            "valid-time": ["2024", "05", "21", "00"],
        },
        "titles": {
            "top_title": "SOUNDERPY PYTEST SYNTHETIC PROFILE",
            "left_title": "VALID: 05/21/2024 00Z",
            "right_title": "TEST | 47.9, -97.0",
        },
    }
    return profile


@pytest.fixture
def second_clean_profile(clean_profile):
    """A slightly perturbed profile for composite-plot testing."""
    profile = copy.deepcopy(clean_profile)
    profile["T"] = profile["T"] + 1.0 * units.delta_degC
    profile["Td"] = profile["Td"] + 0.5 * units.delta_degC
    profile["u"] = profile["u"] + 3.0 * units.kts
    profile["titles"]["top_title"] = "SECOND PYTEST SYNTHETIC PROFILE"
    return profile


@pytest.fixture
def vad_profile():
    """Return a synthetic VAD profile matching ``vad_params`` expectations."""
    z = np.arange(0.0, 9000.0 + 100.0, 100.0)
    z_km = z / 1000.0
    return {
        # vad_params currently expects plain numeric arrays for these fields.
        "z": z,
        "u": 8.0 + 3.0 * z_km,
        "v": 2.0 + 1.5 * z_km + 2.0 * np.sin(z_km / 2.0),
        "site_info": {
            "site-id": "KTEST",
            "site-name": "Synthetic VAD",
            "site-latlon": [47.9, -97.0],
            "source": "PYTEST VAD",
            "valid-time": ["2024", "05", "21", "00"],
        },
    }


# -----------------------------------------------------------------------------
# WINTER PROFILES
# -----------------------------------------------------------------------------

# Idealized profiles that exercise each branch of the winter precip phase trace.
# Temperature control points are (m AGL, deg C); above the last point the
# profile stays isothermal.
#   fzra  warm nose aloft (~99 J/kg) over a shallow cold surface layer (~63 J/kg): freezing rain
#   pl    warm nose aloft (~54 J/kg) over a deep cold layer (~270 J/kg): ice pellets (sleet)
#   rasn  snow falling into a shallow +1.5 deg C surface layer (~12 J/kg): rain/snow mix
#   fzdz  shallow -2 to -6 deg C stratus under a warm, dry inversion: freezing drizzle
WINTER_CASES = {
    "fzra": {"name": "IDEALIZED FREEZING RAIN",
             "t": [(0, -3), (300, -4), (700, 0), (1200, 5), (1800, 0), (2500, -6), (4000, -15), (6000, -30), (10000, -55)]},
    "pl":   {"name": "IDEALIZED SLEET",
             "t": [(0, -8), (1000, -6), (1500, 0), (2000, 3), (2500, 0), (4000, -14), (6000, -30), (10000, -55)]},
    "rasn": {"name": "IDEALIZED RAIN/SNOW MIX",
             "t": [(0, 1.5), (450, 0), (1500, -5), (3000, -13), (6000, -32), (10000, -55)]},
    # shallow supercooled stratus (-2 to -6 deg C) under a warm, dry inversion: no ice from above
    "fzdz": {"name": "IDEALIZED FREEZING DRIZZLE",
             "t": [(0, -2), (1200, -6), (1400, -2), (3000, -10), (6000, -28), (10000, -55)],
             "depression": [(0, 0.3), (1200, 0.3), (1400, 15), (20000, 25)]},
}


def make_winter_profile(case):
    """Return an idealized winter ``clean_data`` dict for one of ``WINTER_CASES``."""
    import metpy.calc as mpcalc

    spec = WINTER_CASES[case]
    p = np.arange(1000, 99, -25.0) * units.hPa
    z_msl = mpcalc.pressure_to_height_std(p).to("m")
    z_agl = (z_msl - z_msl[0]).m

    zs, ts = zip(*spec["t"])
    t = np.interp(z_agl, zs, ts)
    # by default saturated (0.5 deg C depression) through 5 km so the source is high and cold (snow), dry above 6 km
    dz, dd = zip(*spec.get("depression", [(0, 0.5), (5000, 0.5), (6000, 20), (20000, 25)]))
    depression = np.interp(z_agl, dz, dd)
    td = t - depression

    # winds veer from 150 deg at 10 kt to 230 deg at 50 kt by 5 km; rising motion of 5 ub/s from 2 to 5 km
    wdir = np.interp(z_agl, [0, 5000, 20000], [150, 230, 250]) * units.deg
    wspd = np.interp(z_agl, [0, 5000, 12000], [10, 50, 70]) * units.knot
    u, v = mpcalc.wind_components(wspd, wdir)
    omega = np.interp(z_agl, [0, 1500, 2000, 5000, 6000, 20000], [0, 0, -0.5, -0.5, 0, 0]) * units("Pa/s")

    return {"p": p, "z": z_msl, "T": t * units.degC, "Td": td * units.degC, "u": u.to("knot"), "v": v.to("knot"),
            "omega": omega,
            "site_info": {"site-id": "SYN", "site-name": spec["name"], "site-lctn": "TEST",
                          "site-latlon": [45.0, -95.0], "site-elv": int(round(z_msl[0].m)), "source": "SYNTHETIC",
                          "model": "IDEALIZED", "fcst-hour": "F00", "run-time": ["2026", "01", "01", "00"],
                          "valid-time": ["2026", "01", "01", "00"]},
            "titles": {"top_title": f"SYNTHETIC TEST PROFILE | {spec['name']}",
                       "left_title": " IDEALIZED CASE FOR CHECKING THE PRECIP PHASE TRACE",
                       "right_title": "NOT REAL DATA    "}}


@pytest.fixture(params=list(WINTER_CASES))
def winter_profile(request):
    """Each idealized winter profile in turn, as ``(case name, clean_data)``."""
    return request.param, make_winter_profile(request.param)

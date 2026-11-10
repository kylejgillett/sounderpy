"""Synthetic profiles and shared assertions for the SounderPy tests.

Every profile here is built in memory, so no test needs a data service. Test modules import from here
directly (``from helpers import make_winter_profile``); ``conftest.py`` wraps the same builders as fixtures.
"""

from __future__ import annotations

import copy
from pathlib import Path

import metpy.calc as mpcalc
import numpy as np
from metpy.units import units


# -----------------------------------------------------------------------------
# SEVERE / GENERAL PROFILES
# -----------------------------------------------------------------------------

def make_clean_profile():
    """A deep, warm, moist ``clean_data`` dict (heights AGL, above 100 hPa) with curved, deepening winds."""
    z_m = np.arange(0.0, 18000.0 + 250.0, 250.0)
    z_km = z_m / 1000.0
    p_hpa = 970.0 * np.exp(-z_m / 7600.0)
    t_c = 29.0 - 6.2 * z_km
    td_c = t_c - (4.0 + 0.55 * z_km)
    u_kt = 5.0 + 2.25 * z_km
    v_kt = 2.0 + 1.10 * z_km + 3.0 * np.sin(z_km / 3.0)

    return {
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


def make_second_profile():
    """The clean profile, slightly warmer, moister and windier, for composites."""
    profile = copy.deepcopy(make_clean_profile())
    profile["T"] = profile["T"] + 1.0 * units.delta_degC
    profile["Td"] = profile["Td"] + 0.5 * units.delta_degC
    profile["u"] = profile["u"] + 3.0 * units.kts
    profile["titles"]["top_title"] = "SECOND PYTEST SYNTHETIC PROFILE"
    return profile


def make_vad_profile():
    """A VAD profile as ``vad_params`` expects it: plain z/u/v arrays and site info."""
    z = np.arange(0.0, 9000.0 + 100.0, 100.0)
    z_km = z / 1000.0
    return {
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

# Idealized profiles that exercise each branch of the winter precip phase trace. Temperature control
# points are (m AGL, deg C); above the last point the profile stays isothermal.
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
    "fzdz": {"name": "IDEALIZED FREEZING DRIZZLE",
             "t": [(0, -2), (1200, -6), (1400, -2), (3000, -10), (6000, -28), (10000, -55)],
             "depression": [(0, 0.3), (1200, 0.3), (1400, 15), (20000, 25)]},
}


def make_winter_profile(case):
    """An idealized winter ``clean_data`` dict for one of ``WINTER_CASES``, with omega."""
    spec = WINTER_CASES[case]
    p = np.arange(1000, 99, -25.0) * units.hPa
    z_msl = mpcalc.pressure_to_height_std(p).to("m")
    z_agl = (z_msl - z_msl[0]).m

    zs, ts = zip(*spec["t"])
    t = np.interp(z_agl, zs, ts)
    # by default saturated (0.5 deg C depression) through 5 km so the source is high and cold (snow), dry above 6 km
    dz, dd = zip(*spec.get("depression", [(0, 0.5), (5000, 0.5), (6000, 20), (20000, 25)]))
    td = t - np.interp(z_agl, dz, dd)

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


# -----------------------------------------------------------------------------
# ASSERTIONS
# -----------------------------------------------------------------------------

def assert_png(path: Path):
    """The file exists and is a PNG of some size."""
    assert path.exists(), f"expected plot was not created: {path}"
    assert path.stat().st_size > 1_000, f"plot file is unexpectedly small: {path.stat().st_size} bytes"
    with path.open("rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n"


def texts_on(ax):
    """Every string drawn as text on an axis."""
    return [text.get_text() for text in ax.texts]

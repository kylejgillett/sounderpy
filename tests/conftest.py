"""Shared pytest fixtures for SounderPy.

The plotting fixtures are deliberately synthetic and offline.  They provide a
meteorologically reasonable deep profile without depending on NOAA/OU/PSU/IEM
servers, so plotting regressions can be tested independently of data-service
availability.
"""

from __future__ import annotations

import base64
import copy
import io

import matplotlib

# Must be selected before importing pyplot/SounderPy plotting code.  This makes
# the suite safe on GitHub Actions, servers, WSL, and other headless systems.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest
from metpy.units import units

import sounderpy.plot as plot_module


# A valid 1x1 PNG. SounderPy currently retrieves its logo over HTTP while
# building plots; replacing urlopen keeps plot tests offline and deterministic.
_TINY_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
)


@pytest.fixture(autouse=True)
def close_matplotlib_figures():
    """Prevent figures from leaking between tests."""
    yield
    plt.close("all")


@pytest.fixture
def offline_plot_assets(monkeypatch):
    """Prevent plotting tests from downloading the SounderPy logo.

    Plot tests also pass ``map_zoom=0`` for sounding/hodograph plots so radar
    and map retrieval code is not exercised here.  Live data-service behavior
    belongs in the explicitly enabled network tests.
    """

    def fake_urlopen(*args, **kwargs):
        return io.BytesIO(_TINY_PNG)

    monkeypatch.setattr(plot_module, "urlopen", fake_urlopen)


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

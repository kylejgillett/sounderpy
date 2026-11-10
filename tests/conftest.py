"""Shared pytest setup for SounderPy: headless plotting, the --network switch, and profile fixtures."""

from __future__ import annotations

from functools import lru_cache

import matplotlib

# before pyplot or SounderPy's plotting code is imported, so the suite runs headless anywhere
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest

from helpers import WINTER_CASES, make_clean_profile, make_second_profile, make_vad_profile, make_winter_profile


# -----------------------------------------------------------------------------
# LIVE DATA TESTS
# -----------------------------------------------------------------------------

def pytest_addoption(parser):
    parser.addoption("--network", action="store_true", help="also run tests marked 'network', which contact live data services")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--network"):
        return
    skip = pytest.mark.skip(reason="contacts a live data service; run with --network")
    for item in items:
        if "network" in item.keywords:
            item.add_marker(skip)


@pytest.fixture(autouse=True)
def close_figures():
    """Keep figures from leaking between tests."""
    yield
    plt.close("all")


# -----------------------------------------------------------------------------
# PROFILES (a fresh copy for every test, so tests may change them)
# -----------------------------------------------------------------------------

@pytest.fixture
def clean_profile():
    return make_clean_profile()


@pytest.fixture
def second_clean_profile():
    return make_second_profile()


@pytest.fixture
def vad_profile():
    return make_vad_profile()


@pytest.fixture(params=list(WINTER_CASES))
def winter_profile(request):
    """Each idealized winter profile in turn, as ``(case name, clean_data)``."""
    return request.param, make_winter_profile(request.param)


# -----------------------------------------------------------------------------
# CALCULATED PARAMETERS (computed once per session; treat as read-only)
# -----------------------------------------------------------------------------

@pytest.fixture(scope="session")
def severe_params():
    """``sounding_params`` for the clean profile, calculated as the severe figure does."""
    from sounderpy.calc import sounding_params
    return sounding_params(make_clean_profile(), include_all_parcels=True, parcel_depth=3000).calc()


@pytest.fixture(scope="session")
def winter_params_for():
    """A function returning ``winter_params(...).calc()`` for a winter case, calculated once per case."""
    from sounderpy.calc_winter import winter_params

    @lru_cache(maxsize=None)
    def calc(case):
        return winter_params(make_winter_profile(case)).calc()
    return calc

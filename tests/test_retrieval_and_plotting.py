"""Tests for SounderPy data retrieval and plotting.

The suite intentionally has two retrieval layers:

1. Offline routing/API tests, which run on every normal pytest invocation and
   verify SounderPy calls the correct backend with the correct arguments.
2. Live integration tests, disabled by default, which actually contact remote
   archives and validate the returned ``clean_data`` object.

Enable live retrieval tests with::

    SOUNDERPY_RUN_NETWORK_TESTS=1 pytest -m network -v

Plot tests are offline: map/radar retrieval is disabled and the remote logo is
replaced by an in-memory PNG in ``conftest.py``.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

import sounderpy as spy
from sounderpy import sounderpy as public_api


RUN_NETWORK = os.getenv("SOUNDERPY_RUN_NETWORK_TESTS", "0").lower() in {
    "1",
    "true",
    "yes",
    "on",
}

network_only = pytest.mark.skipif(
    not RUN_NETWORK,
    reason="Set SOUNDERPY_RUN_NETWORK_TESTS=1 to run live retrieval tests.",
)


REQUIRED_PROFILE_KEYS = {"p", "z", "T", "Td", "u", "v", "site_info"}
PROFILE_ARRAY_KEYS = ("p", "z", "T", "Td", "u", "v")


def assert_clean_profile(data, *, min_levels=10):
    """Validate the common SounderPy clean-data contract.

    This intentionally checks structural/physical invariants instead of exact
    values, because remote archives can be reprocessed without representing a
    SounderPy regression.
    """
    assert isinstance(data, dict)
    assert REQUIRED_PROFILE_KEYS.issubset(data)

    lengths = []
    for key in PROFILE_ARRAY_KEYS:
        arr = data[key]
        assert hasattr(arr, "units"), f"{key!r} has no units: {type(arr)!r}"
        assert hasattr(arr, "magnitude"), f"{key!r} is not a Pint/MetPy quantity"
        assert np.ndim(arr.magnitude) == 1, f"{key!r} must be one-dimensional"
        lengths.append(len(arr))

    assert len(set(lengths)) == 1, f"profile arrays have different lengths: {lengths}"
    assert lengths[0] >= min_levels

    p = np.asarray(data["p"].to("hPa").magnitude, dtype=float)
    z = np.asarray(data["z"].to("m").magnitude, dtype=float)
    T = np.asarray(data["T"].to("degC").magnitude, dtype=float)
    Td = np.asarray(data["Td"].to("degC").magnitude, dtype=float)

    assert np.count_nonzero(np.isfinite(p)) >= min_levels
    assert np.count_nonzero(np.isfinite(z)) >= min_levels
    assert np.nanmin(p) > 0

    # Allow a few problematic archive levels, but require the overall vertical
    # coordinate ordering expected by SounderPy plotting/calculation routines.
    finite_p = p[np.isfinite(p)]
    finite_z = z[np.isfinite(z)]
    if len(finite_p) > 1:
        assert np.mean(np.diff(finite_p) <= 0) >= 0.90, "pressure is not mostly decreasing"
    if len(finite_z) > 1:
        assert np.mean(np.diff(finite_z) >= 0) >= 0.90, "height is not mostly increasing"

    # Dewpoint should not substantially exceed temperature. Permit small
    # observational/rounding noise rather than enforcing an exact inequality.
    valid_thermo = np.isfinite(T) & np.isfinite(Td)
    if np.any(valid_thermo):
        assert np.nanmax(Td[valid_thermo] - T[valid_thermo]) <= 1.0

    site_info = data["site_info"]
    assert isinstance(site_info, dict)
    assert "source" in site_info
    assert "valid-time" in site_info


# -----------------------------------------------------------------------------
# PUBLIC API / OFFLINE RETRIEVAL ROUTING TESTS
# -----------------------------------------------------------------------------


def test_retrieval_functions_are_public():
    assert callable(spy.get_model_data)
    assert callable(spy.get_obs_data)
    assert callable(spy.get_bufkit_data)
    assert callable(spy.acars_data)


def test_get_obs_data_forwards_arguments(monkeypatch):
    expected = {"backend": "obs"}
    received = {}

    def fake_fetch_obs(station, year, month, day, hour, hush, clean_it):
        received["args"] = (station, year, month, day, hour, hush, clean_it)
        return expected

    monkeypatch.setattr(public_api, "fetch_obs", fake_fetch_obs)

    result = spy.get_obs_data("DTX", "2024", "05", "21", "00", hush=True, clean_it=True)

    assert result is expected
    assert received["args"] == ("DTX", "2024", "05", "21", "00", True, True)


def test_get_bufkit_data_forwards_arguments(monkeypatch):
    expected = {"backend": "bufkit"}
    received = {}

    def fake_fetch_bufkit(model, station, fcst_hour, year, month, day, hour, hush, clean_it):
        received["args"] = (
            model,
            station,
            fcst_hour,
            year,
            month,
            day,
            hour,
            hush,
            clean_it,
        )
        return expected

    monkeypatch.setattr(public_api, "fetch_bufkit", fake_fetch_bufkit)

    result = spy.get_bufkit_data(
        "hrrr",
        "KOUN",
        3,
        "2024",
        "05",
        "21",
        "00",
        hush=True,
        clean_it=True,
    )

    assert result is expected
    assert received["args"] == (
        "hrrr",
        "KOUN",
        3,
        "2024",
        "05",
        "21",
        "00",
        True,
        True,
    )


def test_non_rap_model_routes_to_legacy_model_backend(monkeypatch):
    expected = {"backend": "model"}
    received = {}

    def fake_fetch_model(model, latlon, year, month, day, hour, dataset, box_avg_size, hush, clean_it):
        received["args"] = (
            model,
            latlon,
            year,
            month,
            day,
            hour,
            dataset,
            box_avg_size,
            hush,
            clean_it,
        )
        return expected

    monkeypatch.setattr(public_api, "fetch_model", fake_fetch_model)

    result = spy.get_model_data(
        "era5",
        [44.58, -100.82],
        "2024",
        "08",
        "28",
        "18",
        dataset="era5",
        box_avg_size=0.25,
        hush=True,
        clean_it=True,
    )

    assert result is expected
    assert received["args"] == (
        "era5",
        [44.58, -100.82],
        "2024",
        "08",
        "28",
        "18",
        "era5",
        0.25,
        True,
        True,
    )


@pytest.mark.parametrize("model", ["rap", "ruc", "rap-ruc", "RAP-RUC"])
def test_rap_ruc_models_route_to_grib_retriever(monkeypatch, model):
    expected = {"backend": "rapruc"}
    received = {}

    def fake_get_rapruc_data(latlon, year, month, day, hour, dataset, box_avg_size, hush):
        received["args"] = (latlon, year, month, day, hour, dataset, box_avg_size, hush)
        return expected

    monkeypatch.setattr(public_api, "get_rapruc_data", fake_get_rapruc_data)

    result = spy.get_model_data(
        model,
        [44.58, -100.82],
        "2024",
        "08",
        "28",
        "18",
        box_avg_size=0.25,
        hush=True,
    )

    assert result is expected
    assert received["args"] == (
        [44.58, -100.82],
        "2024",
        "08",
        "28",
        "18",
        None,
        0.25,
        True,
    )


def test_invalid_model_name_is_rejected_without_network():
    with pytest.raises(ValueError, match="not a valid model"):
        spy.get_model_data(
            "definitely-not-a-model",
            [44.58, -100.82],
            "2024",
            "08",
            "28",
            "18",
            hush=True,
        )


def test_acars_connection_object_stores_requested_time():
    conn = spy.acars_data("2024", "05", "21", "00")
    assert conn.year == "2024"
    assert conn.month == "05"
    assert conn.day == "21"
    assert conn.hour == "00"


# -----------------------------------------------------------------------------
# LIVE DATA RETRIEVAL INTEGRATION TESTS
# -----------------------------------------------------------------------------


@pytest.mark.network
@network_only
def test_live_rap_reanalysis_retrieval():
    """Known RAP/AWS-era case used to exercise download + GRIB parsing."""
    data = spy.get_model_data(
        "rap-ruc",
        [44.58, -100.82],
        "2024",
        "08",
        "28",
        "18",
        box_avg_size=0.25,
        hush=True,
    )
    assert_clean_profile(data, min_levels=30)
    assert data["site_info"]["model"].upper() in {"RAP", "RUC"}
    assert data["site_info"]["valid-time"][:4] == ["2024", "08", "28", "18"]
    assert "archive" in data["site_info"]


@pytest.mark.network
@network_only
def test_live_raob_retrieval():
    """Exercise the Siphon/Wyoming observed-sounding path."""
    data = spy.get_obs_data("DTX", "2024", "05", "21", "00", hush=True)
    assert_clean_profile(data, min_levels=15)
    assert data["site_info"]["source"] == "RAOB OBSERVED PROFILE"


@pytest.mark.network
@network_only
def test_live_bufkit_retrieval():
    """Exercise a documented archived GFS BUFKIT sounding."""
    data = spy.get_bufkit_data(
        "gfs",
        "KMOP",
        6,
        "2023",
        "08",
        "05",
        "12",
        hush=True,
    )

    assert_clean_profile(data)
    assert data["site_info"]["source"] == "BUFKIT FORECAST PROFILE"
    assert data["site_info"]["model"] == "GFS"


# -----------------------------------------------------------------------------
# PLOTTING TESTS
# -----------------------------------------------------------------------------


def assert_nonempty_png(path: Path):
    assert path.exists(), f"expected plot was not created: {path}"
    assert path.is_file()
    assert path.stat().st_size > 1_000, f"plot file is unexpectedly small: {path.stat().st_size} bytes"
    with path.open("rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n"


@pytest.mark.slow
def test_build_sounding_saves_rendered_png(clean_profile, offline_plot_assets, tmp_path):
    output = tmp_path / "sounding.png"

    result = spy.build_sounding(
        clean_profile,
        map_zoom=0,      # keeps radar/map retrieval offline
        dpi=50,
        save=True,
        filename=str(output),
    )

    # Current public wrapper saves/shows but does not return the pyplot object.
    assert result is None
    assert_nonempty_png(output)


@pytest.mark.slow
def test_build_sounding_dark_colorblind_variant(clean_profile, offline_plot_assets, tmp_path):
    output = tmp_path / "sounding_dark_cb.png"

    spy.build_sounding(
        clean_profile,
        color_blind=True,
        dark_mode=True,
        map_zoom=0,
        dpi=50,
        save=True,
        filename=str(output),
    )

    assert_nonempty_png(output)


@pytest.mark.slow
def test_build_hodograph_saves_rendered_png(clean_profile, offline_plot_assets, tmp_path):
    output = tmp_path / "hodograph.png"

    result = spy.build_hodograph(
        clean_profile,
        map_zoom=0,
        dpi=50,
        save=True,
        filename=str(output),
    )

    assert result is None
    assert_nonempty_png(output)


@pytest.mark.slow
def test_build_hodograph_storm_relative(clean_profile, offline_plot_assets, tmp_path):
    output = tmp_path / "hodograph_sr.png"

    spy.build_hodograph(
        clean_profile,
        storm_motion="right_moving",
        sr_hodo=True,
        map_zoom=0,
        dpi=50,
        save=True,
        filename=str(output),
    )

    assert_nonempty_png(output)


@pytest.mark.slow
def test_build_composite_saves_rendered_png(
    clean_profile,
    second_clean_profile,
    offline_plot_assets,
    tmp_path,
):
    output = tmp_path / "composite.png"

    result = spy.build_composite(
        [clean_profile, second_clean_profile],
        shade_between=True,
        save=True,
        filename=str(output),
    )

    assert result is None
    assert_nonempty_png(output)


@pytest.mark.slow
def test_build_vad_hodograph_saves_rendered_png(vad_profile, offline_plot_assets, tmp_path):
    output = tmp_path / "vad_hodograph.png"

    result = spy.build_vad_hodograph(
        vad_profile,
        save=True,
        filename=str(output),
    )

    assert result is None
    assert_nonempty_png(output)

"""Data retrieval: offline routing of the public get_* functions, surface observations, and live data tests.

The routing tests mock each backend and check SounderPy calls it with the right arguments. Tests marked
``network`` contact the real archives and only run with ``pytest --network``.
"""

from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pytest

import sounderpy as spy
from sounderpy import sounderpy as public_api
from sounderpy.data import sfc_obs

from helpers import make_winter_profile


def assert_clean_profile(data, min_levels=10):
    """Structural and physical checks on a retrieved ``clean_data`` dict (not exact values, which archives can change)."""
    assert {"p", "z", "T", "Td", "u", "v", "site_info"} <= set(data)

    lengths = set()
    for key in ("p", "z", "T", "Td", "u", "v"):
        assert hasattr(data[key], "units"), f"{key!r} has no units"
        assert np.ndim(data[key].m) == 1
        lengths.add(len(data[key]))
    assert len(lengths) == 1 and lengths.pop() >= min_levels

    p, z = data["p"].to("hPa").m, data["z"].to("m").m
    T, Td = data["T"].to("degC").m, data["Td"].to("degC").m
    assert np.count_nonzero(np.isfinite(p)) >= min_levels and np.nanmin(p) > 0
    # allow a few problem levels, but the profile should mostly run upward
    assert np.mean(np.diff(p[np.isfinite(p)]) <= 0) >= 0.9
    assert np.mean(np.diff(z[np.isfinite(z)]) >= 0) >= 0.9
    ok = np.isfinite(T) & np.isfinite(Td)
    assert np.max(Td[ok] - T[ok]) <= 1.0
    assert {"source", "valid-time"} <= set(data["site_info"])


def capture(monkeypatch, name, result):
    """Replace a backend on the public API module; return the list of argument tuples it receives."""
    calls = []
    monkeypatch.setattr(public_api, name, lambda *args: calls.append(args) or result)
    return calls


# -----------------------------------------------------------------------------
# OFFLINE ROUTING
# -----------------------------------------------------------------------------

def test_get_obs_data_forwards_arguments(monkeypatch):
    calls = capture(monkeypatch, "fetch_obs", {"backend": "obs"})
    assert spy.get_obs_data("DTX", "2024", "05", "21", "00", hush=True, clean_it=True) == {"backend": "obs"}
    assert calls == [("DTX", "2024", "05", "21", "00", True, True)]


def test_get_bufkit_data_forwards_arguments(monkeypatch):
    calls = capture(monkeypatch, "fetch_bufkit", {"backend": "bufkit"})
    result = spy.get_bufkit_data("hrrr", "KOUN", 3, "2024", "05", "21", "00", hush=True, clean_it=True, all_fcst_hours=True)
    assert result == {"backend": "bufkit"}
    assert calls == [("hrrr", "KOUN", 3, "2024", "05", "21", "00", True, True, True)]


def test_non_rap_model_routes_to_model_backend(monkeypatch):
    calls = capture(monkeypatch, "fetch_model", {"backend": "model"})
    result = spy.get_model_data("era5", [44.58, -100.82], "2024", "08", "28", "18", dataset="era5", box_avg_size=0.25,
                                hush=True, clean_it=True)
    assert result == {"backend": "model"}
    assert calls == [("era5", [44.58, -100.82], "2024", "08", "28", "18", "era5", 0.25, True, True)]


@pytest.mark.parametrize("model", ["rap", "ruc", "rap-ruc", "RAP-RUC"])
def test_rap_ruc_models_route_to_grib_retriever(monkeypatch, model):
    calls = capture(monkeypatch, "get_rapruc_data", {"backend": "rapruc"})
    result = spy.get_model_data(model, [44.58, -100.82], "2024", "08", "28", "18", box_avg_size=0.25, hush=True, clean_it=False)
    assert result == {"backend": "rapruc"}
    assert calls == [([44.58, -100.82], "2024", "08", "28", "18", None, 0.25, True, False)]


def test_invalid_model_name_is_rejected_without_network():
    with pytest.raises(ValueError, match="not a valid model"):
        spy.get_model_data("definitely-not-a-model", [44.58, -100.82], "2024", "08", "28", "18", hush=True)


def test_acars_connection_object_stores_requested_time():
    conn = spy.acars_data("2024", "05", "21", "00")
    assert (conn.year, conn.month, conn.day, conn.hour) == ("2024", "05", "21", "00")


# -----------------------------------------------------------------------------
# SURFACE OBSERVATIONS (winter figure)
# -----------------------------------------------------------------------------

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


# -----------------------------------------------------------------------------
# RAP/RUC GRIB PARSING FROM SEVERAL THREADS (small GRIB files written at test time)
# -----------------------------------------------------------------------------

GRIB_LEVELS = {1000: (100, 298), 850: (1500, 290), 700: (3000, 280), 500: (5600, 262), 300: (9200, 232), 200: (11800, 218)}


def write_rapruc_grib(directory, name, offset):
    """Write a small RAP-style pressure file and native file (5x5 grid around 44N, 100W)."""
    eccodes = pytest.importorskip("eccodes")

    def write(f, short, level_type, level, value):
        h = eccodes.codes_grib_new_from_samples("regular_ll_sfc_grib2")
        for key, val in {"Ni": 5, "Nj": 5, "latitudeOfFirstGridPointInDegrees": 46, "longitudeOfFirstGridPointInDegrees": 258,
                         "latitudeOfLastGridPointInDegrees": 42, "longitudeOfLastGridPointInDegrees": 262,
                         "iDirectionIncrementInDegrees": 1, "jDirectionIncrementInDegrees": 1}.items():
            eccodes.codes_set(h, key, val)
        eccodes.codes_set(h, "typeOfLevel", level_type)
        if level is not None:
            eccodes.codes_set(h, "level", level)
        eccodes.codes_set(h, "shortName", short)
        eccodes.codes_set_values(h, np.full(25, float(value)))
        eccodes.codes_write(h, f)
        eccodes.codes_release(h)

    pressure_file, native_file = directory / f"{name}_prs.grib2", directory / f"{name}_nat.grib2"
    with open(pressure_file, "wb") as f:
        for p, (gh, t) in GRIB_LEVELS.items():
            write(f, "gh", "isobaricInhPa", p, gh)
            write(f, "t", "isobaricInhPa", p, t + offset)
            write(f, "r", "isobaricInhPa", p, 60)
            write(f, "w", "isobaricInhPa", p, -0.1)
            write(f, "u", "isobaricInhPa", p, p / 50 + offset)
            write(f, "v", "isobaricInhPa", p, -p / 100 - offset)
    with open(native_file, "wb") as f:
        write(f, "2t", "heightAboveGround", 2, 300 + offset)
        write(f, "2r", "heightAboveGround", 2, 55)
        write(f, "10u", "heightAboveGround", 10, 3 + offset)
        write(f, "10v", "heightAboveGround", 10, 4 - offset)
        write(f, "sp", "surface", None, 98000)
        write(f, "orog", "surface", None, 300)
    return {"pressure_file": str(pressure_file), "native_file": str(native_file)}


@pytest.fixture(scope="module")
def rapruc_grib_cases(tmp_path_factory):
    pytest.importorskip("cfgrib")
    directory = tmp_path_factory.mktemp("rapruc_grib")
    return [write_rapruc_grib(directory, "a", 0.0), write_rapruc_grib(directory, "b", -5.0)]


def run_together(func, cases, repeats=5):
    """Call ``func`` on every case at the same moment from separate threads, ``repeats`` times over."""
    import threading
    from concurrent.futures import ThreadPoolExecutor

    results = []
    with ThreadPoolExecutor(max_workers=len(cases)) as pool:
        for _ in range(repeats):
            barrier = threading.Barrier(len(cases))

            def call(case):
                barrier.wait()
                return func(case)

            results.append(list(pool.map(call, cases)))
    return results


def test_rapruc_grib_parsing_is_thread_safe(rapruc_grib_cases):
    from sounderpy.data.rapruc_grib_parser import parse_rapruc_grib

    def parse(case):
        return parse_rapruc_grib(latlon=[44.0, -100.0], valid_time="2022-05-31 01:00", **case)

    expected = [parse(case) for case in rapruc_grib_cases]
    # the two files differ, so a profile mixed up with the other thread's data would not match
    assert not np.allclose(expected[0]["T"].m, expected[1]["T"].m, equal_nan=True)

    for run in run_together(parse, rapruc_grib_cases):
        for data, reference in zip(run, expected):
            assert_clean_profile(data, min_levels=30)
            for key in ("p", "z", "T", "Td", "u", "v"):
                assert np.all(np.isfinite(data[key].m)), f"{key!r} has missing values"
                np.testing.assert_allclose(data[key].m, reference[key].m)


def test_rapruc_grib_subsetting_is_thread_safe(rapruc_grib_cases):
    import xarray as xr
    from sounderpy.data.rapruc_grib_parser import subset_rapruc_grib

    def subset(case):
        return subset_rapruc_grib(latlon=[44.0, -100.0], **case)

    expected = [subset(case) for case in rapruc_grib_cases]

    for run in run_together(subset, rapruc_grib_cases):
        for raw, reference in zip(run, expected):
            assert set(raw) == {"pressure", "2m", "10m", "surface"}
            assert {"gh", "t", "r", "w", "u", "v"} <= set(raw["pressure"].data_vars)
            for key in raw:
                # values, coordinates, and dims only: the GRIB 'history' attribute holds a timestamp
                xr.testing.assert_equal(raw[key], reference[key])


# -----------------------------------------------------------------------------
# LIVE DATA (pytest --network)
# -----------------------------------------------------------------------------

@pytest.mark.network
def test_live_rap_reanalysis_retrieval():
    data = spy.get_model_data("rap-ruc", [44.58, -100.82], "2024", "08", "28", "18", box_avg_size=0.25, hush=True)
    assert_clean_profile(data, min_levels=30)
    assert data["site_info"]["model"].upper() in {"RAP", "RUC"}
    assert data["site_info"]["valid-time"][:4] == ["2024", "08", "28", "18"]
    assert "archive" in data["site_info"]


@pytest.mark.network
def test_live_raob_retrieval():
    data = spy.get_obs_data("DTX", "2024", "05", "21", "00", hush=True)
    assert_clean_profile(data, min_levels=15)
    assert data["site_info"]["source"] == "RAOB OBSERVED PROFILE"


@pytest.mark.network
def test_live_bufkit_retrieval():
    data = spy.get_bufkit_data("gfs", "KMOP", 6, "2023", "08", "05", "12", hush=True)
    assert_clean_profile(data)
    assert data["site_info"]["source"] == "BUFKIT FORECAST PROFILE"
    assert data["site_info"]["model"] == "GFS"


@pytest.mark.network
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

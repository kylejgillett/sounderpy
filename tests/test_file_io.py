"""Tests for the current SounderPy file I/O API.

This replaces the stale test that referenced ``read_sharppy_file``.  The public
reader in the current package is ``from_sharppy``.
"""

from __future__ import annotations

import numpy as np
import pytest

import sounderpy as spy


REQUIRED_KEYS = {"p", "z", "T", "Td", "u", "v", "wd", "ws", "site_info"}


def test_sharppy_round_trip(clean_profile, tmp_path):
    path = tmp_path / "round_trip.snd"

    spy.to_file("sharppy", clean_profile, filename=str(path))
    assert path.exists()
    assert path.stat().st_size > 0

    loaded = spy.from_sharppy(path)

    assert REQUIRED_KEYS.issubset(loaded)
    assert len(loaded["p"]) == len(clean_profile["p"])
    assert loaded["site_info"]["site-id"] == "TEST"
    assert loaded["site_info"]["valid-time"] == [2024, 5, 21, 0]

    np.testing.assert_allclose(
        loaded["p"].to("hPa").magnitude,
        clean_profile["p"].to("hPa").magnitude,
        rtol=0,
        atol=1e-5,
    )
    np.testing.assert_allclose(
        loaded["T"].to("degC").magnitude,
        clean_profile["T"].to("degC").magnitude,
        rtol=0,
        atol=1e-5,
    )


def test_csv_export_contains_expected_columns(clean_profile, tmp_path):
    path = tmp_path / "profile.csv"
    spy.to_file("csv", clean_profile, filename=str(path))

    first_line = path.read_text(encoding="utf-8").splitlines()[0]
    assert first_line == "p,z,T,Td,u,v"


def test_from_sharppy_rejects_empty_file(tmp_path):
    path = tmp_path / "empty.snd"
    path.write_text("", encoding="utf-8")

    with pytest.raises(ValueError, match="empty"):
        spy.from_sharppy(path)


def test_from_sharppy_rejects_non_sharppy_file(tmp_path):
    path = tmp_path / "not_sharppy.txt"
    path.write_text("this is not a sounding\n", encoding="utf-8")

    with pytest.raises(ValueError, match="does not appear to be a SHARPpy"):
        spy.from_sharppy(path)


# -----------------------------------------------------------------------------
# RADIOSONDE FILE READERS (GRAW, iMet, Windsond)
# trimmed UND launches in test_data/: the first ~300 rows and the end of each flight
# -----------------------------------------------------------------------------

from pathlib import Path

TEST_DATA = Path(__file__).resolve().parents[1] / "test_data"
GRAW_FILE = TEST_DATA / "20260125_2100_graw_sample.txt"
IMET_FILE = TEST_DATA / "20260116_1830_imet_sample.csv"
WINDSOND_FILE = TEST_DATA / "20260118_1550_windsond_sample.csv"


def assert_ascent(data):
    """The reader's output has every key, and pressure falls with height all the way up."""
    assert REQUIRED_KEYS.issubset(data)
    p, z = data["p"].to("hPa").m, data["z"].to("m").m
    assert len(p) > 100
    assert np.all(np.diff(p) < 0)
    assert np.all(np.isfinite(p)) and np.all(np.isfinite(z))
    assert np.all(data["Td"].m <= data["T"].m + 1e-6)


def test_from_graw_reads_the_ascent(capsys):
    data = spy.from_graw(GRAW_FILE)
    assert_ascent(data)
    # the first row is the ground station's surface observation
    assert data["p"][0].m == pytest.approx(1001.1)
    assert data["T"][0].m == pytest.approx(-23.0)
    assert data["Td"][0].m == pytest.approx(-28.4)
    assert data["z"][0].m == pytest.approx(258)
    # the last rows reach the burst
    assert data["p"][-1].m == pytest.approx(80.8)
    # 4.0 m/s from 300 degrees
    assert data["ws"][0].m == pytest.approx(4.0 * 1.943844, rel=1e-3)
    assert data["wd"][0].m == pytest.approx(300)
    assert data["site_info"]["valid-time"] == [2026, 1, 25, 20]
    assert data["site_info"]["site-elv"] == 258
    assert data["site_info"]["site-latlon"] == [47.92, -97.09]
    assert "Still needed: site ID, name, location" in capsys.readouterr().out


def test_from_graw_date_argument_and_missing_date(tmp_path):
    renamed = tmp_path / "launch.txt"
    renamed.write_bytes(GRAW_FILE.read_bytes())
    with pytest.raises(ValueError, match="don't hold the launch date"):
        spy.from_graw(renamed)
    for date in ("2026-01-25", "20260125"):
        assert spy.from_graw(renamed, date=date)["site_info"]["valid-time"] == [2026, 1, 25, 20]


def test_from_graw_drops_levels_where_pressure_does_not_fall():
    raw = [line.split("\t") for line in GRAW_FILE.read_text(encoding="latin-1").splitlines()[1:]]
    pressures = [float(row[2]) for row in raw]
    repeats = sum(1 for a, b in zip(pressures, pressures[1:]) if b >= a)
    assert repeats > 0  # the sample has some to drop
    assert len(spy.from_graw(GRAW_FILE)["p"]) < len(raw)


def test_from_imet_drops_prelaunch_and_descent():
    data = spy.from_imet(IMET_FILE)
    assert_ascent(data)
    # the sample has 40 indoor pre-launch rows (about 21 C) before the launch; the profile starts outside
    assert data["T"][0].m < -10
    # and ends at the top of the flight, not after the burst
    assert data["z"][-1].m == pytest.approx(15844.5, abs=1)
    assert data["site_info"]["site-name"] == "Grand Forks, ND"
    assert data["site_info"]["site-latlon"] == [47.92, -97.09]
    assert data["site_info"]["site-elv"] == 256
    assert data["site_info"]["valid-time"][:3] == [2026, 1, 16]


def test_from_imet_missing_winds_are_nan():
    data = spy.from_imet(IMET_FILE)
    # GPS winds start a little after launch; the 99999 placeholders become NaN, never a 99999 kt wind
    assert np.isnan(data["u"].m[0])
    assert np.nanmax(data["ws"].m) < 200


def test_from_windsond_heights_and_elevation():
    agl = spy.from_windsond(WINDSOND_FILE)
    assert_ascent(agl)
    assert agl["z"][0].m == pytest.approx(20)
    assert agl["site_info"]["site-elv"] is None
    assert agl["site_info"]["site-latlon"] == [47.93, -97.12]
    assert agl["site_info"]["valid-time"] == [2026, 1, 18, 15]

    msl = spy.from_windsond(WINDSOND_FILE, elevation=257)
    np.testing.assert_allclose(msl["z"].m, agl["z"].m + 257)
    assert msl["site_info"]["site-elv"] == 257


@pytest.mark.parametrize("reader, message", [(spy.from_graw, "empty"), (spy.from_imet, "empty"), (spy.from_windsond, "empty")])
def test_radiosonde_readers_reject_empty_files(reader, message, tmp_path):
    path = tmp_path / "20260101_0000.txt"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        reader(path)


@pytest.mark.parametrize("reader, name", [(spy.from_graw, "GRAW"), (spy.from_imet, "iMet"), (spy.from_windsond, "Windsond")])
def test_radiosonde_readers_reject_other_files(reader, name, tmp_path):
    path = tmp_path / "20260101_0000.txt"
    path.write_text("this is not a sounding\nnor is this\n", encoding="utf-8")
    with pytest.raises(ValueError, match=f"does not appear to be a.* {name}"):
        reader(path)


@pytest.mark.slow
@pytest.mark.parametrize("reader, path, kwargs", [(spy.from_graw, GRAW_FILE, {}), (spy.from_imet, IMET_FILE, {}),
                                                   (spy.from_windsond, WINDSOND_FILE, {"elevation": 257})])
def test_radiosonde_files_plot(reader, path, kwargs, tmp_path):
    from helpers import assert_png
    data = reader(path, **kwargs)
    data["site_info"].update({"site-id": "UND", "site-name": "UND", "site-lctn": "ND"})
    data["titles"] = {"top_title": "RADIOSONDE TEST", "left_title": "TEST", "right_title": "UND"}
    out = tmp_path / "radiosonde.png"
    spy.build_sounding(data, mode="winter", save=True, filename=str(out), map_zoom=0, dpi=50)
    assert_png(out)


@pytest.mark.parametrize("reader, path, kwargs, raw_speed_ms, raw_direction", [
    (spy.from_graw, GRAW_FILE, {}, 4.0, 300.0),                    # first row: 4.0 m/s from 300 degrees
    (spy.from_windsond, WINDSOND_FILE, {}, 11.81, 353.1),          # first row: 11.81 m/s from 353.1 degrees
])
def test_radiosonde_winds_are_knots_from_the_direction(reader, path, kwargs, raw_speed_ms, raw_direction):
    data = reader(path, **kwargs)
    for key in ("ws", "u", "v"):
        assert str(data[key].units) == "knot", key
    assert str(data["wd"].units) == "degree"
    speed_kt = raw_speed_ms * 1.943844
    assert data["ws"][0].m == pytest.approx(speed_kt, rel=1e-4)
    # meteorological convention: the direction the wind blows from
    assert data["u"][0].m == pytest.approx(-speed_kt * np.sin(np.radians(raw_direction)), abs=1e-3)
    assert data["v"][0].m == pytest.approx(-speed_kt * np.cos(np.radians(raw_direction)), abs=1e-3)


def test_imet_winds_are_knots():
    data = spy.from_imet(IMET_FILE)
    first = int(np.where(np.isfinite(data["ws"].m))[0][0])
    # the sample's first GPS wind after launch (18:37:38): 18.62 m/s from 333.88 degrees
    assert str(data["ws"].units) == "knot"
    assert data["ws"][first].m == pytest.approx(18.62 * 1.943844, rel=1e-3)
    assert data["wd"][first].m == pytest.approx(333.88, abs=0.01)


# -----------------------------------------------------------------------------
# CM1 AND WRF READERS (from_cm1 / from_wrf, formerly make_cm1_profile / make_wrf_profile)
# -----------------------------------------------------------------------------

CM1_FILE = TEST_DATA / "input_sounding"
CM1_META = {"latlon": [45.10, -100.89], "elev": 358, "top_title": "CM1 TEST"}


def test_from_cm1_reads_input_sounding():
    data = spy.from_cm1(str(CM1_FILE), CM1_META)
    assert {"p", "z", "T", "Td", "u", "v", "site_info", "titles"}.issubset(data)
    assert np.all(np.diff(data["p"].to("hPa").m) < 0)


def test_make_cm1_profile_still_works_with_a_warning():
    with pytest.warns(FutureWarning, match=r"make_cm1_profile\(\) is now from_cm1\(\)"):
        old = spy.make_cm1_profile(str(CM1_FILE), CM1_META)
    new = spy.from_cm1(str(CM1_FILE), CM1_META)
    for key in ("p", "z", "T", "Td", "u", "v"):
        np.testing.assert_allclose(old[key].m, new[key].m)


def test_make_wrf_profile_routes_to_from_wrf(monkeypatch):
    # no WRF output in test_data: check the old name warns and passes every argument through unchanged
    import sounderpy.data.wrf_utils as wrf_utils
    calls = []
    monkeypatch.setattr(wrf_utils, "from_wrf", lambda *args, **kwargs: calls.append((args, kwargs)) or "profile")
    with pytest.warns(FutureWarning, match=r"make_wrf_profile\(\) is now from_wrf\(\)"):
        result = spy.make_wrf_profile("ds", [48.9, -100.8], run_name="RUN", time_idx=2)
    assert result == "profile"
    assert calls == [(("ds", [48.9, -100.8]), {"run_name": "RUN", "time_idx": 2})]
    assert spy.from_wrf is not spy.make_wrf_profile

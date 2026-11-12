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
# small synthetic files in each format, written to tmp_path: a smooth ascent, 10 hPa between levels
# -----------------------------------------------------------------------------

import warnings
from datetime import datetime, timedelta
from pathlib import Path

from sounderpy.utils import SounderPyNotice, SounderPyWarning

TEST_DATA = Path(__file__).resolve().parents[1] / "test_data"
UTF8_BOM = b"\xef\xbb\xbf"


def standard_ascent(n=80, p0=1000.0, dp=10.0):
    """Pressure (hPa), height (m above sea level), and temperature (degC) for n levels dp apart."""
    p = p0 - dp * np.arange(n)
    z = 44330.8 * (1 - (p / 1013.25) ** 0.190263)
    return p, z, 5 - 6.5 * z / 1000


GRAW_HEADER = ["Time [sec]", "Time (UTC)", "P [hPa]", "T [°C]", "Hu [%]", "Ws [m/s]", "Wd [°]", "Alt [m]", "Geopot [m]",
               "Lon [°]", "Lat [°]", "Dew [°C]", "RI", "Elevation [°]", "Azimuth [°]", "Range [m]", "D [kg/m3]"]


def graw_rows(n=80, launch="20:26:35"):
    """GRAW rows: the ground station's surface observation (time 00:00:00), then one level every 10 s from `launch`."""
    p, z, t = standard_ascent(n)
    start = datetime.strptime(launch, "%H:%M:%S")
    rows = []
    for i in range(n):
        time = "00:00:00" if i == 0 else (start + timedelta(seconds=10 * (i - 1))).strftime("%H:%M:%S")
        rows.append([str(10 * i), time, f"{p[i]:.1f}", f"{t[i]:.1f}", "60", "4.0", "300", f"{z[i]:.0f}", f"{z[i]:.0f}",
                     "-97.087270", "47.921420", f"{t[i] - 5:.1f}", "300.0", "10.0", "200", "100", "1.2"])
    return rows


def write_graw(path, rows=None, bom=False, newline="\r\n"):
    text = newline.join("\t".join(row) for row in [GRAW_HEADER] + (graw_rows() if rows is None else rows)) + newline
    path.write_bytes((UTF8_BOM if bom else b"") + text.encode("latin-1"))
    return path


IMET_HEADER = ["date [y-m-d GMT]", "time [h:m:s GMT]", "altitude (from iMet PTU) [km]", "iMet pressure [mb]",
               "iMet air temperature (corrected) [deg C]", "iMet humidity [RH %]", "iMet ascent rate [m/s]",
               "GPS wind speed [m/s]", "GPS wind direction [deg]"]
IMET_META = ["NOAA/GMD balloon flight data file", "", "Station Name:, Grand Forks, ND", "Station Latitude:, 47.92131",
             "Station Longitude:, -97.08764", "Station Altitude:, 256.00 [m]", "Station Pressure:, 99999 [mb]", ""]
IMET_LAUNCH = datetime(2026, 1, 16, 18, 37, 33)


def imet_rows(n=80, prelaunch=5, descent=5, no_wind=0):
    """iMet rows: the sonde sitting indoors (21 C, not rising), the ascent from IMET_LAUNCH, then the fall after the burst.
    The first `no_wind` levels of the ascent have the 99999 placeholder winds."""
    p, z, t = standard_ascent(n, p0=975.0)

    def row(when, z, p, t, ascent, ws="12.0", wd="330.0"):
        return [when.strftime("%Y-%m-%d"), when.strftime("%H:%M:%S"), f"{z / 1000:.4f}", f"{p:.2f}", f"{t:.2f}", "60.0",
                f"{ascent:.2f}", ws, wd]

    rows = [row(IMET_LAUNCH - timedelta(seconds=prelaunch - i), z[0], p[0], 21.0, 0.0, "99999", "99999")
            for i in range(prelaunch)]
    rows += [row(IMET_LAUNCH + timedelta(seconds=10 * i), z[i], p[i], t[i], 5.0,
                 *(("99999", "99999") if i < no_wind else ())) for i in range(n)]
    rows += [row(IMET_LAUNCH + timedelta(seconds=10 * (n + i)), z[-1] - 100 * (i + 1), p[-1] + 5 * (i + 1), t[-1], -10.0)
             for i in range(descent)]
    return rows


def write_imet(path, rows=None, bom=False, newline="\r\n"):
    lines = IMET_META + [", ".join(IMET_HEADER)] + [", ".join(row) for row in (imet_rows() if rows is None else rows)]
    path.write_bytes((UTF8_BOM if bom else b"") + (newline.join(lines) + newline).encode("utf-8"))
    return path


WINDSOND_COMMENT = "# Radiation correction v2.4. Params: lat=47.9293, lon=-97.1152, utc_time=2026-01-18 15:56"
WINDSOND_HEADER = ("Height (m AGL), Pressure (mb), Temperature (C), Relative humidity (%), Wind speed (m/s), "
                   "Wind direction (true deg)")


def windsond_rows(n=80):
    """Windsond rows, heights above ground from 20 m."""
    p, z, t = standard_ascent(n, p0=980.0)
    return [[f"{zi - z[0] + 20:.0f}", f"{pi:.2f}", f"{ti:.2f}", "50.0", "11.81", "353.1"] for pi, zi, ti in zip(p, z, t)]


def write_windsond(path, rows=None, bom=False, newline="\r\n"):
    lines = [WINDSOND_COMMENT, WINDSOND_HEADER] + [", ".join(row) for row in (windsond_rows() if rows is None else rows)]
    path.write_bytes((UTF8_BOM if bom else b"") + (newline.join(lines) + newline).encode("utf-8"))
    return path


def read_quietly(reader, path, **kwargs):
    """A reader's clean_data, with its warnings recorded rather than shown: (data, [warning, ...])."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        data = reader(path, **kwargs)
    return data, caught


def data_warnings(caught):
    return [str(w.message) for w in caught if issubclass(w.category, SounderPyWarning)]


# one entry per format: reader, file writer, row maker, the reader's extra arguments, and the file's name
FORMATS = {
    "graw": (spy.from_graw, write_graw, graw_rows, {}, "20260125_2100.txt"),
    "imet": (spy.from_imet, write_imet, imet_rows, {}, "20260116_1830.csv"),
    "windsond": (spy.from_windsond, write_windsond, windsond_rows, {"elevation": 257}, "20260118_1550.csv"),
}


def read_format(fmt, tmp_path, rows=None, **write_kwargs):
    reader, write, _, kwargs, name = FORMATS[fmt]
    return read_quietly(reader, write(tmp_path / name, rows, **write_kwargs), **kwargs)


def assert_ascent(data):
    """The reader's output has exactly the clean_data keys, and pressure falls with height all the way up."""
    assert set(data) == REQUIRED_KEYS
    p, z = data["p"].to("hPa").m, data["z"].to("m").m
    assert np.all(np.diff(p) < 0)
    assert np.all(np.isfinite(p)) and np.all(np.isfinite(z)) and np.all(np.isfinite(data["T"].m))
    assert np.all(data["Td"].m <= data["T"].m + 1e-6)
    assert {"launch-time", "data-gaps", "missing-sfc-wind"} <= set(data["site_info"])


def test_from_graw_reads_the_ascent(tmp_path):
    data, caught = read_format("graw", tmp_path)
    assert_ascent(data)
    assert len(data["p"]) == 80
    # the first row is the ground station's surface observation
    _, _, t = standard_ascent()
    assert data["p"][0].m == pytest.approx(1000.0)
    assert data["T"][0].m == pytest.approx(round(t[0], 1))
    assert data["Td"][0].m == pytest.approx(round(t[0] - 5, 1))
    # the last row is the burst
    assert data["p"][-1].m == pytest.approx(210.0)
    assert data["wd"][0].m == pytest.approx(300)
    site = data["site_info"]
    assert site["valid-time"] == [2026, 1, 25, 20]
    assert site["launch-time"] == "2026-01-25T20:26:35Z"
    assert site["site-elv"] == round(standard_ascent()[1][0])
    assert site["site-latlon"] == [47.92, -97.09]
    assert site["data-gaps"] == [] and site["missing-sfc-wind"] is False
    assert data_warnings(caught) == []


def test_from_graw_date_argument_and_missing_date(tmp_path):
    path = write_graw(tmp_path / "launch.txt")
    with pytest.raises(ValueError, match="don't hold the launch date"):
        spy.from_graw(path)
    for date in ("2026-01-25", "20260125", datetime(2026, 1, 25).date(), datetime(2026, 1, 25, 21)):
        assert read_quietly(spy.from_graw, path, date=date)[0]["site_info"]["launch-time"] == "2026-01-25T20:26:35Z"
    with pytest.raises(ValueError, match="isn't a date"):
        spy.from_graw(path, date="25 Jan")


@pytest.mark.parametrize("name, launch, expected, moved", [
    ("20260126_0000.txt", "23:40:00", "2026-01-25T23:40:00Z", True),    # the 00Z launch went up before midnight
    ("20260125_2300.txt", "00:20:00", "2026-01-26T00:20:00Z", True),    # a late launch went up after midnight
    ("20260125_2100.txt", "20:26:35", "2026-01-25T20:26:35Z", False),
    ("20260126.txt", "23:40:00", "2026-01-26T23:40:00Z", False),         # no time in the name: the date is the launch date
])
def test_from_graw_launch_near_midnight(tmp_path, name, launch, expected, moved):
    data, caught = read_quietly(spy.from_graw, write_graw(tmp_path / name, graw_rows(launch=launch)))
    assert data["site_info"]["launch-time"] == expected
    assert data["site_info"]["valid-time"] == [int(expected[:4]), int(expected[5:7]), int(expected[8:10]), int(expected[11:13])]
    assert any("more than 12 hours" in message for message in data_warnings(caught)) == moved


def test_from_graw_date_argument_is_not_moved(tmp_path):
    path = write_graw(tmp_path / "20260126_0000.txt", graw_rows(launch="23:40:00"))
    data, caught = read_quietly(spy.from_graw, path, date="2026-01-26")
    assert data["site_info"]["launch-time"] == "2026-01-26T23:40:00Z"
    assert data_warnings(caught) == []


def test_from_graw_drops_levels_where_pressure_does_not_fall(tmp_path):
    rows = graw_rows()
    rows[10:10] = [rows[9][:], rows[9][:]]     # the balloon bobbing: the same pressure twice more
    data, _ = read_format("graw", tmp_path, rows)
    assert_ascent(data)
    assert len(data["p"]) == 80


def test_from_imet_drops_prelaunch_and_descent(tmp_path):
    data, _ = read_format("imet", tmp_path)
    assert_ascent(data)
    # the indoor rows (21 C) before launch and the fall after the burst are gone
    assert len(data["p"]) == 80
    assert data["T"][0].m < 10
    assert data["p"][-1].m == pytest.approx(185.0)
    site = data["site_info"]
    assert site["site-name"] == "Grand Forks, ND"
    assert site["site-latlon"] == [47.92, -97.09]
    assert site["site-elv"] == 256
    assert site["valid-time"] == [2026, 1, 16, 18]
    assert site["launch-time"] == "2026-01-16T18:37:33Z"


def test_from_imet_missing_surface_winds_warn_and_flag(tmp_path):
    data, caught = read_format("imet", tmp_path, imet_rows(no_wind=3))
    # the 99999 placeholders become NaN, never a 99999 kt wind
    assert np.isnan(data["u"].m[:3]).all() and np.isfinite(data["u"].m[3])
    assert np.nanmax(data["ws"].m) < 200
    assert data["site_info"]["missing-sfc-wind"] is True
    assert any("no surface wind: the lowest 3 level(s)" in message for message in data_warnings(caught))


def test_from_windsond_heights_and_elevation(tmp_path):
    data, _ = read_format("windsond", tmp_path)
    assert_ascent(data)
    p, z, _ = standard_ascent(p0=980.0)
    np.testing.assert_allclose(data["z"].m, np.round(z - z[0] + 20) + 257)
    site = data["site_info"]
    assert site["site-elv"] == 257
    assert site["site-latlon"] == [47.93, -97.12]
    assert site["valid-time"] == [2026, 1, 18, 15]
    assert site["launch-time"] == "2026-01-18T15:56:00Z"


def test_from_windsond_requires_an_elevation(tmp_path):
    path = write_windsond(tmp_path / "20260118_1550.csv")
    with pytest.raises(TypeError):
        spy.from_windsond(path)
    for bad in (None, "high", float("nan")):
        with pytest.raises(ValueError, match="elevation"):
            spy.from_windsond(path, elevation=bad)
    assert read_quietly(spy.from_windsond, path, elevation="257")[0]["site_info"]["site-elv"] == 257


@pytest.mark.parametrize("fmt", FORMATS)
def test_data_gaps_warn_and_are_recorded(fmt, tmp_path):
    rows = FORMATS[fmt][2]()
    # a telemetry dropout: 20 levels (200 hPa) missing in the middle of the ascent (after iMet's 5 pre-launch rows)
    start = (5 if fmt == "imet" else 0) + 30
    rows = rows[:start] + rows[start + 20:]
    data, caught = read_format(fmt, tmp_path, rows)
    p = data["p"].m
    i = int(np.argmax(-np.diff(p)))
    assert data["site_info"]["data-gaps"] == [[round(float(p[i]), 1), round(float(p[i + 1]), 1)]]
    assert any(f"data gap from {p[i]:.1f} to {p[i + 1]:.1f} hPa" in message for message in data_warnings(caught))


@pytest.mark.parametrize("fmt", FORMATS)
def test_complete_files_have_no_gaps_or_warnings(fmt, tmp_path):
    data, caught = read_format(fmt, tmp_path)
    assert data["site_info"]["data-gaps"] == []
    assert data_warnings(caught) == []


def test_rows_missing_temperature_are_dropped(tmp_path):
    rows = graw_rows()
    for i in (0, 5, 6, 7):          # the surface observation and three levels aloft
        rows[i][3] = "      "
    data, _ = read_format("graw", tmp_path, rows)
    assert_ascent(data)
    assert len(data["p"]) == 76
    assert data["p"][0].m == pytest.approx(990.0)


def test_a_hole_left_by_missing_temperatures_is_a_gap(tmp_path):
    rows = graw_rows()
    for row in rows[30:50]:
        row[3] = "-----"
    data, caught = read_format("graw", tmp_path, rows)
    assert data["site_info"]["data-gaps"] == [[710.0, 500.0]]
    assert data_warnings(caught)


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("levels", [1, 5, 9])
def test_too_few_levels_raise(fmt, levels, tmp_path):
    rows = FORMATS[fmt][2]()
    rows = rows[:levels + 5] if fmt == "imet" else rows[:levels]
    with pytest.raises(ValueError, match="too short to plot"):
        read_format(fmt, tmp_path, rows)


def test_too_shallow_a_profile_raises(tmp_path):
    # 10 levels, but only 90 hPa deep (the mixed-layer parcel needs 100); one more level makes it
    with pytest.raises(ValueError, match="too short to plot: 10 level"):
        read_format("graw", tmp_path, graw_rows()[:10])
    data, _ = read_format("graw", tmp_path, graw_rows()[:11])
    assert len(data["p"]) == 11


def test_partial_in_flight_file_ends_at_the_last_level(tmp_path):
    # a file saved while the balloon is still rising: the top is just the last level received
    data, caught = read_format("graw", tmp_path, graw_rows()[:40])
    assert data["p"][-1].m == pytest.approx(610.0)
    assert data["site_info"]["data-gaps"] == []
    assert data_warnings(caught) == []


@pytest.mark.parametrize("fmt", FORMATS)
@pytest.mark.parametrize("bom, newline", [(True, "\r\n"), (True, "\n"), (False, "\n")])
def test_byte_order_mark_and_line_endings(fmt, bom, newline, tmp_path):
    plain, _ = read_format(fmt, tmp_path)
    (tmp_path / "other").mkdir()
    other, _ = read_format(fmt, tmp_path / "other", bom=bom, newline=newline)
    np.testing.assert_allclose(other["p"].m, plain["p"].m)
    np.testing.assert_allclose(other["T"].m, plain["T"].m)
    # Windsond's launch lat/lon and time are on its first line, where a byte order mark sits
    assert other["site_info"] == plain["site_info"]


@pytest.mark.parametrize("fmt", FORMATS)
def test_readers_warn_instead_of_printing(fmt, tmp_path, capsys):
    reader, write, _, kwargs, name = FORMATS[fmt]
    path = write(tmp_path / name)
    with pytest.warns(SounderPyNotice, match="Still needed: site ID"):
        reader(path, **kwargs)
    assert capsys.readouterr().out == ""


def test_from_sharppy_warns_instead_of_printing(clean_profile, tmp_path, capsys):
    path = tmp_path / "notice.snd"
    spy.to_file("sharppy", clean_profile, filename=str(path))
    capsys.readouterr()
    with pytest.warns(SounderPyNotice, match="site meta-data"):
        spy.from_sharppy(path)
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("reader, kwargs, message", [(spy.from_graw, {}, "empty"), (spy.from_imet, {}, "empty"),
                                                     (spy.from_windsond, {"elevation": 257}, "empty")])
def test_radiosonde_readers_reject_empty_files(reader, kwargs, message, tmp_path):
    path = tmp_path / "20260101_0000.txt"
    path.write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        reader(path, **kwargs)


@pytest.mark.parametrize("fmt, name", [("graw", "GRAW"), ("imet", "iMet"), ("windsond", "Windsond")])
def test_radiosonde_readers_reject_other_files(fmt, name, tmp_path):
    reader, _, _, kwargs, _ = FORMATS[fmt]
    garbage = tmp_path / "20260101_0000.bin"
    garbage.write_bytes(bytes(range(256)) * 20)
    text = tmp_path / "20260101_0000.txt"
    text.write_text("this is not a sounding\nnor is this\n", encoding="utf-8")
    # and each of the other formats
    others = [FORMATS[other][1](tmp_path / f"20260101_0000_{other}.txt") for other in FORMATS if other != fmt]
    for path in [garbage, text] + others:
        with pytest.raises(ValueError, match=f"does not appear to be a.* {name}"):
            reader(path, **kwargs)


@pytest.mark.slow
@pytest.mark.parametrize("fmt", FORMATS)
def test_radiosonde_files_plot(fmt, tmp_path):
    from helpers import assert_png
    data, _ = read_format(fmt, tmp_path)
    data["site_info"].update({"site-id": "UND", "site-name": "UND", "site-lctn": "ND"})
    data["titles"] = {"top_title": "RADIOSONDE TEST", "left_title": "TEST", "right_title": "UND"}
    out = tmp_path / "radiosonde.png"
    spy.build_sounding(data, mode="winter", save=True, filename=str(out), map_zoom=0, dpi=50)
    assert_png(out)


@pytest.mark.parametrize("fmt, raw_speed_ms, raw_direction", [("graw", 4.0, 300.0), ("imet", 12.0, 330.0),
                                                              ("windsond", 11.81, 353.1)])
def test_radiosonde_winds_are_knots_from_the_direction(fmt, raw_speed_ms, raw_direction, tmp_path):
    data, _ = read_format(fmt, tmp_path)
    for key in ("ws", "u", "v"):
        assert str(data[key].units) == "knot", key
    assert str(data["wd"].units) == "degree"
    speed_kt = raw_speed_ms * 1.943844
    assert data["ws"][0].m == pytest.approx(speed_kt, rel=1e-4)
    # meteorological convention: the direction the wind blows from
    assert data["u"][0].m == pytest.approx(-speed_kt * np.sin(np.radians(raw_direction)), abs=1e-3)
    assert data["v"][0].m == pytest.approx(-speed_kt * np.cos(np.radians(raw_direction)), abs=1e-3)


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

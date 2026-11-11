"""Tests for the SounderPy command-line interface.

These tests are intentionally offline. Retrieval, plotting, and export functions
are mocked so this module tests CLI parsing/routing/output behavior rather than
external data services or scientific calculations.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

import sounderpy.cli as cli


# ---------------------------------------------------------------------------
# Installed/module-level smoke tests
# ---------------------------------------------------------------------------

def test_import_sounderpy_is_silent():
    """Importing the library should not print banners or other text."""
    result = subprocess.run(
        [sys.executable, "-c", "import sounderpy"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert result.stdout == ""


def test_python_module_help_works():
    """``python -m sounderpy`` should expose the CLI."""
    result = subprocess.run(
        [sys.executable, "-m", "sounderpy", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert "usage: sounderpy" in result.stdout
    assert "obs" in result.stdout
    assert "model" in result.stdout
    assert "bufkit" in result.stdout
    assert "acars" in result.stdout


def test_python_module_version_works():
    result = subprocess.run(
        [sys.executable, "-m", "sounderpy", "--version"],
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert result.stdout.startswith("sounderpy ")


# ---------------------------------------------------------------------------
# Argument validation
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("argv", "message"),
    [
        (["obs", "DTX", "2024-13-21", "00"], "not a valid date"),
        (["obs", "DTX", "2024-05-21", "24"], "hour must be from 0 to 23"),
        (
            ["model", "rap-ruc", "91", "-100", "2024-08-28", "18"],
            "latitude must be between -90 and 90",
        ),
        (
            ["model", "rap-ruc", "44", "-181", "2024-08-28", "18"],
            "longitude must be between -180 and 180",
        ),
        (
            ["bufkit", "gfs", "KMOP", "-1"],
            "forecast hour cannot be negative",
        ),
    ],
)
def test_parser_rejects_invalid_arguments(argv, message, capsys):
    """argparse-level validation should fail before retrieval starts."""
    with pytest.raises(SystemExit) as exc:
        cli.main(argv)

    assert exc.value.code == 2
    captured = capsys.readouterr()
    assert message in captured.err


def test_model_rejects_nonpositive_box_size(monkeypatch, capsys):
    """Semantic validation performed after parsing should return CLI status 2."""

    def should_not_run(*args, **kwargs):
        raise AssertionError("retrieval should not have been called")

    monkeypatch.setattr(cli, "get_model_data", should_not_run)

    result = cli.main(
        [
            "model",
            "rap-ruc",
            "44.58",
            "-100.82",
            "2024-08-28",
            "18",
            "--box-size",
            "0",
        ]
    )

    assert result == 2
    assert "--box-size must be greater than zero" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# Retrieval command routing
# ---------------------------------------------------------------------------

def test_obs_command_routes_arguments(monkeypatch, clean_profile, capsys):
    called = {}

    def fake_get_obs_data(
        station,
        year,
        month,
        day,
        hour,
        hush=False,
        clean_it=True,
    ):
        called.update(
            station=station,
            year=year,
            month=month,
            day=day,
            hour=hour,
            hush=hush,
            clean_it=clean_it,
        )
        print("THIS SHOULD BE HIDDEN")
        return clean_profile

    monkeypatch.setattr(cli, "get_obs_data", fake_get_obs_data)

    result = cli.main(["obs", "DTX", "2024-05-21", "0"])

    assert result == 0
    assert called == {
        "station": "DTX",
        "year": "2024",
        "month": "05",
        "day": "21",
        "hour": "00",
        "hush": True,
        "clean_it": True,
    }

    captured = capsys.readouterr()
    assert "SounderPy profile" in captured.out
    assert "THIS SHOULD BE HIDDEN" not in captured.out


def test_obs_verbose_allows_library_output(monkeypatch, clean_profile, capsys):
    def fake_get_obs_data(*args, **kwargs):
        print("VISIBLE RETRIEVAL STATUS")
        return clean_profile

    monkeypatch.setattr(cli, "get_obs_data", fake_get_obs_data)

    result = cli.main(["obs", "DTX", "2024-05-21", "00", "--verbose"])

    assert result == 0
    assert "VISIBLE RETRIEVAL STATUS" in capsys.readouterr().out


def test_model_command_routes_arguments(monkeypatch, clean_profile):
    called = {}

    def fake_get_model_data(
        model,
        latlon,
        year,
        month,
        day,
        hour,
        dataset=None,
        box_avg_size=0.10,
        hush=False,
        clean_it=True,
    ):
        called.update(
            model=model,
            latlon=latlon,
            year=year,
            month=month,
            day=day,
            hour=hour,
            dataset=dataset,
            box_avg_size=box_avg_size,
            hush=hush,
            clean_it=clean_it,
        )
        return clean_profile

    monkeypatch.setattr(cli, "get_model_data", fake_get_model_data)

    result = cli.main(
        [
            "model",
            "rap-ruc",
            "44.58",
            "-100.82",
            "2024-08-28",
            "18",
            "--dataset",
            "aws",
            "--box-size",
            "0.25",
        ]
    )

    assert result == 0
    assert called == {
        "model": "rap-ruc",
        "latlon": [44.58, -100.82],
        "year": "2024",
        "month": "08",
        "day": "28",
        "hour": "18",
        "dataset": "aws",
        "box_avg_size": 0.25,
        "hush": True,
        "clean_it": True,
    }


def test_bufkit_latest_run_routes_none_date(monkeypatch, clean_profile):
    called = {}

    def fake_get_bufkit_data(
        model,
        station,
        forecast_hour,
        run_year=None,
        run_month=None,
        run_day=None,
        run_hour=None,
        hush=False,
        clean_it=True,
    ):
        called.update(
            model=model,
            station=station,
            forecast_hour=forecast_hour,
            run_year=run_year,
            run_month=run_month,
            run_day=run_day,
            run_hour=run_hour,
            hush=hush,
            clean_it=clean_it,
        )
        return clean_profile

    monkeypatch.setattr(cli, "get_bufkit_data", fake_get_bufkit_data)

    result = cli.main(["bufkit", "gfs", "KMOP", "6"])

    assert result == 0
    assert called["model"] == "gfs"
    assert called["station"] == "KMOP"
    assert called["forecast_hour"] == 6
    assert called["run_year"] is None
    assert called["run_month"] is None
    assert called["run_day"] is None
    assert called["run_hour"] is None
    assert called["hush"] is True
    assert called["clean_it"] is True


def test_bufkit_archived_run_routes_date(monkeypatch, clean_profile):
    called = {}

    def fake_get_bufkit_data(*args, **kwargs):
        called["args"] = args
        called["kwargs"] = kwargs
        return clean_profile

    monkeypatch.setattr(cli, "get_bufkit_data", fake_get_bufkit_data)

    result = cli.main(
        [
            "bufkit",
            "gfs",
            "KMOP",
            "6",
            "--run",
            "2023-08-05",
            "12",
        ]
    )

    assert result == 0
    assert called["args"] == ("gfs", "KMOP", 6)
    assert called["kwargs"]["run_year"] == "2023"
    assert called["kwargs"]["run_month"] == "08"
    assert called["kwargs"]["run_day"] == "05"
    assert called["kwargs"]["run_hour"] == "12"


# ---------------------------------------------------------------------------
# JSON behavior
# ---------------------------------------------------------------------------

def test_json_output_is_valid_and_quiet(monkeypatch, clean_profile, capsys):
    """stdout must contain only JSON, even if the underlying retriever prints."""

    def fake_get_obs_data(*args, **kwargs):
        print("NOISY INTERNAL OUTPUT")
        return clean_profile

    monkeypatch.setattr(cli, "get_obs_data", fake_get_obs_data)

    result = cli.main(
        ["obs", "DTX", "2024-05-21", "00", "--json"]
    )

    assert result == 0

    captured = capsys.readouterr()
    parsed = json.loads(captured.out)

    assert "NOISY INTERNAL OUTPUT" not in captured.out
    assert parsed["site_info"]["site-id"] == clean_profile["site_info"]["site-id"]

    # Pint/MetPy quantities should retain both numerical values and units.
    assert set(parsed["p"]) == {"value", "units"}
    assert isinstance(parsed["p"]["value"], list)
    assert parsed["p"]["units"]


def test_json_rejects_interactive_plot(monkeypatch, clean_profile, capsys):
    monkeypatch.setattr(
        cli,
        "get_obs_data",
        lambda *args, **kwargs: clean_profile,
    )

    result = cli.main(
        [
            "obs",
            "DTX",
            "2024-05-21",
            "00",
            "--json",
            "--plot",
            "sounding",
        ]
    )

    assert result == 2
    assert "--json cannot be combined with an interactive plot" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# Export behavior
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("filename", "expected_format"),
    [
        ("profile.csv", "csv"),
        ("profile.snd", "sharppy"),
        ("profile.sharp", "sharppy"),
        ("profile.sharppy", "sharppy"),
        ("input_sounding.cm1", "cm1"),
        ("input_sounding.input", "cm1"),
    ],
)
def test_export_format_is_inferred(
    monkeypatch,
    clean_profile,
    tmp_path,
    filename,
    expected_format,
):
    export_calls = []

    monkeypatch.setattr(
        cli,
        "get_obs_data",
        lambda *args, **kwargs: clean_profile,
    )

    def fake_to_file(file_type, data, filename=None, **kwargs):
        export_calls.append(
            {
                "file_type": file_type,
                "data": data,
                "filename": filename,
            }
        )

    monkeypatch.setattr(cli, "to_file", fake_to_file)

    output = tmp_path / filename

    result = cli.main(
        [
            "obs",
            "DTX",
            "2024-05-21",
            "00",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    assert len(export_calls) == 1
    assert export_calls[0]["file_type"] == expected_format
    assert export_calls[0]["data"] is clean_profile
    assert export_calls[0]["filename"] == str(output)


def test_unknown_export_extension_returns_error(
    monkeypatch,
    clean_profile,
    tmp_path,
    capsys,
):
    monkeypatch.setattr(
        cli,
        "get_obs_data",
        lambda *args, **kwargs: clean_profile,
    )

    output = tmp_path / "profile.unknown"

    result = cli.main(
        [
            "obs",
            "DTX",
            "2024-05-21",
            "00",
            "--output",
            str(output),
        ]
    )

    assert result == 2
    assert "Could not infer export format" in capsys.readouterr().err


def test_explicit_export_format_overrides_extension(
    monkeypatch,
    clean_profile,
    tmp_path,
):
    calls = []

    monkeypatch.setattr(
        cli,
        "get_obs_data",
        lambda *args, **kwargs: clean_profile,
    )

    def fake_to_file(file_type, data, filename=None, **kwargs):
        calls.append((file_type, data, filename))

    monkeypatch.setattr(cli, "to_file", fake_to_file)

    output = tmp_path / "profile.dat"

    result = cli.main(
        [
            "obs",
            "DTX",
            "2024-05-21",
            "00",
            "--output",
            str(output),
            "--format",
            "csv",
        ]
    )

    assert result == 0
    assert calls == [("csv", clean_profile, str(output))]


# ---------------------------------------------------------------------------
# Plot routing
# ---------------------------------------------------------------------------

def test_sounding_plot_routes_options(
    monkeypatch,
    clean_profile,
    tmp_path,
):
    calls = []

    monkeypatch.setattr(
        cli,
        "get_obs_data",
        lambda *args, **kwargs: clean_profile,
    )

    def fake_build_sounding(data, **kwargs):
        calls.append((data, kwargs))

    monkeypatch.setattr(cli, "build_sounding", fake_build_sounding)

    output = tmp_path / "plots" / "sounding.png"

    result = cli.main(
        [
            "obs",
            "DTX",
            "2024-05-21",
            "00",
            "--plot",
            "sounding",
            "--plot-file",
            str(output),
            "--map-zoom",
            "0",
            "--dpi",
            "150",
            "--dark-mode",
            "--color-blind",
            "--show-theta",
        ]
    )

    assert result == 0
    assert len(calls) == 1

    data, kwargs = calls[0]
    assert data is clean_profile
    assert kwargs["save"] is True
    assert kwargs["filename"] == str(output)
    assert kwargs["map_zoom"] == 0
    assert kwargs["dpi"] == 150
    assert kwargs["dark_mode"] is True
    assert kwargs["color_blind"] is True
    assert kwargs["show_theta"] is True
    assert kwargs["mode"] == "severe"
    assert kwargs["show_obs"] is False


def test_sounding_plot_routes_winter_options(
    monkeypatch,
    clean_profile,
    tmp_path,
):
    calls = []

    monkeypatch.setattr(
        cli,
        "get_bufkit_data",
        lambda *args, **kwargs: clean_profile,
    )
    monkeypatch.setattr(
        cli,
        "build_sounding",
        lambda data, **kwargs: calls.append(kwargs),
    )

    result = cli.main(
        [
            "bufkit",
            "hrrr",
            "KGFK",
            "1",
            "--plot",
            "sounding",
            "--plot-file",
            str(tmp_path / "winter.png"),
            "--mode",
            "winter",
            "--obs",
        ]
    )

    assert result == 0
    assert calls[0]["mode"] == "winter"
    assert calls[0]["show_obs"] is True


def test_parser_rejects_unknown_mode(capsys):
    with pytest.raises(SystemExit):
        cli.main(["obs", "DTX", "2024-05-21", "00", "--plot", "sounding", "--mode", "snow"])
    assert "invalid choice: 'snow'" in capsys.readouterr().err


def test_hodograph_plot_routes_storm_relative(
    monkeypatch,
    clean_profile,
    tmp_path,
):
    calls = []

    monkeypatch.setattr(
        cli,
        "get_obs_data",
        lambda *args, **kwargs: clean_profile,
    )

    def fake_build_hodograph(data, **kwargs):
        calls.append((data, kwargs))

    monkeypatch.setattr(cli, "build_hodograph", fake_build_hodograph)

    output = tmp_path / "hodograph.png"

    result = cli.main(
        [
            "obs",
            "DTX",
            "2024-05-21",
            "00",
            "--plot",
            "hodograph",
            "--plot-file",
            str(output),
            "--storm-relative",
            "--map-zoom",
            "0",
        ]
    )

    assert result == 0
    assert len(calls) == 1

    data, kwargs = calls[0]
    assert data is clean_profile
    assert kwargs["save"] is True
    assert kwargs["filename"] == str(output)
    assert kwargs["sr_hodo"] is True
    assert kwargs["map_zoom"] == 0


# ---------------------------------------------------------------------------
# ACARS command routing
# ---------------------------------------------------------------------------

def test_acars_list(monkeypatch, capsys):
    constructor_args = {}

    class FakeACARS:
        def __init__(self, year, month, day, hour):
            constructor_args["values"] = (year, month, day, hour)

        def list_profiles(self):
            print("HIDDEN ACARS STATUS")
            return ["BNA_2320", "DFW_0015"]

    monkeypatch.setattr(cli, "acars_data", FakeACARS)

    result = cli.main(["acars", "list", "2024-05-21", "18"])

    assert result == 0
    assert constructor_args["values"] == ("2024", "05", "21", "18")

    captured = capsys.readouterr()
    assert "HIDDEN ACARS STATUS" not in captured.out
    assert captured.out.splitlines() == ["BNA_2320", "DFW_0015"]


def test_acars_list_json(monkeypatch, capsys):
    class FakeACARS:
        def __init__(self, *args):
            pass

        def list_profiles(self):
            return ["BNA_2320", "DFW_0015"]

    monkeypatch.setattr(cli, "acars_data", FakeACARS)

    result = cli.main(
        ["acars", "list", "2024-05-21", "18", "--json"]
    )

    assert result == 0
    assert json.loads(capsys.readouterr().out) == ["BNA_2320", "DFW_0015"]


def test_acars_get(monkeypatch, clean_profile, capsys):
    calls = {}

    class FakeACARS:
        def __init__(self, year, month, day, hour):
            calls["constructor"] = (year, month, day, hour)

        def get_profile(self, profile, hush=False, clean_it=True):
            calls["get_profile"] = (profile, hush, clean_it)
            return clean_profile

    monkeypatch.setattr(cli, "acars_data", FakeACARS)

    result = cli.main(
        ["acars", "get", "2024-05-21", "18", "BNA_2320"]
    )

    assert result == 0
    assert calls["constructor"] == ("2024", "05", "21", "18")
    assert calls["get_profile"] == ("BNA_2320", True, True)
    assert "SounderPy profile" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

def test_retrieval_exception_becomes_nonzero_cli_exit(monkeypatch, capsys):
    def fake_get_obs_data(*args, **kwargs):
        raise ValueError("synthetic retrieval failure")

    monkeypatch.setattr(cli, "get_obs_data", fake_get_obs_data)

    result = cli.main(["obs", "DTX", "2024-05-21", "00"])

    assert result == 1
    assert "synthetic retrieval failure" in capsys.readouterr().err


def test_legacy_system_exit_is_converted_to_cli_error(monkeypatch, capsys):
    """Temporary regression test while legacy internals still call sys.exit()."""

    def fake_get_obs_data(*args, **kwargs):
        raise SystemExit("legacy retrieval failure")

    monkeypatch.setattr(cli, "get_obs_data", fake_get_obs_data)

    result = cli.main(["obs", "DTX", "2024-05-21", "00"])

    assert result == 1
    assert "legacy retrieval failure" in capsys.readouterr().err


def test_verbose_shows_the_traceback(monkeypatch, capsys):
    def fake_get_obs_data(*args, **kwargs):
        raise IndexError()

    monkeypatch.setattr(cli, "get_obs_data", fake_get_obs_data)

    assert cli.main(["obs", "DTX", "2024-05-21", "00"]) == 1
    err = capsys.readouterr().err
    assert "Traceback" not in err
    assert err.strip() == "sounderpy: error: IndexError"  # an exception without a message is named

    assert cli.main(["obs", "DTX", "2024-05-21", "00", "-v"]) == 1
    assert "Traceback (most recent call last)" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# Options checked before any download
# ---------------------------------------------------------------------------

OBS = ["obs", "DTX", "2024-05-21", "00"]


@pytest.mark.parametrize(
    ("options", "message"),
    [
        (["--json", "--plot", "sounding"], "--json cannot be combined with an interactive plot"),
        (["--output", "profile.unknown"], "Could not infer export format"),
        (["--format", "csv"], "--format needs --output"),
        (["--plot-file", "x.png"], "--plot-file needs --plot sounding or --plot hodograph"),
        (["--dark-mode"], "--dark-mode needs --plot"),
        (["--map-zoom", "0"], "--map-zoom needs --plot"),
        (["--mode", "winter"], "--mode needs --plot"),
        (["--plot", "hodograph", "--plot-file", "h.png", "--color-blind"], "--color-blind only applies to --plot sounding"),
        (["--plot", "hodograph", "--plot-file", "h.png", "--mode", "winter"], "--mode only applies to --plot sounding"),
        (["--plot", "sounding", "--plot-file", "s.png", "--obs"], "--obs only applies to winter soundings"),
        (["--plot", "sounding", "--plot-file", "s.png", "--mode", "severe", "--obs"], "--obs only applies to winter soundings"),
        (["--plot", "sounding", "--plot-file", "s.xyz"], "--plot-file can't be saved as '.xyz'"),
        # the test suite runs Matplotlib's non-interactive Agg backend, as on a server
        (["--plot", "sounding"], "No display is available to open the plot"),
    ],
)
def test_bad_options_fail_before_retrieval(monkeypatch, capsys, tmp_path, options, message):
    def should_not_run(*args, **kwargs):
        raise AssertionError("retrieval should not have been called")

    monkeypatch.setattr(cli, "get_obs_data", should_not_run)
    monkeypatch.chdir(tmp_path)

    assert cli.main(OBS + options) == 2
    assert message in capsys.readouterr().err
    assert not list(tmp_path.iterdir())  # nothing written


@pytest.mark.parametrize(("option", "value"), [("--dpi", "0"), ("--dpi", "high"), ("--map-zoom", "-1")])
def test_parser_rejects_bad_plot_numbers(option, value, capsys):
    with pytest.raises(SystemExit):
        cli.main(OBS + ["--plot", "sounding", option, value])
    assert f"argument {option}" in capsys.readouterr().err


def test_plot_file_without_extension_is_saved_as_png(monkeypatch, clean_profile, tmp_path, capsys):
    calls = []
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)
    monkeypatch.setattr(cli, "build_hodograph", lambda data, **kwargs: calls.append(kwargs))

    assert cli.main(OBS + ["--plot", "hodograph", "--plot-file", str(tmp_path / "hodo"), "-v"]) == 0
    assert calls[0]["filename"] == str(tmp_path / "hodo.png")
    assert f"Saved plot: {tmp_path / 'hodo.png'}" in capsys.readouterr().err


@pytest.mark.parametrize("extension", [".pdf", ".svg", ".JPG"])
def test_plot_file_keeps_other_formats(monkeypatch, clean_profile, tmp_path, extension):
    calls = []
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)
    monkeypatch.setattr(cli, "build_sounding", lambda data, **kwargs: calls.append(kwargs))

    assert cli.main(OBS + ["--plot", "sounding", "--plot-file", str(tmp_path / f"s{extension}")]) == 0
    assert calls[0]["filename"] == str(tmp_path / f"s{extension}")


def test_plot_defaults_are_filled_in(monkeypatch, clean_profile, tmp_path):
    calls = []
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)
    monkeypatch.setattr(cli, "build_sounding", lambda data, **kwargs: calls.append(kwargs))

    assert cli.main(OBS + ["--plot", "sounding", "--plot-file", str(tmp_path / "s.png"), "--mode", "auto", "--obs"]) == 0
    assert (calls[0]["map_zoom"], calls[0]["dpi"], calls[0]["mode"], calls[0]["show_obs"]) == (2, 100, "auto", True)


# ---------------------------------------------------------------------------
# Quiet output and the summary
# ---------------------------------------------------------------------------

def test_library_warnings_are_hidden_unless_verbose(monkeypatch, clean_profile, capsys):
    import logging
    import warnings

    def noisy_get_obs_data(*args, **kwargs):
        warnings.warn("SYNTHETIC LIBRARY WARNING")
        logging.getLogger("matplotlib.axes").warning("SYNTHETIC MATPLOTLIB LOG")
        return clean_profile

    monkeypatch.setattr(cli, "get_obs_data", noisy_get_obs_data)

    assert cli.main(OBS) == 0
    err = capsys.readouterr().err
    assert "SYNTHETIC LIBRARY WARNING" not in err and "SYNTHETIC MATPLOTLIB LOG" not in err

    with pytest.warns(UserWarning, match="SYNTHETIC LIBRARY WARNING"):
        assert cli.main(OBS + ["-v"]) == 0


def test_summary_formats_lat_lon_and_hides_placeholders(monkeypatch, clean_profile, capsys):
    import numpy as np

    clean_profile["site_info"].update({"site-latlon": [np.float64(43.37), np.float64(-84.44)],
                                       "site-lctn": "no-site-location", "model": "no-model"})
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)

    assert cli.main(OBS) == 0
    out = capsys.readouterr().out
    assert "Lat/Lon: 43.37, -84.44" in out
    assert "np.float64" not in out and "no-site-location" not in out and "Model:" not in out


# ---------------------------------------------------------------------------
# Calculation and plot options: --storm-motion, --modify-sfc, --radar, --parcels
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    ("options", "message"),
    [
        (["--storm-motion", "left_moving"], "--storm-motion needs --plot or --params"),
        (["--modify-sfc", "T=30"], "--modify-sfc needs --plot or --params"),
        (["--mode", "winter"], "--mode needs --plot or --params"),
        (["--radar", "single"], "--radar needs --plot"),
        (["--parcels", "simple"], "--parcels needs --plot"),
        (["--params", "--storm-motion", "sideways"], "--storm-motion must be one of"),
        (["--params", "--storm-motion", "270"], "--storm-motion must be one of"),
        (["--params", "--storm-motion", "400", "20"], "direction from 0 to 360"),
        (["--params", "--modify-sfc", "T30"], "--modify-sfc takes KEY=VALUE pairs"),
        (["--params", "--modify-sfc", "P=900"], "--modify-sfc takes KEY=VALUE pairs"),
        (["--params", "--modify-sfc", "T=warm"], "--modify-sfc T must be a number"),
        (["--params", "--modify-sfc", "T=30", "T=31"], "has T more than once"),
        (["--params", "--modify-sfc", "method=blend"], "needs at least one of T, Td"),
        (["--params", "--modify-sfc", "ws=20"], "ws and wd go together"),
        (["--params", "--modify-sfc", "T=30", "method=mixed"], "method must be mix, replace, or blend"),
        (["--params", "--modify-sfc", "T=30", "depth=300"], "depth only applies to method=blend"),
        (["--params", "--modify-sfc", "T=30", "method=blend", "depth=0"], "depth must be a height in m above 0"),
        (["--plot", "hodograph", "--plot-file", "h.png", "--parcels", "simple"], "--parcels only applies to --plot sounding"),
        (["--plot", "hodograph", "--plot-file", "h.png", "--mode", "auto"], "--mode only applies to --plot sounding and --params"),
        (["--plot", "sounding", "--plot-file", "s.png", "--mode", "winter", "--parcels", "simple"], "--parcels doesn't apply to winter soundings"),
        (["--plot", "sounding", "--plot-file", "s.png", "--mode", "winter", "--storm-relative"], "--storm-relative doesn't apply to winter soundings"),
        (["--plot", "sounding", "--plot-file", "s.png", "--parcels", "xx_ps_cape"], "--parcels doesn't know xx_ps_cape"),
        (["--plot", "sounding", "--plot-file", "s.png", "--parcels", "simple", "--background-parcels", "mu_ia_cape"], "can't be combined"),
        (["--plot", "sounding", "--plot-file", "s.png", "--map-zoom", "0", "--radar", "single"], "--radar needs the map"),
    ],
)
def test_bad_calculation_and_plot_options_fail_before_retrieval(monkeypatch, capsys, tmp_path, options, message):
    def should_not_run(*args, **kwargs):
        raise AssertionError("retrieval should not have been called")

    monkeypatch.setattr(cli, "get_obs_data", should_not_run)
    monkeypatch.chdir(tmp_path)

    assert cli.main(OBS + options) == 2
    assert message in capsys.readouterr().err


@pytest.mark.parametrize(
    ("options", "expected"),
    [
        ([], {"storm_motion": "right_moving", "modify_sfc": None, "radar": "mosaic", "radar_time": "sounding",
              "special_parcels": None}),
        (["--storm-motion", "Left_Moving"], {"storm_motion": "left_moving"}),
        (["--storm-motion", "250", "30.5"], {"storm_motion": [250.0, 30.5]}),
        (["--modify-sfc", "T=-6", "Td=-7", "ws=25", "wd=20"], {"modify_sfc": {"T": -6.0, "Td": -7.0, "ws": 25.0, "wd": 20.0}}),
        (["--modify-sfc", "T=30", "method=blend", "depth=300"], {"modify_sfc": {"T": 30.0, "method": "blend", "depth": 300.0}}),
        (["--radar", "single", "--radar-time", "now"], {"radar": "single", "radar_time": "now"}),
        (["--parcels", "simple"], {"special_parcels": "simple"}),
        (["--parcels", "mu_ia_ecape", "sb_ps_cape"], {"special_parcels": [["mu_ia_ecape", "sb_ps_cape"], []]}),
        (["--background-parcels", "ml_ps_cape"], {"special_parcels": [[], ["ml_ps_cape"]]}),
    ],
)
def test_sounding_routes_calculation_and_plot_options(monkeypatch, clean_profile, tmp_path, options, expected):
    calls = []
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)
    monkeypatch.setattr(cli, "build_sounding", lambda data, **kwargs: calls.append(kwargs))

    assert cli.main(OBS + ["--plot", "sounding", "--plot-file", str(tmp_path / "s.png")] + options) == 0
    assert {key: calls[0][key] for key in expected} == expected


def test_hodograph_routes_storm_motion_modify_sfc_and_radar(monkeypatch, clean_profile, tmp_path):
    calls = []
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)
    monkeypatch.setattr(cli, "build_hodograph", lambda data, **kwargs: calls.append(kwargs))

    assert cli.main(OBS + ["--plot", "hodograph", "--plot-file", str(tmp_path / "h.png"), "--storm-motion", "mean_wind",
                           "--modify-sfc", "ws=10", "wd=180", "--radar", "none"]) == 0
    assert calls[0]["storm_motion"] == "mean_wind"
    assert calls[0]["modify_sfc"] == {"ws": 10.0, "wd": 180.0}
    assert calls[0]["radar"] == "none"
    assert "special_parcels" not in calls[0] and "mode" not in calls[0]


# ---------------------------------------------------------------------------
# --params
# ---------------------------------------------------------------------------

def test_params_readout(monkeypatch, clean_profile, capsys):
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)

    assert cli.main(OBS + ["--params", "--storm-motion", "250", "30", "--modify-sfc", "T=35"]) == 0
    out = capsys.readouterr().out
    assert "SounderPy profile" in out and "Parameters" in out
    assert "Modified surface: T 35°C · mixed to 1.0 km" in out
    for line in ("SB parcel: CAPE ", "ML parcel: CAPE ", "MU parcel: CAPE ", "Storm motion (custom): ", "SRH: 0-500 m ",
                 "Bulk shear: 0-1 km ", "STP: "):
        assert line in out
    assert "Precip type" not in out


def test_params_include_winter_values(monkeypatch, capsys):
    from helpers import make_winter_profile

    monkeypatch.setattr(cli, "get_bufkit_data", lambda *args, **kwargs: make_winter_profile("fzra"))

    assert cli.main(["bufkit", "hrrr", "KGFK", "1", "--params", "--mode", "auto"]) == 0
    out = capsys.readouterr().out
    assert "Precip type: FREEZING RAIN (SPC: Freezing Rain)" in out
    assert "Phase from source to surface: SN -> RA -> FZRA" in out


def test_params_json(monkeypatch, clean_profile, capsys):
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)

    assert cli.main(OBS + ["--params", "--json"]) == 0
    parsed = json.loads(capsys.readouterr().out)
    assert {"general", "thermo", "kinem", "intrp"} <= set(parsed["params"])
    assert parsed["params"]["thermo"]["mlcape"] > 0
    assert parsed["site_info"]["site-id"] == "TEST"  # the profile is still at the top level


def test_json_without_params_has_no_params(monkeypatch, clean_profile, capsys):
    monkeypatch.setattr(cli, "get_obs_data", lambda *args, **kwargs: clean_profile)
    assert cli.main(OBS + ["--json"]) == 0
    assert "params" not in json.loads(capsys.readouterr().out)


# ---------------------------------------------------------------------------
# file
# ---------------------------------------------------------------------------

@pytest.fixture
def sharppy_file(clean_profile, tmp_path):
    from sounderpy.data.file_io import to_file

    path = tmp_path / "profile.snd"
    to_file("sharppy", clean_profile, filename=str(path))
    return path


def test_file_summary(sharppy_file, capsys):
    assert cli.main(["file", str(sharppy_file)]) == 0
    out = capsys.readouterr().out
    assert "Source: SHARPPY FILE" in out and "Site: TEST" in out
    assert "Valid time: 2024-05-21 00Z" in out


def test_file_export(sharppy_file, tmp_path):
    output = tmp_path / "profile.csv"
    assert cli.main(["file", str(sharppy_file), "-o", str(output)]) == 0
    assert output.read_text().splitlines()[0] == "p,z,T,Td,u,v"


def test_file_plot_and_params_need_a_location(sharppy_file, capsys):
    assert cli.main(["file", str(sharppy_file), "--params"]) == 2
    assert "--params needs the profile's location; add --latlon LAT LON" in capsys.readouterr().err
    assert cli.main(["file", str(sharppy_file), "--plot", "sounding", "--plot-file", "x.png"]) == 2
    assert "--plot needs the profile's location" in capsys.readouterr().err


def test_file_plot_fills_in_site_details(monkeypatch, sharppy_file, tmp_path):
    calls = []
    monkeypatch.setattr(cli, "build_sounding", lambda data, **kwargs: calls.append(data))

    assert cli.main(["file", str(sharppy_file), "--latlon", "47.9", "-97.0", "--plot", "sounding",
                     "--plot-file", str(tmp_path / "s.png")]) == 0
    site, titles = calls[0]["site_info"], calls[0]["titles"]
    assert site["site-latlon"] == [47.9, -97.0]
    assert site["site-elv"] == pytest.approx(0.0)  # the lowest height in the file
    assert titles == {"top_title": "SHARPPY FILE PROFILE | profile.snd", "left_title": "VALID: 05/21/2024 00Z",
                      "right_title": "TEST | 47.9, -97"}


def test_file_params(sharppy_file, capsys):
    assert cli.main(["file", str(sharppy_file), "--latlon", "47.9", "-97.0", "--params"]) == 0
    assert "SB parcel: CAPE " in capsys.readouterr().out


@pytest.mark.parametrize(("argv", "message"), [
    (["file", "missing.snd"], "No such file: missing.snd"),
    (["file", "x.snd", "--latlon", "95", "0"], "--latlon needs a latitude from -90 to 90"),
])
def test_file_errors(argv, message, capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "x.snd").write_text("")
    assert cli.main(argv) == 2
    assert message in capsys.readouterr().err


def test_file_that_is_not_sharppy(tmp_path, capsys):
    path = tmp_path / "notes.txt"
    path.write_text("not a sounding\n")
    assert cli.main(["file", str(path)]) == 1
    assert "does not appear to be a SHARPpy" in capsys.readouterr().err


def test_sounderpy_warnings_are_shown_without_verbose(monkeypatch, clean_profile, capsys):
    import warnings
    from sounderpy.utils import SounderPyWarning

    def get_obs_data(*args, **kwargs):
        warnings.warn("Radar mosaic could not be plotted: no radar mosaic for X", SounderPyWarning)
        warnings.warn("LIBRARY NOISE")
        return clean_profile

    monkeypatch.setattr(cli, "get_obs_data", get_obs_data)

    assert cli.main(OBS) == 0
    err = capsys.readouterr().err
    assert err.splitlines() == ["sounderpy: warning: Radar mosaic could not be plotted: no radar mosaic for X"]

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

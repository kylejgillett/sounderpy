"""Command-line interface for SounderPy.

The CLI is intentionally a thin wrapper around SounderPy's public Python API.
Scientific calculations, retrieval, plotting, and export logic remain in the
package modules themselves.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
from contextlib import redirect_stdout
from datetime import date as date_type
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .acars_data import acars_data
from .file_io import to_file
from .sounderpy import (
    build_hodograph,
    build_sounding,
    get_bufkit_data,
    get_model_data,
    get_obs_data,
)


class CLIError(Exception):
    """User-facing CLI error."""


def _package_version() -> str:
    """Return the installed SounderPy version."""
    try:
        return version("sounderpy")
    except PackageNotFoundError:
        return "unknown"


def _parse_date(value: str) -> tuple[str, str, str]:
    """Parse YYYY-MM-DD and return zero-padded year, month, and day strings."""
    try:
        parsed = datetime.strptime(value, "%Y-%m-%d")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"{value!r} is not a valid date; use YYYY-MM-DD"
        ) from exc

    return f"{parsed.year:04d}", f"{parsed.month:02d}", f"{parsed.day:02d}"


def _parse_hour(value: str) -> str:
    """Parse a UTC hour and return HH."""
    try:
        hour = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("hour must be an integer from 0 to 23") from exc

    if not 0 <= hour <= 23:
        raise argparse.ArgumentTypeError("hour must be from 0 to 23")

    return f"{hour:02d}"


def _latitude(value: str) -> float:
    try:
        lat = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("latitude must be numeric") from exc
    if not -90.0 <= lat <= 90.0:
        raise argparse.ArgumentTypeError("latitude must be between -90 and 90")
    return lat


def _longitude(value: str) -> float:
    try:
        lon = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("longitude must be numeric") from exc
    if not -180.0 <= lon <= 180.0:
        raise argparse.ArgumentTypeError("longitude must be between -180 and 180")
    return lon


def _forecast_hour(value: str) -> int:
    try:
        hour = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("forecast hour must be an integer") from exc
    if hour < 0:
        raise argparse.ArgumentTypeError("forecast hour cannot be negative")
    return hour


def _run_quietly(
    func: Callable[..., Any],
    *args: Any,
    verbose: bool = False,
    **kwargs: Any,
) -> Any:
    """Run a SounderPy call while hiding library progress output unless requested."""
    if verbose:
        return func(*args, **kwargs)
    with redirect_stdout(io.StringIO()):
        return func(*args, **kwargs)


def _json_ready(value: Any) -> Any:
    """Convert SounderPy/NumPy/Pint objects into JSON-safe structures."""
    if value is None or value is np.ma.masked:
        return None
    if isinstance(value, dict):
        return {str(key): _json_ready(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, set):
        return [_json_ready(item) for item in sorted(value, key=str)]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (datetime, date_type)):
        return value.isoformat()
    if hasattr(value, "magnitude") and hasattr(value, "units"):
        return {
            "value": _json_ready(value.magnitude),
            "units": str(value.units),
        }
    if isinstance(value, np.ma.MaskedArray):
        return _json_ready(value.filled(np.nan))
    if isinstance(value, np.ndarray):
        return [_json_ready(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_ready(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if hasattr(value, "tolist"):
        try:
            return _json_ready(value.tolist())
        except Exception:
            pass
    return str(value)


def _format_valid_time(value: Any) -> str:
    if isinstance(value, (list, tuple)) and len(value) >= 4:
        year, month, day, hour = value[:4]
        return f"{year}-{month}-{day} {hour}Z"
    return str(value)


def _quantity_endpoint(value: Any, index: int) -> str | None:
    """Return a compact formatted endpoint for a quantity/array."""
    try:
        item = value[index]
    except Exception:
        return None
    try:
        magnitude = float(item.magnitude)
        return f"{magnitude:.1f} {item.units:~P}"
    except Exception:
        try:
            return f"{float(item):.1f}"
        except Exception:
            return str(item)


def _print_summary(data: dict[str, Any]) -> None:
    """Print a concise terminal summary of a SounderPy clean-data dictionary."""
    site = data.get("site_info", {}) or {}

    print("SounderPy profile")
    print("-----------------")

    fields = (
        ("Source", site.get("source")),
        ("Site", site.get("site-id")),
        ("Name", site.get("site-name")),
        ("Location", site.get("site-lctn")),
        ("Lat/Lon", site.get("site-latlon")),
        ("Model", site.get("model")),
        ("Forecast hour", site.get("fcst-hour")),
        (
            "Valid time",
            _format_valid_time(site.get("valid-time"))
            if site.get("valid-time") is not None
            else None,
        ),
    )

    for label, value in fields:
        if value not in (None, "", "none", "no-model", "no-fcst-hour"):
            print(f"{label}: {value}")

    try:
        levels = len(data["p"])
    except Exception:
        levels = None

    if levels is not None:
        print(f"Vertical levels: {levels}")

    if "p" in data and levels:
        p0 = _quantity_endpoint(data["p"], 0)
        p1 = _quantity_endpoint(data["p"], -1)
        if p0 and p1:
            print(f"Pressure range: {p0} -> {p1}")

    if "z" in data and levels:
        z0 = _quantity_endpoint(data["z"], 0)
        z1 = _quantity_endpoint(data["z"], -1)
        if z0 and z1:
            print(f"Height range: {z0} -> {z1}")


def _infer_export_format(path: Path) -> str:
    suffix = path.suffix.casefold()
    if suffix == ".csv":
        return "csv"
    if suffix in {".snd", ".sharp", ".sharppy"}:
        return "sharppy"
    if suffix in {".cm1", ".input"}:
        return "cm1"
    raise CLIError(
        "Could not infer export format from the output filename. "
        "Use --format {csv,sharppy,cm1}."
    )


def _export_profile(
    data: dict[str, Any],
    output: str,
    file_format: str | None,
    verbose: bool,
) -> None:
    path = Path(output).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    fmt = file_format or _infer_export_format(path)
    _run_quietly(to_file, fmt, data, filename=str(path), verbose=verbose)
    if verbose:
        print(f"Saved data: {path}", file=sys.stderr)


def _plot_profile(data: dict[str, Any], args: argparse.Namespace) -> None:
    if args.json and args.plot and not args.plot_file:
        raise CLIError(
            "--json cannot be combined with an interactive plot. "
            "Use --plot-file to save the plot instead."
        )
    if not args.plot:
        return

    save = args.plot_file is not None
    filename = str(Path(args.plot_file).expanduser()) if save else None
    if save:
        Path(filename).parent.mkdir(parents=True, exist_ok=True)

    common = {
        "dark_mode": args.dark_mode,
        "map_zoom": args.map_zoom,
        "dpi": args.dpi,
        "save": save,
    }

    if args.plot == "sounding":
        kwargs = {
            **common,
            "color_blind": args.color_blind,
            "show_theta": args.show_theta,
        }
        if save:
            kwargs["filename"] = filename
        _run_quietly(build_sounding, data, verbose=args.verbose, **kwargs)

    elif args.plot == "hodograph":
        kwargs = {**common, "sr_hodo": args.storm_relative}
        if save:
            kwargs["filename"] = filename
        _run_quietly(build_hodograph, data, verbose=args.verbose, **kwargs)

    if save and args.verbose:
        print(f"Saved plot: {filename}", file=sys.stderr)


def _handle_profile_result(data: dict[str, Any], args: argparse.Namespace) -> int:
    if not isinstance(data, dict):
        raise CLIError(
            f"SounderPy returned {type(data).__name__}, but the CLI expected cleaned profile data."
        )
    if args.output:
        _export_profile(data, args.output, args.format, args.verbose)
    _plot_profile(data, args)
    if args.json:
        print(json.dumps(_json_ready(data), indent=2, allow_nan=False))
    else:
        _print_summary(data)
    return 0


def _add_profile_output_options(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("output")
    group.add_argument(
        "--json",
        action="store_true",
        help="Write the cleaned profile to stdout as JSON.",
    )
    group.add_argument(
        "-o", "--output", metavar="PATH", help="Export the retrieved profile to a file."
    )
    group.add_argument(
        "--format",
        choices=("csv", "sharppy", "cm1"),
        help="Export format. If omitted, infer it from --output.",
    )

    plot = parser.add_argument_group("plotting")
    plot.add_argument(
        "--plot",
        choices=("sounding", "hodograph"),
        help="Create a SounderPy plot after retrieval.",
    )
    plot.add_argument(
        "--plot-file",
        metavar="PATH",
        help="Save --plot to this file instead of opening it interactively.",
    )
    plot.add_argument(
        "--dark-mode", action="store_true", help="Use SounderPy's dark plot theme."
    )
    plot.add_argument(
        "--map-zoom",
        type=int,
        default=2,
        metavar="N",
        help="Map inset zoom level; use 0 to disable the map/radar inset. Default: 2.",
    )
    plot.add_argument(
        "--dpi",
        type=int,
        default=100,
        help="Plot resolution in dots per inch. Default: 100.",
    )
    plot.add_argument(
        "--color-blind",
        action="store_true",
        help="Use the color-deficiency-friendly sounding colors.",
    )
    plot.add_argument(
        "--show-theta",
        action="store_true",
        help="Show theta/theta-e on sounding plots.",
    )
    plot.add_argument(
        "--storm-relative",
        action="store_true",
        help="Use a storm-relative hodograph.",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show SounderPy's underlying progress/status output.",
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sounderpy",
        description="Retrieve, export, and plot atmospheric vertical profiles with SounderPy.",
        epilog=(
            "Examples:\n"
            "  sounderpy obs DTX 2024-05-21 00\n"
            "  sounderpy model rap-ruc 44.58 -100.82 2024-08-28 18 --box-size 0.25\n"
            "  sounderpy bufkit gfs KMOP 6 --run 2023-08-05 12\n"
            "  sounderpy acars list 2024-05-21 18\n"
            "  sounderpy obs DTX 2024-05-21 00 --plot sounding --plot-file dtx.png\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {_package_version()}"
    )

    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True)

    obs = subparsers.add_parser(
        "obs",
        help="Retrieve an observed RAOB/IGRA sounding.",
        description="Retrieve an observed RAOB/IGRA vertical profile.",
    )
    obs.add_argument("station", help="RAOB or IGRAv2 station identifier.")
    obs.add_argument("date", type=_parse_date, metavar="YYYY-MM-DD")
    obs.add_argument("hour", type=_parse_hour, metavar="HH")
    _add_profile_output_options(obs)

    model = subparsers.add_parser(
        "model",
        help="Retrieve model/reanalysis profile data.",
        description="Retrieve model or reanalysis vertical profile data.",
    )
    model.add_argument("model", help="Model name, e.g. rap-ruc, era5, or ncep.")
    model.add_argument("lat", type=_latitude, help="Latitude.")
    model.add_argument("lon", type=_longitude, help="Longitude.")
    model.add_argument("date", type=_parse_date, metavar="YYYY-MM-DD")
    model.add_argument("hour", type=_parse_hour, metavar="HH")
    model.add_argument(
        "--dataset",
        help="Target a specific underlying dataset instead of automatic discovery.",
    )
    model.add_argument(
        "--box-size",
        type=float,
        default=0.10,
        metavar="DEGREES",
        help="Area-average box size in degrees. Default: 0.10.",
    )
    _add_profile_output_options(model)

    bufkit = subparsers.add_parser(
        "bufkit",
        help="Retrieve a BUFKIT forecast sounding.",
        description="Retrieve a BUFKIT forecast vertical profile.",
    )
    bufkit.add_argument("model", help="BUFKIT model, e.g. gfs, hrrr, nam.")
    bufkit.add_argument("station", help="BUFKIT station identifier.")
    bufkit.add_argument(
        "forecast_hour", type=_forecast_hour, metavar="FHR", help="Forecast hour."
    )
    bufkit.add_argument(
        "--run",
        nargs=2,
        metavar=("YYYY-MM-DD", "HH"),
        help=(
            "Archived model initialization date/hour. "
            "If omitted, SounderPy requests the most recent run."
        ),
    )
    _add_profile_output_options(bufkit)

    acars = subparsers.add_parser(
        "acars", help="List or retrieve ACARS aircraft profiles."
    )
    acars_sub = acars.add_subparsers(
        dest="acars_command", metavar="ACTION", required=True
    )

    acars_list = acars_sub.add_parser(
        "list", help="List available ACARS profiles for a date/hour."
    )
    acars_list.add_argument("date", type=_parse_date, metavar="YYYY-MM-DD")
    acars_list.add_argument("hour", type=_parse_hour, metavar="HH")
    acars_list.add_argument(
        "--json", action="store_true", help="Write the available profile IDs as JSON."
    )
    acars_list.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show SounderPy's underlying progress/status output.",
    )

    acars_get = acars_sub.add_parser("get", help="Retrieve one ACARS profile.")
    acars_get.add_argument("date", type=_parse_date, metavar="YYYY-MM-DD")
    acars_get.add_argument("hour", type=_parse_hour, metavar="HH")
    acars_get.add_argument("profile", help="ACARS profile ID, e.g. BNA_2320.")
    _add_profile_output_options(acars_get)

    return parser


def _run_command(args: argparse.Namespace) -> int:
    if args.command == "obs":
        year, month, day = args.date
        data = _run_quietly(
            get_obs_data,
            args.station,
            year,
            month,
            day,
            args.hour,
            hush=not args.verbose,
            clean_it=True,
            verbose=args.verbose,
        )
        return _handle_profile_result(data, args)

    if args.command == "model":
        if args.box_size <= 0:
            raise CLIError("--box-size must be greater than zero")
        year, month, day = args.date
        data = _run_quietly(
            get_model_data,
            args.model,
            [args.lat, args.lon],
            year,
            month,
            day,
            args.hour,
            dataset=args.dataset,
            box_avg_size=args.box_size,
            hush=not args.verbose,
            clean_it=True,
            verbose=args.verbose,
        )
        return _handle_profile_result(data, args)

    if args.command == "bufkit":
        run_year = run_month = run_day = run_hour = None
        if args.run:
            try:
                run_year, run_month, run_day = _parse_date(args.run[0])
                run_hour = _parse_hour(args.run[1])
            except argparse.ArgumentTypeError as exc:
                raise CLIError(str(exc)) from exc

        data = _run_quietly(
            get_bufkit_data,
            args.model,
            args.station,
            args.forecast_hour,
            run_year=run_year,
            run_month=run_month,
            run_day=run_day,
            run_hour=run_hour,
            hush=not args.verbose,
            clean_it=True,
            verbose=args.verbose,
        )
        return _handle_profile_result(data, args)

    if args.command == "acars":
        year, month, day = args.date
        connection = acars_data(year, month, day, args.hour)

        if args.acars_command == "list":
            profiles = _run_quietly(connection.list_profiles, verbose=args.verbose)
            if args.json:
                print(json.dumps(_json_ready(profiles), indent=2))
            elif not profiles:
                print("No ACARS profiles found.")
            else:
                for profile in profiles:
                    print(profile)
            return 0

        if args.acars_command == "get":
            data = _run_quietly(
                connection.get_profile,
                args.profile,
                hush=not args.verbose,
                clean_it=True,
                verbose=args.verbose,
            )
            return _handle_profile_result(data, args)

    raise CLIError(f"Unknown command: {args.command}")


def main(argv: list[str] | None = None) -> int:
    """CLI entry point used by the ``sounderpy`` console script."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    try:
        return _run_command(args)
    except KeyboardInterrupt:
        print("\nsounderpy: interrupted", file=sys.stderr)
        return 130
    except CLIError as exc:
        print(f"sounderpy: error: {exc}", file=sys.stderr)
        return 2
    except SystemExit as exc:
        # Some legacy SounderPy internals currently call sys.exit().
        if exc.code not in (None, 0):
            print(f"sounderpy: error: {exc.code}", file=sys.stderr)
            return 1
        return 0
    except Exception as exc:
        print(f"sounderpy: error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

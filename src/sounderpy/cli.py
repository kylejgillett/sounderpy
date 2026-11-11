"""Command-line interface for SounderPy.

The CLI is intentionally a thin wrapper around SounderPy's public Python API.
Scientific calculations, retrieval, plotting, and export logic remain in the
package modules themselves.
"""

from __future__ import annotations

import argparse
import io
import json
import logging
import sys
import traceback
import warnings
from contextlib import contextmanager, redirect_stdout
from datetime import date as date_type
from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Callable

import numpy as np
from matplotlib.backend_bases import FigureCanvasBase

from .data.acars_data import acars_data
from .data.file_io import from_sharppy, to_file
from .utils import SounderPyWarning
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


def _positive_int(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a whole number") from exc
    if number <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return number


def _nonnegative_int(value: str) -> int:
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a whole number") from exc
    if number < 0:
        raise argparse.ArgumentTypeError("cannot be negative")
    return number


# plotting options, and the defaults of those left as None in the parser so the CLI can tell whether they were given
PLOT_OPTIONS = ("plot_file", "dark_mode", "map_zoom", "dpi", "color_blind", "show_theta", "storm_relative", "obs",
                "radar", "radar_time", "parcels", "background_parcels")
SOUNDING_ONLY_OPTIONS = ("color_blind", "show_theta", "obs", "parcels", "background_parcels")
# severe sounding options the winter sounding doesn't use
SEVERE_ONLY_OPTIONS = ("show_theta", "storm_relative", "parcels", "background_parcels")
# options that change the calculation, used by both --plot and --params
CALC_OPTIONS = ("mode", "storm_motion", "modify_sfc")
PLOT_DEFAULTS = {"map_zoom": 2, "dpi": 100, "mode": "severe", "radar": "mosaic", "radar_time": "sounding"}

STORM_MOTIONS = ("right_moving", "left_moving", "mean_wind")
MODIFY_SFC_KEYS = {"T": float, "Td": float, "ws": float, "wd": float, "method": str, "depth": float}


def _storm_motion(values: list[str] | None) -> str | list[float]:
    """--storm-motion as a name, or [direction, speed] from DIR SPD."""
    if values is None:
        return "right_moving"
    if len(values) == 1 and values[0].casefold() in STORM_MOTIONS:
        return values[0].casefold()
    if len(values) == 2:
        try:
            direction, speed = (float(value) for value in values)
        except ValueError:
            pass
        else:
            if not 0 <= direction <= 360 or speed < 0:
                raise CLIError("--storm-motion DIR SPD needs a direction from 0 to 360 degrees and a speed of 0 kt or more")
            return [direction, speed]
    raise CLIError(f"--storm-motion must be one of {', '.join(STORM_MOTIONS)}, or a direction and speed such as 270 25. "
                   f"Got: {' '.join(values)}")


def _modify_sfc(pairs: list[str] | None) -> dict[str, Any] | None:
    """--modify-sfc KEY=VALUE pairs as a modify_sfc dict, checked the way modify_surface needs them."""
    if pairs is None:
        return None
    sfc: dict[str, Any] = {}
    for pair in pairs:
        key, sep, value = pair.partition("=")
        if not sep or key not in MODIFY_SFC_KEYS:
            raise CLIError(f"--modify-sfc takes KEY=VALUE pairs with keys {', '.join(MODIFY_SFC_KEYS)}. Got: {pair!r}")
        if key in sfc:
            raise CLIError(f"--modify-sfc has {key} more than once")
        try:
            sfc[key] = MODIFY_SFC_KEYS[key](value)
        except ValueError as exc:
            raise CLIError(f"--modify-sfc {key} must be a number. Got: {value!r}") from exc

    if not {"T", "Td", "ws", "wd"} & set(sfc):
        raise CLIError("--modify-sfc needs at least one of T, Td, or ws and wd")
    if ("ws" in sfc) != ("wd" in sfc):
        raise CLIError("--modify-sfc ws and wd go together")
    if "wd" in sfc and not 0 <= sfc["wd"] <= 360 or sfc.get("ws", 0) < 0:
        raise CLIError("--modify-sfc wd must be from 0 to 360 degrees and ws 0 kt or more")
    if sfc.get("method", "mix") not in ("mix", "replace", "blend"):
        raise CLIError(f"--modify-sfc method must be mix, replace, or blend. Got: {sfc['method']!r}")
    if "depth" in sfc and sfc.get("method") != "blend":
        raise CLIError("--modify-sfc depth only applies to method=blend")
    if "depth" in sfc and not sfc["depth"] > 0:
        raise CLIError("--modify-sfc depth must be a height in m above 0")
    return sfc


def _special_parcels(parcels: list[str] | None, background: list[str] | None) -> str | list[list[str]] | None:
    """--parcels and --background-parcels as build_sounding's special_parcels."""
    from .plot.components import VALID_PARCEL_CODES

    if parcels == ["simple"]:
        if background:
            raise CLIError("--parcels simple can't be combined with --background-parcels")
        return "simple"
    if parcels is None and background is None:
        return None
    for flag, codes in (("--parcels", parcels or []), ("--background-parcels", background or [])):
        invalid = [code for code in codes if code not in VALID_PARCEL_CODES]
        if invalid:
            raise CLIError(f"{flag} doesn't know {', '.join(invalid)}. Valid parcels: {', '.join(VALID_PARCEL_CODES)}"
                           + (", or 'simple' on its own" if flag == "--parcels" else ""))
    return [list(parcels or []), list(background or [])]


def _flag(name: str) -> str:
    """The command-line flag for an argparse destination, e.g. plot_file -> --plot-file."""
    return "--" + name.replace("_", "-")


def _given(args: argparse.Namespace, name: str) -> bool:
    """Whether an option was given (flags are False and other options None when they aren't; 0 counts)."""
    value = getattr(args, name, None)
    return value is not None and value is not False


def _plot_path(value: str) -> str:
    """A --plot-file path with a format Matplotlib can write; .png is added when there is no extension."""
    path = Path(value).expanduser()
    if not path.suffix:
        path = path.with_name(path.name + ".png")
    supported = FigureCanvasBase.get_supported_filetypes()
    if path.suffix[1:].casefold() not in supported:
        raise CLIError(
            f"--plot-file can't be saved as '{path.suffix}'. "
            f"Use one of: {', '.join('.' + ext for ext in sorted(supported))}."
        )
    return str(path)


def _has_display() -> bool:
    """Whether Matplotlib can open a plot window (its backend isn't a file-only one such as Agg)."""
    import matplotlib.pyplot as plt

    backend = plt.get_backend().casefold()
    try:
        from matplotlib.backends import BackendFilter, backend_registry

        non_interactive = backend_registry.list_builtin(BackendFilter.NON_INTERACTIVE)
    except ImportError:  # Matplotlib < 3.9
        from matplotlib import rcsetup

        non_interactive = [name.casefold() for name in rcsetup.non_interactive_bk]
    return backend not in non_interactive


def _validate_options(args: argparse.Namespace) -> None:
    """Catch option combinations that can't work, before any data is downloaded."""
    if not hasattr(args, "plot"):  # `acars list` has no output or plot options
        return

    if args.format and not args.output:
        raise CLIError("--format needs --output")
    if args.output and not args.format:
        _infer_export_format(Path(args.output).expanduser())

    if not args.plot:
        given = [name for name in PLOT_OPTIONS if _given(args, name)]
        if given:
            raise CLIError(f"{_flag(given[0])} needs --plot sounding or --plot hodograph")
    if not args.plot and not args.params:
        given = [name for name in CALC_OPTIONS if _given(args, name)]
        if given:
            raise CLIError(f"{_flag(given[0])} needs --plot or --params")

    if args.command == "file" and (args.plot or args.params) and args.latlon is None:
        # the calculation uses the latitude (temperature advection), and the plot's map the location
        raise CLIError(f"{'--plot' if args.plot else '--params'} needs the profile's location; add --latlon LAT LON")

    if args.plot:
        if args.plot == "hodograph":
            given = [name for name in SOUNDING_ONLY_OPTIONS if _given(args, name)]
            if args.mode is not None and not args.params:
                given.append("mode")
            if given:
                raise CLIError(f"{_flag(given[0])} only applies to --plot sounding" + (" and --params" if given[0] == "mode" else ""))
        if args.plot == "sounding" and args.mode == "winter":
            given = [name for name in SEVERE_ONLY_OPTIONS if _given(args, name)]
            if given:
                raise CLIError(f"{_flag(given[0])} doesn't apply to winter soundings")
        if args.obs and args.mode not in ("winter", "auto"):
            raise CLIError("--obs only applies to winter soundings; add --mode winter or --mode auto")
        if args.map_zoom == 0 and (_given(args, "radar") or _given(args, "radar_time")):
            raise CLIError(f"{_flag('radar' if _given(args, 'radar') else 'radar_time')} needs the map, which --map-zoom 0 hides")
        if args.plot_file:
            args.plot_file = _plot_path(args.plot_file)
        elif args.json:
            raise CLIError(
                "--json cannot be combined with an interactive plot. "
                "Use --plot-file to save the plot instead."
            )
        elif not _has_display():
            raise CLIError(
                "No display is available to open the plot. "
                "Use --plot-file to save it instead."
            )

    for name, default in PLOT_DEFAULTS.items():
        if getattr(args, name) is None:
            setattr(args, name, default)
    args.storm_motion = _storm_motion(args.storm_motion)
    args.modify_sfc = _modify_sfc(args.modify_sfc)
    args.special_parcels = _special_parcels(args.parcels, args.background_parcels)


@contextmanager
def _hide_warnings():
    """Hide library warnings and Matplotlib's log messages (shown with --verbose).

    SounderPy's own warnings about the figure or data (``SounderPyWarning``, e.g. radar that isn't available)
    are still shown, as one line each.
    """
    mpl_logger = logging.getLogger("matplotlib")
    level = mpl_logger.level
    mpl_logger.setLevel(logging.ERROR)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("ignore")
            warnings.simplefilter("always", SounderPyWarning)
            yield
    finally:
        mpl_logger.setLevel(level)
    for warning in caught:
        if issubclass(warning.category, SounderPyWarning):
            print(f"sounderpy: warning: {warning.message}", file=sys.stderr)


def _run_quietly(
    func: Callable[..., Any],
    *args: Any,
    verbose: bool = False,
    **kwargs: Any,
) -> Any:
    """Run a SounderPy call while hiding library progress output and warnings unless requested.

    ``verbose="stderr"`` shows the output on stderr, keeping stdout free for JSON.
    """
    if verbose == "stderr":
        with redirect_stdout(sys.stderr):
            return func(*args, **kwargs)
    if verbose:
        return func(*args, **kwargs)
    with redirect_stdout(io.StringIO()), _hide_warnings():
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
        try:  # some sources give int parts, others zero-padded strings; ACARS hours can be "HH:MM"
            year, month, day = f"{int(year):04d}", f"{int(month):02d}", f"{int(day):02d}"
            hour = f"{int(hour):02d}" if str(hour).isdigit() else hour
        except (TypeError, ValueError):
            pass
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


def _is_placeholder(value: Any) -> bool:
    """Data sources fill in unknown site info with placeholders such as 'none' or 'no-site-location'."""
    return value is None or (isinstance(value, str) and (value in ("", "none") or value.startswith("no-")))


def _format_latlon(value: Any) -> Any:
    try:
        lat, lon = value
        return f"{float(lat):.2f}, {float(lon):.2f}"
    except (TypeError, ValueError):
        return value


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
        ("Lat/Lon", _format_latlon(site.get("site-latlon"))),
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
        if not _is_placeholder(value):
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
    if not args.plot:
        return

    # options were checked, and --plot-file given its extension, by _validate_options
    save = args.plot_file is not None
    filename = args.plot_file
    if save:
        Path(filename).parent.mkdir(parents=True, exist_ok=True)

    common = {
        "dark_mode": args.dark_mode,
        "map_zoom": args.map_zoom,
        "dpi": args.dpi,
        "save": save,
        "storm_motion": args.storm_motion,
        "modify_sfc": args.modify_sfc,
        "radar": args.radar,
        "radar_time": args.radar_time,
    }

    if args.plot == "sounding":
        kwargs = {
            **common,
            "color_blind": args.color_blind,
            "show_theta": args.show_theta,
            "sr_hodo": args.storm_relative,
            "mode": args.mode,
            "show_obs": args.obs,
            "special_parcels": args.special_parcels,
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


def _calc_params(data: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    """sounding_params for --params, or winter_params when the sounding is a winter one (--mode winter or auto)."""
    from .calc import sounding_params
    from .calc_winter import winter_params
    from .plot import _resolve_sounding_mode

    winter = args.mode == "winter" or (
        args.mode == "auto" and _resolve_sounding_mode("auto", data, args.modify_sfc) == "winter"
    )
    calc = winter_params if winter else sounding_params
    return calc(data, storm_motion=args.storm_motion, modify_sfc=args.modify_sfc).calc()


def _number(value: Any, unit: str = "", decimals: int = 0) -> str:
    """A parameter for the readout: rounded, with its unit, or NA when it's missing."""
    value = getattr(value, "magnitude", value)
    try:
        if np.ma.is_masked(value) or not np.isfinite(float(value)):
            return "NA"
    except (TypeError, ValueError):
        return "NA"
    return f"{float(value):.{decimals}f}{unit}"


def _print_params(params: dict[str, Any], data: dict[str, Any], args: argparse.Namespace) -> None:
    """Print the key calculated parameters under the profile summary."""
    from .plot.components import _describe_surface_modification, _format_motion
    from .plot.style import PTYPE_NAMES

    general, thermo, kinem = params["general"], params["thermo"], params["kinem"]
    jkg, m = " J/kg", " m"

    print()
    print("Parameters")
    print("----------")
    if args.modify_sfc:
        print(f"Modified surface: {_describe_surface_modification(data, args.modify_sfc)}")
    for key in ("sb", "ml", "mu"):
        print(f"{key.upper()} parcel: CAPE {_number(thermo[f'{key}cape'], jkg)} | CIN {_number(thermo[f'{key}cin'], jkg)} | "
              f"LCL {_number(thermo[f'{key}_lcl_z'], m)} | ECAPE {_number(thermo[f'{key}_ecape'], jkg)}")
    print(f"MU 0-3 km CAPE: {_number(thermo['mu3cape'], jkg)} | DCAPE: {_number(thermo['dcape'], jkg)} | "
          f"PWAT: {_number(general['pwat'], ' in', 2)}")
    print(f"Lapse rates: 0-3 km {_number(thermo['lr_03km'], '', 1)} | 3-6 km {_number(thermo['lr_36km'], '', 1)} | "
          f"700-500 hPa {_number(thermo['lr_75'], '', 1)} C/km")
    print(f"Freezing level: {_number(general['frz_pt_z'], m)} | Wet-bulb zero: {_number(general['wb_frz_pt_z'], m)}")

    motion = args.storm_motion if isinstance(args.storm_motion, str) else "custom"
    print(_format_motion(f"Storm motion ({motion.replace('_', ' ')})", kinem["sm_u"], kinem["sm_v"]))
    print(f"SRH: 0-500 m {_number(kinem['srh_0_to_500'])} | 0-1 km {_number(kinem['srh_0_to_1000'])} | "
          f"0-3 km {_number(kinem['srh_0_to_3000'])} | EIL {_number(kinem['srh_eil'])} m2/s2")
    print(f"Bulk shear: 0-1 km {_number(kinem['shear_0_to_1000'])} | 0-3 km {_number(kinem['shear_0_to_3000'])} | "
          f"0-6 km {_number(kinem['shear_0_to_6000'])} kt")
    print(f"STP: {_number(kinem['eil_stp'], '', 1)} | SCP: {_number(kinem['eil_scp'], '', 1)} | "
          f"SHIP: {_number(thermo['ship'], '', 1)} | 0-1 km EHI: {_number(kinem['ehi_0_to_1000'], '', 1)}")

    winter = params.get("winter")
    if winter:
        ptype = PTYPE_NAMES.get(winter["ptype"], winter["ptype"]) or "none expected"
        print(f"Precip type: {ptype} (SPC: {winter['sharppy_ptype'] or 'NA'})")
        phases = " -> ".join(step["phase"] for step in winter["phase_trace"])
        print(f"Phase from source to surface: {phases or 'NA'}")
        slr = "NA" if np.ma.is_masked(winter["slr"]) else f"{winter['slr']:.0f}:1"
        dgz = "NA" if np.ma.is_masked(winter["dgz_depth"]) else f"{winter['dgz_z'][0] / 1000:.1f}-{winter['dgz_z'][1] / 1000:.1f} km"
        print(f"Snow ratio: {slr} | DGZ: {dgz} | Freezing drizzle: {winter['fzdz']['level']} | Freezing fog: {winter['fzfg']['level']}")
        print(f"Surface: T {_number(winter['sfc_T_F'], ' F')} | wind chill {_number(winter['wind_chill'], ' F')}")


def _handle_profile_result(data: dict[str, Any], args: argparse.Namespace) -> int:
    if not isinstance(data, dict):
        raise CLIError(
            f"SounderPy returned {type(data).__name__}, but the CLI expected cleaned profile data."
        )
    if args.output:
        _export_profile(data, args.output, args.format, args.verbose)
    params = _run_quietly(_calc_params, data, args, verbose=args.verbose) if args.params else None
    _plot_profile(data, args)
    if args.json:
        output = data if params is None else {**data, "params": params}
        print(json.dumps(_json_ready(output), indent=2, allow_nan=False))
    else:
        _print_summary(data)
        if params is not None:
            _print_params(params, data, args)
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
    group.add_argument(
        "--params",
        action="store_true",
        help="Print key calculated parameters (CAPE, SRH, shear, ...; precip type with --mode winter). "
             "With --json, add them as \"params\".",
    )

    calc = parser.add_argument_group("calculation (used by --plot and --params)")
    calc.add_argument(
        "--mode",
        choices=("severe", "winter", "auto"),
        help=(
            "severe, winter, or auto (winter when the surface wet-bulb is below 35 F): the sounding plot, "
            "and whether --params includes winter values. Default: severe."
        ),
    )
    calc.add_argument(
        "--storm-motion",
        nargs="+",
        metavar="MOTION",
        help="right_moving, left_moving, mean_wind, or a direction (deg) and speed (kt) such as 270 25. "
             "Default: right_moving.",
    )
    calc.add_argument(
        "--modify-sfc",
        nargs="+",
        metavar="KEY=VALUE",
        help="Change the surface: T and Td (C), ws (kt) with wd (deg), method (mix, replace, blend), and depth "
             "(m, for blend), e.g. T=30 Td=21 or T=-6 method=blend depth=300. Exports (-o) keep the original surface.",
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
        help=(
            "Save --plot to this file instead of opening it in a window. The extension sets the "
            "format (.png, .pdf, .svg, ...); .png is added when there is none."
        ),
    )
    plot.add_argument(
        "--dark-mode", action="store_true", help="Use SounderPy's dark plot theme."
    )
    plot.add_argument(
        "--map-zoom",
        type=_nonnegative_int,
        metavar="N",
        help="Map inset zoom level; use 0 to disable the map/radar inset. Default: 2.",
    )
    plot.add_argument(
        "--dpi",
        type=_positive_int,
        help="Plot resolution in dots per inch. Default: 100.",
    )
    plot.add_argument(
        "--color-blind",
        action="store_true",
        help="Sounding plots: use the color-deficiency-friendly colors.",
    )
    plot.add_argument(
        "--show-theta",
        action="store_true",
        help="Sounding plots: show theta/theta-e.",
    )
    plot.add_argument(
        "--storm-relative",
        action="store_true",
        help="Use a storm-relative hodograph (on its own or on the sounding).",
    )
    plot.add_argument(
        "--radar",
        choices=("mosaic", "single", "none"),
        help="Map radar: the CONUS mosaic (recent dates only), the nearest NEXRAD site, or none. Default: mosaic.",
    )
    plot.add_argument(
        "--radar-time",
        choices=("sounding", "now"),
        help="Radar time: the profile's valid time, or now. Default: sounding.",
    )
    plot.add_argument(
        "--parcels",
        nargs="+",
        metavar="CODE",
        help="Severe sounding plots: highlighted special parcels such as mu_ia_ecape sb_ps_cape, or 'simple' for "
             "only the standard parcels (no MUECAPE trace, faster).",
    )
    plot.add_argument(
        "--background-parcels",
        nargs="+",
        metavar="CODE",
        help="Severe sounding plots: special parcels drawn faded behind --parcels.",
    )
    plot.add_argument(
        "--obs",
        action="store_true",
        help="Winter sounding plots (--mode winter or auto): add the nearest surface observation (IEM archive).",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Show SounderPy's progress output and warnings, and the full traceback on errors.",
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
            "  sounderpy bufkit hrrr KGFK 1 --run 2024-01-24 12 --plot sounding --mode winter --obs\n"
            "  sounderpy bufkit hrrr OUN 3 --params --storm-motion 250 30\n"
            "  sounderpy file oax.snd --plot sounding --plot-file oax.png\n"
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
    model.add_argument(
        "model",
        help="Model: rap-ruc (or rap, ruc), era5 (needs a Copernicus CDS API key in ~/.cdsapirc), or ncep.",
    )
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
    bufkit.add_argument(
        "model",
        help="BUFKIT model: gfs, rap, hrrr, or rrfs (nam and namnest for archived --run dates only).",
    )
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

    file_cmd = subparsers.add_parser(
        "file",
        help="Read a local SHARPpy-format sounding file.",
        description=(
            "Read a local SHARPpy-format sounding file (such as one saved with -o profile.snd). "
            "The file has no location, so --plot and --params need --latlon."
        ),
    )
    file_cmd.add_argument("path", help="SHARPpy-format sounding file.")
    file_cmd.add_argument(
        "--latlon",
        nargs=2,
        type=float,
        metavar=("LAT", "LON"),
        help="The profile's location, needed for --plot and --params (the map, and temperature advection).",
    )
    _add_profile_output_options(file_cmd)

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
    acars_get.add_argument(
        "profile", help="ACARS profile ID from 'sounderpy acars list' for the same date and hour, e.g. BNA_1820."
    )
    _add_profile_output_options(acars_get)

    return parser


def _read_file(args: argparse.Namespace) -> dict[str, Any]:
    """Read a SHARPpy file, filling in the site details and titles the plots need."""
    path = Path(args.path).expanduser()
    if not path.is_file():
        raise CLIError(f"No such file: {path}")
    if args.latlon is not None:
        lat, lon = args.latlon
        if not -90 <= lat <= 90 or not -180 <= lon <= 180:
            raise CLIError("--latlon needs a latitude from -90 to 90 and a longitude from -180 to 180")

    data = _run_quietly(from_sharppy, path, verbose=args.verbose)
    site = data["site_info"]
    site_id = str(site.get("site-id") or path.stem)
    site.update({"site-id": site_id, "site-name": site.get("site-name", path.name), "source": "SHARPPY FILE"})
    # SHARPpy heights are above sea level, so the lowest level is the station elevation
    site.setdefault("site-elv", float(data["z"][0].to("m").m))
    if args.latlon is not None:
        site["site-latlon"] = [float(args.latlon[0]), float(args.latlon[1])]
    time = site.get("valid-time")
    valid = (f"{int(time[1]):02d}/{int(time[2]):02d}/{int(time[0])} {int(time[3]):02d}Z"
             if isinstance(time, (list, tuple)) and len(time) >= 4 else "UNKNOWN TIME")
    where = site_id if args.latlon is None else f"{site_id} | {args.latlon[0]:g}, {args.latlon[1]:g}"
    data.setdefault("titles", {"top_title": f"SHARPPY FILE PROFILE | {path.name}", "left_title": f"VALID: {valid}",
                               "right_title": where})
    return data


def _run_command(args: argparse.Namespace) -> int:
    if getattr(args, "json", False) and args.verbose:
        args.verbose = "stderr"
    _validate_options(args)

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

    if args.command == "file":
        return _handle_profile_result(_read_file(args), args)

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
        if getattr(args, "verbose", False):
            traceback.print_exc()
        print(f"sounderpy: error: {str(exc) or type(exc).__name__}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

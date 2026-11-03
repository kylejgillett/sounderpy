#!/usr/bin/env python3
"""
Generate the example figures in docs/source/_static/images/.

Run this script from the repository root after installing SounderPy in editable
mode:

    pip install -e .
    python docs/scripts/make_docs_figures.py

Examples
--------
Generate the figures the documentation uses:

    python docs/scripts/make_docs_figures.py

Generate every figure, including ones no page currently shows:

    python docs/scripts/make_docs_figures.py --all

Generate specific figures, or write them somewhere else:

    python docs/scripts/make_docs_figures.py \
        --figures example-hodograph_light \
        --output /tmp/figures

Notes
-----
* Every figure requires network access.
* The gallery and tutorial figures are made by make_gallery_figures.py.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

import sounderpy as spy


DEFAULT_OUTPUT = Path("docs/source/_static/images")
DPI = 200


def obs(*args):
    return lambda: spy.get_obs_data(*args)


def bufkit(*args):
    return lambda: spy.get_bufkit_data(*args)


def acars(date_args, profile):
    return lambda: spy.acars_data(*date_args).get_profile(profile)


def reanalysis(*args):
    return lambda: spy.get_model_data(*args)


def composite_light():
    return [spy.get_obs_data("oax", "2014", "06", day, hour) for day, hour in (("16", "12"), ("16", "18"), ("17", "00"))]


def composite_dark():
    return [spy.get_bufkit_data("hrrr", "dtx", 0, "2024", "02", "28", hour, hush=True)
            for hour in ("00", "01", "02", "03", "04", "05", "06")]


def merged_profile():
    # Mirrors the merge_profiles() example in reference/helpertools.rst.
    profile_a = spy.get_obs_data("OAX", "2014", "06", "16", "18")
    profile_b = spy.get_obs_data("OAX", "2014", "06", "17", "00")
    merged = spy.merge_profiles(profile_a, profile_b, weight_a=0.6)
    merged["site_info"] = {
        "site-id": "KOAX", "site-name": "OMAHA VALLEY", "site-lctn": "NE US",
        "site-latlon": [41.32, -96.37], "site-elv": 350, "source": "RAOB OBSERVED PROFILE",
        "model": "none", "fcst-hour": "none", "run-time": "none", "valid-time": "none", "box_area": "none",
    }
    merged["titles"] = {
        "top_title": "COMPOSITE VERTICAL PROFILE",
        "left_title": "Profile-a (60%): 18z OAX 06/16/2014 | Profile-b (40%): 00z OAX 06/16/2014",
        "right_title": "KOAX, OMAHA VALLEY, NE | 41.32, -96.37    ",
    }
    return merged


def smoothed_profile():
    # Mirrors the smooth_profile() example in reference/helpertools.rst.
    return spy.smooth_profile(spy.get_obs_data("OAX", "2014", "06", "16", "18"), 3)


# name: (plot function, data retriever, plot kwargs, shown in the documentation)
FIGURES = {
    "example-sounding_light": (
        "sounding", obs("OAX", "2014", "06", "16", "18"),
        {"color_blind": True, "radar": "single"}, True,
    ),
    "example-sounding_light2": (
        "sounding", bufkit("HRRR", "KMLI", 8, "2023", "03", "31", "12"),
        {"color_blind": False, "radar": "single"}, False,
    ),
    "example-sounding_light3": (
        "sounding", acars(("2024", "04", "28", "03"), "TUL_0320"),
        {"color_blind": False, "radar": "single"}, False,
    ),
    "example-sounding_light4": (
        "sounding", reanalysis("rap", [33.19, -87.46], "2011", "04", "27", "22"),
        {"color_blind": True}, True,
    ),
    "example-sounding_dark": (
        "sounding", obs("OUN", "2024", "05", "07", "03"),
        {"dark_mode": True, "color_blind": True, "radar": "single"}, False,
    ),
    "example-hodograph_light": (
        "hodograph", obs("DTX", "2015", "06", "23", "00"),
        {"radar": "single"}, True,
    ),
    "example-hodograph_dark": (
        "hodograph", obs("ILN", "2024", "03", "15", "00"),
        {"dark_mode": True, "radar": "single"}, False,
    ),
    "example-sounding_dark3": (
        "sounding", obs("OUN", "1999", "05", "04", "00"),
        {"dark_mode": True}, True,
    ),
    "example-composite_light": (
        "composite", composite_light, {}, True,
    ),
    "example-composite_dark": (
        "composite", composite_dark,
        {"dark_mode": True, "linewidths": [4] * 7}, True,
    ),
    "merged_profile_example": (
        "sounding", merged_profile,
        {"radar": None, "special_parcels": "simple"}, True,
    ),
    "smoothed_profile_example": (
        "sounding", smoothed_profile,
        {"radar": None, "special_parcels": "simple"}, True,
    ),
    "example-goofy-acars": (
        "sounding", acars(("2019", "10", "21", "02"), "SMF_0210"),
        {}, True,
    ),
    "test-figure_abr-2024-08-29-00z": (
        "sounding", obs("kabr", "2024", "08", "29", "00"),
        {"color_blind": True, "radar": "single"}, True,
    ),
}

BUILDERS = {
    "sounding": spy.build_sounding,
    "hodograph": spy.build_hodograph,
    "composite": spy.build_composite,
}


def make_figure(name: str, output: Path) -> None:
    kind, get_data, kwargs, _ = FIGURES[name]
    BUILDERS[kind](get_data(), save=True, filename=str(output / f"{name}.png"), dpi=DPI, **kwargs)
    plt.close("all")


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help=f"Output directory. Default: {DEFAULT_OUTPUT}")
    parser.add_argument("--figures", nargs="+", choices=sorted(FIGURES), help="Generate only these figures.")
    parser.add_argument("--all", action="store_true", help="Also generate figures no documentation page shows.")
    parser.add_argument("--strict", action="store_true", help="Stop at the first failure.")
    return parser


def main() -> int:
    args = make_parser().parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    names = args.figures or [name for name, figure in FIGURES.items() if args.all or figure[3]]
    failed = []
    for name in names:
        print(f"Generating {name}")
        try:
            make_figure(name, args.output)
        except Exception:
            if args.strict:
                raise
            traceback.print_exc()
            failed.append(name)

    print(f"{len(names) - len(failed)} of {len(names)} figures written to {args.output}")
    if failed:
        print("Failed:", ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Generate SounderPy figures for the CM1 surface heating sensitivity study.

Run from any directory (the default output folder is found from the
script's own location):

    python docs/scripts/make_cm1_case_study.py

Outputs
-------
docs/source/case_studies/cm1_surface_heating/

    cm1_control.png
    cm1_sensitivity.png

Study design
------------
* control: the CM1 input_sounding in the repository's test data
* +2, +4, +6 °C: the control with a warmer surface, mixed upward along its
  dry adiabat ('method': 'mix'), the surface dewpoint held so its mixing
  ratio mixes through the deeper layer
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt

import sounderpy as spy


# Resolved from this file, so the script works from any directory.
DOCS_DIR = Path(__file__).resolve().parent.parent
OUTPUT_DIR = DOCS_DIR / "source" / "case_studies" / "cm1_surface_heating"

# the CM1 input_sounding in the repository's test data, and the metadata the case study gives it
CM1_SOUNDING = DOCS_DIR.parent / "test_data" / "input_sounding"
CM1_METADATA = {"latlon": [35.2, -97.4], "elev": 380, "top_title": "CM1 INPUT SOUNDING",
                "left_title": "CONTROL BASE STATE", "right_title": "CM1"}

# the control, then each run, in colors that are easy to tell apart (Okabe & Ito, colorblind-safe)
RUN_COLORS = ["black", "#0072B2", "#009E73", "#D55E00"]

DPI = 200


def banner(text: str) -> None:
    print()
    print("=" * 79)
    print(text)
    print("=" * 79)


def close_figures() -> None:
    plt.close("all")


def get_runs() -> dict[str, dict]:
    """The control and each warmer surface, the same calls as the case study page."""
    control = spy.from_cm1(str(CM1_SOUNDING), CM1_METADATA)
    T_sfc, Td_sfc = control["T"][0].to("degC").m, control["Td"][0].to("degC").m
    runs = {"control": control}
    for warming in (2, 4, 6):
        runs[f"+{warming} °C"] = spy.modify_surface(control, {"T": T_sfc + warming, "Td": Td_sfc, "method": "mix"})
    return runs


def print_parameters(runs: dict[str, dict]) -> None:
    """The comparison table printed on the case study page."""
    for name, profile in runs.items():
        params = spy.sounding_params(profile).calc()
        thermo, general = params["thermo"], params["general"]
        print(f"{name:8s}  MLLCL {thermo['ml_lcl_z']:5.0f} m   0-1 km RH {general['rh_0_1000']:3.0f} %   "
              f"PBL depth {thermo['pbl_depth']:5.0f} m   MLCAPE {thermo['mlcape']:5.0f}   MLCIN {thermo['mlcin']:4.0f}")


def save_control(control: dict, filename: Path) -> None:
    spy.build_sounding(control, map_zoom=0, save=True, filename=str(filename), dpi=100)
    close_figures()


def save_sensitivity(runs: dict[str, dict], filename: Path) -> None:
    fig, axes = spy.build_publication_sounding(
        list(runs.values()), labels=list(runs), panels=("skewt", "parameters"),
        colors=RUN_COLORS, ylim=(1000, 500), show_barbs=False,
        parameters=("mllcl", ("0–1 km RH", "general", "rh_0_1000", 1, "%", 0),
                    ("PBL depth", "thermo", "pbl_depth", 0.001, "km", 1), "mlcape", "mlcin",
                    ("DCAPE", "thermo", "dcape", 1, r"J kg$^{-1}$", 0)), dpi=DPI)
    fig.savefig(filename)
    close_figures()


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate SounderPy figures for the CM1 surface heating sensitivity study."
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIR,
        help=f"Output directory. Default: {OUTPUT_DIR}",
    )

    return parser


def main() -> int:
    args = make_parser().parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    runs = get_runs()

    banner("COMPARING THE RUNS")
    print_parameters(runs)

    banner("GENERATING FIGURES")
    save_control(runs["control"], out / "cm1_control.png")
    save_sensitivity(runs, out / "cm1_sensitivity.png")

    banner("CASE STUDY COMPLETE")
    print(out.resolve())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

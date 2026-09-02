#!/usr/bin/env python3
"""
Generate SounderPy figures for the 28 August 2024 Mound City, SD tornado
case study.

Run from the SounderPy repository root:

    python docs/scripts/make_mound_city_case_study.py

Outputs
-------
docs/source/_static/case_studies/mound_city_20240828/

    header.jpg
    abr_12z_sounding.png
    abr_00z_sounding.png
    hrrr_12z_f13_sounding.png
    rap_01z_sounding.png
    hrrr_12z_f13_hodograph.png
    abr_00z_hodograph.png
    rap_01z_hodograph.png
    environment_composite.png
    parameter_comparison.csv

Case design
-----------
* KABR 12Z 28 Aug: morning observed regional environment
* HRRR 12Z run / KMBG F13: morning forecast valid 01Z 29 Aug
* KABR 00Z 29 Aug: evening observed regional environment
* RAP 01Z 29 Aug near 45.68, -100.18: local analysis near tornadogenesis

The first EF2 tornado began around 00:50 UTC 29 August 2024, making the
01 UTC model profiles a useful event-time comparison.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from urllib.request import Request, urlopen

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

import sounderpy as spy


OUTPUT_DIR = Path(
    "../source/case_studies/mound_city_20240828"
)

HEADER_URL = (
    "https://www.weather.gov/images/abr/EventSummaries/"
    "Thunderstorms_Tornadoes/20240828/"
    "PTU%20Kyle%20Gillett%20Twitter%20Messages%20%282%29.jpg"
)

EVENT_LATLON = [45.68, -100.18]

DPI = 160


def banner(text: str) -> None:
    print()
    print("=" * 79)
    print(text)
    print("=" * 79)


def close_figures() -> None:
    plt.close("all")


def download_header(output: Path) -> None:
    """
    Download the Kyle Gillett header photograph from the NWS event page.

    The docs should ultimately keep a local copy so normal Sphinx builds do
    not depend on weather.gov.
    """
    if output.exists():
        print(f"[case study] Header already exists: {output}")
        return

    output.parent.mkdir(parents=True, exist_ok=True)

    request = Request(
        HEADER_URL,
        headers={
            "User-Agent": (
                "SounderPy documentation builder "
                "(https://github.com/kylejgillett/sounderpy)"
            )
        },
    )

    try:
        with urlopen(request, timeout=30) as response:
            output.write_bytes(response.read())
        print(f"[case study] Downloaded header: {output}")
    except Exception as exc:
        print("[case study] Could not download the header photograph.")
        print(f"             {type(exc).__name__}: {exc}")
        print()
        print("Place your original image manually at:")
        print(f"    {output}")
        print()
        print("The case-study figures can still be generated.")


def get_profiles() -> dict[str, dict]:
    profiles: dict[str, dict] = {}

    banner("RETRIEVING 12Z ABR OBSERVATION")
    profiles["abr_12z"] = spy.get_obs_data(
        "KABR",
        "2024", "08", "28", "12",
        hush=True,
    )

    banner("RETRIEVING 00Z ABR OBSERVATION")
    profiles["abr_00z"] = spy.get_obs_data(
        "KABR",
        "2024", "08", "29", "00",
        hush=True,
    )

    banner("RETRIEVING 12Z HRRR F13 KMBG FORECAST")
    try:
        profiles["hrrr_f13"] = spy.get_bufkit_data(
            "hrrr",
            "KMBG",
            13,
            "2024", "08", "28", "12",
            hush=True,
        )
    except (Exception, SystemExit) as exc:
        print(
            "[case study] KMBG HRRR BUFKIT retrieval failed; "
            "trying KABR as a regional fallback."
        )
        print(f"             {type(exc).__name__}: {exc}")

        profiles["hrrr_f13"] = spy.get_bufkit_data(
            "hrrr",
            "KABR",
            13,
            "2024", "08", "28", "12",
            hush=True,
        )

    banner("RETRIEVING 01Z RAP ANALYSIS NEAR MOUND CITY")
    profiles["rap_01z"] = spy.get_model_data(
        "rap-ruc",
        EVENT_LATLON,
        "2024", "08", "29", "01",
        box_avg_size=0.10,
        hush=True,
    )

    return profiles


def save_sounding(data: dict, filename: Path) -> None:
    print(f"[case study] Writing {filename.name}")
    spy.build_sounding(
        data,
        radar=None,
        map_zoom=0,
        special_parcels="simple",
        dpi=DPI,
        save=True,
        filename=str(filename),
    )
    close_figures()


def save_hodograph(data: dict, filename: Path) -> None:
    print(f"[case study] Writing {filename.name}")
    spy.build_hodograph(
        data,
        radar=None,
        map_zoom=0,
        dpi=DPI,
        save=True,
        filename=str(filename),
    )
    close_figures()


def save_composite(profiles: dict[str, dict], filename: Path) -> None:
    print(f"[case study] Writing {filename.name}")

    data_list = [
        profiles["abr_12z"],
        profiles["hrrr_f13"],
        profiles["abr_00z"],
        profiles["rap_01z"],
    ]

    spy.build_composite(
        data_list,
        dark_mode=True,
        shade_between=False,
        colors_to_use=[
            "tab:blue",
            "tab:orange",
            "tab:green",
            "tab:red",
        ],
        lw_to_use=[
            2,
            2,
            3,
            3,
        ],
        alphas_to_use=[
            0.75,
            0.75,
            1.0,
            1.0,
        ],
        save=True,
        filename=str(filename),
    )
    close_figures()


def scalar(value):
    """Convert SounderPy/NumPy/masked values to CSV-friendly scalars."""
    if ma.is_masked(value):
        return ""

    try:
        value = value.magnitude
    except Exception:
        pass

    try:
        arr = np.asarray(value)
        if arr.size == 1:
            value = arr.item()
    except Exception:
        pass

    try:
        value = float(value)
    except Exception:
        return str(value)

    if not math.isfinite(value):
        return ""

    return round(value, 2)


def profile_parameters(data: dict) -> dict[str, object]:
    general, thermo, kinem, _ = spy.sounding_params(
        data,
        storm_motion="right_moving",
    ).calc()

    return {
        "SBCAPE_Jkg": scalar(thermo.get("sbcape")),
        "MLCAPE_Jkg": scalar(thermo.get("mlcape")),
        "MUCAPE_Jkg": scalar(thermo.get("mucape")),
        "MU_ECAPE_Jkg": scalar(thermo.get("mu_ecape")),
        "SBCIN_Jkg": scalar(thermo.get("sbcin")),
        "MLCIN_Jkg": scalar(thermo.get("mlcin")),
        "MUCIN_Jkg": scalar(thermo.get("mucin")),
        "0_1km_SRH_m2s2": scalar(kinem.get("srh_0_to_1000")),
        "0_3km_SRH_m2s2": scalar(kinem.get("srh_0_to_3000")),
        "0_1km_shear_kt": scalar(kinem.get("shear_0_to_1000")),
        "0_6km_shear_kt": scalar(kinem.get("shear_0_to_6000")),
        "effective_STP": scalar(kinem.get("eil_stp")),
        "effective_SCP": scalar(kinem.get("eil_scp")),
        "0_1km_RH_pct": scalar(general.get("rh_0_1000")),
    }


def save_parameter_table(
    profiles: dict[str, dict],
    filename: Path,
) -> None:
    print(f"[case study] Writing {filename.name}")

    labels = {
        "abr_12z": "ABR 12Z observation",
        "hrrr_f13": "HRRR 12Z F13 forecast",
        "abr_00z": "ABR 00Z observation",
        "rap_01z": "RAP 01Z analysis",
    }

    rows = []

    for key in (
        "abr_12z",
        "hrrr_f13",
        "abr_00z",
        "rap_01z",
    ):
        values = profile_parameters(profiles[key])
        rows.append(
            {
                "profile": labels[key],
                **values,
            }
        )

    fieldnames = list(rows[0].keys())

    with filename.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)

    print()
    print("[case study] PARAMETER COMPARISON")
    print("-" * 79)

    for row in rows:
        print(row["profile"])
        print(
            "  MLCAPE:",
            row["MLCAPE_Jkg"],
            "J/kg | 0-1 SRH:",
            row["0_1km_SRH_m2s2"],
            "m2/s2 | 0-6 shear:",
            row["0_6km_shear_kt"],
            "kt",
        )


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Generate SounderPy figures for the 28 August 2024 "
            "Mound City tornado case study."
        )
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_DIR,
        help=f"Output directory. Default: {OUTPUT_DIR}",
    )

    parser.add_argument(
        "--skip-header",
        action="store_true",
        help="Do not attempt to download the NWS-hosted header photograph.",
    )

    return parser


def main() -> int:
    args = make_parser().parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    if not args.skip_header:
        download_header(out / "header.jpg")

    profiles = get_profiles()

    banner("GENERATING SOUNDINGS")
    save_sounding(
        profiles["abr_12z"],
        out / "abr_12z_sounding.png",
    )
    save_sounding(
        profiles["hrrr_f13"],
        out / "hrrr_12z_f13_sounding.png",
    )
    save_sounding(
        profiles["abr_00z"],
        out / "abr_00z_sounding.png",
    )
    save_sounding(
        profiles["rap_01z"],
        out / "rap_01z_sounding.png",
    )

    banner("GENERATING HODOGRAPHS")
    save_hodograph(
        profiles["hrrr_f13"],
        out / "hrrr_12z_f13_hodograph.png",
    )
    save_hodograph(
        profiles["abr_00z"],
        out / "abr_00z_hodograph.png",
    )
    save_hodograph(
        profiles["rap_01z"],
        out / "rap_01z_hodograph.png",
    )

    banner("GENERATING COMPOSITE")
    save_composite(
        profiles,
        out / "environment_composite.png",
    )

    banner("CALCULATING COMPARISON PARAMETERS")
    save_parameter_table(
        profiles,
        out / "parameter_comparison.csv",
    )

    banner("CASE STUDY COMPLETE")
    print(out.resolve())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

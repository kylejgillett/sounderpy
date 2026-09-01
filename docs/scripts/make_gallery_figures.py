#!/usr/bin/env python3
"""
Generate SounderPy documentation/gallery figures.

Run this script from the repository root after installing SounderPy in editable
mode:

    pip install -e .
    python docs/scripts/make_gallery_figures.py

By default, figures are written to:

    docs/source/_static/gallery/

The script is organized into tutorial-sized groups so that the same source file
can be used to regenerate the images shown throughout the documentation.

Examples
--------
Generate everything:

    python docs/scripts/make_gallery_figures.py

Generate only theme/accessibility examples:

    python docs/scripts/make_gallery_figures.py --sections themes accessibility

Generate only radar/map examples:

    python docs/scripts/make_gallery_figures.py --sections map-radar

Use a specific ACARS profile:

    python docs/scripts/make_gallery_figures.py \
        --sections data-sources \
        --acars-profile BNA_2320

Notes
-----
* Most figures use the same OAX 18Z 16 June 2014 RAOB so visual differences
  reflect plotting settings rather than different meteorological environments.
* The archived single-site radar example also uses the OAX/Pilger case.
* The mosaic radar example dynamically searches for a recent OAX RAOB because
  mosaic reflectivity is only available for recent dates.
* RAP/RUC, BUFKIT, ACARS, radar, and recent-observation examples require network
  access.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Callable

import matplotlib.pyplot as plt

import sounderpy as spy


# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_OUTPUT = Path("docs/source/_static/gallery")
DPI = 150

# A deep, severe-weather sounding that works well for theme, parcel, and
# hodograph demonstrations.
BASE_OBS = {
    "station": "OAX",
    "year": "2014",
    "month": "06",
    "day": "16",
    "hour": "18",
}

# Known RAP/RUC retrieval case.
RAP_EXAMPLE = {
    "model": "rap-ruc",
    "latlon": [44.58, -100.82],
    "year": "2024",
    "month": "08",
    "day": "28",
    "hour": "18",
    "box_avg_size": 0.25,
}

# Known archived BUFKIT example.
BUFKIT_EXAMPLE = {
    "model": "gfs",
    "station": "KMOP",
    "fcst_hour": 6,
    "run_year": "2023",
    "run_month": "08",
    "run_day": "05",
    "run_hour": "12",
}

# ACARS is intentionally optional because profile IDs depend on the archive.
# The script can select the first available profile automatically, or you can
# pass --acars-profile explicitly.
ACARS_EXAMPLE = {
    "year": "2024",
    "month": "05",
    "day": "21",
    "hour": "18",
}


# =============================================================================
# HELPERS
# =============================================================================

def banner(text: str) -> None:
    print()
    print("=" * 79)
    print(text)
    print("=" * 79)


def close_figures() -> None:
    """Close all Matplotlib figures after each gallery image."""
    plt.close("all")


def run_example(
    name: str,
    func: Callable[[], None],
    *,
    strict: bool = False,
) -> bool:
    """Run one gallery example without losing the rest if an external source fails."""
    print(f"\n[gallery] {name}")
    try:
        func()
    except (KeyboardInterrupt, SystemExit):
        if strict:
            raise
        print(f"[gallery] SKIPPED/FAILED: {name}")
        close_figures()
        return False
    except Exception as exc:
        if strict:
            raise
        print(f"[gallery] SKIPPED/FAILED: {name}")
        print(f"          {type(exc).__name__}: {exc}")
        close_figures()
        return False

    close_figures()
    print(f"[gallery] COMPLETE: {name}")
    return True


def save_sounding(
    data: dict,
    output: Path,
    **kwargs,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    spy.build_sounding(
        data,
        save=True,
        filename=str(output),
        dpi=DPI,
        **kwargs,
    )


def save_hodograph(
    data: dict,
    output: Path,
    **kwargs,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    spy.build_hodograph(
        data,
        save=True,
        filename=str(output),
        dpi=DPI,
        **kwargs,
    )


# =============================================================================
# DATA RETRIEVAL
# =============================================================================

def get_base_obs() -> dict:
    """Retrieve the common OAX/Pilger profile used by most tutorials."""
    return spy.get_obs_data(
        BASE_OBS["station"],
        BASE_OBS["year"],
        BASE_OBS["month"],
        BASE_OBS["day"],
        BASE_OBS["hour"],
        hush=True,
    )


def get_rap_example() -> dict:
    return spy.get_model_data(
        RAP_EXAMPLE["model"],
        RAP_EXAMPLE["latlon"],
        RAP_EXAMPLE["year"],
        RAP_EXAMPLE["month"],
        RAP_EXAMPLE["day"],
        RAP_EXAMPLE["hour"],
        box_avg_size=RAP_EXAMPLE["box_avg_size"],
        hush=True,
    )


def get_bufkit_example() -> dict:
    return spy.get_bufkit_data(
        BUFKIT_EXAMPLE["model"],
        BUFKIT_EXAMPLE["station"],
        BUFKIT_EXAMPLE["fcst_hour"],
        BUFKIT_EXAMPLE["run_year"],
        BUFKIT_EXAMPLE["run_month"],
        BUFKIT_EXAMPLE["run_day"],
        BUFKIT_EXAMPLE["run_hour"],
        hush=True,
    )


def get_acars_example(profile_id: str | None = None) -> tuple[dict, str]:
    """
    Retrieve one ACARS example.

    If no profile ID is provided, select the first profile returned by the
    archive for ACARS_EXAMPLE.
    """
    connection = spy.acars_data(
        ACARS_EXAMPLE["year"],
        ACARS_EXAMPLE["month"],
        ACARS_EXAMPLE["day"],
        ACARS_EXAMPLE["hour"],
    )

    if profile_id is None:
        profiles = connection.list_profiles()
        if not profiles:
            raise RuntimeError(
                "No ACARS profiles were returned for the configured date/hour."
            )
        profile_id = profiles[0]
        print(f"[gallery] Using ACARS profile: {profile_id}")

    data = connection.get_profile(profile_id, hush=True)
    return data, profile_id


def recent_synoptic_obs(
    station: str = "OAX",
    lookback_cycles: int = 6,
) -> dict:
    """
    Find a recent 00Z/12Z observed sounding for the mosaic-radar example.

    The search begins with the most recent nominal 00Z or 12Z cycle and walks
    backward. This keeps the radar mosaic example within the recent-data window
    without hard-coding a date that will immediately become stale.
    """
    now = datetime.now(timezone.utc)

    current_cycle_hour = 12 if now.hour >= 12 else 0
    candidate = now.replace(
        hour=current_cycle_hour,
        minute=0,
        second=0,
        microsecond=0,
    )

    # Give the latest RAOB a little time to reach upstream archives.
    if (now - candidate) < timedelta(hours=2):
        candidate -= timedelta(hours=12)

    errors: list[str] = []

    for _ in range(lookback_cycles):
        try:
            print(
                "[gallery] Trying recent RAOB:",
                station,
                candidate.strftime("%Y-%m-%d %HZ"),
            )
            return spy.get_obs_data(
                station,
                candidate.strftime("%Y"),
                candidate.strftime("%m"),
                candidate.strftime("%d"),
                candidate.strftime("%H"),
                hush=True,
            )
        except (Exception, SystemExit) as exc:
            errors.append(
                f"{candidate.strftime('%Y-%m-%d %HZ')}: "
                f"{type(exc).__name__}: {exc}"
            )
            candidate -= timedelta(hours=12)

    raise RuntimeError(
        "Could not retrieve a recent synoptic sounding for the mosaic example.\n"
        + "\n".join(errors)
    )


# =============================================================================
# TUTORIAL 1 — LIGHT VS DARK
# =============================================================================

def make_theme_examples(data: dict, out: Path, strict: bool) -> None:
    banner("THEME EXAMPLES — LIGHT VS DARK")

    run_example(
        "Light-mode sounding",
        lambda: save_sounding(
            data,
            out / "themes" / "sounding_light.png",
            dark_mode=False,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    run_example(
        "Dark-mode sounding",
        lambda: save_sounding(
            data,
            out / "themes" / "sounding_dark.png",
            dark_mode=True,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    run_example(
        "Light-mode hodograph",
        lambda: save_hodograph(
            data,
            out / "themes" / "hodograph_light.png",
            dark_mode=False,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    run_example(
        "Dark-mode hodograph",
        lambda: save_hodograph(
            data,
            out / "themes" / "hodograph_dark.png",
            dark_mode=True,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )


# =============================================================================
# TUTORIAL 2 — COLOR-BLIND / ACCESSIBILITY SETTINGS
# =============================================================================

def make_accessibility_examples(data: dict, out: Path, strict: bool) -> None:
    banner("ACCESSIBILITY EXAMPLES — STANDARD VS COLOR-BLIND")

    run_example(
        "Standard sounding colors",
        lambda: save_sounding(
            data,
            out / "accessibility" / "sounding_standard_colors.png",
            color_blind=False,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    run_example(
        "Color-deficiency-friendly sounding",
        lambda: save_sounding(
            data,
            out / "accessibility" / "sounding_colorblind.png",
            color_blind=True,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )


# =============================================================================
# TUTORIAL 3 — MAP AND RADAR SETTINGS
# =============================================================================

def make_map_radar_examples(data: dict, out: Path, strict: bool) -> None:
    banner("MAP / RADAR EXAMPLES")

    # 1. No inset at all.
    run_example(
        "No map / no radar",
        lambda: save_sounding(
            data,
            out / "map_radar" / "map_off.png",
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    # 2. Map inset, but no radar.
    run_example(
        "Map on / radar off",
        lambda: save_sounding(
            data,
            out / "map_radar" / "map_only.png",
            radar=None,
            map_zoom=2,
        ),
        strict=strict,
    )

    # 3. Archived single-site NEXRAD data for the actual sounding time.
    run_example(
        "Single-site radar",
        lambda: save_sounding(
            data,
            out / "map_radar" / "radar_single.png",
            radar="single",
            radar_time="sounding",
            map_zoom=2,
        ),
        strict=strict,
    )

    # 4. Show a wider map extent.
    run_example(
        "Larger map extent",
        lambda: save_sounding(
            data,
            out / "map_radar" / "map_zoom_4.png",
            radar=None,
            map_zoom=4,
        ),
        strict=strict,
    )

    # 5. Mosaic radar requires recent data. Retrieve a recent RAOB independently.
    def mosaic_example() -> None:
        recent = recent_synoptic_obs("OAX")
        save_sounding(
            recent,
            out / "map_radar" / "radar_mosaic.png",
            radar="mosaic",
            radar_time="sounding",
            map_zoom=2,
        )

    run_example(
        "Recent CONUS mosaic radar",
        mosaic_example,
        strict=strict,
    )


# =============================================================================
# TUTORIAL 4 — PARCEL SETTINGS
# =============================================================================

def make_parcel_examples(data: dict, out: Path, strict: bool) -> None:
    banner("PARCEL EXAMPLES")

    # Default SounderPy parcel behavior.
    run_example(
        "Default parcels",
        lambda: save_sounding(
            data,
            out / "parcels" / "parcels_default.png",
            special_parcels=None,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    # Faster/simple plot: common MU/ML/SB CAPE parcels only.
    run_example(
        "Simple parcels",
        lambda: save_sounding(
            data,
            out / "parcels" / "parcels_simple.png",
            special_parcels="simple",
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    # Exact custom example used by the SounderPy parcel documentation:
    # highlighted = surface-based irreversible-adiabatic ECAPE
    # background  = surface-based pseudoadiabatic ECAPE and CAPE
    custom_parcels = [
        ["sb_ia_ecape"],
        ["sb_ps_ecape", "sb_ps_cape"],
    ]

    run_example(
        "Custom highlighted/background parcels",
        lambda: save_sounding(
            data,
            out / "parcels" / "parcels_custom.png",
            special_parcels=custom_parcels,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )


# =============================================================================
# TUTORIAL 5 — OTHER USEFUL PLOT OPTIONS
# =============================================================================

def make_plot_option_examples(data: dict, out: Path, strict: bool) -> None:
    banner("OTHER PLOT OPTION EXAMPLES")

    run_example(
        "Theta / theta-e display",
        lambda: save_sounding(
            data,
            out / "options" / "show_theta.png",
            show_theta=True,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    run_example(
        "Ground-relative hodograph",
        lambda: save_hodograph(
            data,
            out / "options" / "hodograph_ground_relative.png",
            sr_hodo=False,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    run_example(
        "Storm-relative hodograph",
        lambda: save_hodograph(
            data,
            out / "options" / "hodograph_storm_relative.png",
            sr_hodo=True,
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    run_example(
        "Custom storm motion",
        lambda: save_hodograph(
            data,
            out / "options" / "hodograph_custom_storm_motion.png",
            storm_motion=[250.0, 45.0],
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )

    run_example(
        "Hodograph boundary",
        lambda: save_hodograph(
            data,
            out / "options" / "hodograph_boundary.png",
            hodo_boundary={
                "angle": [45],
                "color": ["tab:purple"],
            },
            radar=None,
            map_zoom=0,
        ),
        strict=strict,
    )


# =============================================================================
# TUTORIAL 6 — DIFFERENT DATA SOURCES
# =============================================================================

def make_data_source_examples(
    base_obs: dict,
    out: Path,
    strict: bool,
    acars_profile: str | None,
) -> None:
    banner("DATA-SOURCE EXAMPLES")

    # Observed RAOB
    run_example(
        "Observed RAOB sounding",
        lambda: save_sounding(
            base_obs,
            out / "data_sources" / "data_raob.png",
            radar=None,
            map_zoom=0,
            special_parcels="simple",
        ),
        strict=strict,
    )

    # RAP/RUC reanalysis
    def rap_plot() -> None:
        data = get_rap_example()
        save_sounding(
            data,
            out / "data_sources" / "data_rap_ruc.png",
            radar=None,
            map_zoom=0,
            special_parcels="simple",
        )

    run_example(
        "RAP/RUC reanalysis sounding",
        rap_plot,
        strict=strict,
    )

    # BUFKIT forecast
    def bufkit_plot() -> None:
        data = get_bufkit_example()
        save_sounding(
            data,
            out / "data_sources" / "data_bufkit.png",
            radar=None,
            map_zoom=0,
            special_parcels="simple",
        )

    run_example(
        "BUFKIT forecast sounding",
        bufkit_plot,
        strict=strict,
    )

    # ACARS observation
    def acars_plot() -> None:
        data, selected_profile = get_acars_example(acars_profile)
        print(f"[gallery] ACARS gallery profile: {selected_profile}")
        save_sounding(
            data,
            out / "data_sources" / "data_acars.png",
            radar=None,
            map_zoom=0,
            special_parcels="simple",
        )

    run_example(
        "ACARS aircraft sounding",
        acars_plot,
        strict=strict,
    )


# =============================================================================
# TUTORIAL 7 — COMPOSITE SOUNDINGS
# =============================================================================

def make_composite_examples(out: Path, strict: bool) -> None:
    banner("COMPOSITE SOUNDING EXAMPLES")

    def retrieve_profiles() -> list[dict]:
        profiles = []
        for year, month, day, hour in [
            ("2014", "06", "16", "12"),
            ("2014", "06", "16", "18"),
            ("2014", "06", "17", "00"),
        ]:
            profiles.append(
                spy.get_obs_data("OAX", year, month, day, hour, hush=True)
            )
        return profiles

    def light_composite() -> None:
        profiles = retrieve_profiles()
        output = out / "composites" / "composite_light.png"
        output.parent.mkdir(parents=True, exist_ok=True)

        spy.build_composite(
            profiles,
            shade_between=False,
            dark_mode=False,
            save=True,
            filename=str(output),
        )

    run_example(
        "Light composite — OAX temporal evolution",
        light_composite,
        strict=strict,
    )

    def dark_composite() -> None:
        profiles = retrieve_profiles()
        output = out / "composites" / "composite_dark.png"
        output.parent.mkdir(parents=True, exist_ok=True)

        spy.build_composite(
            profiles,
            shade_between=False,
            dark_mode=True,
            save=True,
            filename=str(output),
        )

    run_example(
        "Dark composite — OAX temporal evolution",
        dark_composite,
        strict=strict,
    )

    def shaded_composite() -> None:
        profiles = retrieve_profiles()
        output = out / "composites" / "composite_shaded.png"
        output.parent.mkdir(parents=True, exist_ok=True)

        spy.build_composite(
            profiles,
            shade_between=True,
            dark_mode=False,
            save=True,
            filename=str(output),
        )

    run_example(
        "Composite with T/Td shading",
        shaded_composite,
        strict=strict,
    )


# =============================================================================
# MAIN
# =============================================================================

VALID_SECTIONS = (
    "themes",
    "accessibility",
    "map-radar",
    "parcels",
    "options",
    "data-sources",
    "composites",
)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate SounderPy documentation gallery figures."
    )

    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Gallery output directory. Default: {DEFAULT_OUTPUT}",
    )

    parser.add_argument(
        "--sections",
        nargs="+",
        choices=VALID_SECTIONS,
        default=list(VALID_SECTIONS),
        help="Only generate selected tutorial/gallery sections.",
    )

    parser.add_argument(
        "--acars-profile",
        help=(
            "Specific ACARS profile ID for the ACARS gallery image. "
            "If omitted, the first available configured profile is used."
        ),
    )

    parser.add_argument(
        "--strict",
        action="store_true",
        help="Stop immediately if any gallery example fails.",
    )

    return parser


def main() -> int:
    args = make_parser().parse_args()

    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    print(f"Gallery output: {out.resolve()}")
    print("Sections:", ", ".join(args.sections))

    # Retrieve the shared profile once if any selected section needs it.
    needs_base = any(
        section in args.sections
        for section in (
            "themes",
            "accessibility",
            "map-radar",
            "parcels",
            "options",
            "data-sources",
        )
    )

    base_obs: dict | None = None
    if needs_base:
        banner("RETRIEVING COMMON OAX PROFILE")
        base_obs = get_base_obs()

    if "themes" in args.sections:
        make_theme_examples(base_obs, out, args.strict)

    if "accessibility" in args.sections:
        make_accessibility_examples(base_obs, out, args.strict)

    if "map-radar" in args.sections:
        make_map_radar_examples(base_obs, out, args.strict)

    if "parcels" in args.sections:
        make_parcel_examples(base_obs, out, args.strict)

    if "options" in args.sections:
        make_plot_option_examples(base_obs, out, args.strict)

    if "data-sources" in args.sections:
        make_data_source_examples(
            base_obs,
            out,
            args.strict,
            args.acars_profile,
        )

    if "composites" in args.sections:
        make_composite_examples(out, args.strict)

    banner("GALLERY GENERATION FINISHED")
    print(f"Images written beneath: {out.resolve()}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

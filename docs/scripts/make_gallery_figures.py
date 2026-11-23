#!/usr/bin/env python3
"""
Generate SounderPy documentation/gallery figures.

Install SounderPy in editable mode, then run this script from any directory
(the default output folder is found from the script's own location):

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

Generate only the winter soundings:

    python docs/scripts/make_gallery_figures.py --sections winter

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

# each profile is downloaded once per run, then copied for every figure that uses it
from profile_cache import get_acars_profile, get_bufkit_data, get_model_data, get_obs_data, summary


# =============================================================================
# CONFIGURATION
# =============================================================================

# Resolved from this file, so the script works from any directory.
DOCS_DIR = Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT = DOCS_DIR / "source" / "_static" / "gallery"
DPI = 200

# A deep, severe-weather sounding that works well for theme, parcel, and
# hodograph demonstrations.
BASE_OBS = {
    "station": "OAX",
    "year": "2014",
    "month": "06",
    "day": "16",
    "hour": "18",
}

# the gallery's custom theme example: a light parchment theme to sit beside the docs' warm, dark palette,
# drawn on the Mound City case (HRRR F03 for Ashley, ND, 21Z 28 August 2024)
PASTORAL_THEME = {
    "background": "#f5edda",        # parchment
    "text": "#241c15",              # dark-brown ink
    "border": "#3b2f24",
    "marker": "#241c15",
    "barb": "#241c15",
    "skew_line": "#3b2f24",
    "shade": "#a08858",             # isotherm shading
    "shade_alpha": 0.09,
    "adiabat": "#8a7348",
    "isotherm": "#5c4a30",
    "temperature": "#a3241a",       # brick
    "dewpoint": "#3c6a1c",          # olive
    "virtual": "#6b1a12",
    "wetbulb": "#2a6595",           # slate blue
    "parcel": {"sb": "#a3241a", "ml": "#cc620c", "mu": "#5a2a73"},
    "hodo_layers": ((0, 1000, "#5a2a73"), (1000, 3000, "#a3241a"), (3000, 6000, "#cc620c"),
                    (6000, 9000, "#c4950f"), (9000, 12000, "#d8b84a")),
    "warm": "#a3241a", "cool": "#2a6595", "eil": "#1f7a72", "dgz": "#2a5596", "hgz": "#3c6a1c",
    "pbl": "#5c4a30", "lapse_max": "#a3241a", "lapse": "#6a3214", "ncape": "#7a3418",
    "moisture": "#3c6a1c", "moisture_strong": "#2a4f0e",
    "kinematics": "#2a6595", "kinematics_strong": "#1a4a70",
    "kinematics_layers": ("#1f7a72", "#22708a", "#2a6595", "#3a5891", "#4a4b88"),
    "composite_index": "#a3241a", "ehi": "#c24e12", "cape": "#cc620c", "cin": "#1f7a72",
    "streamwise": "#cc620c", "theta": "#5a2a73", "low_fill": "#2a6595",
    "u_wind": "#3c6a1c", "v_wind": "#2a6595", "mcs": "#cc620c", "dtm": "#6a3214",
    "watermark": "#7a6440",
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
    return get_obs_data(
        BASE_OBS["station"],
        BASE_OBS["year"],
        BASE_OBS["month"],
        BASE_OBS["day"],
        BASE_OBS["hour"],
        hush=True,
    )


def get_rap_example() -> dict:
    return get_model_data(
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
    return get_bufkit_data(
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

    data = get_acars_profile(ACARS_EXAMPLE["year"], ACARS_EXAMPLE["month"], ACARS_EXAMPLE["day"],
                             ACARS_EXAMPLE["hour"], profile_id, hush=True)
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
            return get_obs_data(
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

    def custom_theme_example() -> None:
        kasy = get_bufkit_data("hrrr", "kasy", 3, "2024", "08", "28", "21", hush=True)
        save_sounding(
            kasy,
            out / "themes" / "sounding_custom_theme.png",
            theme=PASTORAL_THEME,
            radar="single",
            map_zoom=1,
        )

    run_example("Custom-theme sounding", custom_theme_example, strict=strict)


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
            map_zoom=1,
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
# WINTER SOUNDINGS
# =============================================================================

# the gallery's winter cases, each drawn with mode='winter', the nearest surface observation, and single-site
# radar, as the code under each image shows: (output file, how to get the profile, extra build_sounding settings)
WINTER_EXAMPLES = (
    ("winter_bufkit_fzdz.png", "HRRR freezing drizzle at Grand Forks, ND",
     lambda: get_bufkit_data("hrrr", "kgfk", 1, "2024", "01", "24", "12", hush=True), {}),
    ("winter_raob_ice.png", "Detroit, MI RAOB ice storm",
     lambda: get_obs_data("DTX", "2023", "02", "23", "00", hush=True), {}),
    ("winter_raob_dark.png", "Omaha, NE RAOB blizzard (dark mode)",
     lambda: get_obs_data("OAX", "2024", "01", "12", "12", hush=True), {"dark_mode": True}),
    ("winter_acars.png", "Denver, CO ACARS profile",
     lambda: get_acars_profile("2024", "01", "15", "14", "DEN_1420", hush=True), {}),
)


def make_winter_examples(out: Path, strict: bool) -> None:
    banner("WINTER SOUNDING EXAMPLES")

    for filename, label, get_data, settings in WINTER_EXAMPLES:
        run_example(
            f"Winter sounding: {label}",
            lambda get_data=get_data, filename=filename, settings=settings: save_sounding(
                get_data(),
                out / "winter" / filename,
                mode="winter",
                show_obs=True,
                radar="single-site",
                **settings,
            ),
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
                get_obs_data("OAX", year, month, day, hour, hush=True)
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
# PROFILE TOOLS & SURFACE MODIFICATION
# =============================================================================

# the profile tools tutorial's forecast example: the 12Z HRRR run for Mobridge, SD on the day of the Mound City
# tornadoes (see the Mound City case study), through 00Z the next day
PARAMS_TREND_CASE = ("hrrr", "KMBG", 0, "2024", "08", "28", "12")
PARAMS_TREND_HOURS = 13


def save_params_trend(output: Path) -> None:
    """MLCAPE and 0-1 km SRH at each forecast hour, as two panels sharing the forecast-hour axis."""
    forecast = get_bufkit_data(*PARAMS_TREND_CASE, hush=True, all_fcst_hours=True)[:PARAMS_TREND_HOURS]
    hours, mlcape, srh = [], [], []
    for profile in forecast:
        params = spy.sounding_params(profile).calc()
        hours.append(int(profile["site_info"]["fcst-hour"].lstrip("F")))
        mlcape.append(params["thermo"]["mlcape"])
        srh.append(params["kinem"]["srh_0_to_1000"])

    ink, muted, grid, line = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
    fig, axes = plt.subplots(2, 1, figsize=(8, 5.4), sharex=True, dpi=DPI, facecolor="white")
    for ax, values, title, unit in ((axes[0], mlcape, "MLCAPE", "J/kg"),
                                    (axes[1], srh, "0-1 km storm-relative helicity", "m²/s²")):
        ax.plot(hours, values, color=line, linewidth=2, marker="o", markersize=6,
                markeredgecolor="white", markeredgewidth=1.5, zorder=3)
        ax.set_title(f"{title} ({unit})", loc="left", fontsize=11, color=ink, fontweight="bold")
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color=grid, linewidth=0.8)
        ax.tick_params(colors=muted, labelsize=9, length=0)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color(muted)
    axes[1].set_xticks(hours[::2])
    axes[1].set_xticklabels([f"F{h:02d}\n{(12 + h) % 24:02d}Z" for h in hours[::2]])
    axes[1].set_xlabel("forecast hour and valid time (UTC)", fontsize=9, color=muted)
    fig.suptitle("HRRR at Mobridge, SD (KMBG) · run 12Z 28 Aug 2024", x=0.06, ha="left",
                 fontsize=12, color=ink)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")


def make_profile_tool_examples(data: dict, out: Path, strict: bool) -> None:
    banner("PROFILE TOOLS & SURFACE MODIFICATION EXAMPLES")

    # 1. a zoomed skew-t
    run_example(
        "Zoomed skew-t (ylim)",
        lambda: save_sounding(data, out / "profile_tools" / "zoom_ylim.png", ylim=(1000, 500), radar=None, map_zoom=0),
        strict=strict,
    )

    # 2-3. smoothing, with each filter
    import copy

    run_example(
        "Savitzky-Golay smoothing",
        lambda: save_sounding(spy.smooth_profile(copy.deepcopy(data), method="savgol", window_length=15),
                              out / "profile_tools" / "smooth_savgol.png", radar=None, map_zoom=0),
        strict=strict,
    )
    run_example(
        "Gaussian smoothing",
        lambda: save_sounding(spy.smooth_profile(copy.deepcopy(data), sigma=3),
                              out / "profile_tools" / "smooth_gaussian.png", radar=None, map_zoom=0),
        strict=strict,
    )

    # 4. two soundings merged, 60/40
    def merged_example() -> None:
        later = get_obs_data("OAX", "2014", "06", "17", "00", hush=True)
        merged = spy.merge_profiles(data, later, weight_a=0.6)
        merged["site_info"] = dict(data["site_info"], source="MERGED RAOB PROFILES")
        merged["titles"] = {"top_title": "MERGED VERTICAL PROFILE",
                            "left_title": "60% 18Z OAX 06/16/2014 | 40% 00Z OAX 06/17/2014",
                            "right_title": data["titles"]["right_title"]}
        save_sounding(merged, out / "profile_tools" / "merged.png", radar=None, map_zoom=0)

    run_example("Merged profiles", merged_example, strict=strict)

    # 5. fewer levels
    run_example(
        "Thinned profile (max_levels)",
        lambda: save_sounding(data, out / "profile_tools" / "max_levels.png", max_levels=30, radar=None, map_zoom=0),
        strict=strict,
    )

    # parameters through a forecast: one BUFKIT download with every forecast hour, sounding_params at each
    run_example("Parameters through a forecast", lambda: save_params_trend(out / "profile_tools" / "params_trend.png"),
                strict=strict)

    # 6-9. each surface modification method on the 12Z sounding's surface inversion,
    # warmed to the 18Z surface observation and zoomed to the lower atmosphere
    morning = get_obs_data("OAX", "2014", "06", "16", "12", hush=True)
    for name, modify_sfc in (("none", None),
                             ("mix", {"T": 28, "Td": 24}),
                             ("blend", {"T": 28, "Td": 24, "method": "blend"}),
                             ("replace", {"T": 28, "Td": 24, "method": "replace"})):
        run_example(
            f"Surface modification: {name}",
            lambda name=name, modify_sfc=modify_sfc: save_sounding(
                morning, out / "surface_modification" / f"sfc_{name}.png",
                modify_sfc=modify_sfc, ylim=(1000, 600), radar=None, map_zoom=0),
            strict=strict,
        )


# =============================================================================
# TUTORIAL — CUSTOM MAP DATA
# =============================================================================

# the custom map data tutorial uses the Aberdeen, SD RAOB at 00 UTC on 29 August 2024 (see the Mound City case study)
CUSTOM_MAP_CASE = {"station": "ABR", "year": "2024", "month": "08", "day": "29", "hour": "00"}


def make_custom_map_examples(out: Path, strict: bool) -> None:
    banner("CUSTOM MAP DATA EXAMPLES")

    import tempfile

    import cartopy.crs as ccrs
    import fsspec
    import numpy as np
    import pandas as pd
    import xarray as xr
    from metpy.calc import reduce_point_density

    case = CUSTOM_MAP_CASE
    data = get_obs_data(case["station"], case["year"], case["month"], case["day"], case["hour"], hush=True)
    # each example keeps its layers here so the last one can use them all
    layers = {}

    # 1. METAR dewpoints from the Iowa Environmental Mesonet (same code as the tutorial)
    def metar_example() -> None:
        url = ("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?"
               "network=SD_ASOS&network=ND_ASOS&network=NE_ASOS&network=MN_ASOS"
               "&data=dwpf&sts=2024-08-28T23:30Z&ets=2024-08-29T00:30Z&tz=Etc/UTC"
               "&format=onlycomma&latlon=yes&missing=empty&report_type=3")
        metars = pd.read_csv(url, parse_dates=["valid"]).dropna(subset=["dwpf"])
        metars["offset"] = (metars["valid"] - pd.Timestamp("2024-08-29 00:00")).abs()
        metars = metars.sort_values("offset").drop_duplicates("station")

        def plot_metars(ax):
            # keep the stations inside the map (with room for their labels), then thin them so labels don't overlap
            west, east, south, north = ax.get_extent(ccrs.PlateCarree())
            obs = metars[metars["lon"].between(west + 0.2, east - 0.2) & metars["lat"].between(south + 0.3, north - 0.3)]
            obs = obs[reduce_point_density(np.column_stack([obs["lon"], obs["lat"]]), 0.35)]
            ax.scatter(obs["lon"], obs["lat"], c=obs["dwpf"], cmap="BrBG", vmin=40, vmax=75, s=60,
                       edgecolor="black", transform=ccrs.PlateCarree(), zorder=8)
            for _, ob in obs.iterrows():
                ax.text(ob["lon"], ob["lat"] + 0.12, f"{ob['dwpf']:.0f}", ha="center", va="bottom", fontsize=8,
                        weight="bold", transform=ccrs.PlateCarree(), zorder=8)

        layers["metars"] = plot_metars
        save_sounding(data, out / "custom_map" / "metar_dewpoints.png", radar=None, map_zoom=3,
                      map_layers=plot_metars, map_label="METAR DEWPOINTS (F), 00 UTC")

    run_example("METAR dewpoints on the map", metar_example, strict=strict)

    # 2. RAP analysis SBCAPE & MSLP from NOAA's RAP archive on AWS (same code as the tutorial)
    def rap_example() -> None:
        with tempfile.TemporaryDirectory() as tmp:
            grib = f"{tmp}/rap.grib2"
            fs = fsspec.filesystem("s3", anon=True)
            fs.get("s3://noaa-rap-pds/rap.20240829/rap.t00z.awp130pgrbf00.grib2", grib)

            def open_rap(keys):
                return xr.open_dataset(grib, engine="cfgrib", backend_kwargs={"filter_by_keys": keys, "indexpath": ""})

            cape = open_rap({"typeOfLevel": "surface", "stepType": "instant", "shortName": "cape"})["cape"].load()
            mslp = (open_rap({"typeOfLevel": "meanSea"})["mslma"] / 100).load()
            lat, lon = cape.latitude, cape.longitude

            def plot_cape(ax):
                ax.contourf(lon, lat, cape, levels=np.arange(250, 5001, 250), cmap="magma_r", alpha=0.6,
                            transform=ccrs.PlateCarree(), zorder=3)

            def plot_mslp(ax):
                ax.contour(lon, lat, mslp, levels=np.arange(980, 1041, 2), colors="black", linewidths=1,
                           transform=ccrs.PlateCarree(), zorder=7)

            layers["cape"], layers["mslp"] = plot_cape, plot_mslp
            save_sounding(data, out / "custom_map" / "rap_cape_mslp.png", radar=None, map_zoom=3,
                          map_layers=[plot_cape, plot_mslp], map_label="RAP SBCAPE (J/kg) & MSLP (2 hPa)")

    run_example("RAP SBCAPE & MSLP on the map", rap_example, strict=strict)

    # 3. ...and everything at once (the tutorial's closing joke): single-site radar, CAPE, isobars & METARs
    def everything_example() -> None:
        save_sounding(data, out / "custom_map" / "everything.png", radar="single", map_zoom=3,
                      map_layers=[layers["cape"], layers["mslp"], layers["metars"]],
                      map_label="RAP SBCAPE & MSLP + METAR Td + RADAR")

    if {"metars", "cape", "mslp"} <= set(layers):
        run_example("Everything on the map", everything_example, strict=strict)


# =============================================================================
# PUBLICATION & CUSTOM FIGURES
# =============================================================================

# HRRR and RAP 3-hour forecasts for Aberdeen, SD, valid at the custom map case's time (00Z 29 August 2024)
CUSTOM_FIGURES_RUN = ("kabr", 3, "2024", "08", "28", "21")


def make_custom_figure_examples(out: Path, strict: bool) -> None:
    banner("PUBLICATION & CUSTOM FIGURE EXAMPLES")

    import cartopy.crs as ccrs

    folder = out / "custom_figures"
    folder.mkdir(parents=True, exist_ok=True)
    hrrr = get_bufkit_data("hrrr", *CUSTOM_FIGURES_RUN, hush=True)
    rap = get_bufkit_data("rap", *CUSTOM_FIGURES_RUN, hush=True)
    both = dict(profiles=[hrrr, rap], labels=["HRRR F03", "RAP F03"])
    # the observed Aberdeen sounding and the HRRR forecast for Ashley, ND: two places on the map
    raob = get_obs_data(CUSTOM_MAP_CASE["station"], CUSTOM_MAP_CASE["year"], CUSTOM_MAP_CASE["month"],
                            CUSTOM_MAP_CASE["day"], CUSTOM_MAP_CASE["hour"], hush=True)
    ashley = get_bufkit_data("hrrr", "kasy", *CUSTOM_FIGURES_RUN[1:], hush=True)

    def publication(name: str, **kwargs) -> Callable[[], None]:
        def example() -> None:
            kwargs.setdefault("profiles", hrrr)
            fig, _ = spy.build_publication_sounding(kwargs.pop("profiles"), dpi=DPI, **kwargs)
            fig.savefig(folder / name)
            close_figures()
        return example

    # the same calls as the tutorial (tutorials/custom_figures.rst)
    examples = {
        "publication_quickstart.png": dict(profiles=[raob, ashley], labels=["00Z ABR RAOB", "HRRR F03 Ashley, ND"],
                                           panels=("skewt", "hodograph", "map", "parameters"),
                                           parameters=("mllcl", "mlcape", "mlcin", "mucape", "srh01", "srh03", "bwd06", "stp"),
                                           map_radar="single", map_zoom=1),
        "publication_default.png": dict(),
        "publication_pair.png": dict(panels=("skewt", "hodograph")),
        "publication_zoom_parameters.png": dict(**both, panels=("skewt", "parameters"), ylim=(1000, 500), parcel="sb",
                                                parameters=("sbcape", "mlcin", ("0–3 km SBCAPE", "thermo", "sb3cape", 1, r"J kg$^{-1}$", 0), "mllcl")),
        "publication_hodograph.png": dict(wind_units="kt", hodo_heights=(0, 500, 1000, 3000, 6000, 9000), hodo_max_height=9000,
                                          sr_vector_heights=(500, 1000, 3000), parameters=("srh01", "srh03", "bwd06", "scp")),
        "publication_style.png": dict(**both, colors=["black", "darkorange"], show_adiabats=True, parcel="mu", sr_vector_heights=(1000,),
                                      style={"font_family": "serif", "fontsize": 9, "linewidth": 1.4, "hodograph_rings": False,
                                             "sr_vector_color": "purple"}),
        "publication_radar.png": dict(**both, panels=("skewt", "hodograph", "map"), map_radar="single", map_zoom=1),
    }
    for name, kwargs in examples.items():
        run_example(name, publication(name, **kwargs), strict=strict)

    def annotated_example() -> None:
        fig, axes = spy.build_publication_sounding([hrrr, rap], labels=["HRRR F03", "RAP F03"],
                                                   panels=("skewt", "hodograph", "map"), map_zoom=2, dpi=DPI)
        axes["skewt"].ax.axhspan(850, 700, color="gold", alpha=0.25, zorder=0)
        axes["skewt"].ax.text(-35, 775, "EML", fontsize=8, va="center")
        axes["hodograph"].annotate("0–3 km SRH\n342 m$^2$ s$^{-2}$", xy=(18, -6), fontsize=6, ha="center", va="top", style="italic")
        axes["map"].plot(-100.07, 45.73, marker="*", color="black", markersize=7, transform=ccrs.PlateCarree(), zorder=11)
        axes["map"].text(-99.9, 45.73, "Mound City", fontsize=6, va="center", transform=ccrs.PlateCarree(), zorder=11)
        fig.savefig(folder / "publication_annotated.png")
        close_figures()

    def layouts_example() -> None:
        # every set of panels, drawn small and side by side, each titled with the panels that make it
        from io import BytesIO
        from sounderpy.plot.publication import PUBLICATION_LAYOUTS

        images = []
        for panels in PUBLICATION_LAYOUTS:
            fig, _ = spy.build_publication_sounding([hrrr, rap], labels=["HRRR F03", "RAP F03"], panels=panels,
                                                    map_radar="single", map_zoom=1, dpi=120)
            image = BytesIO()
            fig.savefig(image, format="png")
            plt.close(fig)
            image.seek(0)
            images.append((panels, plt.imread(image)))

        # every figure is 6.5 in wide, so rows as tall as their tallest figure draw them all at the same scale
        heights = [max(images[i][1].shape[0], images[i + 1][1].shape[0]) for i in range(0, len(images), 2)]
        sheet, cells = plt.subplots(len(heights), 2, figsize=(10, 5 * sum(heights) / images[0][1].shape[1] + 2), dpi=DPI,
                                    gridspec_kw={"height_ratios": heights, "hspace": 0.12, "wspace": 0.04})
        for cell, (panels, image) in zip(cells.flat, images):
            cell.imshow(image, extent=(0, image.shape[1], image.shape[0], 0))
            cell.set_xlim(0, image.shape[1])
            cell.set_ylim(heights[list(cells.flat).index(cell) // 2], 0)   # anchored to the top of its row
            cell.set_title(f"panels={panels}", fontsize=9, family="monospace")
            cell.set_axis_off()
        sheet.savefig(folder / "publication_layouts.png", bbox_inches="tight", facecolor="white")
        close_figures()

    def return_fig_example() -> None:
        fig, axes = spy.build_sounding(hrrr, radar="single", map_zoom=1, return_fig=True, dpi=100)
        axes["skewt"].ax.axhspan(850, 700, color="gold", alpha=0.2, zorder=0)
        axes["skewt"].ax.text(-28, 775, "EML", fontsize=16, weight="bold", va="center")
        fig.savefig(folder / "return_fig_sounding.png", bbox_inches="tight")
        close_figures()

    def operational_example() -> None:
        spy.build_sounding(hrrr, radar="single", map_zoom=1, save=True, filename=str(folder / "operational_sounding.png"), dpi=100)
        close_figures()

    run_example("operational_sounding.png", operational_example, strict=strict)
    run_example("publication_annotated.png", annotated_example, strict=strict)
    run_example("publication_layouts.png", layouts_example, strict=strict)
    run_example("return_fig_sounding.png", return_fig_example, strict=strict)


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
    "winter",
    "composites",
    "custom-map",
    "custom-figures",
    "profile-tools",
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

    if "winter" in args.sections:
        make_winter_examples(out, args.strict)

    if "composites" in args.sections:
        make_composite_examples(out, args.strict)

    if "custom-map" in args.sections:
        make_custom_map_examples(out, args.strict)

    if "profile-tools" in args.sections:
        make_profile_tool_examples(get_base_obs(), out, args.strict)

    if "custom-figures" in args.sections:
        make_custom_figure_examples(out, args.strict)

    banner("GALLERY GENERATION FINISHED")
    print(f"Images written beneath: {out.resolve()}")

    print(summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

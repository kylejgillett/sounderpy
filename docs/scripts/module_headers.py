#!/usr/bin/env python3
"""
Write the header docstring at the top of every module in src/sounderpy: the module's name, its purpose,
links to its docs, its references (from sounderpy_references.py), and the copyright line.

Each module's old header is replaced wherever it sits, and only the header changes. Run from any
directory (paths are found from this script's own location):

    python docs/scripts/module_headers.py              # rewrite the headers in place
    python docs/scripts/module_headers.py --check      # list modules whose header is out of date
    python docs/scripts/module_headers.py --output DIR # write the rewritten modules to DIR instead
"""

import argparse
import ast
import os
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sounderpy_references import REFS, TAGS  # noqa: E402

PACKAGE = Path(__file__).resolve().parents[2] / "src" / "sounderpy"

D = "https://sounderpy.app/docs/"

SHORT_HEADERS = {
    "__init__.py": ("SOUNDERPY: vertical profile data retrieval and analysis tool for Python", "2026"),
    "__main__.py": ("SOUNDERPY: runs the SounderPy command line tool with `python -m sounderpy`", "2026"),
}

SPECS = {
"calc.py": dict(
    title="SOUNDERPY CALC",
    purpose="""Compute sounding parameters from a SounderPy `clean_data` profile: surface-based, mixed-layer,
and most-unstable parcels (CAPE, CIN, LCL, LFC, EL), entraining CAPE, DCAPE and DCIN, lapse rates,
the PBL, bulk shear, storm-relative helicity and winds, streamwise vorticity, storm motions,
the effective inflow layer, and composite indices. `sounding_params(clean_data).calc()` returns
them in a dict; `vad_params(vad_data).calc()` does the same for VAD wind profiles.""",
    docs=[("sounding parameters", D + "reference/helpertools.html#returning-a-dictionary-of-profile-parameters"),
          ("storm motions", D + "reference/plotsettings.html#storm-motions"),
          ("parcels", D + "reference/plotsettings.html#parcels"),
          ("surface modification", D + "reference/plotsettings.html#surface-modification")],
    refs=["sharppy", "metpy", "doswell_rasmussen1994", "blanchard1998", "peters2023", "emanuel1994", "market_dcape",
          "gilmore_wicker1998", "bunkers2000", "bunkers2014", "corfidi2003", "nixon_allen2021", "davies_jones1984",
          "davies_jones1990", "brandt2022", "thompson2004", "thompson2007", "thompson2012", "hart_korotky1991",
          "davies1993", "weisman_klemp1982", "george1960", "miller1972", "spc_ship", "johnson_sugden2014",
          "cohen2007", "browning_wexler1968"],
    years="2024, 2025, 2026"),
"calc_winter.py": dict(
    title="SOUNDERPY WINTER CALC",
    purpose="""Compute winter-weather sounding parameters: precipitation type, snow growth, and winter
hazards. `winter_params(clean_data).calc()` runs the standard `sounding_params` calc and adds
a `winter` dict to what it returns.

Sources are noted beside each parameter: [peer-reviewed] journal articles, [operational] SPC/NWS
practice that isn't formally published, [informal] methods without a paper, and [SounderPy]
SounderPy's own choices, not from a paper.""",
    docs=[("winter soundings", D + "tutorials/winter_soundings.html"),
          ("winter parameters", D + "reference/helpertools.html#returning-a-dictionary-of-winter-parameters"),
          ("winter values in code", D + "tutorials/winter_soundings.html#winter-values-in-code"),
          ("winter references", D + "tutorials/winter_soundings.html#references")],
    refs=[("sharppy", "SPC precipitation type, OPRH, and PBL"), "metpy", "alduchov1996", "banacos2014", "bernstein2005", "birk2021", "bourgouin2000",
          "gultepe2007", "huffman1988", "libbrecht2005", "mcdonough2004", "kuchera", "nws_windchill",
          "osczevski2005", "rauber2000", "takahashi1991"],
    years="2026"),
"cli.py": dict(
    title="SOUNDERPY CLI",
    purpose="""The `sounderpy` command line tool: get profiles, export them, print their parameters, and plot
soundings and hodographs from a shell. It is a thin wrapper around SounderPy's Python functions,
so the calculations, data retrieval, plotting, and file export all live in the package modules.""",
    docs=[("command line tool", D + "reference/cli.html"),
          ("command line workflows", D + "tutorials/cli_workflows.html")],
    refs=[],
    years="2026"),
"profile_utils.py": dict(
    title="SOUNDERPY PROFILE UTILITIES",
    purpose="""Merge and smooth profiles: `merge_profiles()` blends two `clean_data` profiles by weighted
averaging on a common height axis, and `smooth_profile()` smooths a profile with a Savitzky-Golay
or Gaussian filter.""",
    docs=[("merging profiles", D + "reference/helpertools.html#merging-two-profiles-via-weighted-averaging"),
          ("smoothing a profile", D + "reference/helpertools.html#smooth-a-profile")],
    refs=["savitzky_golay1964", ("scipy", "Gaussian and Savitzky-Golay filters")],
    years="2024, 2025, 2026"),
"sounderpy.py": dict(
    title="SOUNDERPY | Vertical Profile Data Retrieval and Analysis Tool For Python",
    purpose="""SounderPy is an open-source atmospheric science Python package for vertical profile analysis.
This tool is designed to get data, 'clean it up' for simple use, and plot the data on advanced-sounding
plots. SounderPy was developed with the goal in mind to keep the code simple and efficient for users of
all experience levels and for reliability in all use cases.

This module holds the public functions: get_obs_data(), get_bufkit_data(), get_model_data(),
build_sounding(), build_hodograph(), build_composite(), build_vad_hodograph(), and print_variables().

Version: 3.3.0 | November 2026
Code: https://github.com/kylejgillett/sounderpy
PyPI: https://pypi.org/project/sounderpy/
Operational site: https://sounderpy.app/""",
    docs=[("documentation", D),
          ("getting data", D + "reference/gettingdata.html"),
          ("plotting data", D + "reference/plottingdata.html"),
          ("helper tools", D + "reference/helpertools.html")],
    refs=["sounderpy_joss"],
    years="2023, 2024, 2025, 2026"),
"utils.py": dict(
    title="SOUNDERPY UTILITIES",
    purpose="""Helper functions for working with vertical profile data: surface modification
(spy.modify_surface()), interpolation, nearest values, surface-based profiles, station lookups for
the RAOB and BUFKIT station lists, readout formatting, and SounderPy's warning classes. Some are
called by other modules, and some are meant to be called by the user.""",
    docs=[("helper tools", D + "reference/helpertools.html"),
          ("surface modification", D + "reference/plotsettings.html#surface-modification"),
          ("station identifiers", D + "reference/gettingdata.html#station-identifiers")],
    refs=["holzworth1964", "betts1973", "carson1973", "tennekes1973", "metpy"],
    years="2024, 2026"),
"data/__init__.py": dict(
    title="SOUNDERPY DATA",
    purpose="""The modules that get, decode, read, and write profile data:

acars_data.py            : ACARS aircraft profiles from OU
bufkit_data.py           : BUFKIT model forecast profiles from PSU and IEM
ncep_bufr.py             : RRFS, HRRR, and RAP BUFR profiles from NCEP NOMADS
obs_data.py              : RAOB and IGRA observed profiles
model_data.py            : ERA5, NCEP-FNL, and RAP analysis profiles from THREDDS and CDS
rapruc_grib_retriever.py : RAP/RUC GRIB downloads from AWS and NCEI
rapruc_grib_parser.py    : RAP/RUC GRIB parsing
radar_utils.py           : radar mosaic and single-site radar data for the map inset
sfc_obs.py               : the surface observation nearest a sounding, from the IEM ASOS archive
wrf_utils.py             : WRF output profiles
cm1_utils.py             : CM1 input_sounding profiles
file_io.py               : reading and writing SHARPpy, CSV, and CM1 files""",
    docs=[("getting data", D + "reference/gettingdata.html"),
          ("custom data sources", D + "reference/customdatasources.html"),
          ("the clean_data dictionary", D + "reference/clean_data_scheme.html")],
    refs=[],
    years="2026"),
"data/acars_data.py": dict(
    title="SOUNDERPY ACARS DATA",
    purpose="""List and load ACARS aircraft profiles from the University of Oklahoma's archive and parse them into
SounderPy `clean_data` dictionaries. `acars_data(year, month, day, hour)` lists the profiles for an
hour, and `.get_profile(profile_id)` returns one.""",
    docs=[("ACARS data", D + "reference/gettingdata.html#observed-data-acars-aircraft-profiles"),
          ("ACARS tutorial", D + "tutorials/data_sources.html#acars-aircraft-profiles")],
    refs=["ou_acars", "moninger2003", "bs4"],
    years="2024, 2026"),
"data/bufkit_data.py": dict(
    title="SOUNDERPY BUFKIT DATA",
    purpose="""Load and parse BUFKIT model forecast profiles (GFS, NAM, NAM Nest, RAP, HRRR, and RRFS) into
SounderPy `clean_data` dictionaries: the latest runs from Penn State, archived runs from the
Iowa Environmental Mesonet, and NCEP's own BUFR soundings from NOMADS (decoded by ncep_bufr.py)
when the RRFS, HRRR, or RAP file is missing or behind. Called by spy.get_bufkit_data().""",
    docs=[("BUFKIT data", D + "reference/gettingdata.html#model-forecast-data-bufkit-profiles"),
          ("BUFKIT sites", D + "reference/gettingdata.html#available-bufkit-sites"),
          ("BUFKIT tutorial", D + "tutorials/data_sources.html#bufkit-forecast-profiles")],
    refs=["psu_bufkit", "iem_bufkit", "nomads", "benjamin2016", "dowell2022"],
    years="2024, 2026"),
"data/cm1_utils.py": dict(
    title="SOUNDERPY CM1 UTILITIES",
    purpose="""Read CM1 input_sounding files into SounderPy `clean_data` dictionaries with `from_cm1()`, so CM1
initial soundings can be analyzed and plotted like any other profile.""",
    docs=[("CM1 input soundings", D + "reference/customdatasources.html#ingesting-cm1-input-sounding-data")],
    refs=["bryan_fritsch2002", "metpy"],
    years="2024, 2025, 2026"),
"data/file_io.py": dict(
    title="SOUNDERPY FILE READING AND WRITING",
    purpose="""Write profiles to CSV, SHARPpy, and CM1 input_sounding files (`to_file`), and read SHARPpy files and
GRAW, iMet, and Windsond radiosonde files into SounderPy `clean_data` dictionaries (`from_sharppy`,
`from_graw`, `from_imet`, `from_windsond`).""",
    docs=[("saving data to a file", D + "reference/helpertools.html#saving-data-to-a-file"),
          ("exporting data", D + "tutorials/exporting_data.html"),
          ("reading SHARPpy files", D + "reference/customdatasources.html#ingesting-sharppy-files"),
          ("reading radiosonde files", D + "reference/customdatasources.html#ingesting-graw-imet-and-windsond-radiosonde-files")],
    refs=[("sharppy", "SHARPpy file format"), "bryan_fritsch2002", "lofs_sndmod", "metpy"],
    years="2026"),
"data/model_data.py": dict(
    title="SOUNDERPY MODEL REANALYSIS DATA",
    purpose="""Load and parse ERA5, NCEP-FNL, and latest RAP analysis profiles into SounderPy `clean_data`
dictionaries: ERA5 from the Copernicus Climate Data Store, and NCEP-FNL and the latest RAP analysis ('rap-now')
through THREDDS. RAP/RUC reanalysis ('rap-ruc') comes from rapruc_grib_retriever.
Called by spy.get_model_data().""",
    docs=[("reanalysis data", D + "reference/gettingdata.html#model-reanalysis-data-rap-ruc-era5-ncep"),
          ("model keys", D + "reference/gettingdata.html#model-key-names"),
          ("RAP/RUC tutorial", D + "tutorials/data_sources.html#rap-ruc-reanalysis"),
          ("other reanalysis sources", D + "tutorials/data_sources.html#other-model-reanalysis-sources")],
    refs=["hersbach2020", "benjamin2016", "ncep_fnl", ("siphon", "THREDDS NetCDF subsets"), "xarray", "cdsapi"],
    years="2024, 2026"),
"data/ncep_bufr.py": dict(
    title="SOUNDERPY NCEP BUFR SOUNDING DECODER",
    purpose="""Decode NCEP model BUFR sounding files (the 'bufrsnd' station files NCEP builds BUFKIT files from)
and write them out as BUFKIT text, so that bufkit_data.py can read them like any other BUFKIT file.
NCEP BUFR files carry their own descriptor tables (DX tables) in the first few messages, so no outside
BUFR tables or libraries are needed. Also finds and downloads the newest RRFS, HRRR, and RAP soundings
on NOMADS.""",
    docs=[("BUFKIT data (NOMADS backup)", D + "reference/gettingdata.html#model-forecast-data-bufkit-profiles")],
    refs=["wmo_bufr", "nomads", "bolton1980"],
    years="2026"),
"data/obs_data.py": dict(
    title="SOUNDERPY RAOB OBSERVATIONS DATA",
    purpose="""Load and parse observed radiosonde (RAOB) profiles into SounderPy `clean_data` dictionaries, from
the University of Wyoming archive and the Integrated Global Radiosonde Archive (IGRA v2), using
SounderPy's packaged station list to match WMO, ICAO, and IGRA IDs. Called by spy.get_obs_data().""",
    docs=[("RAOB data", D + "reference/gettingdata.html#observed-data-raob-igrav2-profiles"),
          ("RAOB sites", D + "reference/gettingdata.html#available-raob-sites"),
          ("RAOB tutorial", D + "tutorials/data_sources.html#observed-raob-igra-profiles")],
    refs=["uwyo", "durre2018", ("siphon", "Wyoming and IGRA readers")],
    years="2024, 2026"),
"data/radar_utils.py": dict(
    title="SOUNDERPY RADAR UTILITIES",
    purpose="""Radar data for the map inset: the nearest NEXRAD site's lowest-tilt reflectivity from the Unidata
NEXRAD Level II archive on AWS, the CONUS reflectivity mosaic from Unidata's THREDDS server, the
radar valid time, and the reflectivity colormap.""",
    docs=[("map insets & radar", D + "reference/plotsettings.html#map-radar"),
          ("maps & radar tutorial", D + "tutorials/maps_radar.html")],
    refs=["pyart", "nexrad_aws", ("siphon", "THREDDS radar mosaic")],
    years="2025, 2026"),
"data/rapruc_grib_parser.py": dict(
    title="SOUNDERPY RAP/RUC GRIB PARSER",
    purpose="""Parse RAP and RUC GRIB files (pressure-level and native-level) into a SounderPy `clean_data`
profile: open the GRIB files with cfgrib and xarray, average the grid points in a box around the
requested point, and build the profile. Parsing is thread safe, so profiles can be fetched in parallel.

* the development of this module was aided by artificial intelligence""",
    docs=[("reanalysis data", D + "reference/gettingdata.html#model-reanalysis-data-rap-ruc-era5-ncep"),
          ("RAP/RUC tutorial", D + "tutorials/data_sources.html#rap-ruc-reanalysis")],
    refs=["benjamin2016", "cfgrib", "xarray", "metpy"],
    years="2026"),
"data/rapruc_grib_retriever.py": dict(
    title="SOUNDERPY RAP-RUC DATA RETRIEVAL - from NCEI and AWS",
    purpose="""Downloads RAP/RUC GRIB files from AWS and NCEI's archive for the requested date
(NCEI no longer serves the NetCDF subsets earlier versions used).

Part of this script is based on the AWS method used by Chase Archive.""",
    docs=[("reanalysis data", D + "reference/gettingdata.html#model-reanalysis-data-rap-ruc-era5-ncep"),
          ("RAP/RUC tutorial", D + "tutorials/data_sources.html#rap-ruc-reanalysis")],
    refs=["benjamin2016", "noaa_rap_aws", "ncei_rap", "fsspec"],
    years="2026"),
"data/sfc_obs.py": dict(
    title="SOUNDERPY SURFACE OBSERVATION",
    purpose="""Finds the ASOS/AWOS stations nearest a sounding and fetches the METAR closest to its valid time
from the Iowa Environmental Mesonet archive (no API key). The winter sounding draws it as a
station model in the hodograph's upper-right corner (build_sounding(..., mode='winter', show_obs=True)).

    obs = sounding_obs(clean_data)          # None if nothing within reach
    obs["station"], obs["time"], obs["tmpf"], obs["wxcodes"], obs["timing"] ...

Station coordinates come from MetPy's built-in station table (downloaded once to ~/.cache/metpy,
then read locally). IEM rate-limits rapid requests, so callers should treat a failure as "no obs".""",
    docs=[("surface observation", D + "tutorials/winter_soundings.html#add-the-latest-surface-observation")],
    refs=["iem_asos", ("metpy", "station table")],
    years="2026"),
"data/wrf_utils.py": dict(
    title="SOUNDERPY WRF UTILITIES",
    purpose="""Read WRF-ARW output (a netCDF4 dataset of a wrfout file) into a SounderPy `clean_data` dictionary
for a latitude and longitude with `from_wrf()`, so WRF profiles can be analyzed and plotted like any
other profile.""",
    docs=[("WRF output", D + "reference/customdatasources.html#ingesting-wrf-output-data")],
    refs=["skamarock2019", ("bolton1980", "dewpoint from vapor pressure"), "metpy"],
    years="2024, 2026"),
"plot/__init__.py": dict(
    title="SOUNDERPY PLOTTING",
    purpose="""The modules that draw SounderPy's figures:

style.py             : themes, colors, units, layouts, and watermark settings
components.py        : reusable pieces that figures are built from
components_winter.py : the winter sounding's own pieces
operational.py       : full (severe and winter) sounding, hodograph, composite, and VAD figures
publication.py       : journal-style publication figures""",
    docs=[("plotting data", D + "reference/plottingdata.html")],
    refs=[],
    years="2026"),
"plot/components.py": dict(
    title="SOUNDERPY PLOT COMPONENTS",
    purpose="""The reusable pieces that SounderPy figures are built from: the skew-T and its annotations, the
hodograph, side panels, parameter tables, map insets with radar and custom map layers, and figure
text. Figures are put together in operational.py.""",
    docs=[("soundings", D + "reference/plottingdata.html#building-soundings"),
          ("hodographs", D + "reference/plottingdata.html#building-hodographs"),
          ("map insets & radar", D + "reference/plotsettings.html#map-radar"),
          ("custom map data", D + "tutorials/custom_map_data.html")],
    refs=["matplotlib", ("metpy", "skew-T and hodograph"), "cartopy", ("peters2023", "ECAPE parcel traces"), "cohen2007"],
    years="2024, 2026"),
"plot/components_winter.py": dict(
    title="SOUNDERPY WINTER PLOT COMPONENTS",
    purpose="""Drawing tools for the winter sounding: SounderPy's full sounding layout with winter content in each slot.
The skew-t shows melting and refreezing layers, precipitation phase changes, 0 °C levels, clouds,
icing and the DGZ; the hodograph shows PBL and DGZ mean winds, a winds box and the latest surface
observation; the side panels show T & Tw, RHᵢ & RH, wind and theta-e (0-3 km); the tables show precipitation
type, snow growth and winter hazards. SounderPy's own plot components are reused wherever they fit.""",
    docs=[("winter soundings", D + "tutorials/winter_soundings.html"),
          ("a tour of the figure", D + "tutorials/winter_soundings.html#a-tour-of-the-figure")],
    refs=["matplotlib", ("metpy", "skew-T"), "gultepe2007", "huffman1988", "mcdonough2004", "rauber2000", "sczepanski2026"],
    years="2026"),
"plot/operational.py": dict(
    title="SOUNDERPY OPERATIONAL FIGURES",
    purpose="""Build SounderPy's figures from the pieces in components.py and components_winter.py: the full
severe and winter soundings, the hodograph, composite soundings, and the VAD hodograph,
including their map and radar insets. Called by spy.build_sounding(), spy.build_hodograph(),
spy.build_composite(), and spy.build_vad_hodograph().""",
    docs=[("soundings", D + "reference/plottingdata.html#building-soundings"),
          ("hodographs", D + "reference/plottingdata.html#building-hodographs"),
          ("composites", D + "reference/plottingdata.html#building-composite-soundings"),
          ("VAD hodographs", D + "reference/plottingdata.html#building-vad-hodographs"),
          ("winter soundings", D + "tutorials/winter_soundings.html"),
          ("map insets & radar", D + "reference/plotsettings.html#map-radar"),
          ("custom map data", D + "tutorials/custom_map_data.html")],
    refs=["matplotlib", ("metpy", "skew-T and hodograph"), "cartopy"],
    years="2024, 2026"),
"plot/publication.py": dict(
    title="SOUNDERPY PUBLICATION FIGURES",
    purpose="""Clean, journal-style figures of one or more profiles: a skew-t, and optionally a hodograph, a
location map, and a table of parameters, sized for AMS and AGU journals. The settings at the top
set every part of the look, and the figure and its panels are returned for annotating. Called by
spy.build_publication_sounding() (experimental in 3.3.0).""",
    docs=[("publication figures", D + "reference/publicationfigures.html"),
          ("tutorial", D + "tutorials/custom_figures.html"),
          ("map insets & radar", D + "reference/plotsettings.html#map-radar")],
    refs=["matplotlib", ("metpy", "skew-T"), "cartopy"],
    years="2026"),
"plot/style.py": dict(
    title="SOUNDERPY PLOT STYLE",
    purpose="""The look of SounderPy figures in one place: themes, font sizes, units, layouts, winter settings,
and watermark settings. Drawing tools are in components.py (and components_winter.py for the
winter sounding).""",
    docs=[("figure colors & themes", D + "reference/plotsettings.html#figure-colors"),
          ("theme settings", D + "reference/plotsettings.html#theme-settings"),
          ("watermark", D + "reference/plotsettings.html#watermark"),
          ("plot appearance", D + "tutorials/plot_appearance.html")],
    refs=[],
    years="2026"),
}


COPYRIGHT = "(C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, {years}"


def render(spec):
    out = ['"""', f"    {spec['title']}", "", "    Purpose of module:", ""]
    out += [("    " + line).rstrip() for line in spec["purpose"].split("\n")]
    if spec["docs"]:
        width = max(len(label) for label, _ in spec["docs"]) + 1
        out += ["", "    Docs:"]
        out += [f"        {(label + ':').ljust(width)}  {url}" for label, url in spec["docs"]]
    if spec["refs"]:
        out += ["", "    References:"]
        for item in spec["refs"]:
            key, tag = item if isinstance(item, tuple) else (item, TAGS[item])
            text = f"{REFS[key]} ({tag})"
            out += textwrap.wrap(text, width=108, initial_indent=" " * 8, subsequent_indent=" " * 12,
                                 break_long_words=False, break_on_hyphens=False)
    out += ["", "    " + COPYRIGHT.format(years=spec["years"]), '"""']
    return "\n".join(out)


def render_short(title, years):
    return "\n".join(['"""', f"    {title}", "", "    " + COPYRIGHT.format(years=years), '"""'])


def move(src, header):
    """Drop the module's old header string (wherever it sits) and put `header` at the top."""
    tree = ast.parse(src)
    lines = src.split("\n")
    old = next((n for n in tree.body if isinstance(n, ast.Expr) and isinstance(getattr(n, "value", None), ast.Constant)
                and isinstance(n.value.value, str)), None)
    if old is not None:
        start, end = old.lineno - 1, old.end_lineno
        # also drop the blank lines right after the old header
        while end < len(lines) and not lines[end].strip():
            end += 1
        # and blank lines right before it, leaving two
        while start > 0 and not lines[start - 1].strip():
            start -= 1
        keep_blank = 2 if start > 0 and end < len(lines) else 0
        lines = lines[:start] + [""] * keep_blank + lines[end:]
    while lines and not lines[0].strip():
        lines.pop(0)
    return header + "\n\n" + "\n".join(lines)


def rewrite(pkg, output=None, check=False):
    """Rewrite every module's header; return the modules whose header changed."""
    stale = []
    for rel in sorted(list(SPECS) + list(SHORT_HEADERS)):
        path = os.path.join(pkg, rel)
        raw = open(path, newline="", encoding="utf-8").read()
        crlf = "\r\n" in raw
        src = raw.replace("\r\n", "\n")
        header = render(SPECS[rel]) if rel in SPECS else render_short(*SHORT_HEADERS[rel])
        new = move(src, header)
        ast.parse(new)  # still valid Python
        doc = ast.get_docstring(ast.parse(new))
        assert doc and doc.split("\n")[0].strip().startswith("SOUNDERPY"), rel
        if crlf:
            new = new.replace("\n", "\r\n")
        if new != raw:
            stale.append(rel)
        if check:
            continue
        target = Path(output) / rel if output else Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        open(target, "w", newline="", encoding="utf-8").write(new)
    return stale


def main():
    parser = argparse.ArgumentParser(description="Write the header docstring at the top of every SounderPy module.")
    parser.add_argument("--check", action="store_true", help="only list modules whose header is out of date")
    parser.add_argument("--output", type=Path, help="write the rewritten modules here instead of in place")
    args = parser.parse_args()
    stale = rewrite(PACKAGE, args.output, check=args.check)
    if args.check:
        print("\n".join(stale) if stale else "every module header is up to date")
        return 1 if stale else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())

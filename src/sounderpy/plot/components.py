"""
    SOUNDERPY PLOT COMPONENTS

    Purpose of module:

    House the reusable pieces that SounderPy figures are built from: themes,
    layouts, skew-t and hodograph tools, accessory panels, text tables, and
    map insets. Figures are put together in operational.py.

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024, 2026
"""

# METPY OPERATIONS
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT, Hodograph, USCOUNTIES

# MATPLOTLIB OPERATIONS
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.collections import Collection
from matplotlib.lines import Line2D
from matplotlib.text import Text

# NUMPY OPERATIONS
import numpy as np
import numpy.ma as ma
import copy

# STANDARD PYTHON
import warnings
from functools import lru_cache

# ECAPE
from ecape_parcel.calc import calc_ecape_parcel, density_temperature

# FOR MAPPING
from PIL import Image
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import img_tiles
import shapely.geometry as sgeom

# STYLE
from .style import (FONT_SIZES,
                    COMPOSITE_HODOGRAPH_LAYERS,
                    COMPOSITE_PARAMETERS,
                    COMPOSITE_TABLES,
                    COMPOSITE_VIRIDIS_RANGES,
                    HODOGRAPH_KINEMATICS_LAYOUT,
                    HODOGRAPH_PARCEL_LAYOUT,
                    KINEMATICS_TABLE_LAYOUT,
                    MISSING_VALUE,
                    MOISTURE_TABLE_LAYOUT,
                    SKEWT_BOTTOM_PRESSURE,
                    THERMO_TABLE_LAYOUT,
                    UNIT_DEGC,
                    UNIT_GKG,
                    UNIT_JKG,
                    UNIT_KT,
                    UNIT_LAPSE,
                    UNIT_M,
                    UNIT_M2S2,
                    UNIT_PER_S,
                    WATERMARK,
                    WATERMARK_TEXT)

# SOUNDERPY MODULES
from ..utils import modify_surface, find_nearest, mag, mag_round, wind_to_dir
from ..calc import storm_motion_type
from ..data.radar_utils import *


#########################################################################
############################# CONSTANTS #################################
#########################################################################


PARCEL_DETAILS = {"term1": {"sb": "surface_based", "ml": "mixed_layer", "mu": "most_unstable"},
                  "term2": {"ps": [True, "pseudoadiabatic"], "ia": [False, "irrev-adiabatic"]},
                  "term3": {"ecape": True, "cape": False,},
                  "sm": {"rm": "right_moving", "lm": "left_moving", "mw": "mean_wind"},}

# the stepwise CAPE/CIN panel shows up to this height (m AGL), so parcels above it aren't needed
STABILITY_PANEL_DEPTH = 3000





















VALID_PARCEL_CODES = ("mu_ps_cape", "mu_ia_cape", "mu_ps_ecape", "mu_ia_ecape",
                      "ml_ps_cape", "ml_ia_cape", "ml_ps_ecape", "ml_ia_ecape",
                      "sb_ps_cape", "sb_ia_cape", "sb_ps_ecape", "sb_ia_ecape")




#########################################################################
######################### DATA-PREP UTILS ###############################
#########################################################################

def _thin_sounding(data, max_levels):
    # keep at most max_levels evenly spaced levels, always including the surface and top
    n_levels = len(data["p"])
    if max_levels is None or n_levels <= max_levels:
        return data

    keep = np.unique(np.round(np.linspace(0, n_levels - 1, max_levels)).astype(int))
    return {key: value[keep] if getattr(value, "shape", ())[:1] == (n_levels,) else value for key, value in data.items()}


def _prepare_sounding_data(clean_data, modify_sfc=None, max_levels=None):
    plot_data = _thin_sounding(copy.deepcopy(clean_data), max_levels)

    if modify_sfc is not None:
        plot_data = modify_surface(plot_data, modify_sfc)

    return plot_data



###### PARCEL LOGIC UTILS ######
################################
def _validate_special_parcels(special_parcels):
    all_parcels = sum(special_parcels, [])
    invalid_parcels = [parcel for parcel in all_parcels if parcel not in VALID_PARCEL_CODES]
    if invalid_parcels:
        raise ValueError(f"Invalid parcel code(s) provided in 'special_parcels': {invalid_parcels}\n\n"
            f"Valid parcel codes are:\n{list(VALID_PARCEL_CODES)}\n\n Parcel codes include parcel type ('mu', 'ml', 'sb'), adiabatic scheme ('ps', 'ia'), and CAPE type ('cape', 'ecape').")
    return all_parcels

def _calculate_special_parcel(parcel, data, snd_params):
    # calc special ecape parcels with ecape_parcel library 
    kinem = snd_params["kinem"]
    parcel_type, adiabatic_type, cape_type = parcel.split("_")
    trace = calc_ecape_parcel(data["p"], data["z"], data["T"], data["Td"], data["u"], data["v"], False,
        cape_type=PARCEL_DETAILS["term1"][parcel_type], pseudoadiabatic_switch=PARCEL_DETAILS["term2"][adiabatic_type][0],
        entrainment_switch=PARCEL_DETAILS["term3"][cape_type], storm_motion_type="user_defined",
        storm_motion_u=kinem["sm_u"] * units.kts, storm_motion_v=kinem["sm_v"] * units.kts)

    if trace[0][0] is None:
        return trace
    trace_Trho = density_temperature(trace[2], trace[3], trace[4])

    return trace + (trace_Trho,)

def _calculate_special_parcels(special_parcels, data, snd_params):
    # validate and compute special parcels
    all_parcels = _validate_special_parcels(special_parcels)
    parcel_data = {}
    for parcel in all_parcels:
        parcel_data[parcel] = _calculate_special_parcel(parcel, data, snd_params)
    return parcel_data

def _format_special_parcel_name(parcel):
    # combine term components into a name for a special parcel
    parcel_type, adiabatic_type, cape_type = parcel.split("_")
    adiabatic_name = PARCEL_DETAILS["term2"][adiabatic_type][1]
    return (f"{parcel_type} {adiabatic_name} {cape_type}").upper()



##################################################################################
############################ SOUNDING PLOT UTILS #################################
################################################################################## 



# apply theme background and border colors to an axis
def _style_axis(ax, theme):
    ax.set_facecolor(theme["background"])
    for spine in ax.spines.values():
        spine.set_color(theme["border"])


# resize every font and line on a finished figure
def _scale_figure(fig, font_scale=1, line_scale=1):

    if font_scale <= 0 or line_scale <= 0:
        raise ValueError(f"'font_scale' and 'line_scale' must be above 0. Got: {font_scale}, {line_scale}")
    if font_scale == 1 and line_scale == 1:
        return

    # ticks and legends are only complete once the figure has been laid out
    fig.draw_without_rendering()

    for artist in {id(artist): artist for artist in fig.findobj()}.values():
        if isinstance(artist, Text):
            artist.set_fontsize(artist.get_fontsize() * font_scale)
        elif isinstance(artist, Line2D):
            artist.set_linewidth(artist.get_linewidth() * line_scale)
            artist.set_markersize(artist.get_markersize() * line_scale)
        elif isinstance(artist, (Collection, Patch)):
            artist.set_linewidth(np.asarray(artist.get_linewidth()) * line_scale)


# determine skew-t x-axis bounds based on sfc t
def _get_skewt_temperature_bounds(surface_temperature):
    sfc_temp = surface_temperature.to("degC").m
    if sfc_temp <= -20:
        return (-52, 32)
    if sfc_temp <= -10:
        return (-42, 42)
    if sfc_temp >= 20:
        return (-22, 62)
    return (-32, 52)


# create skewt fig
def _setup_skewt(fig, temperature, theme, rect=None, subplot=None, rotation=47, 
                 pressure_limits=(1050, 100), temperature_limits="auto", box_aspect=0.87, show_pressure_ticks=True, pressure_ticks=(1000, 900, 800, 700, 600, 500, 400, 300, 200)):

    skew = SkewT(fig, rotation=rotation, subplot=subplot, rect=rect)
    skew.ax.set_box_aspect(box_aspect)
    skew.ax.zorder = 5
    skew.ax.set_adjustable("datalim")
    skew.ax.set_ylim(*pressure_limits)
    if temperature_limits == "auto":
        xmin, xmax = _get_skewt_temperature_bounds(temperature[0])
    else:
        xmin, xmax = temperature_limits
    skew.ax.set_xlim(xmin, xmax)
    skew.ax.set_xlabel("")
    skew.ax.set_ylabel("")
    skew.ax.tick_params(axis="x", direction="in", pad=-12, colors=theme["text"], labelsize=FONT_SIZES["tick"])
    skew.ax.tick_params(axis="y", direction="in", pad=-10, colors=theme["text"], labelsize=FONT_SIZES["tick"])
    if show_pressure_ticks:
        pressure_ticks = list(pressure_ticks)
        skew.ax.set_yticks(pressure_ticks)
        skew.ax.set_yticklabels(pressure_ticks, color=theme["text"], fontweight="bold")
        for label in skew.ax.get_yticklabels():
            label.set_horizontalalignment("left")
    _style_axis(skew.ax, theme)

    return skew



def _plot_skewt_accessories(skew, theme, shade_isotherms=True, highlight_isotherms=(0, -20), 
                            dry_adiabats=True, moist_adiabats=True, mixing_lines=True):

    if shade_isotherms:
        x1 = np.linspace(-100, 40, 8)
        x2 = np.linspace(-90, 50, 8)
        y = [1200, 50]
        for x1_i, x2_i in zip(x1, x2):
            skew.shade_area(y=y, x1=x1_i, x2=x2_i, color=theme["shade"], alpha=theme["shade_alpha"], zorder=1)

    for temp in highlight_isotherms:
        skew.ax.axvline(temp * units.degC, linestyle="--", color=theme["isotherm"], alpha=0.3)
    if dry_adiabats:
        skew.plot_dry_adiabats(color=theme["adiabat"], linewidth=0.3, alpha=0.8)
    if moist_adiabats:
        skew.plot_moist_adiabats(color=theme["adiabat"], linewidth=0.3, alpha=0.8)
    if mixing_lines:
        skew.plot_mixing_lines(color=theme["skew_line"], linewidth=0.3, alpha=0.8)




def _plot_skewt_profiles(skew, data, snd_params, theme, show_temperature=True, 
                         show_dewpoint=True, show_wetbulb=True, show_virtual=True, linewidth=5):

    p = data["p"]

    if show_wetbulb:
        skew.plot(p, snd_params["general"]["wet_bulb"], theme["wetbulb"], linewidth=2, label="WETBULB TEMP", alpha=0.6,)
    if show_virtual:
        skew.plot(p, snd_params["general"]["virt_temp"], theme["virtual"], linestyle=(0, (1, 1)), linewidth=4, label="VIRTUAL TEMP", alpha=0.7)
    if show_dewpoint:
        skew.plot(p, data["Td"], theme["dewpoint"], linewidth=linewidth,label="DEWPOINT")
    if show_temperature:
        skew.plot(p, data["T"], theme["temperature"], linewidth=linewidth, label="TEMPERATURE")



def _plot_parcel_trace(skew, pressure, temperature, color, label=None, linewidth=2, linestyle="--", alpha=1):
    return skew.plot(pressure, temperature, color=color, label=label, linewidth=linewidth, linestyle=linestyle, alpha=alpha)




def _plot_special_parcels(skew, special_parcels, parcel_data, theme):
    # plot user requested special parcel traces
    primary_parcels = special_parcels[0]
    secondary_parcels = special_parcels[1]

    for parcel, trace in parcel_data.items():

        parcel_name = _format_special_parcel_name(parcel)

        if parcel in primary_parcels:
            if len(trace) == 6:
                _plot_parcel_trace(skew, trace[0], trace[-1], color=theme["special_parcel"], linestyle="--", linewidth=2, alpha=1, label=parcel_name)
        elif parcel in secondary_parcels:
            if len(trace) == 6:
                _plot_parcel_trace(skew, trace[0], trace[-1], color=theme["secondary_parcel"], linestyle="--", linewidth=2, alpha=0.5, label=parcel_name)




def _plot_basic_parcels(skew, snd_params, theme, muecape_trace=None):
    # plot standard MU/SB/ML parcel traces. muecape_trace comes from _calculate_special_parcel
    thermo = snd_params["thermo"]

    # MUECAPE parcel
    if muecape_trace is not None and len(muecape_trace) == 6:
        skew.plot(muecape_trace[0], muecape_trace[-1], linestyle="--", linewidth=3, alpha=1, color=theme["special_parcel"], label="MUECAPE PARCEL")

    for key in ("sb", "ml", "mu"):
        if thermo[f"{key}cape"] > 0:
            skew.plot(thermo[f"{key}P_trace"], thermo[f"{key}T_trace"], linestyle="--", linewidth=2, alpha=0.75,
                      color=theme["parcel"][key], label=f"{key.upper()}CAPE PARCEL")



def _plot_downdraft_parcel(skew, snd_params, theme):
    # plot downdraft parcel trace
    skew.plot(snd_params["thermo"]["dparcel_p"], snd_params["thermo"]["dparcel_T"], 
              linestyle="--", linewidth=0.75, color=theme["downdraft"], alpha=0.8, label="DWNDRFT PARCEL")


def _add_skewt_layer_line(skew, x, p_bottom, p_top, color, alpha=0.7, linewidth=1, linestyle=":"):
    line = mlines.Line2D([x, x], [p_bottom, p_top], color=color,
        linestyle=linestyle, linewidth=linewidth, alpha=alpha, transform=skew.ax.get_yaxis_transform(),)
    skew.ax.add_artist(line)


# km marker heights (m AGL) that suit a skew-t's pressure range
def _get_height_markers(snd_params, ylim):
    intrp = snd_params["intrp"]
    in_view = intrp["zINTRP"][(intrp["pINTRP"] <= ylim[0]) & (intrp["pINTRP"] >= ylim[1])]
    if in_view.size == 0:
        return []
    depth = in_view.max() - in_view.min()
    step = 500 if depth <= 1500 else 1000 if depth <= 6000 else 2000
    return [height for height in np.arange(step, in_view.max() + 1, step) if height >= in_view.min()]


def _plot_skewt_height_annotations(skew, data, snd_params, theme, heights=range(1000, 14000, 2000), xloc=0.048, show_surface=True):
    # plots height markers (km AGL) to left of skew-t
    p = data["p"]
    intrp = snd_params["intrp"]
    trans, _, _ = skew.ax.get_yaxis_text1_transform(0)
    for height in heights:
        # shallow profiles may not reach every height
        if height > np.nanmax(intrp["zINTRP"]):
            continue
        skew.ax.text(xloc, intrp["pINTRP"][find_nearest(intrp["zINTRP"], height)], f"{height / 1000:g}km",
            fontsize=FONT_SIZES["note"], transform=trans, alpha=0.6, weight="bold", color=theme["text"])
    if show_surface:
        # surface elevation
        skew.ax.text(xloc, p[0], f"-SFC ({mag(snd_params['general']['elevation'])}m) -",
            fontsize=FONT_SIZES["note"],transform=trans, alpha=0.6, weight="bold", color=theme["text"])



def _plot_skewt_surface_annotations(skew, data, theme, surface_modified=False, temperature_unit="degF"):
    p = data["p"]
    T = data["T"]
    Td = data["Td"]
    # convert surface T/Td to the requested unit
    unit_label = f"°{temperature_unit[-1]}"
    # nearest whole degree, with halves rounded up (72.5 -> 73)
    T_sfc = int(np.floor(T[0].to(units(temperature_unit)).m + 0.5))
    Td_sfc = int(np.floor(Td[0].to(units(temperature_unit)).m + 0.5))
    # a faded background keeps the surface labels readable where they land on the axis labels
    # surface temperature
    skew.ax.annotate(f"{T_sfc}{unit_label}", (T[0], p[0]), textcoords="offset points", xytext=(16, -15),
                      fontsize=FONT_SIZES["annotation_large"], color=theme["temperature"], weight="bold", alpha=0.9, ha="center", zorder=10,
                      bbox=dict(facecolor=theme["background"], edgecolor="none", alpha=0.6, pad=1.5))
    # surface dewpoint
    skew.ax.annotate(f"{Td_sfc}{unit_label}", (Td[0], p[0]), textcoords="offset points", xytext=(-16, -15),
                      fontsize=FONT_SIZES["annotation_large"], color=theme["dewpoint"], weight="bold", alpha=0.9, ha="center", zorder=10,
                      bbox=dict(facecolor=theme["background"], edgecolor="none", alpha=0.6, pad=1.5))
    # user-modified surface flag
    if surface_modified:
        skew.ax.annotate("USER MODIFIED SFC →", (Td[0], p[0]), textcoords="offset points", xytext=(-50, 0),
                          fontsize=FONT_SIZES["small"], color=theme["text"], weight="bold", alpha=0.5, ha="right", zorder=10)



# space labels at least `gap` pixels apart, keeping each as close to its own pressure as it can
def _spread_label_pressures(ax, pressures, gap):
    to_display, to_data = ax.transData, ax.transData.inverted()
    ys = [to_display.transform((0, p))[1] for p in pressures]
    top = ax.bbox.y1 - gap / 2
    placed = list(ys)
    order = np.argsort(ys)
    # push up from the bottom, then back down from the top of the plot
    for prev, idx in zip(order[:-1], order[1:]):
        placed[idx] = max(placed[idx], placed[prev] + gap)
    placed[order[-1]] = min(placed[order[-1]], top)
    for above, idx in zip(order[::-1][:-1], order[::-1][1:]):
        placed[idx] = min(placed[idx], placed[above] - gap)
    return [to_data.transform((0, y))[1] for y in placed]



def _plot_skewt_parcel_level_annotations(skew, data, snd_params, theme, parcel="mu", levels=("lcl", "lfc", "el"), x=0.85):
    # a column for one parcel: its name at the surface, then its levels up a dotted line. Labels
    # that are too close are spread apart, and an arrow marks each level's true pressure
    thermo = snd_params["thermo"]
    color = theme["parcel"][parcel]
    ax = skew.ax
    transform = ax.get_yaxis_transform()
    p_bottom, p_top = ax.get_ylim()
    fontsize = FONT_SIZES["table"]
    gap = fontsize * 1.2 * ax.figure.dpi / 72

    # with no CAPE, SHARPpy puts the LFC at the profile top and the EL/MPL at the LCL, so only the LCL is labeled
    cape = thermo[f"{parcel}cape"]
    if ma.is_masked(cape) or not cape > 0:
        levels = tuple(level for level in levels if level == "lcl")

    # (text, true pressure); a level above the plot is labeled at the top with an up arrow,
    # and a level below a zoomed plot is left off
    labels = [(parcel.upper(), min(data["p"][0].m, p_bottom))]
    for level in levels:
        pressure = thermo[f"{parcel}_{level}_p"]
        if ma.is_masked(pressure) or pressure > p_bottom:
            continue
        if pressure < p_top:
            labels.append((f"↑{level.upper()}", p_top))
        else:
            labels.append((level.upper(), float(pressure)))
    if len(labels) == 1:
        return

    label_pressures = _spread_label_pressures(ax, [p for _, p in labels], gap)
    renderer = ax.figure.canvas.get_renderer()
    half_width = 0
    for i, ((text, pressure), label_p) in enumerate(zip(labels, label_pressures)):
        is_name = i == 0
        label = ax.text(x, label_p, text, weight="bold", color=color, alpha=0.7 if is_name else 0.95, ha="center", va="center",
                        fontsize=FONT_SIZES["label"] if is_name else fontsize, transform=transform, clip_on=True)
        half_width = max(half_width, label.get_window_extent(renderer).width / 2 / ax.bbox.width)

    # an arrow left of the column at each level's true pressure, since a spread label may sit above it
    for text, pressure in labels[1:]:
        if not text.startswith("↑"):
            arrow = ax.annotate("", xy=(x - half_width - 0.016, pressure), xycoords=transform, xytext=(x - half_width - 0.003, pressure),
                                textcoords=transform, arrowprops=dict(arrowstyle="-|>", color=color, alpha=0.95, linewidth=1.5,
                                                                      mutation_scale=10, shrinkA=0, shrinkB=0))
            arrow.arrow_patch.set_clip_path(ax.patch)

    # dotted line up the column, broken around each label
    ys = sorted(ax.transData.transform((0, p))[1] for p in label_pressures)
    to_data = ax.transData.inverted()
    for lower, upper in zip(ys[:-1], ys[1:]):
        if upper - lower > gap:
            _add_skewt_layer_line(skew, x=x, p_bottom=to_data.transform((0, lower + gap / 2))[1],
                                  p_top=to_data.transform((0, upper - gap / 2))[1], color=color, alpha=0.6)



def _plot_skewt_freezing_annotations(skew, snd_params, theme, xloc=0.745, show_temperature=True, show_wetbulb=True):
    # label the freezing and wet-bulb-zero levels, spreading the labels apart when they're close
    general = snd_params["general"]
    labels = []
    for show, key, text in ((show_temperature, "frz_pt", "←FRZ"), (show_wetbulb, "wb_frz_pt", "←WB0")):
        if show and not ma.is_masked(general[f"{key}_z"]) and general[f"{key}_z"] >= 50 * units.m:
            labels.append((text, getattr(general[f"{key}_p"], "m", general[f"{key}_p"])))
    if not labels:
        return

    gap = FONT_SIZES["label"] * 1.2 * skew.ax.figure.dpi / 72
    for (text, _), pressure in zip(labels, _spread_label_pressures(skew.ax, [p for _, p in labels], gap)):
        skew.ax.text(xloc, pressure, text, weight="bold", color=theme["cool"], alpha=0.6, va="center",
                     fontsize=FONT_SIZES["label"], transform=skew.ax.get_yaxis_transform(), clip_on=True)



# whether there is an effective inflow layer with some depth to draw
def _has_eil(kinem):
    eil = kinem.get("eil", (ma.masked, ma.masked))
    return not ma.is_masked(eil[0]) and not ma.is_masked(eil[1]) and eil[0] > eil[1]


def _plot_skewt_layer_annotations(skew, data, snd_params, theme, layers=("pbl", "lapse_rate", "eil", "growth_zone")):
    # plot PBL, max lapse rate, EIL, and DGZ/HGZ annotations on skew-t
    thermo = snd_params["thermo"]
    kinem = snd_params["kinem"]
    T = data["T"]
    p = data["p"]
    transform = skew.ax.get_yaxis_transform()

    # short horizontal ticks at a layer's bottom and top
    def _plot_layer_ticks(x, bounds, color, alpha, linewidth):
        for bound in bounds:
            if not ma.is_masked(bound):
                skew.ax.add_artist(mlines.Line2D([x - 0.01, x + 0.01], [bound, bound], color=color, alpha=alpha,
                                                 linewidth=linewidth, transform=transform))

    def _plot_labeled_layer(x, bounds, label, color, label_alpha, line_alpha, linewidth=3, ticks=False):
        # sit the label higher when a tick marks the layer top; a multi-line label is anchored at its
        # bottom so the extra lines stack upward instead of running down into the layer line
        va = "bottom" if "\n" in label else "baseline"
        skew.ax.text(x, bounds[1] - (12 if ticks else 8), label, weight="bold", color=color, alpha=label_alpha, ha="center", va=va,
                     linespacing=1.1, fontsize=FONT_SIZES["tick"], transform=transform)
        _add_skewt_layer_line(skew, x=x, p_bottom=bounds[0], p_top=bounds[1], color=color, alpha=line_alpha, linewidth=linewidth)
        if ticks:
            _plot_layer_ticks(x, bounds, color, line_alpha, linewidth - 1)

    if "pbl" in layers:
        skew.ax.text(0.76, thermo["pbl_top"] - 4, "←PBL", weight="bold", color=theme["pbl"], alpha=0.95, fontsize=FONT_SIZES["small"], transform=transform, clip_on=True)
        _add_skewt_layer_line(skew, x=0.78, p_bottom=p[0], p_top=thermo["pbl_top"], color=theme["pbl"], alpha=0.95)
        _plot_layer_ticks(0.78, (p[0], thermo["pbl_top"]), theme["pbl"], 0.95, 1.5)

    if "lapse_rate" in layers:
        # value over unit on two lines, right of the DGZ/HGZ bar so the label doesn't run into it
        _plot_labeled_layer(0.175, thermo["lr_max"][1:3], f"{np.round(thermo['lr_max'][0], 1)}\n°C/km", theme["lapse_max"], 0.8, 0.6)

    if "eil" in layers and _has_eil(kinem):
        _plot_labeled_layer(0.19, kinem["eil"], "EIL", theme["eil"], 1, 1, linewidth=4, ticks=True)

    if "growth_zone" in layers:
        if T[0].m < 5:
            if thermo["dgz"][1] < p[0].m:
                _plot_labeled_layer(0.13, thermo["dgz"], "DGZ", theme["dgz"], 0.6, 0.4, ticks=True)
        elif thermo["hgz"][1] < p[0].m:
            _plot_labeled_layer(0.13, thermo["hgz"], "HGZ", theme["hgz"], 0.6, 0.4, ticks=True)



def _plot_skewt_cape_annotations(skew, snd_params, theme):
    # plot MUCAPE value annotations at 3km, 6km, and the EL along the parcel trace
    thermo = snd_params["thermo"]
    intrp = snd_params["intrp"]
    color = theme["parcel"]["mu"]

    def _plot_cape_label(value, temperature, pressure):
        label = f"← {mag(value)} " + r"$\mathrm{Jkg^{-1}}$"
        skew.ax.annotate(label, xy=(temperature, pressure), xytext=(8, 0), textcoords="offset points",
                         color=color, alpha=0.7, fontsize=FONT_SIZES["annotation"], ha="left", va="center")

    # parcel temperature where the trace reaches a height (m AGL) in the sounding
    def _plot_cape_label_at_height(value, height):
        if height > np.nanmax(intrp["zINTRP"]):
            return
        pressure = intrp["pINTRP"][find_nearest(intrp["zINTRP"], height)]
        trace_p = np.ma.filled(thermo["muP_trace"], np.nan)
        trace_T = np.ma.filled(thermo["muT_trace"], np.nan)
        valid = np.isfinite(trace_p) & np.isfinite(trace_T)
        # the trace may not start until above this height, if the parcel is elevated
        if pressure > trace_p[valid].max():
            return
        temperature = np.interp(-np.log(pressure), -np.log(trace_p[valid]), trace_T[valid])
        _plot_cape_label(value, temperature, pressure)

    # 0-3 km MUCAPE
    if thermo["mu3cape"] > 10:
        _plot_cape_label_at_height(thermo["mu3cape"], 3000)

    # 0-6 km MUCAPE
    if thermo["mu6cape"] > thermo["mu3cape"]:
        _plot_cape_label_at_height(thermo["mu6cape"], 6000)

    # total MUCAPE
    if thermo["mucape"] > thermo["mu6cape"]:
        _plot_cape_label(thermo["mucape"], thermo["mu_el_T"], thermo["mu_el_p"])



def _plot_skewt_rh_annotations(skew, snd_params, theme):
    # plot RH annotations every 3km along the dewpoint trace
    height_indices = [10, 30, 60, 90, 120]
    for idx in height_indices:
        # shallow profiles may not reach every height
        if idx >= len(snd_params["intrp"]["rhINTRP"]):
            continue
        rh = snd_params["intrp"]["rhINTRP"][idx]
        if ma.is_masked(rh) or np.isnan(rh):
            continue
        label = f"{int(rh)}% →"
        dewpoint = mpcalc.dewpoint_from_relative_humidity(snd_params["intrp"]["tINTRP"][idx] * units.degC, rh / 100)
        skew.ax.annotate(label, (dewpoint, snd_params["intrp"]["pINTRP"][idx] * units.hPa), textcoords="offset points", 
                         xytext=(-40, 0), color=theme["dewpoint"], alpha=0.7, fontsize=FONT_SIZES["annotation"], ha="center", clip_on=True)




def _plot_skewt_wind_barbs(skew, data, snd_params, theme, add_sr_barbs=True, 
                           ground_xloc=0.955, storm_xloc=None, n_barbs=30, barb_length=8, sr_alpha=0.4, show_labels=True, ylim=None):

    p = data["p"]
    u = data["u"]
    v = data["v"]
    barb_color = theme["barb"]
    barb_sizes = {"emptybarb": 0.075, "width": 0.18, "height": 0.4}
    transform = skew.ax.get_yaxis_transform()

    # overlay SR barbs on GR barbs unless another location is given
    if storm_xloc is None:
        storm_xloc = ground_xloc
    # resample winds to logarithmically spaced pressure levels
    # and with a pressure range, space them within it
    if ylim is None:
        interval = np.logspace(2.113, 3, n_barbs) * units.hPa
        label_pressure = 115
    else:
        interval = np.logspace(np.log10(ylim[1]) + 0.113 * np.log10(ylim[0] / ylim[1]) / 1.021, np.log10(ylim[0]), n_barbs) * units.hPa
        label_pressure = ylim[1] * (ylim[0] / ylim[1]) ** 0.0594
    idx = mpcalc.resample_nn_1d(p, interval)
    if ylim is not None:
        idx = [i for i in idx if ylim[1] <= p[i].m <= ylim[0]]

    # blank barbs produce origin dots
    blank = np.zeros(len(idx))
    skew.plot_barbs(pressure=p[idx], u=blank, v=blank, xloc=ground_xloc, fill_empty=True, 
                    color=barb_color, sizes=barb_sizes)
    skew.plot_barbs(pressure=p[idx], u=u[idx], v=v[idx], xloc=ground_xloc, fill_empty=True, 
                    color=barb_color, sizes=barb_sizes, length=barb_length)

    plot_sr = (add_sr_barbs and not ma.is_masked(snd_params["kinem"]["sm_u"]) and not ma.is_masked(snd_params["kinem"]["sm_v"]))
    if plot_sr:
        sr_u = u[idx] - (snd_params["kinem"]["sm_u"] * units.kts)
        sr_v = v[idx] - (snd_params["kinem"]["sm_v"] * units.kts)

        # add origin dots if SR barbs use their own column
        if storm_xloc != ground_xloc:
            skew.plot_barbs(pressure=p[idx], u=blank, v=blank, xloc=storm_xloc, fill_empty=True, 
                            color=barb_color, alpha=sr_alpha, sizes=barb_sizes)
        skew.plot_barbs(pressure=p[idx], u=sr_u, v=sr_v, xloc=storm_xloc, fill_empty=True, 
                        color=barb_color, alpha=sr_alpha, sizes=barb_sizes, length=barb_length)

    barb_xlocs = [ground_xloc]
    if plot_sr and storm_xloc != ground_xloc:
        barb_xlocs.append(storm_xloc)
    for xloc in barb_xlocs:
        line = mlines.Line2D([xloc, xloc], [0.01, 0.95], color=barb_color, 
                             linewidth=0.5, transform=skew.ax.transAxes, clip_on=False, zorder=1)
        skew.ax.add_line(line)

    if show_labels:

        # when the two barb sets overlap, offset labels to either side
        if plot_sr and storm_xloc == ground_xloc:
            ground_label_x = ground_xloc - 0.025
            storm_label_x = storm_xloc + 0.025
        # when barbs have separate columns, center labels on each column
        else:
            ground_label_x = ground_xloc
            storm_label_x = storm_xloc

        skew.ax.text(ground_label_x, label_pressure, "Ground\nRelative", weight="bold",
                      color=barb_color, ha="center", alpha=1, fontsize=FONT_SIZES["tiny"], transform=transform, clip_on=True)
        if plot_sr:
            skew.ax.text(storm_label_x, label_pressure, "Storm\nRelative", weight="bold", color=barb_color, 
                         ha="center", alpha=sr_alpha, fontsize=FONT_SIZES["tiny"], transform=transform, clip_on=True)






#################################################################################
############################ HODOGRAPH PLOTTING #################################
################################################################################# 


def _prepare_hodograph_data(snd_params, max_height=12000):
    # use interp u, v, & z for easy hodo plotting
    z = np.asarray(snd_params["intrp"]["zINTRP"])
    u = np.asarray(snd_params["intrp"]["uINTRP"])
    v = np.asarray(snd_params["intrp"]["vINTRP"])
    valid = (np.isfinite(z) & np.isfinite(u) & np.isfinite(v) & (z <= max_height))
    return {"z": z[valid], "u": u[valid], "v": v[valid]}


# determine hodograph bounds based on u/v wind data
def _get_hodograph_bounds(hodo_data, lower_padding=45, upper_padding=30, fallback=(-100, 100), dynamic=True):

    if dynamic:
        u = hodo_data["u"]
        v = hodo_data["v"]

        if len(u) == 0 or len(v) == 0:
            return fallback, fallback

        x_min = min(np.nanmin(u), 0)
        x_max = max(np.nanmax(u), 0)
        y_min = min(np.nanmin(v), 0)
        y_max = max(np.nanmax(v), 0)

        xlim = (x_min - lower_padding, x_max + upper_padding)
        ylim = (y_min - lower_padding, y_max + upper_padding)
        return xlim, ylim
    else:
        return fallback, fallback


# create hodograph axis
def _setup_hodograph(ax, theme, xlim, ylim, 
                     component_range=160, box_aspect=1, major_grid=20, minor_grid=10):

    # create hodograph object
    hodo = Hodograph(ax, component_range=component_range)

    # set bounds
    hodo.ax.set_xlim(*xlim)
    hodo.ax.set_ylim(*ylim)

    # major grid
    if major_grid is not None:
        hodo.add_grid(increment=major_grid, color=theme["text"], linestyle="-", linewidth=1.5, alpha=0.4)

    # minor grid
    if minor_grid is not None:
        hodo.add_grid(increment=minor_grid, color=theme["text"], linestyle="--", linewidth=1, alpha=0.3)

    # background and borders
    _style_axis(hodo.ax, theme)

    # geometry
    hodo.ax.set_box_aspect(box_aspect)

    # remove standard ticks/labels
    hodo.ax.set_xticks([])
    hodo.ax.set_yticks([])
    hodo.ax.set_xlabel("")
    hodo.ax.set_ylabel("")

    return hodo


# add velocity labels to hodograph axes
def _plot_hodograph_velocity_labels(hodo, theme, speeds=range(10, 130, 20), fontsize=FONT_SIZES["label"], alpha=0.5):

    for speed in speeds:
        locations = (( speed, 0), (-speed, 0), (0,  speed), (0, -speed))
        for x, y in locations:
            hodo.ax.annotate(f"{str(speed)}kt", (x, y), xytext=(0, 0), textcoords="offset pixels", clip_on=True,
                fontsize=fontsize, weight="bold", alpha=alpha, zorder=0, color=theme["text"])


# add height markers to hodograph trace
def _plot_hodograph_height_markers(hodo, snd_params, theme,
        heights=(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000), marker_size=40, fontsize=FONT_SIZES["tick"]):

    z = np.asarray(snd_params["intrp"]["zINTRP"])
    u = np.asarray(snd_params["intrp"]["uINTRP"])
    v = np.asarray(snd_params["intrp"]["vINTRP"])
    max_height = np.nanmax(z)

    for height in heights:
        if height > max_height:
            continue
        idx = find_nearest(z, height)

        # label 500 m as ".5"; all others as integer km
        if height == 500:
            label = ".5"
        else:
            label = str(int(height / 1000))

        hodo.ax.plot(u[idx], v[idx], marker=".", color=theme["height_marker"], markeredgecolor=theme["height_marker_edge"], 
                     alpha=1, markersize=marker_size, zorder=5, clip_on=True)
        hodo.ax.annotate(label, (u[idx], v[idx]), weight="bold", fontsize=fontsize, color=theme["height_marker_edge"], 
                         xytext=(0.02, -4), textcoords="offset pixels", horizontalalignment="center", clip_on=True, zorder=6)



# plot hodograph
def _plot_hodograph_trace(hodo, snd_params, theme, layers=None, linewidth=7, zorder=4):

    z = np.asarray(snd_params["intrp"]["zINTRP"])
    u = np.asarray(snd_params["intrp"]["uINTRP"])
    v = np.asarray(snd_params["intrp"]["vINTRP"])

    if layers is None:
        layers = theme["hodo_layers"]

    max_height = np.nanmax(z)
    
    for bottom, top, color in layers:
        if bottom > max_height:
            continue
        top = min(top, max_height)
        idx_bottom = find_nearest(z, bottom)
        idx_top = find_nearest(z, top)

        hodo.ax.plot(u[idx_bottom:idx_top + 1], v[idx_bottom:idx_top + 1], 
                  color=color, linewidth=linewidth, zorder=zorder, clip_on=True)



# format a storm motion's direction and speed for a readout
def _format_motion(label, u, v):
    if not np.isfinite([u, v]).all():
        return f"{label}: ---"
    speed = mpcalc.wind_speed(u * units.kts, v * units.kts).m
    return f"{label}: {wind_to_dir(u, v)} @ {mag(speed)} kt"


# label storm motions at their locations on the hodograph
def _plot_hodograph_motion_labels(hodo, snd_params, storm_motion, theme, show_rm=True, show_lm=True, show_mw=True,
                                  show_dtm=True, show_mcs=True, show_selected_sm=True):
    kinem = snd_params["kinem"]
    text_color = theme["text"]

    # RM/LM AND MEAN WIND - BUNKERS
    if not ma.is_masked(kinem["sm_rm"]):
        for show, key, label in ((show_rm, "sm_rm", "RM"), (show_lm, "sm_lm", "LM"), (show_mw, "sm_mw", "MW")):
            if show:
                hodo.ax.text(kinem[key][0] + 0.5, kinem[key][1] - 0.5, label,
                    weight="bold", ha="left", fontsize=FONT_SIZES["table"], zorder=7, alpha=0.9, color=text_color)

    # MCS MOTION - CORFIDI
    if show_mcs and not ma.is_masked(kinem["mcs"]):
        hodo.ax.text(kinem["mcs"][0], kinem["mcs"][1], "UP", weight="bold", fontsize=FONT_SIZES["label"], color=theme["mcs"], zorder=7, ha="center", alpha=0.5, clip_on=True)
        hodo.ax.text(kinem["mcs"][2], kinem["mcs"][3], "DN", weight="bold", fontsize=FONT_SIZES["label"], color=theme["mcs"], zorder=7, ha="center", alpha=0.5, clip_on=True)

    # DEVIANT TORNADO - NIXON AND ALLEN
    if show_dtm and not ma.is_masked(kinem["dtm"]):
        hodo.ax.text(kinem["dtm"][0], kinem["dtm"][1] + 2, "DTM", weight="bold", fontsize=FONT_SIZES["small"], zorder=7, color=theme["dtm"], ha="center")
        hodo.ax.plot(kinem["dtm"][0], kinem["dtm"][1], marker="v", color=theme["dtm"], markersize=8, zorder=7, alpha=0.8, linestyle="")

    # USER DEFINED STORM MOTION
    valid_sm = not ma.is_masked(kinem["sm_u"]) and not ma.is_masked(kinem["sm_v"])
    if valid_sm and show_selected_sm and isinstance(storm_motion, list):
        hodo.ax.text(kinem["sm_u"] + 0.5, kinem["sm_v"] - 0.5, "SM", weight="bold", ha="left", fontsize=FONT_SIZES["annotation_large"], alpha=0.9, zorder=7, color=text_color)


# draw the selected storm motion's vector from the origin
def _plot_hodograph_sm_vector(hodo, snd_params, theme):
    kinem = snd_params["kinem"]
    if not ma.is_masked(kinem["sm_u"]) and not ma.is_masked(kinem["sm_v"]):
        hodo.ax.arrow(0, 0, kinem["sm_u"] - 0.3, kinem["sm_v"] - 0.3, linewidth=3, color=theme["text"], alpha=0.2, length_includes_head=True, head_width=0.6)


# readout name for a storm motion, matching the names calc() accepts
def _storm_motion_name(storm_motion):
    return {"rm": "RIGHT MOVING", "lm": "LEFT MOVING", "mw": "MEAN WIND", "custom": "USER DEFINED"}[storm_motion_type(storm_motion)]


# list storm motion values on the hodograph, with the selected motion in bold at the top
def _plot_hodograph_motion_readout(hodo, snd_params, storm_motion, theme, show_rm=True, show_lm=True, show_mw=True, show_dtm=True,
                                   show_mcs=True, x=0.085, y=0.99, fontsize=FONT_SIZES["annotation_large"]):
    kinem = snd_params["kinem"]
    motions = []

    if not ma.is_masked(kinem["sm_rm"]):
        for show, key, label in ((show_rm, "sm_rm", "RM"), (show_lm, "sm_lm", "LM"), (show_mw, "sm_mw", "MW")):
            if show:
                motions.append(_format_motion(label, kinem[key][0], kinem[key][1]))
    if show_mcs and not ma.is_masked(kinem["mcs"]):
        motions.append(_format_motion("UP", kinem["mcs"][0], kinem["mcs"][1]))
        motions.append(_format_motion("DN", kinem["mcs"][2], kinem["mcs"][3]))
    if show_dtm and not ma.is_masked(kinem["dtm"]):
        motions.append(_format_motion("DTM", kinem["dtm"][0], kinem["dtm"][1]))

    style = dict(transform=hodo.ax.transAxes, ha="left", va="top", color=theme["text"], fontsize=fontsize, linespacing=1.35, alpha=0.7, zorder=8)

    if not ma.is_masked(kinem["sm_u"]) and not ma.is_masked(kinem["sm_v"]):
        name = _storm_motion_name(storm_motion)
        selected = _format_motion("SM", kinem["sm_u"], kinem["sm_v"]).replace("SM: ", f"SM: {name} | ")
        hodo.ax.text(x, y, selected, weight="bold", **style)
        motions.insert(0, "")
    if motions:
        hodo.ax.text(x, y, "\n".join(motions), **style)


# shift a sounding's winds and storm motions to be relative to the selected storm motion
def _make_storm_relative(snd_params):
    kinem = dict(snd_params["kinem"])
    intrp = dict(snd_params["intrp"])
    sm_u, sm_v = kinem["sm_u"], kinem["sm_v"]

    intrp["uINTRP"] = intrp["uINTRP"] - sm_u
    intrp["vINTRP"] = intrp["vINTRP"] - sm_v
    for key in ("sm_rm", "sm_lm", "sm_mw", "dtm"):
        if not ma.is_masked(kinem[key]):
            kinem[key] = (kinem[key][0] - sm_u, kinem[key][1] - sm_v)
    if not ma.is_masked(kinem["mcs"]):
        kinem["mcs"] = [value - shift for value, shift in zip(kinem["mcs"], (sm_u, sm_v, sm_u, sm_v))]
    kinem["sm_u"], kinem["sm_v"] = 0.0, 0.0

    return {**snd_params, "kinem": kinem, "intrp": intrp}


# add a SRH fill, based on SM or over a layer
def _plot_hodograph_srh(hodo, snd_params, theme, layer="auto", alpha=None, 
                        show_connectors=True, connector_width=2.3, connector_alpha=0.5):

    color = theme["eil"]
    kinem = snd_params["kinem"]
    intrp = snd_params["intrp"]

    # storm motion is required for SRH geometry
    valid_sm = (not ma.is_masked(kinem["sm_u"]) and not ma.is_masked(kinem["sm_v"]))

    if not valid_sm:
        return

    z = intrp["zINTRP"]
    u = intrp["uINTRP"]
    v = intrp["vINTRP"]
    sm_u = kinem["sm_u"]
    sm_v = kinem["sm_v"]

    eil_available = _has_eil(kinem)
    if layer == "auto":
        if eil_available:
            bottom = kinem["eil_z"][0]
            top = kinem["eil_z"][1]
            fill_alpha = 0.3 if alpha is None else alpha
            label = "EIL SRH"
        else:
            bottom = 0
            top = 3000
            fill_alpha = 0.1 if alpha is None else alpha
            label = "0-3 km SRH"

    elif layer == "eil":
        if not eil_available:
            return
        bottom = kinem["eil_z"][0]
        top = kinem["eil_z"][1]
        fill_alpha = 0.3 if alpha is None else alpha
        label = "EIL SRH"

    elif (isinstance(layer, (tuple, list)) and len(layer) == 2):
        bottom, top = layer
        fill_alpha = 0.2 if alpha is None else alpha
        label = (f"{bottom / 1000:g}-{top / 1000:g} km SRH")
    else:
        raise ValueError("'layer' must be 'auto', 'eil', or a (bottom, top) height tuple.")

    bottom_idx = find_nearest(z, bottom)
    top_idx = find_nearest(z, top,)
    # guarantee proper array order
    if bottom_idx > top_idx:
        bottom_idx, top_idx = top_idx, bottom_idx

    if show_connectors:
        hodo.ax.plot([sm_u, u[bottom_idx]], [sm_v, v[bottom_idx]], linestyle="-", linewidth=connector_width, alpha=connector_alpha, zorder=3, color=color)
        hodo.ax.plot([sm_u, u[top_idx]], [sm_v, v[top_idx]], linestyle="-", linewidth=connector_width, alpha=connector_alpha, zorder=3, color=color)

    fill_u = np.append(u[bottom_idx:top_idx + 1], sm_u)
    fill_v = np.append(v[bottom_idx:top_idx + 1], sm_v)
    hodo.ax.fill( fill_u, fill_v, color=color, alpha=fill_alpha, zorder=2, label=label)




def _plot_hodograph_boundary(hodo, hodo_boundary, snd_params, theme, 
                             length=200, show_mw_vector=True, show_lm_vector=False, linestyle="--", linewidth=2, alpha=0.3):

    if hodo_boundary is None:
        return

    kinem = snd_params["kinem"]
    angles = hodo_boundary.get("angle", [])
    colors = hodo_boundary.get("color", [])

    if len(angles) != len(colors):
        raise ValueError("'hodo_boundary' must contain the same number of angles and colors.")
    
    for angle, color in zip(angles, colors):
        # convert meteorological angle
        angle_rad = np.deg2rad(90 - angle)
        delta_u = length * np.cos(angle_rad)
        delta_v = length * np.sin(angle_rad)
        boundary_u = np.array([-delta_u, delta_u])
        boundary_v = np.array([-delta_v, delta_v ])
        hodo.ax.plot(boundary_u, boundary_v, color=color, linestyle=linestyle, linewidth=linewidth, alpha=alpha, label="Boundary Axis")

    for show, key in ((show_mw_vector, "sm_mw"), (show_lm_vector, "sm_lm")):
        if show and not ma.is_masked(kinem[key]):
            hodo.ax.arrow(0, 0, kinem[key][0] - 0.3, kinem[key][1] - 0.3,
                          linewidth=3, color=theme["text"], alpha=0.2, label="SM Vector", length_includes_head=True, head_width=0.6)








#################################################################################
############################## MAP INSET UTILS ##################################
################################################################################# 

# read map shapes once per process
@lru_cache(maxsize=None)
def _load_map_shapes(name):
    feature = {"states": cfeature.STATES.with_scale("10m"), "counties": USCOUNTIES}[name]
    return tuple(feature.geometries()), feature.crs


# get the map shapes within a (west, south, east, north) box
def _get_map_shapes(name, bounds):
    shapes, crs = _load_map_shapes(name)
    box = sgeom.box(*bounds)
    return [shape for shape in shapes if shape.intersects(box)], crs


# find the coarsest tile zoom that still gives the map axis a sharp image
def _get_satellite_zoom(map_ax, map_zoom, oversample=2, tile_size=256):
    pixels = min(map_ax.bbox.width, map_ax.bbox.height) * oversample
    lon_span = 2 * map_zoom
    return int(np.clip(np.ceil(np.log2(pixels * 360 / (tile_size * lon_span))), 1, 18))


def _setup_map_inset(map_ax, latlon, theme, map_zoom=2, show_satellite=True, satellite_zoom="auto", show_states=True, show_counties=True):
    # map_ax needs a cartopy projection (e.g., fig.add_axes(rect, projection=ccrs.PlateCarree()))
    lat, lon = latlon
    west, east, south, north = lon - map_zoom, lon + map_zoom, lat - map_zoom / 1.6, lat + map_zoom * 1.6
    bounds = (west - 0.5, south - 0.5, east + 0.5, north + 0.5)

    # map border
    map_ax.spines["geo"].set_color(theme["border"])

    # map extent
    map_ax.set_extent([west, east, south, north], crs=ccrs.PlateCarree())
    map_ax.set_box_aspect(1)

    if show_states:
        states, crs = _get_map_shapes("states", bounds)
        map_ax.add_geometries(states, crs, color=theme["background"], alpha=0.4, zorder=2)
        map_ax.add_geometries(states, crs, edgecolor=theme["text"], facecolor="none", alpha=0.5, linestyle="-", linewidth=1, zorder=6)
    if show_counties:
        counties, crs = _get_map_shapes("counties", bounds)
        map_ax.add_geometries(counties, crs, facecolor="none", edgecolor=theme["text"], alpha=0.5, linestyle="-", linewidth=0.5, zorder=5)
    if show_satellite:
        if satellite_zoom == "auto":
            satellite_zoom = _get_satellite_zoom(map_ax, map_zoom)
        satellite = img_tiles.GoogleTiles(style="satellite")
        map_ax.add_image(satellite, satellite_zoom, zorder=1)

    return map_ax


# add a boxed label to the bottom of the map
def _add_map_label(map_ax, text, theme):
    label = map_ax.text(0.50, 0.03, text, transform=map_ax.transAxes, ha="center", va="bottom", alpha=0.9, weight="bold",
                        fontsize=FONT_SIZES["small"], zorder=13, color=theme["text"], bbox=dict(facecolor=theme["background"], alpha=0.8, edgecolor=theme["border"], pad=3))
    _fit_text(label, map_ax, 0.92)


def _plot_map_mosaic_radar(map_ax, sounding_data, map_zoom, radar_time, theme, reflectivity_min=10,
                           alpha=1, lon_offset=0.05, lat_offset=0.05, show_timestamp=True):
    try:
        radar_data = get_radar_mosaic(sounding_data, map_zoom, radar_time)
        reflectivity = np.ma.masked_less(np.asarray(radar_data["Reflectivity"])[0], reflectivity_min)

        map_ax.pcolormesh(radar_data["lon"] + lon_offset, radar_data["lat"] + lat_offset, reflectivity, vmin=-32, vmax=95,
                          cmap=rs_expertreflect_cmap, alpha=alpha, zorder=4)

        if show_timestamp:
            scan_time = str(radar_data["time"].values[0])
            _add_map_label(map_ax, f"MOSAIC: {scan_time[5:7]}/{scan_time[8:10]} | {scan_time[11:16]}z", theme)
        return True

    except Exception as exc:
        warnings.warn(f"Radar mosaic could not be plotted: {exc}", Warning)
        if show_timestamp:
            _add_map_label(map_ax, "RADAR UNAVAILABLE", theme)
        return False



def _plot_map_single_site_radar(map_ax, sounding_data, radar_time, theme, reflectivity_min=10, show_timestamp=True):
    try:
        lat, lon = sounding_data["site_info"]["site-latlon"]
        site_ids, lat_radar, lon_radar = parse_nexrad_locs()
        nexrad_site, nexrad_lat, nexrad_lon = find_nearest_station(site_ids, lat_radar, lon_radar, lat, lon)
        radar_obj, scan_timestamp, nexrad_site = get_radar_data(nexrad_site, sounding_data, radar_time)
        radar_obj.longitude["data"][0] = nexrad_lon
        radar_obj.latitude["data"][0] = nexrad_lat

        gate_lat, gate_lon, _ = radar_obj.get_gate_lat_lon_alt(0)
        reflectivity = radar_obj.get_field(0, "reflectivity")
        reflectivity = np.ma.masked_where((reflectivity < reflectivity_min) | (reflectivity == 1) | (reflectivity == 2), reflectivity)

        # drop gates beyond the map, and group those too close together to see on it, keeping the strongest return
        west, east, south, north = map_ax.get_extent(ccrs.PlateCarree())
        in_map = (gate_lon >= west) & (gate_lon <= east) & (gate_lat >= south) & (gate_lat <= north)
        last_gate = np.nonzero(in_map.any(axis=0))[0].max() + 1
        pixel_km = (north - south) * 111 / map_ax.bbox.height
        gate_km = (radar_obj.range["data"][1] - radar_obj.range["data"][0]) / 1000
        step = max(int(pixel_km / gate_km), 1)
        n_groups = last_gate // step
        grouped = reflectivity[:, :n_groups * step].reshape(reflectivity.shape[0], n_groups, step).max(axis=2)

        map_ax.pcolormesh(gate_lon[:, step // 2::step][:, :n_groups], gate_lat[:, step // 2::step][:, :n_groups], grouped,
                          vmin=-32, vmax=95, cmap=rs_expertreflect_cmap, shading="auto", zorder=4)

        if show_timestamp:
            _add_map_label(map_ax, f"{nexrad_site}: {scan_timestamp}", theme)
        return True

    except Exception as exc:
        warnings.warn(f"Single-site radar could not be plotted: {exc}", Warning)
        if show_timestamp:
            _add_map_label(map_ax, "RADAR UNAVAILABLE", theme)
        return False



def _plot_map_profile_location(map_ax, latlon, theme, marker_size=90, center_size=16):

    lat, lon = latlon
    transform = ccrs.PlateCarree()
    map_ax.scatter(lon, lat, s=marker_size + 55, facecolor=theme["background"], edgecolor=theme["background"],
                    linewidth=2, alpha=0.6, transform=transform, zorder=9.9)
    map_ax.scatter(lon, lat, s=marker_size, facecolor="none", edgecolor=theme["marker"], linewidth=2, transform=transform, zorder=10)
    map_ax.scatter(lon, lat, s=center_size, facecolor=theme["marker"], edgecolor="none", transform=transform, zorder=11)




#################################################################################
########################### ACCESSORY PANEL UTILS ###############################
#################################################################################

# set up a log-pressure axis with no ticks, to sit alongside a skew-t
def _setup_pressure_axis(ax, theme, pressure_limits=(1050, 100)):

    ax.set_yscale("log")
    ax.set_ylim(*pressure_limits)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_yticks([], minor=True)
    _style_axis(ax, theme)

    return ax



def _plot_temperature_advection(ax, snd_params, theme, pressure_levels=range(200, 1100, 100), pressure_limits=(1050, 100),
                                x_padding=4, show_values=True):

    advection, bounds = snd_params["thermo"]["temp_adv"]
    advection = np.ma.filled(advection, np.nan)
    bottoms, tops = np.ma.filled(bounds, np.nan).T

    # skip layers without a value
    valid = np.isfinite(advection)
    advection, bottoms, tops = advection[valid], bottoms[valid], tops[valid]

    _setup_pressure_axis(ax, theme, pressure_limits)
    xlim = (advection.min() - x_padding, advection.max() + x_padding) if advection.size else (-x_padding, x_padding)
    ax.set_xlim(*xlim)

    ax.hlines(pressure_levels, *xlim, color=theme["pressure_line"], alpha=0.8, linewidth=1)
    ax.barh(tops, advection, height=bottoms - tops, align="edge", edgecolor=theme["bar_edge"], alpha=0.3,
            color=np.where(advection > 0, theme["warm"], theme["cool"]))
    ax.axvline(0, color=theme["text"], linewidth=1, linestyle="--")

    if show_values:
        for value, top in zip(advection, tops):
            if value != 0:
                ax.text(0.3 if value > 0 else -0.3, top + 10, f"{value:.1f}", ha="left" if value > 0 else "right",
                        color=theme["text"], weight="bold", clip_on=True)

    return ax



# overlay omega on a skew-t. omega (Pa/s) is scaled for plotting simplicity
def _plot_omega(ax, data, theme, x_limits=(40, -300), scale=10, 
                reference_value=10, pressure_limits=(1050, 100)):

    p = data["p"].m
    omega = np.asarray(getattr(data["omega"], "m", data["omega"])) * scale

    _setup_pressure_axis(ax, theme, pressure_limits)
    ax.set_xlim(*x_limits)
    ax.set_axis_off()

    ax.plot(omega, p, linewidth=2, color=theme["text"], alpha=0.1)
    ax.fill_betweenx(p, 0, omega, where=omega < 0, color=theme["ascent"], alpha=0.08, interpolate=True)
    ax.fill_betweenx(p, 0, omega, where=omega > 0, color=theme["cool"], alpha=0.08, interpolate=True)

    ax.axvline(0, color=theme["omega_zero"], alpha=0.2)
    for value in (reference_value, -reference_value):
        ax.axvline(value, color=theme["omega_reference"], linestyle="-.", alpha=0.2)
        ax.text(value, 150, f"{value:+d}", weight="bold", alpha=0.4, color=theme["text"], fontsize=FONT_SIZES["label"], ha="center", clip_on=True)

    return ax




# configure the x-axis of a low-level panel
def _set_panel_x_axis(ax, theme, limits, ticks, tick_labels, tick_pad=-12, rotation=0):

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, weight="bold", alpha=0.5, fontstyle="italic", color=theme["text"], rotation=rotation)
    ax.set_xlim(*limits)
    ax.tick_params(axis="x", direction="in", pad=tick_pad)


# set up a height axis for a low-level panel
def _setup_height_panel(ax, theme, x_limits, x_ticks, x_tick_labels, z_limits=(0, 3000), z_interval=500, tick_pad=-12, rotation=0):

    ax.set_yticks(np.arange(0, z_limits[1] + 1, z_interval))
    ax.set_ylim(*z_limits)
    ax.tick_params(axis="y", length=0, labelleft=False)
    ax.grid(True, axis="y")
    _set_panel_x_axis(ax, theme, x_limits, x_ticks, x_tick_labels, tick_pad=tick_pad, rotation=rotation)
    _style_axis(ax, theme)

    return ax


def _add_panel_title(ax, title, theme, fontsize=FONT_SIZES["label"]):
    ax.text(0.5, 0.98, title, transform=ax.transAxes, ha="center", va="top", weight="bold", color=theme["text"], fontsize=fontsize, alpha=0.9)


# shade the effective inflow layer, when one exists
def _shade_eil(ax, snd_params, theme, alpha=0.2):
    eil_z = snd_params["kinem"].get("eil_z", (ma.masked, ma.masked))
    if _has_eil(snd_params["kinem"]):
        ax.axhspan(eil_z[0], eil_z[1], color=theme["eil"], alpha=alpha)


# plot a profile using the 0-1 km and 1-3 km hodograph colors, with a fill beneath it
def _plot_low_level_trace(ax, values, z, theme, linewidth=3, alpha=1, fill_alpha=0.2):

    idx_1km = find_nearest(z, 1000)
    idx_3km = find_nearest(z, 3000)

    ax.plot(values[:idx_1km + 1], z[:idx_1km + 1], color=theme["hodo_layers"][0][2], linewidth=linewidth, alpha=alpha)
    ax.plot(values[idx_1km:idx_3km + 1], z[idx_1km:idx_3km + 1], color=theme["hodo_layers"][1][2], linewidth=linewidth, alpha=alpha)
    ax.fill_betweenx(z[:idx_3km + 1], values[:idx_3km + 1], color=theme["low_fill"], linewidth=0, alpha=fill_alpha)



def _plot_streamwiseness_panel(ax, snd_params, theme,
                               height_labels=(500, 1000, 1500, 2000), label_x=47):

    kinem = snd_params["kinem"]
    _setup_height_panel(ax, theme, (40, 102), [50, 70, 90], [50, 70, 90])
    _add_panel_title(ax, "Streamwiseness\n of ω (%)", theme, fontsize=FONT_SIZES["small"])

    for height in height_labels:
        ax.text(label_x, height + 2, f"{height / 1000:g} km".lstrip("0"), fontsize=FONT_SIZES["small"], weight="bold", color=theme["text"])

    if ma.is_masked(kinem["sm_u"]):
        warnings.warn("Streamwiseness could not be plotted (no valid storm motion/not enough data)", Warning)
        return ax

    _plot_low_level_trace(ax, kinem["swv_perc"], snd_params["intrp"]["zINTRP"], theme)
    _shade_eil(ax, snd_params, theme)

    return ax



def _plot_vorticity_panel(ax, snd_params, theme):

    streamwise_color = theme["streamwise"]
    kinem = snd_params["kinem"]
    _setup_height_panel(ax, theme, (0, 0.06), [0.01, 0.03, 0.05], [".01", ".03", ".05"])
    _add_panel_title(ax, "Total ω &\n Streamwise ω\n(/sec)", theme)

    if ma.is_masked(kinem["sm_u"]):
        warnings.warn("Total Vorticity could not be plotted (no valid storm motion/not enough data)", Warning)
        return ax

    z = snd_params["intrp"]["zINTRP"]
    idx_3km = find_nearest(z, 3000) + 1
    z = z[:idx_3km]
    streamwise = kinem["swv"][:idx_3km]
    total = kinem["vort"][:idx_3km]

    ax.plot(streamwise, z, color=streamwise_color, linewidth=3, alpha=0.8, label="SW ζ")
    ax.plot(total, z, color=theme["text"], linewidth=4, alpha=0.4, label="Total ζ")
    ax.fill_betweenx(z, total, where=total > streamwise, color=theme["text"], linewidth=0, alpha=0.1)
    ax.fill_betweenx(z, streamwise, color=streamwise_color, linewidth=0, alpha=0.2)
    _shade_eil(ax, snd_params, theme)

    return ax



def _plot_sr_wind_panel(ax, snd_params, theme, label_x=35):

    kinem = snd_params["kinem"]
    _setup_height_panel(ax, theme, (10, 50), [20, 30, 40], [20, 30, 40])
    _add_panel_title(ax, "Storm Relative\nWind (kts)", theme)

    # LCL and LFC labels, for data with thermodynamics
    if "thermo" in snd_params:
        thermo = snd_params["thermo"]
        ax.text(label_x, thermo["sb_lcl_z"], "-LCL-", fontsize=FONT_SIZES["small"], weight="bold", color=theme["text"], clip_on=True)
        if thermo["mu_lfc_z"] < 2500:
            ax.text(label_x, thermo["mu_lfc_z"], "-LFC-", fontsize=FONT_SIZES["small"], weight="bold", color=theme["text"], clip_on=True)

    if ma.is_masked(kinem["sm_u"]):
        warnings.warn("Storm Relative Wind could not be plotted (no valid storm motion/not enough data)", Warning)
        return ax

    _plot_low_level_trace(ax, kinem["srw"], snd_params["intrp"]["zINTRP"], theme, alpha=0.8)
    _shade_eil(ax, snd_params, theme)

    return ax



def _plot_theta_panel(ax, snd_params, theme):

    intrp = snd_params["intrp"]
    idx_3km = find_nearest(intrp["zINTRP"], 3000) + 1
    z = intrp["zINTRP"]
    theta = intrp["thetaINTRP"]
    theta_e = intrp["thetaeINTRP"]

    theta_min = np.nanmin(theta[:idx_3km])
    theta_max = np.nanmax(theta_e[:idx_3km])
    _setup_height_panel(ax, theme, (theta_min - 5, theta_max + 5), [theta_min, theta_max],
                        [int(theta_min), int(theta_max)], z_limits=(z[0], 3000))
    _add_panel_title(ax, "Theta-e &\nTheta (K)", theme)

    ax.plot(theta, z, color=theme["theta"], linewidth=3.5, alpha=0.5)
    ax.plot(theta_e, z, color=theme["theta"], linewidth=3.5, alpha=0.8)
    _shade_eil(ax, snd_params, theme)

    return ax



def _plot_wind_component_panel(ax, snd_params, theme):

    intrp = snd_params["intrp"]
    z = intrp["zINTRP"]
    idx_3km = find_nearest(z, 3000) + 1

    low = min(np.nanmin(intrp["uINTRP"][:idx_3km]), np.nanmin(intrp["vINTRP"][:idx_3km]))
    high = max(np.nanmax(intrp["uINTRP"][:idx_3km]), np.nanmax(intrp["vINTRP"][:idx_3km]))
    _setup_height_panel(ax, theme, (low - 10, high + 5), [low, high], [int(low), int(high)], z_limits=(z[0], 3000))
    _add_panel_title(ax, "U-Wind (green)\nV-Wind (blue)", theme)

    ax.plot(intrp["uINTRP"], z, color=theme["u_wind"], linewidth=3.5, alpha=0.5)
    ax.plot(intrp["vINTRP"], z, color=theme["v_wind"], linewidth=3.5, alpha=0.8)

    return ax



# CIN and CAPE share a height axis, but have separate x-axes
def _plot_cape_cin_panel(ax, snd_params, theme):

    intrp = snd_params["intrp"]
    z = intrp["zINTRP"]

    _setup_height_panel(ax, theme, (-300, 300), [-200, -100], [-200, -100],
                        z_limits=(z[0], 3000), tick_pad=-25, rotation=45)
    ax.axvline(0, linestyle=":", alpha=0.8, color=theme["zero_line"])
    ax.fill_betweenx(z, intrp["cin_profileINTRP"], color=theme["cin"], linewidth=0, alpha=0.5)
    _shade_eil(ax, snd_params, theme)

    cape_ax = ax.twiny()
    _set_panel_x_axis(cape_ax, theme, (-3000, 3000), [0, 1000, 2000], ["0 ", "1k", "2k"], tick_pad=-20, rotation=45)
    cape_ax.xaxis.set_ticks_position("bottom")
    _style_axis(cape_ax, theme)
    cape_ax.fill_betweenx(z, intrp["cape_profileINTRP"], color=theme["cape"], linewidth=0, alpha=0.8)
    # title sits on the top axis so the fills can't cover it
    _add_panel_title(cape_ax, "Stepwise\n CIN & CAPE\n (J/kg)", theme)

    return ax




#################################################################################
############################### TEXT TABLE UTILS ################################
#################################################################################

# set up a bordered box for text. positions inside it run from 0 to 1, so its contents move with it
def _setup_text_box(ax, theme):

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    _style_axis(ax, theme)
    for spine in ax.spines.values():
        spine.set_linewidth(1)

    return ax


def _draw_text(ax, x, y, text, color, fontsize=FONT_SIZES["table"], alpha=1, ha="left"):
    ax.text(x, y, text, color=color, fontsize=fontsize, alpha=alpha, ha=ha, weight="bold")


# shrink text that is wider than a fraction of an axis or figure, keeping it readable in its space
def _fit_text(text, target, max_width):
    if hasattr(target, "apply_aspect"):
        target.apply_aspect()
    width = text.get_window_extent(target.get_figure().canvas.get_renderer()).width / target.bbox.width
    if width > max_width:
        text.set_fontsize(text.get_fontsize() * max_width / width)


# format a value and its unit for a readout, showing MISSING_VALUE for missing values
def _format_value(value, unit="", decimals=None):
    if ma.is_masked(value) or np.isnan(getattr(value, "m", value)):
        return MISSING_VALUE
    text = mag(value) if decimals is None else round(float(getattr(value, "m", value)), decimals)
    return f"{text}{unit}"


# draw a table of header, row labels, and row values. each row sets its own color and alpha
def _draw_table(ax, table, layout, theme, fontsize=FONT_SIZES["table"]):

    for x, header in zip(layout["header_columns"], table["headers"]):
        _draw_text(ax, x, layout["header_y"], header, theme["text"], fontsize, layout.get("header_alpha", 1))

    for y, row in zip(layout["rows_y"], table["rows"]):
        _draw_text(ax, layout["label_x"], y, row["label"], theme["text"], fontsize, layout.get("label_alpha", 1))
        for x, value in zip(layout["columns"], row["values"]):
            _draw_text(ax, x, y, value, row["color"], fontsize, row.get("alpha", 1))



def _get_parcel_table(snd_params, theme):

    thermo = snd_params["thermo"]
    parcels = (("SB:", "sb"), ("ML:", "ml"), ("MU:", "mu"))

    rows = [{"label": label, "color": theme["parcel"][key], "values": [
                _format_value(thermo[f'{key}_ecape'], f" {UNIT_JKG}"),
                _format_value(thermo[f'{key}cape'], f" {UNIT_JKG}"),
                _format_value(thermo[f'{key}6cape'], f" {UNIT_JKG}"),
                _format_value(thermo[f'{key}3cape'], f" {UNIT_JKG}"),
                _format_value(thermo[f'{key}cin'], f" {UNIT_JKG}"),
                _format_value(thermo[f'{key}_lcl_z'], f" {UNIT_M}")]}
            for label, key in parcels]

    return {"headers": ["SR-ECAPE", "CAPE", "6CAPE", "3CAPE", "CIN", "LCL"], "rows": rows}



def _get_moisture_table(snd_params, theme):

    general = snd_params["general"]
    layers = (("0-.5ₖₘ:", "0_500", 1), ("0-1ₖₘ:", "0_1000", 0.9), ("1-3ₖₘ:", "1_3000", 0.8),
              ("3-6ₖₘ:", "3_6000", 0.7), ("6-9ₖₘ:", "6_9000", 0.64))

    rows = [{"label": label, "color": theme["moisture"], "alpha": alpha, "values": [
                _format_value(general[f'rh_{key}'], "%"),
                _format_value(general[f'w_{key}'], UNIT_GKG, decimals=1)]}
            for label, key, alpha in layers]

    return {"headers": ["RelHum", "MixRto"], "rows": rows}



def _get_kinematics_table(snd_params, theme, include_eil=True, layer_colors=False):

    kinem = snd_params["kinem"]
    layers = (("0-.5ₖₘ:", "0_to_500"), ("0-1ₖₘ:", "0_to_1000"), ("1-3ₖₘ:", "1_to_3000"),
              ("3-6ₖₘ:", "3_to_6000"), ("6-9ₖₘ:", "6_to_9000"))
    if include_eil:
        layers += (("EIL:", "eil"),)

    # blues can deepen with height, with the EIL in the standard blue
    colors = (*theme["kinematics_layers"], theme["kinematics"]) if layer_colors else [theme["kinematics"]] * len(layers)

    rows = [{"label": label, "color": color, "values": [
                _format_value(kinem[f'shear_{key}'], UNIT_KT),
                _format_value(kinem[f'srh_{key}'], UNIT_M2S2),
                _format_value(kinem[f'srw_{key}'], UNIT_KT),
                _format_value(kinem[f'swv_perc_{key}'], "%"),
                _format_value(kinem[f'swv_{key}'], UNIT_PER_S, decimals=3)]}
            for (label, key), color in zip(layers, colors)]

    return {"headers": ["BWD", "SRH", "SRW", "SWω%", "SWω"], "rows": rows}



def _plot_thermo_table(ax, snd_params, theme):

    thermo = snd_params["thermo"]
    text = theme["text"]

    _setup_text_box(ax, theme)
    _draw_table(ax, _get_parcel_table(snd_params, theme), THERMO_TABLE_LAYOUT, theme)
    _draw_text(ax, 0.015, THERMO_TABLE_LAYOUT["rows_y"][1] - 0.04, "100 hPa", text, fontsize=FONT_SIZES["fine"], alpha=0.5)

    # single-value thermo parameters, as a 2x4 grid:
    # downdraft | MUNCAPE & 700-500 hPa lapse rate | low & mid lapse rates | K-index & total totals
    top, bottom = 0.310, 0.175
    _draw_text(ax, 0.015, top, "DCAPE:", text)
    _draw_text(ax, 0.126, top, _format_value(thermo['dcape'], f" {UNIT_JKG}"), theme["cool"], alpha=0.7)
    _draw_text(ax, 0.015, bottom, "DCIN:", text)
    _draw_text(ax, 0.126, bottom, _format_value(thermo['dcin'], f" {UNIT_JKG}"), theme["cool"], alpha=0.7)

    _draw_text(ax, 0.300, top, "MUNCAPE:", text)
    _draw_text(ax, 0.450, top, _format_value(thermo['mu_ncape'], decimals=1), theme["ncape"], alpha=0.7)
    _draw_text(ax, 0.300, bottom, "Γ₇₀₀₋₅₀₀:", text, fontsize=FONT_SIZES["table_large"])
    _draw_text(ax, 0.435, bottom, _format_value(thermo['lr_75'], f" {UNIT_LAPSE}", decimals=1), theme["lapse"])

    _draw_text(ax, 0.555, top, "Γ₀₋₃:", text, fontsize=FONT_SIZES["table_large"])
    _draw_text(ax, 0.635, top, _format_value(thermo['lr_03km'], f" {UNIT_LAPSE}", decimals=1), theme["lapse"])
    _draw_text(ax, 0.555, bottom, "Γ₃₋₆:", text, fontsize=FONT_SIZES["table_large"])
    _draw_text(ax, 0.635, bottom, _format_value(thermo['lr_36km'], f" {UNIT_LAPSE}", decimals=1), theme["lapse"])

    for y, label, key in ((top, "KI:", 'k_index'), (bottom, "TT:", 't_totals')):
        value_text = _format_value(thermo[key], decimals=0)
        _draw_text(ax, 0.800, y, label, text)
        _draw_text(ax, 0.865, y, value_text.removesuffix(".0"), theme["ncape"], alpha=0.7)

    # composite indices, set off below a divider. all are unitless, so six fit on one row
    kinem = snd_params["kinem"]
    ax.axhline(0.125, xmin=0.012, xmax=0.988, color=text, alpha=0.45, linewidth=1.2)
    composites = (("SCP", kinem['eil_scp'], 1, theme["composite_index"]),
                  ("STP", kinem['eil_stp'], 1, theme["composite_index"]),
                  ("EHI₀₋₁", kinem['ehi_0_to_1000'], 1, theme["ehi"]),
                  ("EHI₀₋₃", kinem['ehi_0_to_3000'], 1, theme["ehi"]),
                  ("SHIP", thermo['ship'], 1, theme["composite_index"]),
                  ("MLBRN", thermo['ml_brn'], 0, theme["composite_index"]))
    # each label/value pair is centered in one of six equal slots
    slot = 0.97 / len(composites)
    for i, (label, value, decimals, color) in enumerate(composites):
        x = 0.015 + (i + 0.5) * slot
        value_text = _format_value(value, decimals=decimals)
        _draw_text(ax, x + 0.004, 0.035, label, text, fontsize=FONT_SIZES["tick"], ha="right")
        _draw_text(ax, x + 0.014, 0.035, value_text.removesuffix(".0") if decimals == 0 else value_text, color,
                   fontsize=FONT_SIZES["tick"])

    return ax



def _plot_moisture_table(ax, snd_params, theme):

    general = snd_params["general"]
    text = theme["text"]
    table = _get_moisture_table(snd_params, theme)

    _setup_text_box(ax, theme)
    _draw_table(ax, table, MOISTURE_TABLE_LAYOUT, theme)

    # line above each layer label, fading with height
    for y, row in zip(MOISTURE_TABLE_LAYOUT["rows_y"], table["rows"]):
        _draw_text(ax, 0.073, y + 0.087, "___", text, alpha=row["alpha"])

    # single-level values
    _draw_text(ax, 0.067, 0.122, "PWAT:", text, fontsize=FONT_SIZES["label"])
    _draw_text(ax, 0.400, 0.122, _format_value(general['pwat'], "in", decimals=3), theme["moisture_strong"], fontsize=FONT_SIZES["label"], ha="center")
    _draw_text(ax, 0.633, 0.122, "WB:", text, fontsize=FONT_SIZES["label"], ha="center")
    _draw_text(ax, 0.833, 0.122, _format_value(general['wet_bulb'][0], UNIT_DEGC), theme["moisture_strong"], fontsize=FONT_SIZES["label"], ha="center")
    _draw_text(ax, 0.067, 0.035, "FRZ:", text, fontsize=FONT_SIZES["label"])
    _draw_text(ax, 0.400, 0.035, _format_value(general['frz_pt_z'], "m"), theme["cool"], fontsize=FONT_SIZES["label"], ha="center")
    _draw_text(ax, 0.633, 0.035, "WB0:", text, fontsize=FONT_SIZES["label"], ha="center")
    _draw_text(ax, 0.840, 0.035, _format_value(general['wb_frz_pt_z'], "m"), theme["cool"], fontsize=FONT_SIZES["label"], ha="center")

    return ax



def _plot_kinematics_table(ax, snd_params, theme):
    _setup_text_box(ax, theme)
    _draw_table(ax, _get_kinematics_table(snd_params, theme), KINEMATICS_TABLE_LAYOUT, theme)

    return ax



def _plot_indices_box(ax, snd_params, theme):

    kinem = snd_params["kinem"]
    text = theme["text"]
    _setup_text_box(ax, theme)

    def _draw_stat(x, y_label, y_value, label, value, color):
        _draw_text(ax, x, y_label, label, text, fontsize=FONT_SIZES["tick"], alpha=0.9, ha="center")
        _draw_text(ax, x, y_value, value, color, fontsize=FONT_SIZES["tick"], alpha=0.9, ha="center")

    # SRH and bulk shear (composite indices are in the thermo table)
    _draw_stat(0.30, 0.800, 0.581, "  0-3ₖₘ SRH,", _format_value(kinem['srh_0_to_3000'], UNIT_M2S2), theme["kinematics_strong"])
    _draw_stat(0.77, 0.800, 0.581, "BWD", _format_value(kinem['shear_0_to_3000'], UNIT_KT), theme["kinematics_strong"])
    _draw_stat(0.30, 0.361, 0.161, "  0-6ₖₘ SRH,", _format_value(kinem['srh_0_to_6000'], UNIT_M2S2), theme["kinematics_strong"])
    _draw_stat(0.77, 0.361, 0.161, "BWD", _format_value(kinem['shear_0_to_6000'], UNIT_KT), theme["kinematics_strong"])

    return ax




# the kinematics table on its own, for wind-only data
def _plot_wind_parameter_box(ax, snd_params, theme):

    _setup_text_box(ax, theme)
    _draw_table(ax, _get_kinematics_table(snd_params, theme, include_eil=False, layer_colors=True), HODOGRAPH_KINEMATICS_LAYOUT, theme)

    return ax


# the kinematics and parcel tables and key values, all in one box
def _plot_hodograph_parameter_box(ax, snd_params, theme):

    thermo = snd_params["thermo"]
    text = theme["text"]

    _setup_text_box(ax, theme)
    _draw_table(ax, _get_kinematics_table(snd_params, theme, layer_colors=True), HODOGRAPH_KINEMATICS_LAYOUT, theme)
    _draw_table(ax, _get_parcel_table(snd_params, theme), HODOGRAPH_PARCEL_LAYOUT, theme)

    _draw_text(ax, 0.155, 0.133, "MUNCAPE:", text)
    _draw_text(ax, 0.374, 0.133, " " + _format_value(thermo["mu_ncape"], decimals=1), theme["ncape"], alpha=0.7)
    _draw_text(ax, 0.549, 0.133, "DCAPE:", text)
    _draw_text(ax, 0.733, 0.133, _format_value(thermo["dcape"], f" {UNIT_JKG}"), theme["cool"], alpha=0.7)
    _draw_text(ax, 0.155, 0.067, "Γ₀₋₃:", text, fontsize=FONT_SIZES["table_large"])
    _draw_text(ax, 0.374, 0.067, _format_value(thermo["lr_03km"], f" {UNIT_LAPSE}", decimals=1), theme["lapse"])
    _draw_text(ax, 0.549, 0.067, "Γ₃₋₆:", text, fontsize=FONT_SIZES["table_large"])
    _draw_text(ax, 0.733, 0.067, _format_value(thermo["lr_36km"], f" {UNIT_LAPSE}", decimals=1), theme["lapse"])

    return ax




#################################################################################
############################### COMPOSITE UTILS #################################
#################################################################################

# the letter that labels a profile: A-Z, then AA, AB, ...
def _get_profile_letter(index):
    letter = ""
    index += 1
    while index > 0:
        index, remainder = divmod(index - 1, 26)
        letter = chr(ord("A") + remainder) + letter
    return letter


# a user's list of styles for each profile, or the default for every profile
def _get_profile_styles(styles, default, n_profiles):
    if styles is None or (isinstance(styles, str) and styles == "none"):
        return [default] * n_profiles
    return list(styles)


# a user's list of colors, or colors evenly spaced along a colormap
def _get_profile_colors(colors, cmap, n_profiles, dark_mode=False):
    if colors is not None and not (isinstance(colors, str) and colors == "none"):
        return list(colors)
    start, end = (0, 1)
    if cmap == "viridis":
        start, end = COMPOSITE_VIRIDIS_RANGES["dark" if dark_mode else "light"]
    colormap = plt.get_cmap(cmap)
    return [colormap(value) for value in np.linspace(start, end, n_profiles)]


# data sources fill in unknown site info with placeholders
def _is_known(value):
    return value not in (None, "", "none") and not str(value).startswith("no-")


# describe what a profile is, from any data source
def _get_profile_type(site):
    source = str(site.get("source", ""))
    if "RAOB" in source:
        return "RAOB"
    if "ACARS" in source:
        return "ACARS"
    if _is_known(site.get("model")):
        hour = site.get("fcst-hour")
        return f"{site['model']} {hour}".upper() if _is_known(hour) else str(site["model"]).upper()
    return source if _is_known(source) else "PROFILE"


# describe where a profile is, from any data source
def _get_profile_location(site, short=False):
    if short and _is_known(site.get("site-id")):
        return str(site["site-id"])
    names = [str(site[key]) for key in ("site-id", "site-name") if _is_known(site.get(key))]
    if names:
        return ", ".join(names)
    latlon = site.get("site-latlon")
    return f"{latlon[0]}, {latlon[1]}" if latlon is not None else MISSING_VALUE


# custom data may not have a time
def _has_profile_time(time_parts):
    return isinstance(time_parts, (list, tuple)) and len(time_parts) >= 4


def _format_profile_time(time_parts):
    if not _has_profile_time(time_parts):
        return MISSING_VALUE
    return f"{time_parts[3]}Z | {time_parts[1]}-{time_parts[2]}-{time_parts[0]}"


# sort key for a profile's valid time, so int (2024, 4, 26, 9) and zero-padded ('2024', '04', '26', '09') parts match
def _profile_time_key(time_parts):
    try:
        return (0, tuple(int(float(part)) for part in time_parts[:4]))
    except (TypeError, ValueError):
        return (1, tuple(str(part) for part in time_parts[:4]))


# describe a set of profiles for the figure's subtitles
def _get_composite_titles(profiles, style="full", sr_hodo=False):
    sites = [profile["site_info"] for profile in profiles]
    times = {}
    for site in sites:
        if _has_profile_time(site.get("valid-time")):
            times.setdefault(_profile_time_key(site["valid-time"]), site["valid-time"])
    times = [times[key] for key in sorted(times)]
    locations = sorted({_get_profile_location(site) for site in sites})

    subtitle = f"{len(profiles)} PROFILE{'S' if len(profiles) > 1 else ''}"
    if times:
        subtitle += f" | VALID: {_format_profile_time(times[0])}"
    if len(times) > 1:
        subtitle += f" - {_format_profile_time(times[-1])}"

    top_titles = {"full": "COMPOSITE VERTICAL PROFILE", "sounding": "COMPOSITE SOUNDING",
                  "hodograph": "COMPOSITE STORM RELATIVE HODOGRAPH" if sr_hodo else "COMPOSITE HODOGRAPH"}
    return {"top_title": top_titles[style],
            "left_title": subtitle,
            "right_title": locations[0] if len(locations) == 1 else f"{len(locations)} LOCATIONS"}


# a skew-t's pressure range from its top pressure (hPa). the bottom is always the same
def _get_pressure_limits(ptop):
    if ptop is None:
        return None
    if isinstance(ptop, bool) or not isinstance(ptop, (int, float, np.integer, np.floating)) or not 0 < ptop < SKEWT_BOTTOM_PRESSURE:
        raise ValueError(f"'ptop' must be a pressure in hPa between 0 and {SKEWT_BOTTOM_PRESSURE}, such as 700. Got: {ptop}")
    return (SKEWT_BOTTOM_PRESSURE, ptop)


# check a skew-t's pressure range, as (bottom, top) in hPa
def _validate_pressure_limits(ylim):
    if ylim is None:
        return None
    if not isinstance(ylim, (list, tuple)) or len(ylim) != 2 or not all(isinstance(value, (int, float, np.integer, np.floating)) and not isinstance(value, bool) for value in ylim) \
            or ylim[0] <= ylim[1] or ylim[1] <= 0:
        raise ValueError(f"'ylim' must be (bottom, top) in hPa with bottom above top, such as (1000, 700). Got: {ylim}")
    return tuple(ylim)


# pressure ticks that suit a skew-t's pressure range
def _get_pressure_ticks(ylim):
    bottom, top = ylim
    span = bottom - top
    step = 100 if span > 400 else 50 if span > 150 else 25 if span > 75 else 10
    return [tick for tick in range(step * int(top // step + 1), int(bottom) + 1, step)]


# find the temperature range that fits every profile between two pressures on a skew-t.
# isotherms lean to the right with height, so the range is measured along the lean
def _get_zoomed_temperature_limits(profiles, ylim, skew_per_axis=0.933, margin=0.04, right_margin=0.04, extra_traces=(), min_span=12):
    # extra_traces are (pressure in hPa, temperature in degC) arrays to keep in view, such as parcel traces
    bottom, top = ylim
    traces = [(profile["p"].to("hPa").m, profile[key].to("degC").m) for profile in profiles for key in ("T", "Td")]
    traces += [(np.ma.filled(np.asarray(p, dtype=float), np.nan), np.ma.filled(np.asarray(t, dtype=float), np.nan)) for p, t in extra_traces]

    temperatures, lean = [], []
    for p, values in traces:
        in_view = np.isfinite(p) & np.isfinite(values) & (p <= bottom) & (p >= top)
        temperatures.append(values[in_view])
        lean.append(skew_per_axis * np.log(bottom / p[in_view]) / np.log(bottom / top))
    temperatures, lean = np.concatenate(temperatures), np.concatenate(lean)

    # no data in view
    if temperatures.size == 0:
        return (-32, 52)

    for span in np.arange(min_span, 300, 2):
        position = temperatures / span + lean
        if position.max() - position.min() <= 1 - margin - right_margin:
            break
    center = (position.max() + position.min()) / 2
    left = (center - (1 + margin - right_margin) / 2) * span
    return (left, left + span)


# keep a skew-t's temperature range when it shows only part of the atmosphere
# temperature ticks that suit a skew-t's temperature range
def _get_temperature_ticks(temperature_limits, target=6):
    span = temperature_limits[1] - temperature_limits[0]
    step = next(step for step in (2, 5, 10, 20) if span / step <= target + 2)
    ticks = np.arange(step * np.ceil(temperature_limits[0] / step), temperature_limits[1] + 0.01, step)
    # leave the left edge to the pressure labels
    return ticks[(ticks - temperature_limits[0]) / span > 0.07]


def _fit_skewt_aspect(skew, temperature_limits, ylim, box_aspect=0.87):
    span = temperature_limits[1] - temperature_limits[0]
    skew.ax.set_aspect(box_aspect * span / np.log10(ylim[0] / ylim[1]), adjustable="datalim")


# create a skew-t for a pressure range, fitting its temperature axis and pressure ticks to the data in it
def _setup_zoomed_skewt(fig, profiles, theme, ylim, rect, right_margin=0.04, extra_traces=()):
    temperature_limits = _get_zoomed_temperature_limits(profiles, ylim, right_margin=right_margin, extra_traces=extra_traces)
    skew = _setup_skewt(fig, profiles[0]["T"], theme, rect=rect, pressure_limits=ylim,
                        temperature_limits=temperature_limits, pressure_ticks=_get_pressure_ticks(ylim))
    skew.ax.set_xticks(_get_temperature_ticks(temperature_limits))
    skew.ax.set_xlim(*temperature_limits)
    _fit_skewt_aspect(skew, temperature_limits, ylim)
    return skew


# hide text outside a skew-t's axis, which would otherwise stretch a saved figure
def _clip_skewt_text(skew):
    for text in skew.ax.texts:
        text.set_clip_on(True)


def _plot_composite_skewt_profiles(skew, profiles, colors, styles, widths, alphas, shade_between=True):
    for i, profile in enumerate(profiles):
        line = dict(color=colors[i], linestyle=styles[i], linewidth=widths[i], alpha=alphas[i])
        skew.plot(profile["p"], profile["Td"], **line)
        skew.plot(profile["p"], profile["T"], label=f"PROFILE {_get_profile_letter(i)}", **line)
        if shade_between:
            skew.shade_area(profile["p"], x1=profile["Td"], x2=profile["T"], color=colors[i], alpha=0.05, zorder=1)


# interpolate a profile's winds to regular heights (m AGL)
def _get_profile_winds(profile, step=100, max_height=9000):
    z, u, v = profile["z"].m, profile["u"].m, profile["v"].m
    valid = np.isfinite(z) & np.isfinite(u) & np.isfinite(v)
    z, u, v = z[valid] - z[valid][0], u[valid], v[valid]
    levels = np.arange(0, min(z.max(), max_height) + step, step)
    return {"z": levels, "u": np.interp(levels, z, u), "v": np.interp(levels, z, v)}


# the letters of profiles without a storm motion
def _get_missing_storm_motions(snd_params):
    return [_get_profile_letter(i) for i, params in enumerate(snd_params)
            if params is None or ma.is_masked(params["kinem"]["sm_u"]) or np.isnan(params["kinem"]["sm_u"])]


# shift each profile's winds by its own storm motion
def _make_winds_storm_relative(winds, snd_params):
    return [dict(profile_winds, u=profile_winds["u"] - params["kinem"]["sm_u"], v=profile_winds["v"] - params["kinem"]["sm_v"])
            for profile_winds, params in zip(winds, snd_params)]


# draw each profile's hodograph trace, fading with height
def _plot_composite_hodograph_traces(hodo, winds, colors, linewidth=5):
    for profile_winds, color in zip(winds, colors):
        z, u, v = profile_winds["z"], profile_winds["u"], profile_winds["v"]
        for bottom, top, alpha in COMPOSITE_HODOGRAPH_LAYERS:
            layer = (z >= bottom) & (z <= top)
            hodo.ax.plot(u[layer], v[layer], color=color, alpha=alpha, linewidth=linewidth, clip_on=True)


# label what each hodograph alpha means
def _plot_composite_hodograph_key(hodo, theme, x=0.02, y=0.02, spacing=0.085, fontsize=FONT_SIZES["table"]):
    for i, (bottom, top, alpha) in enumerate(COMPOSITE_HODOGRAPH_LAYERS):
        hodo.ax.text(x + i * spacing, y, f"{bottom // 1000}-{top // 1000}km", transform=hodo.ax.transAxes,
                     color=theme["text"], weight="bold", fontsize=fontsize, alpha=alpha, zorder=8)


def _get_composite_parameter(snd_params, group, key, unit, decimals):
    if snd_params is None:
        return MISSING_VALUE
    unit = unit.format(jkg=UNIT_JKG, m=UNIT_M, m2s2=UNIT_M2S2, kt=UNIT_KT)
    return _format_value(snd_params[group][key], unit, decimals=decimals)


# one row per profile, in that profile's color. snd_params is None without parameter columns
def _plot_composite_table(ax, profiles, snd_params, colors, theme, style="full", fontsize=FONT_SIZES["table"]):
    _setup_text_box(ax, theme)
    columns = COMPOSITE_TABLES[style]
    fontsize = columns.get("fontsize", fontsize)
    parameters = [next(parameter for parameter in COMPOSITE_PARAMETERS if parameter[0] == header) for header in columns["parameters"]]
    n_rows = len(profiles) + 1
    y_of = lambda row: 1 - (row + 0.7) / (n_rows + 0.2)

    # the other columns can use the parameter columns' space when there aren't any
    header_y = y_of(0)
    _draw_text(ax, columns["profile"], header_y, "PROFILE", theme["text"], fontsize)
    _draw_text(ax, columns["location"], header_y, "LOCATION", theme["text"], fontsize)
    _draw_text(ax, columns["valid"], header_y, "VALID", theme["text"], fontsize)
    if snd_params is not None:
        for x, parameter in zip(columns["columns"], parameters):
            _draw_text(ax, x, header_y, parameter[0], theme["text"], fontsize)

    # keep text within its column
    def _draw_fitted(x, y, text, color, right):
        label = ax.text(x, y, text, color=color, fontsize=fontsize, weight="bold")
        _fit_text(label, ax, right - x - 0.006)

    for i, profile in enumerate(profiles):
        y = y_of(i + 1)
        site = profile["site_info"]
        _draw_text(ax, columns["letter"], y, _get_profile_letter(i), colors[i], fontsize)
        _draw_fitted(columns["profile"], y, _get_profile_type(site), colors[i], columns["location"])
        _draw_fitted(columns["location"], y, _get_profile_location(site, columns.get("short_location", False)), colors[i], columns["valid"])
        _draw_fitted(columns["valid"], y, _format_profile_time(site.get("valid-time")), colors[i], columns["columns"][0])
        if snd_params is not None:
            for x, parameter in zip(columns["columns"], parameters):
                _draw_text(ax, x, y, _get_composite_parameter(snd_params[i], *parameter[1:]), colors[i], fontsize)

    return ax




#################################################################################
############################### FIGURE EXTRAS UTILS #############################
#################################################################################

# add an empty area that a tight save has to include, setting the whitespace around the figure
def _add_figure_margin(fig, bounds):
    fig.add_artist(Rectangle(bounds[:2], bounds[2], bounds[3], transform=fig.transFigure, fill=False, linewidth=0, clip_on=False))


def _plot_figure_titles(fig, titles, theme, left_x=0.1237, right_x=1.206, top_y=0.985, sub_y=0.959, top_fontsize=FONT_SIZES["title"], sub_fontsize=FONT_SIZES["subtitle"]):
    fig.text(left_x, top_y, titles["top_title"], weight="bold", ha="left", fontsize=top_fontsize, color=theme["text"])
    fig.text(left_x, sub_y, titles["left_title"], ha="left", fontsize=sub_fontsize, color=theme["text"])
    fig.text(right_x, sub_y, titles["right_title"].strip(), ha="right", fontsize=sub_fontsize, color=theme["text"])


def _plot_skewt_legend(skew, theme):
    return skew.ax.legend(loc="upper left", framealpha=0.5, labelcolor=theme["text"], facecolor=theme["background"])




# add the watermark text to a figure (in figure fractions) or an axis (in axis fractions, behind its data)
def _plot_credits(target, theme, x=0.168, title_y=-0.105, subtitle_y=-0.125, style="full", title_fontsize=FONT_SIZES["watermark"], subtitle_fontsize=FONT_SIZES["watermark_sub"], max_x=0.525):
    if not WATERMARK["show"]:
        return
    default_title, default_subtitle = WATERMARK_TEXT[style]
    title = default_title if WATERMARK["title"] is None else WATERMARK["title"]
    subtitle = default_subtitle if WATERMARK["subtitle"] is None else WATERMARK["subtitle"]

    text_style = dict(ha="left", color=theme["watermark"] if WATERMARK["color"] is None else WATERMARK["color"], weight="bold", alpha=WATERMARK["alpha"])
    if hasattr(target, "transAxes"):
        text_style.update(transform=target.transAxes, zorder=1)
    texts = [target.text(x, title_y, title, fontsize=title_fontsize, **text_style),
             target.text(x, subtitle_y, subtitle, fontsize=subtitle_fontsize, **text_style)]

    # shrink long custom text to stay within its space
    for text in texts:
        _fit_text(text, target, max_x - x)


# add the sounderpy logo
def _plot_logo(ax):
    ax.axis("off")
    if WATERMARK["show"] and WATERMARK["logo"] is not None:
        ax.imshow(Image.open(WATERMARK["logo"]))
        ax.axis("off")

    return ax

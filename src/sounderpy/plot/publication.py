"""
    SOUNDERPY PUBLICATION FIGURES

    Purpose of module:

    House tools for plotting clean, journal-style sounding figures. These
    figures are experimental, and are meant to be restyled: every style setting
    can be changed, and each figure returns its Matplotlib objects.

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""

# METPY OPERATIONS
import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units

# MATPLOTLIB OPERATIONS
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import NullFormatter, NullLocator
from matplotlib.transforms import blended_transform_factory

# NUMPY OPERATIONS
import numpy as np
import numpy.ma as ma

# SOUNDERPY MODULES
from ..calc import sounding_params
from .style import MISSING_VALUE, SKEWT_BOTTOM_PRESSURE
from .components import _format_profile_time, _get_pressure_limits, _get_profile_location, _get_profile_type, _get_profile_winds, _prepare_sounding_data


#########################################################################
############################# CONSTANTS #################################
#########################################################################

KTS_TO_MS = 0.514444

PUBLICATION_STYLE = {"font_family": "sans-serif", "fontsize": 8, "panel_label_fontsize": 10,
                     "colors": ("red", "#1a7fd4", "black", "darkorange", "purple", "forestgreen"),
                     "linewidth": 1.1, "parcel_linewidth": 0.9, "parcel_linestyle": (0, (4, 2)),
                     "grid_color": "#bdbdbd", "grid_linewidth": 0.6, "axis_linewidth": 0.8,
                     "barb_length": 4.5, "barb_linewidth": 0.5, "barb_spacing": 0.13,
                     "hodograph_linewidth": 1.1, "marker_size": 4, "height_label_weight": "bold", "height_label_style": "italic",
                     "storm_motion_color": "black", "storm_motion_linewidth": 1.2, "sr_vector_color": "limegreen", "sr_vector_linewidth": 1.2}

# parameter name: label, group of sounding_params, key, multiplier, unit, and decimals
PUBLICATION_PARAMETERS = {"mllcl": ("MLLCL", "thermo", "ml_lcl_z", 0.001, "km", 1),
                          "mlcape": ("MLCAPE", "thermo", "mlcape", 1, r"J kg$^{-1}$", 0),
                          "mlcin": ("MLCIN", "thermo", "mlcin", -1, r"J kg$^{-1}$", 0),
                          "sbcape": ("SBCAPE", "thermo", "sbcape", 1, r"J kg$^{-1}$", 0),
                          "mucape": ("MUCAPE", "thermo", "mucape", 1, r"J kg$^{-1}$", 0),
                          "muecape": ("MU SR-ECAPE", "thermo", "mu_ecape", 1, r"J kg$^{-1}$", 0),
                          "srh01": ("0–1 km SRH", "kinem", "srh_0_to_1000", 1, r"m$^2$ s$^{-2}$", 0),
                          "srh03": ("0–3 km SRH", "kinem", "srh_0_to_3000", 1, r"m$^2$ s$^{-2}$", 0),
                          "bwd06": ("0–6 km BWD", "kinem", "shear_0_to_6000", KTS_TO_MS, r"m s$^{-1}$", 0),
                          "stp": ("STP", "kinem", "eil_stp", 1, "", 1),
                          "scp": ("SCP", "kinem", "eil_scp", 1, "", 1)}



#########################################################################
########################## PUBLICATION UTILS ############################
#########################################################################

# skew-t with isobars and isotherms, and its axis labels outside of the plot
def _setup_publication_skewt(fig, subplot, style, ylim, temperature_limits=(-40, 40), rotation=45, show_adiabats=False):

    skew = SkewT(fig, rotation=rotation, subplot=subplot)
    bottom, top = ylim
    skew.ax.set_ylim(bottom, top)
    skew.ax.set_xlim(*temperature_limits)

    ticks = [tick for tick in (200, 300, 400, 500, 600, 700, 800, 900, 1000) if top <= tick <= bottom]
    skew.ax.yaxis.set_major_locator(plt.FixedLocator(ticks))
    skew.ax.yaxis.set_major_formatter(plt.FixedFormatter([str(tick) for tick in ticks]))
    skew.ax.yaxis.set_minor_locator(NullLocator())
    skew.ax.yaxis.set_minor_formatter(NullFormatter())
    skew.ax.set_xticks(np.arange(temperature_limits[0], temperature_limits[1] + 1, 10))

    skew.ax.grid(True, color=style["grid_color"], linewidth=style["grid_linewidth"])
    skew.ax.tick_params(labelsize=style["fontsize"], width=style["axis_linewidth"], length=2.5, pad=2)
    skew.ax.set_xlabel("temperature (°C)", fontsize=style["fontsize"], labelpad=1)
    skew.ax.set_ylabel("pressure (mb)", fontsize=style["fontsize"], labelpad=1)
    for spine in skew.ax.spines.values():
        spine.set_linewidth(style["axis_linewidth"])

    if show_adiabats:
        skew.plot_dry_adiabats(colors=style["grid_color"], linewidth=style["grid_linewidth"], alpha=0.6)
        skew.plot_moist_adiabats(colors=style["grid_color"], linewidth=style["grid_linewidth"], alpha=0.6)

    return skew


def _plot_publication_profile(skew, data, snd_params, color, style, parcel="ml"):

    skew.plot(data["p"], data["T"], color=color, linewidth=style["linewidth"])
    skew.plot(data["p"], data["Td"], color=color, linewidth=style["linewidth"])

    if parcel is not None and snd_params["thermo"][f"{parcel}cape"] > 0:
        thermo = snd_params["thermo"]
        skew.plot(thermo[f"{parcel}P_trace"], thermo[f"{parcel}T_trace"], color=color,
                  linewidth=style["parcel_linewidth"], linestyle=style["parcel_linestyle"])


# a column of wind barbs beside the skew-t, in the profile's color
def _plot_publication_barbs(skew, data, color, column, style, ylim, n_barbs=14, first_x=1.09):

    p, u, v = data["p"], data["u"], data["v"]
    interval = np.logspace(np.log10(ylim[1]) + 0.02, np.log10(min(ylim[0], p[0].m)), n_barbs) * units.hPa
    idx = [i for i in mpcalc.resample_nn_1d(p, interval) if ylim[1] <= p[i].m <= ylim[0]]
    skew.plot_barbs(pressure=p[idx], u=u[idx], v=v[idx], xloc=first_x + column * style["barb_spacing"], color=color,
                    length=style["barb_length"], linewidth=style["barb_linewidth"], sizes={"emptybarb": 0.1})


# an open hodograph axis with its two axes crossing at the origin
def _setup_publication_hodograph(ax, limits, wind_unit, style, tick_spacing=None):

    xlim, ylim = limits
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")

    for name in ("top", "right"):
        ax.spines[name].set_visible(False)
    for name in ("left", "bottom"):
        ax.spines[name].set_position(("data", 0))
        ax.spines[name].set_linewidth(style["axis_linewidth"])

    if tick_spacing is None:
        tick_spacing = 10 if max(xlim[1] - xlim[0], ylim[1] - ylim[0]) < 60 else 20
    xticks = [tick for tick in np.arange(tick_spacing * np.ceil(xlim[0] / tick_spacing), xlim[1], tick_spacing) if tick != 0]
    yticks = [tick for tick in np.arange(tick_spacing * np.ceil(ylim[0] / tick_spacing), ylim[1], tick_spacing) if tick != 0]
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([f"{tick:g} {wind_unit}" if tick == xticks[-1] else f"{tick:g}" for tick in xticks])
    ax.tick_params(labelsize=style["fontsize"], width=style["axis_linewidth"], length=2.5, pad=2)
    ax.patch.set_alpha(0)

    ax.text(1.0, 0, r"$\overline{u}$", transform=ax.get_yaxis_transform(), ha="left", va="bottom", fontsize=style["fontsize"] + 1)
    ax.text(0, 1.0, r"$\overline{v}$", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=style["fontsize"] + 1)


def _format_height_label(height):
    return "sfc" if height == 0 else f"{height / 1000:g} km"


def _plot_publication_hodograph_trace(ax, winds, scale, color, style, heights=(0, 1000, 3000, 6000, 12000)):

    u, v, z = winds["u"] * scale, winds["v"] * scale, winds["z"]
    ax.plot(u, v, color=color, linewidth=style["hodograph_linewidth"], clip_on=False, zorder=3)

    for height in heights:
        if height > z.max():
            continue
        idx = np.abs(z - height).argmin()
        ax.plot(u[idx], v[idx], marker="o", color=color, markersize=style["marker_size"], clip_on=False, zorder=4)
        ax.annotate(_format_height_label(height), (u[idx], v[idx]), xytext=(0, 5), textcoords="offset points", ha="center", color=color,
                    fontsize=style["fontsize"], weight=style["height_label_weight"], fontstyle=style["height_label_style"], annotation_clip=False, zorder=5,
                    path_effects=[patheffects.withStroke(linewidth=1.5, foreground="white")])


def _draw_publication_arrow(ax, start, end, color, linewidth):
    ax.annotate("", xy=end, xytext=start, annotation_clip=False, zorder=2,
                arrowprops=dict(arrowstyle="-|>", color=color, linewidth=linewidth, shrinkA=0, shrinkB=0, mutation_scale=7))


def _plot_publication_storm_motion(ax, snd_params, winds, scale, style, sr_vector_heights=()):

    kinem = snd_params["kinem"]
    if ma.is_masked(kinem["sm_u"]) or ma.is_masked(kinem["sm_v"]):
        return
    motion = (kinem["sm_u"] * scale, kinem["sm_v"] * scale)
    _draw_publication_arrow(ax, (0, 0), motion, style["storm_motion_color"], style["storm_motion_linewidth"])

    # storm-relative wind vectors
    for height in sr_vector_heights:
        if height <= winds["z"].max():
            idx = np.abs(winds["z"] - height).argmin()
            _draw_publication_arrow(ax, motion, (winds["u"][idx] * scale, winds["v"][idx] * scale), style["sr_vector_color"], style["sr_vector_linewidth"])


def _format_publication_parameter(snd_params, parameter):

    label, group, key, multiplier, unit, decimals = PUBLICATION_PARAMETERS[parameter] if isinstance(parameter, str) else parameter
    value = snd_params[group][key]
    value = getattr(value, "m", value)
    if ma.is_masked(value) or not np.isfinite(value):
        return f"{label} = {MISSING_VALUE}"
    return f"{label} = {value * multiplier:.{decimals}f} {unit}".rstrip()


# list each profile's parameters in its color, in side-by-side columns
def _plot_publication_parameters(ax, snd_params_list, parameters, colors, style):

    ax.set_axis_off()
    for column, (snd_params, color) in enumerate(zip(snd_params_list, colors)):
        text = "\n".join(_format_publication_parameter(snd_params, parameter) for parameter in parameters)
        ax.text(column / max(len(snd_params_list), 1), 1, text, transform=ax.transAxes, ha="left", va="top", color=color,
                fontsize=style["fontsize"], linespacing=1.35)


# hodograph limits (in the plotted unit) with room for labels, always including the origin
def _get_publication_hodograph_limits(all_winds, scale, snd_params_list, padding=0.12):

    u = np.concatenate([winds["u"] for winds in all_winds]) * scale
    v = np.concatenate([winds["v"] for winds in all_winds]) * scale
    for snd_params in snd_params_list:
        if not ma.is_masked(snd_params["kinem"]["sm_u"]):
            u = np.append(u, snd_params["kinem"]["sm_u"] * scale)
            v = np.append(v, snd_params["kinem"]["sm_v"] * scale)

    limits = []
    for values in (u, v):
        low, high = min(np.nanmin(values), 0), max(np.nanmax(values), 0)
        pad = padding * max(high - low, 10)
        limits.append((low - pad, high + pad))
    return limits



#########################################################################
############################## FIGURES ##################################
#########################################################################

def build_publication_sounding(profiles, colors=None, storm_motion='right_moving', parcel='ml', ptop=None, show_adiabats=False,
                               show_barbs=True, hodo_heights=(0, 1000, 3000, 6000, 12000), hodo_max_height=12000, wind_units='m/s',
                               hodo_limits=None, sr_vector_heights=(), storm_motion_profile=0,
                               parameters=('mllcl', 'mlcape', 'mlcin', 'srh01', 'srh03', 'stp'), parameter_profiles=None,
                               panel_labels=('(a)', '(b)'), figsize=(7.5, 3.8), dpi=300, style=None, max_levels=250,
                               save=False, filename='sounderpy_publication.pdf'):
    '''
    EXPERIMENTAL: Return a clean, journal-style figure of a skew-t and hodograph for one or more profiles

    :param profiles: one data dictionary, or a list of them
    :type profiles: dict or list of dicts, required
    :param colors: a color for each profile. Default is ``None``, which uses the style's colors
    :type colors: list of strings, optional
    :param storm_motion: the storm motion used for every profile's calculations. Default is ``right_moving``
    :type storm_motion: str or list of floats, optional
    :param parcel: the parcel trace to draw: ``'sb'``, ``'ml'``, ``'mu'``, or ``None``. Default is ``'ml'``
    :type parcel: str or None, optional
    :param ptop: the pressure (hPa) at the top of the skew-t. Default is ``None``, which shows up to 100 hPa
    :type ptop: int or float, optional
    :param show_adiabats: whether to draw dry and moist adiabats. Default is ``False``
    :type show_adiabats: bool, optional
    :param show_barbs: whether to draw wind barb columns beside the skew-t. Default is ``True``
    :type show_barbs: bool, optional
    :param hodo_heights: heights (m AGL) to mark and label on the hodograph. Default is ``(0, 1000, 3000, 6000, 12000)``
    :type hodo_heights: tuple of numbers, optional
    :param hodo_max_height: the height (m AGL) that the hodograph stops at. Default is ``12000``
    :type hodo_max_height: int or float, optional
    :param wind_units: the hodograph's units: ``'m/s'`` or ``'kt'``. Default is ``'m/s'``
    :type wind_units: str, optional
    :param hodo_limits: hodograph limits as ``((u_min, u_max), (v_min, v_max))``. Default is ``None``, which fits the data
    :type hodo_limits: tuple, optional
    :param sr_vector_heights: heights (m AGL) to draw storm-relative wind vectors to. Default is none
    :type sr_vector_heights: tuple of numbers, optional
    :param storm_motion_profile: which profile's storm motion to draw, or ``None`` for none. Default is ``0``
    :type storm_motion_profile: int or None, optional
    :param parameters: parameters to list: names from ``PUBLICATION_PARAMETERS``, or ``(label, group, key, multiplier, unit, decimals)`` tuples
    :type parameters: list, optional
    :param parameter_profiles: which profiles to list parameters for. Default is ``None``, which lists all of them
    :type parameter_profiles: list of ints, optional
    :param panel_labels: the labels of the two panels, or ``None`` for no labels. Default is ``('(a)', '(b)')``
    :type panel_labels: tuple of strings, optional
    :param figsize: the figure's size in inches. Default is ``(7.5, 3.8)``
    :type figsize: tuple of numbers, optional
    :param dpi: dots per inch resolution of the figure. Default is ``300``
    :type dpi: int, optional
    :param style: changes to ``PUBLICATION_STYLE``, such as ``{'fontsize': 9, 'linewidth': 1.5}``
    :type style: dict, optional
    :param max_levels: the most vertical levels used for each profile. Default is ``250``
    :type max_levels: int or None, optional
    :param save: whether to save the figure to ``filename``. Default is ``False``
    :type save: bool, optional
    :param filename: the file to save to. Its extension sets the format, such as ``.pdf``, ``.svg``, or ``.png``
    :type filename: str, optional
    :return: a dict of the figure's Matplotlib objects, for further changes: ``fig``, ``skew``, ``hodo_ax``, and ``parameter_ax``
    :rtype: dict
    '''

    if wind_units not in ('m/s', 'kt'):
        raise ValueError(f"'wind_units' must be 'm/s' or 'kt'. Got: {wind_units}")

    style = {**PUBLICATION_STYLE, **(style or {})}
    profiles = [profiles] if isinstance(profiles, dict) else list(profiles)
    colors = list(style["colors"] if colors is None else colors)
    if len(colors) < len(profiles):
        raise ValueError(f"{len(profiles)} profiles need at least that many colors. Got {len(colors)}")
    scale = KTS_TO_MS if wind_units == 'm/s' else 1
    wind_label = r"m s$^{-1}$" if wind_units == 'm/s' else "kt"

    ylim = _get_pressure_limits(ptop) or (SKEWT_BOTTOM_PRESSURE, 100)

    # data
    profiles = [_prepare_sounding_data(profile, None, max_levels) for profile in profiles]
    snd_params = [sounding_params(profile, storm_motion).calc() for profile in profiles]
    winds = [_get_profile_winds(profile, max_height=hodo_max_height) for profile in profiles]

    with plt.rc_context({"font.family": style["font_family"], "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}):

        # layout
        fig = plt.figure(figsize=figsize, dpi=dpi, facecolor="white")
        grid = GridSpec(1, 2, figure=fig, width_ratios=[1.15, 1], left=0.075, right=0.98, top=0.95, bottom=0.11, wspace=0.42)
        right = GridSpecFromSubplotSpec(2, 1, subplot_spec=grid[1], height_ratios=[1.7, 1], hspace=0.12)

        # skew-t
        skew = _setup_publication_skewt(fig, grid[0], style, ylim, show_adiabats=show_adiabats)
        for i, (profile, params) in enumerate(zip(profiles, snd_params)):
            _plot_publication_profile(skew, profile, params, colors[i], style, parcel=parcel)
            if show_barbs:
                _plot_publication_barbs(skew, profile, colors[i], i, style, ylim)

        # hodograph
        hodo_ax = fig.add_subplot(right[0], anchor="NW")
        limits = hodo_limits or _get_publication_hodograph_limits(winds, scale, snd_params)
        _setup_publication_hodograph(hodo_ax, limits, wind_label, style)
        for i in range(len(profiles)):
            _plot_publication_hodograph_trace(hodo_ax, winds[i], scale, colors[i], style, heights=hodo_heights)
        if storm_motion_profile is not None:
            _plot_publication_storm_motion(hodo_ax, snd_params[storm_motion_profile], winds[storm_motion_profile], scale, style, sr_vector_heights)

        # parameters
        parameter_ax = fig.add_subplot(right[1])
        listed = range(len(profiles)) if parameter_profiles is None else parameter_profiles
        _plot_publication_parameters(parameter_ax, [snd_params[i] for i in listed], parameters, [colors[i] for i in listed], style)

        # panel labels
        if panel_labels is not None:
            skew.ax.text(-0.11, 1.0, panel_labels[0], transform=skew.ax.transAxes, va="top", fontsize=style["panel_label_fontsize"])
            hodo_ax.text(-0.10, 1.0, panel_labels[1], transform=blended_transform_factory(hodo_ax.transAxes, skew.ax.transAxes), va="top", fontsize=style["panel_label_fontsize"])

        if save:
            fig.savefig(filename, bbox_inches='tight', pad_inches=0.03)

    return {"fig": fig, "skew": skew, "hodo_ax": hodo_ax, "parameter_ax": parameter_ax}



#########################################################################
########################### EDITORIAL FIGURE ############################
#########################################################################

EDITORIAL_THEMES = {"light": {"background": "white", "text": "#1f2430", "muted": "#7a8190", "faint": "#e7e9ee", "halo": "white",
                              "temperature": "#d1344b", "dewpoint": "#0f8f7e", "parcel": "#3a3f4b", "cape": "#f0883e", "cin": "#5b8def"},
                    "dark": {"background": "#14161c", "text": "#eef0f4", "muted": "#8d94a3", "faint": "#2a2e38", "halo": "#14161c",
                             "temperature": "#ff5d73", "dewpoint": "#35d0b6", "parcel": "#d6d9e0", "cape": "#ff9d57", "cin": "#7aa5ff"}}

# colors up through the atmosphere, shared by the hodograph and wind barbs
EDITORIAL_HEIGHT_COLORS = ("#5b2a86", "#c2257a", "#ef5b3c", "#f6a821", "#7fb069")

EDITORIAL_STATS = (("MLCAPE", "thermo", "mlcape", 1, r"J kg$^{-1}$", 0), ("SR-ECAPE", "thermo", "mu_ecape", 1, r"J kg$^{-1}$", 0),
                   ("MLCIN", "thermo", "mlcin", 1, r"J kg$^{-1}$", 0), ("ML LCL", "thermo", "ml_lcl_z", 1, "m", 0),
                   ("0–1 km SRH", "kinem", "srh_0_to_1000", 1, r"m$^2$ s$^{-2}$", 0), ("0–3 km SRH", "kinem", "srh_0_to_3000", 1, r"m$^2$ s$^{-2}$", 0),
                   ("0–6 km shear", "kinem", "shear_0_to_6000", 1, "kt", 0), ("STP", "kinem", "eil_stp", 1, "", 1))


# a station's ID stays in capitals, and its name doesn't
def _get_editorial_title(site):
    parts = _get_profile_location(site).split(", ", 1)
    if len(parts) == 2 and parts[0].isalnum() and parts[0].isupper():
        return f"{parts[0]}  ·  {parts[1].title()}"
    return ", ".join(parts)


def _smooth_trace(values, window):
    # light smoothing for display. the ends stay fixed
    if window <= 1 or len(values) <= window:
        return values
    pad = window // 2
    padded = np.r_[np.full(pad, values[0]), values, np.full(pad, values[-1])]
    smoothed = np.convolve(padded, np.ones(window) / window, mode="valid")
    smoothed[0], smoothed[-1] = values[0], values[-1]
    return smoothed


def _get_trace_arrays(pressure, temperature):
    pressure, temperature = np.ma.filled(pressure, np.nan), np.ma.filled(temperature, np.nan)
    valid = np.isfinite(pressure) & np.isfinite(temperature)
    return pressure[valid], temperature[valid]


def _setup_editorial_skewt(fig, rect, data, theme, ylim):

    skew = SkewT(fig, rotation=45, rect=rect)
    surface = data["T"][0].to("degC").m
    left = 10 * np.floor((surface - 62) / 10)
    skew.ax.set_ylim(*ylim)
    skew.ax.set_xlim(left, left + 80)
    skew.ax.set_facecolor(theme["background"])

    ticks = [tick for tick in (1000, 850, 700, 500, 400, 300, 200, 150) if ylim[1] <= tick <= ylim[0]]
    skew.ax.yaxis.set_major_locator(plt.FixedLocator(ticks))
    skew.ax.yaxis.set_major_formatter(plt.FixedFormatter([str(tick) for tick in ticks]))
    skew.ax.yaxis.set_minor_locator(NullLocator())
    skew.ax.yaxis.set_minor_formatter(NullFormatter())
    skew.ax.set_xticks(np.arange(left, left + 81, 20))
    skew.ax.set_xlabel("")
    skew.ax.set_ylabel("")

    skew.ax.grid(True, color=theme["faint"], linewidth=0.7)
    skew.ax.tick_params(length=0, labelsize=8, labelcolor=theme["muted"], pad=3)
    skew.ax.tick_params(axis="x", pad=18)
    for spine in skew.ax.spines.values():
        spine.set_visible(False)
    skew.plot_dry_adiabats(colors=theme["faint"], linewidth=0.5, alpha=0.8)

    return skew


# soft fills where a parcel is warmer (CAPE) or cooler (CIN) than the environment
def _plot_editorial_buoyancy(skew, data, snd_params, theme, parcel):

    thermo = snd_params["thermo"]
    parcel_p, parcel_T = _get_trace_arrays(thermo[f"{parcel}P_trace"], thermo[f"{parcel}T_trace"])
    if parcel_p.size == 0:
        return
    env_p, env_T = _get_trace_arrays(data["p"].m, np.asarray(getattr(snd_params["general"]["virt_temp"], "m", snd_params["general"]["virt_temp"])))
    env_on_parcel = np.interp(-np.log(parcel_p), -np.log(env_p), env_T)

    # nothing is shaded above the equilibrium level, where a parcel is no longer rising
    el_p = thermo[f"{parcel}_el_p"]
    top_p = parcel_p.min() if ma.is_masked(el_p) or not np.isfinite(el_p) else el_p
    below_el = parcel_p >= top_p

    skew.ax.fill_betweenx(parcel_p, env_on_parcel, parcel_T, where=below_el & (parcel_T >= env_on_parcel), color=theme["cape"], alpha=0.22, linewidth=0, interpolate=True, zorder=1.5)
    skew.ax.fill_betweenx(parcel_p, env_on_parcel, parcel_T, where=below_el & (parcel_T < env_on_parcel), color=theme["cin"], alpha=0.18, linewidth=0, interpolate=True, zorder=1.5)
    skew.plot(parcel_p, parcel_T, color=theme["parcel"], linewidth=0.9, linestyle=(0, (1, 1.8)), zorder=3)

    # a small key, instead of a legend
    for y, text, color, weight in ((0.995, f"{parcel.upper()} parcel", theme["muted"], "normal"), (0.97, "CAPE", theme["cape"], "bold"), (0.945, "CIN", theme["cin"], "bold")):
        skew.ax.text(0.0, y, text, transform=skew.ax.transAxes, va="top", fontsize=7.5, color=color, weight=weight, zorder=6)


def _plot_editorial_profiles(skew, data, theme):

    for key, color in (("Td", theme["dewpoint"]), ("T", theme["temperature"])):
        line = skew.plot(data["p"], data[key], color=color, linewidth=2.0, zorder=4, solid_capstyle="round", solid_joinstyle="round")[0]
        line.set_path_effects([patheffects.withStroke(linewidth=3.4, foreground=theme["halo"], alpha=0.8), patheffects.Normal()])

    # label each line at its surface end, instead of using a legend
    for key, color, ha, offset, name in (("Td", theme["dewpoint"], "right", -6, "dewpoint"), ("T", theme["temperature"], "left", 6, "temperature")):
        value = data[key][0].to("degC").m
        skew.ax.annotate(f"{value:.0f}°", (value, data["p"][0].m), xytext=(offset, -3), textcoords="offset points", ha=ha, va="top", color=color, fontsize=13, weight="bold", annotation_clip=False)
        skew.ax.annotate(name, (value, data["p"][0].m), xytext=(offset, -18), textcoords="offset points", ha=ha, va="top", color=color, fontsize=7.5, annotation_clip=False)


# small marks for a parcel's key levels, along the right of the skew-t
def _plot_editorial_levels(skew, snd_params, theme, parcel, ylim, min_gap=0.022):

    thermo = snd_params["thermo"]
    transform = skew.ax.get_yaxis_transform()
    levels = []
    for name in ("lcl", "lfc", "el"):
        pressure = thermo[f"{parcel}_{name}_p"]
        if not ma.is_masked(pressure) and np.isfinite(pressure) and ylim[1] < pressure < ylim[0]:
            levels.append((float(pressure), name.upper()))

    # spread the labels when levels sit close together
    label_log = None
    for pressure, name in sorted(levels, reverse=True):
        log_p = np.log10(pressure)
        label_log = log_p if label_log is None else min(log_p, label_log - min_gap)
        skew.ax.plot([0.955, 0.985], [pressure, pressure], color=theme["muted"], linewidth=0.9, transform=transform, clip_on=False, zorder=5, solid_capstyle="round")
        skew.ax.text(0.95, 10 ** label_log, name, transform=transform, ha="right", va="center", fontsize=7.5, color=theme["muted"], weight="bold",
                     path_effects=[patheffects.withStroke(linewidth=2, foreground=theme["halo"])], zorder=5)


def _get_height_color(height, colormap, max_height=12000):
    return colormap(min(height / max_height, 1))


def _plot_editorial_barbs(skew, data, snd_params, colormap, theme, ylim, n_barbs=16, xloc=1.05):

    p, u, v = data["p"], data["u"], data["v"]
    height = data["z"].m - data["z"].m[0]
    interval = np.logspace(np.log10(ylim[1]) + 0.03, np.log10(p[0].m), n_barbs) * units.hPa
    for i in mpcalc.resample_nn_1d(p, interval):
        if ylim[1] <= p[i].m <= ylim[0] and np.isfinite(u[i].m):
            skew.plot_barbs(pressure=p[[i]], u=u[[i]], v=v[[i]], xloc=xloc, color=_get_height_color(height[i], colormap),
                            length=5, linewidth=0.8, sizes={"emptybarb": 0.08, "spacing": 0.14, "height": 0.42})


def _plot_editorial_hodograph(ax, snd_params, colormap, theme, smooth=5, max_height=12000, heights=(1000, 3000, 6000, 9000, 12000)):

    intrp = snd_params["intrp"]
    z = np.asarray(intrp["zINTRP"])
    keep = np.isfinite(z) & np.isfinite(intrp["uINTRP"]) & np.isfinite(intrp["vINTRP"]) & (z <= max_height)
    z, u, v = z[keep], _smooth_trace(np.asarray(intrp["uINTRP"])[keep], smooth), _smooth_trace(np.asarray(intrp["vINTRP"])[keep], smooth)
    kinem = snd_params["kinem"]

    # limits: a square window around the data and the origin
    has_motion = not ma.is_masked(kinem["sm_u"])
    all_u, all_v = np.r_[u, 0, kinem["sm_u"] if has_motion else 0], np.r_[v, 0, kinem["sm_v"] if has_motion else 0]
    half = max(np.ptp(all_u), np.ptp(all_v), 40) / 2 * 1.16
    center_u, center_v = (all_u.max() + all_u.min()) / 2, (all_v.max() + all_v.min()) / 2
    ax.set_xlim(center_u - half, center_u + half)
    ax.set_ylim(center_v - half, center_v + half)
    ax.set_aspect("equal")
    ax.set_axis_off()

    # speed rings
    corner = np.hypot(half + abs(center_u), half + abs(center_v))
    step = 20 if half < 70 else 40
    for speed in np.arange(step, corner, step):
        ax.add_patch(Circle((0, 0), speed, fill=False, edgecolor=theme["faint"], linewidth=0.8, zorder=1))
        x = speed * np.cos(np.deg2rad(-45)); y = speed * np.sin(np.deg2rad(-45))
        if abs(x - center_u) < half * 0.94 and abs(y - center_v) < half * 0.94:
            ax.text(x, y, f"{speed:.0f}", ha="center", va="center", fontsize=6.5, color=theme["muted"], zorder=1.2,
                    bbox=dict(facecolor=theme["background"], edgecolor="none", pad=0.8))
    ax.plot([center_u - half, center_u + half], [0, 0], color=theme["faint"], linewidth=0.8, zorder=1)
    ax.plot([0, 0], [center_v - half, center_v + half], color=theme["faint"], linewidth=0.8, zorder=1)

    # 0-3 km SRH, softly
    if has_motion:
        idx_3km = np.abs(z - 3000).argmin() + 1
        ax.fill(np.r_[u[:idx_3km], kinem["sm_u"]], np.r_[v[:idx_3km], kinem["sm_v"]], color=theme["cape"], alpha=0.13, linewidth=0, zorder=1.5)

    # the trace, with color changing smoothly with height
    points = np.column_stack([u, v]).reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    ax.plot(u, v, color=theme["halo"], linewidth=5, solid_capstyle="round", zorder=2)
    trace = LineCollection(segments, colors=[_get_height_color(height, colormap) for height in z[:-1]], linewidths=3, capstyle="round", zorder=3)
    ax.add_collection(trace)

    # height marks
    for height in heights:
        if height <= z.max():
            idx = np.abs(z - height).argmin()
            ax.plot(u[idx], v[idx], marker="o", markersize=5.5, markerfacecolor=theme["background"], markeredgecolor=_get_height_color(height, colormap), markeredgewidth=1.5, zorder=4)
            ax.annotate(f"{height // 1000}", (u[idx], v[idx]), xytext=(0, 6.5), textcoords="offset points", ha="center", fontsize=8, weight="bold", color=_get_height_color(height, colormap), zorder=5,
                        path_effects=[patheffects.withStroke(linewidth=2, foreground=theme["halo"])])

    # storm motions
    if not ma.is_masked(kinem["sm_rm"]):
        for key, name, offset in (("sm_rm", "RM", (6, -4)), ("sm_lm", "LM", (6, 4))):
            ax.plot(kinem[key][0], kinem[key][1], marker="o", markersize=3.2, color=theme["text"], zorder=5)
            ax.annotate(name, (kinem[key][0], kinem[key][1]), xytext=offset, textcoords="offset points", fontsize=7.5, weight="bold", color=theme["text"], va="center", zorder=5,
                        path_effects=[patheffects.withStroke(linewidth=2, foreground=theme["halo"])])

    ax.text(0.0, 1.0, "HODOGRAPH", transform=ax.transAxes, fontsize=7.5, weight="bold", color=theme["muted"], va="top")
    ax.text(0.0, 0.965, "knots, labeled by height in km", transform=ax.transAxes, fontsize=7, color=theme["muted"], va="top")


def _plot_editorial_stats(fig, snd_params, theme, stats, left, bottom, width, height, columns=4):

    rows = int(np.ceil(len(stats) / columns))
    for i, stat in enumerate(stats):
        label, group, key, multiplier, unit, decimals = stat
        value = snd_params[group][key]
        value = getattr(value, "m", value)
        missing = ma.is_masked(value) or not np.isfinite(value)
        x = left + (i % columns) * width / columns
        y = bottom + height - (i // columns) * height / rows

        fig.text(x, y, label, fontsize=7.5, color=theme["muted"], weight="bold", va="top")
        number = fig.text(x, y - 0.024, MISSING_VALUE if missing else f"{value * multiplier:,.{decimals}f}", fontsize=17, color=theme["text"], weight="bold", va="top")
        if unit and not missing:
            fig.text(x, y - 0.064, unit, fontsize=7, color=theme["muted"], va="top")


def build_editorial_sounding(data, dark_mode=False, storm_motion='right_moving', parcel='ml', ptop=None, smooth=5, stats=EDITORIAL_STATS,
                             height_colors=EDITORIAL_HEIGHT_COLORS, theme=None, title=None, subtitle=None, show_credit=True,
                             figsize=(11, 6.6), dpi=200, max_levels=250, save=False, filename='sounderpy_editorial.png'):
    '''
    EXPERIMENTAL: Return an open, modern figure of a sounding with no boxes, soft buoyancy shading, and winds colored by height

    :param data: one data dictionary
    :type data: dict, required
    :param dark_mode: whether to use the dark colors. Default is ``False``
    :type dark_mode: bool, optional
    :param storm_motion: the storm motion used for calculations. Default is ``right_moving``
    :type storm_motion: str or list of floats, optional
    :param parcel: the parcel to shade and mark: ``'sb'``, ``'ml'``, or ``'mu'``. Default is ``'ml'``
    :type parcel: str, optional
    :param ptop: the pressure (hPa) at the top of the skew-t. Default is ``None``, which shows up to 100 hPa
    :type ptop: int or float, optional
    :param smooth: the window (in 100 m levels) used to smooth the hodograph for display, or ``1`` for none. Default is ``5``. Calculations always use unsmoothed data
    :type smooth: int, optional
    :param stats: the values to show: ``(label, group, key, multiplier, unit, decimals)`` tuples
    :type stats: tuple of tuples, optional
    :param height_colors: colors from the surface to 12 km, for the hodograph and wind barbs
    :type height_colors: tuple of strings, optional
    :param theme: changes to the figure's colors, such as ``{'temperature': 'black'}``
    :type theme: dict, optional
    :param title: the figure's title. Default is ``None``, which uses the data's location
    :type title: str, optional
    :param subtitle: the figure's subtitle. Default is ``None``, which uses the data's type and time
    :type subtitle: str, optional
    :param show_credit: whether to show a small SounderPy credit. Default is ``True``
    :type show_credit: bool, optional
    :param figsize: the figure's size in inches. Default is ``(11, 6.6)``
    :type figsize: tuple of numbers, optional
    :param dpi: dots per inch resolution of the figure. Default is ``200``
    :type dpi: int, optional
    :param max_levels: the most vertical levels used. Default is ``250``
    :type max_levels: int or None, optional
    :param save: whether to save the figure to ``filename``. Default is ``False``
    :type save: bool, optional
    :param filename: the file to save to. Its extension sets the format
    :type filename: str, optional
    :return: a dict of the figure's Matplotlib objects: ``fig``, ``skew``, and ``hodo_ax``
    :rtype: dict
    '''

    colors = {**EDITORIAL_THEMES["dark" if dark_mode else "light"], **(theme or {})}
    colormap = LinearSegmentedColormap.from_list("height", height_colors)
    ylim = _get_pressure_limits(ptop) or (SKEWT_BOTTOM_PRESSURE, 100)

    profile = _prepare_sounding_data(data, None, max_levels)
    snd_params = sounding_params(profile, storm_motion).calc()
    site = profile["site_info"]

    with plt.rc_context({"font.family": "sans-serif", "pdf.fonttype": 42, "svg.fonttype": "none"}):
        fig = plt.figure(figsize=figsize, dpi=dpi, facecolor=colors["background"])

        # heading
        fig.text(0.045, 0.955, _get_editorial_title(site) if title is None else title, fontsize=21, weight="bold", color=colors["text"], va="top")
        fig.text(0.045, 0.892, f"{_get_profile_type(site)}  ·  {_format_profile_time(site.get('valid-time'))}" if subtitle is None else subtitle, fontsize=10.5, color=colors["muted"], va="top")

        # skew-t, winds, and hodograph
        skew = _setup_editorial_skewt(fig, (0.055, 0.085, 0.47, 0.77), profile, colors, ylim)
        _plot_editorial_buoyancy(skew, profile, snd_params, colors, parcel)
        _plot_editorial_profiles(skew, profile, colors)
        _plot_editorial_levels(skew, snd_params, colors, parcel, ylim)
        _plot_editorial_barbs(skew, profile, snd_params, colormap, colors, ylim)

        hodo_ax = fig.add_axes((0.60, 0.325, 0.37, 0.61))
        hodo_ax.set_facecolor(colors["background"])
        _plot_editorial_hodograph(hodo_ax, snd_params, colormap, colors, smooth=smooth)

        _plot_editorial_stats(fig, snd_params, colors, stats, 0.61, 0.075, 0.37, 0.215)

        if show_credit:
            fig.text(0.98, 0.022, "SounderPy", fontsize=7, color=colors["muted"], ha="right", weight="bold")

        if save:
            fig.savefig(filename, facecolor=colors["background"])

    return {"fig": fig, "skew": skew, "hodo_ax": hodo_ax}

"""
    SOUNDERPY PUBLICATION FIGURES

    Purpose of module:

    Clean, journal-style figures of one or more profiles: a skew-t, and optionally a hodograph, a
    location map, and a table of parameters, sized for AMS and AGU journals. The settings at the top
    set every part of the look, and the figure and its panels are returned for annotating. Called by
    spy.build_publication_sounding() (experimental in 3.3.0).

    Docs:
        publication figures:  https://sounderpy.app/docs/reference/publicationfigures.html
        tutorial:             https://sounderpy.app/docs/tutorials/custom_figures.html
        map insets & radar:   https://sounderpy.app/docs/reference/plotsettings.html#map-radar

    References:
        Hunter, J. D., 2007: Matplotlib: A 2D graphics environment. Computing in Science & Engineering, 9
            (3), 90-95 (figures)
        May, R. M., S. C. Arms, P. Marsh, E. Bruning, J. R. Leeman, K. Goebbert, J. E. Thielen, Z. S.
            Bruick, and M. D. Camron, 2023: MetPy: A Python package for meteorological data. Unidata,
            https://doi.org/10.5065/D6WW7G29 (skew-T)
        Met Office, 2010-2015: Cartopy: a cartographic Python library with a Matplotlib interface. Exeter,
            Devon, https://cartopy.readthedocs.io (map insets)

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""

# METPY OPERATIONS
import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units

# MATPLOTLIB OPERATIONS
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.patches import Circle
from matplotlib.ticker import NullFormatter, NullLocator

# NUMPY OPERATIONS
import numpy as np
import numpy.ma as ma

# STANDARD PYTHON
import warnings

# FOR MAPPING
import cartopy.crs as ccrs
from cartopy.mpl.feature_artist import FeatureArtist
from PIL import Image

# SOUNDERPY MODULES
from ..calc import sounding_params
from ..utils import SounderPyWarning
from .style import MISSING_VALUE, SKEWT_BOTTOM_PRESSURE, WATERMARK, WATERMARK_TEXT
from .components import (_get_profile_winds, _get_zoomed_temperature_limits, _plot_map_mosaic_radar,
                         _plot_map_single_site_radar, _prepare_sounding_data, _setup_map_inset, _validate_pressure_limits)



#########################################################################
############################# CONSTANTS #################################
#########################################################################
# Everything about how a publication figure looks is set here. Change a value to restyle every figure, or pass
# build_publication_sounding(style={...}) to change it for one figure.

PUBLICATION_STYLE = {
    # text: AMS and most AGU journals ask for about 8 pt, GRL for 6 pt
    "font_family": "sans-serif", "fontsize": 8, "panel_label_fontsize": 10,
    # one color per profile, in order
    "colors": ("red", "#1a7fd4", "black", "darkorange", "purple", "forestgreen"),
    # skew-t
    "linewidth": 1.1, "parcel_linewidth": 0.9, "parcel_linestyle": (0, (4, 2)), "temperature_range": (-40, 40),
    "grid_color": "#bdbdbd", "grid_linewidth": 0.6, "axis_linewidth": 0.8,
    "barb_length": 4.5, "barb_linewidth": 0.5, "barb_spacing": 0.17,   # barb_spacing: between profiles' columns (in)
    # hodograph: height markers are see-through white circles edged in the profile's color, with "S" (the surface)
    # or the height in km inside (for GRL, use a height_label_fontsize of 6 and a height_marker_size of 8)
    "hodograph_linewidth": 1.1, "height_marker_size": 7, "height_label_fontsize": 5.5, "height_label_weight": "bold",
    "height_marker_fill": (1, 1, 1, 0.8), "surface_label": "S",
    "hodograph_rings": True, "ring_color": "#d4d4d4", "ring_linewidth": 0.5,
    "storm_motion_color": "black", "storm_motion_linewidth": 1.2, "sr_vector_color": "gray", "sr_vector_linewidth": 1.2,
    # map: its layers, bottom to top (zorder), are radar 4, county lines 5, state lines 6, and profile markers 10
    "map_line_color": "#555555", "map_marker_size": 30,
    "radar_colorbar": True, "radar_range": (10, 75),
    # watermark (only drawn when asked for)
    "watermark_fontsize": 6, "watermark_color": "#9a9a9a",
}

# figure widths (in), from the journals' author guides (2026). "single" and "double" fit both AMS and AGU:
#   AMS: 19 / 27 / 33 / 39 picas = 3.2 / 4.5 / 5.5 / 6.5 in (one column / two-thirds page / two columns / full width)
#   AGU (Wiley): one column up to 85 mm (3.35 in), two columns up to 170 mm (6.7 in)
PUBLICATION_WIDTHS = {"single": 3.2, "double": 6.5,
                      "ams_single": 3.2, "ams_two_thirds": 4.5, "ams_two_column": 5.5, "ams_full": 6.5,
                      "agu_single": 3.35, "agu_double": 6.7}

# the panels, most important first: they're labeled (a), (b), ... in this order
PUBLICATION_PANELS = ("skewt", "hodograph", "map", "parameters")

# the layout for each set of panels: each panel's space as (left, top, width, height) in inches from the figure's top
# left, at the "double" width (6.5 in). Other widths scale the whole layout. The skew-t, hodograph, and map fill their
# spaces whatever their shape, the skew-t's space includes its wind barbs, and the map's includes its radar colorbar
PUBLICATION_LAYOUTS = {
    ("skewt",): {"skewt": (0.45, 0.25, 5.9, 4.6)},
    ("skewt", "hodograph"): {"skewt": (0.45, 0.25, 3.35, 2.8), "hodograph": (4.35, 0.25, 2.05, 2.8)},
    ("skewt", "map"): {"skewt": (0.45, 0.25, 3.35, 2.8), "map": (4.35, 0.25, 2.05, 2.8)},
    ("skewt", "parameters"): {"skewt": (0.45, 0.25, 5.9, 3.7), "parameters": (0.45, 4.6, 5.9, 0.7)},
    ("skewt", "hodograph", "map"): {"skewt": (0.45, 0.25, 3.45, 3.9), "hodograph": (4.45, 0.25, 1.95, 1.9),
                                    "map": (4.45, 2.55, 1.95, 1.85)},
    ("skewt", "hodograph", "parameters"): {"skewt": (0.45, 0.25, 3.35, 2.8), "hodograph": (4.35, 0.25, 2.05, 2.8),
                                           "parameters": (0.45, 3.7, 5.9, 0.7)},
    ("skewt", "map", "parameters"): {"skewt": (0.45, 0.25, 3.35, 2.8), "map": (4.35, 0.25, 2.05, 2.8),
                                     "parameters": (0.45, 3.7, 5.9, 0.7)},
    ("skewt", "hodograph", "map", "parameters"): {"skewt": (0.45, 0.25, 3.45, 3.9), "hodograph": (4.45, 0.25, 1.95, 1.9),
                                                  "map": (4.45, 2.55, 1.95, 1.85), "parameters": (0.45, 4.85, 5.9, 0.7)},
}

# the room (in) for the radar colorbar under the map, and the smallest panel (in) before a figure gives a warning
RADAR_COLORBAR_ROOM = 0.5
PUBLICATION_MIN_PANEL = 1.4

KTS_TO_MS = 0.514444

# parameter name: label, sounding_params group, key, multiplier, unit, and decimals
PUBLICATION_PARAMETERS = {"mllcl": ("MLLCL", "thermo", "ml_lcl_z", 0.001, "km", 1),
                          "mlcape": ("MLCAPE", "thermo", "mlcape", 1, r"J kg$^{-1}$", 0),
                          "mlcin": ("MLCIN", "thermo", "mlcin", 1, r"J kg$^{-1}$", 0),
                          "sbcape": ("SBCAPE", "thermo", "sbcape", 1, r"J kg$^{-1}$", 0),
                          "mucape": ("MUCAPE", "thermo", "mucape", 1, r"J kg$^{-1}$", 0),
                          "muecape": ("MU SR-ECAPE", "thermo", "mu_ecape", 1, r"J kg$^{-1}$", 0),
                          "srh01": ("0–1 km SRH", "kinem", "srh_0_to_1000", 1, r"m$^2$ s$^{-2}$", 0),
                          "srh03": ("0–3 km SRH", "kinem", "srh_0_to_3000", 1, r"m$^2$ s$^{-2}$", 0),
                          "bwd06": ("0–6 km BWD", "kinem", "shear_0_to_6000", KTS_TO_MS, r"m s$^{-1}$", 0),
                          "stp": ("STP", "kinem", "eil_stp", 1, "", 1),
                          "scp": ("SCP", "kinem", "eil_scp", 1, "", 1)}

# the skew-t's full pressure range (hPa)
FULL_PRESSURE_RANGE = (SKEWT_BOTTOM_PRESSURE, 100)

# fonts kept as text in PDF, PS, and SVG files, as journals ask
FONT_TYPES = {"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"}



#########################################################################
############################ FIGURE UTILS ###############################
#########################################################################

class PublicationFigure(Figure):
    """A Matplotlib figure that keeps its fonts as text whenever it's saved, including after it's returned."""

    def savefig(self, *args, **kwargs):
        with plt.rc_context(FONT_TYPES):
            return super().savefig(*args, **kwargs)


def _check_publication_options(profiles, panels, width, colors, labels, parcel, wind_units, map_radar,
                               storm_motion_profile, parameter_profiles, parameters):
    """Raise a ValueError that says what's wrong with an option."""

    n = len(profiles)
    checks = [
        (set(panels) <= set(PUBLICATION_PANELS) and len(set(panels)) == len(panels) and "skewt" in panels,
         f"'panels' must name the skew-t and any of the others, each once: {PUBLICATION_PANELS}. Got: {panels!r}"),
        (isinstance(PUBLICATION_WIDTHS.get(width, width), (int, float)) and PUBLICATION_WIDTHS.get(width, width) > 0,
         f"'width' must be one of {list(PUBLICATION_WIDTHS)} or a width in inches. Got: {width!r}"),
        (len(colors) >= n, f"{n} profiles need at least that many colors. Got {len(colors)}"),
        (labels is None or len(labels) >= n, f"{n} profiles need at least that many labels. Got {len(labels or ())}"),
        (parcel in (None, "sb", "ml", "mu"), f"'parcel' must be 'sb', 'ml', 'mu', or None. Got: {parcel!r}"),
        (wind_units in ("m/s", "kt"), f"'wind_units' must be 'm/s' or 'kt'. Got: {wind_units}"),
        (map_radar in (None, "mosaic", "single"), f"'map_radar' must be 'mosaic', 'single', or None. Got: {map_radar!r}"),
        ("map" in panels or not map_radar, "'map_radar' needs the 'map' panel"),
        ("map" not in panels or all(profile["site_info"].get("site-latlon") is not None for profile in profiles),
         "the 'map' panel needs each profile's location in site_info['site-latlon']"),
        ("parameters" not in panels or parameters, "the 'parameters' panel needs at least one parameter"),
        (storm_motion_profile is None or storm_motion_profile in range(n),
         f"'storm_motion_profile' must be a profile's number, 0 to {n - 1}, or None. Got: {storm_motion_profile!r}"),
        (parameter_profiles is None or all(number in range(n) for number in parameter_profiles),
         f"'parameter_profiles' must list profile numbers, 0 to {n - 1}. Got: {parameter_profiles}"),
    ]
    for ok, message in checks:
        if not ok:
            raise ValueError(message)
    for parameter in parameters:
        _get_parameter_spec(parameter)


def _get_panel_rects(panels, width, n_profiles, style, colorbar):
    """Each panel's axes as (left, top, width, height) in inches from the figure's top left, from its space in
    PUBLICATION_LAYOUTS scaled to the figure's width: the skew-t leaves room on its right for each profile's wind
    barbs, and the map leaves room under it for a radar colorbar."""

    scale = width / PUBLICATION_WIDTHS["double"]
    rects = {panel: [value * scale for value in space] for panel, space in PUBLICATION_LAYOUTS[tuple(panels)].items()}
    rects["skewt"][2] -= 0.12 + n_profiles * style["barb_spacing"]
    if "map" in rects and colorbar:
        rects["map"][3] -= RADAR_COLORBAR_ROOM
    return rects


def _position_panels(fig, axes, rects, shift, scale):
    """Put each panel's axes at its place (`rects`, in inches from the figure's top left), moved by `shift` (in,
    right and down) and scaled by `scale`."""

    width, height = fig.get_size_inches()
    for panel, (left, top, w, h) in rects.items():
        ax = axes[panel].ax if panel == "skewt" else axes[panel]
        x, y = shift[0] + left * scale, shift[1] + top * scale
        ax.set_position([x / width, 1 - (y + h * scale) / height, w * scale / width, h * scale / height])


def _fit_figure(fig, axes, rects, legend=None, pad=0.03, tries=8):
    """Fit the figure to everything drawn, with `pad` inches to spare on each side.

    The panels keep their places, and the figure's height changes to fit them; if anything reaches past the figure's
    sides, the panels shrink until it doesn't, so the figure stays the width asked for. `legend` is (skew, style) when
    there's a legend, which is placed again under the panels each time; a watermark is placed again under everything.
    """

    width = fig.get_figwidth()
    renderer = fig.canvas.get_renderer()
    panels_width = max(left + w for left, top, w, h in rects.values()) - min(left for left, top, w, h in rects.values())
    shift, scale = [0.0, 0.0], 1.0
    for _ in range(tries):
        _position_panels(fig, axes, rects, shift, scale)
        if legend:
            _add_publication_legend(fig, *legend)
        _place_publication_watermark(fig)
        drawn = fig.get_tightbbox(renderer)   # inches, from the bottom left
        height = fig.get_figheight()
        overhang = drawn.width - panels_width * scale   # labels and barbs past the panels' sides
        new_scale = scale * min(1.0, (width - 2 * pad - overhang) / (panels_width * scale))
        new_shift = [shift[0] + pad - drawn.x0, shift[1] + pad - (height - drawn.y1)]
        new_height = drawn.height + 2 * pad
        settled = abs(new_scale - scale) < 1e-3 and abs(new_height - height) < 0.005
        if settled and max(map(abs, np.subtract(new_shift, shift))) < 0.005:
            break
        shift, scale = new_shift, new_scale
        fig.set_size_inches(width, new_height, forward=False)


def _add_panel_labels(axes, panel_labels, style):
    """Label the panels "(a)", "(b)", ... just above their top left corners, in order."""
    size = style["panel_label_fontsize"]
    for (panel, ax), label in zip(axes.items(), panel_labels or ()):
        ax = ax.ax if panel == "skewt" else ax
        ax.annotate(label, xy=(0, 1), xycoords="axes fraction", xytext=(-0.4 * size, 6), textcoords="offset points",
                    ha="right", va="bottom", fontsize=size, gid="panel label")


def _add_publication_legend(fig, skew, style):
    """The profiles' names in one row centered under the panels, across the bottom of the figure (more rows only if
    one row is wider than the figure). Placed again each time the figure is fitted."""

    renderer = fig.canvas.get_renderer()
    for legend in list(fig.legends):
        legend.remove()
    panels = [ax for ax in fig.axes if ax.get_label() != "watermark logo"]
    lowest = min(ax.get_tightbbox(renderer).y0 for ax in panels) / fig.bbox.height
    top = lowest - style["fontsize"] * 0.6 / 72 / fig.get_figheight()
    handles, labels = skew.ax.get_legend_handles_labels()
    for columns in range(len(labels), 0, -1):
        legend = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, top), ncol=columns, frameon=False,
                            fontsize=style["fontsize"], borderaxespad=0, handlelength=1.5, columnspacing=1.5)
        if columns == 1 or legend.get_window_extent(renderer).width <= fig.bbox.width * 0.98:
            break
        legend.remove()


def _add_publication_watermark(fig, settings, style):
    """Small credit text (and a logo, if one is given) in a row of its own across the bottom of the figure, under the
    panels and the legend: the text at the right and the logo at the left.

    `settings` uses set_watermark()'s names; a title, subtitle, color, or alpha left out comes from the session's
    set_watermark(). The logo is drawn only when given in `settings`.
    """

    settings = {} if settings is True else dict(settings)
    title, subtitle = WATERMARK_TEXT["full"]
    title = settings.get("title", WATERMARK["title"] or title)
    subtitle = settings.get("subtitle", WATERMARK["subtitle"] or subtitle)
    text = "\n".join(line for line in (title, subtitle) if line)
    fig.text(0.995, 0.005, text, ha="right", va="top", fontsize=style["watermark_fontsize"], linespacing=1.2,
             color=settings.get("color", WATERMARK["color"] or style["watermark_color"]),
             alpha=settings.get("alpha", WATERMARK["alpha"]), gid="watermark")
    if settings.get("logo"):
        logo_ax = fig.add_axes((0.005, 0.005, 0.06, 0.06 * fig.get_figwidth() / fig.get_figheight()), label="watermark logo")
        logo_ax.imshow(Image.open(settings["logo"]))
        logo_ax.axis("off")
    _place_publication_watermark(fig)


def _place_publication_watermark(fig):
    """Move the watermark just under the lowest panel or legend. Placed again each time the figure is fitted."""

    text = next((t for t in fig.texts if t.get_gid() == "watermark"), None)
    logo = next((ax for ax in fig.axes if ax.get_label() == "watermark logo"), None)
    if text is None and logo is None:
        return
    renderer = fig.canvas.get_renderer()
    drawn = [ax.get_tightbbox(renderer) for ax in fig.axes if ax is not logo]
    drawn += [legend.get_window_extent(renderer) for legend in fig.legends]
    gap = 4 / 72 * fig.dpi   # 4 points
    top = (min(box.y0 for box in drawn if box is not None) - gap) / fig.bbox.height
    if text is not None:
        text.set_y(top)
    if logo is not None:
        size = 0.06 * fig.get_figwidth()   # inches, square
        logo.set_position((0.005, top - size / fig.get_figheight(), size / fig.get_figwidth(), size / fig.get_figheight()))



#########################################################################
############################# SKEW-T UTILS ##############################
#########################################################################

def _setup_publication_skewt(fig, style, ylim, temperature_range, show_adiabats=False):
    """Create an open skew-t with light isobars and isotherms, and its tick and axis labels outside of it."""

    # placed in its space once the figure is laid out, and filling it whatever its shape (the isotherms stay at 45°)
    skew = SkewT(fig, rotation=45, rect=(0, 0, 1, 1), aspect="auto")
    skew.ax.set_ylim(*ylim)
    skew.ax.set_xlim(*temperature_range)

    ticks = [tick for tick in range(100, 1100, 100) if ylim[1] <= tick <= ylim[0]]
    skew.ax.yaxis.set_major_locator(plt.FixedLocator(ticks))
    skew.ax.yaxis.set_major_formatter(plt.FixedFormatter([str(tick) for tick in ticks]))
    skew.ax.yaxis.set_minor_locator(NullLocator())
    skew.ax.yaxis.set_minor_formatter(NullFormatter())
    _set_temperature_ticks(skew, style)

    # open on the top and right
    for name in ("top", "right"):
        skew.ax.spines[name].set_visible(False)
    for spine in skew.ax.spines.values():
        spine.set_linewidth(style["axis_linewidth"])
    skew.ax.grid(True, color=style["grid_color"], linewidth=style["grid_linewidth"])
    skew.ax.tick_params(labelsize=style["fontsize"], width=style["axis_linewidth"], length=2.5, pad=2)
    skew.ax.set_xlabel("temperature (°C)", fontsize=style["fontsize"], labelpad=1)
    skew.ax.set_ylabel("pressure (hPa)", fontsize=style["fontsize"], labelpad=1)

    if show_adiabats:
        skew.plot_dry_adiabats(colors=style["grid_color"], linewidth=style["grid_linewidth"], alpha=0.6)
        skew.plot_moist_adiabats(colors=style["grid_color"], linewidth=style["grid_linewidth"], alpha=0.6)
    return skew


def _set_temperature_ticks(skew, style):
    """Temperature ticks every 10 °C (5 °C on a narrow range), or every 20 (10) °C when the skew-t is too narrow
    for their labels. Returns whether the ticks changed."""

    ax = skew.ax
    axis_width = ax.get_window_extent(ax.figure.canvas.get_renderer()).width * 72 / ax.figure.dpi   # points
    low, high = ax.get_xlim()
    step = 10 if high - low > 40 else 5
    while True:
        ticks = np.arange(step * np.ceil(low / step), high + 0.01, step)
        # each label is about 0.6 of the font size wide per character, with half the font size between labels
        needed = sum((len(f"{tick:g}") * 0.6 + 0.5) * style["fontsize"] for tick in ticks)
        if needed <= axis_width or len(ticks) <= 3:
            break
        step *= 2
    changed = not np.array_equal(ticks, ax.get_xticks())
    ax.set_xticks(ticks)
    return changed


def _plot_publication_profile(skew, data, snd_params, color, style, parcel="ml", label=None):
    """Plot one profile's temperature and dewpoint, and its parcel trace when the parcel has CAPE."""
    skew.plot(data["p"], data["T"], color=color, linewidth=style["linewidth"], label=label)
    skew.plot(data["p"], data["Td"], color=color, linewidth=style["linewidth"])

    thermo = snd_params["thermo"]
    if parcel is not None and thermo[f"{parcel}cape"] > 0:
        skew.plot(thermo[f"{parcel}P_trace"], thermo[f"{parcel}T_trace"], color=color,
                  linewidth=style["parcel_linewidth"], linestyle=style["parcel_linestyle"])


def _plot_publication_barbs(skew, data, color, column, style, ylim, size, n_barbs=14):
    """Plot a column of wind barbs beside the skew-t; each profile gets the next column over, barb_spacing apart.
    `size` is the skew-t's (width, height) in inches, to place the columns in inches."""

    p, u, v = data["p"], data["u"], data["v"]
    axis_width, axis_height = size
    # the lowest barb sits a barb's length above the bottom of the axis, so a barb pointing down stays inside it
    room = min(style["barb_length"] * 2.2 / (axis_height * 72), 0.5)
    lowest = 10 ** (np.log10(ylim[0]) - room * np.log10(ylim[0] / ylim[1]))
    levels = np.logspace(np.log10(ylim[1]) + 0.02, np.log10(min(lowest, p[0].m)), n_barbs) * units.hPa
    idx = [i for i in mpcalc.resample_nn_1d(p, levels) if ylim[1] <= p[i].m <= ylim[0]]
    xloc = 1 + (0.06 + (column + 0.5) * style["barb_spacing"]) / axis_width
    skew.plot_barbs(pressure=p[idx], u=u[idx], v=v[idx], xloc=xloc, color=color,
                    length=style["barb_length"], linewidth=style["barb_linewidth"], sizes={"emptybarb": 0.1})



#########################################################################
########################### HODOGRAPH UTILS #############################
#########################################################################

def _get_publication_hodograph_limits(winds, scale, snd_params, shape, padding=0.12):
    """Hodograph limits that fit every profile and storm motion, always including the origin, in the panel's
    `shape` (height / width)."""

    u = [profile["u"] * scale for profile in winds]
    v = [profile["v"] * scale for profile in winds]
    for params in snd_params:
        if not ma.is_masked(params["kinem"]["sm_u"]):
            u.append([params["kinem"]["sm_u"] * scale])
            v.append([params["kinem"]["sm_v"] * scale])

    limits = []
    for values in (np.concatenate(u), np.concatenate(v)):
        low, high = min(np.nanmin(values), 0), max(np.nanmax(values), 0)
        pad = padding * max(high - low, 10)
        limits.append((low - pad, high + pad))
    return _fit_limits(limits, shape)


def _fit_limits(limits, shape):
    """Widen the u or v limits around their middle, so the hodograph (one unit the same length both ways) fills a
    panel of `shape` (height / width)."""
    (u_low, u_high), (v_low, v_high) = limits
    u_span, v_span = max(u_high - u_low, (v_high - v_low) / shape), max(v_high - v_low, (u_high - u_low) * shape)
    return [((u_low + u_high - u_span) / 2, (u_low + u_high + u_span) / 2),
            ((v_low + v_high - v_span) / 2, (v_low + v_high + v_span) / 2)]


def _setup_publication_hodograph(ax, limits, wind_unit, style):
    """Set up an open hodograph, its axes crossing at the origin, with light range rings."""

    (xmin, xmax), (ymin, ymax) = limits
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_aspect("equal", adjustable="box")
    ax.patch.set_alpha(0)

    for name in ("top", "right"):
        ax.spines[name].set_visible(False)
    for name in ("left", "bottom"):
        ax.spines[name].set_position(("data", 0))
        ax.spines[name].set_linewidth(style["axis_linewidth"])

    spacing = 10 if min(xmax - xmin, ymax - ymin) < 60 else 20
    xticks = [tick for tick in np.arange(spacing * np.ceil(xmin / spacing), xmax, spacing) if tick != 0]
    ax.set_xticks(xticks, [f"{tick:g}" for tick in xticks])
    ax.set_yticks([tick for tick in np.arange(spacing * np.ceil(ymin / spacing), ymax, spacing) if tick != 0])
    ax.tick_params(labelsize=style["fontsize"], width=style["axis_linewidth"], length=2.5, pad=2)

    if style["hodograph_rings"]:
        corner = max(np.hypot(x, y) for x in (xmin, xmax) for y in (ymin, ymax))
        for radius in np.arange(spacing, corner + spacing, spacing):
            ax.add_patch(Circle((0, 0), radius, fill=False, edgecolor=style["ring_color"], linewidth=style["ring_linewidth"],
                                zorder=0))

    ax.text(1.0, 0, r"$\overline{u}$", transform=ax.get_yaxis_transform(), ha="left", va="bottom",
            fontsize=style["fontsize"] + 1)
    ax.text(0, 1.0, r"$\overline{v}$", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
            fontsize=style["fontsize"] + 1)
    # the unit beside the v label at the top, where it can't run into a neighboring panel
    ax.annotate(f"({wind_unit})", xy=(0, 1.0), xycoords=ax.get_xaxis_transform(), xytext=(6, 0), textcoords="offset points",
                ha="left", va="bottom", fontsize=style["fontsize"])


def _place_height_marker(ax, point, placed, style):
    """Where to draw a height marker so it doesn't cover one already drawn (`placed`, in display pixels).

    That's its own point when it's clear, otherwise the nearest clear spot inside the axis on rings around it.
    Returns the spot (in data units) and whether it moved.
    """

    size = style["height_marker_size"] * ax.figure.dpi / 72 * 1.1
    box = ax.bbox
    true = ax.transData.transform(point)

    def clear(spot):
        inside = box.x0 + size / 2 <= spot[0] <= box.x1 - size / 2 and box.y0 + size / 2 <= spot[1] <= box.y1 - size / 2
        return inside and all(np.hypot(*(spot - other)) >= size for other in placed)

    spot = true
    if any(np.hypot(*(true - other)) < size for other in placed):
        candidates = (true + radius * size * np.array([np.cos(angle), np.sin(angle)])
                      for radius in (1, 1.5, 2, 2.5, 3) for angle in np.radians(np.arange(0, 360, 30)))
        spot = next((candidate for candidate in candidates if clear(candidate)), true)
    placed.append(spot)
    return ax.transData.inverted().transform(spot), spot is not true


def _plot_publication_height_markers(ax, winds, scale, color, style, heights, placed):
    """Plot a labeled marker at each of `heights` (m AGL) on one profile's hodograph.

    "S" marks the surface and the rest show the height in km. A marker that would cover another (from any profile,
    shared through `placed`) moves aside, with a thin line to its point.
    """

    u, v, z = winds["u"] * scale, winds["v"] * scale, winds["z"]
    for height in heights:
        if height > z.max():
            continue
        idx = np.abs(z - height).argmin()
        (x, y), moved = _place_height_marker(ax, (u[idx], v[idx]), placed, style)
        if moved:
            ax.plot([u[idx], x], [v[idx], y], color=color, linewidth=0.6, clip_on=False, zorder=4)
        ax.plot(x, y, marker="o", markersize=style["height_marker_size"], markerfacecolor=style["height_marker_fill"],
                markeredgecolor=color, markeredgewidth=style["hodograph_linewidth"] * 0.8, clip_on=False, zorder=4,
                gid="height_marker")
        label = style["surface_label"] if height == 0 else f"{height / 1000:g}".lstrip("0")
        ax.text(x, y, label, ha="center", va="center_baseline", color=color, fontsize=style["height_label_fontsize"],
                weight=style["height_label_weight"], clip_on=False, zorder=5)


def _draw_publication_arrow(ax, start, end, color, linewidth, zorder=2, stop=0):
    """Draw an arrow from `start` to `end`, its head stopping `stop` points short of `end`."""
    ax.annotate("", xy=end, xytext=start, annotation_clip=False, zorder=zorder,
                arrowprops=dict(arrowstyle="-|>", color=color, linewidth=linewidth, shrinkA=0, shrinkB=stop, mutation_scale=7))


def _plot_publication_storm_motion(ax, snd_params, winds, scale, style, sr_vector_heights=(), marker_heights=()):
    """Draw the storm motion arrow, and storm-relative wind vectors to `sr_vector_heights` (m AGL).

    A storm-relative vector is drawn over the hodograph traces, from the storm motion to the wind at its height,
    and stops at the edge of a height marker there; one too short to show outside the marker is left out.
    """

    kinem = snd_params["kinem"]
    if ma.is_masked(kinem["sm_u"]) or ma.is_masked(kinem["sm_v"]):
        return
    motion = (kinem["sm_u"] * scale, kinem["sm_v"] * scale)
    _draw_publication_arrow(ax, (0, 0), motion, style["storm_motion_color"], style["storm_motion_linewidth"])

    for height in sr_vector_heights:
        if height > winds["z"].max():
            continue
        idx = np.abs(winds["z"] - height).argmin()
        end = (winds["u"][idx] * scale, winds["v"][idx] * scale)
        stop = style["height_marker_size"] / 2 + 0.5 if height in marker_heights else 0
        length = np.hypot(*(ax.transData.transform(end) - ax.transData.transform(motion))) * 72 / ax.figure.dpi   # points
        if not stop or length >= stop + 4:
            _draw_publication_arrow(ax, motion, end, style["sr_vector_color"], style["sr_vector_linewidth"], zorder=3.5, stop=stop)



#########################################################################
############################## MAP UTILS ################################
#########################################################################

def _map_colors(style):
    """The colors SounderPy's map helpers take (as a theme), from the publication style."""
    return {"border": style["map_line_color"], "background": "white", "text": style["map_line_color"]}


def _plot_publication_map(ax, profiles, colors, style, shape, zoom=1.5, show_counties=True):
    """Draw a clean map centered on the profiles, with a dot for each in its color, filling a panel of `shape`
    (height / width).

    Its shorter side shows +/- `zoom` degrees, grown to fit every profile, with no satellite image and no state fill.
    """

    latlons = [profile["site_info"]["site-latlon"] for profile in profiles]
    lats, lons = zip(*latlons)
    zoom = max(zoom, (max(lons) - min(lons)) / 2 + 0.5, (max(lats) - min(lats)) / 2 + 0.5)
    center = ((max(lats) + min(lats)) / 2, (max(lons) + min(lons)) / 2)

    # the map helper frames its point low, so draw it wider (loading the lines for the whole view), then center the
    # view in the panel's shape and take off the state fill (it's there to dim satellite images)
    across, up = zoom * max(1, 1 / shape), zoom * max(1, shape)
    _setup_map_inset(ax, center, _map_colors(style), map_zoom=max(across, up) * 1.6, show_satellite=False,
                     show_counties=show_counties)
    ax.set_box_aspect(shape)
    ax.set_extent([center[1] - across, center[1] + across, center[0] - up, center[0] + up], crs=ccrs.PlateCarree())
    for artist in list(ax.get_children()):
        if isinstance(artist, FeatureArtist) and artist.get_zorder() == 2:
            artist.remove()

    for (lat, lon), color in zip(latlons, colors):
        ax.scatter(lon, lat, s=style["map_marker_size"], color=color, edgecolor="white", linewidth=0.8,
                   transform=ccrs.PlateCarree(), zorder=10)
    return center, zoom


def _plot_publication_radar(ax, profile, center, zoom, radar, radar_time, style, height):
    """Draw radar reflectivity on the map, with a dBZ colorbar under it labeled with the scan time. `height` is the
    map's height (in), to place the colorbar in inches.

    "mosaic" is the CONUS mosaic (recent dates only), centered on the map; "single" is the NEXRAD site nearest
    `profile`. Radar that isn't found leaves the map as it is (the map helpers give a SounderPyWarning saying why).
    """

    if radar == "mosaic":
        centered = {**profile, "site_info": {**profile["site_info"], "site-latlon": list(center)}}
        found = _plot_map_mosaic_radar(ax, centered, zoom, radar_time, _map_colors(style))
    else:
        found = _plot_map_single_site_radar(ax, profile, radar_time, _map_colors(style))

    # the map helpers put the scan time (such as "KTLX: 2024-05-06 | 21:04") in a box on the map; take it off the map
    captions = [label.get_text().replace(": ", " ", 1).replace(" | ", " ") for label in ax.texts]
    for label in list(ax.texts):
        label.remove()
    caption = captions[0] if captions else ""

    if found and style["radar_colorbar"]:
        low, high = style["radar_range"]
        cax = ax.inset_axes((0, -0.15 / height, 1, 0.07 / height))   # 0.15 in under the map, 0.07 in tall
        bar = ax.figure.colorbar(ax.collections[-1], cax=cax, orientation="horizontal",
                                 boundaries=np.arange(low, high + 1, 1), ticks=np.arange(low, high + 1, 15))
        bar.ax.tick_params(labelsize=style["fontsize"] - 1, length=2, pad=1)
        bar.set_label("\n".join(line for line in ("reflectivity (dBZ)", caption) if line), fontsize=style["fontsize"] - 1,
                      labelpad=1, linespacing=1.2)
        bar.outline.set_linewidth(style["axis_linewidth"])
    elif found and caption:
        ax.text(0, 1.02, caption, transform=ax.transAxes, ha="left", va="bottom", fontsize=style["fontsize"] - 1.5)



#########################################################################
########################### PARAMETER UTILS #############################
#########################################################################

def _get_parameter_spec(parameter):
    """A parameter's (label, group, key, multiplier, unit, decimals), from its name or its own tuple."""
    if isinstance(parameter, str):
        if parameter not in PUBLICATION_PARAMETERS:
            raise ValueError(f"Unknown parameter {parameter!r}. Use one of {list(PUBLICATION_PARAMETERS)}, or a "
                             f"(label, group, key, multiplier, unit, decimals) tuple")
        return PUBLICATION_PARAMETERS[parameter]
    if len(parameter) != 6:
        raise ValueError(f"A parameter tuple needs (label, group, key, multiplier, unit, decimals). Got: {parameter!r}")
    return parameter


def _format_publication_parameter(snd_params, parameter):
    """A parameter's (label, unit, value as text), with MISSING_VALUE when it's missing."""
    label, group, key, multiplier, unit, decimals = _get_parameter_spec(parameter)
    value = snd_params[group][key]
    value = getattr(value, "m", value)
    if ma.is_masked(value) or not np.isfinite(value):
        return label, unit, MISSING_VALUE
    return label, unit, f"{value * multiplier + 0.0:.{decimals}f}"   # + 0.0 so a zero never prints as -0


def _plot_publication_parameters(ax, snd_params, parameters, colors, style, width):
    """Plot the parameter table as a strip `width` inches across: a column for each parameter, with its name and unit
    on top and each profile's value under them in its color. Each column is as wide as its widest entry, and
    the room left over is shared evenly between them."""

    ax.set_axis_off()
    renderer = ax.figure.canvas.get_renderer()
    line, gap = style["fontsize"] * 1.5, style["fontsize"]   # a row's height and the least room between columns (pt)
    columns = []
    for parameter in parameters:
        entries = [_format_publication_parameter(params, parameter) for params in snd_params]
        rows = [(entries[0][0], None), (entries[0][1], None)] + [(entry[2], color) for entry, color in zip(entries, colors)]
        columns.append([ax.annotate(text, xy=(0, 1), xycoords="axes fraction", xytext=(0, -(row + 0.2 * (row > 1)) * line),
                                    textcoords="offset points", ha="center", va="top", fontsize=style["fontsize"],
                                    color=color, annotation_clip=False) for row, (text, color) in enumerate(rows)])
    widths = [max(text.get_window_extent(renderer).width for text in column) * 72 / ax.figure.dpi for column in columns]
    spare = (width * 72 - sum(widths)) / len(columns)
    if spare < gap:
        raise ValueError(f"The {len(parameters)} parameters are too wide for the table at this width and font size; "
                         f"list fewer, or use a wider figure.")
    x = 0
    for column, column_width in zip(columns, widths):
        for text in column:
            text.xyann = (x + (spare + column_width) / 2, text.xyann[1])
        x += spare + column_width



#########################################################################
########################### PUBLICATION FIGURE ##########################
#########################################################################

def build_publication_sounding(profiles, panels=('skewt', 'hodograph', 'parameters'),
                               parameters=('mllcl', 'mlcape', 'mlcin', 'srh01', 'srh03', 'stp'),
                               width='double', panel_labels=('(a)', '(b)', '(c)', '(d)'),
                               colors=None, labels=None, storm_motion='right_moving', max_levels=250,
                               ylim=None, parcel='ml', show_adiabats=False, show_barbs=True,
                               hodo_heights=(0, 1000, 3000, 6000, 12000), hodo_max_height=12000, wind_units='m/s',
                               hodo_limits=None, storm_motion_profile=0, sr_vector_heights=(),
                               map_zoom=1.5, map_counties=True, map_radar=None, radar_time='sounding', parameter_profiles=None,
                               style=None, watermark=None, dpi=300, save=False, filename='sounderpy_publication.pdf'):
    '''
    EXPERIMENTAL: Return a clean, journal-style figure of one or more profiles: a skew-t, and any of a hodograph, a
    location map, and a table of parameters, sized for AMS and AGU journals. Each set of panels has its own predefined
    layout, where every panel fills a set space: the skew-t on the left, the hodograph and map beside it, and the table
    in a strip across the bottom. Every part of its look can be changed, and the figure and its panels
    are returned for annotating.

    :param profiles: one SounderPy data dictionary, or a list of them, drawn in ``colors``
    :type profiles: dict or list of dicts, required
    :param panels: the panels to draw: ``'skewt'`` and any of ``'hodograph'``, ``'map'`` (each profile needs ``site-latlon``), and ``'parameters'``. They're labeled in that order, whatever order they're listed in. Default is ``('skewt', 'hodograph', 'parameters')``
    :type panels: tuple of strings, optional
    :param parameters: the parameter table's columns (with the ``'parameters'`` panel): names from ``PUBLICATION_PARAMETERS``, or ``(label, group, key, multiplier, unit, decimals)`` tuples for any ``sounding_params`` value. The table is a strip across the bottom of the figure with a column for each parameter; about 8 fit at the default width and font. Default is ``('mllcl', 'mlcape', 'mlcin', 'srh01', 'srh03', 'stp')``
    :type parameters: list or None, optional
    :param width: the figure's width: ``'single'`` (3.2 in) or ``'double'`` (6.5 in), which fit AMS and AGU journals; a journal's own width from ``PUBLICATION_WIDTHS``; or inches. The height follows from the panels. Other widths scale the layout, and ``'single'`` suits a skew-t alone. Default is ``'double'``
    :type width: str or float, optional
    :param panel_labels: labels for the panels, in the order skew-t, hodograph, map, and table, or ``None`` for no labels. Default is ``('(a)', '(b)', '(c)', '(d)')``
    :type panel_labels: tuple of strings, optional
    :param colors: a color for each profile. Default is ``None``, which uses the style's colors
    :type colors: list of strings, optional
    :param labels: a legend name for each profile, such as ``['12Z OUN RAOB', 'HRRR F03']``. Default is ``None``, for no legend
    :type labels: list of strings, optional
    :param storm_motion: the storm motion used for every profile's parameters. Default is ``'right_moving'``
    :type storm_motion: str or list of floats, optional
    :param max_levels: the most vertical levels used from each profile. Default is ``250``
    :type max_levels: int or None, optional
    :param ylim: the skew-t's pressure range as ``(bottom, top)`` in hPa, such as ``(1000, 500)``. Default is ``None``, for 1050 to 100 hPa
    :type ylim: tuple of two numbers, optional
    :param parcel: the parcel trace to draw: ``'sb'``, ``'ml'``, ``'mu'``, or ``None``. Default is ``'ml'``
    :type parcel: str or None, optional
    :param show_adiabats: whether to draw dry and moist adiabats. Default is ``False``
    :type show_adiabats: bool, optional
    :param show_barbs: whether to draw a column of wind barbs for each profile beside the skew-t. Default is ``True``
    :type show_barbs: bool, optional
    :param hodo_heights: heights (m AGL) to mark on the hodograph. Default is ``(0, 1000, 3000, 6000, 12000)``
    :type hodo_heights: tuple of numbers, optional
    :param hodo_max_height: the height (m AGL) the hodograph stops at. Default is ``12000``
    :type hodo_max_height: int or float, optional
    :param wind_units: the hodograph's units: ``'m/s'`` or ``'kt'``. Default is ``'m/s'``
    :type wind_units: str, optional
    :param hodo_limits: the hodograph's limits as ``((u_min, u_max), (v_min, v_max))``. Default is ``None``, which fits the data
    :type hodo_limits: tuple, optional
    :param storm_motion_profile: which profile's storm motion to draw, or ``None`` for none. Default is ``0``
    :type storm_motion_profile: int or None, optional
    :param sr_vector_heights: heights (m AGL) to draw storm-relative wind vectors to. Default is none
    :type sr_vector_heights: tuple of numbers, optional
    :param map_zoom: half the map's shorter side, in degrees (it grows to fit every profile). Default is ``1.5``
    :type map_zoom: float, optional
    :param map_counties: whether to draw county lines on the map. Default is ``True``
    :type map_counties: bool, optional
    :param map_radar: radar reflectivity on the map: ``'single'`` (the NEXRAD site nearest the first profile), ``'mosaic'`` (the CONUS mosaic; only for about the last month), or ``None``. Default is ``None``
    :type map_radar: str or None, optional
    :param radar_time: the radar's time: ``'sounding'`` (the first profile's valid time) or ``'now'``. Default is ``'sounding'``
    :type radar_time: str, optional
    :param parameter_profiles: which profiles get a row in the parameter table. Default is ``None``, for all of them
    :type parameter_profiles: list of ints, optional
    :param style: changes to ``PUBLICATION_STYLE`` for this figure, such as ``{'fontsize': 6}`` for GRL. Default is ``None``
    :type style: dict, optional
    :param watermark: ``None`` for no watermark, ``True`` for the session's watermark text, or settings such as ``{'title': 'UND Atmospheric Sciences', 'subtitle': '', 'logo': 'logo.png'}``. Default is ``None``
    :type watermark: None, bool, or dict, optional
    :param dpi: dots per inch resolution of the figure. Default is ``300``
    :type dpi: int, optional
    :param save: whether to save the figure to ``filename``. Default is ``False``
    :type save: bool, optional
    :param filename: the file to save to; its extension sets the format, such as ``.pdf``, ``.svg``, or ``.png``. Default is ``'sounderpy_publication.pdf'``
    :type filename: str, optional
    :return: ``(fig, axes)``: the Matplotlib figure, and a dict of its panels: ``'skewt'`` (a MetPy ``SkewT``), ``'hodograph'``, ``'map'`` (a cartopy ``GeoAxes``), and ``'parameters'``, each when drawn
    :rtype: tuple
    '''

    #################################################################
    ### CHECK THE OPTIONS ###
    #################################################################
    user_style = style or {}
    style = {**PUBLICATION_STYLE, **user_style}
    profiles = [profiles] if isinstance(profiles, dict) else list(profiles)
    colors = list(style["colors"] if colors is None else colors)
    parameters = list(parameters or [])
    _check_publication_options(profiles, list(panels), width, colors, labels, parcel, wind_units, map_radar,
                               storm_motion_profile, parameter_profiles, parameters)
    panels = [panel for panel in PUBLICATION_PANELS if panel in panels]   # in order: skew-t, hodograph, map, table
    width = PUBLICATION_WIDTHS.get(width, width)
    ylim = _validate_pressure_limits(ylim) or FULL_PRESSURE_RANGE
    scale, wind_label = (KTS_TO_MS, r"m s$^{-1}$") if wind_units == "m/s" else (1, "kt")


    #################################################################
    ### SET UP THE DATA ###
    #################################################################
    data = [_prepare_sounding_data(profile, None, max_levels) for profile in profiles]
    snd_params = [sounding_params(profile, storm_motion).calc() for profile in data]
    for parameter in parameters:
        label, group, key = _get_parameter_spec(parameter)[:3]
        if key not in snd_params[0].get(group, {}):
            raise ValueError(f"Parameter {label!r}: sounding_params has no {group!r} value named {key!r}. "
                             f"Its groups are {list(snd_params[0])}")
    winds = [_get_profile_winds(profile, max_height=hodo_max_height) for profile in data]
    # a zoomed skew-t fits its temperature range to the profiles, unless the style sets one
    if ylim != FULL_PRESSURE_RANGE and "temperature_range" not in user_style:
        temperature_range = _get_zoomed_temperature_limits(data, ylim)
    else:
        temperature_range = style["temperature_range"]

    # each panel's place, from the layout for these panels, and the hodograph's limits to fill its panel
    listed = range(len(data)) if parameter_profiles is None else parameter_profiles
    rects = _get_panel_rects(panels, width, len(profiles), style, colorbar=bool(map_radar) and style["radar_colorbar"])
    size = min(min(w, h) for panel, (left, top, w, h) in rects.items() if panel != "parameters")
    if size < PUBLICATION_MIN_PANEL:
        warnings.warn(f"A panel is only {size:.1f} in across at this width; use a wider figure or fewer panels.",
                      SounderPyWarning)
    if "hodograph" in panels:
        shape = rects["hodograph"][3] / rects["hodograph"][2]
        limits = (_fit_limits(hodo_limits, shape) if hodo_limits else
                  _get_publication_hodograph_limits(winds, scale, snd_params, shape))


    with plt.rc_context({"font.family": style["font_family"]}):

        #################################################################
        ### BUILD THE FIGURE ###
        #################################################################
        fig = plt.figure(figsize=(width, max(top + h for left, top, w, h in rects.values())), dpi=dpi, facecolor="white", FigureClass=PublicationFigure)
        axes = {}


        #################################################################
        ### SKEW-T ###
        #################################################################
        skew = _setup_publication_skewt(fig, style, ylim, temperature_range, show_adiabats=show_adiabats)
        axes["skewt"] = skew
        for i, (profile, params) in enumerate(zip(data, snd_params)):
            _plot_publication_profile(skew, profile, params, colors[i], style, parcel=parcel,
                                      label=labels[i] if labels is not None else None)
            if show_barbs:
                _plot_publication_barbs(skew, profile, colors[i], i, style, ylim, size=rects["skewt"][2:])


        #################################################################
        ### HODOGRAPH ###
        #################################################################
        # the traces here; the storm motion and height markers once the figure is laid out (below)
        if "hodograph" in panels:
            axes["hodograph"] = fig.add_axes((0, 0, 1, 1))
            _setup_publication_hodograph(axes["hodograph"], limits, wind_label, style)
            for profile_winds, color in zip(winds, colors):
                axes["hodograph"].plot(profile_winds["u"] * scale, profile_winds["v"] * scale, color=color,
                                       linewidth=style["hodograph_linewidth"], clip_on=False, zorder=3)


        #################################################################
        ### MAP ###
        #################################################################
        if "map" in panels:
            axes["map"] = fig.add_axes((0, 0, 1, 1), projection=ccrs.PlateCarree())
            shape = rects["map"][3] / rects["map"][2]
            center, zoom = _plot_publication_map(axes["map"], data, colors, style, shape, zoom=map_zoom,
                                                 show_counties=map_counties)
            if map_radar:
                _plot_publication_radar(axes["map"], data[0], center, zoom, map_radar, radar_time, style,
                                        height=rects["map"][3])


        #################################################################
        ### PARAMETER TABLE ###
        #################################################################
        if "parameters" in panels:
            axes["parameters"] = fig.add_axes((0, 0, 1, 1))
            _plot_publication_parameters(axes["parameters"], [snd_params[i] for i in listed], parameters,
                                         [colors[i] for i in listed], style, width=rects["parameters"][2])


        #################################################################
        ### LAY OUT THE FIGURE ###
        #################################################################
        # the panels labeled; then everything fitted inside the figure at the journal width, with the legend at the
        # bottom and fewer temperature ticks if the skew-t is too narrow for them
        _add_panel_labels(axes, panel_labels, style)
        if watermark:
            _add_publication_watermark(fig, watermark, style)
        legend = (skew, style) if labels is not None else None
        _fit_figure(fig, axes, rects, legend)
        if _set_temperature_ticks(skew, style):
            _fit_figure(fig, axes, rects, legend)

        # the storm motion and height markers, measured at the hodograph's final size
        if "hodograph" in panels:
            axes["hodograph"].apply_aspect()
            if storm_motion_profile is not None:
                _plot_publication_storm_motion(axes["hodograph"], snd_params[storm_motion_profile], winds[storm_motion_profile],
                                               scale, style, sr_vector_heights, marker_heights=hodo_heights)
            placed = []   # markers drawn so far, so every profile's markers stay clear of each other
            for profile_winds, color in zip(winds, colors):
                _plot_publication_height_markers(axes["hodograph"], profile_winds, scale, color, style, hodo_heights, placed)

        # saved at the figure's exact size
        if save:
            fig.savefig(filename)

    return fig, axes

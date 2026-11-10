"""
    SOUNDERPY WINTER PLOT COMPONENTS

    Purpose of module:

    Drawing tools for the winter sounding: SounderPy's full sounding layout with winter content in each slot.
    The skew-t shows melting and refreezing layers, precipitation phase changes, 0 °C levels, clouds,
    icing and the DGZ; the hodograph shows ML and DGZ mean winds, a winds box and the latest surface
    observation; the side panels show T & Tw, RHᵢ & RH, wind and theta-e (0-3 km); the tables show precipitation
    type, snow growth and winter hazards. SounderPy's own plot components are reused wherever they fit.

    (C) KYLE J GILLETT, 2026
"""

# NUMPY OPERATIONS
import numpy as np
import numpy.ma as ma

# MATPLOTLIB OPERATIONS
from matplotlib.patches import Rectangle

# METPY OPERATIONS
import metpy.calc as mpcalc
from metpy.units import units

# SOUNDERPY MODULES
from ..utils import wind_to_dir
from .style import (FONT_SIZES, HAZARD_LEVEL_COLORS, HAZARDS_LAYOUT, LAYER_TABLE_LAYOUT, PANEL_TOP,
                    PARCEL_TABLE_LAYOUT, PTYPE_COLORS, PTYPE_HEADER_LAYOUT, PTYPE_NAMES, SHARPPY_TO_CODE,
                    SKEWT_COLUMNS, SKY_COVER_CODES, SNOW_GROWTH_LAYOUT, STATION_BOX, UNIT_JKG, WINDS_BOX_LAYOUT)
from .components import _add_skewt_layer_line, _draw_text, _format_value, _setup_height_panel, _setup_text_box




#########################################################################
############################## HELPERS ##################################
#########################################################################

def _p_at(z, p, height):
    """Pressure at a height (m AGL) on the data levels."""
    return float(np.interp(height, z, p))


def _fmt(value, unit="", decimals=0):
    """A value with its unit, dropping a trailing ".0" when no decimals are asked for."""
    text = _format_value(value, unit, decimals=decimals)
    return text.replace(".0" + unit, unit) if decimals == 0 else text


def _fmt_signed(value, unit="", decimals=1):
    """A change, always with its sign: +6.0°C, -2.3°C."""
    if ma.is_masked(value) or not np.isfinite(float(value)):
        return "NA"
    return f"{float(value):+.{decimals}f}{unit}"


def _wind_text(u, v):
    """A wind as "NNW @ 21 kt", named for the direction it blows from."""
    # uses SounderPy's wind_to_dir() like the severe figure; it names the direction a vector points toward
    # (right for storm motions), so the wind is negated
    if ma.is_masked(u) or ma.is_masked(v):
        return "NA"
    return f"{wind_to_dir(-float(u), -float(v))} @ {np.hypot(u, v):.0f} kt"


def _slr_text(winter):
    """Kuchera snow ratio, e.g. "12:1"."""
    if ma.is_masked(winter["slr"]):
        return "NA"
    return f"{winter['slr']:.0f}:1"


def _fit_width(text, ax, max_width):
    """Shrink a text to fit a width (fraction of its axis); returns the text."""
    width = text.get_window_extent(ax.figure.canvas.get_renderer()).transformed(ax.transAxes.inverted()).width
    if width > max_width:
        text.set_fontsize(text.get_fontsize() * max_width / width)
    return text


def _shade_dgz(ax, winter, theme):
    """Shade the DGZ on a height panel, and mark the mixed-layer top."""
    if not ma.is_masked(winter["dgz_depth"]):
        ax.axhspan(*winter["dgz_z"], color=theme["dgz_fill"], alpha=0.3, linewidth=0)
    _mark_ml_top(ax, winter, theme)


def _mark_ml_top(ax, winter, theme):
    """Mark the mixed-layer top on a height panel: a dashed line labeled "ML" at the right edge."""
    depth = winter.get("ml_depth", ma.masked)
    if ma.is_masked(depth) or not 0 < depth < PANEL_TOP:
        return
    ax.axhline(depth, color=theme["lapse"], linewidth=1.5, linestyle="--", alpha=0.85, zorder=5)
    # under the line, clear of the km labels above it (above the line when the layer is too shallow)
    below = depth > 200
    ax.text(0.97, depth - 25 if below else depth + 25, "ML", transform=ax.get_yaxis_transform(), ha="right",
            va="top" if below else "bottom", color=theme["lapse"], weight="bold", fontsize=FONT_SIZES["small"],
            alpha=0.9, zorder=6)


def _section_title(ax, text, theme, x, y):
    """A faded section title in a table."""
    _draw_text(ax, x, y, text, theme["text"], fontsize=FONT_SIZES["small"], alpha=0.55)




#########################################################################
############################## SKEW-T ###################################
#########################################################################

def _plot_melt_refreeze(skew, data, winter, theme):
    """Shade melting (Tw > 0 °C) and refreezing layers on the skew-t, labeling each layer's energy."""
    # each layer's energy is labeled just right of the temperature trace, level with the middle of the layer
    layers = winter["layers"]
    if not layers:
        return
    p = data["p"].m
    z = data["z"].m - data["z"].m[0]
    T = data["T"].to("degC").m
    Tw = np.asarray(winter["Tw_levels"], dtype=float)

    below_source = z <= max(lyr["top"] for lyr in layers) + 1
    lowest_warm_aloft = min((lyr["bottom"] for lyr in layers if lyr["kind"] == "warm" and not lyr["surface"]), default=None)
    skew.ax.fill_betweenx(p, Tw, 0, where=below_source & (Tw > 0), interpolate=True, color=theme["melt"], alpha=0.28,
                          linewidth=0, zorder=2)
    if lowest_warm_aloft is not None:
        refreeze = below_source & (Tw < 0) & (z <= lowest_warm_aloft + 1)
        skew.ax.fill_betweenx(p, Tw, 0, where=refreeze, interpolate=True, color=theme["refreeze"], alpha=0.28,
                              linewidth=0, zorder=2)

    for lyr in layers:
        # only melting layers and the refreezing layers under them
        if lyr["kind"] == "cold" and (lowest_warm_aloft is None or lyr["bottom"] >= lowest_warm_aloft):
            continue
        if abs(lyr["energy"]) < 0.5:
            continue
        # the trace leans right with height, so find its rightmost point in the layer in display space
        p_mid = _p_at(z, p, (lyr["bottom"] + lyr["top"]) / 2)
        in_layer = (z >= lyr["bottom"] - 1) & (z <= lyr["top"] + 1)
        p_layer = np.append(p[in_layer], [_p_at(z, p, lyr["bottom"]), _p_at(z, p, lyr["top"])])
        T_layer = np.interp(-np.log(p_layer), -np.log(p), T)
        right_px = skew.ax.transData.transform(np.column_stack([T_layer, p_layer]))[:, 0].max()
        anchor = skew.ax.transData.inverted().transform((right_px, skew.ax.transData.transform((0, p_mid))[1]))
        skew.ax.annotate(f"{lyr['energy']:+.0f} J/kg", xy=anchor, xytext=(10, 0), textcoords="offset points",
                         color=theme["melt"] if lyr["kind"] == "warm" else theme["refreeze"], weight="bold",
                         va="center", ha="left", fontsize=FONT_SIZES["annotation"], zorder=8, clip_on=True)


def _plot_freezing_levels(skew, data, winter, theme):
    """Label every 0 °C crossing of T (FRZ) and of Tw (WB0), in their own columns."""
    p = data["p"].m
    z = data["z"].m - data["z"].m[0]
    transform = skew.ax.get_yaxis_transform()
    for key, label, x, color, alpha in (("frz_pts_z", "←FRZ", SKEWT_COLUMNS["frz"], theme["frz"], 0.8),
                                        ("wb_frz_pts_z", "←WB0", SKEWT_COLUMNS["wb_frz"], theme["wb_frz"], 0.7)):
        for height in winter[key]:
            skew.ax.text(x, _p_at(z, p, height), label, color=color, weight="bold", alpha=alpha, va="center",
                         fontsize=FONT_SIZES["annotation"], transform=transform, clip_on=True)


def _plot_labeled_bar(skew, x, p_bottom, p_top, label, color, alpha, label_alpha, linewidth, label_offset, label_va, fontsize):
    """A labeled bar with end ticks between two pressures, in SounderPy's layer-bar style."""
    transform = skew.ax.get_yaxis_transform()
    skew.ax.text(x, p_top - label_offset, label, weight="bold", color=color, alpha=label_alpha, ha="center",
                 va=label_va, fontsize=fontsize, transform=transform, clip_on=True)
    _add_skewt_layer_line(skew, x=x, p_bottom=p_bottom, p_top=p_top, color=color, alpha=alpha, linewidth=linewidth)
    for pressure in (p_bottom, p_top):
        skew.ax.plot([x - 0.01, x + 0.01], [pressure, pressure], color=color, alpha=alpha, linewidth=max(linewidth - 1, 1.5),
                     transform=transform, clip_on=True)


def _plot_dgz_bar(skew, winter, snd, theme):
    """Draw the DGZ bar on the skew-t."""
    # always the DGZ: SounderPy's growth-zone helper shows the hail growth zone when the surface is >= 5 °C
    if not ma.is_masked(winter["dgz_depth"]):
        p_bottom, p_top = winter["dgz_p"]
        _plot_labeled_bar(skew, SKEWT_COLUMNS["dgz"], p_bottom, p_top, "DGZ", theme["dgz"], alpha=0.4, label_alpha=0.6,
                          linewidth=3, label_offset=12, label_va="baseline", fontsize=FONT_SIZES["tick"])


def _plot_pbl(skew, data, snd, theme):
    """Draw the PBL bar, labeled above its top tick, in its own column."""
    pbl_top = snd["thermo"]["pbl_top"]
    if not ma.is_masked(pbl_top):
        _plot_labeled_bar(skew, SKEWT_COLUMNS["pbl"], data["p"][0].m, pbl_top, "PBL", theme["pbl"], alpha=0.95,
                          label_alpha=0.95, linewidth=1, label_offset=6, label_va="bottom", fontsize=FONT_SIZES["small"])


def _plot_ml_bar(skew, data, winter, theme):
    """Draw the mixed-layer bar, labeled like the PBL bar, in its own column beside it."""
    if not ma.is_masked(winter["ml_top_p"]):
        _plot_labeled_bar(skew, SKEWT_COLUMNS["ml"], data["p"][0].m, winter["ml_top_p"], "ML", theme["lapse"], alpha=0.95,
                          label_alpha=0.95, linewidth=1, label_offset=6, label_va="bottom", fontsize=FONT_SIZES["small"])


def _plot_winter_parcels(skew, snd, theme):
    """Plot the SB and MU parcel traces, when they have CAPE."""
    thermo = snd["thermo"]
    for key in ("sb", "mu"):
        cape = thermo[f"{key}cape"]
        if not ma.is_masked(cape) and cape > 0:
            skew.plot(thermo[f"{key}P_trace"], thermo[f"{key}T_trace"], linestyle="--", linewidth=2, alpha=0.75,
                      color=theme["parcel"][key], label=f"{key.upper()}CAPE PARCEL")


def _plot_warm_nose(skew, winter, theme):
    """Mark the warm nose on the temperature trace."""
    if not ma.is_masked(winter["warm_nose_T"]):
        skew.ax.plot(winter["warm_nose_T"], winter["warm_nose_p"], marker="o", markersize=9, color=theme["warm_nose"],
                     markeredgecolor=theme["background"], zorder=9)


def _plot_phase_annotations(skew, data, winter, theme):
    """Label the precipitation phase on the temperature trace wherever it changes on the way down."""
    # at the source, under each layer that changes it, and at the surface when it freezes on contact.
    # snow carries its ratio
    p = data["p"].m
    T = data["T"].to("degC").m
    for step in winter["phase_trace"]:
        # the surface label sits a little above the ground, clear of the surface T/Td labels
        pressure = p[0] - 12 if step["z"] <= 1 else step["p"]
        temperature = float(np.interp(-np.log(pressure), -np.log(p), T))
        phase = step["phase"]
        label = f"← {phase}" + (f" {_slr_text(winter)}" if phase == "SN" and not ma.is_masked(winter["slr"]) else "")
        color = PTYPE_COLORS.get(phase, theme["text"])
        skew.ax.plot(temperature, pressure, marker="o", markersize=7, color=color, markeredgecolor=theme["background"], zorder=10)
        skew.ax.annotate(label, xy=(temperature, pressure), xytext=(9, 0), textcoords="offset points", color=color,
                         weight="bold", fontsize=FONT_SIZES["annotation_large"], ha="left", va="center", zorder=10,
                         clip_on=True)


def _plot_rh_ice_annotations(skew, snd, winter, theme, heights=(1000, 3000, 6000, 9000, 12000)):
    """Label RH with respect to ice (RH over water at or above 0 °C) at fixed heights, left of the dewpoint trace."""
    for height in heights:
        i = int(height // 100)
        if i >= len(winter["z"]) or not np.isfinite(winter["rh"][i]):
            continue
        below_freezing = winter["T"][i] < 0
        value = winter["rh_ice"][i] if below_freezing else winter["rh"][i]
        dewpoint = mpcalc.dewpoint_from_relative_humidity(winter["T"][i] * units.degC, max(winter["rh"][i], 1) / 100).m
        # two short lines (RHᵢ over the value) rather than one long one
        skew.ax.annotate(f"{'RHᵢ' if below_freezing else 'RH'}\n{value:.0f}% →", (dewpoint, snd["intrp"]["pINTRP"][i]),
                         xytext=(-8, 0), textcoords="offset points", alpha=0.8, weight="bold", multialignment="center",
                         linespacing=1.0,
                         color=theme["rh_ice"] if below_freezing else theme["dewpoint"],
                         fontsize=FONT_SIZES["annotation"], ha="right", va="center", clip_on=True)


def _plot_layer_bars(skew, layers, x, color, alpha, label, side):
    """Draw a bar for each layer, each labeled along its side and centered on its depth."""
    for layer in layers:
        _add_skewt_layer_line(skew, x=x, p_bottom=layer["bottom_p"], p_top=layer["top_p"], color=color,
                              alpha=alpha, linewidth=5, linestyle="-")
        skew.ax.annotate(label, xy=(x, np.sqrt(layer["bottom_p"] * layer["top_p"])), xycoords=skew.ax.get_yaxis_transform(),
                         xytext=(5 if side == "right" else -5, 0), textcoords="offset points", rotation=90,
                         ha="left" if side == "right" else "right", va="center", color=color, weight="bold",
                         alpha=0.85, fontsize=FONT_SIZES["fine"] + 1, annotation_clip=True, zorder=9)


def _plot_clouds(skew, winter, theme):
    """Draw cloud bars labeled on their left, with the top temperature of the highest cloud."""
    x = SKEWT_COLUMNS["cloud"]
    _plot_layer_bars(skew, winter["clouds"], x, theme["cloud"], 0.7, "CLOUD", side="left")
    if winter["clouds"]:
        top = winter["clouds"][-1]
        skew.ax.text(x, top["top_p"] - 8, f"CLD TOP\n{top['top_T']:.0f}°C", color=theme["cloud"], weight="bold",
                     ha="center", va="bottom", linespacing=1.1, fontsize=FONT_SIZES["tick"],
                     transform=skew.ax.get_yaxis_transform(), clip_on=True, zorder=9)


def _plot_icing_layers(skew, winter, theme):
    """Draw icing-layer bars beside the cloud bars, each labeled on its right, with one legend entry."""
    if winter["icing_layers"]:
        _plot_layer_bars(skew, winter["icing_layers"], SKEWT_COLUMNS["icing"], theme["icing"], 0.45, "ICING", side="right")
        skew.ax.plot([], [], color=theme["icing"], alpha=0.45, linewidth=5, label="ICING LAYER")




#########################################################################
############################## HODOGRAPH ################################
#########################################################################

def _plot_mean_wind_vectors(hodo, winter, theme):
    """Draw the ML and DGZ mean-wind vectors from the hodograph origin."""
    for key, label in (("mean_wind_ml", "ML"), ("mean_wind_dgz", "DGZ")):
        u, v = winter[key]
        if ma.is_masked(u) or ma.is_masked(v):
            continue
        hodo.ax.annotate("", xy=(u, v), xytext=(0, 0), zorder=6,
                         arrowprops=dict(arrowstyle="-|>", color=theme["mean_wind"], linewidth=3, alpha=0.85,
                                         mutation_scale=22, shrinkA=0, shrinkB=0))
        hodo.ax.annotate(label, xy=(u, v), xytext=(6, 6), textcoords="offset points", color=theme["text"],
                         weight="bold", fontsize=FONT_SIZES["annotation"], zorder=7, clip_on=True)


def _plot_winds_box(ax, winter, theme):
    """Draw the winds box in the hodograph's lower-right corner."""
    layout = WINDS_BOX_LAYOUT
    _setup_text_box(ax, theme)
    _section_title(ax, "WINDS", theme, layout["label_x"], layout["title_y"])
    rows = [("10 m:", _wind_text(*winter["sfc_wind_uv"])),
            ("ML mean:", _wind_text(*winter["mean_wind_ml"])),
            ("DGZ mean:", _wind_text(*winter["mean_wind_dgz"])),
            ("Max 0-2 km:", f"{_fmt(winter['max_wind_0_to_2000'])} kt")]
    # round(...) + 0 keeps a -0.4 °C from printing as "-0"
    rows += [(f"{level} hPa:", f"{_wind_text(w['u'], w['v'])}, {round(w['T']) + 0:.0f}°C")
             for level, w in winter["level_winds"].items()]
    for y, (label, value) in zip(np.linspace(layout["rows_top"], layout["rows_bottom"], len(rows)), rows):
        _draw_text(ax, layout["label_x"], y, label, theme["text"], fontsize=FONT_SIZES["label"], alpha=0.9)
        _draw_text(ax, layout["value_x"], y, value, theme["kinematics_strong"], fontsize=FONT_SIZES["label"], ha="right")
    return ax


def _plot_station_obs(hodo, obs, theme):
    """Draw a compact station model for the latest surface observation in the hodograph's upper-right corner."""
    # T and Td (left), coded SLP and gust (right), present weather and visibility, sky cover and wind barb.
    # Plain Matplotlib text plus MetPy's symbol font: MetPy 1.7.1's StationPlot text fails under Matplotlib 3.11
    from metpy.plots import current_weather, sky_cover, wx_code_to_numeric
    from metpy.plots.wx_symbols import wx_symbol_font
    title = f"Latest Observation\n{obs['station']} ({obs['time']:%H%MZ}) · {obs['distance_km']:.0f} km"

    # as wide as the title and STATION_BOX["aspect"] times as tall, in screen pixels. Measured against the
    # hodograph's final size (its fixed aspect shrinks it at draw time), right edge flush with the figure's edge
    hodo.ax.apply_aspect()
    probe = hodo.ax.text(0, 0, title, fontsize=FONT_SIZES["small"], weight="bold", linespacing=1.25)
    side_px = max(probe.get_window_extent(hodo.ax.figure.canvas.get_renderer()).width + 2 * STATION_BOX["pad_px"],
                  STATION_BOX["min_px"])
    probe.remove()
    width, height = side_px / hodo.ax.bbox.width, side_px * STATION_BOX["aspect"] / hodo.ax.bbox.height
    position = hodo.ax.get_position()
    left = (STATION_BOX["right"] - position.x0) / position.width - width
    ax = hodo.ax.inset_axes((left, STATION_BOX["top"] - height, width, height))
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor(theme["background"])
    ax.patch.set_alpha(0.92)
    for spine in ax.spines.values():
        spine.set_color(theme["text"])
        spine.set_alpha(0.35)
    ax.set_zorder(9)
    ax.text(0, 0.785, title, ha="center", va="center", fontsize=FONT_SIZES["small"], weight="bold", linespacing=1.25,
            color=theme["text"], alpha=0.8)
    ax.plot([-0.86, 0.86], [0.575, 0.575], color=theme["text"], alpha=0.2, linewidth=1)

    # the station circle sits centered in the space under the title (room for a barb up or down) and a little
    # right of center (room for weather and visibility on its left); dy keeps T/Td the same on-screen distance
    # from the circle as in a square box
    cx, cy, dx, dy = 0.14, -0.21, 0.19, 0.3 / STATION_BOX["aspect"]
    size = FONT_SIZES["annotation"]
    text = dict(weight="bold", fontsize=size, va="center", zorder=3)
    symbol = dict(fontproperties=wx_symbol_font, va="center", ha="center", zorder=3)
    # a value the report doesn't have reads NA (no gust or weather is a real "none", so those stay blank)
    ax.text(cx - dx, cy + dy, "NA" if obs["tmpf"] is None else f"{obs['tmpf']:.0f}", color=theme["temperature"],
            ha="right", **text)
    ax.text(cx - dx, cy - dy, "NA" if obs["dwpf"] is None else f"{obs['dwpf']:.0f}", color=theme["dewpoint"],
            ha="right", **text)
    ax.text(cx + dx, cy + dy, "NA" if obs["mslp"] is None else f"{int(round(obs['mslp'] * 10)) % 1000:03d}",
            color=theme["text"], ha="left", **text)
    if obs["gust"] is not None:
        ax.text(cx + dx, cy - dy, f"G{obs['gust']:.0f}", color=theme["kinematics_strong"], ha="left", **text)

    x_left = cx - dx
    wx = obs["wxcodes"] if isinstance(obs["wxcodes"], str) else ""
    code = int(wx_code_to_numeric([wx])[0]) if wx else 0
    if code:
        ax.text(x_left - 0.07, cy, current_weather(code), color=theme["text"], size=size + 3, **symbol)
        x_left -= 0.22
    ax.text(x_left, cy, "NA" if obs["vsby"] is None else f"{obs['vsby']:g}", color=theme["text"], ha="right", **text)

    # sky cover in the circle; the barb draws over the values, as on a standard station model
    covers = [SKY_COVER_CODES[str(c).upper()] for c in (obs["skyc1"], obs["skyc2"], obs["skyc3"])
              if c is not None and str(c).upper() in SKY_COVER_CODES]
    ax.text(cx, cy, sky_cover(max(covers) if covers else 10), color=theme["text"], size=size + 2, **symbol)
    if obs["sknt"] is not None and obs["drct"] is not None and obs["sknt"] > 0:
        direction = np.radians(obs["drct"])
        ax.barbs([cx], [cy], [-obs["sknt"] * np.sin(direction)], [-obs["sknt"] * np.cos(direction)], color=theme["text"],
                 length=9, linewidth=1.5, zorder=4)
    return ax




#########################################################################
############################ SIDE PANELS ################################
#########################################################################

def _setup_winter_panel(ax, theme, x_limits, ticks, title, tick_labels=None):
    """Set up a 0-3 km height panel with its title and km labels every 500 m."""
    tick_labels = tick_labels or [f"{v:.0f}" for v in ticks]
    _setup_height_panel(ax, theme, x_limits, ticks, tick_labels, z_limits=(0, PANEL_TOP), z_interval=500)
    ax.text(0.5, 0.98, title, transform=ax.transAxes, ha="center", va="top", weight="bold", color=theme["text"],
            fontsize=FONT_SIZES["small"], alpha=0.9, zorder=20)
    for height in range(500, PANEL_TOP, 500):
        ax.text(0.03, height + 20, f"{height / 1000:g} km".lstrip("0"), transform=ax.get_yaxis_transform(),
                fontsize=FONT_SIZES["small"], weight="bold", color=theme["text"], alpha=0.8)


def _plot_temperature_panel(ax, winter, theme):
    """Plot the T and Tw panel, with melting shaded and an axis that fits the data."""
    sel = winter["z"] <= PANEL_TOP
    z, T, Tw = winter["z"][sel], winter["T"][sel], winter["Tw"][sel]
    low = min(-25.0, np.floor(np.nanmin(np.r_[T, Tw]) / 5) * 5 - 2)
    high = max(15.0, np.ceil(np.nanmax(T) / 5) * 5 + 2)
    _setup_winter_panel(ax, theme, (low, high), [v for v in (-30, -15, 0, 10, 20) if low < v < high], "T & Tw\n(°C)")
    ax.axvline(0, color=theme["frz"], linewidth=1.5, alpha=0.7)
    ax.fill_betweenx(z, T, 0, where=T > 0, interpolate=True, color=theme["melt"], alpha=0.25, linewidth=0)
    ax.plot(T, z, color=theme["temperature"], linewidth=3)
    ax.plot(Tw, z, color=theme["wetbulb"], linewidth=2.5)
    _shade_dgz(ax, winter, theme)


def _plot_rh_panel(ax, winter, theme):
    """Plot the RHᵢ (solid) and RH (dashed) panel, with clouds shaded."""
    sel = winter["z"] <= PANEL_TOP
    _setup_winter_panel(ax, theme, (0, 125), [25, 50, 75, 100], "RHᵢ & RH\n(%)")
    ax.axvline(100, color=theme["text"], linewidth=1, alpha=0.4)
    ax.axvline(90, color=theme["cloud"], linewidth=1, alpha=0.6, linestyle="--")
    for cloud in winter["clouds"]:
        ax.axhspan(cloud["bottom"], min(cloud["top"], PANEL_TOP), color=theme["cloud"], alpha=0.15, linewidth=0)
    ax.plot(winter["rh"][sel], winter["z"][sel], color=theme["moisture"], linewidth=1.8, linestyle="--", alpha=0.8, zorder=3)
    ax.plot(np.clip(winter["rh_ice"][sel], 0, 125), winter["z"][sel], color=theme["rh_ice"], linewidth=3, zorder=4)
    _shade_dgz(ax, winter, theme)


def _plot_wind_panel(ax, winter, theme):
    """Plot the wind speed panel (kt), with barbs every 500 m along the right edge."""
    sel = winter["z"] <= PANEL_TOP
    z, u, v = winter["z"][sel], winter["u"][sel], winter["v"][sel]
    speed = np.hypot(u, v)
    top = max(30, int(np.ceil((np.nanmax(speed) + 5) / 10) * 10))
    top += 10 if top <= 40 else 20  # room for one tick past the strongest wind
    high = top * 1.45  # room on the right for the barbs
    # at most two labeled ticks
    step = next(s for s in (10, 20, 25, 50, 100) if top // s <= 2)
    ticks = list(range(step, top + 1, step))
    _setup_winter_panel(ax, theme, (0, high), ticks, "Wind\n(kt)")
    ax.plot(speed, z, color=theme["kinematics_strong"], linewidth=2.5, zorder=3)
    levels = np.arange(250, PANEL_TOP - 500, 500)
    ax.barbs(np.full(levels.size, high * 0.83), levels, np.interp(levels, z, u), np.interp(levels, z, v),
             color=theme["text"], length=5.5, linewidth=1, zorder=4, clip_on=False)
    _shade_dgz(ax, winter, theme)


def _plot_thetae_panel(ax, winter, theme):
    """Plot the theta-e panel."""
    sel = winter["z"] <= PANEL_TOP
    thetae = winter["thetae"][sel]
    low, high = np.floor(np.nanmin(thetae) - 2), np.ceil(np.nanmax(thetae) + 2)
    _setup_winter_panel(ax, theme, (low, high), np.linspace(low, high, 4)[1:-1].round(), "θe\n(K)")
    ax.plot(thetae, winter["z"][sel], color=theme["theta"], linewidth=2.5)
    _shade_dgz(ax, winter, theme)




#########################################################################
############################### TABLES ##################################
#########################################################################

def _plot_ptype_header(ax, winter, theme):
    """Draw the precipitation-type band across the top of the thermo box, tinted in the precip-type color."""
    # four centered columns of two lines:  Estimated SFC PTYPE | type (Bourgouin/Birk) | type (SPC) | SOURCE hPa · km · °C
    layout = PTYPE_HEADER_LAYOUT
    color = PTYPE_COLORS.get(winter["ptype"], theme["text"])
    ax.add_patch(Rectangle((0, layout["band_bottom"]), 1, 1 - layout["band_bottom"], transform=ax.transAxes,
                           facecolor=color, alpha=0.10, linewidth=0, zorder=0))
    ax.axhline(layout["band_bottom"], color=color, alpha=0.8, linewidth=2)

    sharppy = winter["sharppy_ptype"]
    if ma.is_masked(winter["source_p"]):
        source = "none saturated"
    else:
        source = f"{winter['source_p']:.0f} hPa · {winter['source_z'] / 1000:.1f} km · {winter['source_T']:+.0f}°C"
    # (top line, bottom line) per column, each as (text, color, font size, alpha)
    columns = [(("Estimated", theme["text"], FONT_SIZES["tick"], 0.75), ("SFC PTYPE", theme["text"], FONT_SIZES["tick"], 0.75)),
               ((PTYPE_NAMES.get(winter["ptype"], winter["ptype"]), color, FONT_SIZES["table_large"] + 3, 1.0),
                ("(Bourgouin/Birk)", theme["text"], FONT_SIZES["small"], 0.6)),
               ((sharppy.upper(), PTYPE_COLORS.get(SHARPPY_TO_CODE.get(sharppy, ""), theme["text"]), FONT_SIZES["table"] + 1, 1.0),
                ("(SPC)", theme["text"], FONT_SIZES["small"], 0.6)),
               (("SOURCE", theme["text"], FONT_SIZES["tick"], 0.75), (source, theme["source"], FONT_SIZES["small"] + 1, 1.0))]
    texts = []
    for x, lines in zip(layout["columns"], columns):
        for y, (text, text_color, size, alpha), max_width in zip((layout["top_y"], layout["bottom_y"]), lines, (0.23, 0.26)):
            # centered in its column, shrunk if it would spill into the next one
            texts.append(_fit_width(ax.text(x, y, text, color=text_color, fontsize=size, alpha=alpha, weight="bold",
                                            ha="center", va="baseline", transform=ax.transAxes), ax, max_width))
    # the two precip types (Bourgouin, SPC) share one size: the smaller of the two after fitting
    bourgouin, spc = texts[2], texts[4]
    size = min(bourgouin.get_fontsize(), spc.get_fontsize())
    bourgouin.set_fontsize(size)
    spc.set_fontsize(size)
    for x in layout["dividers"]:
        ax.plot([x, x], [layout["band_bottom"] + 0.03, 0.97], color=theme["text"], alpha=0.2, linewidth=1, transform=ax.transAxes)


def _plot_layer_table(ax, winter, theme):
    """Draw the melting and refreezing layers below the source (wet-bulb based), top first."""
    layout = LAYER_TABLE_LAYOUT
    text, small = theme["text"], FONT_SIZES["tick"]
    for x, header in zip(layout["columns"], ["LAYER", "BASE", "TOP", "DEPTH", "MAX/MIN Tw", "ENERGY (Tw)"]):
        _draw_text(ax, x, layout["header_y"], header, text, fontsize=FONT_SIZES["small"], alpha=0.7)

    if ma.is_masked(winter["source_p"]):
        _draw_text(ax, 0.015, layout["rows_y"][0], "no saturated source (no precip expected)", text, fontsize=small, alpha=0.6)
        return
    if not winter["layers"]:
        _draw_text(ax, 0.015, layout["rows_y"][0], "no melting or refreezing layer below the source", text, fontsize=small, alpha=0.6)
        return

    highest_warm = max((lyr["top"] for lyr in winter["layers"] if lyr["kind"] == "warm"), default=None)
    for y, lyr in zip(layout["rows_y"], winter["layers"]):
        # a cold layer above every warm layer is where snow forms, not where melted precip refreezes
        if lyr["kind"] == "warm":
            name, color = ("Sfc melt" if lyr["surface"] else "Melt"), theme["melt"]
        elif highest_warm is None or lyr["bottom"] >= highest_warm:
            name, color = ("Sfc cold" if lyr["surface"] else "Cold aloft"), text
        else:
            name, color = ("Sfc refrz" if lyr["surface"] else "Refreeze"), theme["refreeze"]
        values = [name, f"{lyr['bottom']:.0f} m", f"{lyr['top']:.0f} m", f"{lyr['depth']:.0f} m",
                  f"{lyr['extreme']:+.1f}°C", f"{lyr['energy']:+.0f} J/kg"]
        for x, value in zip(layout["columns"], values):
            _draw_text(ax, x, y, value, color, fontsize=small)


def _plot_parcel_table(ax, snd, theme):
    """Draw the SB and MU parcel table: CAPE, CIN, LCL, LFC, and EL."""
    layout = PARCEL_TABLE_LAYOUT
    thermo = snd["thermo"]
    ax.axhline(layout["divider_y"], xmin=0.012, xmax=0.988, color=theme["text"], alpha=0.45, linewidth=1.2)
    for x, header in zip(layout["columns"], ["CAPE", "CIN", "LCL", "LFC", "EL"]):
        _draw_text(ax, x, layout["header_y"], header, theme["text"], fontsize=FONT_SIZES["small"], alpha=0.7)
    for y, key in zip(layout["rows_y"], ("sb", "mu")):
        cape = thermo[f"{key}cape"]
        has_cape = not ma.is_masked(cape) and cape > 0
        values = [_format_value(cape, f" {UNIT_JKG}"), _format_value(thermo[f"{key}cin"], f" {UNIT_JKG}"),
                  _format_value(thermo[f"{key}_lcl_z"], " m"),
                  _format_value(thermo[f"{key}_lfc_z"], " m") if has_cape else "NA",
                  _format_value(thermo[f"{key}_el_z"], " m") if has_cape else "NA"]
        _draw_text(ax, layout["label_x"], y, f"{key.upper()}:", theme["text"])
        for x, value in zip(layout["columns"], values):
            _draw_text(ax, x, y, value, theme["parcel"][key])


def _plot_precip_table(ax, winter, snd, theme):
    """Draw the thermo box: precipitation-type band, melting and refreezing layers, then SB and MU parcels."""
    _setup_text_box(ax, theme)
    _plot_ptype_header(ax, winter, theme)
    _plot_layer_table(ax, winter, theme)
    _plot_parcel_table(ax, snd, theme)
    return ax


def _plot_snow_growth_table(ax, winter, snd, theme):
    """Draw the snow growth (DGZ) box, then every 0 °C crossing of T and Tw."""
    layout = SNOW_GROWTH_LAYOUT
    _setup_text_box(ax, theme)
    _section_title(ax, "SNOW GROWTH (DGZ)", theme, layout["label_x"], layout["title_y"])
    has_dgz = not ma.is_masked(winter["dgz_depth"])
    rows = [("DGZ:", f"{winter['dgz_z'][0] / 1000:.1f}-{winter['dgz_z'][1] / 1000:.1f} km" if has_dgz else "NA", theme["dgz"]),
            ("Depth:", _fmt(winter["dgz_depth"], " m"), theme["dgz"]),
            ("RH / RHᵢ:", f"{winter['dgz_rh']:.0f} / {winter['dgz_rh_ice']:.0f}%" if not ma.is_masked(winter["dgz_rh"]) else "NA", theme["rh_ice"]),
            ("Omega:", _fmt_signed(winter["dgz_omega"], " μb/s"), theme["ascent"]),
            ("OPRH:", _fmt(winter["oprh"], decimals=2), theme["composite_index"]),
            ("Kuchera SLR:", _slr_text(winter), PTYPE_COLORS["SN"]),
            ("PWAT:", _format_value(snd["general"]["pwat"], " in", decimals=2), theme["moisture_strong"])]
    for y, (label, value, color) in zip(np.linspace(layout["rows_top"], layout["rows_bottom"], len(rows)), rows):
        _draw_text(ax, layout["label_x"], y, label, theme["text"], fontsize=FONT_SIZES["annotation"], alpha=0.9)
        _draw_text(ax, layout["value_x"], y, value, color, fontsize=FONT_SIZES["annotation"], ha="right")

    ax.axhline(layout["divider_y"], xmin=0.03, xmax=0.97, color=theme["text"], alpha=0.25, linewidth=1)
    for y, label, key, color in zip(layout["levels_y"], ("FRZ:", "WB0:"), ("frz_pts_z", "wb_frz_pts_z"), (theme["frz"], theme["wb_frz"])):
        _draw_text(ax, layout["label_x"], y, label, theme["text"], fontsize=FONT_SIZES["annotation"], alpha=0.9)
        levels = ", ".join(f"{height:.0f} m" for height in winter[key]) or "NA"
        _fit_width(ax.text(layout["value_x"], y, levels, color=color, fontsize=FONT_SIZES["annotation"], weight="bold",
                           ha="right", transform=ax.transAxes), ax, layout["levels_width"])
    return ax


def _draw_even_column(ax, lines, theme, x_left, x_right):
    """Draw a column of title, row, check, note, and gap lines spread evenly from top to bottom."""
    # each line is one of
    #   ("title", text, level)          section title, with an optional HIGH / MOD / LOW / NONE level on the right
    #   ("row", label, value, color)    label on the left, value on the right
    #   ("check", met, label, value)    label, value, then ✓ / ✗ on the right
    #   ("note", text, n_lines)         faded citation taking n_lines line slots
    #   ("gap",)                        an empty line slot between sections
    top, bottom = HAZARDS_LAYOUT["top"], HAZARDS_LAYOUT["bottom"]
    slots = sum(line[2] if line[0] == "note" else 1 for line in lines)
    step = (top - bottom) / max(slots - 1, 1)
    small = FONT_SIZES["small"] + 0.5
    y = top
    for kind, *line in lines:
        if kind == "title":
            _section_title(ax, line[0], theme, x_left, y)
            if line[1]:
                _draw_text(ax, x_right, y, line[1], HAZARD_LEVEL_COLORS[line[1]], fontsize=FONT_SIZES["label"], ha="right")
        elif kind == "row":
            _draw_text(ax, x_left, y, line[0], theme["text"], fontsize=FONT_SIZES["label"], alpha=0.9)
            _draw_text(ax, x_right, y, line[1], line[2], fontsize=FONT_SIZES["label"], ha="right")
        elif kind == "check":
            met, label, value = line
            _draw_text(ax, x_left, y, label, theme["text"], fontsize=small, alpha=0.85)
            _draw_text(ax, x_right - 0.04, y, value, theme["text"], fontsize=small, alpha=0.85, ha="right")
            _draw_text(ax, x_right, y, "✓" if met else "✗", theme["moisture"] if met else theme["warm"],
                       fontsize=FONT_SIZES["label"], ha="right")
        elif kind == "note":
            # anchored by its top, so extra lines grow down into the slots below
            _fit_width(ax.text(x_left, y + 0.03, line[0], color=theme["text"], fontsize=FONT_SIZES["fine"] + 0.5,
                               alpha=0.45, weight="bold", va="top", linespacing=1.15, transform=ax.transAxes),
                       ax, x_right - x_left)
            y -= step * (line[1] - 1)
        y -= step


def _fzfg_lines(winter):
    """Freezing fog hazard lines: its level, four checks, and the source of its criteria."""
    fzfg = winter["fzfg"]
    ok, values = fzfg["ok"], fzfg["values"]
    return [("title", "FREEZING FOG", fzfg["level"]),
            ("check", ok["cold"], "Sfc T < 0°C", f"{values['cold']:+.1f}°C"),
            ("check", ok["sat"], "Sfc RH ≥ 95%", f"{values['sat']:.0f}%"),
            ("check", ok["calm"], "10 m wind ≤ 6 kt", f"{values['calm']:.0f} kt"),
            ("check", ok["inv"], "Sfc inversion",
             f"{_fmt_signed(values['inv'], '°C', 0)} / {winter['sfc_inv_depth'] / 1000:.1f} km" if values["inv"] > 0 else "none"),
            ("note", "Based on Gultepe et al. 2007", 1)]


def _fzdz_lines(winter):
    """Freezing drizzle hazard lines: its level, four checks (five with omega), and the sources of its criteria."""
    fzdz = winter["fzdz"]
    ok, cloud = fzdz["ok"], fzdz["cloud"]
    missing = "none" if cloud is None else "NA"
    lines = [("title", "FREEZING DRIZZLE", fzdz["level"]),
             ("check", ok["sfc"], "Sfc Tw ≤ 0°C", f"{winter['sfc_Tw']:+.1f}°C"),
             ("check", ok["top"], "Cloud top > −10°C", f"{cloud['top_T']:+.0f}°C" if cloud else missing),
             ("check", ok["dry"], "Dry 1 km above", f"RHᵢ {fzdz['rh_above']:.0f}%" if np.isfinite(fzdz["rh_above"]) else missing),
             ("check", ok["deep"], "Cloud ≥ 500 m deep", f"{fzdz['depth']:.0f} m" if np.isfinite(fzdz["depth"]) else missing)]
    if fzdz["lift"] is not None:
        lines.append(("check", fzdz["lift"], "Lift in cloud", f"ω {fzdz['omega']:+.1f}" if np.isfinite(fzdz["omega"]) else missing))
    lines.append(("note", "Based on Huffman & Norman 1988,\nRauber et al. 2000 & McDonough et al. 2004", 2))
    return lines


def _plot_hazards_table(ax, winter, theme):
    """Draw the hazards box: cold, blowing snow, and other on the left; freezing fog and drizzle checklists on the right."""
    layout = HAZARDS_LAYOUT
    _setup_text_box(ax, theme)

    left = [("title", "COLD", None),
            ("row", "Sfc T:", _fmt(winter["sfc_T_F"], "°F"), theme["temperature"]),
            ("row", "Wind chill:", _fmt(winter["wind_chill"], "°F"), theme["cool"]),
            ("row", "10 m wind:", _wind_text(*winter["sfc_wind_uv"]), theme["kinematics_strong"]),
            ("title", "BLOWING SNOW", None),
            ("row", "PBL Max Gust:", f"{_fmt(winter['gust_potential'])} kt", theme["kinematics_strong"]),
            ("row", "ML mean:", _wind_text(*winter["mean_wind_ml"]), theme["kinematics_strong"]),
            ("row", "ML depth:", _fmt(winter["ml_depth"], " m"), theme["lapse"]),
            ("row", "ML LR:", _fmt(winter["ml_lapse_rate"], "°C/km", 1), theme["lapse"]),
            ("note", "Based on Sczepanski et al. 2026", 1),
            ("title", "OTHER", None),
            ("row", "Snow squall:", _fmt(winter["snsq"], decimals=1), theme["composite_index"])]
    _draw_even_column(ax, left, theme, *layout["left"])
    ax.plot([layout["divider_x"]] * 2, [0.03, 0.97], color=theme["text"], alpha=0.25, linewidth=1, transform=ax.transAxes)
    _draw_even_column(ax, _fzfg_lines(winter) + [("gap",)] + _fzdz_lines(winter), theme, *layout["right"])
    return ax

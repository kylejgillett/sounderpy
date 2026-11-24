"""
    SOUNDERPY READOUTS

    Purpose of module:

    EXPERIMENTAL: interactive readouts for SounderPy figures on the web. A figure saved with ``readout=True`` gets a JSON
    file next to its PNG that maps the image's pixels to the data: where the skew-t and hodograph sit in
    the image, how to turn a pixel into pressure and temperature (or u and v), and the plotted profile.
    The browser script ``web/sounderpy-readout.js`` reads that file to show hover readouts and to let
    users drag the profile, and ``apply_profile_edits()`` turns an edited profile back into ``clean_data``
    so SounderPy can redraw it.

    readout_data()      : the readout dict for a finished figure
    save_readout()      : write the readout JSON for a saved figure
    apply_profile_edits(): rebuild clean_data from a profile edited in the browser
    readout_js_path()   : where the browser script is installed

    Docs:
        web readouts:   https://sounderpy.app/docs/developers/interactive_figures.html

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""

import copy
import json
from importlib.metadata import PackageNotFoundError, version as _package_version
from pathlib import Path

import numpy as np
import matplotlib as mpl
import metpy.calc as mpcalc
from metpy.units import units


# the readout JSON's format version; the browser script checks it
READOUT_FORMAT = "sounderpy-readout"
READOUT_VERSION = 1

# the hodograph is drawn up to 12 km AGL
HODOGRAPH_TOP_M = 12000

# a value that changes by less than this between the original and edited profile isn't counted as an edit
EDIT_TOLERANCE = 0.05


def readout_js_path():
    """Path to the browser script (``sounderpy-readout.js``) that reads readout files.

    Copy it next to your web page, or serve it from here, e.g. ``shutil.copy(spy.readout_js_path(), 'site/')``.

    :return: the script's path
    :rtype: pathlib.Path
    """
    return Path(__file__).parent / "web" / "sounderpy-readout.js"


def _sounderpy_version():
    try:
        return _package_version("sounderpy")
    except PackageNotFoundError:
        return None


def _plain(values, unit=None):
    """A 1-D list of floats (None for missing) from a Quantity, masked array, or list."""
    if values is None:
        return None
    if hasattr(values, "to") and unit is not None:
        values = values.to(unit)
    values = getattr(values, "magnitude", values)
    values = np.ma.filled(np.ma.asarray(values, dtype=float), np.nan).ravel()
    return [round(float(v), 2) if np.isfinite(v) else None for v in values]


def save_bbox(fig, pad_inches=None):
    """The region of the figure that ``savefig(bbox_inches='tight')`` keeps, in inches.

    Saving with this bbox gives the same image as 'tight', and tells the readout where the image starts.
    """
    if pad_inches is None:
        pad_inches = mpl.rcParams["savefig.pad_inches"]
    fig.draw_without_rendering()
    return fig.get_tightbbox(fig.canvas.get_renderer()).padded(pad_inches)


def save_dpi(fig):
    """The resolution savefig uses for this figure."""
    dpi = mpl.rcParams["savefig.dpi"]
    return fig.dpi if dpi == "figure" else float(dpi)


def image_size(bbox, dpi):
    """The saved image's (width, height) in pixels; Matplotlib drops partial pixels."""
    return int(bbox.width * dpi), int(bbox.height * dpi)


def _image_mapping(fig, bbox, dpi, height):
    """A function from display coordinates (Matplotlib pixels, origin bottom left) to saved-image pixels
    (origin top left)."""
    # display coordinates are at the figure's dpi; the saved image may use another
    scale = dpi / fig.dpi
    x0, y0 = bbox.x0 * fig.dpi, bbox.y0 * fig.dpi
    # rows count down from the top of the saved image

    def to_image(points):
        points = np.asarray(points, dtype=float)
        return np.column_stack([(points[:, 0] - x0) * scale, height - (points[:, 1] - y0) * scale])

    return to_image


def _fit_affine(forward, samples, to_image):
    """The 2x3 affine matrix taking data-space points to image pixels, fitted to sample points.

    The skew-t is affine in (temperature, log10 pressure) and the hodograph in (u, v), so the fit is exact.
    """
    samples = np.asarray(samples, dtype=float)
    pixels = to_image(forward(samples))
    design = np.column_stack([samples, np.ones(len(samples))])
    coef, *_ = np.linalg.lstsq(design, pixels, rcond=None)
    error = float(np.abs(design @ coef - pixels).max())
    if error > 0.5:
        raise ValueError(f"The panel's pixel mapping isn't affine (error {error:.2f} px); the readout would be misaligned.")
    return [[round(float(v), 6) for v in row] for row in coef.T]


def _bbox_pixels(ax, to_image):
    (x0, y0), (x1, y1) = to_image(ax.get_window_extent().get_points())
    return [round(float(min(x0, x1)), 2), round(float(min(y0, y1)), 2), round(float(max(x0, x1)), 2), round(float(max(y0, y1)), 2)]


def _skewt_panel(skew, to_image):
    ax = skew.ax
    bottom, top = ax.get_ylim()
    pressures = np.geomspace(max(min(bottom, top), 10), max(bottom, top), 4)
    samples = np.array([(t, np.log10(p)) for t in (-60, 0, 60) for p in pressures])
    forward = lambda pts: ax.transData.transform(np.column_stack([pts[:, 0], 10 ** pts[:, 1]]))
    return {"bbox": _bbox_pixels(ax, to_image), "transform": _fit_affine(forward, samples, to_image),
            "coords": ["T_degC", "log10_p_hPa"], "p_limits": [round(float(bottom), 2), round(float(top), 2)]}


def _hodograph_panel(hodo, to_image, offset):
    ax = hodo.ax
    samples = np.array([(u, v) for u in (-50, 0, 50) for v in (-50, 0, 50)], dtype=float)
    return {"bbox": _bbox_pixels(ax, to_image), "transform": _fit_affine(ax.transData.transform, samples, to_image),
            "coords": ["u_kt", "v_kt"], "offset": [round(float(offset[0]), 3), round(float(offset[1]), 3)],
            "max_height": HODOGRAPH_TOP_M}


def _profile(data, snd_params):
    """The plotted profile and the level-by-level values a readout shows."""
    p, z, T, Td = data["p"], data["z"], data["T"], data["Td"]
    u, v = data["u"].to("kt"), data["v"].to("kt")
    general = (snd_params or {}).get("general", {})

    def level_values(key, compute):
        value = general.get(key)
        if value is not None and np.shape(value) == np.shape(p):
            return value
        try:
            return compute()
        except Exception:
            return None

    # relative humidity in percent, as sounding_params stores it
    rh = level_values("rel_humidity", lambda: mpcalc.relative_humidity_from_dewpoint(T, Td) * 100)
    profile = {
        "p": _plain(p, "hPa"),
        "z": _plain(z, "m"),
        "z_agl": _plain(z - z[0], "m"),
        "T": _plain(T, "degC"),
        "Td": _plain(Td, "degC"),
        "u": _plain(u, "kt"),
        "v": _plain(v, "kt"),
        "wdir": _plain(mpcalc.wind_direction(u, v), "degree"),
        "wspd": _plain(mpcalc.wind_speed(u, v), "kt"),
        "Tw": _plain(level_values("wet_bulb", lambda: mpcalc.wet_bulb_temperature(p, T, Td)), "degC"),
        "Tv": _plain(level_values("virt_temp", lambda: mpcalc.virtual_temperature_from_dewpoint(p, T, Td)), "degC"),
        "rh": _plain(rh.to("dimensionless").m if rh is not None else None),
        "mr": _plain(level_values("mix_ratio", lambda: mpcalc.mixing_ratio_from_relative_humidity(p, T, rh / 100)), "g/kg"),
        "theta": _plain(level_values("theta", lambda: mpcalc.potential_temperature(p, T)), "K"),
        "thetae": _plain(level_values("theta_e", lambda: mpcalc.equivalent_potential_temperature(p, T, Td)), "K"),
    }
    if "omega" in data and np.shape(data["omega"]) == np.shape(p):
        profile["omega"] = _plain(data["omega"], "Pa/s")
    return {key: value for key, value in profile.items() if value is not None}


def _hodograph_trace(snd_params):
    """The ground-relative wind every 100 m up to the hodograph's top, for matching the drawn trace."""
    intrp = (snd_params or {}).get("intrp")
    if not intrp:
        return None
    z = np.asarray(intrp["zINTRP"], dtype=float)
    u = np.asarray(np.ma.filled(intrp["uINTRP"], np.nan), dtype=float)
    v = np.asarray(np.ma.filled(intrp["vINTRP"], np.nan), dtype=float)
    keep = np.isfinite(z) & np.isfinite(u) & np.isfinite(v) & (z <= HODOGRAPH_TOP_M)
    return {"z_agl": _plain(z[keep]), "u": _plain(u[keep]), "v": _plain(v[keep])}


# ---------------------------------------------------------------------------------------------------
# PARAMETERS IN THE FIGURE'S TABLES
# ---------------------------------------------------------------------------------------------------
PARCEL_NAMES = {"sb": "surface-based", "ml": "mixed-layer", "mu": "most-unstable"}
# table layers: key -> (bottom, top) in m AGL
LAYERS = {"0_to_500": (0, 500), "0_to_1000": (0, 1000), "1_to_3000": (1000, 3000), "3_to_6000": (3000, 6000),
          "6_to_9000": (6000, 9000), "0_to_3000": (0, 3000), "0_to_6000": (0, 6000), "0_to_10000": (0, 10000),
          "0_500": (0, 500), "0_1000": (0, 1000), "1_3000": (1000, 3000), "3_6000": (3000, 6000), "6_9000": (6000, 9000)}
KINEMATIC_PARAMS = {
    "shear": ("Bulk wind difference", "kt", "The magnitude of the vector difference between the winds at the layer's top and bottom."),
    "srh": ("Storm-relative helicity", "m²/s²", "Streamwise vorticity in the storm's inflow over the layer, for the selected storm motion."),
    "srw": ("Storm-relative wind", "kt", "Mean storm-relative wind speed over the layer."),
    "swv_perc": ("Streamwiseness of vorticity", "%", "The share of the layer's horizontal vorticity that is streamwise to the storm-relative inflow."),
    "swv": ("Streamwise vorticity", "/s", "Mean streamwise vorticity over the layer."),
}
SINGLE_PARAMS = {
    "thermo.dcape": ("DCAPE", "J/kg", "Downdraft CAPE: the energy available to a saturated downdraft descending to the surface.", None),
    "thermo.dcin": ("DCIN", "J/kg", "Downdraft CIN: the negative buoyancy a downdraft meets near the surface.", None),
    "thermo.mu_ncape": ("MU normalized CAPE", "m/s²", "Most-unstable CAPE divided by the depth of its positive area; how fast an updraft accelerates.", None),
    "kinem.eil_depth": ("Effective inflow layer depth", "m", "Depth of the layer with at least 100 J/kg CAPE and no more than 250 J/kg CIN.", "eil"),
    "thermo.k_index": ("K-index", "", "A thunderstorm potential index from 850-500 hPa lapse rate and 850-700 hPa moisture.", {"p": [850, 500]}),
    "thermo.t_totals": ("Total totals", "", "A severe weather index from the 850-500 hPa lapse rate and 850 hPa dewpoint.", {"p": [850, 500]}),
    "thermo.lr_03km": ("0-3 km lapse rate", "°C/km", "The temperature lapse rate from the surface to 3 km AGL.", {"z": [0, 3000]}),
    "thermo.lr_36km": ("3-6 km lapse rate", "°C/km", "The temperature lapse rate from 3 to 6 km AGL.", {"z": [3000, 6000]}),
    "thermo.lr_75": ("700-500 hPa lapse rate", "°C/km", "The mid-level temperature lapse rate.", {"p": [700, 500]}),
    "kinem.eil_scp": ("Supercell composite parameter", "", "Effective-layer SCP (Thompson et al. 2004): MUCAPE, effective SRH, and effective bulk shear.", "eil"),
    "kinem.eil_stp": ("Significant tornado parameter", "", "Effective-layer STP (Thompson et al. 2012): MLCAPE, MLLCL, effective SRH and shear, and MLCIN.", "eil"),
    "kinem.ehi_0_to_1000": ("0-1 km energy-helicity index", "", "CAPE and 0-1 km SRH combined.", {"z": [0, 1000]}),
    "kinem.ehi_0_to_3000": ("0-3 km energy-helicity index", "", "CAPE and 0-3 km SRH combined.", {"z": [0, 3000]}),
    "thermo.ship": ("Significant hail parameter", "", "SHIP: MUCAPE, mixing ratio, lapse rates, 500 hPa temperature, and shear.", None),
    "thermo.ml_brn": ("ML bulk Richardson number", "", "Mixed-layer CAPE over the 0-6 km shear term; low values favor supercells.", {"z": [0, 6000]}),
    "general.pwat": ("Precipitable water", "in", "The depth of water if all the column's vapor condensed.", None),
    "general.wet_bulb[0]": ("Surface wet-bulb temperature", "°C", "The temperature air cools to by evaporating water into it, at the surface.", None),
    "general.frz_pt_z": ("Freezing level", "m", "The lowest height (AGL) where the temperature reaches 0 °C.", "point"),
    "general.wb_frz_pt_z": ("Wet-bulb zero height", "m", "The lowest height (AGL) where the wet-bulb temperature reaches 0 °C.", "point"),
    "thermo.thetae_diff": ("θe difference", "K", "The drop in θe from its low-level maximum to its mid-level minimum (Cohen et al. 2007); larger favors strong outflow.", None),
    "thermo.pbl_depth": ("PBL depth", "m", "Depth of the planetary boundary layer.", "pbl"),
    "kinem.mean_wind_pbl": ("PBL mean wind", "kt", "The mean wind through the planetary boundary layer (the transport wind), from the direction it blows.", "pbl"),
    # winter sounding
    "winter.ptype": ("Surface precipitation type (Bourgouin/Birk)", "", "The precipitation type at the surface, from the melting and refreezing energy of the wet-bulb profile below the precipitation source (Bourgouin 2000; Birk et al. 2021).", None),
    "winter.sharppy_ptype": ("SPC best-guess precipitation type", "", "SHARPpy's best-guess precipitation type (Blumberg et al. 2017).", None),
    "winter.source_p": ("Precipitation source", "hPa", "Where precipitation is assumed to start: the highest saturated (RH > 80%) 50 hPa layer below 5 km AGL, in rising air when omega is available.", None),
    "winter.layer.kind": ("Layer", "", "A melting (warm) or refreezing (cold) layer below the precipitation source, from the wet-bulb profile.", None),
    "winter.layer.bottom": ("Layer base", "m", "The bottom of the layer (AGL).", None),
    "winter.layer.top": ("Layer top", "m", "The top of the layer (AGL).", None),
    "winter.layer.depth": ("Layer depth", "m", "The depth of the layer.", None),
    "winter.layer.extreme": ("Layer max/min wet-bulb", "°C", "The warmest wet-bulb temperature in a melting layer, or the coldest in a refreezing layer.", None),
    "winter.layer.energy": ("Layer energy", "J/kg", "The area between the wet-bulb profile and 0 °C: melting energy (positive) or refreezing energy (negative) (Bourgouin 2000).", None),
    "winter.dgz_z": ("Dendritic growth zone", "", "The -12 to -17 °C layer in cloud, where snow crystals grow fastest (Takahashi et al. 1991).", None),
    "winter.dgz_depth": ("DGZ depth", "m", "The depth of the dendritic growth zone.", None),
    "winter.dgz_rh_ice": ("DGZ humidity", "", "Mean relative humidity over water and over ice in the dendritic growth zone.", None),
    "winter.dgz_omega": ("DGZ omega", "µb/s", "Mean vertical motion in the dendritic growth zone; negative is rising air.", None),
    "winter.oprh": ("OPRH", "", "Omega × precipitable water × RH in the dendritic growth zone; negative favors snow (Blumberg et al. 2017).", None),
    "winter.slr": ("Kuchera snow-to-liquid ratio", ":1", "The Kuchera snow ratio, when snow reaches the surface.", None),
    "winter.frz_pts_z": ("Freezing levels", "", "Every 0 °C crossing of the temperature profile.", None),
    "winter.wb_frz_pts_z": ("Wet-bulb zero levels", "", "Every 0 °C crossing of the wet-bulb profile.", None),
    "winter.sfc_T_F": ("Surface temperature", "°F", "The temperature at the surface.", None),
    "winter.wind_chill": ("Wind chill", "°F", "Wind chill, at or below 50 °F with at least 3 mph of wind (NWS 2001).", None),
    "winter.sfc_wind_uv": ("10 m wind", "kt", "The surface wind: the lowest level with a wind, within 30 m of the surface, from the direction it blows.", None),
    "winter.gust_potential": ("PBL gust potential", "kt", "The strongest wind between the surface and the PBL top: what mixing could bring down.", "pbl"),
    "winter.mean_wind_pbl": ("PBL mean wind", "kt", "The mean wind through the planetary boundary layer, from the direction it blows.", "pbl"),
    "winter.mean_wind_dgz": ("DGZ mean wind", "kt", "The mean wind through the dendritic growth zone, which carries the growing snow.", "dgz"),
    "winter.max_wind_0_to_2000": ("0-2 km max wind", "kt", "The strongest wind in the lowest 2 km.", {"z": [0, 2000]}),
    "winter.pbl_depth": ("PBL depth", "m", "Depth of the planetary boundary layer, for blowing snow and mixing.", "pbl"),
    "winter.pbl_lapse_rate": ("PBL lapse rate", "°C/km", "The temperature lapse rate through the planetary boundary layer (Sczepanski et al. 2026).", "pbl"),
    "winter.snsq": ("Snow squall parameter", "", "Snow squall potential from low-level moisture, instability, and wind (Banacos et al. 2014).", None),
    "winter.fzfg": ("Freezing fog potential", "", "An ingredient check (Gultepe et al. 2007): HIGH when every check is met, MOD when one is short, LOW otherwise.", None),
    "winter.fzfg_check": ("Freezing fog check", "", "One ingredient of the freezing fog potential.", None),
    "winter.fzdz": ("Freezing drizzle potential", "", "An ingredient check on the lowest cloud (Huffman & Norman 1988; Rauber et al. 2000; McDonough et al. 2004): HIGH when every check is met, MOD when one is short, LOW otherwise.", None),
    "winter.fzdz_check": ("Freezing drizzle check", "", "One ingredient of the freezing drizzle potential.", None),
    "winter.level_winds": ("Wind", "kt", "The wind and temperature at a standard pressure level.", None),
}


def _describe_param(key):
    """(name, units, about, layer) for a parameter key, such as 'kinem.srh_0_to_3000'."""
    if key in SINGLE_PARAMS:
        return SINGLE_PARAMS[key]
    group, name = key.split(".", 1)
    if group == "thermo":
        for parcel, parcel_name in PARCEL_NAMES.items():
            label = parcel.upper()
            if name == f"{parcel}_ecape":
                return (f"{label} storm-relative ECAPE", "J/kg", f"Entraining CAPE of the {parcel_name} parcel, from the storm-relative inflow (Peters et al. 2023).", None)
            if name == f"{parcel}cape":
                return (f"{label}CAPE", "J/kg", f"Convective available potential energy of the {parcel_name} parcel.", None)
            if name == f"{parcel}6cape":
                return (f"{label} 0-6 km CAPE", "J/kg", f"CAPE of the {parcel_name} parcel below 6 km AGL.", {"z": [0, 6000]})
            if name == f"{parcel}3cape":
                return (f"{label} 0-3 km CAPE", "J/kg", f"CAPE of the {parcel_name} parcel below 3 km AGL.", {"z": [0, 3000]})
            if name == f"{parcel}cin":
                return (f"{label}CIN", "J/kg", f"Convective inhibition of the {parcel_name} parcel.", None)
            if name == f"{parcel}_lcl_z":
                return (f"{label} LCL height", "m", f"The lifting condensation level of the {parcel_name} parcel (AGL).", "point")
            if name == f"{parcel}_lfc_z":
                return (f"{label} LFC height", "m", f"The level of free convection of the {parcel_name} parcel (AGL).", "point")
            if name == f"{parcel}_el_z":
                return (f"{label} EL height", "m", f"The equilibrium level of the {parcel_name} parcel (AGL).", "point")
    if group == "general" and name[:3] in ("rh_", "w_0", "w_1", "w_3", "w_6"):
        layer = name.split("_", 1)[1]
        bottom, top = LAYERS.get(layer, (None, None))
        depth = f"{bottom / 1000:g}-{top / 1000:g} km"
        if name.startswith("rh_"):
            return (f"{depth} relative humidity", "%", f"Mean relative humidity from {depth} AGL.", {"z": [bottom, top]})
        return (f"{depth} mixing ratio", "g/kg", f"Mean water vapor mixing ratio from {depth} AGL.", {"z": [bottom, top]})
    if group == "kinem":
        for prefix in sorted(KINEMATIC_PARAMS, key=len, reverse=True):
            if name.startswith(prefix + "_"):
                layer = name[len(prefix) + 1:]
                title, units, about = KINEMATIC_PARAMS[prefix]
                if layer == "eil":
                    return (f"Effective {title[0].lower()}{title[1:]}", units, about + " Over the effective inflow layer.", "eil")
                bottom, top = LAYERS.get(layer, (None, None))
                if bottom is None:
                    break
                return (f"{bottom / 1000:g}-{top / 1000:g} km {title[0].lower()}{title[1:]}", units, about, {"z": [bottom, top]})
    return (name, "", "", None)


def _resolve_param(params, key):
    """A parameter's value from sounding_params by key, such as 'thermo.mlcape' or 'general.wet_bulb[0]'."""
    group, name = key.split(".", 1)
    index = None
    if name.endswith("]"):
        name, index = name[:-1].split("[")
    value = params[group][name]
    return value[int(index)] if index is not None else value


def _finite(value):
    try:
        value = float(np.ma.filled(np.ma.asarray(getattr(value, "magnitude", value), dtype=float), np.nan))
    except (TypeError, ValueError):
        return None
    return round(value, 4) if np.isfinite(value) else None


def _layer_heights(layer, key, params):
    """A parameter's layer as {'z': [bottom, top]} (m AGL), {'p': [bottom, top]} (hPa), or {'z_at': height}."""
    if layer is None or isinstance(layer, dict):
        return layer
    kinem, thermo = params.get("kinem", {}), params.get("thermo", {})
    if layer == "eil":
        bottom, top = (_finite(v) for v in kinem.get("eil_z", (None, None)))
        return {"z": [bottom, top]} if bottom is not None and top is not None else None
    if layer == "pbl":
        depth = _finite(thermo.get("pbl_depth"))
        return {"z": [0, depth]} if depth is not None else None
    if layer == "dgz":
        bottom, top = (_finite(v) for v in params.get("winter", {}).get("dgz_z", (None, None)))
        return {"z": [bottom, top]} if bottom is not None and top is not None else None
    if layer == "point":
        height = _finite(_resolve_param(params, key))
        return {"z_at": height} if height is not None else None
    return None


def _clean_layer(layer):
    """A layer from a tag with plain numbers, or None when any part is missing."""
    if not layer:
        return None
    clean = {}
    for kind, values in layer.items():
        values = [_finite(v) for v in np.atleast_1d(values)]
        if not values or any(v is None for v in values):
            return None
        clean[kind] = values if kind != "z_at" or len(values) > 1 else values[0]
    return clean


def _table_params(fig, params, to_image, renderer):
    """Every tagged value in the figure's tables: where it is in the image and what it means."""
    found = []
    for text in fig.findobj(lambda artist: hasattr(artist, "sounderpy_param")):
        if not text.get_visible() or not text.get_text().strip():
            continue
        # a key, or a dict with the key and anything the key can't look up (value, name, units, layer, extra rows)
        tag = text.sounderpy_param if isinstance(text.sounderpy_param, dict) else {"key": text.sounderpy_param}
        key = tag["key"]
        name, units, about, layer = _describe_param(key)
        try:
            raw = tag["value"] if "value" in tag else _resolve_param(params, key)
        except (KeyError, IndexError, TypeError, ValueError):
            continue
        extent = text.get_window_extent(renderer)
        (x0, y0), (x1, y1) = to_image(extent.get_points())
        entry = {"key": key, "name": tag.get("name", name), "units": tag.get("units", units), "about": tag.get("about", about),
                 "bbox": [round(float(min(x0, x1)), 1), round(float(min(y0, y1)), 1), round(float(max(x0, x1)), 1), round(float(max(y0, y1)), 1)]}
        if isinstance(raw, str):
            entry["value"], entry["kind"] = raw, "text"
        elif isinstance(raw, (tuple, list)) and len(raw) == 2:
            # a wind, as (u, v) in kt
            u, v = (_finite(c) for c in raw)
            entry["value"] = None if u is None or v is None else [u, v]
            entry["kind"] = "wind"
        else:
            entry["value"] = _finite(raw)
        try:
            entry["layer"] = _clean_layer(tag["layer"]) if "layer" in tag else _layer_heights(layer, key, params)
        except (KeyError, IndexError, TypeError, ValueError):
            entry["layer"] = None
        if tag.get("extra"):
            entry["extra"] = [[str(label), str(value)] for label, value in tag["extra"]]
        found.append(entry)
    return found


def _match_length(values, length):
    """Values as floats, cut or padded with NaN to ``length``."""
    values = np.ma.filled(np.ma.asarray(getattr(values, "magnitude", values), dtype=float), np.nan).ravel()[:length]
    return np.concatenate([values, np.full(length - len(values), np.nan)])


def _height_panels(fig, to_image):
    """The figure's 0-3 km height panels: where each is, its height axis, and its values every 100 m."""
    panels = []
    for ax in fig.axes:
        info = getattr(ax, "sounderpy_panel", None)
        if not info:
            continue
        z = np.ma.filled(np.ma.asarray(info["z"], dtype=float), np.nan)
        bottom, top = ax.get_ylim()
        keep = np.isfinite(z) & (z <= max(bottom, top) + 1)
        samples = np.array([(x, h) for x in ax.get_xlim() for h in (bottom, top)], dtype=float)
        panels.append({
            "title": info["title"],
            "bbox": _bbox_pixels(ax, to_image),
            "transform": _fit_affine(ax.transData.transform, samples, to_image),
            "coords": ["x", "z_agl_m"],
            "z_limits": [round(float(bottom), 1), round(float(top), 1)],
            "z": _plain(z[keep]),
            "series": [{"label": label, "units": units, "values": _plain(_match_length(values, len(z))[keep])}
                       for label, units, values in info["series"]],
        })
    return panels


def readout_data(fig, axes, source, bbox, dpi=None, size=None):
    """The readout dict for a finished SounderPy figure.

    :param fig: the figure
    :param axes: the figure's panels, as returned with ``return_fig=True``
    :param source: what the figure plotted: ``{'kind', 'data', 'params', 'hodo_offset'}``
    :param bbox: the saved region of the figure in inches (see ``save_bbox``)
    :param dpi: the saved image's dpi. Default is the figure's savefig dpi
    :param size: the saved image's (width, height) in pixels. Default is the size Matplotlib should give it
    :return: a JSON-ready dict
    """
    dpi = save_dpi(fig) if dpi is None else dpi
    size = image_size(bbox, dpi) if size is None else tuple(size)
    to_image = _image_mapping(fig, bbox, dpi, size[1])
    data, params = source["data"], source.get("params")
    info = data.get("site_info", {})

    panels = {}
    if "skewt" in axes:
        panels["skewt"] = _skewt_panel(axes["skewt"], to_image)
    if "hodograph" in axes:
        panels["hodograph"] = _hodograph_panel(axes["hodograph"], to_image, source.get("hodo_offset", (0.0, 0.0)))

    readout = {
        "format": READOUT_FORMAT,
        "version": READOUT_VERSION,
        "sounderpy": _sounderpy_version(),
        "figure": source["kind"],
        "image": {"width": int(size[0]), "height": int(size[1])},
        # a dark-mode figure gets a dark readout
        "dark": bool(np.dot(mpl.colors.to_rgb(fig.get_facecolor()), (0.299, 0.587, 0.114)) < 0.5),
        "site": {"id": info.get("site-id"), "name": info.get("site-name"), "source": info.get("source"),
                 "valid": info.get("valid-time"), "latlon": _plain(info.get("site-latlon")) if info.get("site-latlon") is not None else None},
        "edited": info.get("edited"),
        "units": {"p": "hPa", "z": "m", "z_agl": "m", "T": "degC", "Td": "degC", "Tw": "degC", "Tv": "degC",
                  "rh": "%", "mr": "g/kg", "theta": "K", "thetae": "K", "u": "kt", "v": "kt", "wdir": "deg",
                  "wspd": "kt", "omega": "Pa/s"},
        "panels": panels,
        "profile": _profile(data, params),
    }
    trace = _hodograph_trace(params)
    if trace is not None and "hodograph" in panels:
        readout["hodograph_trace"] = trace
    # the storm motion (ground relative, kt) the figure's SRH and storm-relative values use
    kinem = (params or {}).get("kinem", {})
    storm_motion = [_finite(kinem.get("sm_u")), _finite(kinem.get("sm_v"))]
    if None not in storm_motion:
        readout["storm_motion"] = storm_motion
    if params:
        renderer = fig.canvas.get_renderer()
        readout["params"] = _table_params(fig, params, to_image, renderer)
        readout["height_panels"] = _height_panels(fig, to_image)
    return readout


def readout_filename(filename):
    """The readout file's name for a figure saved as ``filename``: ``kokc.png`` -> ``kokc.readout.json``."""
    path = Path(str(filename))
    stem = path.with_suffix("") if path.suffix.lower() in (".png", ".jpg", ".jpeg", ".svg", ".pdf", ".webp") else path
    return Path(f"{stem}.readout.json")


def saved_image_size(filename):
    """The (width, height) of a saved raster figure, or None for a vector one (or a file that can't be read)."""
    path = Path(str(filename))
    if not path.suffix:
        # savefig adds the default format's extension to a name without one
        path = path.with_name(f"{path.name}.{mpl.rcParams['savefig.format']}")
    try:
        from PIL import Image
        with Image.open(path) as image:
            return image.size
    except Exception:
        return None


def save_readout(readout, filename):
    """Write a readout dict as JSON next to the figure saved as ``filename``; returns the JSON's path."""
    path = readout_filename(filename)
    with open(path, "w", encoding="utf-8") as file:
        json.dump(readout, file, separators=(",", ":"), allow_nan=False)
    return path


def _as_array(values, name, length):
    if values is None:
        raise ValueError(f"The edited profile needs '{name}'.")
    array = np.array([np.nan if v is None else v for v in values], dtype=float)
    if array.shape != (length,):
        raise ValueError(f"The edited profile's '{name}' has {array.size} values; 'p' has {length}.")
    return array


def _describe_edits(original, edited, p):
    """A short note of what changed, e.g. 'T AT 2 LEVELS, WIND AT 1 LEVEL (1000-850 hPa)'."""
    changed_levels = np.zeros(len(p), dtype=bool)
    parts = []
    for label, keys in (("T", ("T",)), ("TD", ("Td",)), ("WIND", ("u", "v"))):
        changed = np.zeros(len(p), dtype=bool)
        for key in keys:
            a, b = original[key], edited[key]
            changed |= (np.abs(a - b) > EDIT_TOLERANCE) | (np.isfinite(a) != np.isfinite(b))
        if changed.any():
            count = int(changed.sum())
            parts.append(f"{label} AT {count} LEVEL{'S' if count > 1 else ''}")
        changed_levels |= changed
    if not parts:
        return None
    layer = p[changed_levels]
    bottom, top = float(np.nanmax(layer)), float(np.nanmin(layer))
    where = f"{bottom:.0f} hPa" if round(bottom) == round(top) else f"{bottom:.0f}-{top:.0f} hPa"
    return f"{', '.join(parts)} ({where})"


def apply_profile_edits(clean_data, edits):
    """Rebuild ``clean_data`` from a profile edited in the browser.

    ``edits`` is the profile the readout script returns from ``getProfile()``: lists of ``p`` (hPa), ``z`` (m),
    ``T`` and ``Td`` (°C), and ``u`` and ``v`` (kt), one value per plotted level. The result has those levels
    and values, the original's site information and titles, and any other per-level values (such as omega)
    interpolated to the edited levels. Its ``site_info['edited']`` describes what changed, and SounderPy's
    figures show that in an 'EDITED PROFILE' note.

    :param clean_data: the profile the readout was made from
    :type clean_data: dict, required
    :param edits: the edited profile, from the browser
    :type edits: dict, required
    :return: the edited profile, as ``clean_data``
    :rtype: dict
    """
    p = np.array([np.nan if v is None else v for v in edits.get("p", [])], dtype=float)
    if p.ndim != 1 or p.size < 3 or not np.all(np.isfinite(p)) or np.any(np.diff(p) >= 0):
        raise ValueError("The edited profile's pressures must be at least 3 values that decrease with height.")
    n = p.size
    new = {key: _as_array(edits.get(key), key, n) for key in ("z", "T", "Td", "u", "v")}
    if not np.all(np.isfinite(new["z"])) or np.any(np.diff(new["z"]) <= 0):
        raise ValueError("The edited profile's heights must increase with height.")
    # a dewpoint above the temperature isn't physical
    new["Td"] = np.fmin(new["Td"], new["T"])

    # the original's values at the edited levels, to describe what changed and to carry other fields over
    orig_p = np.asarray(clean_data["p"].to("hPa").m, dtype=float)
    order = np.argsort(np.log(orig_p))
    log_p = np.log(p)

    def at_edited_levels(values, unit=None):
        values = values.to(unit).m if unit else getattr(values, "m", values)
        values = np.ma.filled(np.ma.asarray(values, dtype=float), np.nan)
        result = np.interp(log_p, np.log(orig_p)[order], values[order], left=np.nan, right=np.nan)
        # levels that match an original level keep its exact value
        for i, level in enumerate(p):
            # the browser's pressures are rounded to 0.01 hPa
            match = np.flatnonzero(np.isclose(orig_p, level, rtol=0, atol=0.006))
            if match.size:
                result[i] = values[match[0]]
        return result

    original = {"T": at_edited_levels(clean_data["T"], "degC"), "Td": at_edited_levels(clean_data["Td"], "degC"),
                "u": at_edited_levels(clean_data["u"], "kt"), "v": at_edited_levels(clean_data["v"], "kt")}

    edited = {key: value for key, value in clean_data.items() if key not in ("p", "z", "T", "Td", "u", "v")}
    edited = copy.deepcopy(edited)
    edited.update({"p": p * units.hPa, "z": new["z"] * units.m, "T": new["T"] * units.degC, "Td": new["Td"] * units.degC,
                   "u": new["u"] * units.kt, "v": new["v"] * units.kt})
    # other per-level values follow the profile to its edited levels
    for key, values in clean_data.items():
        if key in edited and key not in ("p", "z", "T", "Td", "u", "v") and np.shape(values) == np.shape(orig_p):
            unit = getattr(values, "units", None)
            moved = at_edited_levels(values)
            edited[key] = moved * unit if unit is not None else moved

    info = dict(edited.get("site_info") or {})
    note = _describe_edits(original, new, p)
    if note:
        info["edited"] = note
    edited["site_info"] = info
    return edited

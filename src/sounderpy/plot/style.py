"""
    SOUNDERPY PLOT STYLE

    Purpose of module:

    House the look of SounderPy figures in one place: themes, font sizes, units,
    layouts, and watermark settings. Drawing tools are in components.py.

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""

# STANDARD PYTHON
from importlib.metadata import PackageNotFoundError, version as package_version
from pathlib import Path



#########################################################################
############################### THEMES ##################################
#########################################################################

PARCEL_COLORS = {"sb": "#c62828",
                 "ml": "#d95f02",
                 "mu": "#8828c3"}

# hodograph trace layers: bottom (m), top (m), and color
HODOGRAPH_LAYERS = ((0, 1000, "purple"),
                    (1000, 3000, "red"),
                    (3000, 6000, "darkorange"),
                    (6000, 9000, "gold"),
                    (9000, 12000, "#fff09f"))

# colors shared by the light and dark themes
BASE_THEME = {# skew-t profiles
              "temperature": "red",
              "dewpoint": "green",
              "virtual": "darkred",
              "wetbulb": "#3d8aff",

              # skew-t background
              "adiabat": "cornflowerblue",
              "isotherm": "blue",
              "shade": "gray",

              # parcels
              "parcel": PARCEL_COLORS,
              "special_parcel": "k",
              "secondary_parcel": "#808080",
              "downdraft": "purple",

              # hodograph
              "hodo_layers": HODOGRAPH_LAYERS,
              "height_marker": "white",
              "height_marker_edge": "black",
              "mcs": "orange",
              "dtm": "brown",

              # skew-t annotations
              "pbl": "gray",
              "eil": "lightblue",
              "dgz": "blue",
              "hgz": "green",
              "lapse_max": "firebrick",

              # shared warm and cool values
              "warm": "red",
              "cool": "cornflowerblue",
              "ascent": "firebrick",

              # tables
              "moisture": "forestgreen",
              "moisture_strong": "darkgreen",
              "kinematics": "#6495ED",
              "kinematics_strong": "#4169E1",
              "kinematics_layers": ("#1f9ad6", "#3a86e0", "#4f74e0", "#6463dc", "#7558d0"),
              "ncape": "brown",
              "lapse": "saddlebrown",
              "composite_index": "firebrick",
              "ehi": "orangered",

              # accessory panels
              "bar_edge": "black",
              "pressure_line": "gray",
              "low_fill": "cornflowerblue",
              "omega_zero": "b",
              "omega_reference": "purple",
              "zero_line": "C0",
              "streamwise": "orange",
              "theta": "purple",
              "u_wind": "forestgreen",
              "v_wind": "cornflowerblue",
              "cin": "teal",
              "cape": "tab:orange",

              # watermark
              "watermark": "cornflowerblue"}

LIGHT_THEME = {**BASE_THEME,
               "text": "black",
               "background": "white",
               "border": "black",
               "marker": "black",
               "barb": "black",
               "skew_line": "black",
               "shade_alpha": 0.02}

DARK_THEME = {**BASE_THEME,
              "text": "white",
              "background": "#161616",
              "border": "white",
              "marker": "white",
              "barb": "white",
              "skew_line": "white",
              "shade_alpha": 0.04,

              # the browns in BASE_THEME are too dark for this background
              "ncape": "#e08a6a",
              "lapse": "#d9a066"}


def _get_plot_theme(dark_mode=False, colorblind=False, overrides=None):
    # build a figure's theme: light or dark, then any changes from the user

    theme = dict(DARK_THEME if dark_mode else LIGHT_THEME)

    if colorblind:
        theme["dewpoint"] = "cornflowerblue"

    if overrides is not None:
        if not isinstance(overrides, dict):
            raise ValueError(f"'theme' must be a dict of colors to change, such as {{'border': 'k'}}. Got: {overrides}")

        invalid = set(overrides) - set(theme)
        if invalid:
            raise ValueError(f"Invalid theme setting(s): {sorted(invalid)}. Valid settings are: {sorted(theme)}")

        # parcel colors can be changed one at a time
        parcel = {**theme["parcel"], **overrides.get("parcel", {})}
        theme.update(overrides)
        theme["parcel"] = parcel

    return theme



#########################################################################
############################# FONT SIZES ################################
#########################################################################

# the size of each kind of text. use font_scale to resize a whole figure
FONT_SIZES = {"fine": 7,
              "tiny": 8,
              "small": 10,
              "note": 11,
              "label": 12,
              "tick": 13,
              "annotation": 13.5,
              "annotation_large": 14,
              "table": 15,
              "table_large": 16,
              "watermark_sub": 13,
              "watermark": 17,
              "title_small": 22,
              "subtitle": 23,
              "title": 30}



#########################################################################
################################ UNITS ##################################
#########################################################################

UNIT_JKG = r"$\mathrm{\mathbf{Jkg^{-1}}}$"
UNIT_M = r"$\mathrm{\mathbf{m}}$"
UNIT_KT = r"$\mathrm{\mathbf{kt}}$"
UNIT_M2S2 = r"$\mathrm{\mathbf{m^2s^{-2}}}$"
UNIT_PER_S = r"$\mathrm{\mathbf{s^{-1}}}$"
UNIT_GKG = r"$\mathrm{\mathbf{gkg^{-1}}}$"
UNIT_DEGC = r"$\mathrm{\mathbf{^{\circ}C}}$"
UNIT_LAPSE = r"$\mathrm{\mathbf{^{\circ}C}}$"



#########################################################################
########################## FIGURE LAYOUTS ###############################
#########################################################################

# the pressure (hPa) at the bottom of every skew-t
SKEWT_BOTTOM_PRESSURE = 1050

# axis positions for the full sounding figure, as figure fractions
FULL_SOUNDING_LAYOUT = {"skewt": (0.1124, 0.1005, 0.60, 0.85),
                        "hodograph": (0.53, 0.1003, 0.85, 0.85),
                        "map": (0.7005, 0.10, 0.20, 0.20),
                        "temp_adv": (0.701, 0.1003, 0.04, 0.85),
                        "omega": (0.127, 0.1003, 0.3, 0.85),
                        "streamwiseness": (0.945, -0.13, 0.065, 0.23),
                        "vorticity": (1.01, -0.13, 0.066, 0.23),
                        "sr_wind": (1.0755, -0.13, 0.066, 0.23),
                        "stability": (1.141, -0.13, 0.0653, 0.23),
                        "thermo_table": (0.124, -0.08, 0.406, 0.18),
                        "moisture_table": (0.53, -0.13, 0.15, 0.23),
                        "kinematics_table": (0.68, -0.13, 0.265, 0.23),
                        "indices_box": (1.106, 0.10, 0.10, 0.105),
                        "logo": (0.115, -0.135, 0.05, 0.05),
                        "margin": (0.1136, -0.1454, 1.1064, 1.1639)}

# axis positions for the hodograph and VAD figures, as figure fractions
FULL_HODOGRAPH_LAYOUT = {"hodograph": (0.13, 0.11, 0.77, 0.77),
                         "map": (0.2263, 0.11, 0.12, 0.16),
                         "streamwiseness": (0.8015, 0.11, 0.095, 0.32),
                         "vorticity": (0.8965, 0.11, 0.095, 0.32),
                         "sr_wind": (0.9915, 0.11, 0.095, 0.32),
                         "stability": (1.0865, 0.11, 0.095, 0.32),
                         "parameter_box": (0.8015, 0.43, 0.38, 0.45),
                         "logo": (0.352, 0.115, 0.042, 0.048),
                         "logo_no_map": (0.2303, 0.115, 0.042, 0.048),
                         "margin": (0.215, 0.095, 0.97, 0.865)}

# positions for each style of composite figure, as figure fractions. the table and everything below it
# grow down from table_top with the number of profiles, and "logo" is its (left, width)
COMPOSITE_LAYOUTS = {"full": {"figsize": (22, 13),
                              "skewt": (0.1124, 0.1005, 0.60, 0.85),
                              "hodograph": (0.53, 0.1003, 0.85, 0.85),
                              "table_left": 0.124,
                              "table_width": 1.082,
                              "table_top": 0.1005,
                              "table_row_height": 0.03,
                              "title_x": (0.1237, 1.206),
                              "subtitle_fontsize": FONT_SIZES["subtitle"],
                              "credits_x": (0.168, 0.525),
                              "logo": (0.115, 0.05),
                              "margin_x": (0.1136, 1.22)},
                     "sounding": {"figsize": (13.4, 13),
                                  "skewt": (0.0299, 0.1005, 0.9478, 0.85),
                                  "table_left": 0.0299,
                                  "table_width": 0.9478,
                                  "table_top": 0.1005,
                                  "table_row_height": 0.03,
                                  "title_x": (0.0299, 0.9777),
                                  "subtitle_fontsize": FONT_SIZES["table_large"],
                                  "credits_x": (0.1022, 0.97),
                                  "logo": (0.0157, 0.0821),
                                  "margin_x": (0.0127, 1.0)},
                     "hodograph": {"figsize": (11.85, 13),
                                   "hodograph": (0.0338, 0.1005, 0.9325, 0.85),
                                   "table_left": 0.0338,
                                   "table_width": 0.9325,
                                   "table_top": 0.1005,
                                   "table_row_height": 0.03,
                                   "title_x": (0.0338, 0.9662),
                                   "subtitle_fontsize": FONT_SIZES["table_large"],
                                   "credits_x": (0.1156, 0.96),
                                   "logo": (0.0177, 0.0928),
                                   "margin_x": (0.0143, 0.9915)}}



#########################################################################
########################### TABLE LAYOUTS ###############################
#########################################################################

# shown in tables for values that are missing
MISSING_VALUE = "NA"

# text positions inside each table, as fractions of the table
THERMO_TABLE_LAYOUT = {"label_x": 0.015,
                       "columns": [0.113, 0.273, 0.433, 0.591, 0.741, 0.889],
                       "header_columns": [0.113, 0.296, 0.433, 0.588, 0.750, 0.890],
                       "header_y": 0.875,
                       "rows_y": [0.735, 0.595, 0.455]}

MOISTURE_TABLE_LAYOUT = {"label_x": 0.067,
                         "columns": [0.367, 0.600],
                         "header_columns": [0.282, 0.667],
                         "header_y": 0.870,
                         "rows_y": [0.739, 0.609, 0.478, 0.348, 0.217]}

KINEMATICS_TABLE_LAYOUT = {"label_x": 0.034,
                           "label_alpha": 0.9,
                           "header_alpha": 0.8,
                           "columns": [0.196, 0.336, 0.555, 0.694, 0.811],
                           "header_columns": [0.208, 0.386, 0.554, 0.682, 0.863],
                           "header_y": 0.870,
                           "rows_y": [0.739, 0.609, 0.478, 0.348, 0.217, 0.087]}

HODOGRAPH_KINEMATICS_LAYOUT = {"label_x": 0.036,
                               "label_alpha": 0.9,
                               "header_alpha": 0.8,
                               "columns": [0.193, 0.336, 0.539, 0.688, 0.800],
                               "header_columns": [0.200, 0.370, 0.540, 0.670, 0.825],
                               "header_y": 0.911,
                               "rows_y": [0.844, 0.778, 0.711, 0.644, 0.578, 0.511]}

# the first four columns of the parcel table
HODOGRAPH_PARCEL_LAYOUT = {"label_x": 0.036,
                           "columns": [0.155, 0.391, 0.614, 0.828],
                           "header_columns": [0.155, 0.420, 0.614, 0.828],
                           "header_y": 0.430,
                           "rows_y": [0.363, 0.296, 0.230]}



#########################################################################
############################## COMPOSITE ################################
#########################################################################

# the part of viridis, the default colormap, that stands out from each background (at least 3:1 contrast)
COMPOSITE_VIRIDIS_RANGES = {"light": (0, 0.62),
                            "dark": (0.35, 1)}

# hodograph trace layers: bottom (m), top (m), and alpha
COMPOSITE_HODOGRAPH_LAYERS = ((0, 1000, 1),
                              (1000, 3000, 0.8),
                              (3000, 6000, 0.6),
                              (6000, 9000, 0.4))

# composite table columns: header, group of sounding_params, key, unit, and decimals
COMPOSITE_PARAMETERS = (("MLCAPE", "thermo", "mlcape", " {jkg}", None),
                        ("MUCAPE", "thermo", "mucape", " {jkg}", None),
                        ("SR-ECAPE", "thermo", "mu_ecape", " {jkg}", None),
                        ("MLCIN", "thermo", "mlcin", " {jkg}", None),
                        ("ML LCL", "thermo", "ml_lcl_z", " {m}", None),
                        ("0-1 SRH", "kinem", "srh_0_to_1000", "{m2s2}", None),
                        ("0-3 SRH", "kinem", "srh_0_to_3000", "{m2s2}", None),
                        ("0-6 BWD", "kinem", "shear_0_to_6000", "{kt}", None),
                        ("PWAT", "general", "pwat", "in", 2))

# each composite style's table: the headers of its parameter columns, and column positions as fractions of the table width.
# a short location is the site id alone, when there is one, and fontsize replaces the usual table font
COMPOSITE_TABLES = {"full": {"parameters": ("MLCAPE", "MUCAPE", "SR-ECAPE", "MLCIN", "ML LCL", "0-1 SRH", "0-3 SRH", "0-6 BWD", "PWAT"),
                             "letter": 0.006,
                             "profile": 0.024,
                             "location": 0.135,
                             "valid": 0.305,
                             "columns": (0.415, 0.482, 0.549, 0.618, 0.682, 0.742, 0.802, 0.868, 0.928)},
                    "sounding": {"parameters": ("MLCAPE", "MUCAPE", "SR-ECAPE", "MLCIN", "ML LCL", "PWAT"),
                                 "short_location": True,
                                 "fontsize": FONT_SIZES["tick"],
                                 "letter": 0.008,
                                 "profile": 0.03,
                                 "location": 0.115,
                                 "valid": 0.205,
                                 "columns": (0.37, 0.475, 0.58, 0.685, 0.79, 0.895)},
                    "hodograph": {"parameters": ("0-1 SRH", "0-3 SRH", "0-6 BWD"),
                                  "letter": 0.008,
                                  "profile": 0.035,
                                  "location": 0.145,
                                  "valid": 0.42,
                                  "columns": (0.64, 0.765, 0.885)}}



#########################################################################
############################## TIMESERIES ###############################
#########################################################################

# preset panels: the parameters each one plots, by their timeseries_params() names
TIMESERIES_PANELS = {"thermodynamics": ("sbcape", "mlcape", "mucape", "ml3cape", "mlcin"),
                     "kinematics": ("srh_0_to_1000", "srh_0_to_3000", "srh_eil", "shear_0_to_1000", "shear_0_to_3000", "shear_0_to_6000"),
                     "composites": ("eil_stp", "ehi_0_to_1000", "eil_scp"),
                     "low_level": ("swv_0_to_500", "srw_0_to_500"),
                     "surface": ("sfc_ws", "sfc_T", "sfc_Td")}
TIMESERIES_DEFAULT_PANELS = ("thermodynamics", "kinematics", "composites", "surface")

# timeseries figure layout. the figure grows with the number of panels, so spacing above and below
# the panels is in inches: titles sit the same distance above the panels as on the hodograph figure,
# and the logo and watermark sit below the time axis, lined up with the panels' left edge
TIMESERIES_LAYOUT = {"width": 16,               # figure width (in)
                     "panel_height": 3.4,       # height of each panel (in)
                     "panel_x": (0.08, 0.91),   # left and right edges of the panels (figure fraction)
                     "panel_gap": 0.12,         # space between panels (fraction of a panel)
                     "top": 1.25,               # space above the panels for the titles and padding (in)
                     "bottom": 1.45,            # space below the panels for the time axis, watermark, and padding (in)
                     "title_y": (0.6, 0.2),     # title and subtitle baselines above the panels (in)
                     "logo": (0.22, 0.6),       # logo bottom and height (in)
                     "credits_y": (0.56, 0.3),  # watermark title and subtitle baselines (in)
                     "credits_fontsize": ("table_large", "label")}  # watermark title and subtitle FONT_SIZES

# how each parameter is drawn: label, color (a theme setting or a color), style ("-", "--", ":", or "bar"),
# and optionally which axis it shares, for parameters with the same unit but very different ranges
# ("low_cape" puts CIN and 0-3 km CAPE on one axis, apart from total CAPE).
# other parameters are drawn as lines, labeled by name, in the next color of the default cycle
TIMESERIES_PARAMS = {"sbcape": ("SBCAPE", "parcel.sb", "-"),
                     "mlcape": ("MLCAPE", "parcel.ml", "-"),
                     "mucape": ("MUCAPE", "parcel.mu", "-"),
                     "sb_ecape": ("SB-ECAPE", "parcel.sb", "--"),
                     "ml_ecape": ("ML-ECAPE", "parcel.ml", "--"),
                     "mu_ecape": ("MU-ECAPE", "parcel.mu", "--"),
                     "sb3cape": ("0-3 km SBCAPE", "parcel.sb", ":", "low_cape"),
                     "ml3cape": ("0-3 km MLCAPE", "parcel.ml", ":", "low_cape"),
                     "mu3cape": ("0-3 km MUCAPE", "parcel.mu", ":", "low_cape"),
                     "sbcin": ("SBCIN", "parcel.sb", "bar", "low_cape"),
                     "mlcin": ("MLCIN", "parcel.ml", "bar", "low_cape"),
                     "mucin": ("MUCIN", "parcel.mu", "bar", "low_cape"),
                     "srh_0_to_500": ("0-500 m SRH", "mediumpurple", "-"),
                     "srh_0_to_1000": ("0-1 km SRH", "purple", "-"),
                     "srh_0_to_3000": ("0-3 km SRH", "red", "-"),
                     "srh_eil": ("Effective SRH", "kinematics_strong", "-"),
                     "shear_0_to_1000": ("0-1 km BWD", "purple", "--"),
                     "shear_0_to_3000": ("0-3 km BWD", "red", "--"),
                     "shear_0_to_6000": ("0-6 km BWD", "darkorange", "--"),
                     "ebwd": ("Effective BWD", "kinematics", "--"),
                     "srw_0_to_500": ("0-500 m SR wind", "dodgerblue", "-"),
                     "swv_0_to_500": ("0-500 m streamwise vorticity", "cool", "-"),
                     "eil_stp": ("STP", "composite_index", "-"),
                     "eil_scp": ("SCP", "theta", "-", "scp"),
                     "ehi_0_to_1000": ("0-1 km EHI", "ehi", "--"),
                     "ehi_0_to_3000": ("0-3 km EHI", "ehi", ":"),
                     "lr_03km": ("0-3 km lapse rate", "lapse", "-"),
                     "lr_36km": ("3-6 km lapse rate", "lapse", "--"),
                     "pwat": ("PWAT", "moisture", "-"),
                     "sfc_ws": ("Surface wind speed", "#90ee90", "bar"),
                     "sfc_T": ("Surface temperature", "temperature", "-"),
                     "sfc_Td": ("Surface dewpoint", "dewpoint", "-")}

# axis labels for the table's units
TIMESERIES_UNITS = {"J/kg": UNIT_JKG, "m2/s2": UNIT_M2S2, "kt": UNIT_KT, "1/s": UNIT_PER_S, "g/kg": UNIT_GKG,
                    "degC": UNIT_DEGC, "degC/km": r"$\mathrm{\mathbf{^{\circ}C\,km^{-1}}}$", "m": UNIT_M,
                    "hPa": r"$\mathrm{\mathbf{hPa}}$", "in": r"$\mathrm{\mathbf{in}}$", "%": "%", "deg": r"$\mathrm{\mathbf{^{\circ}}}$"}




#########################################################################
############################## WATERMARK ################################
#########################################################################

LOGO_PATH = Path(__file__).with_name("sounderpy_logo.png")

# the installed package version, which comes from pyproject.toml
try:
    SOUNDERPY_VERSION = f"v{package_version('sounderpy')}"
except PackageNotFoundError:
    SOUNDERPY_VERSION = None

# default watermark text for each figure style
WATERMARK_TEXT = {"full": ("SOUNDERPY VERTICAL PROFILE ANALYSIS TOOL",
                           " · ".join(part for part in (f"Built with {SOUNDERPY_VERSION}", "Make your own at sounderpy.app", "(Gillett, 2024)") if part)),
                  "compact": ("SOUNDERPY", "sounderpy.app")}

# a title or subtitle of None uses the default text for that figure's style
# a color of None uses the color of the figure's theme
WATERMARK_DEFAULTS = {"show": True,
                      "logo": LOGO_PATH,
                      "color": None,
                      "alpha": 0.6,
                      "title": None,
                      "subtitle": None}

WATERMARK = dict(WATERMARK_DEFAULTS)


def set_watermark(**settings):
    r"""Change the watermark on SounderPy figures, for this Python session

    Calling this with no arguments restores the default watermark.

    :param show: whether to show the watermark text and logo
    :type show: bool, optional
    :param title: the top line of text, on every figure. By default, each figure uses its own text
    :type title: str, optional
    :param subtitle: the bottom line of text, on every figure. By default, each figure uses its own text
    :type subtitle: str, optional
    :param logo: a path to an image to use as the logo, or ``None`` for no logo
    :type logo: str or None, optional
    :param color: the color of the text
    :type color: str, optional
    :param alpha: the transparency of the text
    :type alpha: float, optional
    """
    invalid = set(settings) - set(WATERMARK_DEFAULTS)
    if invalid:
        raise ValueError(f"Invalid watermark setting(s): {sorted(invalid)}. Valid settings are: {list(WATERMARK_DEFAULTS)}")
    WATERMARK.update(settings if settings else WATERMARK_DEFAULTS)
